"""Null-model controls transform one declared input and stay comparable to their source arm."""
from collections import Counter
import json

import pytest

from virtualtissue_xenium.v2.cli import compare, main
from virtualtissue_xenium.v2.controls import CONTROLS, apply_control
from virtualtissue_xenium.v2.inputs import prepare, validate_prepared
from virtualtissue_xenium.v2.simulation import run


def context():
    return dict(schema_version="2.0", tissue="skin", organism="human", disease="melanoma", assay="xenium",
                sampling=dict(design="subsampled", description="Test fixture", missing_populations=[]),
                perturbations=[dict(type="treatment", name="drug", targets=[], mechanism_assumptions=[])],
                recognition=True, time_step_hours=1,
                context_modifiers=[dict(action="ACTIVATE", lineages=["T_CD8"], multiplier=2.0,
                                        reason="test", when=dict(perturbation="drug"))])


def cells():
    out = []
    for i in range(12):
        lineage, label = [("malignant", "mal"), ("T_CD8", "cd8"), ("myeloid", "mac")][i % 3]
        out.append(dict(cell_id=f"c{i:02d}", cell_type=label, lineage=lineage, sample_id="s",
                        x=float(i * 7), y=float(i % 2), expression={"MKI67": float(i), "GZMB": float(12 - i), "ACTB": 3.0}))
    return out


@pytest.fixture
def bundle():
    return prepare(cells(), context(), assayed_genes=["MKI67", "GZMB", "ACTB"])


def test_every_control_yields_a_valid_bundle_with_a_control_record(bundle):
    for kind in CONTROLS:
        out = apply_control(bundle, kind, seed=3)
        validate_prepared(out)
        record = out["guidebook"]["control"]
        assert record["kind"] == kind and record["source_fingerprint"] == bundle["fingerprint"]
        assert record["source_input_sha256"] == bundle["guidebook"]["provenance"]["input_sha256"]
        assert out["guidebook"]["provenance"]["gene_sets_sha256"] == bundle["guidebook"]["provenance"]["gene_sets_sha256"]
    with pytest.raises(ValueError):
        apply_control(bundle, "beta_zero")


def test_shuffles_change_only_the_declared_field(bundle):
    by_id = {c["cell_id"]: c for c in bundle["cells"]}
    within = apply_control(bundle, "shuffle_expression_within_type", seed=1)
    for c in within["cells"]:
        src = by_id[c["cell_id"]]
        assert (c["cell_type"], c["lineage"], c["x"], c["y"]) == (src["cell_type"], src["lineage"], src["x"], src["y"])
    # Within-type: the multiset of expression vectors per type is preserved and at least one cell moved.
    for label in ("mal", "cd8", "mac"):
        before = Counter(json.dumps(c["expression"], sort_keys=True) for c in bundle["cells"] if c["cell_type"] == label)
        after = Counter(json.dumps(c["expression"], sort_keys=True) for c in within["cells"] if c["cell_type"] == label)
        assert before == after
    assert any(c["expression"] != by_id[c["cell_id"]]["expression"] for c in within["cells"])
    labels = apply_control(bundle, "label_permutation", seed=1)
    assert Counter(c["lineage"] for c in labels["cells"]) == Counter(c["lineage"] for c in bundle["cells"])
    assert all(c["expression"] == by_id[c["cell_id"]]["expression"] for c in labels["cells"])
    coords = apply_control(bundle, "coordinate_shuffle", seed=1)
    assert sorted((c["x"], c["y"]) for c in coords["cells"]) == sorted((c["x"], c["y"]) for c in bundle["cells"])
    assert any((c["x"], c["y"]) != (by_id[c["cell_id"]]["x"], by_id[c["cell_id"]]["y"]) for c in coords["cells"])
    assert apply_control(bundle, "no_modifier_treatment")["context"]["context_modifiers"] == []
    assert apply_control(bundle, "well_mixed_equivalent")["context"]["interaction_mode"] == "well_mixed"
    assert apply_control(bundle, "shuffle_expression_global", seed=2) == apply_control(bundle, "shuffle_expression_global", seed=2)
    assert apply_control(bundle, "shuffle_expression_global", seed=2) != apply_control(bundle, "shuffle_expression_global", seed=3)


def test_compare_accepts_a_declared_control_arm_and_rejects_an_undeclared_input_change(bundle, tmp_path):
    run(bundle, tmp_path / "real", steps=2, seed=5)
    run(apply_control(bundle, "shuffle_expression_within_type", seed=1), tmp_path / "shuffled", steps=2, seed=5)
    run(apply_control(bundle, "well_mixed_equivalent"), tmp_path / "mixed", steps=2, seed=5)
    report = compare([str(tmp_path / "real"), str(tmp_path / "shuffled"), str(tmp_path / "mixed")])
    assert [r["control"] and r["control"]["kind"] for r in report["runs"]] == [None, "shuffle_expression_within_type", "well_mixed_equivalent"]
    other = cells()
    other[0]["expression"]["ACTB"] = 99.0
    run(prepare(other, context(), assayed_genes=["MKI67", "GZMB", "ACTB"]), tmp_path / "other", steps=2, seed=5)
    with pytest.raises(ValueError, match="no declared control"):
        compare([str(tmp_path / "real"), str(tmp_path / "other")])


def test_cli_control_writes_once_and_lists_kinds(bundle, tmp_path, capsys):
    source = tmp_path / "b.json"
    source.write_text(json.dumps(bundle))
    assert main(["control", str(source), "--kind", "list"]) == 0
    assert "coordinate_shuffle" in capsys.readouterr().out
    out = tmp_path / "ctl.json"
    assert main(["control", str(source), str(out), "--kind", "label_permutation", "--seed", "4"]) == 0
    assert json.loads(out.read_text())["guidebook"]["control"]["seed"] == 4
    assert main(["control", str(source), str(out), "--kind", "label_permutation"]) == 1
    assert main(["control", str(source), "--kind", "label_permutation"]) == 1
