"""Integration regressions found by independent review of the v3 changes."""
import json
from collections import Counter

import pytest

from virtualtissue_xenium.v2.context import validate_context
from virtualtissue_xenium.v2.inputs import load_cells, prepare
from virtualtissue_xenium.v2.provider import prompt_version
from virtualtissue_xenium.v2.simulation import run
from virtualtissue_xenium.v2.cli import compare
from virtualtissue_xenium.v2.controls import apply_control
from test_v3_contracts import cells, context


class UniformHiddenProvider:
    name = "mock_jev"
    model = "jev-1.14.0"
    prior_mode = "hidden"
    prompt_version = prompt_version("hidden", "3.0")
    max_questions = 20
    sampling_parameters = {"temperature": None, "provider_seed": None, "status": "not exposed"}

    def choose(self, rows, context):
        return {row["cell_id"]: dict(probabilities={a: 1 / len(row["actions"]) for a in row["actions"]},
                                    choice=next(iter(row["actions"])), confidence=.5) for row in rows}

    def usage(self):
        return dict(http_attempts=0, input_tokens=0, output_tokens=0, attempts_with_unknown_usage=0)


def test_hidden_prior_direct_run_is_not_pooled_with_half_prior(tmp_path):
    bundle = prepare(cells(1), context())
    result = run(bundle, tmp_path / "hidden", steps=1, mode="jev", provider=UniformHiddenProvider())
    assert result["completed"]
    config = json.loads((tmp_path / "hidden/config.json").read_text())
    assert config["jev_weight"] == 1.0
    decisions = [json.loads(line) for line in (tmp_path / "hidden/decisions.jsonl").read_text().splitlines()]
    assert decisions
    for decision in decisions:
        assert decision["applied_probabilities"] == decision["jev"]["probabilities"]
        assert decision["applied_probabilities"] != decision["prior"]


@pytest.mark.parametrize("platform", ["visium_hd", "stereo_seq"])
def test_ambiguous_spatial_platform_does_not_default_to_segmented_cells(platform):
    with pytest.raises(ValueError, match="unit|lattice|measurement"):
        prepare(cells(1), context(assay=platform))


@pytest.mark.parametrize("key,value", [("kinetics", None), ("kinetics", []), ("scoring", None), ("scoring", [])])
def test_malformed_v3_nested_settings_raise_user_facing_value_error(key, value):
    with pytest.raises(ValueError):
        validate_context(context(**{key: value}))


def test_h5ad_observation_bound_precedes_eager_layer_loading(tmp_path, monkeypatch):
    import anndata
    import h5py
    path = tmp_path / "large.h5ad"
    with h5py.File(path, "w") as handle:
        obs = handle.create_group("obs")
        obs.attrs["_index"] = "named_index"
        obs.create_dataset("named_index", data=[b"a", b"b", b"c"])

    def must_not_open(*args, **kwargs):
        raise AssertionError("read_h5ad called before enforcing max_cells; backed layers may materialize")

    monkeypatch.setattr(anndata, "read_h5ad", must_not_open)
    with pytest.raises(ValueError, match="exceeds max_cells"):
        load_cells(path, max_cells=1, layer="counts")


def test_compare_rejects_different_measured_unit_semantics(tmp_path):
    population = cells(1)
    first = prepare(population, context(), measurement={"unit": "segmented_cell"})
    second = prepare(population, context(), measurement={"unit": "nucleus"})
    run(first, tmp_path / "cells", steps=1, seed=3)
    run(second, tmp_path / "nuclei", steps=1, seed=3)
    with pytest.raises(ValueError, match="measurement|unit"):
        compare([tmp_path / "cells", tmp_path / "nuclei"])


def test_lattice_label_permutation_does_not_destroy_labels_via_purity_gate():
    population = [dict(cell_id=f"c{i}", sample_id="s", cell_type="tumour", lineage="malignant",
                       purity=p, x=i * 8., y=0., expression={"MKI67": 40.})
                  for i, p in enumerate([.9, .4])]
    bundle = prepare(population, context(assay="visium_hd", kernel="lattice"),
        measurement=dict(unit="binned_spot", unit_size_um=8., deconvolution_method="synthetic method"))
    expected = Counter(c["lineage"] for c in bundle["cells"])
    assert expected == {"malignant": 1, "unknown": 1}
    for seed in range(20):
        control = apply_control(bundle, "label_permutation", seed=seed)
        assert Counter(c["lineage"] for c in control["cells"]) == expected
