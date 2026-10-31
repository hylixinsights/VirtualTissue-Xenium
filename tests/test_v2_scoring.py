"""Lane B: depth-conditioned within-type signature evidence and the beta effect size."""
import math
import random
from statistics import mean

import pytest

from virtualtissue_xenium.v2.cli import compare
from virtualtissue_xenium.v2.controls import apply_control
from virtualtissue_xenium.v2.guidebook import DEFAULT_GENE_SETS
from virtualtissue_xenium.v2.inputs import prepare
from virtualtissue_xenium.v2.scoring import depth_bins, recovery_auc
from virtualtissue_xenium.v2.simulation import run

HOUSEKEEPING = [f"HK{i:02d}" for i in range(60)]
PANEL = sorted({g for genes in DEFAULT_GENE_SETS.values() for g in genes} | set(HOUSEKEEPING))
CYCLE = DEFAULT_GENE_SETS["cell_cycle"]


def context(**scoring):
    s = dict(method="depth_conditioned_auc", depth_bins=2, min_reference_cells=30)
    s.update(scoring)
    return dict(schema_version="2.0", tissue="skin", organism="human", disease="melanoma", assay="scrna",
                interaction_mode="well_mixed", sampling=dict(design="subsampled", description="fixture",
                missing_populations=[]), perturbations=[], recognition=False, time_step_hours=0.5,
                kinetics=dict(mode="hazard_per_hour"), scoring=s)


def population(n_per_type=90, seed=7, small_type=0):
    """Three types; depth varies 10-fold; a third of malignant cells express cell-cycle genes."""
    rng = random.Random(seed)
    cells = []
    spec = [("mal", "malignant"), ("cd8", "T_CD8"), ("mac", "myeloid")] + ([("rare", "fibroblast")] if small_type else [])
    for label, lineage in spec:
        for i in range(small_type if label == "rare" else n_per_type):
            depth = rng.choice([1, 3, 10])
            expr = {g: float(rng.randint(1, 4) * depth) for g in rng.sample(HOUSEKEEPING, 20 + 3 * depth)}
            if label == "mal" and i % 3 == 0:
                expr.update({g: float(rng.randint(3, 8) * depth) for g in CYCLE})
            elif rng.random() < .3:
                expr[rng.choice(CYCLE)] = float(depth)
            cells.append(dict(cell_id=f"{label}{i:03d}", cell_type=label, lineage=lineage, sample_id="s", expression=expr))
    return cells


def hazards(bundle, action, label):
    rows = [v for v in bundle["guidebook"]["cells"].values() if v["cell_type"] == label and action in v["actions"]]
    return [(r["actions"][action]["hazard_per_hour"], r["actions"][action]["prior_hazard_per_hour"], r) for r in rows]


def test_recovery_auc_is_deterministic_rank_based_and_scale_invariant():
    panel = {"A", "B", "C", "D", "E"}
    # Ranks: A=1, B=2, C=3; window 3; signature {A, C}: area (3)+(1) = 4; best (3)+(2) = 5.
    assert recovery_auc({"A": 9, "B": 5, "C": 2}, panel, ["A", "C"], 3) == pytest.approx(0.8)
    assert recovery_auc({"A": 90, "B": 50, "C": 20}, panel, ["A", "C"], 3) == pytest.approx(0.8)
    # Ties share the average rank: A and B tied at 1.5.
    assert recovery_auc({"A": 5, "B": 5}, panel, ["A"], 2) == pytest.approx(1.5 / 2)
    assert recovery_auc({"B": 5}, panel, ["A", "C"], 3) == 0.0
    bins = depth_bins([dict(cell_id=f"c{i}", sample_id="s", expression={"X": float(i)}) for i in range(10)], 2)
    assert [bins[f"c{i}"] for i in range(10)] == [0] * 5 + [1] * 5


def test_reference_hierarchy_and_unscored_cells_keep_multiplier_one():
    cells = population(n_per_type=90, small_type=12)
    by_bin = prepare(cells, context(depth_bins=2), assayed_genes=PANEL)["guidebook"]
    statuses = {v["evidence"]["cell_cycle"]["reference"] for v in by_bin["cells"].values() if v["cell_type"] == "mal"}
    assert statuses == {"cell_type_x_depth_bin"}
    fine = prepare(cells, context(depth_bins=10), assayed_genes=PANEL)["guidebook"]
    assert {v["evidence"]["cell_cycle"]["reference"] for v in fine["cells"].values() if v["cell_type"] == "mal"} == {"cell_type_no_depth"}
    rare = [v for v in by_bin["cells"].values() if v["cell_type"] == "rare"]
    assert all(v["evidence"]["cell_cycle"]["status"] == "insufficient_reference" for v in rare)
    assert all(v["evidence"]["cell_cycle"]["z"] is None for v in rare)
    assert all(v["actions"]["PROLIFERATE"]["hazard_per_hour"] == v["actions"]["PROLIFERATE"]["prior_hazard_per_hour"] for v in rare)
    coverage = by_bin["provenance"]["scoring_coverage"]["cell_cycle"]
    assert coverage["insufficient_reference"] == 12 and coverage["available"] == 270
    # Signatures with too little panel overlap stay explicitly unavailable.
    narrow = prepare(cells, context(min_overlap=9), assayed_genes=PANEL)["guidebook"]
    assert {v["evidence"]["cell_cycle"]["status"] for v in narrow["cells"].values()} == {"insufficient_panel_overlap"}


def test_multipliers_redistribute_but_preserve_the_reference_mean_and_follow_expression():
    bundle = prepare(population(), context(), assayed_genes=PANEL)
    rows = hazards(bundle, "PROLIFERATE", "mal")
    assert mean(h / prior for h, prior, _ in rows) == pytest.approx(1.0, abs=1e-9)
    for ref in {r["evidence"]["cell_cycle"]["reference_id"] for _, _, r in rows}:
        group = [h / prior for h, prior, r in rows if r["evidence"]["cell_cycle"]["reference_id"] == ref]
        assert mean(group) == pytest.approx(1.0, abs=1e-9)
    cycling = [h for h, _, r in rows if len(r["evidence"]["cell_cycle"]["detected_genes"]) == len(CYCLE)]
    quiet = [h for h, _, r in rows if not r["evidence"]["cell_cycle"]["detected_genes"]]
    assert cycling and quiet and min(cycling) > max(quiet)
    audit = [a for _, _, r in rows for a in r["weight_adjustments"] if a.get("action") == "PROLIFERATE"]
    assert {a["source"] for a in audit} == {"depth_conditioned_auc"} and all(a["beta"] == 1.0 for a in audit)


def test_beta_zero_and_per_action_beta():
    cells = population()
    zero = prepare(cells, context(beta=0), assayed_genes=PANEL)
    assert all(a.get("hazard_per_hour") == a.get("prior_hazard_per_hour")
               for v in zero["guidebook"]["cells"].values() for k, a in v["actions"].items() if k != "REST")
    only = prepare(cells, context(beta={"PROLIFERATE": 1.5}), assayed_genes=PANEL)
    cells_out = only["guidebook"]["cells"].values()
    assert any(v["actions"]["PROLIFERATE"]["hazard_per_hour"] != v["actions"]["PROLIFERATE"]["prior_hazard_per_hour"]
               for v in cells_out if v["cell_type"] == "mal")
    assert all(v["actions"]["APOPTOSIS"]["hazard_per_hour"] == v["actions"]["APOPTOSIS"]["prior_hazard_per_hour"] for v in cells_out)
    control = apply_control(prepare(cells, context(), assayed_genes=PANEL), "beta_zero")
    assert control["context"]["scoring"]["beta"] == 0.0
    assert control["guidebook"]["cells"] == zero["guidebook"]["cells"]


def test_depth_bins_remove_a_depth_confound_that_a_single_reference_keeps():
    rng = random.Random(3)
    cells = []
    for i in range(120):
        deep = i % 2 == 0
        expr = {g: float(rng.randint(1, 3)) for g in rng.sample(HOUSEKEEPING, 50 if deep else 15)}
        # Deep cells detect more cycle genes for purely technical reasons.
        for g in CYCLE[: (6 if deep else 2)]:
            expr[g] = float(rng.randint(1, 3) * (10 if deep else 1))
        if deep:
            expr = {g: v * 10 for g, v in expr.items()}
        cells.append(dict(cell_id=f"m{i:03d}", cell_type="mal", lineage="malignant", sample_id="s", expression=expr))
    def mean_percentile(n_bins):
        g = prepare(cells, context(depth_bins=n_bins), assayed_genes=PANEL)["guidebook"]["cells"]
        deep = [g[c["cell_id"]]["evidence"]["cell_cycle"]["percentile"] for c in cells if int(c["cell_id"][1:]) % 2 == 0]
        shallow = [g[c["cell_id"]]["evidence"]["cell_cycle"]["percentile"] for c in cells if int(c["cell_id"][1:]) % 2 == 1]
        return mean(deep), mean(shallow)
    deep1, shallow1 = mean_percentile(1)
    deep2, shallow2 = mean_percentile(2)
    assert deep1 - shallow1 > 0.4          # one reference: depth masquerades as cycling
    assert abs(deep2 - shallow2) < 0.05    # depth bins: compared only with similar depth


def test_qc_gate_makes_low_count_cells_rest_only_and_excludes_them_from_references():
    cells = population()
    cells[0]["expression"] = {"HK00": 1.0}
    bundle = prepare(cells, context(min_transcripts=5), assayed_genes=PANEL)
    entry = bundle["guidebook"]["cells"][cells[0]["cell_id"]]
    assert entry["qc"] == "below_min_transcripts" and set(entry["actions"]) == {"REST"}
    assert entry["evidence"]["cell_cycle"]["status"] == "qc_fail_below_min_transcripts"
    assert bundle["guidebook"]["provenance"]["qc_fail_cells"] == 1
    assert all(cells[0]["cell_id"] not in v["evidence"]["cell_cycle"].get("reference_id", "")
               for v in bundle["guidebook"]["cells"].values())
    assert all(v["evidence"]["cell_cycle"]["reference_size"] <= 89 for v in bundle["guidebook"]["cells"].values()
               if v["cell_type"] == "mal" and v["evidence"]["cell_cycle"]["status"] == "available"
               and v["evidence"]["cell_cycle"]["reference"] == "cell_type_no_depth")


def test_qc_fail_cell_only_rests_in_a_run(tmp_path):
    cells = population()
    cells[0]["expression"] = {"HK00": 1.0}
    report = run(prepare(cells, context(min_transcripts=5), assayed_genes=PANEL), tmp_path / "qc", steps=6, seed=1)
    assert report["completed"]
    lines = (tmp_path / "qc" / "decisions.jsonl").read_text().splitlines()
    import json
    mine = [json.loads(l) for l in lines if json.loads(l)["cell_id"] == cells[0]["cell_id"]]
    assert mine and all(d["action"] == "REST" and set(d["applied_probabilities"]) == {"REST"} for d in mine)


def test_invalid_scoring_settings_are_rejected():
    for bad in (dict(beta=4), dict(beta={"REST": 1}), dict(depth_bins=0), dict(auc_top_fraction=2), dict(unknown=1)):
        with pytest.raises(ValueError):
            prepare(population(n_per_type=30), context(**bad), assayed_genes=PANEL)
    legacy = context()
    legacy["scoring"] = dict(method="rank_legacy", beta=1)
    with pytest.raises(ValueError, match="rank_legacy"):
        prepare(population(n_per_type=30), legacy, assayed_genes=PANEL)


def test_shuffled_expression_arm_is_a_real_control_and_compare_accepts_it(tmp_path):
    bundle = prepare(population(), context(), assayed_genes=PANEL)
    shuffled = apply_control(bundle, "shuffle_expression_within_type", seed=1)
    real_h = {cid: v["actions"]["PROLIFERATE"]["hazard_per_hour"] for cid, v in bundle["guidebook"]["cells"].items() if v["cell_type"] == "mal"}
    shuf_h = {cid: v["actions"]["PROLIFERATE"]["hazard_per_hour"] for cid, v in shuffled["guidebook"]["cells"].items() if v["cell_type"] == "mal"}
    assert real_h != shuf_h                                    # assignment to cells changes
    assert mean(real_h.values()) == pytest.approx(mean(shuf_h.values()), rel=1e-9)   # type mean preserved
    run(bundle, tmp_path / "real", steps=4, seed=2)
    run(shuffled, tmp_path / "shuffled", steps=4, seed=2)
    run(apply_control(bundle, "beta_zero"), tmp_path / "beta0", steps=4, seed=2)
    out = compare([str(tmp_path / "real"), str(tmp_path / "shuffled"), str(tmp_path / "beta0")])
    assert [r["control"] and r["control"]["kind"] for r in out["runs"]] == [None, "shuffle_expression_within_type", "beta_zero"]
