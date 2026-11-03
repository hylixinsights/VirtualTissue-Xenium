"""Rare types retain bounded panel evidence without invented reference cells."""
from copy import deepcopy
from itertools import combinations
from statistics import mean, pvariance

import pytest

from virtualtissue_xenium.v2.scoring import (
    SCORING_DEFAULTS, expression_multiplier, panel_null_moments,
    recovery_auc, score_population,
)

PANEL = set("ABCDEF")
SIGNATURES = {"cell_cycle": ["A", "B"]}


def settings(**changes):
    value = dict(SCORING_DEFAULTS, method="tiered_panel_auc", reference_scope="sample",
                 reference_fallback="cell_type", min_overlap=2, auc_min_rank=2,
                 depth_bins=2, panel_max_multiplier=2)
    value.update(changes)
    return value


def cell(cid="rare", *, sample="s", label="rare", expression=None):
    return dict(cell_id=cid, sample_id=sample, cell_type=label, lineage="malignant",
                expression=expression or dict(A=20, B=10, C=3, D=2, E=1, F=1), x=0, y=0)


def test_panel_null_moments_match_all_gene_subsets_with_ties_and_zeros():
    expression = dict(A=6, B=4, C=4, D=1, E=0, F=0)
    actual = panel_null_moments(expression, PANEL, 2, 4)
    all_aucs = [recovery_auc(expression, PANEL, genes, 4) for genes in combinations(sorted(PANEL), 2)]
    assert actual["mean_auc"] == pytest.approx(mean(all_aucs))
    assert actual["sd_auc"]**2 == pytest.approx(pvariance(all_aucs))
    assert "not biological calibration" in actual["assumption"]


def test_singleton_panel_evidence_is_bounded_and_not_divided_by_itself():
    population, refs = score_population([cell()], PANEL, SIGNATURES, settings())
    evidence = population["rare"]["cell_cycle"]
    assert evidence["status"] == "available_panel_heuristic"
    assert evidence["reference_size"] == 0 and evidence["candidate_cell_reference_size"] == 1
    assert evidence["percentile"] is None and evidence["calibration"] == "weak_panel_calibration"
    assert evidence["auc"] == 1 and evidence["z"] > 0
    multiplier, normalizer = expression_multiplier(evidence, refs["cell_cycle"], 1)
    assert 1 < multiplier <= 2 and normalizer == 1
    assert expression_multiplier(evidence, refs["cell_cycle"], 0) == (1, 1)
    assert expression_multiplier(evidence, refs["cell_cycle"], -3)[0] >= .5
    assert refs["cell_cycle"] == {}  # No fabricated reference-cell observations.


def test_hierarchy_partitions_empirical_groups_and_preserves_their_mean():
    cells = [cell(f"large{i}", sample="large", label="large", expression=dict(A=i+1, B=3, C=20, D=2, E=1, F=1)) for i in range(60)]
    cells += [cell(f"mid{i}", sample="mid", label="mid", expression=dict(A=i+1, B=3, C=20, D=2, E=1, F=1)) for i in range(35)]
    cells += [cell()]
    scores, refs = score_population(cells, PANEL, SIGNATURES, settings())
    large = [scores[f"large{i}"]["cell_cycle"] for i in range(60)]
    mid = [scores[f"mid{i}"]["cell_cycle"] for i in range(35)]
    assert {r["reference"] for r in large} == {"cell_type_x_depth_bin"}
    assert {r["reference_size"] for r in large} == {30}
    assert {r["reference"] for r in mid} == {"cell_type_no_depth"}
    assert {r["reference_size"] for r in mid} == {35}
    for group in (large[:30], large[30:], mid):
        # Depth ordering is by count, matching these monotonically increasing cells.
        assert mean(expression_multiplier(r, refs["cell_cycle"], 1)[0] for r in group) == pytest.approx(1)
    assert scores["rare"]["cell_cycle"]["reference"] == "panel_gene_set_null"


def test_rare_panel_score_does_not_borrow_other_types_or_samples():
    first, _ = score_population([cell()], PANEL, SIGNATURES, settings())
    others = [cell(f"other{i}", sample="other", label="other", expression=dict(A=1, B=1, C=100)) for i in range(100)]
    second, _ = score_population([cell(), *others], PANEL, SIGNATURES, settings())
    a, b = first["rare"]["cell_cycle"], second["rare"]["cell_cycle"]
    assert a == b
    with pytest.raises(ValueError, match="sample-scoped"):
        score_population([cell()], PANEL, SIGNATURES, dict(settings(), reference_scope="pooled"))


@pytest.mark.parametrize("kind", ["overlap", "undetected", "constant_panel", "qc"])
def test_unsupported_evidence_remains_neutral(kind):
    source, genes, scoring = cell(), SIGNATURES, settings()
    if kind == "overlap":
        genes = {"cell_cycle": ["NOT_IN_PANEL"]}
    elif kind == "undetected":
        source["expression"] = dict(C=30, D=20)
    elif kind == "constant_panel":
        source["expression"] = {g: 10 for g in PANEL}
    else:
        scoring["_qc_fail"] = {source["cell_id"]}
    scores, refs = score_population([source], PANEL, genes, scoring)
    record = scores["rare"]["cell_cycle"]
    assert record["z"] is None and record["informative"] is False
    assert expression_multiplier(record, refs["cell_cycle"], 1) == (1, 1)


def test_legacy_depth_mode_is_unchanged_for_rare_cells():
    original = dict(settings(), method="depth_conditioned_auc", reference_fallback="null")
    original.pop("panel_max_multiplier")
    scores, refs = score_population([cell()], PANEL, SIGNATURES, original)
    assert scores["rare"]["cell_cycle"]["status"] == "insufficient_reference"
    assert expression_multiplier(scores["rare"]["cell_cycle"], refs["cell_cycle"], 1) == (1, 1)


def test_expression_audit_attributes_fixed_state_local_changes_without_claiming_jev():
    from virtualtissue_xenium.v2.inputs import prepare
    from virtualtissue_xenium.v3.expression_audit import audit_expression
    context = dict(schema_version="3.0", tissue="skin", organism="human", disease="fixture",
        assay="xenium", recognition=False, perturbations=[], time_step_hours=.1,
        sampling=dict(design="subsampled", description="expression audit fixture", missing_populations=[]),
        scoring=dict(method="tiered_panel_auc", min_overlap=2, auc_min_rank=2, min_transcripts=0, beta=1))
    prepared = prepare([cell()], context, assayed_genes=sorted(PANEL), gene_sets=SIGNATURES)
    unchanged = deepcopy(prepared)
    audit = audit_expression(prepared, seed=4)
    assert prepared == unchanged
    assert audit == audit_expression(prepared, seed=4)
    beta = audit["comparisons"]["beta_zero"]
    assert beta["summary"]["n_probability_changed"] == 1
    assert beta["summary"]["n_evidence_changed"] == 0
    assert beta["summary"]["n_menu_changed"] == 0
    assert audit["comparisons"]["shuffle_expression_within_type"]["summary"]["n_probability_changed"] == 0
    assert audit["jev_response_effect_measured"] is False
    action = beta["cells"][0]["action_probabilities"]["PROLIFERATE"]
    assert action["reference_hazard_per_hour"] > action["counterfactual_hazard_per_hour"]
    assert action["local_path"] == "expression_hazard_or_qc"
    assert beta["cells"][0]["action_probabilities"]["REST"]["local_path"] == "competing_risk_redistribution"


@pytest.mark.parametrize("case", ["unknown_lineage", "unassayed_signature", "unit_bound"])
def test_audit_accepts_legitimate_zero_influence_without_forcing_positive_results(case):
    from virtualtissue_xenium.v2.inputs import prepare
    from virtualtissue_xenium.v3.expression_audit import audit_expression
    source = cell()
    genes = SIGNATURES
    scoring = dict(method="tiered_panel_auc", min_overlap=2, auc_min_rank=2, min_transcripts=0, beta=1)
    if case == "unknown_lineage":
        source["lineage"] = "unknown"
    elif case == "unassayed_signature":
        genes = {"cell_cycle": ["NOT1", "NOT2"]}
    else:
        scoring["panel_max_multiplier"] = 1
    context = dict(schema_version="3.0", tissue="skin", organism="human", disease="fixture",
        assay="xenium", recognition=False, perturbations=[], time_step_hours=.1,
        sampling=dict(design="subsampled", description="neutral audit fixture", missing_populations=[]), scoring=scoring)
    prepared = prepare([source], context, assayed_genes=sorted(PANEL), gene_sets=genes)
    audit = audit_expression(prepared)
    assert audit["comparisons"]["beta_zero"]["summary"]["n_probability_changed"] == 0
    if case == "unknown_lineage":
        assert audit["comparisons"]["beta_zero"]["summary"]["n_baseline_multiple_actions"] == 0
    if case == "unassayed_signature":
        support = audit["arms"]["baseline"]["coverage_by_type"][0]["signatures"]["cell_cycle"]
        assert support["n_without_score"] == 1 and support["n_with_z"] == 0
