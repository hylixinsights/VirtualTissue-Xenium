"""Evidence policies change which evidence modulates a hazard, never the menu, gates or law."""
from copy import deepcopy
import json
import math

import pytest

from virtualtissue_xenium.v2.context import validate_context
from virtualtissue_xenium.v2.controls import apply_control
from virtualtissue_xenium.v2.inputs import digest, prepare, read_json, validate_prepared
from virtualtissue_xenium.v2.simulation import initialize, preflight, run
from virtualtissue_xenium.v3.expression_audit import audit_expression
from virtualtissue_xenium.v3.ontology import lineage_rules
from virtualtissue_xenium.v5.evidence import (POLICY_VERSION, V4_RULES, action_evidence, v4_equivalent_policy,
                                              validate_measured_evidence, validate_policy)
from v5_helpers import PANEL, bundle, cells, context, manifest, policy, rule


def hazards(prepared):
    return {cid: {a: spec["hazard_per_hour"] for a, spec in entry["actions"].items()}
            for cid, entry in prepared["guidebook"]["cells"].items()}


def adjustment(prepared, cell_id, action):
    return next(a for a in prepared["guidebook"]["cells"][cell_id]["weight_adjustments"] if a.get("action") == action)


# ---- no policy: v4 exactly

def test_context_without_policy_gains_no_key_and_builds_the_v4_guidebook():
    plain = prepare(cells(scores=False), context(), assayed_genes=PANEL)
    assert "evidence_policy" not in plain["context"]
    assert plain["guidebook"]["provenance"]["guidebook_version"] == "v3-guidebook-1"
    assert "evidence_policy_sha256" not in plain["guidebook"]["provenance"]
    assert "measured_evidence" not in plain
    validate_prepared(plain)


def test_policy_restating_the_v4_table_reproduces_v4_hazards_and_trajectory(tmp_path):
    plain = prepare(cells(scores=False), context(), assayed_genes=PANEL)
    restated = prepare(cells(scores=False), context(evidence_policy=v4_equivalent_policy()), assayed_genes=PANEL)
    assert [(r["action"], r["evidence"]) for r in restated["context"]["evidence_policy"]["rules"]] == list(V4_RULES)
    assert hazards(restated) == hazards(plain)
    assert any(a.get("multiplier", 1) != 1 for e in plain["guidebook"]["cells"].values() for a in e["weight_adjustments"])
    first = run(plain, tmp_path / "plain", steps=6, seed=2)
    second = run(restated, tmp_path / "restated", steps=6, seed=2)
    assert (tmp_path / "plain/frames.jsonl").read_bytes() == (tmp_path / "restated/frames.jsonl").read_bytes()
    assert first["events"] == second["events"]


# ---- policy validation

@pytest.mark.parametrize("change,match", [
    (dict(policy_version="other"), "policy_version"),
    (dict(rules="none"), "rules must be a list"),
    (dict(rules=[rule("REST", "S_CYCLE")]), "non-REST"),
    (dict(rules=[rule("FLY", "S_CYCLE")]), "not permitted for any lineage"),
    (dict(rules=[rule("KILL", "S_CYCLE", lineages=["malignant"])]), "cannot add an action"),
    (dict(rules=[rule("PROLIFERATE", "S_CYCLE", lineages=["martian"])]), "unsupported lineages"),
    (dict(rules=[rule("PROLIFERATE", "S_CYCLE", kind="vibes")]), "kind"),
    (dict(rules=[rule("PROLIFERATE", "S_CYCLE", direction="sideways")]), "direction"),
    (dict(rules=[rule("PROLIFERATE", "S_CYCLE"), rule("PROLIFERATE", "S_ACT", lineages=["malignant"])]), "one evidence source"),
    (dict(rules=[dict(rule("PROLIFERATE", "S_CYCLE"), rationale="")]), "rationale"),
    (dict(provider_evidence=["S_CYCLE", "S_CYCLE"]), "provider_evidence"),
    (dict(unexpected=1), "Unknown evidence_policy"),
])
def test_invalid_policies_are_rejected(change, match):
    value = dict(policy(rule("PROLIFERATE", "S_CYCLE")), **change)
    with pytest.raises(ValueError, match=match):
        validate_policy(value, lineage_rules())


def test_policy_is_schema3_only_and_normalisation_is_idempotent():
    with pytest.raises(ValueError, match="schema 3.0"):
        validate_context(dict(context(schema_version="2.0", evidence_policy=policy()), scoring={"method": "rank_legacy"},
                              dynamics=None))
    once = validate_context(context(evidence_policy=policy(rule("ACTIVATE", "S_ACT", lineages=["T_CD8", "B"]))))
    assert once["evidence_policy"]["rules"][0]["lineages"] == ["B", "T_CD8"]
    assert once["evidence_policy"]["rules"][0]["review_status"] == "requires_biological_review"
    assert validate_context(once) == once


def test_rules_must_reference_declared_evidence():
    with pytest.raises(ValueError, match="not declared for this bundle"):
        prepare(cells(), context(evidence_policy=policy(rule("PROLIFERATE", "S_MISSING"))), assayed_genes=PANEL,
                measured_evidence=manifest())
    with pytest.raises(ValueError, match="not declared for this bundle"):
        prepare(cells(), context(evidence_policy=policy(rule("PROLIFERATE", "S_CYCLE", kind="gene_set"))),
                assayed_genes=PANEL, measured_evidence=manifest())
    with pytest.raises(ValueError, match="provider_evidence names undeclared"):
        prepare(cells(), context(evidence_policy=policy(rule("PROLIFERATE", "S_CYCLE"), provider_evidence=["nope"])),
                assayed_genes=PANEL, measured_evidence=manifest())


# ---- measured scores

def test_measured_scores_need_a_policy_and_a_manifest():
    with pytest.raises(ValueError, match="evidence_policy"):
        prepare(cells(), context(), assayed_genes=PANEL, measured_evidence=manifest())
    with pytest.raises(ValueError, match="without a measured_evidence manifest"):
        prepare(cells(), context(evidence_policy=policy()), assayed_genes=PANEL)
    incomplete = manifest()
    del incomplete["scores"]["S_ACT"]
    with pytest.raises(ValueError, match="without a manifest entry"):
        prepare(cells(), context(evidence_policy=policy()), assayed_genes=PANEL, measured_evidence=incomplete)
    with pytest.raises(ValueError, match="method"):
        validate_measured_evidence(dict(scores=dict(S=dict(description="x"))), [])
    bad = cells()
    bad[0]["measured_scores"]["S_CYCLE"] = float("nan")
    with pytest.raises(ValueError, match="finite numbers, or null"):
        prepare(bad, context(evidence_policy=policy()), assayed_genes=PANEL, measured_evidence=manifest())


def test_measured_score_redistributes_the_hazard_and_preserves_the_group_mean():
    prepared = bundle()
    validate_prepared(prepared)
    tumour = {cid: h["PROLIFERATE"] for cid, h in hazards(prepared).items() if cid.startswith("tumour")}
    ordered = [tumour[f"tumour{i:02}"] for i in range(12)]
    assert ordered == sorted(ordered) and ordered[0] < ordered[-1]  # higher score, higher hazard
    base = prepared["guidebook"]["provenance"]["default_hazards_per_hour"]["malignant"]["PROLIFERATE"]
    assert sum(ordered) / len(ordered) == pytest.approx(base)  # mean of the reference group is unchanged
    record = adjustment(prepared, "tumour05", "PROLIFERATE")
    assert record["source"] == "measured_score" and record["signature"] == "S_CYCLE"
    assert record["evidence_rule"]["review_status"] == "requires_biological_review"
    evidence = prepared["guidebook"]["cells"]["tumour05"]["evidence"]["S_CYCLE"]
    assert evidence["kind"] == "measured_score" and evidence["measured_value"] == pytest.approx(5 / 12)
    assert evidence["reference"] == "cell_type_x_depth_bin" and evidence["reference_size"] == 12


def test_missing_score_stays_missing_and_neutral():
    prepared = bundle()
    evidence = prepared["guidebook"]["cells"]["tcell11"]["evidence"]["S_ACT"]
    assert evidence["status"] == "not_measured" and evidence["z"] is None and evidence["measured_value"] is None
    assert adjustment(prepared, "tcell11", "ACTIVATE")["multiplier"] == 1.0
    base = prepared["guidebook"]["provenance"]["default_hazards_per_hour"]["T_CD8"]["ACTIVATE"]
    assert hazards(prepared)["tcell11"]["ACTIVATE"] == base
    # The unscored cell is not a reference for the others.
    assert prepared["guidebook"]["cells"]["tcell00"]["evidence"]["S_ACT"]["reference_size"] == 11
    assert prepared["guidebook"]["provenance"]["measured_coverage"]["S_ACT"]["not_measured"] == 1


def test_small_reference_is_neutral_not_invented():
    prepared = prepare(cells(4), context(evidence_policy=policy(rule("PROLIFERATE", "S_CYCLE"))), assayed_genes=PANEL,
                       measured_evidence=manifest())
    statuses = {e["evidence"]["S_CYCLE"]["status"] for e in prepared["guidebook"]["cells"].values()}
    assert statuses == {"insufficient_reference"}
    assert all(a.get("multiplier", 1) == 1 for e in prepared["guidebook"]["cells"].values() for a in e["weight_adjustments"])


def test_opposing_rule_inverts_the_direction_with_the_same_law():
    supports = bundle(rule("PROLIFERATE", "S_CYCLE"))
    opposes = bundle(rule("PROLIFERATE", "S_CYCLE", direction="opposes"))
    up, down = hazards(supports), hazards(opposes)
    assert up["tumour11"]["PROLIFERATE"] > up["tumour00"]["PROLIFERATE"]
    assert down["tumour11"]["PROLIFERATE"] < down["tumour00"]["PROLIFERATE"]
    assert adjustment(opposes, "tumour03", "PROLIFERATE")["beta"] == -adjustment(supports, "tumour03", "PROLIFERATE")["beta"]


def test_policy_never_changes_the_menu_or_the_gates():
    plain = prepare(cells(scores=False), context(), assayed_genes=PANEL)
    with_policy = bundle()
    assert {cid: sorted(e["actions"]) for cid, e in with_policy["guidebook"]["cells"].items()} == \
           {cid: sorted(e["actions"]) for cid, e in plain["guidebook"]["cells"].items()}
    assert all("KILL" not in e["actions"] for e in with_policy["guidebook"]["cells"].values())  # recognition is false
    assert with_policy["guidebook"]["provenance"]["kinetics_sha256"] == plain["guidebook"]["provenance"]["kinetics_sha256"]


def test_lineage_scoped_rule_leaves_other_lineages_label_only():
    prepared = bundle(rule("ACTIVATE", "S_ACT", lineages=["T_CD8"]))
    assert all(a.get("action") != "PROLIFERATE" or "signature" not in a
               for e in prepared["guidebook"]["cells"].values() for a in e["weight_adjustments"])
    summary = action_evidence(prepared)
    assert "PROLIFERATE" in summary["label_only_actions"] and "ACTIVATE" not in summary["label_only_actions"]
    assert summary["by_action"]["PROLIFERATE"]["status_counts"] == {"no_rule_declared": 24}
    assert summary["policy"] == "test_policy"


def test_preflight_warns_when_declared_evidence_is_unusable():
    # The packaged activation signature has too few assayed genes on this panel.
    restated = prepare(cells(scores=False), context(evidence_policy=v4_equivalent_policy()),
                       assayed_genes=[g for g in PANEL if g not in ("IFNG", "TNFRSF9", "JUN")])
    report = preflight(restated)["action_evidence"]
    assert "ACTIVATE" in report["label_only_actions"]
    assert any("ACTIVATE" in w and "insufficient_panel_overlap" in w for w in report["warnings"])
    plain = preflight(prepare(cells(scores=False), context(), assayed_genes=PANEL))
    assert plain["action_evidence"]["policy"] == "v4_default_table"


# ---- bundle integrity, controls, state

def test_prepared_bundle_binds_manifest_and_policy():
    prepared = bundle()
    assert prepared["guidebook"]["provenance"]["evidence_policy_sha256"] == digest(prepared["context"]["evidence_policy"])
    assert prepared["guidebook"]["provenance"]["guidebook_version"] == "v5-guidebook-1"
    edited = deepcopy(prepared)
    edited["measured_evidence"]["scores"]["S_CYCLE"]["method"] = "something else"
    edited["fingerprint"] = digest({k: v for k, v in edited.items() if k != "fingerprint"})
    with pytest.raises(ValueError, match="manifest"):
        validate_prepared(edited)


def test_expression_shuffle_moves_measured_scores_with_the_counts():
    prepared = bundle()
    shuffled = apply_control(prepared, "shuffle_expression_within_type", seed=3)
    validate_prepared(shuffled)
    original = {c["cell_id"]: (digest(c["expression"]), digest(c["measured_scores"])) for c in prepared["cells"]}
    pairs = set(original.values())
    moved = 0
    for cell in shuffled["cells"]:
        pair = (digest(cell["expression"]), digest(cell["measured_scores"]))
        assert pair in pairs  # counts and scores travel as one unit
        moved += pair != original[cell["cell_id"]]
    assert moved > 0
    assert shuffled["guidebook"]["control"]["measured_scores"].startswith("permuted together")


def test_beta_zero_neutralises_measured_evidence():
    prepared = bundle()
    label_only = apply_control(prepared, "beta_zero")
    assert len({h["PROLIFERATE"] for cid, h in hazards(label_only).items() if cid.startswith("tumour")}) == 1
    audit = audit_expression(prepared, seed=1)
    assert audit["comparisons"]["beta_zero"]["summary"]["n_probability_changed"] > 0
    declared = {row["action_probabilities"]["PROLIFERATE"]["declared_signature"]
                for row in audit["comparisons"]["beta_zero"]["cells"] if "PROLIFERATE" in row["action_probabilities"]}
    assert declared == {"S_CYCLE"}


def test_measured_scores_are_not_copied_into_the_evolving_state(tmp_path):
    prepared = bundle()
    assert all("measured_scores" not in state and "expression" not in state for state in initialize(prepared).values())
    report = run(prepared, tmp_path / "run", steps=3, seed=1)
    assert report["completed"]
    assert "measured_scores" not in (tmp_path / "run/frames.jsonl").read_text()


# ---- inputs

def test_json_object_carries_scores_and_manifest_through_the_cli(tmp_path, capsys):
    from virtualtissue_xenium.v2.cli import main
    (tmp_path / "object.json").write_text(json.dumps(dict(cells=cells(), assayed_genes=PANEL, measured_evidence=manifest())))
    (tmp_path / "context.json").write_text(json.dumps(context()))
    (tmp_path / "policy.json").write_text(json.dumps(policy(rule("PROLIFERATE", "S_CYCLE"))))
    assert main(["prepare", str(tmp_path / "object.json"), str(tmp_path / "context.json"), str(tmp_path / "prepared.json"),
                 "--evidence-policy", str(tmp_path / "policy.json")]) == 0
    output = json.loads(capsys.readouterr().out)
    assert output["preflight"]["action_evidence"]["by_action"]["PROLIFERATE"]["cells_with_changed_hazard"] > 0
    prepared = read_json(tmp_path / "prepared.json")
    assert sorted(prepared["measured_evidence"]["scores"]) == ["S_ACT", "S_CYCLE"]
    assert main(["evidence-policy", str(tmp_path / "template.json")]) == 0
    assert read_json(tmp_path / "template.json")["policy_version"] == POLICY_VERSION


def test_anndata_score_columns_keep_nan_as_null():
    anndata = pytest.importorskip("anndata")
    numpy = pytest.importorskip("numpy")
    pandas = pytest.importorskip("pandas")
    from virtualtissue_xenium.v2.inputs import from_anndata
    obs = pandas.DataFrame(dict(cell_type=["t", "t"], lineage=["malignant", "malignant"], score=[0.25, float("nan")]),
                           index=["c1", "c2"])
    adata = anndata.AnnData(X=numpy.array([[1., 0.], [0., 2.]]), obs=obs, var=pandas.DataFrame(index=["A", "B"]))
    loaded, _ = from_anndata(adata, score_keys=["score"])
    assert [c["measured_scores"]["score"] for c in loaded] == [0.25, None]
    with pytest.raises(ValueError, match="Missing measured-score column"):
        from_anndata(adata, score_keys=["absent"])
    assert all(not math.isnan(v) for c in loaded for v in c["measured_scores"].values() if v is not None)


def test_roi_selection_keeps_measured_scores_and_applies_the_policy(tmp_path, capsys):
    from virtualtissue_xenium.v2.cli import main
    from virtualtissue_xenium.v3.roi_catalog import build_roi_catalog
    source = dict(cells=cells(), assayed_genes=PANEL, measured_evidence=manifest())
    (tmp_path / "object.json").write_text(json.dumps(source))
    (tmp_path / "context.json").write_text(json.dumps(context()))
    (tmp_path / "policy.json").write_text(json.dumps(policy(rule("PROLIFERATE", "S_CYCLE"))))
    selection = build_roi_catalog(source, cells_per_roi=20, grid_um=60)["candidates"][0]["selection"]
    (tmp_path / "selection.json").write_text(json.dumps(selection))
    assert main(["select-roi", str(tmp_path / "object.json"), str(tmp_path / "selection.json"), str(tmp_path / "roi.json"),
                 "--context", str(tmp_path / "context.json"), "--evidence-policy", str(tmp_path / "policy.json")]) == 0
    capsys.readouterr()
    roi = read_json(tmp_path / "roi.json")
    validate_prepared(roi)
    assert len(roi["cells"]) == 20 and all("measured_scores" in c for c in roi["cells"])
    assert sorted(roi["measured_evidence"]["scores"]) == ["S_ACT", "S_CYCLE"]
    assert roi["context"]["evidence_policy"]["name"] == "test_policy"
    assert any(a.get("source") == "measured_score" for e in roi["guidebook"]["cells"].values() for a in e["weight_adjustments"])


def test_selected_decisions_carry_more_declared_evidence_than_offered_ones(tmp_path, capsys):
    from virtualtissue_xenium.v2.cli import main
    from virtualtissue_xenium.v5.evidence import decision_evidence_association
    prepared = bundle(n=40)
    run(prepared, tmp_path / "baseline", steps=60, seed=5)
    run(apply_control(prepared, "beta_zero"), tmp_path / "label_only", steps=60, seed=5)
    def difference(folder):
        rows = {(r["lineage"], r["action"]): r for r in decision_evidence_association(tmp_path / folder)["rows"]}
        row = rows[("T_CD8", "ACTIVATE")]
        assert row["evidence"] == "S_ACT" and row["decisions_selected"] > 20
        return row["difference"], row["standard_error_selected"]
    with_evidence, error = difference("baseline")
    label_only, _ = difference("label_only")
    assert with_evidence > 3 * error and with_evidence > label_only + .3
    assert abs(label_only) < 3 * error
    assert main(["evidence-decisions", str(tmp_path / "baseline"), str(tmp_path / "label_only"),
                 "--output", str(tmp_path / "association.json")]) == 0
    capsys.readouterr()
    saved = read_json(tmp_path / "association.json")
    assert [r["control"] for r in saved["runs"]] == [None, "beta_zero"]


def test_unscored_cell_does_not_shrink_the_reference_of_scored_cells():
    """With type-level fallback, a stratum made small only by an unscored cell must still fall back to the type."""
    population = [c for c in cells(10) if c["lineage"] == "malignant"]
    population[9]["measured_scores"]["S_CYCLE"] = None  # the deepest cell was not scored
    prepared = prepare(population, context(evidence_policy=policy(rule("PROLIFERATE", "S_CYCLE")),
                                           scoring=dict(method="depth_conditioned_auc", depth_bins=2, min_reference_cells=5,
                                                        min_transcripts=0, auc_top_fraction=1.0, reference_fallback="cell_type")),
                       assayed_genes=PANEL, measured_evidence=manifest())
    records = {cid: e["evidence"]["S_CYCLE"] for cid, e in prepared["guidebook"]["cells"].items()}
    assert records["tumour09"]["status"] == "not_measured"
    scored = [r for cid, r in records.items() if cid != "tumour09"]
    assert {r["status"] for r in scored} == {"available"}
    assert {r["reference"] for r in scored} == {"cell_type_no_depth"} and {r["reference_size"] for r in scored} == {9}


def test_anndata_missing_and_boolean_score_columns():
    anndata = pytest.importorskip("anndata")
    numpy = pytest.importorskip("numpy")
    pandas = pytest.importorskip("pandas")
    from virtualtissue_xenium.v2.inputs import from_anndata
    obs = pandas.DataFrame(dict(cell_type=["t", "t"], lineage=["malignant", "malignant"],
                                nullable=pandas.array([0.5, pandas.NA], dtype="Float64"), flag=[True, False]), index=["c1", "c2"])
    adata = anndata.AnnData(X=numpy.array([[1., 0.], [0., 2.]]), obs=obs, var=pandas.DataFrame(index=["A", "B"]))
    loaded, _ = from_anndata(adata, score_keys=["nullable"])
    assert [c["measured_scores"]["nullable"] for c in loaded] == [0.5, None]
    with pytest.raises(ValueError, match="not boolean"):
        from_anndata(adata, score_keys=["flag"])


def test_manifest_has_one_source_and_malformed_logs_are_user_errors(tmp_path, capsys):
    from virtualtissue_xenium.v2.cli import main
    from virtualtissue_xenium.v2.simulation import read_decisions
    (tmp_path / "object.json").write_text(json.dumps(dict(cells=cells(), assayed_genes=PANEL, measured_evidence=manifest())))
    (tmp_path / "context.json").write_text(json.dumps(context(evidence_policy=policy(rule("PROLIFERATE", "S_CYCLE")))))
    (tmp_path / "manifest.json").write_text(json.dumps(manifest()))
    assert main(["prepare", str(tmp_path / "object.json"), str(tmp_path / "context.json"), str(tmp_path / "out.json"),
                 "--measured-evidence", str(tmp_path / "manifest.json")]) == 1
    assert "use one source" in capsys.readouterr().err
    (tmp_path / "decisions.jsonl").write_text("[1, 2]\n")
    with pytest.raises(ValueError, match="not a decision record"):
        list(read_decisions(tmp_path))
    reference = dict(file="guidebook.json", guidebook_id="a", sha256="0" * 64)
    (tmp_path / "decisions.jsonl").write_text(json.dumps(dict(observation=dict(t0_evidence_ref=reference))) + "\n")
    for content in ("[]", json.dumps(dict(cells=dict(a="text")))):
        (tmp_path / "guidebook.json").write_text(content)
        with pytest.raises(ValueError, match="not this run's guidebook"):
            list(read_decisions(tmp_path))
