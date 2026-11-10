"""Frozen-state probes attribute a provider's answers; these tests use an in-process provider only."""
import json

import pytest

from virtualtissue_xenium.v2.cli import main
from virtualtissue_xenium.v2.inputs import digest, read_json, write_json
from virtualtissue_xenium.v5 import probe
from v5_helpers import bundle, context


def responsive(body):
    """Answers that follow the evidence shown (when any) and are otherwise identical across repeats."""
    answers = {}
    for qid, question in body["questions"].items():
        observation = question["instructions"]["observation"]
        actions = sorted(question["criteria"])
        z = observation["t0_evidence"].get("S_CYCLE") if isinstance(observation.get("t0_evidence"), dict) else None
        lift = 0.0 if z is None else max(-.15, min(.15, .1 * z))
        weights = {a: 1.0 for a in actions}
        weights["REST"] = 6.0
        if "PROLIFERATE" in weights:
            weights["PROLIFERATE"] = 1.0 + 10 * (lift + .15)
        total = sum(weights.values())
        probabilities = {a: w / total for a, w in weights.items()}
        answers[qid] = dict(type="choice", choice=max(actions, key=lambda a: (probabilities[a], a)),
                            probabilities=probabilities, confidence=.5)
    return dict(model=body["model"], answers=answers, usage=dict(input_tokens=5, output_tokens=3))


@pytest.fixture()
def planned(tmp_path):
    prepared = bundle()
    result = probe.plan(prepared, tmp_path / "probe", cells_per_lineage=6, seed=4, repeats=2, max_questions=5)
    return prepared, tmp_path / "probe", result


def test_plan_sends_nothing_and_fixes_cells_menus_and_batches(planned):
    prepared, directory, result = planned
    assert not (directory / "journals").exists()
    assert len(result["cells"]) == 12 and {c["lineage"] for c in result["cells"]} == {"malignant", "T_CD8"}
    assert result["estimate"]["http_requests"] == len(result["requests"]) * 2 == 4 * 2 * len(result["groups"]) * 2
    payloads = read_json(directory / "payloads.json")
    baseline = payloads["baseline"]["hidden"]
    for arm in ("evidence_masked", "evidence_exchanged", "neighbourhood_exchanged"):
        for mine, theirs in zip(payloads[arm]["hidden"], baseline):
            assert list(mine["questions"]) == list(theirs["questions"])
            for qid in mine["questions"]:
                assert mine["questions"][qid]["criteria"] == theirs["questions"][qid]["criteria"]  # same menu
                a, b = mine["questions"][qid]["instructions"]["observation"], theirs["questions"][qid]["instructions"]["observation"]
                assert a["cell_id"] == b["cell_id"] and a["state"] == b["state"]
    masked = payloads["evidence_masked"]["hidden"][0]
    assert all(q["instructions"]["observation"]["t0_evidence"] == {} for q in masked["questions"].values())
    assert "masked_by_probe" in json.dumps(masked["state"]["evidence_unavailable"])
    shown = payloads["evidence_masked"]["shown"][0]["questions"]["q0"]["instructions"]["guidebook_prior"]
    assert shown == payloads["baseline"]["shown"][0]["questions"]["q0"]["instructions"]["guidebook_prior"]
    exchanged = result["exchanges"]["evidence_exchanged"]
    assert all(donor is not None and donor != cell for cell, donor in exchanged.items())
    with pytest.raises(ValueError, match="not empty"):
        probe.plan(prepared, directory)


def test_run_refuses_a_budget_below_the_plan_and_is_resumable(planned):
    prepared, directory, result = planned
    needed = result["estimate"]["http_requests"]
    with pytest.raises(ValueError, match="Nothing sent"):
        probe.run_probe(directory, "", max_http_attempts=needed - 1, bundle=prepared, transport=responsive)
    assert not (directory / "journals").exists()
    sent = []
    first = probe.run_probe(directory, "", max_http_attempts=needed, bundle=prepared,
                            transport=lambda body: sent.append(1) or responsive(body))
    assert first["usage"]["http_attempts"] == needed == len(sent)
    again = probe.run_probe(directory, "", max_http_attempts=needed, bundle=prepared,
                            transport=lambda body: pytest.fail("A completed probe request was resent"))
    assert again["usage"]["cache_hits"] == needed


def test_analysis_separates_arm_effects_from_replicate_noise(planned):
    prepared, directory, result = planned
    probe.run_probe(directory, "", max_http_attempts=result["estimate"]["http_requests"], bundle=prepared, transport=responsive)
    analysis = probe.analyze(directory, bundle=prepared, bootstrap_replicates=200)
    hidden = analysis["by_prior"]["hidden"]
    assert hidden["replicate_noise"]["mean"] == 0 and hidden["replicate_noise"]["modal_action_agreement"] == 1
    masked, swapped = hidden["arms"]["evidence_masked"], hidden["arms"]["evidence_exchanged"]
    assert masked["js_excess_over_replicates"]["mean"] > 0 and swapped["js_excess_over_replicates"]["mean"] > 0
    low, high = masked["js_excess_over_replicates"]["bootstrap_ci95_over_request_groups"]
    assert 0 <= low <= masked["js_excess_over_replicates"]["mean"] <= high
    # Both agreement figures are shares of repeat pairs, so they are on one scale.
    assert hidden["arms"]["neighbourhood_exchanged"]["modal_action_agreement_with_baseline"] == 1
    assert masked["baseline_event_mass_same_cells"] == pytest.approx(hidden["arms"]["evidence_exchanged"]["baseline_event_mass_same_cells"])
    assert analysis["request_groups_excluded"] == [] and analysis["n_cells"] == analysis["n_cells_planned"] == 12
    # The mock ignores the neighbourhood, so that arm must show no effect at all.
    assert hidden["arms"]["neighbourhood_exchanged"]["js_between_arm_and_baseline"]["mean"] == 0
    tracking = {(r["lineage"], r["action"]): r for r in hidden["evidence_tracking"]["rows"]}
    row = tracking[("malignant", "PROLIFERATE")]
    assert row["evidence"] == "S_CYCLE" and row["spearman_jev_share_vs_z"] > .9 and row["spearman_guidebook_share_vs_z"] > .9
    assert analysis["fingerprint"] == digest({k: v for k, v in analysis.items() if k != "fingerprint"})
    with pytest.raises(ValueError, match="exists"):
        probe.analyze(directory, bundle=prepared)


def test_edited_plan_or_other_bundle_is_refused(planned):
    prepared, directory, result = planned
    with pytest.raises(ValueError, match="differs from the one the probe was planned on"):
        probe.run_probe(directory, "", max_http_attempts=10 ** 6, bundle=bundle(n=13), transport=responsive)
    edited = read_json(directory / "plan.json")
    edited["settings"]["repeats"] = 1
    write_json(directory / "plan.json", edited)
    with pytest.raises(ValueError, match="unmodified probe plan"):
        probe.analyze(directory, bundle=prepared)


def test_context_exchange_needs_a_different_donor(tmp_path):
    prepared = bundle()
    with pytest.raises(ValueError, match="requires a donor context"):
        probe.plan(prepared, tmp_path / "a", arms=["baseline", "context_exchanged"])
    with pytest.raises(ValueError, match="change nothing"):
        probe.plan(prepared, tmp_path / "b", arms=["baseline", "context_exchanged"], donor_context=context())
    result = probe.plan(prepared, tmp_path / "c", arms=["baseline", "context_exchanged"], prior_modes=["hidden"],
                        donor_context=context(tissue="skin", disease="melanoma"), repeats=1)
    bodies = read_json(tmp_path / "c/payloads.json")
    assert bodies["context_exchanged"]["hidden"][0]["state"]["context"]["tissue"] == "skin"
    assert bodies["baseline"]["hidden"][0]["state"]["context"]["tissue"] == "lung"
    assert result["settings"]["donor_context_sha256"]


def test_cli_plan_reports_cost_and_run_requires_key_and_budget(tmp_path, capsys, monkeypatch):
    prepared = bundle()
    write_json(tmp_path / "prepared.json", prepared)
    assert main(["jev-probe", "plan", str(tmp_path / "p"), "--prepared", str(tmp_path / "prepared.json"),
                 "--cells-per-lineage", "3", "--repeats", "2"]) == 0
    output = json.loads(capsys.readouterr().out)
    assert output["requests_sent"] == 0 and output["estimate"]["http_requests"] > 0 and output["cells"] == 6
    monkeypatch.delenv("TYPESAFE_API_KEY", raising=False)
    monkeypatch.delenv("JEV_API_KEY", raising=False)
    assert main(["jev-probe", "run", str(tmp_path / "p")]) == 1
    assert "TYPESAFE_API_KEY" in capsys.readouterr().err
    monkeypatch.setenv("JEV_API_KEY", "placeholder-not-a-real-key")
    assert main(["jev-probe", "run", str(tmp_path / "p")]) == 1
    assert "--max-http-attempts" in capsys.readouterr().err
    assert not (tmp_path / "p/journals").exists()


def test_exchange_that_changes_nothing_is_not_counted(tmp_path):
    """On a panel where no evidence is usable, a donor's evidence is the same request as the baseline."""
    from virtualtissue_xenium.v2.inputs import prepare
    from v5_helpers import cells
    panel = ["ACTB", "GAPDH", "MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "CDK1", "CD69", "IL2RA", "FOS"]
    plain = prepare(cells(scores=False), context(scoring=dict(method="depth_conditioned_auc", depth_bins=1,
                                                               min_reference_cells=5, min_transcripts=0, min_overlap=9)),
                    assayed_genes=panel)
    result = probe.plan(plain, tmp_path / "probe", arms=["baseline", "evidence_exchanged"], prior_modes=["hidden"],
                        repeats=2, cells_per_lineage=4)
    assert set(result["exchanges"]["evidence_exchanged"].values()) == {None}
    assert result["cells_unchanged_by_design"] == {"evidence_exchanged": 8}
    probe.run_probe(tmp_path / "probe", "", max_http_attempts=result["estimate"]["http_requests"], bundle=plain, transport=responsive)
    arm = probe.analyze(tmp_path / "probe", bundle=plain, bootstrap_replicates=100)["by_prior"]["hidden"]["arms"]["evidence_exchanged"]
    assert arm["n_cells"] == 0 and arm["n_cells_unchanged_by_design"] == 8
    assert arm["js_between_arm_and_baseline"]["mean"] is None and arm["js_excess_over_replicates"]["mean"] is None


def test_failed_request_is_reported_not_resent_and_its_group_is_excluded(planned):
    prepared, directory, result = planned
    needed, calls = result["estimate"]["http_requests"], []

    def flaky(body):
        calls.append(1)
        if len(calls) == 3:
            raise RuntimeError("network down")
        return responsive(body)

    first = probe.run_probe(directory, "", max_http_attempts=needed, bundle=prepared, transport=flaky)
    assert first["completed"] is False and len(first["failed_requests"]) == 1 and len(calls) == needed
    failed = first["failed_requests"][0]
    again = probe.run_probe(directory, "", max_http_attempts=needed, bundle=prepared,
                            transport=lambda body: pytest.fail("A probe request was resent"))
    assert again["completed"] is False and [f["group"] for f in again["failed_requests"]] == [failed["group"]]
    analysis = probe.analyze(directory, bundle=prepared, bootstrap_replicates=100)
    excluded = analysis["request_groups_excluded"]
    assert [g["group"] for g in excluded] == [failed["group"]] and "failed request" in excluded[0]["reason"]
    assert analysis["n_cells"] == analysis["n_cells_planned"] - excluded[0]["cells"]
    # The excluded group is absent from every arm, so arms are still compared on the same cells.
    for section in analysis["by_prior"].values():
        assert section["replicate_noise"]["n_cells"] == analysis["n_cells"]
        assert all(arm["n_cells"] + arm["n_cells_unchanged_by_design"] == analysis["n_cells"] for arm in section["arms"].values())


def test_analysis_rejects_copied_journals_and_bad_bootstrap_settings(planned):
    import shutil
    prepared, directory, result = planned
    probe.run_probe(directory, "", max_http_attempts=result["estimate"]["http_requests"], bundle=prepared, transport=responsive)
    with pytest.raises(ValueError, match="bootstrap_replicates"):
        probe.analyze(directory, bundle=prepared, bootstrap_replicates=0)
    source = directory / "journals/repeat_0/baseline/hidden"
    target = directory / "journals/repeat_1/baseline/hidden"
    shutil.rmtree(target)
    shutil.copytree(source, target)  # a copied repeat is not an independent repeat
    with pytest.raises(ValueError, match="another probe or repeat"):
        probe.analyze(directory, bundle=prepared, bootstrap_replicates=100)


def test_tracking_follows_the_rule_direction_and_flags_withheld_evidence(tmp_path):
    from v5_helpers import rule
    prepared = bundle(rule("PROLIFERATE", "S_CYCLE", direction="opposes"))
    result = probe.plan(prepared, tmp_path / "probe", arms=["baseline", "evidence_masked"], prior_modes=["hidden"], repeats=2)
    probe.run_probe(tmp_path / "probe", "", max_http_attempts=result["estimate"]["http_requests"], bundle=prepared, transport=responsive)
    rows = probe.analyze(tmp_path / "probe", bundle=prepared, bootstrap_replicates=100)["by_prior"]["hidden"]["evidence_tracking"]["rows"]
    row = next(r for r in rows if (r["lineage"], r["action"]) == ("malignant", "PROLIFERATE"))
    # The mock raises PROLIFERATE with the score; the rule says the score opposes it.
    assert row["direction"] == "opposes" and row["spearman_jev_share_vs_z"] < -.9 and row["spearman_guidebook_share_vs_z"] > .9
    assert row["evidence_shown_to_provider"] is True
