"""Regression coverage for rounded choices, declared tissue and t0 diagnostics."""
import copy
import importlib.util
import json
import os
from pathlib import Path
from unittest.mock import patch

import pytest

import virtualtissue_xenium as vtx
from virtualtissue_xenium import results
from jev_payload import payload, v2_payload, build_batches, declared_tissue
from jev_provider import JevProvider, ProviderError, validate_response
from recording_analysis import tumour_cycle_diagnostics
from xenium_runtime.engine import validate_answer
from test_provider import rows, response, MODEL, PARAMS
from test_release_012 import fixture

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ablation", ROOT / "tools/jev_prompt_ablation.py")
ablation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ablation)


@pytest.mark.parametrize("p,choice", [
    ({"A": .50, "B": .49}, "B"),
    ({"A": .51, "B": .50}, "B"),
    ({"A": .34, "B": .33, "C": .33}, "B"),
    ({"A": .34, "B": .33, "C": .32}, "B"),
    ({"A": .35, "B": .34, "C": .32}, "B"),
    ({"A": .4950004975, "B": .4949995025}, "B"),
])
def test_raw_values_survive_rounding_repair_and_kernel_accepts(p, choice):
    body = dict(model=MODEL, questions={"c1": dict(criteria={k: "" for k in p})})
    raw = dict(model=MODEL, answers={"c1": dict(type="choice", choice=choice, probabilities=p, confidence=.6)})
    original = copy.deepcopy(raw)
    answers, repairs = validate_response(raw, body)
    answer = answers["c1"]
    validate_answer(answer, body["questions"]["c1"]["criteria"])
    assert answer["reported_probabilities"] == p
    assert answer["choice"] == choice
    assert answer["probabilities"][choice] == max(answer["probabilities"].values())
    assert answer["probability_repairs"] == repairs["c1"]
    assert raw == original


def test_repaired_round_replays_without_another_call(tmp_path):
    calls = []
    def send(body):
        calls.append(body)
        raw = response(body)
        raw["answers"]["c0"].update(choice="ENTER_CYCLE", probabilities=dict(WAIT=.50, ENTER_CYCLE=.49))
        return raw
    def provider():
        return JevProvider("FAKE_SECRET", MODEL, PARAMS, 1, tmp_path / "journal", {}, transport=send, pause_s=0)
    expected = provider().choose(rows(1))
    validate_answer(expected[0], rows(1)[0]["options"])
    assert provider().choose(rows(1)) == expected
    assert len(calls) == 1


def test_rejected_answers_are_auditable_and_redacted(tmp_path):
    def invalid(body):
        raw = response(body)
        raw.update(secret="FAKE_SECRET")
        raw["answers"]["c0"].update(choice="FAKE_SECRET")
        return raw
    journal = tmp_path / "journal"
    p = JevProvider("FAKE_SECRET", MODEL, PARAMS, 1, journal, {}, transport=invalid, pause_s=0)
    with pytest.raises(ProviderError):
        p.choose(rows(1))
    entry = json.loads(journal.read_text().splitlines()[-1])
    assert entry["status"] == "failed"
    assert entry["usage"]["input_tokens"] == 100
    assert set(entry["response"]) == {"model", "answers", "usage"}
    assert "FAKE_SECRET" not in journal.read_text()


@pytest.mark.parametrize("key", ["", "FAKE_SECRET"])
def test_nonfinite_response_still_leaves_a_failed_journal(tmp_path, key):
    def invalid(body):
        raw = response(body)
        raw["answers"]["c0"]["probabilities"]["WAIT"] = float("nan")
        return raw
    journal = tmp_path / "journal"
    p = JevProvider(key, MODEL, PARAMS, 1, journal, {}, transport=invalid, pause_s=0)
    with pytest.raises(ProviderError):
        p.choose(rows(1))
    assert json.loads(journal.read_text().splitlines()[-1])["status"] == "failed"


def test_prompt_snapshots_and_isolation():
    requests = rows(3)
    requests[0]["observation"].update(recognition_sim=True, excluded={"x": ["recognition_declared:true"]})
    requests[1]["observation"]["primary_label"] = "T_CD8_cytotoxic"
    for name, body in {"v2": v2_payload(requests, MODEL),
                       "v3": payload(requests, MODEL, tissue="human breast carcinoma")}.items():
        assert body == json.loads((ROOT / "tests/fixtures" / f"prompt_{name}.json").read_text())
    other = copy.deepcopy(requests)
    other[0]["observation"].update(recognition_sim=False, excluded={"x": ["recognition_declared:false"]})
    assert payload(requests, MODEL, tissue="breast") == payload(other, MODEL, tissue="breast")
    batches = build_batches(requests, MODEL, PARAMS, tissue="breast")
    assert [[r["cell_id"] for r in b] for b in batches] == [[0, 2], [1]]
    assert "LUAD" not in payload(requests, MODEL)["state"]["setting"]
    with pytest.raises(ValueError):
        declared_tissue("x" * 301)


def test_tissue_changes_journal_and_comparison_identity(tmp_path):
    journal = tmp_path / "journal"
    JevProvider("", MODEL, PARAMS, 1, journal, {}, tissue="breast", transport=response, pause_s=0).choose(rows(1))
    with pytest.raises(ProviderError):
        JevProvider("", MODEL, PARAMS, 1, journal, {}, tissue="lung", transport=response, pause_s=0)
    a, _ = fixture()
    b = copy.deepcopy(a)
    a["config"]["prompt_context"] = dict(tissue="breast")
    b["config"]["prompt_context"] = dict(tissue="lung")
    with patch.object(results, "load_results", side_effect=[a, b]), pytest.raises(ValueError, match="not matched"):
        results.compare_runs(["a", "b"])


def test_context_is_opt_in_and_custom_text_does_not_fall_back_to_metadata(tmp_path):
    bundle = dict(folder=tmp_path)
    (tmp_path / "SOURCE.json").write_text(json.dumps(dict(tissue="human lung adenocarcinoma")))
    wf = vtx.workflow
    default = wf.prompt_context(bundle, dict(model=MODEL))
    assert default["mode"] == "none" and default["tissue"] is None
    assert "lung" not in default["state"]["setting"]
    custom = wf.prompt_context(bundle, dict(model=MODEL, tissue="human breast carcinoma"))
    assert custom["mode"] == "custom" and custom["tissue_source"] == "settings"
    assert "human breast carcinoma" in custom["state"]["setting"]
    source = wf.prompt_context(bundle, dict(model=MODEL, context_mode="bundle"))
    assert source["tissue"] == "human lung adenocarcinoma" and source["tissue_source"] == "SOURCE.json"
    (tmp_path / "SOURCE.json").unlink()
    assert wf.prompt_context(bundle, dict(model=MODEL)) == default
    assert wf.prompt_context(bundle, dict(model=MODEL, tissue="human breast carcinoma")) == custom


@pytest.mark.parametrize("settings", [
    dict(context_mode="custom"), dict(context_mode="custom", tissue=" "),
    dict(context_mode="none", tissue="breast"), dict(context_mode="bundle", tissue="breast"),
    dict(context_mode="bundle"), dict(context_mode="invalid"),
])
def test_ambiguous_or_empty_context_choices_fail_before_requests(tmp_path, settings):
    (tmp_path / "SOURCE.json").write_text("{}")
    with pytest.raises(ValueError):
        vtx.workflow.prompt_context(dict(folder=tmp_path), dict(model=MODEL, **settings))


def test_cli_exposes_context_choice_and_keeps_tissue_alias(tmp_path):
    from virtualtissue_xenium import cli
    parser = cli.build_parser()
    bundle = dict(folder=tmp_path)
    for option in ("--context", "--tissue"):
        args = parser.parse_args(["check", "bundle.zip", option, "human breast carcinoma"])
        context = vtx.workflow.prompt_context(bundle, cli._settings(args))
        assert context["mode"] == "custom" and context["tissue"] == "human breast carcinoma"
    args = parser.parse_args(["check", "bundle.zip"])
    assert vtx.workflow.prompt_context(bundle, cli._settings(args))["mode"] == "none"


def test_cycle_diagnostics_use_reported_values_and_binary_denominator():
    def decision(value, choice, time=0, other=None):
        p = dict(REMAIN_QUIESCENT=1-value, ENTER_CYCLE=value)
        if other:
            p.update(other)
        return dict(time_min=time, cell_id=1, status="validated", observation=dict(primary_label="Tumor",
                    t0_aucell_percentile=dict(MITOTIC_CELL_CYCLE=.9)),
                    answer=dict(choice=choice, probabilities=p))
    a = decision(.5, "ENTER_CYCLE")
    a["answer"].update(reported_probabilities=dict(REMAIN_QUIESCENT=.50, ENTER_CYCLE=.49),
                       probability_repairs=["choice_tied_with_maximum_within_rounding"])
    episode = dict(summary=dict(time_min=60), decisions=[a, decision(.03, "REMAIN_QUIESCENT"),
                   decision(.2, "REMAIN_QUIESCENT", other=dict(X=.1)), decision(.9, "ENTER_CYCLE", time=80)])
    report = tumour_cycle_diagnostics(episode)
    assert report["decisions"] == 3
    assert report["probability_mean"] == pytest.approx((.49+.03+.2)/3)
    assert report["binary_fraction_within_0_1_of_half"] == .5
    assert report["repaired_decisions"] == 1
    assert report["by_mitotic_quartile"]["Q4"]["enter_cycle_choices"] == 1
    assert sum(report["histogram_0_1_steps"]) == 3
    assert tumour_cycle_diagnostics(dict(summary=dict(time_min=0), decisions=[]))["probability_mean"] is None
    saved = dict(episode=dict(episode, bodies=[dict(id=1, label="Tumor")]))
    assert results.decision_table(saved).iloc[0].probability == .49
    history = results.cell_history(saved, 1)
    assert history.iloc[0].probabilities["ENTER_CYCLE"] == .49
    assert history.iloc[0].probability_repairs == a["answer"]["probability_repairs"]


def test_ablation_limits_isolation_redaction_and_no_resend(tmp_path):
    makers = ablation.variants("breast")
    requests = rows(5)
    groups = ablation.groups(requests, MODEL, list(makers.values()), 100000, 2)
    assert list(map(len, groups)) == [2, 2, 1]
    with pytest.raises(ValueError):
        ablation.groups(requests, MODEL, list(makers.values()), 1, 2)
    v2 = makers["v2"](requests, MODEL)
    changed = makers["v2_legacy_semantics"](requests, MODEL)
    changed["state"]["semantics"] = v2["state"]["semantics"]
    assert changed == v2  # a single-factor ablation, including prompt_version
    bodies = {k: [make(g, MODEL) for g in groups[:2]] for k, make in makers.items()}
    calls = []
    def send(body):
        calls.append(body)
        raw = response(body)
        raw["secret"] = "FAKE_SECRET"
        return raw
    with pytest.raises(ValueError):
        ablation.send_plan(bodies, 2, tmp_path, "FAKE_SECRET", 23, 100000, send)
    assert not calls and not (tmp_path / "ablation_journal.jsonl").exists()
    _, sent, failed = ablation.send_plan(bodies, 2, tmp_path, "FAKE_SECRET", 24, 100000, send)
    assert sent == len(calls) == 24 and not failed
    assert "FAKE_SECRET" not in (tmp_path / "ablation_journal.jsonl").read_text()
    with pytest.raises(FileExistsError):
        ablation.send_plan(bodies, 2, tmp_path, "FAKE_SECRET", 24, 100000, send)
    assert len(calls) == 24


def test_ablation_interrupted_intent_is_not_resent(tmp_path):
    bodies = {"v3": [payload(rows(1), MODEL)]}
    def interrupted(body):
        raise KeyboardInterrupt()
    with pytest.raises(KeyboardInterrupt):
        ablation.send_plan(bodies, 1, tmp_path, "", 1, 100000, interrupted)
    assert json.loads((tmp_path / "ablation_journal.jsonl").read_text())["phase"] == "intent"
    with pytest.raises(FileExistsError):
        ablation.send_plan(bodies, 1, tmp_path, "", 1, 100000, response)


@pytest.mark.skipif(not os.environ.get("VTX_TEST_BUNDLE"), reason="set VTX_TEST_BUNDLE")
def test_real_bundle_preflight_does_not_confuse_scheduler_with_decision_cap(tmp_path):
    wf = vtx.workflow
    bundle = wf.extract_bundle(os.environ["VTX_TEST_BUNDLE"], tmp_path)
    settings = dict(scenario="declared_recognition", seed=21, horizon_hours=1,
                    max_population=1000, max_decisions=1, model=MODEL)
    facts, body = wf.preflight(bundle, settings)
    activity = facts["t0_activity"]
    assert sum(v["asked_at_t0"] for v in activity["by_type"].values()) == 1
    assert sum(v["eligible_at_t0"] for v in activity["by_type"].values()) == 280
    assert sum(v["deferred_by_decision_cap"] for v in activity["by_type"].values()) == 279
    assert activity["cd8"]["touching_tumour"] == 2
    assert activity["cd8"]["without_question_at_t0"] == 0
    assert facts["prompt_context"]["mode"] == "none"
    assert facts["prompt_context"]["tissue"] is None
    assert body["state"] == facts["prompt_context"]["state"]
    source_facts, _ = wf.preflight(bundle, dict(settings, context_mode="bundle"))
    assert source_facts["prompt_context"]["tissue_source"] == "SOURCE.json"
    settings["tissue"] = "human breast carcinoma"
    sim = wf.create_sim(bundle, settings)
    assert sim.episode()["prompt_context"]["tissue"] == "human breast carcinoma"
    facts, body = wf.preflight(bundle, settings)
    assert "human breast carcinoma" in body["state"]["setting"]
    assert facts["prompt_context"]["tissue_source"] == "settings"
