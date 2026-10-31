"""Paid-provider boundaries, bounded concurrency and crash-safe cache behavior."""
from copy import deepcopy
import json
import threading

import pytest

from virtualtissue_xenium.v2.inputs import canonical, digest
from virtualtissue_xenium.v2.provider import (
    BudgetExhausted, JevProvider, ProviderError, batches, durable_write,
    payload, validate_response,
)

MODEL = "jev-1.14.0"


def rows(n=4):
    return [dict(cell_id=f"c{i}", actions={"REST": "Remain at rest", "ACTIVATE": "Activate"},
                 prior={"REST": .75, "ACTIVATE": .25},
                 observation={"cell_id": f"c{i}", "step": 0, "state": {"activated": False}})
            for i in range(n)]


def context():
    return dict(schema_version="2.0", tissue="brain", organism="human", disease="control",
                assay="xenium", sampling=dict(design="subsampled", description="test region",
                missing_populations=[]), perturbations=[], recognition=False,
                time_step_hours=1, interaction_mode="spatial")


def response(body):
    return dict(model=body["model"], answers={qid: dict(type="choice", choice="REST",
                probabilities={"REST": .75, "ACTIVATE": .25}, confidence=.7)
                for qid in body["questions"]}, usage=dict(input_tokens=10, output_tokens=5))


def test_full_5k_panel_is_hashed_without_overflowing_the_prompt():
    c = context()
    c["assayed_genes"] = [f"HUMAN_PANEL_FEATURE_{i:05d}" for i in range(5001)]
    body = payload(rows(1), c, MODEL)
    assert len(canonical(c)) > 40000
    assert len(canonical(body)) <= 40000
    assert "assayed_genes" not in body["state"]["context"]
    assert body["state"]["full_context_sha256"] == digest(c)
    assert len(batches(rows(1), c, MODEL)) == 1
    changed = deepcopy(c)
    changed["assayed_genes"][0] = "DIFFERENT_FEATURE"
    assert payload(rows(1), changed, MODEL) != body


def test_requests_run_concurrently_then_replay_without_another_paid_call(tmp_path):
    barrier, lock, calls = threading.Barrier(2), threading.Lock(), []

    def transport(body):
        with lock:
            calls.append(body)
        barrier.wait(timeout=10)
        return response(body)

    p = JevProvider("FAKE_KEY", tmp_path, {"run": 1}, model=MODEL,
                    max_attempts=2, workers=2, max_questions=2, transport=transport)
    answers = p.choose(rows(), context())
    assert set(answers) == {"c0", "c1", "c2", "c3"}
    assert len(calls) == 2
    assert p.usage()["http_attempts"] == 2
    assert p.usage()["input_tokens"] == 20

    def must_not_send(body):
        pytest.fail("A cached, successful request was sent again")

    resumed = JevProvider("FAKE_KEY", tmp_path, {"run": 1}, model=MODEL,
                          max_attempts=2, max_questions=2, transport=must_not_send)
    assert resumed.choose(rows(), context()) == answers
    assert resumed.usage()["http_attempts"] == 2
    assert resumed.usage()["cache_hits"] == 2


def test_complete_round_budget_prevents_any_partial_dispatch(tmp_path):
    calls = []
    p = JevProvider("FAKE_KEY", tmp_path, {}, model=MODEL, max_attempts=1,
                    max_questions=2, transport=lambda body: calls.append(body))
    with pytest.raises(BudgetExhausted, match="complete synchronous round"):
        p.choose(rows(), context())
    assert calls == []
    assert list(tmp_path.glob("*.intent.json")) == []
    assert p.usage()["http_attempts"] == 0


def test_unresolved_intent_is_not_resent(tmp_path):
    p = JevProvider("FAKE_KEY", tmp_path, {}, model=MODEL, max_attempts=2,
                    transport=lambda body: pytest.fail("Unresolved attempt was resent"))
    body = payload(rows(1), context(), MODEL)
    h = digest(body)
    durable_write(tmp_path / f"{h}.intent.json",
                  dict(identity=p.identity, request_hash=h, request=body))
    resumed = JevProvider("FAKE_KEY", tmp_path, {}, model=MODEL, max_attempts=2,
                          transport=lambda body: pytest.fail("Unresolved attempt was resent"))
    with pytest.raises(ProviderError, match="Unresolved Jev intent"):
        resumed.choose(rows(1), context())
    assert resumed.usage()["http_attempts"] == 1
    assert resumed.usage()["attempts_with_unknown_usage"] == 1


def test_failed_attempt_is_redacted_and_never_silently_retried(tmp_path):
    key = "TEST_SECRET_DO_NOT_PERSIST"

    def transport(body):
        raise ProviderError(f"Provider rejected {key}")

    p = JevProvider(key, tmp_path, {}, model=MODEL, max_attempts=2, transport=transport)
    with pytest.raises(ProviderError, match="redacted"):
        p.choose(rows(1), context())
    assert key not in "".join(path.read_text() for path in tmp_path.glob("*.json"))
    resumed = JevProvider(key, tmp_path, {}, model=MODEL, max_attempts=2,
                          transport=lambda body: pytest.fail("Failed request was retried"))
    with pytest.raises(ProviderError, match="previously failed"):
        resumed.choose(rows(1), context())


def test_changed_configuration_and_tampered_cache_are_rejected(tmp_path):
    p = JevProvider("", tmp_path, {"seed": 1}, model=MODEL, transport=response)
    p.choose(rows(1), context())
    with pytest.raises(ProviderError, match="another run"):
        JevProvider("", tmp_path, {"seed": 2}, model=MODEL, transport=response)
    path = next(tmp_path.glob("*.result.json"))
    record = json.loads(path.read_text())
    record["response"]["answers"]["q0"]["confidence"] = .1
    path.write_text(json.dumps(record))
    with pytest.raises(ProviderError, match="result changed"):
        p.choose(rows(1), context())


@pytest.mark.parametrize("fault", ["model", "questions", "options", "nan", "choice", "confidence"])
def test_invalid_provider_answer_fails_closed(fault):
    body = payload(rows(1), context(), MODEL)
    raw = response(body)
    answer = raw["answers"]["q0"]
    if fault == "model":
        raw["model"] = "jev-0.0.0"
    elif fault == "questions":
        raw["answers"]["unexpected"] = deepcopy(answer)
    elif fault == "options":
        answer["probabilities"]["PROLIFERATE"] = 0
    elif fault == "nan":
        answer["probabilities"]["REST"] = float("nan")
    elif fault == "choice":
        answer["choice"] = "ACTIVATE"
    else:
        answer["confidence"] = True
    with pytest.raises(ProviderError):
        validate_response(raw, body)


def test_hidden_prior_mode_omits_prior_and_keeps_its_own_journal(tmp_path):
    shown = payload(rows(1), context(), MODEL, "shown")
    hidden = payload(rows(1), context(), MODEL, "hidden")
    assert "guidebook_prior" in shown["questions"]["q0"]["instructions"]
    assert "guidebook_prior" not in hidden["questions"]["q0"]["instructions"]
    assert hidden["state"]["prompt_version"] != shown["state"]["prompt_version"]
    assert "pooled" not in hidden["state"]["semantics"]
    with pytest.raises(ValueError):
        payload(rows(1), context(), MODEL, "anchored")
    sent = []
    p = JevProvider("FAKE_KEY", tmp_path / "a", {"run": 1}, model=MODEL, max_attempts=4,
                    transport=lambda body: sent.append(body) or response(body), prior_mode="hidden")
    p.choose(rows(2), context())
    assert all("guidebook_prior" not in q["instructions"] for q in sent[0]["questions"].values())
    # The same journal directory cannot be reused by the other prompt mode.
    with pytest.raises(ProviderError):
        JevProvider("FAKE_KEY", tmp_path / "a", {"run": 1}, model=MODEL, max_attempts=4,
                    transport=response, prior_mode="shown")


def test_cost_estimate_reports_budget_coverage_before_any_call():
    from virtualtissue_xenium.v2.provider import estimate_cost
    many = rows(50)
    estimate = estimate_cost(many, context(), model=MODEL, steps=12, max_questions=20, max_http_attempts=10)
    assert estimate["http_calls_per_round"] == 3 and estimate["http_calls_total"] == 36
    assert estimate["steps_covered_by_budget"] == 3 and estimate["budget_sufficient"] is False
    assert estimate["input_tokens_total_estimate"] > 0
    assert estimate_cost(many, context(), steps=12, max_http_attempts=36)["budget_sufficient"] is True
    hidden = estimate_cost(many, context(), steps=12, prior_mode="hidden")
    assert hidden["bytes_per_round"] < estimate["bytes_per_round"]
