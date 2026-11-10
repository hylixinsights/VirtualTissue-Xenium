"""The v5 prompt and the opt-in pooling rule are versioned and leave the v4 defaults alone."""
import json
import math

import pytest

from virtualtissue_xenium.v2.inputs import canonical, prepare, read_json
from virtualtissue_xenium.v2.provider import (EVIDENCE_PROMPT_VERSION, JevProvider, ProviderError, payload,
                                              prompt_version)
from virtualtissue_xenium.v2.simulation import decision_rows, distribution, initialize, run
from v5_helpers import PANEL, FixedProvider, bundle, cells, context, rule


def rows_of(prepared):
    rows, _ = decision_rows(initialize(prepared), prepared["guidebook"], prepared["context"], seed=1, step=0)
    return [r for r in rows if len(r["actions"]) > 1]


# ---- prompt identity

def test_v4_prompt_identities_are_unchanged():
    assert prompt_version("shown", "2.0") == "vtx-context-guidebook-2.0"
    assert prompt_version("shown", "3.0") == "vtx-context-guidebook-3.1-evidence-basis"
    assert prompt_version("hidden", "3.0") == "vtx-context-guidebook-3.1-evidence-basis-hidden-prior"
    assert prompt_version("shown", "3.0", evidence_policy=True) == EVIDENCE_PROMPT_VERSION
    assert prompt_version("hidden", "3.0", evidence_policy=True, pooling="event_share").endswith("-event-share-hidden-prior")
    with pytest.raises(ValueError):
        prompt_version("shown", "2.0", evidence_policy=True)


def test_bundle_without_policy_gets_the_v4_request():
    plain = prepare(cells(scores=False), context(), assayed_genes=PANEL)
    body = payload(rows_of(plain), plain["context"], "jev-1.14.0")
    assert body["state"]["prompt_version"] == "vtx-context-guidebook-3.1-evidence-basis"
    assert "signature_definitions" in body["state"] and "evidence_definitions" not in body["state"]
    question = body["questions"]["q0"]["instructions"]["observation"]
    assert set(question["t0_evidence"]) == set(question["evidence_status"]) == set(question["evidence_basis"])


# ---- v5 evidence prompt

def test_policy_prompt_states_support_definitions_and_unavailable_evidence():
    prepared = bundle()
    rows = rows_of(prepared)
    body = payload(rows, prepared["context"], "jev-1.14.0")
    state = body["state"]
    assert state["prompt_version"] == EVIDENCE_PROMPT_VERSION
    definition = state["evidence_definitions"]["S_ACT"]
    assert definition == dict(kind="measured_score", genes_assayed=["CD69", "IL2RA"], set_size=12,
                              description="Activation proxy", caveat="Two genes only.")
    assert state["measured_scores"] == dict(method=["test-rank-score 1.0"], caveat="Synthetic scores.")
    assert state["evidence_definitions"]["cell_cycle"]["kind"] == "gene_set"
    # Signatures with no assayed genes on this panel are stated once, not repeated as nulls per cell.
    assert state["evidence_unavailable"]["hypoxia"] == ["insufficient_panel_overlap"]
    by_cell = {q["instructions"]["observation"]["cell_id"]: q["instructions"]["observation"] for q in body["questions"].values()}
    scored = by_cell["tumour05"]
    assert "hypoxia" not in scored["t0_evidence"] and "evidence_status" not in scored
    assert all(isinstance(z, float) and round(z, 2) == z for z in scored["t0_evidence"].values())
    assert scored["evidence_support"]["basis"] == "cell_type_x_depth_bin" and scored["evidence_support"]["n_ref"] == 12
    assert "cell_cycle" in scored["evidence_support"]["n_detected"]
    unscored = by_cell["tcell11"]
    assert "S_ACT" not in unscored["t0_evidence"] and unscored["evidence_unavailable"] == {"S_ACT": "not_measured"}
    assert "measured, not simulated" not in state["evidence_semantics"] and "not measurements" in state["evidence_semantics"]


def test_policy_prompt_does_not_reveal_the_rules():
    prepared = bundle()
    text = canonical(payload(rows_of(prepared), prepared["context"], "jev-1.14.0", "hidden")).decode()
    assert "Test rule." not in text and "evidence_policy" not in text and "requires_biological_review" not in text
    assert "guidebook_prior" not in text


def test_provider_evidence_can_restrict_what_is_sent():
    from v5_helpers import policy
    prepared = prepare(cells(), context(evidence_policy=policy(rule("PROLIFERATE", "S_CYCLE"), provider_evidence=["S_CYCLE"])),
                       assayed_genes=PANEL, measured_evidence=__import__("v5_helpers").manifest())
    body = payload(rows_of(prepared), prepared["context"], "jev-1.14.0")
    assert set(body["state"]["evidence_definitions"]) == {"S_CYCLE"}
    assert all(set(q["instructions"]["observation"]["t0_evidence"]) <= {"S_CYCLE"} for q in body["questions"].values())


def test_provider_and_context_must_agree_about_the_prompt(tmp_path):
    prepared = bundle()
    stale = JevProvider("", tmp_path / "a", {}, transport=lambda body: {}, schema_version="3.0")
    with pytest.raises(ProviderError, match="evidence policy"):
        stale.choose(rows_of(prepared), prepared["context"])
    matching = JevProvider("", tmp_path / "b", {}, transport=lambda body: {}, schema_version="3.0", evidence_policy=True)
    assert matching.prompt_version == EVIDENCE_PROMPT_VERSION


# ---- pooling

def test_default_pooling_is_the_v4_rule():
    prior, answer = {"REST": .9, "ACTIVATE": .1}, {"probabilities": {"REST": .8, "ACTIVATE": .2}}
    assert distribution(prior, answer, 1.0) == answer["probabilities"]
    pooled = distribution(prior, answer, .5)
    assert pooled["ACTIVATE"] == pytest.approx(math.sqrt(.02) / (math.sqrt(.02) + math.sqrt(.72)))
    assert distribution(prior, answer, .5, "full_distribution") == pooled
    with pytest.raises(ValueError):
        distribution(prior, answer, .5, "other")


@pytest.mark.parametrize("alpha", [.25, .5, 1.0])
def test_event_share_keeps_the_kernel_event_mass_at_every_time_step(alpha):
    answer = {"probabilities": {"REST": .6, "ACTIVATE": .3, "SECRETE": .1}}
    rates = []
    for dt in (.25, .5, 1.0):
        hazard = {"ACTIVATE": .08, "SECRETE": .02}
        mass = 1 - math.exp(-sum(hazard.values()) * dt)
        prior = {"REST": 1 - mass, **{a: mass * k / sum(hazard.values()) for a, k in hazard.items()}}
        applied = distribution(prior, answer, alpha, "event_share")
        assert applied["REST"] == prior["REST"] and sum(applied.values()) == pytest.approx(1)
        rates.append((1 - applied["REST"]) / dt)
        v4 = distribution(prior, answer, alpha)
        assert (1 - v4["REST"]) / dt > rates[-1]  # the v4 rule inflates events for this stated preference
        if alpha == 1:
            assert applied["ACTIVATE"] / applied["SECRETE"] == pytest.approx(3)  # Jev alone decides which event
    assert rates[0] == pytest.approx((1 - math.exp(-.025)) / .25)
    assert max(rates) / min(rates) < 1.06  # only the first-event discretisation remains
    v4_rates = [(1 - distribution({"REST": math.exp(-.1 * dt), "ACTIVATE": 1 - math.exp(-.1 * dt)},
                                  {"probabilities": {"REST": .6, "ACTIVATE": .4}}, 1.0)["REST"]) / dt for dt in (.25, 1.0)]
    assert v4_rates[0] == pytest.approx(4 * v4_rates[1])  # v4 hidden-prior rate scales as 1/dt


def test_event_share_falls_back_to_the_prior_when_the_provider_names_no_event():
    prior = {"REST": .9, "ACTIVATE": .06, "SECRETE": .04}
    assert distribution(prior, {"probabilities": {"REST": 1.0, "ACTIVATE": 0.0, "SECRETE": 0.0}}, 1.0, "event_share") == prior
    assert distribution(prior, {"probabilities": {"REST": .2, "ACTIVATE": .8, "SECRETE": 0.0}}, 0, "event_share") == prior
    only = distribution(prior, {"probabilities": {"REST": .2, "ACTIVATE": .8, "SECRETE": 0.0}}, 1.0, "event_share")
    assert only == {"REST": .9, "ACTIVATE": pytest.approx(.1), "SECRETE": 0.0}


def test_runs_report_event_mass_and_record_pooling_only_when_opted_in(tmp_path):
    plain = prepare(cells(scores=False), context(), assayed_genes=PANEL)
    v4 = run(plain, tmp_path / "v4", steps=2, seed=1, mode="jev", provider=FixedProvider(rest=.6))
    config = read_json(tmp_path / "v4/config.json")
    assert "jev_pooling" not in config and v4["jev_pooling"] == "full_distribution"
    mass = v4["jev_event_mass"]["by_lineage"]["T_CD8"]
    assert mass["jev"] == pytest.approx(.4) and mass["applied"] == pytest.approx(.4)
    assert mass["applied_over_prior"] > 5  # hidden prior: the stated REST preference set the event frequency

    shared = run(plain, tmp_path / "share", steps=2, seed=1, mode="jev", jev_pooling="event_share",
                 provider=FixedProvider(rest=.6, pooling="event_share"))
    config = read_json(tmp_path / "share/config.json")
    assert config["jev_pooling"] == "event_share" and "event mass" in config["probability_semantics"]
    mass = shared["jev_event_mass"]["by_lineage"]["T_CD8"]
    assert mass["applied"] == pytest.approx(mass["prior"]) and mass["applied_over_prior"] == pytest.approx(1)
    assert shared["config_hash"] != v4["config_hash"]
    for line in (tmp_path / "share/decisions.jsonl").read_text().splitlines():
        decision = json.loads(line)
        if decision["jev"]:
            assert decision["applied_probabilities"]["REST"] == decision["prior"]["REST"]


def test_pooling_requires_jev_mode_and_a_matching_provider(tmp_path):
    plain = prepare(cells(scores=False), context(), assayed_genes=PANEL)
    with pytest.raises(ValueError, match="only to Jev mode"):
        run(plain, tmp_path / "a", steps=1, jev_pooling="event_share")
    with pytest.raises(ValueError, match="disagree about jev_pooling"):
        run(plain, tmp_path / "b", steps=1, mode="jev", jev_pooling="event_share", provider=FixedProvider())
    with pytest.raises(ValueError, match="jev_pooling"):
        run(plain, tmp_path / "c", steps=1, mode="jev", jev_pooling="magic", provider=FixedProvider())


def test_event_share_prompt_tells_the_provider_what_its_answer_controls():
    plain = prepare(cells(scores=False), context(), assayed_genes=PANEL)
    rows = rows_of(plain)
    default = payload(rows, plain["context"], "jev-1.14.0", "hidden")
    shared = payload(rows, plain["context"], "jev-1.14.0", "hidden", "event_share")
    assert "simulator's clock" not in default["state"]["semantics"]
    assert "simulator's clock" in shared["state"]["semantics"]
    assert shared["state"]["prompt_version"].endswith("-event-share-hidden-prior")
    assert shared["questions"] == default["questions"]


def test_vtx5_is_an_alias_of_the_same_workflow():
    import tomllib
    from pathlib import Path
    import virtualtissue_xenium
    project = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())["project"]
    assert project["scripts"]["vtx5"] == project["scripts"]["vtx4"] == project["scripts"]["vtx2"]
    assert project["version"] == virtualtissue_xenium.__version__ == "5.0.0"
