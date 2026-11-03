"""State clocks, spatial motion and streamed snapshots have explicit v3 semantics."""
from copy import deepcopy
import json
import math
import re

import pytest

from virtualtissue_xenium.v2.dynamics import validate_dynamics
from virtualtissue_xenium.v2.frame_store import iter_frames
from virtualtissue_xenium.v2.inputs import prepare
from virtualtissue_xenium.v2.simulation import commit, decision_rows, initialize, interaction_graph, run


def context(**overrides):
    value = dict(schema_version="3.0", tissue="skin", organism="human", disease="engineering fixture",
                 assay="xenium", sampling=dict(design="subsampled", description="synthetic fixture", missing_populations=[]),
                 perturbations=[], recognition=True, time_step_hours=.5)
    value.update(overrides)
    return value


def cells():
    return [dict(cell_id=cid, cell_type=lineage, lineage=lineage, sample_id="s", expression={"ACTB": 100},
                 x=x, y=0, initial_state={"activated": activated})
            for cid, lineage, x, activated in [("tum", "malignant", 0, False), ("cd8", "T_CD8", 15, True),
                                               ("mac", "myeloid", 20, False)]]


def bundle(**overrides):
    return prepare(cells(), context(**overrides), assayed_genes=["ACTB"])


def proposal(cid, action):
    return dict(cell_id=cid, action=action, targets=[], corpses=[])


def test_dynamic_parameters_reject_unknown_nonfinite_or_boolean_numbers():
    assert validate_dynamics()["corpse_lifetime_hours"] == 72
    for kwargs in ({"persistence_hours": 0}, {"speed_um_per_hour": float("nan")}, {"speed_um_per_hour": True},
                   {"n_max": 0}, {"movement_enabled": "yes"}, {"unrecognised": 1}):
        with pytest.raises(ValueError):
            validate_dynamics(kwargs)


@pytest.mark.parametrize("mode", ["spatial", "well_mixed"])
def test_live_graph_is_symmetric_bounded_and_corpses_do_not_consume_slots(mode):
    b = bundle(max_neighbors=1, interaction_mode=mode)
    state = initialize(b)
    state["corpse1"] = dict(state["tum"], cell_id="corpse1", status="dead", x=5)
    state["corpse2"] = dict(state["tum"], cell_id="corpse2", status="dead", x=7)
    state["other_sample"] = dict(state["tum"], cell_id="other_sample", sample_id="other")
    graph = interaction_graph(state, b["context"], seed=3)
    assert graph["other_sample"] == []
    for cid, neighbors in graph.items():
        assert all(cid in graph[nid] for nid in neighbors)
        if state[cid]["status"] == "alive":
            assert sum(state[n]["status"] == "alive" for n in neighbors) <= 1
    assert any(any(state[n]["status"] == "alive" for n in ns) and any(state[n]["status"] == "dead" for n in ns)
               for cid, ns in graph.items() if state[cid]["status"] == "alive")


def test_corpses_carry_no_signal_and_degrade_after_declared_hours():
    b = bundle(dynamics={"corpse_lifetime_hours": 1})
    state = initialize(b)
    state["tum"].update(status="dead", dead_since=0, signal=1000)
    state["mac"]["signal"] = 4
    rows, _ = decision_rows(state, b["guidebook"], b["context"], seed=1, step=0)
    cd8 = next(r for r in rows if r["cell_id"] == "cd8")
    assert cd8["observation"]["local"]["incoming_signal"] == 4
    updated, _ = commit(state, [], max_population=10, seed=1, step=0, context=b["context"])
    assert updated["tum"]["status"] == "dead" and updated["tum"]["signal"] == 0
    updated, events = commit(updated, [], max_population=10, seed=1, step=1, context=b["context"])
    assert updated["tum"]["status"] == "removed"
    assert any(e["kind"] == "corpse_degraded" for e in events)


def test_activation_expires_and_daughters_inherit_clocks_instead_of_resetting():
    b = bundle()
    state = initialize(b)
    state["tum"].update(activated=True, activated_since=0, activated_until=.75, exhausted=True, exhausted_since=.1)
    updated, events = commit(state, [proposal("tum", "PROLIFERATE")], max_population=10, seed=4, step=0, context=b["context"])
    daughter_id = next(e["daughter"] for e in events if e["kind"] == "division")
    assert updated[daughter_id]["activated_until"] == .75
    assert updated[daughter_id]["exhausted_since"] == .1
    assert math.hypot(updated[daughter_id]["x"] - state["tum"]["x"], updated[daughter_id]["y"] - state["tum"]["y"]) >= 9.9
    final, events = commit(updated, [], max_population=10, seed=4, step=1, context=b["context"])
    assert not final["tum"]["activated"] and not final[daughter_id]["activated"]
    assert {e["cell_id"] for e in events if e["kind"] == "deactivate"} >= {"tum", daughter_id}


def test_terminal_exhaustion_blocks_recovery_and_activation_can_restimulate():
    b = bundle(dynamics={"terminal_exhaustion_hours": 1})
    state = initialize(b)
    state["cd8"].update(exhausted=True, exhausted_since=0, activated_until=100)
    early, _ = decision_rows(state, b["guidebook"], b["context"], seed=2, step=1)
    late, _ = decision_rows(state, b["guidebook"], b["context"], seed=2, step=2)
    assert "RECOVER" in next(r for r in early if r["cell_id"] == "cd8")["actions"]
    c = next(r for r in late if r["cell_id"] == "cd8")
    assert "RECOVER" not in c["actions"] and "terminal" in c["excluded"]["RECOVER"]
    assert "ACTIVATE" in c["actions"]


def test_density_gate_and_inactive_lymphocyte_proliferation():
    b = bundle(dynamics={"n_max": 1})
    state = initialize(b)
    state["cd8"].update(x=1, activated=False)
    rows, _ = decision_rows(state, b["guidebook"], b["context"], seed=3, step=0)
    tumor = next(r for r in rows if r["cell_id"] == "tum")
    cd8 = next(r for r in rows if r["cell_id"] == "cd8")
    assert "PROLIFERATE" not in tumor["actions"] and tumor["excluded"]["PROLIFERATE"] == "local density gate"
    assert "PROLIFERATE" not in cd8["actions"]


def test_movement_integrates_hours_and_is_not_multiplied_by_time_step_twice():
    b = bundle(dynamics={"speed_um_per_hour": 10})
    state = initialize(b)
    state["cd8"].update(moving_until=1.0, heading_radians=0.0, activated_until=100)
    results = []
    for dt, n in [(.5, 2), (.25, 4)]:
        s = deepcopy(state)
        c = dict(b["context"], time_step_hours=dt)
        for step in range(n):
            s, _ = commit(s, [proposal("cd8", "REST")], max_population=10, seed=9, step=step, context=c)
        results.append(s["cd8"]["x"])
    assert results == [25, 25]
    disabled = deepcopy(b["context"])
    disabled["dynamics"]["movement_enabled"] = False
    s, events = commit(state, [proposal("cd8", "MOVE")], max_population=10, seed=9, step=0, context=disabled)
    assert s["cd8"]["x"] == state["cd8"]["x"] and not any(e["kind"] == "move" for e in events)


def test_move_starts_seeded_episode_and_can_bring_a_target_into_range():
    b = bundle(dynamics={"speed_um_per_hour": 20, "persistence_hours": 10})
    state = initialize(b)
    state["cd8"].update(x=-40, moving_until=10, heading_radians=0, activated_until=100)
    rows, _ = decision_rows(state, b["guidebook"], b["context"], seed=9, step=0)
    assert "KILL" not in next(r for r in rows if r["cell_id"] == "cd8")["actions"]
    moved, events = commit(state, [proposal("cd8", "REST")], max_population=10, seed=9, step=0, context=b["context"])
    rows, _ = decision_rows(moved, b["guidebook"], b["context"], seed=9, step=1)
    assert "KILL" in next(r for r in rows if r["cell_id"] == "cd8")["actions"]
    assert any(e["kind"] == "move" and e["distance_um"] == 10 for e in events)
    first = commit(initialize(b), [proposal("cd8", "MOVE")], max_population=10, seed=9, step=0, context=b["context"])
    second = commit(initialize(b), [proposal("cd8", "MOVE")], max_population=10, seed=9, step=0, context=b["context"])
    assert first == second and any(e["kind"] == "move" for e in first[1])


def test_spatial_target_ties_are_seeded_without_cell_id_bias():
    b = bundle()
    state = initialize(b)
    state["cd8"].update(x=0, y=0)
    state["tum"].update(x=5)
    state["tum2"] = dict(state["tum"], cell_id="tum2", x=-5)
    choices = []
    for seed in range(1000):
        rows, _ = decision_rows(state, b["guidebook"], b["context"], seed=seed, step=0)
        choices.append(next(r for r in rows if r["cell_id"] == "cd8")["targets"][0])
    assert 450 < choices.count("tum") < 550


def test_frames_stream_equal_array_and_replay_and_repeatability(tmp_path):
    b = bundle()
    first, second = tmp_path / "first", tmp_path / "second"
    report = run(b, first, steps=6, seed=21)
    repeated = run(b, second, steps=6, seed=21)
    assert report == repeated and report["completed"]
    frames = list(iter_frames(first / "frames.jsonl"))
    assert frames == json.loads((first / "frames.json").read_text())
    embedded = re.search(r'<script id="run-data" type="application/json">(.*?)</script>', (first / "replay.html").read_text(), re.S)
    replay_data = json.loads(embedded.group(1))
    assert all("cells" not in frame for frame in replay_data["frames"])
    loaded = []
    for chunk in replay_data["replay"]["chunks"]:
        content = (first / chunk["src"]).read_text()
        loaded.extend(json.loads(re.search(r"window\.__vtxReplayChunk\(\d+,(\[.*\])\);", content, re.S).group(1)))
    assert loaded == frames
    assert (first / "frames.jsonl").read_bytes() == (second / "frames.jsonl").read_bytes()
    assert not (first / ".running.lock").exists()
    assert not (first / "checkpoint.json").exists()


def test_horizon_guard_rejects_before_running(tmp_path):
    b = bundle(dynamics={"max_horizon_hours": 1})
    with pytest.raises(ValueError, match="horizon"):
        run(b, tmp_path, steps=3)
    assert not (tmp_path / "frames.jsonl").exists()


def test_expanding_clone_does_not_stack_daughters_at_two_micrometres():
    b = bundle(dynamics={"n_max": 100})
    state = {"tum": initialize(b)["tum"]}
    for step in range(4):
        state, _ = commit(state, [proposal(cid, "PROLIFERATE") for cid in sorted(state)],
                          max_population=30, seed=31, step=step, context=b["context"])
    assert len(state) == 16
    width = max(c["x"] for c in state.values()) - min(c["x"] for c in state.values())
    height = max(c["y"] for c in state.values()) - min(c["y"] for c in state.values())
    assert width * height > 400


def test_suppressive_secretion_does_not_act_as_inflammatory_activation():
    obj = cells() + [dict(cell_id="treg", cell_type="Treg", lineage="T_reg", sample_id="s", x=12, y=0, expression={"ACTB": 100})]
    b = prepare(obj, context(), assayed_genes=["ACTB"])
    state = initialize(b)
    base, _ = decision_rows(state, b["guidebook"], b["context"], seed=1, step=0)
    updated, _ = commit(state, [proposal("treg", "SECRETE")], max_population=10, seed=1, step=0, context=b["context"])
    assert updated["treg"]["signal"] == 0 and updated["treg"]["suppressive_signal"] == 1
    after, _ = decision_rows(updated, b["guidebook"], b["context"], seed=1, step=1)
    a = next(r for r in base if r["cell_id"] == "cd8")
    z = next(r for r in after if r["cell_id"] == "cd8")
    assert z["effective_hazards_per_hour"]["ACTIVATE"] < a["effective_hazards_per_hour"]["ACTIVATE"]


def test_report_rates_use_gated_t0_hazards_and_first_event_probability(tmp_path):
    b = bundle()
    state = initialize(b)
    rows, _ = decision_rows(state, b["guidebook"], b["context"], seed=3, step=0)
    report = run(b, tmp_path, steps=1, seed=3)
    rates = report["preflight"]["kinetics"]
    for row in rows:
        lineage = row["observation"]["lineage"]
        assert rates["gated_hazard_per_hour_by_lineage"][lineage] == row["effective_hazards_per_hour"]
        for action, probability in row["prior"].items():
            if action != "REST":
                assert math.isclose(report["expected_events_per_hour_by_lineage"][lineage][action], probability / .5)


@pytest.mark.parametrize("version", ["2.0", "3.0"])
def test_same_target_kill_and_apoptosis_commit_one_death_for_any_order(version):
    b = bundle() if version == "3.0" else prepare(cells(), context(schema_version="2.0"), assayed_genes=["ACTB"])
    state = initialize(b)
    proposals = [dict(cell_id="tum", action="APOPTOSIS", targets=[], corpses=[]),
                 dict(cell_id="cd8", action="KILL", targets=["tum"], corpses=[])]
    causes = set()
    for seed in range(30):
        updated, events = commit(state, proposals, max_population=10, seed=seed, step=0, context=b["context"])
        assert updated["tum"]["status"] == "dead"
        lethal = [e for e in events if e["kind"] in {"apoptosis", "kill"}]
        assert len(lethal) == 1
        causes.add(lethal[0]["kind"])
        assert sum(e["kind"] == "conflict_blocked" for e in events) == 1
    assert causes == {"kill", "apoptosis"}


def test_legacy_snapshot_remains_only_for_v2_compatibility(tmp_path):
    b = prepare(cells(), context(schema_version="2.0"), assayed_genes=["ACTB"])
    run(b, tmp_path, steps=1)
    assert (tmp_path / "checkpoint.json").exists()
