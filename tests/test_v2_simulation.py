"""End-to-end v2 state, synchronous update and provenance boundaries."""
from copy import deepcopy
import json

import pytest

from virtualtissue_xenium.v2.inputs import digest, prepare, validate_prepared
from virtualtissue_xenium.v2.simulation import (
    commit, decision_rows, initialize, interaction_graph, run,
)


def context(assay="xenium"):
    return dict(schema_version="2.0", tissue="skin", organism="human", disease="melanoma",
                assay=assay, sampling=dict(design="subsampled", description="Test fixture region",
                missing_populations=["outside this region"]), perturbations=[], recognition=True,
                time_step_hours=1, max_neighbors=20)


def cells():
    return [dict(cell_id=cid, cell_type=label, lineage=lineage, expression={"ACTB": 10},
                 sample_id="sample_A", x=x, y=0, **extra)
            for cid, label, lineage, x, extra in [
                ("tumor", "reviewed_malignant", "malignant", 0, {}),
                ("cd8", "reviewed_CD8", "T_CD8", 5, {"initial_state": {"activated": True}}),
                ("macro", "reviewed_macrophage", "myeloid", 10, {}),
                ("neuron", "reviewed_neuron", "neuron", 12, {}),
            ]]


def rehash(bundle):
    bundle["fingerprint"] = digest({k: v for k, v in bundle.items() if k != "fingerprint"})
    return bundle


def test_seeded_trajectories_are_identical_after_cell_storage_reordering(tmp_path):
    a = prepare(cells(), context(), assayed_genes=["ACTB"])
    b = prepare(list(reversed(cells())), context(), assayed_genes=["ACTB"])
    assert a["fingerprint"] == b["fingerprint"]
    first = run(a, tmp_path / "first", steps=8, seed=21, max_population=20)
    second = run(b, tmp_path / "second", steps=8, seed=21, max_population=20)
    assert first == second
    for filename in ("frames.json", "decisions.jsonl"):
        assert (tmp_path / "first" / filename).read_bytes() == (tmp_path / "second" / filename).read_bytes()
    decisions = [json.loads(line) for line in (tmp_path / "first" / "decisions.jsonl").read_text().splitlines()]
    assert any(d["action"] != "REST" for d in decisions)
    assert not (tmp_path / "first" / ".running.lock").exists()


@pytest.mark.parametrize("mode", ["spatial", "well_mixed"])
def test_neighborhoods_never_cross_samples(mode):
    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    state = initialize(bundle)
    state["tumor"]["sample_id"] = "other_sample"
    c = deepcopy(bundle["context"])
    c["interaction_mode"] = mode
    graph = interaction_graph(state, c, seed=7, step=0)
    assert graph["tumor"] == []
    assert all("tumor" not in neighbors for neighbors in graph.values())


def test_nuclei_without_coordinates_use_explicit_nonspatial_interactions(tmp_path):
    obj = cells()
    for cell in obj:
        cell.pop("x")
        cell.pop("y")
    bundle = prepare(obj, context("snrna"), assayed_genes=["ACTB"])
    assert bundle["context"]["interaction_mode"] == "well_mixed"
    run(bundle, tmp_path, steps=2, seed=7, max_population=10)
    frames = json.loads((tmp_path / "frames.json").read_text())
    assert all("x" not in c and "y" not in c for f in frames for c in f["cells"].values())
    spatial = deepcopy(bundle["context"])
    spatial["interaction_mode"] = "spatial"
    with pytest.raises(ValueError, match="requires measured x/y"):
        interaction_graph(initialize(bundle), spatial)


def test_rehashed_guidebook_cannot_license_neuron_division():
    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    bundle["guidebook"]["cells"]["neuron"]["actions"]["PROLIFERATE"] = dict(
        weight=1, description="Improperly added division action")
    with pytest.raises(ValueError, match="lineage-incompatible"):
        validate_prepared(rehash(bundle))


def test_local_gates_and_next_snapshot_signal_response():
    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    state = initialize(bundle)
    rows, _ = decision_rows(state, bundle["guidebook"], bundle["context"], seed=21, step=0)
    cd8 = next(r for r in rows if r["cell_id"] == "cd8")
    assert "KILL" in cd8["actions"]
    assert "CLEAR_CORPSE" not in next(r for r in rows if r["cell_id"] == "macro")["actions"]
    baseline_activation = next(r for r in rows if r["cell_id"] == "neuron")["prior"]["ACTIVATE"]
    state["macro"]["signal"] = 4
    state["cd8"]["exhausted"] = True
    changed, _ = decision_rows(state, bundle["guidebook"], bundle["context"], seed=21, step=1)
    assert "KILL" not in next(r for r in changed if r["cell_id"] == "cd8")["actions"]
    assert next(r for r in changed if r["cell_id"] == "neuron")["prior"]["ACTIVATE"] > baseline_activation
    state["tumor"]["status"] = "dead"
    final, _ = decision_rows(state, bundle["guidebook"], bundle["context"], seed=21, step=2)
    assert "tumor" not in {r["cell_id"] for r in final}
    assert "CLEAR_CORPSE" in next(r for r in final if r["cell_id"] == "macro")["actions"]


def test_same_round_death_blocks_division_and_conflicts_are_reproducible():
    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    state = initialize(bundle)
    original = deepcopy(state)
    proposals = [dict(cell_id="tumor", action="PROLIFERATE", targets=[], corpses=[]),
                 dict(cell_id="cd8", action="KILL", targets=["tumor"], corpses=[])]
    updated, events = commit(state, proposals, max_population=10, seed=4, step=0)
    assert state == original
    assert updated["tumor"]["status"] == "dead"
    assert len(updated) == len(state)
    assert not any(e["kind"] == "division" for e in events)
    assert any(e["kind"] == "actor_died_same_step" for e in events)
    assert commit(state, list(reversed(proposals)), max_population=10, seed=4, step=0) == (updated, events)


def test_two_phagocytes_cannot_clear_the_same_snapshot_corpse_twice():
    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    state = initialize(bundle)
    state["tumor"]["status"] = "dead"
    state["macro2"] = dict(state["macro"], cell_id="macro2")
    proposals = [dict(cell_id=cid, action="CLEAR_CORPSE", targets=[], corpses=["tumor"])
                 for cid in ("macro", "macro2")]
    updated, events = commit(state, proposals, max_population=10, seed=4, step=0)
    assert updated["tumor"]["status"] == "removed"
    assert sum(e["kind"] == "clear_corpse" for e in events) == 1
    assert sum(e["kind"] == "conflict_blocked" for e in events) == 1


def test_capacity_counts_living_cells_and_can_be_reused_after_death():
    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    state = initialize(bundle)
    proposals = [dict(cell_id="tumor", action="PROLIFERATE", targets=[], corpses=[]),
                 dict(cell_id="neuron", action="APOPTOSIS", targets=[], corpses=[])]
    updated, events = commit(state, proposals, max_population=len(state), seed=4, step=0)
    assert len(updated) == len(state)+1
    assert sum(c["status"] == "alive" for c in updated.values()) == len(state)
    assert any(e["kind"] == "division" for e in events)


def test_exclusive_lock_and_changed_context_prevent_unsafe_resume(tmp_path):
    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    lock = tmp_path / ".running.lock"
    lock.write_text("another process")
    with pytest.raises(ValueError, match="Run lock"):
        run(bundle, tmp_path, steps=1)
    assert lock.read_text() == "another process"
    lock.unlink()
    run(bundle, tmp_path, steps=1)
    changed_context = context()
    changed_context["disease"] = "different declared disease"
    changed = prepare(cells(), changed_context, assayed_genes=["ACTB"])
    with pytest.raises(ValueError, match="identical settings"):
        run(changed, tmp_path, steps=1, resume=True)
    assert not lock.exists()


def test_missing_jev_answers_do_not_silently_apply_local_choices(tmp_path):
    class BrokenProvider:
        model = "jev-1.14.0"

        def choose(self, rows, context):
            return {}

        def usage(self):
            return {"http_attempts": 0}

    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    report = run(bundle, tmp_path, steps=2, mode="jev", provider=BrokenProvider())
    assert report["completed"] is False
    assert report["completed_steps"] == 0
    assert report["stop_reason"] == "provider_stopped_before_round_commit"
    assert json.loads((tmp_path / "frames.json").read_text())[0]["cells"] == initialize(bundle)
    assert (tmp_path / "decisions.jsonl").read_text() == ""
