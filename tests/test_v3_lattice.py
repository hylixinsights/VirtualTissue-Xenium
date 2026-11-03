"""Binned units use fixed occupancy, never secretly become moving cell agents."""
from copy import deepcopy
import json

import pytest

from virtualtissue_xenium.v2.inputs import prepare
from virtualtissue_xenium.v2.lattice import validate_lattice
from virtualtissue_xenium.v2.simulation import commit, decision_rows, initialize, interaction_graph, run


def prepared(*, lattice=None, cell_updates=None):
    context = dict(schema_version="3.0", tissue="skin", organism="human", disease="fixture", assay="visium_hd",
                   kernel="lattice", interaction_mode="spatial", sampling=dict(design="subsampled", description="test bins", missing_populations=[]),
                   perturbations=[], recognition=False, time_step_hours=.5, lattice=lattice or {})
    cells = [dict(cell_id="a", cell_type="tumor", lineage="malignant", sample_id="s", x=0, y=0,
                  expression={"ACTB": 100}, purity=1, occupancy=1),
             dict(cell_id="b", cell_type="empty", lineage="unknown", sample_id="s", x=8, y=0,
                  expression={"ACTB": 100}, purity=1, occupancy=0),
             dict(cell_id="c", cell_type="ambiguous", lineage="malignant", sample_id="s", x=16, y=0,
                  expression={"ACTB": 100}, purity=.2)]
    for cid, changes in (cell_updates or {}).items():
        next(c for c in cells if c["cell_id"] == cid).update(changes)
    return prepare(cells, context, assayed_genes=["ACTB"],
                   measurement=dict(unit="binned_spot", unit_size_um=8, deconvolution_method="synthetic fixture"))


def test_unknown_measured_bin_is_occupied_while_only_explicit_zero_is_empty():
    bundle = prepared()
    cells = initialize(bundle)
    assert cells["b"]["status"] == "removed" and cells["b"]["occupancy"] == 0
    assert cells["c"]["lineage"] == "unknown" and cells["c"]["status"] == "alive" and cells["c"]["occupancy"] == 1
    rows, graph = decision_rows(cells, bundle["guidebook"], bundle["context"], seed=2, step=0)
    assert {k: set(v) for k, v in graph.items()} == {"a": {"b"}, "b": {"a", "c"}, "c": {"b"}}
    assert next(r for r in rows if r["cell_id"] == "a")["division_targets"] == ["b"]
    assert set(next(r for r in rows if r["cell_id"] == "c")["actions"]) == {"REST"}


def test_growth_populates_existing_bin_without_creating_or_moving_ids():
    bundle = prepared()
    cells = initialize(bundle)
    rows, _ = decision_rows(cells, bundle["guidebook"], bundle["context"], seed=2, step=0)
    a = next(r for r in rows if r["cell_id"] == "a")
    proposal = dict(cell_id="a", action="PROLIFERATE", targets=[], corpses=[], division_targets=a["division_targets"])
    updated, events = commit(cells, [proposal], max_population=3, seed=2, step=0, context=bundle["context"])
    assert set(updated) == set(cells)
    assert all((updated[c]["x"], updated[c]["y"]) == (cells[c]["x"], cells[c]["y"]) for c in cells)
    assert updated["b"]["lineage"] == "malignant" and updated["b"]["guidebook_id"] == "a"
    assert updated["b"]["bin_lineage"] == "unknown" and updated["b"]["occupancy"] == 1
    assert events == [dict(kind="occupancy_gain", cell_id="a", target="b", lineage="malignant", added_occupancy=1)]
    rows, _ = decision_rows(updated, bundle["guidebook"], bundle["context"], seed=2, step=1)
    assert all("PROLIFERATE" not in r["actions"] for r in rows)


def test_loss_empties_bin_and_clocks_and_growth_conflicts_respect_capacity():
    bundle = prepared(cell_updates={"c": dict(cell_type="tumor", lineage="malignant", purity=1)})
    state = initialize(bundle)
    proposals = [dict(cell_id=cid, action="PROLIFERATE", targets=[], corpses=[], division_targets=["b"]) for cid in ("a", "c")]
    updated, events = commit(state, proposals, max_population=3, seed=2, step=0, context=bundle["context"])
    assert updated["b"]["occupancy"] == 1
    assert sum(e["kind"] == "occupancy_gain" for e in events) == 1
    assert sum(e["kind"] == "occupancy_capacity_blocked" for e in events) == 1
    emptied, events = commit(updated, [dict(cell_id="b", action="APOPTOSIS", targets=[], corpses=[])],
                             max_population=3, seed=2, step=1, context=bundle["context"])
    assert emptied["b"]["status"] == "removed" and emptied["b"]["occupancy"] == 0 and emptied["b"]["signal"] == 0
    assert events[0]["kind"] == "occupancy_loss" and events[0]["cause"] == "apoptosis"


def test_partial_same_lineage_occupancy_only_with_declared_increment_and_threshold():
    bundle = prepared(lattice=dict(division_increment=.5, low_density_threshold=.5),
                      cell_updates={"b": dict(cell_type="tumor", lineage="malignant", occupancy=.25)})
    state = initialize(bundle)
    rows, _ = decision_rows(state, bundle["guidebook"], bundle["context"], seed=2, step=0)
    assert next(r for r in rows if r["cell_id"] == "a")["division_targets"] == ["b"]
    updated, _ = commit(state, [dict(cell_id="a", action="PROLIFERATE", division_targets=["b"], targets=[], corpses=[])],
                        max_population=3, seed=2, step=0, context=bundle["context"])
    assert updated["b"]["occupancy"] == .75


def test_lattice_run_reports_occupancy_separately_and_has_streamed_frames(tmp_path):
    bundle = prepared()
    report = run(bundle, tmp_path, steps=4, seed=2)
    assert report["kernel"] == "lattice" and report["completed"]
    assert "occupancy" in report and "not cells" in report["interpretation"]
    frames = json.loads((tmp_path / "frames.json").read_text())
    assert all(set(f["cells"]) == {"a", "b", "c"} for f in frames)
    assert all(all((c["x"], c["y"]) == (frames[0]["cells"][cid]["x"], frames[0]["cells"][cid]["y"])
                   for cid, c in f["cells"].items()) for f in frames)
    assert not any(e["kind"] in {"division", "move", "apoptosis", "kill"} for f in frames for e in f["events"])


def test_lattice_rejects_invalid_capacity_or_duplicate_coordinates():
    for settings in (dict(capacity=0), dict(capacity=1, division_increment=2), dict(low_density_threshold=-1)):
        with pytest.raises(ValueError):
            validate_lattice(settings)
    bundle = prepared(cell_updates={"b": dict(x=0)})
    with pytest.raises(ValueError, match="unique centroid"):
        interaction_graph(initialize(bundle), bundle["context"], seed=1)
