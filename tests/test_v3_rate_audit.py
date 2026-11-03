"""Scientific metadata corrections preserve rates and describe actual event units."""
from copy import deepcopy
import hashlib
import json

import pytest

from virtualtissue_xenium.v2.inputs import prepare
from virtualtissue_xenium.v2.simulation import commit, decision_rows, initialize
from virtualtissue_xenium.v3.ontology import load_ontology


def context(**overrides):
    value = dict(schema_version="3.0", tissue="skin", organism="human", disease="rate audit fixture",
                 assay="xenium", recognition=True, perturbations=[], time_step_hours=.1,
                 sampling=dict(design="subsampled", description="synthetic rate audit", missing_populations=[]),
                 scoring=dict(beta=0), dynamics=dict(crowding_enabled=False))
    value.update(overrides)
    return value


def cell(cid, lineage, x, **extra):
    return dict(cell_id=cid, cell_type=lineage, lineage=lineage, sample_id="s", x=x, y=0,
                expression={"ACTB": 100}, **extra)


def rows(bundle, state, step=0):
    result, _ = decision_rows(state, bundle["guidebook"], bundle["context"], seed=1, step=step)
    return {r["cell_id"]: r for r in result}


def test_scientific_metadata_audit_did_not_change_any_numeric_default():
    table = load_ontology()["kinetics"]
    projection = dict(hazards={lineage: {action: {k: row[k] for k in ("hazard_per_hour", "beta")}
                                       for action, row in records.items()}
                               for lineage, records in table["hazards"].items()},
                      beta_by_action=table["beta_by_action"])
    # Frozen from v3-kinetics-2 before editing the scientific annotations. This
    # independent baseline detects any accidental parameter drift in the audit.
    assert hashlib.sha256(json.dumps(projection, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest() == (
        "31d02e512cc7ab52dab2277ba27666268140a852fd0465e0d9dc6b8edaf28a5c")
    records = [r for values in table["hazards"].values() for r in values.values()]
    assert len(records) == 91
    assert all(r["confidence"] == r["beta_confidence"] == "assumed" for r in records)
    assert all(not r["evidence_mapping"]["parameter_estimated_from_cited_data"] for r in records)
    assert all(set(r["evidence_mapping"]["anchors"]) <= table["anchors"].keys() for r in records)


def test_cd8_activation_is_a_gate_not_a_twentyfold_proliferation_multiplier():
    bundle = prepare([cell("cd8", "T_CD8", 0), cell("nk", "NK", 5), cell("tumor", "malignant", 10)],
                     context(), assayed_genes=["ACTB"])
    state = initialize(bundle)
    inactive = rows(bundle, state)
    assert "PROLIFERATE" not in inactive["cd8"]["actions"]
    assert inactive["nk"]["effective_hazards_per_hour"]["PROLIFERATE"] == .004
    for cid in ("cd8", "nk"):
        assert "KILL" not in inactive[cid]["actions"]
        state[cid].update(activated=True, activated_since=0, activated_until=200)
    active = rows(bundle, state)
    for cid in ("cd8", "nk"):
        assert active[cid]["effective_hazards_per_hour"]["PROLIFERATE"] == .004
        assert active[cid]["effective_hazards_per_hour"]["KILL"] == .2
        state[cid].update(exhausted=True, exhausted_since=0)
    terminal = rows(bundle, state, step=1680)
    for cid in ("cd8", "nk"):
        assert "RECOVER" not in terminal[cid]["actions"]
        assert "KILL" not in terminal[cid]["actions"]
        assert {"PROLIFERATE", "ACTIVATE", "SECRETE", "MOVE"} <= terminal[cid]["actions"].keys()
        assert terminal[cid]["effective_hazards_per_hour"]["PROLIFERATE"] == .004


def test_lattice_hazard_is_per_occupied_bin_and_loss_removes_fractional_occupancy():
    bundle = prepare([cell("bin", "malignant", 0, purity=1, occupancy=.25)],
                     context(assay="visium_hd", kernel="lattice", recognition=False),
                     assayed_genes=["ACTB"],
                     measurement=dict(unit="binned_spot", unit_size_um=8,
                                      deconvolution_method="synthetic audit fixture"))
    partial = initialize(bundle)
    full = deepcopy(partial)
    full["bin"]["occupancy"] = 1
    a, b = rows(bundle, partial)["bin"], rows(bundle, full)["bin"]
    assert a["effective_hazards_per_hour"]["APOPTOSIS"] == .002
    assert a["prior"] == b["prior"]  # Neither hazard nor menu is occupancy-scaled.
    updated, events = commit(partial, [dict(cell_id="bin", action="APOPTOSIS", targets=[], corpses=[])],
                             max_population=10, seed=1, step=0, context=bundle["context"])
    assert updated["bin"]["occupancy"] == 0
    assert updated["bin"]["status"] == "removed"
    assert events == [dict(kind="occupancy_loss", cell_id="bin", target="bin", lineage="malignant",
                           removed_occupancy=.25, cause="apoptosis")]
