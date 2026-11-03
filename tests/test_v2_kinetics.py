"""Lane A: per-hour hazards, first-event discretisation and time-step invariance."""
from copy import deepcopy
import math

import pytest

from virtualtissue_xenium.v2 import simulation
from virtualtissue_xenium.v2.guidebook import KineticsStepError, load_kinetics
from virtualtissue_xenium.v2.inputs import prepare
from virtualtissue_xenium.v2.simulation import decision_rows, initialize, preflight, run

PANEL = ["ACTB", "GAPDH", "B2M"]   # no action-linked signature reaches the 3-gene overlap, so multipliers are 1


def context(dt=1.0, mode="hazard_per_hour", **extra):
    c = dict(schema_version="2.0", tissue="skin", organism="human", disease="melanoma", assay="scrna",
             interaction_mode="well_mixed", sampling=dict(design="subsampled", description="fixture",
             missing_populations=[]), perturbations=[], recognition=True, time_step_hours=dt,
             kinetics=dict(mode=mode))
    c.update(extra)
    return c


def population(n=200, lineage="malignant", label="mal"):
    return [dict(cell_id=f"c{i:03d}", cell_type=label, lineage=lineage, sample_id="s", expression={"ACTB": 5.0})
            for i in range(n)]


def living_at(bundle, seed, steps, tmp_path, name):
    report = run(bundle, tmp_path / name, steps=steps, seed=seed)
    assert report["completed"]
    return report["living"]


def test_hazard_table_matches_licensed_actions_and_converts_to_first_event_probabilities():
    table, hazards = load_kinetics()
    bundle = prepare(population(1), context(dt=2.0), assayed_genes=PANEL)
    actions = bundle["guidebook"]["cells"]["c000"]["actions"]
    k = hazards["malignant"]
    burden = sum(k.values()) * 2.0
    assert math.isclose(actions["REST"]["probability"], math.exp(-burden))
    assert math.isclose(actions["PROLIFERATE"]["probability"], (1 - math.exp(-burden)) * k["PROLIFERATE"] / sum(k.values()))
    assert actions["PROLIFERATE"]["hazard_per_hour"] == k["PROLIFERATE"]
    assert math.isclose(sum(a["probability"] for a in actions.values()), 1.0)
    prov = bundle["guidebook"]["provenance"]
    assert prov["kinetics_mode"] == "hazard_per_hour" and prov["kinetics_sha256"] and prov["kinetics_version"]
    assert "DRAFT" in prov["kinetics_review_status"]
    assert all(set(rows) <= {"hazard_per_hour", "confidence", "anchor", "note"} and rows["confidence"] in table["confidence_levels"]
               for lineage in table["hazards"].values() for rows in lineage.values())
    assert all(rows.get("anchor", "ctl_serial_killing") in table["anchors"] for lineage in table["hazards"].values() for rows in lineage.values())


def test_time_step_changes_hazard_mode_but_not_legacy_mode():
    fine = prepare(population(1), context(dt=0.5), assayed_genes=PANEL)["guidebook"]["cells"]["c000"]["actions"]
    coarse = prepare(population(1), context(dt=2.0), assayed_genes=PANEL)["guidebook"]["cells"]["c000"]["actions"]
    assert coarse["PROLIFERATE"]["probability"] > fine["PROLIFERATE"]["probability"]
    legacy_fine = prepare(population(1), context(dt=0.5, mode="legacy_per_decision"), assayed_genes=PANEL)
    legacy_coarse = prepare(population(1), context(dt=2.0, mode="legacy_per_decision"), assayed_genes=PANEL)
    assert legacy_fine["guidebook"]["cells"] == legacy_coarse["guidebook"]["cells"]
    assert legacy_fine["guidebook"]["provenance"]["kinetics_mode"] == "legacy_per_decision"


def test_coarse_step_is_refused_with_a_recommended_time_step():
    with pytest.raises(ValueError, match="time_step_hours <="):
        prepare(population(1, "myeloid", "mac"), context(dt=4.0), assayed_genes=PANEL)
    with pytest.raises(ValueError, match="REST is residual"):
        prepare(population(1), context(dt=1.0, action_priors={"malignant": {"REST": 2.0}}), assayed_genes=PANEL)
    with pytest.raises(ValueError, match="kinetics.mode"):
        prepare(population(1), context(dt=1.0, kinetics={"mode": "per_minute"}), assayed_genes=PANEL)
    with pytest.raises(ValueError, match="hazard_per_hour"):
        prepare(population(1), context(dt=1.0, kinetics={"mode": "legacy_per_decision", "max_step_probability": 0.3}), assayed_genes=PANEL)
    facts = preflight(prepare(population(1, "myeloid", "mac"), context(dt=0.5), assayed_genes=PANEL))
    kin = facts["kinetics"]
    assert kin["mode"] == "hazard_per_hour" and kin["gated_worst_case_step_burden"] > kin["max_cell_event_probability_t0"]
    assert 0 < kin["recommended_max_time_step_hours"] < 1
    assert math.isclose(preflight(prepare(population(1), context(dt=1.0), assayed_genes=PANEL))["kinetics"]
                        ["by_lineage"]["malignant"]["implied_doubling_hours_ignoring_death"], math.log(2) / 0.006)


def test_gated_actions_return_their_share_to_rest_and_the_kernel_guards_the_live_menu():
    cells = [dict(cell_id="cd8", cell_type="cd8", lineage="T_CD8", sample_id="s", expression={"ACTB": 1.0},
                  initial_state={"activated": True}),
             dict(cell_id="tum", cell_type="mal", lineage="malignant", sample_id="s", expression={"ACTB": 1.0})]
    bundle = prepare(cells, context(dt=0.5), assayed_genes=PANEL)
    state = initialize(bundle)
    rows, _ = decision_rows(state, bundle["guidebook"], bundle["context"], seed=1, step=0)
    cd8 = next(r for r in rows if r["cell_id"] == "cd8")
    available = {a: bundle["guidebook"]["cells"]["cd8"]["actions"][a]["hazard_per_hour"] for a in cd8["actions"] if a != "REST"}
    assert "KILL" in available and "ACTIVATE" not in available
    assert math.isclose(cd8["prior"]["REST"], math.exp(-sum(available.values()) * 0.5))
    # Remove the target: KILL leaves the menu and its hazard goes back to REST, not to the other actions.
    state["tum"]["status"] = "removed"
    rows, _ = decision_rows(state, bundle["guidebook"], bundle["context"], seed=1, step=0)
    cd8 = next(r for r in rows if r["cell_id"] == "cd8")
    assert "KILL" not in cd8["actions"]
    assert math.isclose(cd8["prior"]["REST"], math.exp(-(sum(available.values()) - available["KILL"]) * 0.5))
    tight = prepare(cells, context(dt=0.5, kinetics=dict(mode="hazard_per_hour", max_step_probability=0.1)), assayed_genes=PANEL)
    with pytest.raises(KineticsStepError):
        decision_rows(initialize(tight), tight["guidebook"], tight["context"], seed=1, step=0)


def test_run_stops_with_a_report_when_a_step_becomes_too_coarse(tmp_path, monkeypatch):
    bundle = prepare(population(3), context(dt=1.0), assayed_genes=PANEL)
    original = simulation.decision_rows

    def failing(cells, guidebook, context, *, seed, step):
        if step == 1:
            raise KineticsStepError("Cell 'c000' at step 1: summed hazard x time step = 0.9 exceeds kinetics.max_step_probability 0.25")
        return original(cells, guidebook, context, seed=seed, step=step)

    monkeypatch.setattr(simulation, "decision_rows", failing)
    report = run(bundle, tmp_path / "coarse", steps=3, seed=1)
    assert report["completed"] is False and report["stop_reason"] == "kinetics_step_too_coarse"
    assert report["completed_steps"] == 1 and "time step" in report["error"]


def test_time_step_invariance_in_hazard_mode_and_step_dependence_in_legacy_mode(tmp_path):
    horizon, seeds = 48.0, range(8)
    means = {}
    for dt in (0.5, 1.0, 2.0, 4.0):
        bundle = prepare(population(200), context(dt=dt), assayed_genes=PANEL)
        steps = int(horizon / dt)
        means[dt] = sum(living_at(bundle, s, steps, tmp_path, f"h{dt}_{s}") for s in seeds) / len(seeds)
    _, hazards = load_kinetics()
    k = hazards["malignant"]
    expected = 200 * math.exp((k["PROLIFERATE"] - k["APOPTOSIS"]) * horizon)
    for dt, mean in means.items():
        # One event per cell per step under-counts rates by about burden / 2 (documented); allow that plus noise.
        bias = sum(k.values()) * dt / 2
        assert abs(mean - expected) / expected < bias + 0.03, (dt, mean, expected, means)
    assert abs(means[1.0] - means[0.5]) / means[0.5] < 0.04, means
    # Legacy mode: the same number of steps gives the same trajectory whatever time_step_hours says.
    legacy = {dt: living_at(prepare(population(60), context(dt=dt, mode="legacy_per_decision"), assayed_genes=PANEL),
                            3, 12, tmp_path, f"l{dt}") for dt in (0.5, 4.0)}
    assert legacy[0.5] == legacy[4.0]
