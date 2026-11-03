"""Seed-level analyses must preserve missingness, pairing and failed arms."""
from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from virtualtissue_xenium.v2 import ensemble
from virtualtissue_xenium.v2.inputs import prepare, write_json
from virtualtissue_xenium.v2.readouts import extract_readouts, load_readout_spec, spatial_readouts


def fixture_bundle(dt=.5):
    context = dict(schema_version="2.0", tissue="skin", organism="human", disease="melanoma",
        assay="xenium", sampling=dict(design="subsampled", description="ensemble engineering fixture",
        missing_populations=[]), perturbations=[], recognition=False, time_step_hours=dt,
        kinetics=dict(mode="hazard_per_hour"), scoring=dict(method="depth_conditioned_auc"))
    cells = [dict(cell_id="m", lineage="malignant", cell_type="mal", sample_id="s",
                  expression={"ACTB": 30}, x=0, y=0),
             dict(cell_id="t", lineage="T_CD8", cell_type="cd8", sample_id="s",
                  expression={"ACTB": 30}, x=4, y=3, initial_state=dict(activated=True))]
    return prepare(cells, context, assayed_genes=["ACTB"])


def read_csv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle))


def test_packaged_and_public_readout_contract_are_identical():
    spec, packaged_hash = load_readout_spec()
    public, public_hash = load_readout_spec(Path(__file__).parents[1] / "analysis/readouts.yaml")
    assert spec == public and public_hash == packaged_hash


def test_seed_quantile_envelope_is_distinct_from_median_uncertainty_and_missingness():
    summary = ensemble.summarize_values([*range(101), None], bootstrap_replicates=200, seed=7)
    assert summary["median"] == 50
    assert summary["seed_q025"] == 2.5 and summary["seed_q975"] == 97.5
    assert summary["median_ci95_low"] > summary["seed_q025"]
    assert summary["median_ci95_high"] < summary["seed_q975"]
    assert summary["n"] == 101 and summary["n_missing"] == 1
    one = ensemble.summarize_values([3, None], bootstrap_replicates=100)
    assert one["median"] == 3 and one["median_ci95_low"] is None
    assert ensemble.summarize_values([None], bootstrap_replicates=100)["median"] is None


def test_paired_standardization_uses_difference_sd_and_does_not_hide_degeneracy():
    # Mean=2, sample SD=1, hence dz=2; an independent-arm SD would differ.
    result = ensemble.standardized_paired_effect([1, 2, 3, None], bootstrap_replicates=200, seed=7)
    assert result["paired_dz"] == 2
    assert result["paired_mean_difference"] == 2 and result["paired_difference_sd"] == 1
    assert result["paired_dz_n"] == 3 and result["paired_dz_n_missing"] == 1
    assert result["paired_dz_bootstrap_undefined"] > 0
    assert result["paired_dz_ci95_low"] is None and result["paired_dz_ci_status"] == "undefined_bootstrap_ratios"
    for values in ([0, 0, 0], [4, 4, 4]):
        constant = ensemble.standardized_paired_effect(values, bootstrap_replicates=100)
        assert constant["paired_difference_sd"] == 0
        assert constant["paired_dz"] is None and constant["paired_dz_ci95_low"] is None
    one = ensemble.standardized_paired_effect([3, None], bootstrap_replicates=100)
    assert one["paired_dz"] is None and one["paired_dz_status"] == "insufficient_seeds"


def test_standardized_bootstrap_is_reproducible_and_sign_scale_equivariant():
    values = list(range(-10, 11))
    first = ensemble.standardized_paired_effect(values, bootstrap_replicates=200, seed=42)
    assert first == ensemble.standardized_paired_effect(values, bootstrap_replicates=200, seed=42)
    assert first["paired_dz"] == 0 and first["paired_dz_ci95_low"] < 0 < first["paired_dz_ci95_high"]
    assert first["paired_dz_bootstrap_valid"] == 200
    negative = ensemble.standardized_paired_effect([-3*v for v in values], bootstrap_replicates=200, seed=42)
    assert negative["paired_dz_ci95_low"] == pytest.approx(-first["paired_dz_ci95_high"])
    assert negative["paired_dz_ci95_high"] == pytest.approx(-first["paired_dz_ci95_low"])


def test_aggregate_pairs_before_standardizing_and_preserves_raw_medians():
    records = []
    for arm, values in (("baseline", [100, 1000, 10000]), ("contrast", [101, 1002, 10003])):
        for seed, value in enumerate(values):
            records.append(dict(arm=arm, seed=seed, completed=True,
                readouts=[dict(time_hours=1, readout="living_cells", lineage="ALL", value=value)]))
    result = ensemble.aggregate_readouts(records, bootstrap_replicates=100)["paired_differences"][0]
    assert result["median"] == 2 and result["paired_dz"] == 2
    assert result["paired_difference_sd"] == 1


def test_spatial_readouts_never_cross_samples_and_use_declared_denominators():
    def cell(lineage, sample, x, y, activated=False):
        return dict(lineage=lineage, sample_id=sample, x=x, y=y, status="alive", activated=activated)
    cells = {"m": cell("malignant", "s1", 0, 0), "t": cell("T_CD8", "s1", 3, 4, True),
             "other_m": cell("malignant", "s2", 3, 4), "lonely_t": cell("T_CD8", "s3", 0, 0, True)}
    result = spatial_readouts(cells, dict(interaction_mode="spatial", neighbor_radius_um=6))
    assert result["malignant_activated_cd8_contact_fraction"] == (.5, 2, "ok")
    assert result["mean_cd8_to_malignant_distance_um"] == (5, 1, "ok")
    mixed = spatial_readouts(cells, dict(interaction_mode="well_mixed"))
    assert all(row[0] is None and row[2] == "nonspatial" for row in mixed.values())


def test_readouts_count_committed_events_by_actor_and_target_and_stream_frames(tmp_path):
    initial = {"m": dict(lineage="malignant", sample_id="s", status="alive", activated=False, exhausted=False, x=0, y=0),
               "t": dict(lineage="T_CD8", sample_id="s", status="alive", activated=True, exhausted=False, x=1, y=0)}
    final = deepcopy(initial)
    final["m"]["status"] = "dead"
    frames = [dict(step=0, time_hours=0, cells=initial, events=[]),
              dict(step=1, time_hours=.5, cells=final, events=[dict(kind="kill", cell_id="t", target="m"),
                    dict(kind="conflict_blocked", cell_id="t", action="KILL")])]
    (tmp_path/"frames.jsonl").write_text("".join(json.dumps(frame)+"\n" for frame in frames))
    write_json(tmp_path/"report.json", dict(completed_steps=1))
    (tmp_path/"decisions.jsonl").write_text(json.dumps(dict(step=0, jev_prior_js_bits=.4,
        observation=dict(lineage="T_CD8")))+"\n")
    rows = extract_readouts(tmp_path, context=dict(interaction_mode="spatial", neighbor_radius_um=30))
    index = {(row["step"], row["readout"], row["lineage"]): row for row in rows}
    assert index[(1, "cumulative_kills", "T_CD8")]["value"] == 1
    assert index[(1, "cumulative_kills", "malignant")]["value"] == 0
    assert index[(1, "cumulative_killed_targets", "malignant")]["value"] == 1
    assert index[(1, "activated_fraction", "malignant")]["value"] is None
    assert index[(1, "jev_prior_js_bits", "T_CD8")]["value"] == .4
    assert index[(0, "jev_prior_js_bits", "T_CD8")]["value"] is None


def test_ensemble_runs_seed_arm_grid_reproducibly_and_pairs_by_seed(tmp_path):
    bundle = fixture_bundle()
    kwargs = dict(seeds=[3, 8], steps=3, max_population=20, bootstrap_replicates=100)
    first = ensemble.run_ensemble({"control": bundle, "same": bundle}, tmp_path/"first", **kwargs)
    second = ensemble.run_ensemble({"control": bundle, "same": bundle}, tmp_path/"second", **kwargs)
    assert first == second
    assert first["completed"] and first["n_completed_runs"] == 4
    assert (tmp_path/"first/summary.csv").read_bytes() == (tmp_path/"second/summary.csv").read_bytes()
    differences = read_csv(tmp_path/"first/paired_differences.csv")
    assert differences and all(float(row["median"]) == 0 for row in differences if row["median"])
    assert any(row["median"] == "" for row in differences if row["readout"] == "jev_prior_js_bits")
    assert first["warnings"]
    with pytest.raises(FileExistsError):
        ensemble.run_ensemble({"control": bundle}, tmp_path/"first", **kwargs)


def test_different_dt_requires_shared_horizon_and_no_interpolation(tmp_path):
    arms = {"fine": fixture_bundle(.5), "coarse": fixture_bundle(1)}
    with pytest.raises(ValueError, match="horizon_hours"):
        ensemble.run_ensemble(arms, tmp_path/"bad", steps=2)
    report = ensemble.run_ensemble(arms, tmp_path/"ok", horizon_hours=2, seeds=[1, 2], bootstrap_replicates=100)
    assert report["completed"]
    pairs = read_csv(tmp_path/"ok/paired_differences.csv")
    assert {float(row["time_hours"]) for row in pairs} == {0, 1, 2}
    with pytest.raises(ValueError, match="whole time steps"):
        ensemble.run_ensemble(arms, tmp_path/"bad_horizon", horizon_hours=.75)


def test_incomplete_seed_excludes_the_entire_arm_but_preserves_records(tmp_path, monkeypatch):
    real_run = ensemble.run
    def incomplete(bundle, output, **kwargs):
        report = real_run(bundle, output, **kwargs)
        if "fail" in Path(output).parts and kwargs["seed"] == 2:
            report.update(completed=False, stop_reason="provider_stopped_before_round_commit", error="fixture")
            write_json(Path(output)/"report.json", report)
        return report
    monkeypatch.setattr(ensemble, "run", incomplete)
    report = ensemble.run_ensemble({"control": fixture_bundle(), "fail": fixture_bundle()}, tmp_path/"runs",
                                  steps=1, seeds=[1, 2], bootstrap_replicates=100)
    assert not report["completed"] and report["excluded_arms"] == ["fail"]
    assert len(report["failures"]) == 1 and len(report["runs"]) == 4
    assert {row["arm"] for row in read_csv(tmp_path/"runs/summary.csv")} == {"control"}
    assert read_csv(tmp_path/"runs/paired_differences.csv") == []


def test_unrelated_inputs_duplicate_seeds_and_invalid_names_fail_before_output(tmp_path):
    other = fixture_bundle()
    other["cells"][0]["expression"]["ACTB"] = 40
    other = prepare(other["cells"], other["context"], assayed_genes=["ACTB"])
    with pytest.raises(ValueError, match="unrelated measured input"):
        ensemble.run_ensemble({"a": fixture_bundle(), "b": other}, tmp_path/"wrong", steps=1)
    with pytest.raises(ValueError, match="unique integer seeds"):
        ensemble.run_ensemble({"a": fixture_bundle()}, tmp_path/"repeat", steps=1, seeds=[1, 1])
    with pytest.raises(ValueError, match="Arm names"):
        ensemble.run_ensemble({"../escape": fixture_bundle()}, tmp_path/"escape", steps=1)
    assert not list(tmp_path.iterdir())


def test_provider_factory_is_per_run_and_never_used_for_local_baseline(tmp_path):
    class MockProvider:
        model = "mock"
        def choose(self, rows, context):
            return {row["cell_id"]: dict(probabilities=row["prior"]) for row in rows}
        def usage(self):
            return dict(http_attempts=0)
    calls = []
    def factory(**kwargs):
        calls.append(kwargs)
        return MockProvider()
    report = ensemble.run_ensemble({"local": fixture_bundle(), "mock": fixture_bundle()}, tmp_path/"ensemble",
        steps=1, seeds=[1, 2], provider_factory=factory, arm_options={"mock": dict(mode="jev", jev_weight=1)},
        bootstrap_replicates=100)
    assert report["completed"] and len(calls) == 2
    assert {call["arm"] for call in calls} == {"mock"}
    assert {call["seed"] for call in calls} == {1, 2}


def test_modified_readout_contract_is_rejected(tmp_path):
    spec, _ = load_readout_spec()
    spec["replication_unit"] = "cells"
    write_json(tmp_path/"changed.json", spec)
    with pytest.raises(ValueError, match="differs"):
        load_readout_spec(tmp_path/"changed.json")


def test_lattice_readouts_use_occupancy_not_cell_units(tmp_path):
    initial = {"a": dict(lineage="malignant", status="alive", occupancy=.5)}
    final = {"a": dict(lineage="malignant", status="alive", occupancy=1.)}
    frames = [dict(step=0, time_hours=0, cells=initial, events=[]),
              dict(step=1, time_hours=1, cells=final, events=[dict(kind="occupancy_gain", cell_id="a",
                        target="a", lineage="malignant", added_occupancy=.5)])]
    write_json(tmp_path/"frames.json", frames)
    write_json(tmp_path/"report.json", dict(completed_steps=1))
    rows = extract_readouts(tmp_path, context=dict(kernel="lattice", interaction_mode="spatial"))
    index = {(r["step"], r["readout"], r["lineage"]): r for r in rows}
    assert index[(1, "occupied_bins", "malignant")]["value"] == 1
    assert index[(1, "total_occupancy", "malignant")]["value"] == 1
    assert index[(1, "cumulative_occupancy_gained", "malignant")]["value"] == .5
    assert not {r["readout"] for r in rows} & {"living_cells", "cumulative_divisions", "activated_fraction"}


def test_ensemble_rejects_cell_nucleus_unit_mismatch_even_with_same_input(tmp_path):
    from test_v3_contracts import context, cells
    source = cells()
    a = prepare(source, context(), measurement={"unit": "segmented_cell"})
    b = prepare(source, context(), measurement={"unit": "nucleus"})
    assert a["guidebook"]["provenance"]["input_sha256"] == b["guidebook"]["provenance"]["input_sha256"]
    with pytest.raises(ValueError, match="measurement semantics"):
        ensemble.run_ensemble({"cells": a, "nuclei": b}, tmp_path/"mismatched", steps=1)
    assert not (tmp_path/"mismatched").exists()


def test_ensemble_rejects_lattice_capacity_scale_mismatch(tmp_path):
    from test_v3_lattice import prepared
    a, b = prepared(lattice=dict(capacity=1)), prepared(lattice=dict(capacity=2))
    with pytest.raises(ValueError, match="same occupancy capacity"):
        ensemble.run_ensemble({"a": a, "b": b}, tmp_path/"different_units", steps=1)
    assert not (tmp_path/"different_units").exists()
