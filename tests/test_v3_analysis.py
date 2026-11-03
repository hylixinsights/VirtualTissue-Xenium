"""Analytical fixtures for sensitivity and explicit-window spatial estimates."""
from copy import deepcopy
import math

import pytest

from virtualtissue_xenium.v2.sensitivity import analyze_sobol, make_design, run_sensitivity
from virtualtissue_xenium.v2.spatial_stats import cluster_sizes, cross_l, held_out_envelope
from test_v3_ensemble import fixture_bundle


def test_pick_freeze_indices_recover_known_additive_variance_and_shift_invariance():
    design = make_design({"x": [0, 1], "y": [0, 1]}, n=2048, seed=77)
    values = {p["id"]: p["parameters"]["x"] + 2*p["parameters"]["y"] for p in design["points"]}
    result = analyze_sobol(design, values, bootstrap_replicates=100, seed=5)
    shifted = analyze_sobol(design, {k: v+1000 for k, v in values.items()}, bootstrap_replicates=100, seed=5)
    for name, expected in (("x", .2), ("y", .8)):
        assert result["indices"][name]["S1"] == pytest.approx(expected, abs=.06)
        assert result["indices"][name]["ST"] == pytest.approx(expected, abs=.06)
        assert result["indices"][name]["S1"] == pytest.approx(shifted["indices"][name]["S1"], abs=1e-10)
    assert design["n_evaluations"] == 2048*4


def test_lhs_strata_integer_distribution_and_no_invalid_sobol_label():
    design = make_design({"x": [0, 1], "n": {"low": 1, "high": 4, "distribution": "uniform_integer"}},
                         n=20, method="lhs", seed=4)
    assert {int(p["parameters"]["x"]*20) for p in design["points"]} == set(range(20))
    assert all(type(p["parameters"]["n"]) is int and 1 <= p["parameters"]["n"] <= 4 for p in design["points"])
    with pytest.raises(ValueError, match="not LHS"):
        analyze_sobol(design, {p["id"]: 1 for p in design["points"]}, bootstrap_replicates=100)


def test_constant_or_missing_sobol_response_is_not_misreported():
    design = make_design({"x": [0, 1]}, n=3)
    values = {p["id"]: 1 for p in design["points"]}
    result = analyze_sobol(design, values, bootstrap_replicates=100)
    assert result["indices"]["x"]["S1"] is None
    assert result["indices"]["x"]["status"] == "constant_response_undefined"
    values.pop(next(iter(values)))
    with pytest.raises(ValueError, match="every design point"):
        analyze_sobol(design, values, bootstrap_replicates=100)


def test_local_sensitivity_records_complete_parameter_grid_and_rejects_pooling(tmp_path):
    report = run_sensitivity(fixture_bundle(), {"action_priors.malignant.PROLIFERATE": [.001, .004]},
        tmp_path/"sensitivity", horizon_hours=.5, seeds=[1, 2], n=2, method="lhs", bootstrap_replicates=100)
    assert report["completed"] and report["n_parameter_points"] == 2
    assert not report["sobol_indices_written"]
    with pytest.raises(ValueError, match="Unsupported local sensitivity"):
        run_sensitivity(fixture_bundle(), {"jev_weight": [0, 1]}, tmp_path/"wrong",
                        horizon_hours=.5, n=2, method="lhs")
    assert not (tmp_path/"wrong").exists()


def spatial_cells():
    return [dict(cell_id="t", sample_id="a", lineage="T_CD8", x=1., y=1.),
            dict(cell_id="m", sample_id="a", lineage="malignant", x=4., y=5.),
            dict(cell_id="other", sample_id="b", lineage="malignant", x=1., y=1.)]


def test_cross_l_uses_rectangular_translation_weight_and_same_sample():
    result = cross_l(spatial_cells(), [0, 4, 5], window=[0, 0, 10, 10], sample_id="a")
    assert result["n_source"] == result["n_target"] == 1
    assert result["curve"][1]["K_um2"] == 0
    expected_k = 100**2/((10-3)*(10-4))
    assert result["curve"][2]["K_um2"] == pytest.approx(expected_k)
    assert result["curve"][2]["L_um"] == pytest.approx(math.sqrt(expected_k/math.pi))
    shifted = deepcopy(spatial_cells())
    for c in shifted:
        c["x"] += 1
        c["y"] += 1
    assert cross_l(shifted, [0, 4, 5], window=[0, 0, 10, 10], sample_id="a")["curve"] == result["curve"]


def test_cross_l_rejects_implicit_or_unstable_windows_and_keeps_missing_null():
    with pytest.raises(ValueError, match="outside"):
        cross_l(spatial_cells(), [0, 1], window=[0, 0, 3, 3], sample_id="a")
    with pytest.raises(ValueError, match="half"):
        cross_l(spatial_cells(), [0, 6], window=[0, 0, 10, 10], sample_id="a")
    missing = cross_l(spatial_cells(), [0, 5], window=[0, 0, 10, 10], sample_id="b")
    assert missing["status"] == "missing_lineage" and missing["curve"][0]["K_um2"] is None


def test_clusters_and_held_out_envelopes_keep_explicit_semantics():
    cells = spatial_cells()+[dict(cell_id="m2", sample_id="a", lineage="malignant", x=4.1, y=5.)]
    assert cluster_sizes(cells, radius_um=1, sample_id="a")["sizes"] == [2]
    observed = cross_l(cells, [0, 4], window=[0, 0, 10, 10], sample_id="a")
    comparison = held_out_envelope(observed, [observed, deepcopy(observed)])
    assert all(row["inside_pointwise_envelope"] for row in comparison["curve"])
    with pytest.raises(ValueError, match="two"):
        held_out_envelope(observed, [observed])
