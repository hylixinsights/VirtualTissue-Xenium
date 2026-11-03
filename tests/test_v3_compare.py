"""Historical ontology contrasts stay explicit and lattice/token units stay honest."""
from copy import deepcopy

import pytest

from virtualtissue_xenium.v2.cli import compare, parser
from virtualtissue_xenium.v2.inputs import prepare, read_json, write_json
from virtualtissue_xenium.v2.simulation import run


def bundle(*, lattice=False):
    context = dict(schema_version="3.0", tissue="skin", organism="human", disease="fixture",
        assay="visium_hd" if lattice else "xenium", recognition=False, perturbations=[],
        sampling=dict(design="subsampled", description="compare fixture", missing_populations=[]),
        time_step_hours=.1)
    cells = [dict(cell_id="a", cell_type="unresolved", lineage="unknown", sample_id="s",
                  x=0, y=0, expression={"ACTB": 100})]
    measurement = None
    if lattice:
        context.update(kernel="lattice", lattice=dict(capacity=2))
        cells[0]["occupancy"] = .75
        cells.append(dict(cells[0], cell_id="b", x=8, occupancy=0))
        measurement = dict(unit="binned_spot", unit_size_um=8, deconvolution_method="fixture")
    return prepare(cells, context, assayed_genes=["ACTB"], measurement=measurement)


@pytest.fixture
def paired(tmp_path):
    prepared = bundle()
    paths = [tmp_path / "a", tmp_path / "b"]
    for path, steps in zip(paths, [1, 2]):
        run(prepared, path, steps=steps, seed=4)
    return paths


@pytest.mark.parametrize("key", ["lineages_sha256", "kinetics_sha256", "dataset_registry_sha256", "data_version_sha256"])
def test_asset_hash_changes_require_override_and_are_reported_exactly(paired, key):
    path = paired[1] / "guidebook.json"
    book = read_json(path)
    old = book["provenance"]["ontology_hashes"][key]
    book["provenance"]["ontology_hashes"][key] = "historical-other-hash"
    write_json(path, book)
    before = path.read_bytes()
    with pytest.raises(ValueError, match="allow-ontology-diff"):
        compare(paired)
    result = compare(paired, allow_ontology_diff=True)
    diff = result["ontology_differences"][0]
    assert diff["changed_asset_hashes"] == [dict(asset=key, reference=old, candidate="historical-other-hash")]
    assert diff["lineage_rows_available"] is True
    assert diff["recorded_lineage_row_differences"] == []
    assert "not stored" in diff["unavailable_rows"]
    assert "not a matched experiment" in result["interpretation"]
    assert path.read_bytes() == before


def test_override_reports_frozen_lineage_rows_and_kinetics_metadata(paired):
    path = paired[1] / "guidebook.json"
    book = read_json(path)
    before = deepcopy(book["provenance"]["lineage_rules"]["unknown"])
    book["provenance"]["lineage_rules"]["unknown"]["mean_radius_um"] = 9
    book["provenance"]["kinetics_sha256"] = "historic-table"
    book["provenance"]["ontology_version"] = "historical"
    book["context"]["ontology_version"] = "historical"
    write_json(path, book)
    result = compare(paired, allow_ontology_diff=True)
    diff = result["ontology_differences"][0]
    assert diff["recorded_lineage_row_differences"] == [dict(lineage="unknown", reference=before,
        candidate=book["provenance"]["lineage_rules"]["unknown"])]
    assert {item["field"] for item in diff["changed_metadata"]} >= {"kinetics_sha256", "ontology_version", "context.ontology_version"}


def test_missing_historical_rows_are_not_reconstructed(paired):
    for folder in paired:
        book = read_json(folder / "guidebook.json")
        book["provenance"].pop("lineage_rules")
        if folder == paired[1]:
            book["provenance"]["ontology_hashes"]["lineages_sha256"] = "old"
        write_json(folder / "guidebook.json", book)
    diff = compare(paired, allow_ontology_diff=True)["ontology_differences"][0]
    assert diff["lineage_rows_available"] is False
    assert diff["recorded_lineage_row_differences"] == []


def test_forensic_override_does_not_allow_unmatched_execution_or_measurements(paired):
    path = paired[1] / "config.json"
    config = read_json(path)
    config["seed"] = 7
    write_json(path, config)
    with pytest.raises(ValueError, match="Unmatched run setting: seed"):
        compare(paired, allow_ontology_diff=True)


def test_lattice_compare_uses_raw_occupancy_at_common_time_and_streams(tmp_path):
    prepared = bundle(lattice=True)
    paths = [tmp_path / "a", tmp_path / "b"]
    for folder, steps in zip(paths, [1, 2]):
        run(prepared, folder, steps=steps, seed=4)
        (folder / "frames.json").unlink()  # JSONL is sufficient; no full history load.
    result = compare(paths)
    assert result["shared_step"] == 1 and result["time_hours"] == .1
    for row in result["runs"]:
        assert not ({"living", "dead", "activated", "exhausted"} & row.keys())
        assert row["occupied_bins"] == 1 and row["empty_bins"] == 1
        assert row["occupancy"]["total"] == .75
        assert row["occupancy"]["capacity_per_bin"] == 2
        assert row["cumulative_occupancy_gained"] == 0
        assert row["cumulative_occupancy_lost"] == 0


def test_actual_tokens_include_unknown_attempts_and_full_run_scope(paired):
    for folder, inp, out, unknown in zip(paired, [10, 20], [4, 6], [0, 1]):
        config = read_json(folder / "config.json")
        config["mode"] = "jev"
        write_json(folder / "config.json", config)
        report = read_json(folder / "report.json")
        report["usage"] = dict(http_attempts=2, cache_hits=1, input_tokens=inp, output_tokens=out,
                               attempts_with_unknown_usage=unknown)
        write_json(folder / "report.json", report)
    ledger = compare(paired)["cost_ledger"]
    assert ledger["known_totals"] == dict(http_attempts=4, cache_hits=2, input_tokens=30, output_tokens=10,
                                         attempts_with_unknown_usage=1)
    assert ledger["token_totals_complete"] is False
    assert "entire runs" in ledger["scope"]
    assert ledger["currency_cost"] is None
    report = read_json(paired[1] / "report.json")
    report.pop("usage")
    write_json(paired[1] / "report.json", report)
    ledger = compare(paired)["cost_ledger"]
    assert ledger["known_totals"]["input_tokens"] == 10
    assert ledger["runs_missing_fields"]["input_tokens"] == 1
    assert ledger["by_run"][1]["input_tokens"] is None


def test_duplicate_paths_and_missing_common_frame_fail_closed(paired):
    with pytest.raises(ValueError, match="distinct run directories"):
        compare([paired[0], paired[0]])
    (paired[1] / "frames.jsonl").write_text("")
    with pytest.raises(ValueError, match="missing the common completed frame"):
        compare(paired)


def test_parser_forensic_override_is_explicit():
    assert parser().parse_args(["compare", "a", "b"]).allow_ontology_diff is False
    assert parser().parse_args(["compare", "a", "b", "--allow-ontology-diff"]).allow_ontology_diff is True
