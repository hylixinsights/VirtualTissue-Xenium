"""Collaborator CLI workflow, matched endpoints and overwrite protection."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from virtualtissue_xenium.v2.cli import compare, main, parser
from virtualtissue_xenium.v2.inputs import prepare, read_json, write_json
from virtualtissue_xenium.v2.simulation import run


def context():
    return dict(schema_version="2.0", tissue="brain", organism="human", disease="control",
                assay="snrna", sampling=dict(design="enriched", description="Nuclei subset",
                missing_populations=["uncaptured cells"]), perturbations=[], recognition=False,
                time_step_hours=2, action_priors={"neuron": {"ACTIVATE": 10, "APOPTOSIS": 0}})


def cells():
    return [dict(cell_id="n1", cell_type="reviewed_mature_neuron", lineage="neuron",
                 expression={"ACTB": 5}, sample_id="sample1")]


@pytest.fixture
def paired_runs(tmp_path):
    bundle = prepare(cells(), context(), assayed_genes=["ACTB"])
    paths = [tmp_path / "short", tmp_path / "long"]
    run(bundle, paths[0], steps=1, seed=4)
    run(bundle, paths[1], steps=3, seed=4)
    return paths


def snapshot(directory):
    return {p.relative_to(directory).as_posix(): p.read_bytes()
            for p in Path(directory).rglob("*") if p.is_file()}


def test_compare_uses_latest_common_completed_endpoint(paired_runs):
    result = compare(paired_runs)
    assert result["shared_step"] == 1
    assert all(r["shared_step"] == 1 and r["time_hours"] == 2 for r in result["runs"])
    for actual, path in zip(result["runs"], paired_runs):
        frame = read_json(path / "frames.json")[1]
        assert actual["living"] == sum(c["status"] == "alive" for c in frame["cells"].values())
        assert actual["activated"] == sum(c["activated"] for c in frame["cells"].values())
    assert "no causal treatment efficacy claim" in result["interpretation"]


@pytest.mark.parametrize("key,value", [
    ("seed", 5), ("mode", "jev"), ("model", "jev-9.0.0"),
    ("provider_name", "different_provider"), ("provider_prompt_version", "changed-prompt"),
    ("max_questions", 1),
])
def test_compare_rejects_mismatched_execution_settings(paired_runs, key, value):
    path = paired_runs[1] / "config.json"
    config = read_json(path)
    config[key] = value
    write_json(path, config)
    with pytest.raises(ValueError, match=f"Unmatched run setting: {key}"):
        compare(paired_runs)


@pytest.mark.parametrize("key,value", [
    ("assayed_genes", ["ACTB", "UNMEASURED_EXTRA_FEATURE"]),
    ("score_method", "different_scoring_method"),
    ("gene_sets_sha256", "different_gene_set_hash"),
])
def test_compare_rejects_panel_and_scoring_changes(paired_runs, key, value):
    path = paired_runs[1] / "guidebook.json"
    guidebook = read_json(path)
    guidebook["provenance"][key] = value
    write_json(path, guidebook)
    with pytest.raises(ValueError, match="scoring assumption|different gene sets"):
        compare(paired_runs)


def test_compare_rejects_different_sampling_and_measured_expression(tmp_path, paired_runs):
    altered_sampling = context()
    altered_sampling["sampling"]["design"] = "whole_tissue"
    bundle = prepare(cells(), altered_sampling, assayed_genes=["ACTB"])
    sampling_run = tmp_path / "different_sampling"
    run(bundle, sampling_run, steps=1, seed=4)
    with pytest.raises(ValueError, match="context setting: sampling"):
        compare([paired_runs[0], sampling_run])

    altered_cells = cells()
    altered_cells[0]["expression"]["ACTB"] = 6
    bundle = prepare(altered_cells, context(), assayed_genes=["ACTB"])
    measured_run = tmp_path / "different_measurement"
    run(bundle, measured_run, steps=1, seed=4)
    with pytest.raises(ValueError, match="different measured input"):
        compare([paired_runs[0], measured_run])


def test_prompt_only_perturbation_is_an_accepted_declared_contrast(tmp_path, paired_runs):
    treated = context()
    treated["perturbations"] = [dict(type="treatment", name="hypothetical inflammation modulation",
        targets=["declared target"], mechanism_assumptions=["No measured treatment response is available."])]
    bundle = prepare(cells(), treated, assayed_genes=["ACTB"])
    output = tmp_path / "perturbed"
    run(bundle, output, steps=1, seed=4)
    result = compare([paired_runs[0], output])
    assert result["runs"][0]["context"]["perturbations"] == []
    assert result["runs"][1]["context"]["perturbations"][0]["status"] == "hypothetical"
    # Local mode does not create a drug response merely from a changed name.
    assert result["runs"][0]["living"] == result["runs"][1]["living"]
    assert result["runs"][0]["activated"] == result["runs"][1]["activated"]


def test_parser_keeps_network_and_resume_options_explicit():
    p = parser()
    defaults = p.parse_args(["run", "prepared.json", "run"])
    assert defaults.mode == "local"
    assert defaults.resume is False
    args = p.parse_args(["run", "prepared.json", "run", "--mode", "jev", "--resume",
        "--model", "jev-1.14.0", "--workers", "2", "--max-questions", "1", "--max-http-attempts", "3"])
    assert (args.mode, args.resume, args.max_questions, args.max_http_attempts) == ("jev", True, 1, 3)
    xenium = p.parse_args(["prepare-xenium", "raw", "context.json", "prepared.json", "--annotations", "reviewed.csv"])
    assert xenium.annotations == "reviewed.csv"


def test_prepare_inspect_run_and_existing_destinations_preserve_outputs(tmp_path, capsys):
    input_path, context_path, prepared_path = tmp_path / "cells.json", tmp_path / "context.json", tmp_path / "prepared.json"
    write_json(input_path, dict(cells=cells(), assayed_genes=["ACTB"]))
    write_json(context_path, context())
    command = ["prepare", str(input_path), str(context_path), str(prepared_path)]
    assert main(command) == 0
    assert read_json(prepared_path)["schema_version"] == "2.0"
    capsys.readouterr()
    original_bundle = prepared_path.read_bytes()
    original_guidebook = prepared_path.with_suffix(".guidebook.json").read_bytes()
    assert main(command) == 1
    assert prepared_path.read_bytes() == original_bundle
    assert prepared_path.with_suffix(".guidebook.json").read_bytes() == original_guidebook
    assert "exists" in capsys.readouterr().err.lower()
    assert main(["inspect", str(prepared_path)]) == 0
    assert json.loads(capsys.readouterr().out)["interaction_mode"] == "well_mixed"

    output = tmp_path / "run"
    command = ["run", str(prepared_path), str(output), "--steps", "2", "--seed", "4"]
    assert main(command) == 0
    assert json.loads(capsys.readouterr().out)["completed"] is True
    original = snapshot(output)
    assert main(command) == 1
    assert snapshot(output) == original
    assert "already has a run" in capsys.readouterr().err


def test_prepare_checks_sidecar_before_creating_any_output(tmp_path, capsys):
    input_path, context_path, output = tmp_path / "cells.json", tmp_path / "context.json", tmp_path / "prepared.json"
    write_json(input_path, dict(cells=cells(), assayed_genes=["ACTB"]))
    write_json(context_path, context())
    sidecar = output.with_suffix(".guidebook.json")
    sidecar.write_bytes(b"existing collaborator file\n")
    assert main(["prepare", str(input_path), str(context_path), str(output)]) == 1
    assert not output.exists()
    assert sidecar.read_bytes() == b"existing collaborator file\n"
    assert "exists" in capsys.readouterr().err.lower()


def test_run_preserves_preexisting_artifacts_without_a_config(tmp_path, capsys):
    prepared_path = tmp_path / "prepared.json"
    write_json(prepared_path, prepare(cells(), context(), assayed_genes=["ACTB"]))
    output = tmp_path / "prior_output"
    output.mkdir()
    (output / "frames.json").write_bytes(b"preexisting collaborator artifact\n")
    original = snapshot(output)
    assert main(["run", str(prepared_path), str(output), "--steps", "1"]) == 1
    assert snapshot(output) == original
    assert capsys.readouterr().err


def test_compare_requires_two_runs_without_an_unhandled_traceback(paired_runs, capsys):
    assert main(["compare", str(paired_runs[0])]) == 1
    assert "at least two" in capsys.readouterr().err.lower()
