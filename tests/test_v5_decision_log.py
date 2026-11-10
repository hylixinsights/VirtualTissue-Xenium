"""The by-reference decision log is a lossless storage change, never a scientific one."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from virtualtissue_xenium.v2.cli import main
from virtualtissue_xenium.v2.inputs import read_json, write_json
from virtualtissue_xenium.v2.simulation import DECISION_LOG_FORMAT, read_decisions, run

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """One synthetic 270-cell bundle run twice: by reference (default) and inline (v4 layout)."""
    root = tmp_path_factory.mktemp("decision_log")
    subprocess.run([sys.executable, str(ROOT / "examples/v3/make_demo.py"), "--output", str(root / "demo")],
                   check=True, capture_output=True)
    bundle = read_json(root / "demo/control.json")
    reports = dict(reference=run(bundle, root / "reference", steps=12, seed=3),
                   inline=run(bundle, root / "inline", steps=12, seed=3, decision_log="inline"))
    return root, reports


def test_expanded_log_is_byte_identical_to_the_inline_log(runs):
    root, reports = runs
    assert reports["reference"]["decision_log"]["format"] == DECISION_LOG_FORMAT
    assert reports["inline"]["decision_log"]["format"] == "inline"
    inline = (root / "inline/decisions.jsonl").read_text().splitlines()
    expanded = [json.dumps(d, allow_nan=False) for d in read_decisions(root / "reference")]
    assert expanded == inline  # every step, every cell, key order included
    assert reports["reference"]["events"].get("division", 0) > 0  # daughters are part of the comparison
    daughters = [json.loads(line) for line in inline if "/d" in json.loads(line)["cell_id"]]
    assert daughters and all(d["observation"]["t0_evidence"] for d in daughters)


def test_storage_choice_changes_no_scientific_output(runs):
    root, reports = runs
    for name in ("frames.jsonl", "config.json", "guidebook.json"):
        assert sha(root / "reference" / name) == sha(root / "inline" / name)
    assert reports["reference"]["config_hash"] == reports["inline"]["config_hash"]
    assert "decision_log" not in read_json(root / "reference/config.json")
    assert (root / "reference/decisions.jsonl").stat().st_size * 5 < (root / "inline/decisions.jsonl").stat().st_size


def test_stored_records_reference_the_run_guidebook_with_a_hash(runs):
    root, _ = runs
    raw = [json.loads(line) for line in (root / "reference/decisions.jsonl").read_text().splitlines()]
    guidebook = read_json(root / "reference/guidebook.json")["cells"]
    assert raw and all("t0_evidence" not in d["observation"] for d in raw)
    for decision in raw:
        reference = decision["observation"]["t0_evidence_ref"]
        assert reference["file"] == "guidebook.json" and reference["guidebook_id"] in guidebook
        assert len(reference["sha256"]) == 64
    assert list(read_decisions(root / "reference", expand=False)) == raw


def test_expanded_records_do_not_share_evidence_objects(runs):
    root, _ = runs
    records = [d for d in read_decisions(root / "reference") if d["cell_id"] == "malignant_000"][:2]
    records[0]["observation"]["t0_evidence"]["cell_cycle"]["z"] = "edited"
    assert records[1]["observation"]["t0_evidence"]["cell_cycle"]["z"] != "edited"


def _copy(root, tmp_path):
    for name in ("decisions.jsonl", "guidebook.json"):
        (tmp_path / name).write_bytes((root / "reference" / name).read_bytes())
    return tmp_path


def test_edited_guidebook_is_detected(runs, tmp_path):
    folder = _copy(runs[0], tmp_path)
    guidebook = read_json(folder / "guidebook.json")
    guidebook["cells"]["malignant_000"]["evidence"]["cell_cycle"]["z"] = 2.5
    write_json(folder / "guidebook.json", guidebook)
    with pytest.raises(ValueError, match="does not match the hash"):
        list(read_decisions(folder))


def test_missing_or_foreign_guidebook_is_detected(runs, tmp_path):
    folder = _copy(runs[0], tmp_path)
    guidebook = read_json(folder / "guidebook.json")
    del guidebook["cells"]["malignant_000"]
    write_json(folder / "guidebook.json", guidebook)
    with pytest.raises(ValueError, match="no entry"):
        list(read_decisions(folder))
    (folder / "guidebook.json").unlink()
    with pytest.raises(ValueError, match="requires guidebook.json"):
        list(read_decisions(folder))


@pytest.mark.parametrize("file", ["../guidebook.json", "/etc/hosts", "other.json", None])
def test_reference_cannot_name_another_file(runs, tmp_path, file):
    folder = _copy(runs[0], tmp_path)
    first = json.loads((folder / "decisions.jsonl").read_text().splitlines()[0])
    first["observation"]["t0_evidence_ref"]["file"] = file
    (folder / "decisions.jsonl").write_text(json.dumps(first) + "\n")
    with pytest.raises(ValueError, match="invalid t0 evidence reference"):
        list(read_decisions(folder))


def test_reader_passes_inline_and_schema2_logs_through(tmp_path):
    line = dict(step=0, cell_id="a", observation=dict(lineage="x", t0_evidence=dict(s=dict(z=0.0))))
    (tmp_path / "decisions.jsonl").write_text(json.dumps(line) + "\n")
    assert list(read_decisions(tmp_path)) == [line]


def test_cli_exports_the_inline_layout_without_touching_the_run(runs, tmp_path, capsys):
    root, _ = runs
    before = {p.name: sha(p) for p in (root / "reference").iterdir() if p.is_file()}
    output = tmp_path / "expanded.jsonl"
    assert main(["decisions", str(root / "reference"), "--output", str(output)]) == 0
    assert output.read_bytes() == (root / "inline/decisions.jsonl").read_bytes()
    assert {p.name: sha(p) for p in (root / "reference").iterdir() if p.is_file()} == before
    capsys.readouterr()
    assert main(["decisions", str(root / "reference"), "--output", str(output)]) == 1  # refuses to overwrite
    assert main(["decisions", str(root / "reference"), "--output", str(root / "reference/decisions.jsonl")]) == 1


def test_unknown_log_mode_is_rejected(runs, tmp_path):
    bundle = read_json(runs[0] / "demo/control.json")
    with pytest.raises(ValueError, match="decision_log"):
        run(bundle, tmp_path / "bad", steps=1, decision_log="compressed")
