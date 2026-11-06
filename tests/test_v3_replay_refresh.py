"""Refreshing presentation never changes the saved scientific trajectory."""
import json
import re

from virtualtissue_xenium.v2.cli import main


def saved_run(directory, schema_version):
    directory.mkdir()
    frames = [dict(step=i, time_hours=i * .09,
                   cells={"a": dict(cell_id="a", lineage="malignant", x=i, y=0, status="alive")},
                   events=[]) for i in range(3)]
    report = dict(schema_version=schema_version, completed=True, completed_steps=2,
                  seed=4, mode="local", preflight=dict(interaction_mode="spatial"), events={}, usage={})
    (directory / "report.json").write_text(json.dumps(report))
    (directory / "frames.json").write_text(json.dumps(frames))
    if schema_version == "3.0":
        (directory / "frames.jsonl").write_text("".join(json.dumps(f) + "\n" for f in frames))
    for name in ("decisions.jsonl", "readouts.csv", "guidebook.json", "config.json"):
        (directory / name).write_text("saved scientific evidence\n")
    (directory / "replay.html").write_text("old replay")
    return frames


def evidence(directory):
    return {path.name: path.read_bytes() for path in directory.iterdir()
            if path.is_file() and path.name != "replay.html"}


def payload(path):
    return json.loads(re.search(r'<script id="run-data" type="application/json">(.*?)</script>',
                                path.read_text(), re.S).group(1))


def prohibit_simulation(*args, **kwargs):
    raise AssertionError("Replay refresh must not rerun the simulation.")


def test_v3_refresh_overwrites_only_presentation_and_retains_streaming(tmp_path, monkeypatch):
    directory = tmp_path / "run"
    frames = saved_run(directory, "3.0")
    before = evidence(directory)
    monkeypatch.setattr("virtualtissue_xenium.v2.cli.run", prohibit_simulation)
    # A duplicate JSON export may be absent: v3 refresh must consume JSONL.
    (directory / "frames.json").unlink()
    before.pop("frames.json")
    assert main(["replay", str(directory)]) == 0
    assert evidence(directory) == before
    data = payload(directory / "replay.html")
    assert len(data["frames"]) == len(frames)
    assert all("cells" not in stub for stub in data["frames"])
    assert data["replay"]["cache_chunks"] == 2
    assert (directory / "replay_frames/manifest.json").is_file()


def test_v2_refresh_can_write_a_separate_self_contained_copy(tmp_path, monkeypatch):
    directory = tmp_path / "run"
    frames = saved_run(directory, "2.0")
    before = evidence(directory)
    monkeypatch.setattr("virtualtissue_xenium.v2.cli.run", prohibit_simulation)
    output = tmp_path / "separate" / "replay.html"
    assert main(["replay", str(directory), "--output", str(output)]) == 0
    assert evidence(directory) == before
    assert (directory / "replay.html").read_text() == "old replay"
    assert payload(output)["frames"] == frames
    assert not output.with_name("replay_frames").exists()


def test_refresh_refuses_a_scientific_export_destination(tmp_path):
    directory = tmp_path / "run"
    saved_run(directory, "3.0")
    before = evidence(directory)
    assert main(["replay", str(directory), "--output", str(directory / "frames.jsonl")]) == 1
    assert evidence(directory) == before
    assert (directory / "replay.html").read_text() == "old replay"


def test_v3_refresh_requires_saved_jsonl_before_writing_any_presentation(tmp_path):
    directory = tmp_path / "run"
    saved_run(directory, "3.0")
    (directory / "frames.jsonl").unlink()
    assert main(["replay", str(directory)]) == 1
    assert (directory / "replay.html").read_text() == "old replay"
    assert not (directory / "replay_frames").exists()


def test_refresh_refuses_a_symbolic_link_to_saved_evidence(tmp_path):
    directory = tmp_path / "run"
    saved_run(directory, "3.0")
    before = evidence(directory)
    output = tmp_path / "linked.html"
    output.symlink_to(directory / "report.json")
    assert main(["replay", str(directory), "--output", str(output)]) == 1
    assert evidence(directory) == before


def test_refresh_refuses_linked_sidecars_before_writing(tmp_path):
    directory = tmp_path / "run"
    saved_run(directory, "3.0")
    before = evidence(directory)
    sidecars = directory / "replay_frames"
    sidecars.mkdir()
    (sidecars / "chunk_00000.js").symlink_to(directory / "report.json")
    assert main(["replay", str(directory)]) == 1
    assert evidence(directory) == before
    assert (directory / "replay.html").read_text() == "old replay"


def test_refresh_refuses_malformed_report_before_writing(tmp_path):
    directory = tmp_path / "run"
    saved_run(directory, "3.0")
    (directory / "report.json").write_text("[]")
    before = evidence(directory)
    assert main(["replay", str(directory)]) == 1
    assert evidence(directory) == before
    assert (directory / "replay.html").read_text() == "old replay"
