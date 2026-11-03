"""Offline lazy replay preserves data while bounding the browser's frame cache."""
import json
from pathlib import Path
import re
import shutil
import subprocess

import pytest

from virtualtissue_xenium.v2.replay import LAZY_SCRIPT, write_replay, write_replay_stream


def payload(path):
    return json.loads(re.search(r'<script id="run-data" type="application/json">(.*?)</script>', path.read_text(), re.S).group(1))


def fixtures(count):
    return [dict(step=i, time_hours=i*.5,
                 cells={"cell": dict(cell_id="cell", lineage="malignant", x=i, y=-i, status="alive")}, events=[])
            for i in range(count)]


def test_chunks_are_bounded_lightweight_and_roundtrip_exactly(tmp_path):
    frames = fixtures(127)
    source = tmp_path / "frames.jsonl"
    source.write_text("".join(json.dumps(frame)+"\n" for frame in frames))
    write_replay_stream(tmp_path / "replay.html", source, dict(schema_version="3.0"))
    data = payload(tmp_path / "replay.html")
    assert len(data["frames"]) == 127 and all("cells" not in f for f in data["frames"])
    assert [c["count"] for c in data["replay"]["chunks"]] == [50, 50, 27]
    assert data["replay"]["cache_chunks"] == 2
    assert data["replay"]["bounds"] == [0, 126, -126, 0]
    assert data["replay"]["cell_order"] == ["cell"]
    reconstructed = []
    for chunk in data["replay"]["chunks"]:
        script = (tmp_path / chunk["src"]).read_text()
        match = re.fullmatch(r'window\.__vtxReplayChunk\((\d+),(\[.*\])\);\s*', script, re.S)
        assert match
        reconstructed.extend(json.loads(match.group(2)))
    assert reconstructed == frames
    assert json.loads((tmp_path / "replay_frames/manifest.json").read_text()) == data


def test_v2_still_uses_self_contained_legacy_replay(tmp_path):
    frames = fixtures(2)
    write_replay(tmp_path / "v2.html", frames, dict(schema_version="2.0"))
    data = payload(tmp_path / "v2.html")
    assert data["frames"] == frames and "replay" not in data
    assert not (tmp_path / "v2_frames").exists()


def test_untrusted_text_cannot_terminate_data_or_chunk_scripts(tmp_path):
    marker = '</script><script>globalThis.unwanted=true</script>'
    frame = dict(step=0, time_hours=0, cells={marker: dict(cell_id=marker, lineage=marker, x=0, y=0, status="alive")}, events=[])
    source = tmp_path / "frames.jsonl"
    source.write_text(json.dumps(frame)+"\n")
    write_replay_stream(tmp_path / "replay.html", source, dict(schema_version="3.0", note=marker))
    html = (tmp_path / "replay.html").read_text()
    script = (tmp_path / "replay_frames/chunk_00000.js").read_text()
    assert marker not in html and marker not in script
    assert html.count('<script') == 2
    assert payload(tmp_path / "replay.html")["report"]["note"] == marker
    assert r'\u003c/script>' in script


def test_lazy_player_javascript_parses(tmp_path):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is optional; CI syntax check runs when available.")
    script = tmp_path / "lazy-player.js"
    script.write_text(LAZY_SCRIPT)
    result = subprocess.run([node, "--check", str(script)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_invalid_chunk_size_or_empty_trajectory_is_refused(tmp_path):
    source = tmp_path / "empty.jsonl"
    source.write_text("")
    with pytest.raises(ValueError, match="chunk_size"):
        write_replay_stream(tmp_path / "replay.html", source, {}, chunk_size=0)
    with pytest.raises(ValueError, match="at least one"):
        write_replay_stream(tmp_path / "replay.html", source, {})


def test_actual_player_loader_seek_cache_play_and_selection_with_dom_double(tmp_path):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is optional; functional DOM-double test runs when available.")
    frames = fixtures(151)
    for frame in frames[1:]:
        i = frame["step"]
        frame["events"] = [dict(kind="move", cell_id="cell", from_xy=[i-1, 1-i], to_xy=[i, -i])]
    source = tmp_path / "frames.jsonl"
    source.write_text("".join(json.dumps(f)+"\n" for f in frames))
    report = dict(schema_version="3.0", kernel="cell_agent", preflight=dict(interaction_mode="spatial"),
                  mode="local", completed=True, completed_steps=150, seed=1, events={}, usage={})
    write_replay_stream(tmp_path / "replay.html", source, report)
    harness = Path(__file__).with_name("test_v3_replay_loader.cjs")
    checked = subprocess.run([node, str(harness), str(tmp_path)], capture_output=True, text=True, timeout=15)
    assert checked.returncode == 0, checked.stderr
    evidence = json.loads(checked.stdout)
    assert evidence["passed"] and evidence["cacheChunks"] <= 2
