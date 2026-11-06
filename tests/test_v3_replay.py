"""Offline lazy replay preserves data while bounding the browser's frame cache."""
from copy import deepcopy
import json
from pathlib import Path
import re
import shutil
import subprocess

import pytest

from virtualtissue_xenium import __version__
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


def test_sidebar_history_counts_recorded_events_and_populations_across_chunks(tmp_path):
    frames = fixtures(103)
    for frame in frames:
        step = frame["step"]
        if step >= 50:
            frame["cells"]["daughter"] = dict(cell_id="daughter", lineage="malignant", x=0, y=0,
                                               status="alive" if step < 100 else "dead")
        if step >= 101:
            frame["cells"]["daughter"]["status"] = "removed"
        if step == 50:
            frame["events"] = [dict(kind="division", cell_id="cell", daughter="daughter")]
        if step == 100:
            frame["events"] = [dict(kind="kill", cell_id="cd8", target="daughter")]
        if step == 101:
            frame["events"] = [dict(kind="clear_corpse", cell_id="macrophage", target="daughter")]
    source = tmp_path / "frames.jsonl"
    source.write_text("".join(json.dumps(f) + "\n" for f in frames))
    write_replay_stream(tmp_path / "replay.html", source, dict(schema_version="3.0"))
    stubs = payload(tmp_path / "replay.html")["frames"]
    assert all("cells" not in frame for frame in stubs)
    assert [stubs[i]["metrics"]["malignant"] for i in [0, 49, 50, 99, 100, 101]] == [1, 1, 2, 2, 1, 1]
    assert [stubs[i]["metrics"]["dead"] for i in [99, 100, 101]] == [0, 1, 0]
    assert stubs[101]["metrics"]["removed"] == 1
    assert stubs[49]["metrics"]["events"] == {}
    assert stubs[50]["metrics"]["events"] == {"division": 1}
    assert stubs[102]["metrics"]["events"] == {"division": 1, "kill": 1, "clear_corpse": 1}


def test_v2_still_uses_self_contained_legacy_replay(tmp_path):
    frames = fixtures(2)
    write_replay(tmp_path / "v2.html", frames, dict(schema_version="2.0"))
    data = payload(tmp_path / "v2.html")
    assert data["frames"] == frames and "replay" not in data
    assert not (tmp_path / "v2_frames").exists()


def test_embedded_and_streamed_players_restore_the_same_dark_dashboard(tmp_path):
    frames = fixtures(2)
    source = tmp_path / "frames.jsonl"
    source.write_text("".join(json.dumps(frame) + "\n" for frame in frames))
    for writer, argument, name in [(write_replay, frames, "v2.html"),
                                   (write_replay_stream, source, "v3.html")]:
        target = tmp_path / name
        writer(target, argument, dict(schema_version="3.0" if name == "v3.html" else "2.0"))
        html = target.read_text()
        assert payload(target)["product_version"] == __version__
        assert "--bg:#0a0e1a" in html.replace(" ", "")
        for element in ["tTumor", "tLiving", "tDiv", "tDeath", "chart", "inspector", "feed", "slider"]:
            assert f'id="{element}"' in html
        assert html.count("<script") == 2
        node = shutil.which("node")
        if node:
            harness = Path(__file__).with_name("test_v3_replay_loader.cjs")
            checked = subprocess.run([node, str(harness), str(tmp_path), "version", name],
                                     capture_output=True, text=True, timeout=15)
            assert checked.returncode == 0, checked.stderr
            actual = json.loads(checked.stdout)
            assert actual["productVersion"] == __version__
            assert actual["schemaVersion"] == ("3.0" if name == "v3.html" else "2.0")
            assert actual["title"] == f"VirtualTissue v{__version__} · research replay"


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
        if i >= 50:
            frame["cells"]["daughter"] = dict(cell_id="daughter", lineage="malignant", x=i+10, y=-i,
                                               status="alive" if i < 110 else "dead")
        if i == 50:
            frame["events"].append(dict(kind="division", cell_id="cell", daughter="daughter"))
        if i == 110:
            frame["events"].append(dict(kind="kill", cell_id="cd8", target="daughter"))
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


def test_default_cell_bodies_preserve_organic_shapes_nuclei_radii_and_supplied_outlines(tmp_path):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is optional for renderer DOM tests.")
    from virtualtissue_xenium.v3.geometry import normalize_geometry

    def cell(cid, x, y, lineage="malignant", radius=None, polygon=None):
        value = dict(cell_id=cid, cell_type=lineage, lineage=lineage, x=x, y=y,
                     status="alive", expression={}, sample_id="fixture")
        if radius is not None:
            value["radius_um"] = radius
        if polygon is not None:
            value["polygon"] = polygon
        return value

    cells = [cell("cell", 0, 0, radius=8), cell("cd8", 30, 0, "T_CD8", radius=5),
             cell("pair_a", 0, 30, radius=8), cell("pair_b", 10, 30, radius=8),
             cell("provided", 40, 30, polygon=[[36, 26], [44, 26], [44, 34], [36, 34]]),
             cell("default_radius", 60, 0, "unknown")]
    normalized, _ = normalize_geometry(cells)
    frame = dict(step=0, time_hours=0, cells={c["cell_id"]: c for c in normalized}, events=[])
    source = tmp_path / "frames.jsonl"
    frames = []
    for i, time in enumerate([0, 1/6, 1/3, 1]):
        saved = deepcopy(frame)
        saved.update(step=i, time_hours=time)
        saved["cells"]["provided"].update(x=40+4*i, y=30+2*i)
        frames.append(saved)
    source.write_text("".join(json.dumps(saved) + "\n" for saved in frames))
    original = source.read_bytes()
    report = dict(schema_version="3.0", kernel="cell_agent", preflight=dict(interaction_mode="spatial"),
                  mode="local", completed=True, completed_steps=3, seed=1, events={}, usage={})
    write_replay_stream(tmp_path / "replay.html", source, report)
    harness = Path(__file__).with_name("test_v3_replay_loader.cjs")
    checked = subprocess.run([node, str(harness), str(tmp_path), "shapes"],
                             capture_output=True, text=True, timeout=15)
    assert checked.returncode == 0, checked.stderr
    assert json.loads(checked.stdout)["passed"]
    assert source.read_bytes() == original
