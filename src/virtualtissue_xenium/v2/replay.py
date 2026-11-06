"""Original dark replay presentation, with bounded offline v3 frame loading.

The presentation reads committed snapshots only. It never changes simulation
rates, decisions, coordinates, geometry evidence, or scientific exports.
"""
from collections import Counter
from pathlib import Path
import json
import hashlib


WEB = Path(__file__).with_name("web")
LAZY_SCRIPT = (WEB / "replay_dark.js").read_text(encoding="utf8")
TEMPLATE = (WEB / "replay_dark.html").read_text(encoding="utf8").replace("__SCRIPT__", LAZY_SCRIPT)
LAZY_TEMPLATE = TEMPLATE


def _display_radii(lineages):
    """Freeze the declared radius fallback for the presentation, not the model."""
    from ..v3.ontology import lineage_rules
    rules = lineage_rules()
    return {name: rules[name]["mean_radius_um"] for name in sorted(lineages) if name in rules}


def _frame_metrics(frame, cumulative):
    """Tiny summaries make the sidebar/plot independent of the frame cache."""
    cells = frame["cells"].values()
    statuses = Counter(c.get("status", "alive") for c in cells)
    lineages = Counter(c["lineage"] for c in cells if c.get("status", "alive") == "alive")
    cumulative.update(e["kind"] for e in frame.get("events", []))
    return dict(living=statuses["alive"], dead=statuses["dead"], removed=statuses["removed"],
                malignant=lineages["malignant"],
                occupied=sum(c.get("occupancy", 0) > 0 for c in cells),
                empty=sum(c.get("occupancy") == 0 for c in cells),
                event_count=len(frame.get("events", [])), events=dict(cumulative), lineages=dict(lineages))


def write_replay(path, frames, report):
    cumulative = Counter()
    lineages = {c["lineage"] for frame in frames for c in frame["cells"].values()}
    value = dict(frames=frames, report=report, metrics=[_frame_metrics(f, cumulative) for f in frames],
                 radii_um=_display_radii(lineages))
    data = json.dumps(value, allow_nan=False).replace("<", "\\u003c")
    Path(path).write_text(TEMPLATE.replace("__DATA__", data), encoding="utf8")


def write_replay_stream(path, frames_jsonl, report, *, chunk_size=50):
    """Write bounded offline chunks; immutable t0 shapes are shared within a chunk.

    The player rehydrates geometry references without copying the shared shape.
    At most two chunks are retained. Full JSONL/JSON snapshots remain the audit
    exports; the sidecars are a lossless presentation representation.
    """
    from ..v3.geometry import projected_polygon
    if type(chunk_size) is not int or not 1 <= chunk_size <= 500:
        raise ValueError("Replay chunk_size must be an integer from 1 to 500.")
    path = Path(path)
    sidecars = path.with_name(path.stem + "_frames")
    sidecars.mkdir(parents=True, exist_ok=True)
    stubs, chunks, lineages, cell_order, seen_ids = [], [], set(), [], set()
    bounds = [None, None, None, None]
    cumulative = Counter()
    stream, shapes, has_polygons = None, {}, False
    def close_chunk():
        if stream is not None and not stream.closed:
            stream.write("]")
            if shapes:
                stream.write("," + json.dumps(shapes,allow_nan=False).replace("<","\\u003c"))
            stream.write(");\n")
            stream.close()
    try:
        with Path(frames_jsonl).open(encoding="utf8") as source:
            for raw in source:
                if not raw.strip():
                    continue
                frame = json.loads(raw)
                index = len(stubs)
                chunk_id, offset = divmod(index, chunk_size)
                if offset == 0:
                    close_chunk()
                    shapes = {}
                    filename = f"chunk_{chunk_id:05d}.js"
                    stream = (sidecars / filename).open("w", encoding="utf8")
                    stream.write(f"window.__vtxReplayChunk({chunk_id},[\n")
                    chunks.append(dict(src=sidecars.name + "/" + filename, start=index, count=0))
                displayed = dict(frame,cells={})
                for cid, cell in frame["cells"].items():
                    lineages.add(cell["lineage"])
                    if cid not in seen_ids:
                        seen_ids.add(cid)
                        cell_order.append(cid)
                    contour = projected_polygon(cell)
                    has_polygons = has_polygons or contour is not None
                    positions = (contour or []) + ([[cell["x"],cell["y"]]] if "x" in cell and "y" in cell else [])
                    for x,y in positions:
                        bounds[0] = x if bounds[0] is None else min(bounds[0], x)
                        bounds[1] = x if bounds[1] is None else max(bounds[1], x)
                        bounds[2] = y if bounds[2] is None else min(bounds[2], y)
                        bounds[3] = y if bounds[3] is None else max(bounds[3], y)
                    value = dict(cell)
                    if "t0_geometry" in value:
                        geometry = value.pop("t0_geometry")
                        geometry_key = hashlib.sha256(json.dumps(geometry,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
                        shapes[geometry_key] = geometry
                        value["t0_geometry_ref"] = geometry_key
                    displayed["cells"][cid] = value
                if offset:
                    stream.write(",\n")
                stream.write(json.dumps(displayed, allow_nan=False).replace("<", "\\u003c"))
                chunks[-1]["count"] += 1
                stubs.append(dict(step=frame["step"], time_hours=frame["time_hours"], chunk=chunk_id, offset=offset,
                                  metrics=_frame_metrics(frame, cumulative)))
    finally:
        close_chunk()
    if not stubs:
        raise ValueError("Replay requires at least one committed frame.")
    metadata = dict(format="vtx-offline-chunks-2", chunk_size=chunk_size, cache_chunks=2,
                    chunks=chunks, lineages=sorted(lineages), cell_order=cell_order, has_polygons=has_polygons,
                    radii_um=_display_radii(lineages),
                    geometry_encoding="immutable t0 records shared once per chunk; runtime contours are rigid translations",
                    bounds=bounds if bounds[0] is not None else [0, 1, 0, 1])
    value = dict(frames=stubs, report=report, replay=metadata)
    encoded = json.dumps(value, allow_nan=False).replace("<", "\\u003c")
    path.write_text(LAZY_TEMPLATE.replace("__DATA__", encoded), encoding="utf8")
    (sidecars / "manifest.json").write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf8")
