"""Frozen, explicitly denominated readouts from committed simulation snapshots.

Missing values stay null. Spatial measurements never cross samples and are not
defined for well-mixed runs. Event counts use committed events, not proposals.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path


SPEC_PATH = Path(__file__).with_name("readouts_v3.json")
EVENTS = {"kill": "cumulative_kills", "apoptosis": "cumulative_apoptoses",
          "clear_corpse": "cumulative_clearances", "division": "cumulative_divisions"}


def load_readout_spec(path=None):
    """The public .yaml is JSON-compatible YAML; no YAML dependency is needed."""
    source = Path(path) if path is not None else SPEC_PATH
    raw = source.read_bytes()
    spec = json.loads(raw)
    if spec.get("schema_version") != "3.0" or spec.get("readout_version") != "1":
        raise ValueError("Unsupported readout specification version.")
    # A changed prose definition cannot silently drive the existing extractor.
    if spec != json.loads(SPEC_PATH.read_text()):
        raise ValueError("Readout specification differs from the implemented version; version the extractor with it.")
    return spec, hashlib.sha256(raw).hexdigest()


def iter_frames(directory):
    """Read v3 streams without retaining historical cell populations in memory."""
    directory = Path(directory)
    stream = directory / "frames.jsonl"
    if stream.exists():
        with stream.open(encoding="utf8") as handle:
            for line in handle:
                if line.strip():
                    yield json.loads(line)
    else:
        yield from json.loads((directory / "frames.json").read_text())


def _js_by_step(directory):
    result = defaultdict(lambda: [0.0, 0])
    path = Path(directory) / "decisions.jsonl"
    if not path.exists():
        return result
    with path.open(encoding="utf8") as handle:
        for line in handle:
            if not line.strip():
                continue
            decision = json.loads(line)
            value = decision.get("jev_prior_js_bits")
            if value is None:
                continue
            if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1 + 1e-12:
                raise ValueError("Invalid Jev-prior divergence in decision log.")
            lineage = decision["observation"]["lineage"]
            for key in ((decision["step"] + 1, lineage), (decision["step"] + 1, "ALL")):
                result[key][0] += value
                result[key][1] += 1
    return result


def spatial_readouts(cells, context):
    """Same-sample centroid distances, independent of the capped action graph."""
    if context.get("interaction_mode") != "spatial":
        return {"malignant_activated_cd8_contact_fraction": (None, 0, "nonspatial"),
                "mean_cd8_to_malignant_distance_um": (None, 0, "nonspatial")}
    living = [c for c in cells.values() if c["status"] == "alive"]
    if any("x" not in c or "y" not in c for c in living):
        raise ValueError("Spatial readouts require measured coordinates for every living cell.")
    malignant = [c for c in living if c["lineage"] == "malignant"]
    cd8 = [c for c in living if c["lineage"] == "T_CD8"]
    radius2 = context["neighbor_radius_um"] ** 2
    def distance2(a, b):
        return (a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2
    contacts = sum(any(t["sample_id"] == m["sample_id"] and t.get("activated", False)
                       and distance2(t, m) <= radius2 for t in cd8) for m in malignant)
    distances = []
    for t in cd8:
        candidates = [distance2(t, m) for m in malignant if t["sample_id"] == m["sample_id"]]
        if candidates:
            distances.append(math.sqrt(min(candidates)))
    return {
        "malignant_activated_cd8_contact_fraction":
            (contacts / len(malignant), len(malignant), "ok") if malignant else (None, 0, "no_malignant_cells"),
        "mean_cd8_to_malignant_distance_um":
            (sum(distances) / len(distances), len(distances), "ok") if distances else (None, 0, "no_eligible_cd8_target_pairs"),
    }


def extract_readouts(directory, *, context=None, readouts_path=None):
    """Return tidy records, one value per time × readout × lineage.

Event lineage denotes the actor (including clearance and kill); killed-target
lineages are additionally available as cumulative_killed_targets. JS at a
snapshot describes the decision round that produced that snapshot.
"""
    load_readout_spec(readouts_path)
    directory = Path(directory)
    if context is None:
        context = json.loads((directory / "guidebook.json").read_text())["context"]
    lattice = context.get("kernel") == "lattice"
    report = json.loads((directory / "report.json").read_text())
    js = _js_by_step(directory)
    counts, identity = Counter(), {}
    records, last_step, last_time, lineages = [], -1, -1.0, set()
    for frame in iter_frames(directory):
        step, time = frame["step"], frame["time_hours"]
        if type(step) is not int or step != last_step + 1 or not math.isfinite(time) or time <= last_time:
            raise ValueError("Frames must start at step zero and have contiguous steps and increasing times.")
        if step == 0 and time != 0:
            raise ValueError("Initial frame must be at zero hours.")
        cells = frame["cells"]
        for cid, cell in cells.items():
            identity[cid] = cell["lineage"]
            lineages.add(cell["lineage"])
        for event in frame["events"]:
            if lattice and event["kind"] in ("occupancy_gain", "occupancy_loss"):
                gain = event["kind"] == "occupancy_gain"
                metric = "cumulative_occupancy_gained" if gain else "cumulative_occupancy_lost"
                amount = event["added_occupancy" if gain else "removed_occupancy"]
                lineage = event["lineage"]
                if type(amount) not in (int, float) or not math.isfinite(amount) or amount < 0:
                    raise ValueError("Occupancy events require a finite nonnegative amount.")
                counts[(metric, lineage)] += amount
                counts[(metric, "ALL")] += amount
            if event["kind"] in EVENTS:
                actor = identity.get(event["cell_id"])
                if actor is None:
                    raise ValueError("Committed event actor cannot be linked to a lineage.")
                metric = EVENTS[event["kind"]]
                counts[(metric, actor)] += 1
                counts[(metric, "ALL")] += 1
            if event["kind"] == "kill":
                target = identity.get(event["target"])
                if target is None:
                    raise ValueError("Killed target cannot be linked to a lineage.")
                counts[("cumulative_killed_targets", target)] += 1
                counts[("cumulative_killed_targets", "ALL")] += 1
        def emit(metric, lineage, value, n, status="ok"):
            records.append(dict(step=step, time_hours=float(time), readout=metric,
                                lineage=lineage, value=value, n=n, status=status))
        for lineage in ["ALL", *sorted(lineages)]:
            selected = [c for c in cells.values() if lineage == "ALL" or c["lineage"] == lineage]
            living = [c for c in selected if c["status"] == "alive"]
            if lattice:
                occupied = [c for c in living if c.get("occupancy", 1.0) > 0]
                emit("occupied_bins", lineage, len(occupied), len(selected))
                emit("total_occupancy", lineage, sum(c.get("occupancy", 1.0) for c in occupied), len(occupied))
                for metric in ("cumulative_occupancy_gained", "cumulative_occupancy_lost"):
                    emit(metric, lineage, counts[(metric, lineage)], counts[(metric, lineage)])
                total, n = js[(step, lineage)]
                emit("jev_prior_js_bits", lineage, total / n if n else None, n,
                     "ok" if n else "no_provider_decisions")
                continue
            emit("living_cells", lineage, len(living), len(selected))
            for metric, state in (("activated_fraction", "activated"), ("exhausted_fraction", "exhausted")):
                value = sum(bool(c.get(state, False)) for c in living) / len(living) if living else None
                emit(metric, lineage, value, len(living), "ok" if living else "no_living_cells")
            for metric in (*EVENTS.values(), "cumulative_killed_targets"):
                emit(metric, lineage, counts[(metric, lineage)], counts[(metric, lineage)])
            total, n = js[(step, lineage)]
            emit("jev_prior_js_bits", lineage, total / n if n else None, n,
                 "ok" if n else "no_provider_decisions")
        if not lattice:
            for metric, (value, n, status) in spatial_readouts(cells, context).items():
                emit(metric, "ALL", value, n, status)
        last_step, last_time = step, time
    if last_step != report["completed_steps"]:
        raise ValueError("Frame stream and report disagree about completed steps.")
    return records
