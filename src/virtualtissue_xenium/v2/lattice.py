"""Fixed-bin occupancy dynamics, explicitly separate from the cell-agent kernel.

Occupancy is an abstract amount per measured bin, not a count of cells. Only
existing bin IDs can be colonised, and unknown annotations are not empty bins.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import math

from . import dynamics

DEFAULTS = dict(capacity=1.0, division_increment=1.0, low_density_threshold=0.0,
                neighbor_radius_multiplier=1.01)
DESCRIPTIONS = {
    "REST": "Maintain current fixed-bin occupancy and abstract occupant state.",
    "PROLIFERATE": "Increase occupancy in one compatible neighboring fixed bin; no new cell or bin is created.",
    "APOPTOSIS": "Empty this fixed bin under the declared occupancy model; this is not a counted cell death.",
    "KILL": "Empty one recognized malignant neighboring bin under the declared occupancy model.",
    "ACTIVATE": "Activate the abstract state attached to this bin's occupancy.",
    "EXHAUST": "Enter the abstract exhaustion state attached to this bin's occupancy.",
    "RECOVER": "Recover the abstract exhaustion state attached to this bin's occupancy.",
    "SECRETE": "Emit an abstract local signal from this occupied bin; no molecular concentration is implied.",
}


def validate_lattice(raw=None):
    if raw is None:
        raw = {}
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError(f"lattice accepts only {sorted(DEFAULTS)}.")
    value = dict(DEFAULTS, **raw)
    for key, number in value.items():
        if type(number) not in (int, float) or not math.isfinite(number) or number < 0:
            raise ValueError(f"lattice.{key} must be finite and nonnegative.")
        if key != "low_density_threshold" and number == 0:
            raise ValueError(f"lattice.{key} must be positive.")
        value[key] = float(number)
    if value["division_increment"] > value["capacity"] or value["low_density_threshold"] > value["capacity"]:
        raise ValueError("Lattice division increment and low-density threshold cannot exceed bin capacity.")
    return value


def enabled(context):
    return context.get("schema_version") == "3.0" and context.get("kernel") == "lattice"


def initialize_bin(cell, measurement, context):
    if measurement.get("unit") != "binned_spot":
        raise ValueError("kernel=lattice requires measurement.unit=binned_spot.")
    size = measurement.get("unit_size_um")
    if type(size) not in (int, float) or not math.isfinite(size) or size <= 0:
        raise ValueError("Lattice bins require a positive measurement.unit_size_um.")
    if "x" not in cell or "y" not in cell:
        raise ValueError("Lattice bins require measured x/y in micrometres.")
    occupancy = cell.get("occupancy", 1.0)
    if type(occupancy) not in (int, float) or not math.isfinite(occupancy) or not 0 <= occupancy <= context["lattice"]["capacity"]:
        raise ValueError("Bin occupancy must be finite and lie between zero and lattice.capacity.")
    cell.update(occupancy=float(occupancy), bin_size_um=float(size), bin_lineage=cell["lineage"],
                bin_cell_type=cell["cell_type"], occupancy_source_id=cell["cell_id"],
                status="alive" if occupancy > 0 else "removed")
    if occupancy == 0:
        cell.update(activated=False, exhausted=False, signal=0.0, suppressive_signal=0.0, activated_until=None,
                    activated_since=None, exhausted_since=None)


def graph(cells, context, *, seed, step, digest):
    if context.get("interaction_mode") != "spatial":
        raise ValueError("Lattice dynamics require measured spatial bin coordinates, not well_mixed mode.")
    result = {cid: [] for cid in sorted(cells)}
    # Fixed lattice includes empty bins; max_neighbors is a cell-agent setting.
    grid = defaultdict(list)
    seen = set()
    for cid in sorted(cells):
        c = cells[cid]
        location = (c["sample_id"], c["x"], c["y"])
        if location in seen:
            raise ValueError("Lattice bins must have unique centroid coordinates within each sample.")
        seen.add(location)
        width = c["bin_size_um"] * context["lattice"]["neighbor_radius_multiplier"]
        key = (c["sample_id"], width, math.floor(c["x"] / width), math.floor(c["y"] / width))
        grid[key].append(cid)
    for cid in sorted(cells):
        c = cells[cid]
        width = c["bin_size_um"] * context["lattice"]["neighbor_radius_multiplier"]
        gx, gy = math.floor(c["x"] / width), math.floor(c["y"] / width)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for nid in grid.get((c["sample_id"], width, gx + dx, gy + dy), []):
                    if nid == cid:
                        continue
                    other = cells[nid]
                    if math.hypot(c["x"] - other["x"], c["y"] - other["y"]) <= width:
                        result[cid].append(nid)
        result[cid].sort(key=lambda n: ((c["x"] - cells[n]["x"]) ** 2 + (c["y"] - cells[n]["y"]) ** 2,
                                      digest([seed, "lattice_target", step, cid, n])))
    return result


def expansion_targets(cell, neighbors, cells, context):
    settings = context["lattice"]
    out = []
    for cid in neighbors:
        target = cells[cid]
        amount = target["occupancy"]
        compatible = amount == 0 or (target["lineage"] == cell["lineage"] and target["cell_type"] == cell["cell_type"])
        if compatible and amount <= settings["low_density_threshold"] and amount + settings["division_increment"] <= settings["capacity"]:
            out.append(cid)
    return out


def commit(cells, proposals, context, *, seed, step, digest, seeded_rng):
    next_cells = deepcopy(cells)
    events = []
    dt = context["time_step_hours"]
    now = (step + 1) * dt
    retain = math.exp(-math.log(2) * dt / context["dynamics"]["signal_half_life_hours"])
    for c in next_cells.values():
        c["signal"] = c["signal"] * retain if c["occupancy"] > 0 else 0.0
        c["suppressive_signal"] = c.get("suppressive_signal", 0.0) * retain if c["occupancy"] > 0 else 0.0
    order = sorted(proposals, key=lambda p: digest([seed, "commit", step, p["cell_id"]]))
    emptied = set()
    for p in order:
        cid, action = p["cell_id"], p["action"]
        if action not in {"APOPTOSIS", "KILL"}:
            continue
        targets = [cid] if action == "APOPTOSIS" else p.get("targets", [])
        target_id = next((t for t in targets if t not in emptied and next_cells[t]["occupancy"] > 0), None)
        if target_id is None:
            events.append(dict(kind="conflict_blocked", cell_id=cid, action=action))
            continue
        target = next_cells[target_id]
        events.append(dict(kind="occupancy_loss", cell_id=cid, target=target_id, lineage=target["lineage"],
                           removed_occupancy=target["occupancy"], cause="apoptosis" if action == "APOPTOSIS" else "kill"))
        target.update(occupancy=0.0, status="removed", signal=0.0, suppressive_signal=0.0, activated=False, exhausted=False,
                      activated_since=None, activated_until=None, exhausted_since=None)
        emptied.add(target_id)
    for p in order:
        cid, action = p["cell_id"], p["action"]
        c = next_cells[cid]
        if c["occupancy"] == 0:
            continue
        if action == "PROLIFERATE":
            # Targets come from the old snapshot; deaths never open a new target
            # within the same simultaneous round.
            candidates = expansion_targets(c, p.get("division_targets", []), next_cells, context)
            if not candidates:
                events.append(dict(kind="occupancy_capacity_blocked", cell_id=cid))
                continue
            target_id = candidates[0]
            target = next_cells[target_id]
            amount = context["lattice"]["division_increment"]
            if target["occupancy"] == 0:
                for key in ("lineage", "cell_type", "guidebook_id", "activated", "exhausted", "activated_since",
                            "activated_until", "exhausted_since"):
                    target[key] = deepcopy(c.get(key))
                target.update(generation=c["generation"] + 1, occupancy_source_id=cid, signal=0.0, suppressive_signal=0.0)
            target.update(occupancy=target["occupancy"] + amount, status="alive")
            events.append(dict(kind="occupancy_gain", cell_id=cid, target=target_id, lineage=c["lineage"], added_occupancy=amount))
        elif action == "ACTIVATE":
            c["activated"] = True
            duration = seeded_rng(seed, "activation_clock", step, cid).expovariate(math.log(2) / context["dynamics"]["activation_half_life_hours"])
            c.update(activated_since=now, activated_until=now + duration)
            events.append(dict(kind="activate", cell_id=cid, activated_until=c["activated_until"]))
        elif action == "EXHAUST":
            c.update(exhausted=True, exhausted_since=now)
            events.append(dict(kind="exhaust", cell_id=cid))
        elif action == "RECOVER":
            c.update(exhausted=False, exhausted_since=None)
            events.append(dict(kind="recover", cell_id=cid))
        elif action == "SECRETE":
            channel = "suppressive_signal" if dynamics.lineage_rule(c, context).get("secretion_kind") == "suppressive" else "signal"
            c[channel] = c.get(channel, 0.0) + 1.0
    dynamics.apply_clocks(next_cells, context, now, events)
    return next_cells, events


def summary(cells, context):
    sums = defaultdict(float)
    for c in cells.values():
        if c["occupancy"] > 0:
            sums[c["lineage"]] += c["occupancy"]
    return dict(total=sum(sums.values()), by_lineage=dict(sorted(sums.items())),
                unit="abstract occupancy units; the declared capacity represents a full measured bin",
                capacity_per_bin=context["lattice"]["capacity"],
                occupied_bins=sum(c["occupancy"] > 0 for c in cells.values()))
