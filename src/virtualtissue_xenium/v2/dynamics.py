"""Declared v3 state clocks and spatial rules; no changes to the v2 kernel.

Defaults are modelling assumptions, not literature-calibrated measurements.
MOVE starts a persistent episode; movement then integrates over elapsed hours,
so a per-hour MOVE hazard is not incorrectly multiplied by dt twice.
"""
from __future__ import annotations

import math

DEFAULTS = {
    "activation_half_life_hours": 48.0,
    "terminal_exhaustion_hours": 168.0,
    "corpse_lifetime_hours": 72.0,
    "signal_half_life_hours": 0.5,
    "movement_enabled": True,
    "speed_um_per_hour": 6.0,
    "persistence_hours": 1.0,
    "chemotaxis_gain": 0.0,
    "crowding_enabled": True,
    "default_radius_um": 5.0,
    "n_max": 12,
    "repulsion_iterations": 3,
    "max_horizon_hours": 8760.0,
}


def validate_dynamics(raw=None):
    if raw is None:
        raw = {}
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError(f"dynamics accepts only {sorted(DEFAULTS)}.")
    result = dict(DEFAULTS, **raw)
    for name in ("movement_enabled", "crowding_enabled"):
        if type(result[name]) is not bool:
            raise ValueError(f"dynamics.{name} must be boolean.")
    for name in ("n_max", "repulsion_iterations"):
        value = result[name]
        if type(value) is not int or not 1 <= value <= (10000 if name == "n_max" else 20):
            raise ValueError(f"dynamics.{name} must be a positive bounded integer.")
    for name in set(DEFAULTS) - {"movement_enabled", "crowding_enabled", "n_max", "repulsion_iterations"}:
        value = result[name]
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError(f"dynamics.{name} must be finite and nonnegative.")
        if name not in {"speed_um_per_hour", "chemotaxis_gain"} and value == 0:
            raise ValueError(f"dynamics.{name} must be positive.")
        result[name] = float(value)
    return result


def enabled(context):
    return context.get("schema_version") == "3.0"


def lineage_rule(cell, context):
    from ..v3.ontology import lineage_rule as get_lineage_rule
    return get_lineage_rule(cell["lineage"], context.get("ontology_version", "3.0"))


def radius(cell, context):
    return float(lineage_rule(cell, context).get("mean_radius_um", context["dynamics"]["default_radius_um"]))


def initialise_clocks(cell, context, rng):
    cell.update(activated_since=0.0 if cell["activated"] else None,
                activated_until=(rng.expovariate(math.log(2) / context["dynamics"]["activation_half_life_hours"])
                                 if cell["activated"] else None),
                exhausted_since=0.0 if cell["exhausted"] else None,
                dead_since=None, moving_until=0.0, heading_radians=None, suppressive_signal=0.0)


def density_factor(cell, cells, context):
    if not context["dynamics"]["crowding_enabled"] or "x" not in cell:
        return 1.0, 0
    cutoff2 = (2 * radius(cell, context)) ** 2
    n = sum(other["cell_id"] != cell["cell_id"] and other["status"] == "alive"
            and other["sample_id"] == cell["sample_id"] and "x" in other
            and (other["x"] - cell["x"]) ** 2 + (other["y"] - cell["y"]) ** 2 <= cutoff2
            for other in cells.values())
    return max(0.0, 1 - n / context["dynamics"]["n_max"]), n


def apply_clocks(cells, context, now, events):
    """Implicit transitions at the end of a committed interval."""
    dynamics = context["dynamics"]
    for cid in sorted(cells):
        c = cells[cid]
        if c["status"] != "alive":
            c["signal"] = 0.0
            c["suppressive_signal"] = 0.0
            if c["status"] == "dead" and c.get("dead_since") is not None and now - c["dead_since"] >= dynamics["corpse_lifetime_hours"]:
                c["status"] = "removed"
                events.append(dict(kind="corpse_degraded", cell_id=cid, time_hours=now))
            continue
        if c["activated"] and c.get("activated_until") is not None and c["activated_until"] <= now:
            c["activated"] = False
            c["activated_until"] = None
            c["activated_since"] = None
            events.append(dict(kind="deactivate", cell_id=cid, time_hours=now))


def move(cell, snapshot, context, *, start, end, selected, rng):
    """Return a displacement event, integrating the active part of a moving episode."""
    d = context["dynamics"]
    if not d["movement_enabled"] or "x" not in cell or not lineage_rule(cell, context)["motile"]:
        return None
    if selected:
        previous = cell.get("heading_radians")
        direction = rng.uniform(-math.pi, math.pi)
        if previous is not None and cell.get("moving_until", 0) > start:
            retain = math.exp(-(end - start) / d["persistence_hours"])
            direction = previous + (1 - retain) * ((direction - previous + math.pi) % math.tau - math.pi)
        gx = gy = 0.0
        if d["chemotaxis_gain"]:
            for n in snapshot.values():
                if n["status"] != "alive" or n["sample_id"] != cell["sample_id"] or n["cell_id"] == cell["cell_id"] or "x" not in n:
                    continue
                dx, dy = n["x"] - cell["x"], n["y"] - cell["y"]
                distance = math.hypot(dx, dy)
                if 0 < distance <= context["neighbor_radius_um"]:
                    gx += n["signal"] * dx / distance
                    gy += n["signal"] * dy / distance
            direction = math.atan2(math.sin(direction) + d["chemotaxis_gain"] * gy,
                                   math.cos(direction) + d["chemotaxis_gain"] * gx)
        cell["heading_radians"] = direction
        cell["moving_until"] = start + rng.expovariate(1 / d["persistence_hours"])
    duration = max(0.0, min(end, cell.get("moving_until", 0.0)) - start)
    if not duration or cell.get("heading_radians") is None:
        return None
    old = [cell["x"], cell["y"]]
    distance = d["speed_um_per_hour"] * duration
    cell["x"] += distance * math.cos(cell["heading_radians"])
    cell["y"] += distance * math.sin(cell["heading_radians"])
    return dict(kind="move", cell_id=cell["cell_id"], from_xy=old, to_xy=[cell["x"], cell["y"]],
                distance_um=distance, active_hours=duration, moving_until=cell["moving_until"])


def place_daughter(parent, daughter, cells, context, rng):
    """Contact-distance placement followed by three declared local soft-repulsion passes.

    Only the newborn moves during relaxation, preserving synchronous positions of
    existing actors. This is a local overlap heuristic, not a mechanics solver.
    """
    if "x" not in daughter:
        return []
    angle = rng.random() * math.tau
    r = radius(daughter, context)
    separation = r + radius(parent, context)
    daughter["x"] = parent["x"] + math.cos(angle) * separation
    daughter["y"] = parent["y"] + math.sin(angle) * separation
    for _ in range(context["dynamics"]["repulsion_iterations"]):
        dx_total = dy_total = 0.0
        for cid in sorted(cells):
            other = cells[cid]
            if other["status"] != "alive" or other["sample_id"] != daughter["sample_id"] or "x" not in other:
                continue
            dx, dy = daughter["x"] - other["x"], daughter["y"] - other["y"]
            distance = math.hypot(dx, dy)
            minimum = r + radius(other, context)
            if distance < minimum:
                if distance == 0:
                    direction = rng.random() * math.tau
                    ux, uy = math.cos(direction), math.sin(direction)
                else:
                    ux, uy = dx / distance, dy / distance
                dx_total += .5 * (minimum - distance) * ux
                dy_total += .5 * (minimum - distance) * uy
        daughter["x"] += dx_total
        daughter["y"] += dy_total
    return [daughter["x"], daughter["y"]]
