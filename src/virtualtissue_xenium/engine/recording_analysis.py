"""Read-only accounting shared by result tables and the replay exporter."""
from collections import Counter


def blocked_reasons(event, contract):
    if event.get("block_reasons"):
        return event["block_reasons"]
    kind = event["kind"]
    for cause in ("population_ceiling", "local_pressure", "compression"):
        if kind == "division_blocked_" + cause:
            return [cause]
    if kind == "division_blocked_crowding":
        pressure = event.get("local_pressure_sim")
        threshold = contract.get("division", {}).get("division_blocked_local_pressure")
        if pressure is not None and threshold is not None and pressure <= threshold:
            return ["population_ceiling"]  # the other side of the legacy OR was false
        return ["legacy_crowding_unspecified"]
    return []


def event_totals(episode, until=None):
    until = episode["summary"]["time_min"] if until is None else until
    events = [e for e in episode.get("events", []) if e["time_min"] <= until]
    counts = Counter(e["kind"] for e in events)
    deaths = {e["cell_id"] for e in events
              if e["kind"] in ("apoptosis_committed", "hypoxic_necrosis")}
    causes = Counter(c for e in events if e["kind"].startswith("division_blocked_")
                     for c in blocked_reasons(e, episode.get("contract", {})))
    return dict(divisions=counts["division"], apoptosis_committed=counts["apoptosis_committed"],
                hypoxic_necrosis=counts["hypoxic_necrosis"], unique_deaths=len(deaths),
                cytotoxic_pulses=counts["cytotoxic_pulse"], corpses_cleared=counts["clearance"],
                corpses_degraded=counts["corpse_degraded"],
                division_blocks_population_ceiling=causes["population_ceiling"],
                division_blocks_local_pressure=causes["local_pressure"],
                division_blocks_compression=causes["compression"],
                division_blocks_legacy_unspecified=causes["legacy_crowding_unspecified"],
                growth_blocks_compression=counts["growth_blocked_compression"])


def population_limit(episode, until=None):
    until = episode["summary"]["time_min"] if until is None else until
    maximum = episode.get("contract", {}).get("biology", {}).get("max_population")
    bodies = episode.get("bodies", [])
    initial = sum(b["origin"] not in ("born", "recruited") for b in bodies)
    births = sorted(b["born_at_min"] for b in bodies if b["origin"] in ("born", "recruited"))
    reached = None
    if maximum is not None:
        headroom = maximum - initial
        if headroom <= 0:
            reached = 0
        elif len(births) >= headroom:
            reached = births[headroom - 1]
        if reached is not None and reached > until:
            reached = None
    return dict(maximum=maximum, initial_bodies=initial,
                initial_headroom=None if maximum is None else maximum-initial,
                allocated_bodies=initial+sum(t <= until for t in births), reached_min=reached)
