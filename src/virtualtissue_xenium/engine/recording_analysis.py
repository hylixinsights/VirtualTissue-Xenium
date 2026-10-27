"""Read-only accounting shared by result tables and the replay exporter."""
from collections import Counter
from statistics import mean, median


def tumour_cycle_diagnostics(episode, until=None):
    """Decision-level evidence; use Jev's reported values before compatibility repairs."""
    until = episode["summary"]["time_min"] if until is None else until
    labels = {b["id"]: b["label"] for b in episode.get("bodies", [])}
    groups, probabilities, margins = {}, [], []
    selected = binary = near_half = repaired = 0
    for d in episode.get("decisions", []):
        observation = d.get("observation") or {}
        label = observation.get("primary_label") or labels.get(d["cell_id"], "")
        answer = d.get("answer") or {}
        p = answer.get("reported_probabilities") or answer.get("probabilities") or {}
        if (d["time_min"] > until or d.get("status") != "validated" or not label.startswith("Tumor")
                or "ENTER_CYCLE" not in p):
            continue
        value = p["ENTER_CYCLE"]
        probabilities.append(value)
        selected += answer.get("choice") == "ENTER_CYCLE"
        repaired += bool(answer.get("probability_repairs"))
        alternatives = [v for k, v in p.items() if k != "ENTER_CYCLE"]
        if alternatives:
            margins.append(value - max(alternatives))
        if set(p) == {"ENTER_CYCLE", "REMAIN_QUIESCENT"}:
            binary += 1
            near_half += abs(value - .5) <= .1 + 1e-9
        mitotic = (observation.get("t0_aucell_percentile") or {}).get("MITOTIC_CELL_CYCLE")
        quartile = "missing" if mitotic is None else f"Q{min(4, int(mitotic * 4) + 1)}"
        row = groups.setdefault(quartile, dict(values=[], enter_cycle_choices=0))
        row["values"].append(value)
        row["enter_cycle_choices"] += answer.get("choice") == "ENTER_CYCLE"
    return dict(decisions=len(probabilities), enter_cycle_choices=selected,
                probability_source="reported_probabilities when available; otherwise probabilities",
                probability_min=min(probabilities) if probabilities else None,
                probability_median=median(probabilities) if probabilities else None,
                probability_mean=mean(probabilities) if probabilities else None,
                probability_max=max(probabilities) if probabilities else None,
                histogram_0_1_steps=[sum(min(9, int(p * 10)) == i for p in probabilities) for i in range(10)],
                binary_decisions=binary, binary_within_0_1_of_half=near_half,
                binary_fraction_within_0_1_of_half=near_half / binary if binary else None,
                mean_enter_minus_best_alternative=mean(margins) if margins else None,
                repaired_decisions=repaired,
                by_mitotic_quartile={k: dict(decisions=len(v["values"]), probability_mean=mean(v["values"]),
                                              enter_cycle_choices=v["enter_cycle_choices"])
                                      for k, v in sorted(groups.items())})


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
