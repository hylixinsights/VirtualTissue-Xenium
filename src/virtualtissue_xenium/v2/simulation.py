"""Simultaneous, seeded state transitions on an explicitly declared interaction graph.

This small v2 kernel models discrete states and local signalling. It does not
reuse the calibrated-looking physical constants of the lung-specific engine.
There is no learned gene-expression evolution, mechanics, or pharmacokinetics.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
import json
import math
import os
from pathlib import Path
import random

from .guidebook import KineticsStepError, competing_risk_probabilities
from .inputs import digest, validate_prepared, write_json
from .provider import BudgetExhausted, ProviderError, batches
from . import dynamics, lattice
from .frame_store import FrameStore
from ..v3 import geometry


def seeded_rng(seed, *parts):
    return random.Random(int(digest([seed, *parts]), 16))


def interaction_graph(cells, context, *, seed=21, step=0):
    """No edges across samples; nonspatial neighborhoods are explicit assumptions."""
    if lattice.enabled(context):
        return lattice.graph(cells, context, seed=seed, step=step, digest=digest)
    if dynamics.enabled(context):
        return _v3_interaction_graph(cells, context, seed=seed, step=step)
    ids = sorted(cid for cid, c in cells.items() if c["status"] != "removed")
    graph = {cid: [] for cid in ids}
    mode = context.get("interaction_mode", "spatial" if context["assay"] == "xenium" else "well_mixed")
    limit = context.get("max_neighbors", 20)
    if type(limit) is not int or not 1 <= limit <= 100:
        raise ValueError("max_neighbors must be an integer from 1 to 100.")
    by_sample = {}
    for cid in ids:
        by_sample.setdefault(cells[cid]["sample_id"], []).append(cid)
    if mode == "well_mixed":
        for cid in ids:
            candidates = [x for x in by_sample[cells[cid]["sample_id"]] if x != cid]
            graph[cid] = sorted(seeded_rng(seed, "neighbors", step, cid).sample(candidates, min(limit, len(candidates))))
        return graph
    if mode != "spatial":
        raise ValueError("interaction_mode must be spatial or well_mixed.")
    radius = context.get("neighbor_radius_um", 30.0)
    if type(radius) not in (int, float) or not math.isfinite(radius) or radius <= 0:
        raise ValueError("neighbor_radius_um must be finite and positive.")
    if any("x" not in cells[cid] or "y" not in cells[cid] for cid in ids):
        raise ValueError("Spatial mode requires measured x/y for every cell. Choose well_mixed explicitly for sc/snRNA-seq.")
    # Grid hash avoids quadratic all-pairs distances and needs no scientific stack.
    grid = {}
    for cid in ids:
        c = cells[cid]
        key = (c["sample_id"], math.floor(c["x"] / radius), math.floor(c["y"] / radius))
        grid.setdefault(key, []).append(cid)
    for cid in ids:
        c = cells[cid]
        gx, gy = math.floor(c["x"] / radius), math.floor(c["y"] / radius)
        candidates = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for nid in grid.get((c["sample_id"], gx + dx, gy + dy), []):
                    if nid == cid:
                        continue
                    n = cells[nid]
                    d2 = (c["x"] - n["x"]) ** 2 + (c["y"] - n["y"]) ** 2
                    if d2 <= radius ** 2:
                        candidates.append((d2, nid))
        graph[cid] = [nid for _, nid in sorted(candidates)[:limit]]
    return graph


def _v3_interaction_graph(cells, context, *, seed, step):
    """Symmetric live graph with a degree bound; corpses have independent edges.

    Spatial edges are admitted in increasing distance order while both live
    endpoints have capacity. Well-mixed candidate edges are sampled, then
    admitted in seeded order. Corpse edges never consume the live-neighbor cap.
    """
    ids = sorted(cid for cid, c in cells.items() if c["status"] != "removed")
    graph = {cid: [] for cid in ids}
    limit = context.get("max_neighbors", 20)
    if type(limit) is not int or not 1 <= limit <= 100:
        raise ValueError("max_neighbors must be an integer from 1 to 100.")
    mode = context.get("interaction_mode", "spatial")
    if mode not in {"spatial", "well_mixed"}:
        raise ValueError("interaction_mode must be spatial or well_mixed.")
    radius = context.get("neighbor_radius_um", 30.0)
    by_sample = {}
    for cid in ids:
        by_sample.setdefault(cells[cid]["sample_id"], []).append(cid)
    edges = []
    corpse_edges = []
    distances = {}
    metric = context.get("geometry", {}).get("interaction_distance", "centroid")
    if mode == "spatial":
        if type(radius) not in (int, float) or not math.isfinite(radius) or radius <= 0:
            raise ValueError("neighbor_radius_um must be finite and positive.")
        if any("x" not in cells[cid] or "y" not in cells[cid] for cid in ids):
            raise ValueError("Spatial mode requires measured x/y for every cell.")
        if metric == "boundary":
            for gap, first, second, _ in geometry.boundary_candidates(cells, ids, radius):
                cid,nid = sorted((first,second))
                distances[(cid,nid)] = gap ** 2
                edge = (gap ** 2, digest([seed,"edge_tie",step,cid,nid]), cid, nid)
                if cells[cid]["status"] == cells[nid]["status"] == "alive":
                    edges.append(edge)
                elif (cells[cid]["status"] == "alive") != (cells[nid]["status"] == "alive"):
                    corpse_edges.append(edge)
        else:
            grid = {}
            for cid in ids:
                c = cells[cid]
                key = (c["sample_id"], math.floor(c["x"] / radius), math.floor(c["y"] / radius))
                grid.setdefault(key, []).append(cid)
            for cid in ids:
                c = cells[cid]
                gx, gy = math.floor(c["x"] / radius), math.floor(c["y"] / radius)
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for nid in grid.get((c["sample_id"], gx + dx, gy + dy), []):
                            if nid <= cid:
                                continue
                            n = cells[nid]
                            d2 = (c["x"] - n["x"]) ** 2 + (c["y"] - n["y"]) ** 2
                            if d2 <= radius ** 2:
                                edge = (d2, digest([seed, "edge_tie", step, cid, nid]), cid, nid)
                                if c["status"] == n["status"] == "alive":
                                    edges.append(edge)
                                elif (c["status"] == "alive") != (n["status"] == "alive"):
                                    corpse_edges.append(edge)
    else:
        seen = set()
        for sample, members in sorted(by_sample.items()):
            live = [cid for cid in members if cells[cid]["status"] == "alive"]
            dead = [cid for cid in members if cells[cid]["status"] == "dead"]
            for cid in live:
                candidates = [n for n in live if n != cid]
                for nid in seeded_rng(seed, "neighbors", step, cid).sample(candidates, min(limit, len(candidates))):
                    a, b = sorted((cid, nid))
                    if (a, b) not in seen:
                        seen.add((a, b))
                        edges.append((0, digest([seed, "edge", step, a, b]), a, b))
                for nid in seeded_rng(seed, "corpse_neighbors", step, cid).sample(dead, min(limit, len(dead))):
                    corpse_edges.append((0, "", cid, nid))
    for _, _, cid, nid in sorted(edges):
        if len(graph[cid]) < limit and len(graph[nid]) < limit:
            graph[cid].append(nid)
            graph[nid].append(cid)
    for _, _, cid, nid in corpse_edges:
        graph[cid].append(nid)
        graph[nid].append(cid)
    for cid in graph:
        if mode == "spatial":
            c = cells[cid]
            graph[cid].sort(key=lambda nid: (distances[tuple(sorted((cid,nid)))] if metric == "boundary" else
                                           (c["x"] - cells[nid]["x"]) ** 2 + (c["y"] - cells[nid]["y"]) ** 2,
                                           digest([seed, "target_tie", step, cid, nid])))
        else:
            graph[cid].sort()
    return graph


def initialize(bundle, *, seed=21):
    result = {}
    for c in bundle["cells"]:
        # Frozen t0 measurements (counts, imported scores) live in the bundle and
        # guidebook; the evolving state carries neither.
        state = {k: deepcopy(v) for k, v in c.items() if k not in ("expression", "measured_scores")}
        state.update(status="alive", activated=c.get("initial_state", {}).get("activated", False),
                     exhausted=c.get("initial_state", {}).get("exhausted", False), signal=0.0,
                     guidebook_id=c["cell_id"], generation=0)
        if dynamics.enabled(bundle["context"]):
            dynamics.initialise_clocks(state, bundle["context"], seeded_rng(seed, "initial_clocks", c["cell_id"]))
        if dynamics.enabled(bundle["context"]) and "t0_geometry" in state:
            state["geometry_model"] = ("provided_t0_rigid_translation" if state["t0_geometry"]["status"] == "valid" else "centroid_fallback")
        if lattice.enabled(bundle["context"]):
            lattice.initialize_bin(state, bundle["measurement"], bundle["context"])
        result[c["cell_id"]] = state
    return result


def decision_rows(cells, guidebook, context, *, seed, step):
    graph = interaction_graph(cells, context, seed=seed, step=step)
    spatial = context.get("interaction_mode", "spatial" if context["assay"] == "xenium" else "well_mixed") == "spatial"
    kinetics = context.get("kinetics", {"mode": "legacy_per_decision"})
    hazard_mode, dt = kinetics["mode"] == "hazard_per_hour", context["time_step_hours"]
    v3 = dynamics.enabled(context)
    lattice_mode = lattice.enabled(context)
    # v5 policy bundles only: one shared record describing imported scores for the provider prompt.
    definitions = guidebook.get("provenance", {}).get("measured_definitions") if v3 else None
    rows = []
    for cid in sorted(cells):
        cell = cells[cid]
        if cell["status"] != "alive":
            continue
        neighbors = [cells[n] for n in graph[cid]]
        signal_neighbors = [n for n in neighbors if n["status"] == "alive"] if v3 else neighbors
        incoming_signal = sum(n["signal"] for n in signal_neighbors)/max(1, len(signal_neighbors))
        incoming_suppressive_signal = sum(n.get("suppressive_signal", 0.0) for n in signal_neighbors)/max(1, len(signal_neighbors)) if v3 else 0.0
        rule = dynamics.lineage_rule(cell, context) if v3 else {}
        density, local_density = dynamics.density_factor(cell, cells, context) if v3 and not lattice_mode else (1.0, 0)
        # Target order decides which neighbor a KILL or CLEAR_CORPSE acts on first.
        # Spatial graphs are nearest-first (interaction_graph sorts by distance,
        # then ID); well-mixed neighborhoods have no distance, so the order is a
        # seeded per-cell draw. Lexicographic cell-ID order would be a hidden bias.
        if spatial:
            ordered, target_rule = list(graph[cid]), ("nearest_centroid_seeded_tie" if v3 else "nearest_centroid_then_id")
            if v3 and not lattice_mode and context.get("geometry", {}).get("interaction_distance") == "boundary":
                target_rule = "nearest_boundary_with_centroid_fallback_seeded_tie"
        else:
            ordered, target_rule = list(graph[cid]), "seeded_uniform_draw"
            seeded_rng(seed, "target_order", step, cid).shuffle(ordered)
        licensed_targets = [n for n in ordered if cells[n]["status"] == "alive" and
                            (dynamics.lineage_rule(cells[n], context)["can_be_killed_target"] if v3 else cells[n]["lineage"] == "malignant")]
        corpses = [n for n in ordered if cells[n]["status"] == "dead"]
        division_targets = lattice.expansion_targets(cell, ordered, cells, context) if lattice_mode else []
        entry = guidebook["cells"][cell["guidebook_id"]]
        actions, prior, excluded = {}, {}, {}
        for action, spec in sorted(entry["actions"].items()):
            reason = None
            if spec["weight"] <= 0:
                reason = "zero guidebook weight"
            elif action == "KILL" and (not (rule.get("can_kill") if v3 else cell["lineage"] == "T_CD8") or not context.get("recognition", False) or
                                       not cell["activated"] or cell["exhausted"] or not licensed_targets):
                reason = "requires CD8, declared recognition, activated/nonexhausted state and a local live malignant target"
            elif lattice_mode and action in {"MOVE", "CLEAR_CORPSE"}:
                reason = "lattice bins neither move nor leave cell corpses"
            elif lattice_mode and action == "PROLIFERATE" and not division_targets:
                reason = "no compatible neighboring bin has declared occupancy capacity"
            elif action == "CLEAR_CORPSE" and not corpses:
                reason = "no local corpse"
            elif action == "ACTIVATE" and cell["activated"] and not v3:
                reason = "already activated"
            elif action == "EXHAUST" and (cell["exhausted"] or not cell["activated"]):
                reason = "requires activated, nonexhausted state"
            elif action == "RECOVER" and not cell["exhausted"]:
                reason = "requires exhausted state"
            elif action == "RECOVER" and v3 and cell.get("exhausted_since") is not None and step * dt - cell["exhausted_since"] >= context["dynamics"]["terminal_exhaustion_hours"]:
                reason = "terminal exhaustion duration reached"
            elif action == "MOVE" and v3 and (not rule.get("motile") or not context["dynamics"]["movement_enabled"] or not spatial):
                reason = "requires enabled motility and measured spatial coordinates"
            elif action == "PROLIFERATE" and v3 and rule.get("proliferation_requires_activation", cell["lineage"] == "T_CD8") and not cell["activated"]:
                reason = "declared lineage proliferation requires activation"
            elif action == "PROLIFERATE" and density <= 0:
                reason = "local density gate"
            if reason:
                excluded[action] = reason
            else:
                actions[action] = lattice.DESCRIPTIONS.get(action, spec["description"]) if lattice_mode else spec["description"]
                # A declared toy communication rule; no cytokine concentration or
                # inferred receptor specificity is implied by this abstract signal.
                multiplier = (1 + incoming_signal)/(1 + incoming_suppressive_signal) if action == "ACTIVATE" else 1.0
                if action == "PROLIFERATE":
                    multiplier *= density
                if hazard_mode and action != "REST":
                    prior[action] = spec["hazard_per_hour"] * multiplier  # hazard, converted below
                else:
                    prior[action] = spec["weight"] * multiplier
        if not actions or "REST" not in actions:
            raise ValueError("Every living cell must have a positive REST prior.")
        hazards = {}
        if hazard_mode:
            # An unavailable action's hazard is zero; its share returns to REST rather
            # than being redistributed among the remaining actions.
            hazards = {a: k for a, k in prior.items() if a != "REST"}
            prior, _ = competing_risk_probabilities(hazards, dt, kinetics["max_step_probability"],
                                                    where=f"Cell {cid!r} at step {step}")
        else:
            total = sum(prior.values())
            prior = {a: v/total for a, v in prior.items()}
        observation = dict(cell_id=cid, cell_type=cell["cell_type"], lineage=cell["lineage"],
            step=step, time_hours=step*context["time_step_hours"],
            state={k: cell[k] for k in ("activated", "exhausted", "signal", "generation")},
            t0_evidence=entry.get("evidence", {}),
            local=dict(neighbor_counts=dict(Counter(n["lineage"] for n in neighbors if n["status"] == "alive")),
                       live_malignant_targets=len(licensed_targets), corpses=len(corpses),
                       incoming_signal=incoming_signal, **(dict(local_density_count=local_density, proliferation_density_factor=density, incoming_suppressive_signal=incoming_suppressive_signal) if v3 else {})))
        if v3 and context.get("geometry", {}).get("interaction_distance") == "boundary" and not lattice_mode:
            observation["local"].update(geometry_distance="boundary_when_both_valid_else_centroid",
                centroid_fallback_neighbors=sum(cell.get("t0_geometry", {}).get("status") != "valid" or n.get("t0_geometry", {}).get("status") != "valid" for n in neighbors))
        if lattice_mode:
            observation.update(unit="binned_spot")
            observation["state"]["occupancy"] = cell["occupancy"]
            observation["local"]["counts_unit"] = "occupied bins, not cells"
        if v3:
            observation["state"].update({k: cell.get(k) for k in ("activated_since", "activated_until", "exhausted_since", "moving_until", "suppressive_signal")})
        rows.append(dict(cell_id=cid, actions=actions, prior=prior, observation=observation,
                         excluded=excluded, targets=licensed_targets, corpses=corpses, target_rule=target_rule,
                         **(dict(division_targets=division_targets) if lattice_mode else {}),
                         **(dict(effective_hazards_per_hour=hazards) if v3 and hazard_mode else {}),
                         **(dict(evidence_definitions=definitions) if definitions else {})))
    return rows, graph


DECISION_LOG_FORMAT = "t0-evidence-by-reference-1"
DECISION_LOG_MODES = ("reference", "inline")
_EVIDENCE_FILE = "guidebook.json"


def evidence_reference(guidebook_id, evidence):
    """Pointer to one guidebook entry's frozen t0 evidence, with its content hash."""
    return dict(file=_EVIDENCE_FILE, guidebook_id=guidebook_id, sha256=digest(evidence))


def logged_observation(observation, reference):
    """Decision-log view of an observation: t0 evidence by reference, same key order.

    Frozen t0 evidence is identical at every step and for every descendant of a
    cell, and ``guidebook.json`` in the run folder already holds it. Only the
    stored log changes; the provider payload is built from the full observation.
    """
    return {("t0_evidence_ref" if key == "t0_evidence" else key): (reference if key == "t0_evidence" else value)
            for key, value in observation.items()}


def read_decisions(directory, *, expand=True):
    """Yield the decisions of a run folder; ``expand`` restores the inline v4 form.

    A by-reference record is rebuilt from ``guidebook.json`` in the same folder
    and verified against the hash stored with the reference, so a missing
    guidebook, or one whose evidence differs from what the run recorded, fails
    instead of reconstructing other evidence. This is an integrity check against
    mix-ups and accidental edits; it does not detect someone rewriting both files.
    Each expanded record owns an independent copy of its evidence, in the
    original key order. Inline logs (schema 2.0, v4 runs, ``--decision-log
    inline``) are yielded unchanged.
    """
    directory = Path(directory)
    entries, verified = None, {}
    with (directory / "decisions.jsonl").open(encoding="utf8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            decision = json.loads(line)
            if not isinstance(decision, dict):
                raise ValueError(f"decisions.jsonl line {number} is not a decision record.")
            observation = decision.get("observation")
            if expand and isinstance(observation, dict) and "t0_evidence_ref" in observation:
                reference = observation["t0_evidence_ref"]
                if "t0_evidence" in observation:
                    raise ValueError(f"decisions.jsonl line {number} holds both inline and referenced t0 evidence.")
                if (not isinstance(reference, dict) or reference.get("file") != _EVIDENCE_FILE
                        or not isinstance(reference.get("guidebook_id"), str)
                        or not isinstance(reference.get("sha256"), str)):
                    raise ValueError(f"decisions.jsonl line {number} has an invalid t0 evidence reference; "
                                     f"only {_EVIDENCE_FILE} in the run folder is accepted.")
                key = reference["guidebook_id"]
                if key not in verified:
                    if entries is None:
                        path = directory / _EVIDENCE_FILE
                        if not path.is_file():
                            raise ValueError(f"Expanding this decision log requires {_EVIDENCE_FILE} beside decisions.jsonl.")
                        document = json.loads(path.read_text(encoding="utf8"))
                        entries = document.get("cells") if isinstance(document, dict) else None
                        if not isinstance(entries, dict):
                            raise ValueError(f"{_EVIDENCE_FILE} has no cells map; it is not this run's guidebook.")
                    if not isinstance(entries.get(key), dict):
                        raise ValueError(f"{_EVIDENCE_FILE} has no entry {key!r}; it is not this run's guidebook.")
                    evidence = entries[key].get("evidence", {})
                    verified[key] = (digest(evidence), json.dumps(evidence, allow_nan=False))
                checksum, serialized = verified[key]
                if checksum != reference["sha256"]:
                    raise ValueError(f"t0 evidence of {key!r} in {_EVIDENCE_FILE} does not match the hash recorded "
                                     f"in decisions.jsonl line {number}; the guidebook was replaced or edited.")
                decision["observation"] = {
                    ("t0_evidence" if name == "t0_evidence_ref" else name):
                        (json.loads(serialized) if name == "t0_evidence_ref" else value)
                    for name, value in observation.items()}
            yield decision


JEV_POOLING_MODES = ("full_distribution", "event_share")


def distribution(prior, answer=None, jev_weight=.5, pooling="full_distribution"):
    """Explicit log-linear preference pooling, not Bayesian biological inference.

    ``full_distribution`` (the v4 rule, default) pools every action including
    REST, so Jev's REST preference also sets how often events happen and the
    event rate depends on the time step. ``event_share`` (opt-in) keeps the
    kernel's event mass, 1 - prior[REST] = 1 - exp(-sum(k) dt), and pools only
    the shares among non-REST actions; Jev then decides which event, not when.
    """
    if answer is None or jev_weight == 0:
        return dict(prior)
    p = answer["probabilities"]
    if set(p) != set(prior):
        raise ValueError("Provider actions differ from the gated menu.")
    if pooling == "event_share":
        return _event_share(prior, p, jev_weight)
    if pooling != "full_distribution":
        raise ValueError(f"jev_pooling must be one of {JEV_POOLING_MODES}.")
    if jev_weight == 1:
        return dict(p)
    weights = {a: prior[a] ** (1-jev_weight) * p[a] ** jev_weight for a in prior}
    total = sum(weights.values())
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Combined distribution is empty or invalid.")
    return {a: w/total for a, w in weights.items()}


def _event_share(prior, p, jev_weight):
    """P(a) = (1 - prior[REST]) * share_a; REST keeps the kernel's exp(-sum(k) dt).

    share_a is the normalised pool of the prior share and Jev's share among the
    non-REST actions. If Jev puts no mass on any event, or the pool is empty,
    the prior shares are kept: the provider cannot veto the kernel's clock.
    """
    events = [a for a in prior if a != "REST"]
    event_mass = 1.0 - prior["REST"]
    prior_events, jev_events = sum(prior[a] for a in events), sum(p[a] for a in events)
    if not events or event_mass <= 0 or prior_events <= 0 or jev_events <= 0:
        return dict(prior)
    pooled = {a: (prior[a] / prior_events) ** (1 - jev_weight) * (p[a] / jev_events) ** jev_weight for a in events}
    total = sum(pooled.values())
    if not math.isfinite(total) or total <= 0:
        return dict(prior)
    return {a: (prior["REST"] if a == "REST" else event_mass * pooled[a] / total) for a in prior}


def js_divergence(p, q):
    """Jensen-Shannon divergence in bits between two distributions on the same menu."""
    keys = sorted(set(p) | set(q))
    m = {a: (p.get(a, 0.0) + q.get(a, 0.0)) / 2 for a in keys}
    def kl(x):
        return sum(x.get(a, 0.0) * math.log2(x.get(a, 0.0) / m[a]) for a in keys if x.get(a, 0.0) > 0 and m[a] > 0)
    return (kl(p) + kl(q)) / 2


def draw(probabilities, rng):
    u, cumulative = rng.random(), 0.0
    for action, p in sorted(probabilities.items()):
        cumulative += p
        if u < cumulative:
            return action
    return sorted(probabilities)[-1]


def commit(cells, proposals, *, max_population, seed, step, context=None):
    """Deaths precede divisions. Target conflicts resolve by a seeded stable order.

    Each proposal's ``targets``/``corpses`` arrive already ordered by the rule
    recorded in ``target_rule`` (nearest centroid in spatial mode, a seeded draw
    in well-mixed mode); the first still-available entry is acted on.
    """
    context = context or {}
    if lattice.enabled(context):
        return lattice.commit(cells, proposals, context, seed=seed, step=step, digest=digest, seeded_rng=seeded_rng)
    v3 = dynamics.enabled(context)
    dt = context.get("time_step_hours", 1.0)
    start, now = step * dt, (step + 1) * dt
    next_cells = deepcopy(cells)
    events, doomed, cleared = [], set(), set()
    for cell in next_cells.values():
        retention = math.exp(-math.log(2) * dt / context["dynamics"]["signal_half_life_hours"]) if v3 else .5
        cell["signal"] *= retention
        if v3:
            cell["suppressive_signal"] = cell.get("suppressive_signal", 0.0) * retention
        if v3 and cell["status"] != "alive":
            cell["signal"] = 0.0
            cell["suppressive_signal"] = 0.0
    order = sorted(proposals, key=lambda p: digest([seed, "commit", step, p["cell_id"]]))
    for p in order:
        cid, action = p["cell_id"], p["action"]
        if action == "APOPTOSIS":
            if cid in doomed:
                events.append(dict(kind="conflict_blocked", cell_id=cid, action=action, reason="already_doomed_this_round"))
            else:
                doomed.add(cid)
                events.append(dict(kind="apoptosis", cell_id=cid))
        elif action == "KILL":
            targets = [t for t in p["targets"] if t not in doomed]
            if targets:
                target = targets[0]
                doomed.add(target)
                events.append(dict(kind="kill", cell_id=cid, target=target, target_rule=p.get("target_rule")))
            else:
                events.append(dict(kind="conflict_blocked", cell_id=cid, action=action))
        elif action == "CLEAR_CORPSE":
            targets = [t for t in p["corpses"] if t not in cleared]
            if targets:
                cleared.add(targets[0])
                events.append(dict(kind="clear_corpse", cell_id=cid, target=targets[0], target_rule=p.get("target_rule")))
            else:
                events.append(dict(kind="conflict_blocked", cell_id=cid, action=action))
    for cid in doomed:
        next_cells[cid]["status"] = "dead"
        if v3:
            next_cells[cid].update(dead_since=now, signal=0.0, suppressive_signal=0.0, moving_until=now)
    for cid in cleared:
        next_cells[cid]["status"] = "removed"
        if v3:
            next_cells[cid]["signal"] = 0.0
            next_cells[cid]["suppressive_signal"] = 0.0
    living = sum(c["status"] == "alive" for c in next_cells.values())
    for p in order:
        cid, action = p["cell_id"], p["action"]
        c = next_cells[cid]
        if c["status"] != "alive":
            if action not in ("APOPTOSIS",):
                events.append(dict(kind="actor_died_same_step", cell_id=cid, action=action))
            continue
        if action == "PROLIFERATE":
            if living >= max_population:
                events.append(dict(kind="population_cap_blocked", cell_id=cid))
                continue
            daughter_id = f"{cid}/d{step+1}"
            while daughter_id in next_cells:
                daughter_id += "d"
            daughter = deepcopy(c)
            daughter.update(cell_id=daughter_id, generation=c["generation"]+1)
            if v3:
                dynamics.place_daughter(c, daughter, next_cells, context, seeded_rng(seed, "daughter", step, cid))
                if daughter.get("t0_geometry", {}).get("status") == "valid":
                    daughter["geometry_model"] = "synthetic_inherited_rigid_translation"
                    daughter["shape_inherited_from"] = cid
            elif "x" in daughter:
                angle = seeded_rng(seed, "daughter", step, cid).random() * math.tau
                daughter["x"] += math.cos(angle) * 2
                daughter["y"] += math.sin(angle) * 2
            next_cells[daughter_id] = daughter
            living += 1
            events.append(dict(kind="division", cell_id=cid, daughter=daughter_id,
                               **(dict(daughter_shape="synthetic_inherited_rigid_translation", shape_source_cell_id=daughter["t0_geometry"]["source_cell_id"]) if v3 and daughter.get("t0_geometry", {}).get("status") == "valid" else {})))
        elif action == "ACTIVATE":
            c["activated"] = True
            if v3:
                duration = seeded_rng(seed, "activation_clock", step, cid).expovariate(math.log(2) / context["dynamics"]["activation_half_life_hours"])
                c.update(activated_since=now, activated_until=now + duration)
                events.append(dict(kind="activate", cell_id=cid, activated_until=c["activated_until"]))
        elif action == "EXHAUST":
            c["exhausted"] = True
            if v3:
                c["exhausted_since"] = now
                events.append(dict(kind="exhaust", cell_id=cid))
        elif action == "RECOVER":
            c["exhausted"] = False
            if v3:
                c["exhausted_since"] = None
                events.append(dict(kind="recover", cell_id=cid))
        elif action == "SECRETE":
            channel = "suppressive_signal" if v3 and dynamics.lineage_rule(c, context).get("secretion_kind") == "suppressive" else "signal"
            c[channel] = c.get(channel, 0.0) + 1.0
    if v3:
        selected_moves = {p["cell_id"] for p in proposals if p["action"] == "MOVE"}
        # Use the shared old snapshot for chemotaxis, independent of commit order.
        for cid in sorted(next_cells):
            c = next_cells[cid]
            if c["status"] == "alive" and cid in cells:
                event = dynamics.move(c, cells, context, start=start, end=now, selected=cid in selected_moves,
                                      rng=seeded_rng(seed, "move", step, cid))
                if event is not None:
                    events.append(event)
        dynamics.apply_clocks(next_cells, context, now, events)
    return next_cells, events


def kinetics_summary(bundle, rows):
    """Implied event rates at t0 so the user can judge the time scale before running.

    In hazard mode the per-lineage mean hazard and the implied population
    doubling time (ln 2 / mean PROLIFERATE hazard, ignoring death) are reported
    in hours. In legacy mode only per-step probabilities exist, and the
    doubling time is per step, not per hour.
    """
    context = bundle["context"]
    mode, dt = context["kinetics"]["mode"], context["time_step_hours"]
    guidebook = bundle["guidebook"]["cells"]
    per_lineage = {}
    for row in rows:
        lineage = row["observation"]["lineage"]
        entry = guidebook[row["observation"]["cell_id"]]
        acc = per_lineage.setdefault(lineage, dict(n=0, hazard_sum={}, probability_sum={}))
        acc["n"] += 1
        for action, spec in entry["actions"].items():
            if action == "REST":
                continue
            acc["probability_sum"][action] = acc["probability_sum"].get(action, 0.0) + spec["probability"]
            if mode == "hazard_per_hour":
                acc["hazard_sum"][action] = acc["hazard_sum"].get(action, 0.0) + spec["hazard_per_hour"]
    out = {}
    for lineage, acc in sorted(per_lineage.items()):
        n = acc["n"]
        item = dict(cells=n, mean_step_probability={a: v / n for a, v in sorted(acc["probability_sum"].items())})
        if mode == "hazard_per_hour":
            item["mean_hazard_per_hour"] = {a: v / n for a, v in sorted(acc["hazard_sum"].items())}
            k = item["mean_hazard_per_hour"].get("PROLIFERATE", 0.0)
            item["implied_doubling_hours_ignoring_death"] = (math.log(2) / k) if k > 0 else None
        else:
            p = item["mean_step_probability"].get("PROLIFERATE", 0.0)
            item["implied_doubling_steps_ignoring_death"] = (math.log(2) / math.log1p(p)) if p > 0 else None
        out[lineage] = item
    max_event_probability = max((1 - r["prior"]["REST"] for r in rows), default=0.0)
    max_burden = -math.log(1 - max_event_probability) if max_event_probability < 1 else float("inf")
    provenance = bundle["guidebook"]["provenance"]
    effective_rates, expected_rates = {}, {}
    if dynamics.enabled(context) and mode == "hazard_per_hour":
        for row in rows:
            lineage = row["observation"]["lineage"]
            rates = effective_rates.setdefault(lineage, {})
            expected = expected_rates.setdefault(lineage, {})
            for action, value in row.get("effective_hazards_per_hour", {}).items():
                rates[action] = rates.get(action, 0.0) + value
            for action, value in row["prior"].items():
                if action != "REST":
                    expected[action] = expected.get(action, 0.0) + value / dt
    result = dict(mode=mode, time_step_hours=dt, by_lineage=out, max_cell_event_probability_t0=max_event_probability,
                max_cell_step_burden_t0=max_burden,
                max_step_probability=context["kinetics"].get("max_step_probability"),
                gated_worst_case_step_burden=provenance.get("kinetics_worst_case_step_burden"),
                recommended_max_time_step_hours=provenance.get("kinetics_recommended_max_time_step_hours"),
                expected_rate_undercount_fraction=(max_burden / 2 if mode == "hazard_per_hour" else None),
                discretisation="at most one event per cell per step; rates are under-counted by about burden/2",
                kinetics_review_status=provenance.get("kinetics_review_status"),
                note=("Hazards are per hour; time_step_hours changes the trajectory." if mode == "hazard_per_hour" else
                      "Legacy mode: weights are per decision, time_step_hours is a label only."))
    if dynamics.enabled(context) and mode == "hazard_per_hour":
        result.update(gated_hazard_per_hour_by_lineage=effective_rates,
                      expected_events_per_hour_by_lineage=expected_rates,
                      expected_events_semantics="t0 sum of gated first-event probabilities divided by dt; local prior before Jev; successful commits may be reduced by conflicts; lattice events act on occupancy, not cells")
    return result


def _action_evidence(bundle):
    from ..v5.evidence import action_evidence
    return action_evidence(bundle)


def preflight(bundle, *, model="jev-1.14.0", max_questions=20, seed=21, prior_mode="shown"):
    validate_prepared(bundle)
    cells = initialize(bundle, seed=seed)
    rows, graph = decision_rows(cells, bundle["guidebook"], bundle["context"], seed=seed, step=0)
    active_rows = [r for r in rows if len(r["actions"]) > 1]
    groups = batches(active_rows, bundle["context"], model, max_questions, prior_mode=prior_mode)
    result = dict(schema_version=bundle["schema_version"], cells=len(cells), by_lineage=dict(Counter(c["lineage"] for c in cells.values())),
        initial_questions=len(active_rows), initial_http_batches=len(groups),
        kinetics=kinetics_summary(bundle, rows),
        scoring=dict(method=bundle["context"]["scoring"]["method"],
                     coverage=bundle["guidebook"]["provenance"].get("scoring_coverage"),
                     qc_fail_cells=bundle["guidebook"]["provenance"].get("qc_fail_cells", 0),
                     note="Statuses 'available' and 'available_panel_heuristic' can support expression modulation; the latter is weak panel calibration. Unavailable evidence keeps multiplier 1; available evidence may still be neutral."),
        isolated_cells=sum(not ns for ns in graph.values()), interaction_mode=bundle["context"].get("interaction_mode"),
        **(dict(action_evidence=_action_evidence(bundle)) if bundle["schema_version"] == "3.0" else {}),
        limitations=[("Hazards and state clocks have hour units but remain uncalibrated." if bundle["context"]["kinetics"]["mode"] == "hazard_per_hour" else "Legacy per-decision weights use time_step_hours as labels only."),
                    "t0 expression stays fixed; state transitions do not predict a new transcriptome.",
                    ("Spatial edges use declared outline gap with centroid fallback; segmentation proximity does not establish physical contact."
                     if bundle["context"].get("geometry", {}).get("interaction_distance") == "boundary" and not lattice.enabled(bundle["context"])
                     else "Spatial edges use centroid distance; they do not establish physical contact."),
                    "Enrichment/subsampling does not reconstruct absent populations or population frequencies."])
    if dynamics.enabled(bundle["context"]) and ("geometry" in bundle["context"] or any("t0_geometry" in c for c in bundle["cells"])):
        result["geometry"] = dict(geometry.summarize_geometry(bundle["cells"]),
            interaction_distance=bundle["context"].get("geometry", {}).get("interaction_distance", "centroid"),
            lattice_note="Lattice uses fixed-bin adjacency; polygon display does not change lattice topology." if lattice.enabled(bundle["context"]) else None)
    return result


def run(bundle, output_dir, **kwargs):
    """Acquire an exclusive run lock; unresolved process exits require inspection."""
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    lock = directory / ".running.lock"
    try:
        with lock.open("x", encoding="utf8") as f:
            f.write(str(os.getpid()))
    except FileExistsError:
        raise ValueError("Run lock exists. Check that no process is running before manually removing .running.lock.") from None
    try:
        return _run(bundle, output_dir, **kwargs)
    finally:
        lock.unlink()


def _run(bundle, output_dir, *, steps=12, seed=21, mode="local", provider=None, jev_weight=.5,
        max_population=2000, resume=False, decision_log="reference", jev_pooling="full_distribution"):
    validate_prepared(bundle)
    if decision_log not in DECISION_LOG_MODES:
        raise ValueError(f"decision_log must be one of {DECISION_LOG_MODES}.")
    if jev_pooling not in JEV_POOLING_MODES:
        raise ValueError(f"jev_pooling must be one of {JEV_POOLING_MODES}.")
    if jev_pooling != "full_distribution" and mode != "jev":
        raise ValueError("jev_pooling applies only to Jev mode.")
    if jev_pooling != "full_distribution" and bundle["schema_version"] != "3.0":
        raise ValueError("event_share keeps the hazard kernel's event probability and requires schema 3.0.")
    if mode == "jev" and getattr(provider, "pooling", "full_distribution") != jev_pooling:
        raise ValueError("Provider and run disagree about jev_pooling; the prompt must describe the applied rule.")
    context = bundle["context"]
    v3 = dynamics.enabled(context)
    if type(steps) is not int or not 1 <= steps <= 10000 or type(seed) is not int:
        raise ValueError("Use 1–10000 steps and an integer seed.")
    if v3 and steps * context["time_step_hours"] > context["dynamics"]["max_horizon_hours"]:
        raise ValueError("Requested horizon exceeds dynamics.max_horizon_hours; reduce steps or declare a reviewed longer horizon.")
    if mode not in ("local", "jev") or (mode == "jev" and provider is None):
        raise ValueError("mode must be local, or jev with a provider.")
    if type(jev_weight) not in (int, float) or not math.isfinite(jev_weight) or not 0 <= jev_weight <= 1:
        raise ValueError("jev_weight must be in [0,1].")
    if mode == "jev" and getattr(provider, "prior_mode", None) == "hidden":
        # Hiding the prompt prior must also remove prior pooling; otherwise the
        # supposedly prior-blind control remains anchored in the same guidebook.
        jev_weight = 1.0
    if type(max_population) is not int or max_population < len(bundle["cells"]):
        raise ValueError("max_population must cover all starting cells.")
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    config = dict(schema_version=bundle["schema_version"], bundle_fingerprint=bundle["fingerprint"], steps=steps, seed=seed,
        mode=mode, jev_weight=jev_weight, max_population=max_population,
        model=provider.model if provider else None,
        provider_name=getattr(provider, "name", None), provider_prompt_version=getattr(provider, "prompt_version", None),
        provider_sampling_parameters=deepcopy(getattr(provider, "sampling_parameters", None)),
        jev_prior_mode=getattr(provider, "prior_mode", None),
        max_questions=getattr(provider, "max_questions", None),
        probability_semantics="guidebook prior (local) or normalized prior^(1-alpha)*Jev^alpha; seeded categorical sample",
        kinetics=bundle["context"]["kinetics"],
        kernel="vtx-discrete-context-2.0.1", assumptions=dict(signal_retention_per_step=.5, signal_increment=1.0,
            signal_activation_multiplier="1 + mean neighbor signal from previous snapshot",
            kill="one licensed KILL commits one target death per step; target = nearest live malignant neighbor (spatial) or seeded draw (well_mixed)",
            clearance="target = nearest corpse (spatial) or seeded draw (well_mixed)",
            division="one daughter per selected action per step",
            geometry="fixed centroids; daughter offset 2 um; no mechanics; capped nearest-neighbor graph",
            simultaneous="snapshot proposals; death/clearance conflicts first; surviving actors divide/change state"))
    if v3:
        config.update(kernel="vtx-cell-agent-3.0.1", dynamics=context["dynamics"],
                      ontology_version=context["ontology_version"], measurement=bundle.get("measurement"))
        config["assumptions"].update(
            signal_activation_multiplier="(1 + mean live-neighbor inflammatory signal)/(1 + mean live-neighbor suppressive signal); abstract channels without molecular receptor inference",
            signal_retention_per_step=math.exp(-math.log(2) * context["time_step_hours"] / context["dynamics"]["signal_half_life_hours"]),
            geometry="symmetric bounded-degree live graph; corpse edges independent; ontology radii; local newborn soft repulsion",
            motility="MOVE starts a persistent episode; speed integrated over active hours; abstract signal bias optional",
            activation="exponential lifetime with declared half-life; daughters inherit clocks",
            exhaustion="RECOVER gated after declared terminal duration",
            crowding="PROLIFERATE hazard multiplied by max(0,1-local_live_neighbors/n_max)",
            corpse="zero signal; automatic removal after declared lifetime",
            frame_storage="committed JSONL/JSON array; lazy offline replay with 50-frame sidecars and a two-chunk cache",
            resume="reconstruct from t0 with the configuration-bound provider journal; no checkpoint file")
    if v3 and ("geometry" in context or any("t0_geometry" in c for c in bundle["cells"])):
        config["geometry"] = geometry.validate_geometry_config(context.get("geometry"))
        config["assumptions"].update(segmentation="provided t0 polygons with QC; invalid/missing contours use centroids; no contour repair or inferred histology",
            shape_evolution="rigid translation; daughters inherit an explicitly synthetic copy of the ancestral t0 shape",
            interaction_distance=config["geometry"]["interaction_distance"],
            distance_cost="bounding-box candidate sweep; exact polygon gap worst-case O(vertices_a*vertices_b) per pair")
    if jev_pooling != "full_distribution":
        # Recorded only when it departs from the v4 rule, so default configurations
        # and their hashes stay identical to v4.
        config.update(jev_pooling=jev_pooling,
            probability_semantics="kernel event mass 1-exp(-sum(k)dt) kept; Jev pooled only among non-REST shares: "
                                  "share_a ~ prior_share^(1-alpha)*Jev_share^alpha; seeded categorical sample")
    if lattice.enabled(context):
        config.update(kernel="vtx-lattice-occupancy-3.0", lattice=context["lattice"])
        config["assumptions"].update(geometry="fixed measured bin centroids; occupancy changes on existing bins only",
            division="increment compatible empty/low-density neighboring bin occupancy; no daughter IDs",
            kill="declared killing empties one eligible target bin; occupancy loss, not a cell death count",
            motility="disabled for fixed bins", crowding="declared per-bin capacity; no cell-radius density gate",
            corpse="empty bins have zero occupancy and signal, no cell corpses")
    identity = digest(config)
    old_config = directory / "config.json"
    if old_config.exists():
        if not resume or json.loads(old_config.read_text()) != config:
            raise ValueError("Output already has a run. Use a new directory or --resume with identical settings.")
    else:
        permitted = {".running.lock"}
        if provider is not None and getattr(provider, "directory", None) == directory / "jev_journal":
            permitted.add("jev_journal")
        if any(p.name not in permitted for p in directory.iterdir()):
            raise ValueError("Output directory contains unrelated files; use a new empty run directory.")
    write_json(old_config, config)
    write_json(directory / "guidebook.json", bundle["guidebook"])
    cells = initialize(bundle, seed=seed)
    frames = FrameStore(directory) if v3 else []
    frames.append(dict(step=0, time_hours=0, cells=cells if v3 else deepcopy(cells), events=[]))
    decisions, stop, error, divergences = [], "completed", None, {}
    by_reference, evidence_refs, event_mass = v3 and decision_log == "reference", {}, {}
    decision_file = (directory / "decisions.jsonl").open("w", encoding="utf8") if v3 else None
    completed_steps, event_counts = 0, Counter()
    try:
        for step in range(steps):
            try:
                rows, graph = decision_rows(cells, bundle["guidebook"], bundle["context"], seed=seed, step=step)
                if not rows and not (v3 and any(c["status"] == "dead" for c in cells.values())):
                    stop = "no_living_cells"
                    break
                active_rows = [r for r in rows if len(r["actions"]) > 1]
                answers = provider.choose(active_rows, bundle["context"]) if mode == "jev" else {}
                if mode == "jev" and set(answers) != {r["cell_id"] for r in active_rows}:
                    raise ProviderError("Provider omitted or added cell answers; no state committed.")
            except (ProviderError, BudgetExhausted) as e:
                stop, error = "provider_stopped_before_round_commit", str(e)
                break
            except KineticsStepError as e:
                stop, error = "kinetics_step_too_coarse", str(e)
                break
            proposals = []
            if v3:
                decisions = []
            for row in rows:
                cid = row["cell_id"]
                answer = answers.get(cid)
                probabilities = distribution(row["prior"], answer, jev_weight, jev_pooling)
                action = draw(probabilities, seeded_rng(seed, "decision", step, cid))
                divergence = js_divergence(row["prior"], answer["probabilities"]) if answer else None
                if divergence is not None:
                    stats = divergences.setdefault(row["observation"]["lineage"], dict(n=0, total=0.0, maximum=0.0))
                    stats["n"] += 1
                    stats["total"] += divergence
                    stats["maximum"] = max(stats["maximum"], divergence)
                    # How much event probability the kernel declared, Jev stated and the run applied.
                    mass = event_mass.setdefault(row["observation"]["lineage"], dict(n=0, prior=0.0, jev=0.0, applied=0.0))
                    mass["n"] += 1
                    mass["prior"] += 1 - row["prior"]["REST"]
                    mass["jev"] += 1 - answer["probabilities"]["REST"]
                    mass["applied"] += 1 - probabilities["REST"]
                observation = row["observation"]
                if by_reference:
                    key = cells[cid]["guidebook_id"]
                    if key not in evidence_refs:
                        evidence_refs[key] = evidence_reference(key, observation["t0_evidence"])
                    observation = logged_observation(observation, evidence_refs[key])
                proposals.append(dict(cell_id=cid, action=action, targets=row["targets"], corpses=row["corpses"],
                                      target_rule=row["target_rule"],
                                      **(dict(division_targets=row["division_targets"]) if "division_targets" in row else {})))
                decisions.append(dict(step=step, time_hours=step*bundle["context"]["time_step_hours"], cell_id=cid,
                    action=action, prior=row["prior"], jev=answer, applied_probabilities=probabilities,
                    jev_prior_js_bits=divergence,
                    observation=observation, excluded=row["excluded"],
                    provenance="local_guidebook" if mode == "local" else ("jev_pooled_preferences" if answer else "single_permitted_action")))
            cells, events = commit(cells, proposals, max_population=max_population, seed=seed, step=step, context=context)
            frames.append(dict(step=step+1, time_hours=(step+1)*context["time_step_hours"], cells=cells if v3 else deepcopy(cells), events=events))
            completed_steps = step + 1
            event_counts.update(e["kind"] for e in events)
            if v3:
                for decision in decisions:
                    decision_file.write(json.dumps(decision, allow_nan=False) + "\n")
                decision_file.flush()
            # v2 compatibility artifact only: resume reconstructs from t0 using
            # its durable provider journal, never by reading this snapshot.
            if not v3:
                write_json(directory / "checkpoint.json", dict(config_hash=identity, completed_steps=step+1, cells=cells))
    finally:
        if v3:
            frames.close()
            decision_file.close()
    report = dict(schema_version=bundle["schema_version"], completed=completed_steps == steps, stop_reason=stop, error=error,
        completed_steps=completed_steps, requested_steps=steps, seed=seed, mode=mode,
        living=sum(c["status"] == "alive" for c in cells.values()), dead=sum(c["status"] == "dead" for c in cells.values()),
        removed=sum(c["status"] == "removed" for c in cells.values()), allocated=len(cells),
        events=dict(event_counts),
        bundle_fingerprint=bundle["fingerprint"], config_hash=identity,
        usage=provider.usage() if provider else dict(http_attempts=0),
        jev_prior_mode=getattr(provider, "prior_mode", None),
        jev_prior_divergence=dict(unit="Jensen-Shannon divergence, bits, Jev vs gated guidebook prior",
            by_lineage={lineage: dict(n=v["n"], mean=v["total"]/v["n"], max=v["maximum"]) for lineage, v in sorted(divergences.items())},
            interpretation="Near zero with prior shown and larger with prior hidden indicates anchoring; compare against replicate variance."),
        biological_validation="not established", preflight=preflight(bundle, seed=seed,
            model=getattr(provider, "model", "jev-1.14.0"), max_questions=getattr(provider, "max_questions", 20),
            prior_mode=getattr(provider, "prior_mode", "shown")))
    if v3:
        report["decision_log"] = (dict(format=DECISION_LOG_FORMAT, t0_evidence_file=_EVIDENCE_FILE,
                t0_evidence="by reference: cells[observation.t0_evidence_ref.guidebook_id].evidence, verified by sha256",
                reader="virtualtissue_xenium.v2.simulation.read_decisions; CLI: decisions RUN_DIR --output FILE",
                note="Keep guidebook.json with decisions.jsonl; the reader restores the inline v4 records exactly "
                     "and fails if the guidebook's evidence differs from the recorded hashes.")
            if by_reference else dict(format="inline", t0_evidence="copied into every decision (v4 layout)"))
        report["expected_events_per_hour_by_lineage"] = report["preflight"]["kinetics"]["expected_events_per_hour_by_lineage"]
    if mode == "jev" and v3:
        dt = context["time_step_hours"]
        report["jev_pooling"] = jev_pooling
        report["jev_event_mass"] = dict(
            unit="mean probability per decision that the cell does anything but REST; divide by time_step_hours for a rate",
            time_step_hours=dt,
            by_lineage={lineage: dict(n=v["n"], prior=v["prior"]/v["n"], jev=v["jev"]/v["n"], applied=v["applied"]/v["n"],
                                      applied_over_prior=(v["applied"]/v["prior"] if v["prior"] > 0 else None),
                                      prior_events_per_hour=v["prior"]/v["n"]/dt, applied_events_per_hour=v["applied"]/v["n"]/dt)
                        for lineage, v in sorted(event_mass.items())},
            interpretation=("full_distribution pools REST too, so Jev's REST preference sets event frequency and "
                            "applied_over_prior changes with time_step_hours; compare arms only at one time step. "
                            "event_share keeps applied equal to prior by construction."))
    if lattice.enabled(context):
        report.update(kernel="lattice", occupancy=lattice.summary(cells, context),
                      occupied_bins=sum(c["occupancy"] > 0 for c in cells.values()),
                      interpretation="Fixed-bin occupancy model; living/removed legacy fields count occupied/empty bins, not cells.")
    elif v3:
        report["kernel"] = "cell_agent"
    if not v3:
        write_json(directory / "frames.json", frames)
        with (directory / "decisions.jsonl").open("w", encoding="utf8") as f:
            for d in decisions:
                f.write(json.dumps(d, allow_nan=False) + "\n")
    write_json(directory / "report.json", report)
    from .replay import write_replay, write_replay_stream
    if v3:
        write_replay_stream(directory / "replay.html", directory / "frames.jsonl", report)
    else:
        write_replay(directory / "replay.html", frames, report)
    return report
