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

from .inputs import digest, validate_prepared, write_json
from .provider import BudgetExhausted, ProviderError, batches


def seeded_rng(seed, *parts):
    return random.Random(int(digest([seed, *parts]), 16))


def interaction_graph(cells, context, *, seed=21, step=0):
    """No edges across samples; nonspatial neighborhoods are explicit assumptions."""
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


def initialize(bundle):
    result = {}
    for c in bundle["cells"]:
        state = {k: deepcopy(v) for k, v in c.items() if k != "expression"}
        state.update(status="alive", activated=c.get("initial_state", {}).get("activated", False),
                     exhausted=c.get("initial_state", {}).get("exhausted", False), signal=0.0,
                     guidebook_id=c["cell_id"], generation=0)
        result[c["cell_id"]] = state
    return result


def decision_rows(cells, guidebook, context, *, seed, step):
    graph = interaction_graph(cells, context, seed=seed, step=step)
    spatial = context.get("interaction_mode", "spatial" if context["assay"] == "xenium" else "well_mixed") == "spatial"
    rows = []
    for cid in sorted(cells):
        cell = cells[cid]
        if cell["status"] != "alive":
            continue
        neighbors = [cells[n] for n in graph[cid]]
        incoming_signal = sum(n["signal"] for n in neighbors)/max(1, len(neighbors))
        # Target order decides which neighbor a KILL or CLEAR_CORPSE acts on first.
        # Spatial graphs are nearest-first (interaction_graph sorts by distance,
        # then ID); well-mixed neighborhoods have no distance, so the order is a
        # seeded per-cell draw. Lexicographic cell-ID order would be a hidden bias.
        if spatial:
            ordered, target_rule = list(graph[cid]), "nearest_centroid_then_id"
        else:
            ordered, target_rule = list(graph[cid]), "seeded_uniform_draw"
            seeded_rng(seed, "target_order", step, cid).shuffle(ordered)
        licensed_targets = [n for n in ordered if cells[n]["status"] == "alive" and cells[n]["lineage"] == "malignant"]
        corpses = [n for n in ordered if cells[n]["status"] == "dead"]
        entry = guidebook["cells"][cell["guidebook_id"]]
        actions, prior, excluded = {}, {}, {}
        for action, spec in sorted(entry["actions"].items()):
            reason = None
            if spec["weight"] <= 0:
                reason = "zero guidebook weight"
            elif action == "KILL" and (cell["lineage"] != "T_CD8" or not context.get("recognition", False) or
                                       not cell["activated"] or cell["exhausted"] or not licensed_targets):
                reason = "requires CD8, declared recognition, activated/nonexhausted state and a local live malignant target"
            elif action == "CLEAR_CORPSE" and not corpses:
                reason = "no local corpse"
            elif action == "ACTIVATE" and cell["activated"]:
                reason = "already activated"
            elif action == "EXHAUST" and (cell["exhausted"] or not cell["activated"]):
                reason = "requires activated, nonexhausted state"
            elif action == "RECOVER" and not cell["exhausted"]:
                reason = "requires exhausted state"
            if reason:
                excluded[action] = reason
            else:
                actions[action] = spec["description"]
                prior[action] = spec["weight"]
                # A declared toy communication rule; no cytokine concentration or
                # inferred receptor specificity is implied by this abstract signal.
                if action == "ACTIVATE":
                    prior[action] *= 1 + incoming_signal
        if not actions or "REST" not in actions:
            raise ValueError("Every living cell must have a positive REST prior.")
        total = sum(prior.values())
        prior = {a: v/total for a, v in prior.items()}
        observation = dict(cell_id=cid, cell_type=cell["cell_type"], lineage=cell["lineage"],
            step=step, time_hours=step*context["time_step_hours"],
            state={k: cell[k] for k in ("activated", "exhausted", "signal", "generation")},
            t0_evidence=entry.get("evidence", {}),
            local=dict(neighbor_counts=dict(Counter(n["lineage"] for n in neighbors if n["status"] == "alive")),
                       live_malignant_targets=len(licensed_targets), corpses=len(corpses),
                       incoming_signal=incoming_signal))
        rows.append(dict(cell_id=cid, actions=actions, prior=prior, observation=observation,
                         excluded=excluded, targets=licensed_targets, corpses=corpses, target_rule=target_rule))
    return rows, graph


def distribution(prior, answer=None, jev_weight=.5):
    """Explicit log-linear preference pooling, not Bayesian biological inference."""
    if answer is None or jev_weight == 0:
        return dict(prior)
    p = answer["probabilities"]
    if set(p) != set(prior):
        raise ValueError("Provider actions differ from the gated menu.")
    if jev_weight == 1:
        return dict(p)
    weights = {a: prior[a] ** (1-jev_weight) * p[a] ** jev_weight for a in prior}
    total = sum(weights.values())
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Combined distribution is empty or invalid.")
    return {a: w/total for a, w in weights.items()}


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


def commit(cells, proposals, *, max_population, seed, step):
    """Deaths precede divisions. Target conflicts resolve by a seeded stable order.

    Each proposal's ``targets``/``corpses`` arrive already ordered by the rule
    recorded in ``target_rule`` (nearest centroid in spatial mode, a seeded draw
    in well-mixed mode); the first still-available entry is acted on.
    """
    next_cells = deepcopy(cells)
    events, doomed, cleared = [], set(), set()
    for cell in next_cells.values():
        cell["signal"] *= .5  # Explicit per-step toy signal decay.
    order = sorted(proposals, key=lambda p: digest([seed, "commit", step, p["cell_id"]]))
    for p in order:
        cid, action = p["cell_id"], p["action"]
        if action == "APOPTOSIS":
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
    for cid in cleared:
        next_cells[cid]["status"] = "removed"
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
            if "x" in daughter:
                angle = seeded_rng(seed, "daughter", step, cid).random() * math.tau
                daughter["x"] += math.cos(angle) * 2
                daughter["y"] += math.sin(angle) * 2
            next_cells[daughter_id] = daughter
            living += 1
            events.append(dict(kind="division", cell_id=cid, daughter=daughter_id))
        elif action == "ACTIVATE":
            c["activated"] = True
        elif action == "EXHAUST":
            c["exhausted"] = True
        elif action == "RECOVER":
            c["exhausted"] = False
        elif action == "SECRETE":
            c["signal"] += 1.0
    return next_cells, events


def preflight(bundle, *, model="jev-1.14.0", max_questions=20, seed=21, prior_mode="shown"):
    validate_prepared(bundle)
    cells = initialize(bundle)
    rows, graph = decision_rows(cells, bundle["guidebook"], bundle["context"], seed=seed, step=0)
    active_rows = [r for r in rows if len(r["actions"]) > 1]
    groups = batches(active_rows, bundle["context"], model, max_questions, prior_mode=prior_mode)
    return dict(schema_version="2.0", cells=len(cells), by_lineage=dict(Counter(c["lineage"] for c in cells.values())),
        initial_questions=len(active_rows), initial_http_batches=len(groups),
        isolated_cells=sum(not ns for ns in graph.values()), interaction_mode=bundle["context"].get("interaction_mode"),
        limitations=["Preferences are uncalibrated; simulation hours are declared step labels, not fitted biological time.",
                    "t0 expression stays fixed; state transitions do not predict a new transcriptome.",
                    "Spatial edges use centroid distance; they do not establish physical contact.",
                    "Enrichment/subsampling does not reconstruct absent populations or population frequencies."])


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
        max_population=2000, resume=False):
    validate_prepared(bundle)
    if type(steps) is not int or not 1 <= steps <= 10000 or type(seed) is not int:
        raise ValueError("Use 1–10000 steps and an integer seed.")
    if mode not in ("local", "jev") or (mode == "jev" and provider is None):
        raise ValueError("mode must be local, or jev with a provider.")
    if type(jev_weight) not in (int, float) or not math.isfinite(jev_weight) or not 0 <= jev_weight <= 1:
        raise ValueError("jev_weight must be in [0,1].")
    if type(max_population) is not int or max_population < len(bundle["cells"]):
        raise ValueError("max_population must cover all starting cells.")
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    config = dict(schema_version="2.0", bundle_fingerprint=bundle["fingerprint"], steps=steps, seed=seed,
        mode=mode, jev_weight=jev_weight, max_population=max_population,
        model=provider.model if provider else None,
        provider_name=getattr(provider, "name", None), provider_prompt_version=getattr(provider, "prompt_version", None),
        jev_prior_mode=getattr(provider, "prior_mode", None),
        max_questions=getattr(provider, "max_questions", None),
        probability_semantics="guidebook prior (local) or normalized prior^(1-alpha)*Jev^alpha; seeded categorical sample",
        kernel="vtx-discrete-context-2.0", assumptions=dict(signal_retention_per_step=.5, signal_increment=1.0,
            signal_activation_multiplier="1 + mean neighbor signal from previous snapshot",
            kill="one licensed KILL commits one target death per step; target = nearest live malignant neighbor (spatial) or seeded draw (well_mixed)",
            clearance="target = nearest corpse (spatial) or seeded draw (well_mixed)",
            division="one daughter per selected action per step",
            geometry="fixed centroids; daughter offset 2 um; no mechanics; capped nearest-neighbor graph",
            simultaneous="snapshot proposals; death/clearance conflicts first; surviving actors divide/change state"))
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
    cells = initialize(bundle)
    frames = [dict(step=0, time_hours=0, cells=deepcopy(cells), events=[])]
    decisions, stop, error, divergences = [], "completed", None, {}
    for step in range(steps):
        rows, graph = decision_rows(cells, bundle["guidebook"], bundle["context"], seed=seed, step=step)
        if not rows:
            stop = "no_living_cells"
            break
        active_rows = [r for r in rows if len(r["actions"]) > 1]
        try:
            answers = provider.choose(active_rows, bundle["context"]) if mode == "jev" else {}
            if mode == "jev" and set(answers) != {r["cell_id"] for r in active_rows}:
                raise ProviderError("Provider omitted or added cell answers; no state committed.")
        except (ProviderError, BudgetExhausted) as e:
            stop, error = "provider_stopped_before_round_commit", str(e)
            break
        proposals = []
        for row in rows:
            cid = row["cell_id"]
            answer = answers.get(cid)
            probabilities = distribution(row["prior"], answer, jev_weight)
            action = draw(probabilities, seeded_rng(seed, "decision", step, cid))
            divergence = js_divergence(row["prior"], answer["probabilities"]) if answer else None
            if divergence is not None:
                divergences.setdefault(row["observation"]["lineage"], []).append(divergence)
            proposals.append(dict(cell_id=cid, action=action, targets=row["targets"], corpses=row["corpses"],
                                  target_rule=row["target_rule"]))
            decisions.append(dict(step=step, time_hours=step*bundle["context"]["time_step_hours"], cell_id=cid,
                action=action, prior=row["prior"], jev=answer, applied_probabilities=probabilities,
                jev_prior_js_bits=divergence,
                observation=row["observation"], excluded=row["excluded"],
                provenance="local_guidebook" if mode == "local" else ("jev_pooled_preferences" if answer else "single_permitted_action")))
        cells, events = commit(cells, proposals, max_population=max_population, seed=seed, step=step)
        frames.append(dict(step=step+1, time_hours=(step+1)*bundle["context"]["time_step_hours"], cells=deepcopy(cells), events=events))
        # Checkpoint after each complete round; never expose half-applied decisions.
        write_json(directory / "checkpoint.json", dict(config_hash=identity, completed_steps=step+1, cells=cells))
    report = dict(schema_version="2.0", completed=len(frames)-1 == steps, stop_reason=stop, error=error,
        completed_steps=len(frames)-1, requested_steps=steps, seed=seed, mode=mode,
        living=sum(c["status"] == "alive" for c in cells.values()), dead=sum(c["status"] == "dead" for c in cells.values()),
        removed=sum(c["status"] == "removed" for c in cells.values()), allocated=len(cells),
        events=dict(Counter(e["kind"] for f in frames for e in f["events"])),
        bundle_fingerprint=bundle["fingerprint"], config_hash=identity,
        usage=provider.usage() if provider else dict(http_attempts=0),
        jev_prior_mode=getattr(provider, "prior_mode", None),
        jev_prior_divergence=dict(unit="Jensen-Shannon divergence, bits, Jev vs gated guidebook prior",
            by_lineage={lineage: dict(n=len(v), mean=sum(v)/len(v), max=max(v)) for lineage, v in sorted(divergences.items())},
            interpretation="Near zero with prior shown and larger with prior hidden indicates anchoring; compare against replicate variance."),
        biological_validation="not established", preflight=preflight(bundle, seed=seed,
            model=getattr(provider, "model", "jev-1.14.0"), max_questions=getattr(provider, "max_questions", 20),
            prior_mode=getattr(provider, "prior_mode", "shown")))
    write_json(directory / "frames.json", frames)
    with (directory / "decisions.jsonl").open("w", encoding="utf8") as f:
        for d in decisions:
            f.write(json.dumps(d, allow_nan=False) + "\n")
    write_json(directory / "report.json", report)
    from .replay import write_replay
    write_replay(directory / "replay.html", frames, report)
    return report
