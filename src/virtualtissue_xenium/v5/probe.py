"""Frozen-state provider probes: what moves Jev's answer, beyond replicate noise?

A simulated trajectory cannot attribute Jev's behaviour: after one step two arms
visit different states and menus. A probe asks about one frozen t0 state instead.
The same cells, the same permitted menus and the same request batches are sent
under several arms, each repeated into independent journals:

- ``baseline``: the questions exactly as the simulation asks them at t0;
- ``evidence_masked``: every t0 evidence value withheld (with a policy, its
  definition too);
- ``evidence_exchanged``: a cell receives the evidence of another cell of the
  same declared type;
- ``neighbourhood_exchanged``: a cell receives the local summary of another cell
  of the same lineage and menu;
- ``context_exchanged``: tissue, disease and perturbation text from a donor
  context (needs ``donor_context``).

``plan`` builds the question bank and the cost estimate without any request.
``run_probe`` is the only step that sends requests; it refuses a budget that
cannot cover the whole plan and never resends a request. ``analyze`` reads the
journals offline and reports divergence between arms next to divergence between
identical repeats; a request group with a failed request is left out of every arm.

With the prior shown, the guidebook prior stays the baseline prior in every arm,
so an evidence arm then isolates the evidence shown directly to the model. No
threshold here certifies a result; the numbers describe one model version on one
question bank.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from itertools import combinations
from pathlib import Path
import random
from statistics import mean

from ..v2.inputs import canonical, digest, read_json, validate_prepared, write_json
from ..v2.provider import (TOKENS_PER_BYTE, JevProvider, ProviderError, payload, prompt_version, validate_response)
from ..v2.simulation import decision_rows, initialize, js_divergence
from .evidence import declared_rule

PROBE_VERSION = "jev-probe-1"
ARMS = {
    "baseline": "Questions exactly as the simulation asks them at t0.",
    "evidence_masked": "Every t0 evidence value withheld (status masked_by_probe), together with its definition when an "
                       "evidence policy is in use; state, neighbourhood and menu unchanged.",
    "evidence_exchanged": "Each cell receives the t0 evidence of another cell of the same declared type; everything else unchanged.",
    "neighbourhood_exchanged": "Each cell receives the local neighbourhood summary of another cell of the same lineage and menu.",
    "context_exchanged": "Tissue, disease and perturbation text taken from a declared donor context; numeric settings unchanged.",
}
MAX_BYTES = 40000


def _order(seed, purpose, cell_id):
    return digest([seed, PROBE_VERSION, purpose, cell_id])


def _donor(row, candidates, seed, purpose):
    """Deterministic other cell from ``candidates``; None when the cell has no counterpart."""
    others = [c for c in candidates if c["cell_id"] != row["cell_id"]]
    return min(others, key=lambda c: _order(seed, purpose, [row["cell_id"], c["cell_id"]])) if others else None


def _arm_rows(arm, bank, active, seed):
    """Question rows for one arm, plus {cell_id: donor cell_id or None} for exchange arms."""
    rows, exchanges = deepcopy(bank), {}
    if arm == "evidence_masked":
        for row in rows:
            row["observation"]["t0_evidence"] = {
                name: dict(spec, z=None, status="masked_by_probe", reference=None)
                for name, spec in row["observation"]["t0_evidence"].items()}
    elif arm == "evidence_exchanged":
        for row in rows:
            same = [c for c in active if c["observation"]["cell_type"] == row["observation"]["cell_type"]]
            donor = _donor(row, same, seed, arm)
            exchanges[row["cell_id"]] = donor["cell_id"] if donor else None
            if donor:
                row["observation"]["t0_evidence"] = deepcopy(donor["observation"]["t0_evidence"])
    elif arm == "neighbourhood_exchanged":
        for row in rows:
            same = [c for c in active if c["observation"]["lineage"] == row["observation"]["lineage"]
                    and set(c["actions"]) == set(row["actions"])]
            donor = _donor(row, same, seed, arm)
            exchanges[row["cell_id"]] = donor["cell_id"] if donor else None
            if donor:
                row["observation"]["local"] = deepcopy(donor["observation"]["local"])
    return rows, exchanges


def _arm_context(arm, context, donor_context):
    if arm != "context_exchanged":
        return context
    from ..v2.context import validate_context
    donor = validate_context(donor_context)
    if donor["organism"].casefold() != context["organism"].casefold():
        raise ValueError("The donor context must preserve the organism.")
    changed = deepcopy(context)
    for field in ("tissue", "disease", "perturbations"):
        changed[field] = deepcopy(donor[field])
    changed.pop("sample_context", None)
    if "sample_context" in donor:
        changed["sample_context"] = deepcopy(donor["sample_context"])
    if all(changed.get(field) == context.get(field) for field in ("tissue", "disease", "perturbations", "sample_context")):
        raise ValueError("The donor context has the same tissue, disease, perturbations and sample details; the arm would change nothing.")
    return changed


def _build(bundle, settings, donor_context):
    """Deterministic question bank, arm rows, contexts and fixed groups from a plan's settings."""
    validate_prepared(bundle)
    if bundle["schema_version"] != "3.0":
        raise ValueError("Provider probes require a schema 3.0 prepared bundle.")
    seed = settings["seed"]
    rows, _ = decision_rows(initialize(bundle, seed=seed), bundle["guidebook"], bundle["context"], seed=seed, step=0)
    active = [r for r in rows if len(r["actions"]) > 1]
    by_lineage = defaultdict(list)
    for row in active:
        by_lineage[row["observation"]["lineage"]].append(row)
    bank = sorted((row for members in by_lineage.values()
                   for row in sorted(members, key=lambda r: _order(seed, "sample", r["cell_id"]))[:settings["cells_per_lineage"]]),
                  key=lambda r: r["cell_id"])
    if not bank:
        raise ValueError("No cell has more than one permitted action at t0; there is nothing to ask.")
    arms, contexts, exchanges = {}, {}, {}
    for arm in settings["arms"]:
        arms[arm], changed = _arm_rows(arm, bank, active, seed)
        contexts[arm] = _arm_context(arm, bundle["context"], donor_context)
        if changed:
            exchanges[arm] = changed
    policy = bundle["context"].get("evidence_policy") is not None

    def question(arm, index, prior):
        return payload([arms[arm][index]], contexts[arm], settings["model"], prior, settings["pooling"])

    # A donor can carry the same visible content (for example when no evidence is
    # usable on the panel). The request is then the baseline request: nothing was
    # manipulated, so the cell must not count as evidence that the arm had no effect.
    for arm, changed in exchanges.items():
        for index, row in enumerate(bank):
            if changed[row["cell_id"]] is not None and all(
                    question(arm, index, prior) == question("baseline", index, prior) for prior in settings["prior_modes"]):
                changed[row["cell_id"]] = None

    def fits(indices):
        return len(indices) <= settings["max_questions"] and all(
            len(canonical(payload([arms[arm][i] for i in indices], contexts[arm], settings["model"], prior, settings["pooling"]))) <= MAX_BYTES
            for arm in settings["arms"] for prior in settings["prior_modes"])

    groups, current = [], []
    for index in range(len(bank)):
        if not fits([index]):
            raise ValueError("A single probe question exceeds the provider request limit.")
        if current and not fits(current + [index]):
            groups.append(current)
            current = []
        current.append(index)
    groups.append(current)
    return bank, arms, contexts, exchanges, groups, policy


def _settings(bundle, *, cells_per_lineage, seed, arms, prior_modes, repeats, model, max_questions, pooling, donor_context):
    arms = list(dict.fromkeys(arms))
    if not arms or arms[0] != "baseline" or set(arms) - set(ARMS):
        raise ValueError(f"Arms must start with baseline and be chosen from {sorted(ARMS)}.")
    if ("context_exchanged" in arms) != (donor_context is not None):
        raise ValueError("context_exchanged requires a donor context, and a donor context is used only by that arm.")
    prior_modes = list(dict.fromkeys(prior_modes))
    if not prior_modes or set(prior_modes) - {"shown", "hidden"}:
        raise ValueError("prior_modes must be shown and/or hidden.")
    for value, name, low, high in ((cells_per_lineage, "cells_per_lineage", 1, 200), (repeats, "repeats", 1, 20),
                                   (max_questions, "max_questions", 1, 20)):
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f"{name} must be an integer from {low} to {high}.")
    if type(seed) is not int:
        raise ValueError("seed must be an integer.")
    return dict(bundle_fingerprint=bundle["fingerprint"], cells_per_lineage=cells_per_lineage, seed=seed, arms=arms,
                prior_modes=prior_modes, repeats=repeats, model=model, max_questions=max_questions, pooling=pooling,
                donor_context_sha256=digest(donor_context) if donor_context is not None else None)


def plan(bundle, output_dir, *, cells_per_lineage=6, seed=0, arms=("baseline", "evidence_masked", "evidence_exchanged",
         "neighbourhood_exchanged"), prior_modes=("shown", "hidden"), repeats=3, donor_context=None,
         model="jev-1.14.0", max_questions=10, pooling="full_distribution", prepared_path=None):
    """Write plan.json and payloads.json; send nothing."""
    settings = _settings(bundle, cells_per_lineage=cells_per_lineage, seed=seed, arms=arms, prior_modes=prior_modes,
                         repeats=repeats, model=model, max_questions=max_questions, pooling=pooling, donor_context=donor_context)
    directory = Path(output_dir)
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Probe directory is not empty; choose a new one.")
    bank, arm_rows, contexts, exchanges, groups, policy = _build(bundle, settings, donor_context)
    requests, bodies = [], {}
    for arm in settings["arms"]:
        for prior in settings["prior_modes"]:
            for index, group in enumerate(groups):
                body = payload([arm_rows[arm][i] for i in group], contexts[arm], model, prior, pooling)
                bodies.setdefault(arm, {}).setdefault(prior, []).append(body)
                requests.append(dict(arm=arm, prior=prior, group=index, request_sha256=digest(body), bytes=len(canonical(body))))
    planned = len(requests) * repeats
    total_bytes = sum(r["bytes"] for r in requests) * repeats
    result = dict(
        probe_version=PROBE_VERSION, settings=settings, prepared_path=str(Path(prepared_path).resolve()) if prepared_path else None,
        arms={arm: ARMS[arm] for arm in settings["arms"]},
        prompt_versions={prior: prompt_version(prior, "3.0", evidence_policy=policy, pooling=pooling) for prior in settings["prior_modes"]},
        cells=[dict(cell_id=row["cell_id"], lineage=row["observation"]["lineage"], cell_type=row["observation"]["cell_type"],
                    actions=sorted(row["actions"])) for row in bank],
        groups=[[bank[i]["cell_id"] for i in group] for group in groups], exchanges=exchanges,
        cells_unchanged_by_design={arm: sum(donor is None for donor in changed.values()) for arm, changed in exchanges.items()},
        requests=requests,
        estimate=dict(http_requests=planned, questions=len(bank) * len(settings["arms"]) * len(settings["prior_modes"]) * repeats,
                      input_bytes=total_bytes, input_tokens_estimate=int(round(total_bytes * TOKENS_PER_BYTE)),
                      tokens_per_byte_assumption=TOKENS_PER_BYTE, note="Output tokens are not estimated; this is a plan, not a tariff."),
        design="One frozen t0 state; identical cells, menus and batches in every arm; repeats go to independent journals.",
        limits=["With the prior shown, the guidebook prior is the baseline prior in every arm.",
                "An exchange arm leaves a cell unchanged when it has no counterpart or when its counterpart's visible "
                "content is identical; those cells (exchanges[arm][cell] = null) are excluded from that arm's statistics.",
                "With an evidence policy, evidence_masked withholds the values together with their definitions.",
                "context_exchanged also changes the opaque full-context hash carried in the request.",
                "A failed request is never resent; its request group is then excluded from the analysis in every arm.",
                "A probe describes one model version and one question bank at t0; it is not a trajectory result."])
    result["fingerprint"] = digest(result)
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory / "plan.json", result)
    write_json(directory / "payloads.json", bodies)
    if donor_context is not None:
        write_json(directory / "donor_context.json", donor_context)
    return result


def _load(directory, bundle=None):
    directory = Path(directory)
    plan_record = read_json(directory / "plan.json")
    content = {k: v for k, v in plan_record.items() if k != "fingerprint"}
    if plan_record.get("probe_version") != PROBE_VERSION or digest(content) != plan_record.get("fingerprint"):
        raise ValueError("plan.json is not an unmodified probe plan.")
    if bundle is None:
        if not plan_record.get("prepared_path") or not Path(plan_record["prepared_path"]).is_file():
            raise ValueError("The prepared bundle recorded in the plan is unavailable; pass it explicitly.")
        bundle = read_json(plan_record["prepared_path"])
    if not isinstance(bundle, dict) or bundle.get("fingerprint") != plan_record["settings"]["bundle_fingerprint"]:
        raise ValueError("The prepared bundle differs from the one the probe was planned on.")
    validate_prepared(bundle)  # the fingerprint field alone is not proof of the content
    donor_path = directory / "donor_context.json"
    donor = read_json(donor_path) if donor_path.exists() else None
    return directory, plan_record, bundle, donor


def _journal(directory, plan_record, repeat, arm, prior):
    """Journal folder of one repeat/arm/prior and the identity its records must carry."""
    identity = digest(dict(run_identity=dict(probe=plan_record["fingerprint"], arm=arm, prior=prior, repeat=repeat),
                           model=plan_record["settings"]["model"], prompt=plan_record["prompt_versions"][prior]))
    return directory / "journals" / f"repeat_{repeat}" / arm / prior, identity


def run_probe(output_dir, key, *, max_http_attempts, bundle=None, transport=None, workers=1):
    """Send the planned requests; a budget below the plan sends nothing.

    Completed requests are never sent again, so an interrupted run can be
    continued. A request that failed is not resent either (its usage is
    unknown): it is reported, the remaining requests are still sent, and
    ``analyze`` leaves that request group out of every arm.
    """
    directory, plan_record, bundle, donor = _load(output_dir, bundle)
    settings = plan_record["settings"]
    planned = plan_record["estimate"]["http_requests"]
    if type(max_http_attempts) is not int or max_http_attempts < planned:
        raise ValueError(f"The plan needs {planned} HTTP requests; --max-http-attempts {max_http_attempts} cannot cover it. Nothing sent.")
    bank, arm_rows, contexts, _, groups, policy = _build(bundle, settings, donor)
    expected = {(r["arm"], r["prior"], r["group"]): r["request_sha256"] for r in plan_record["requests"]}
    for arm in settings["arms"]:
        for prior in settings["prior_modes"]:
            for index, group in enumerate(groups):
                body = payload([arm_rows[arm][i] for i in group], contexts[arm], settings["model"], prior, settings["pooling"])
                if digest(body) != expected[(arm, prior, index)]:
                    raise ValueError("Rebuilt requests differ from the plan; the package or bundle changed since planning.")
    usage = dict(http_attempts=0, cache_hits=0, input_tokens=0, output_tokens=0, attempts_with_unknown_usage=0)
    failures = []
    for repeat in range(settings["repeats"]):
        for arm in settings["arms"]:
            for prior in settings["prior_modes"]:
                folder, _ = _journal(directory, plan_record, repeat, arm, prior)
                provider = JevProvider(key, folder, dict(probe=plan_record["fingerprint"], arm=arm, prior=prior, repeat=repeat),
                                       model=settings["model"], max_attempts=len(groups), workers=workers,
                                       max_questions=settings["max_questions"], transport=transport, prior_mode=prior,
                                       schema_version="3.0", evidence_policy=policy, pooling=settings["pooling"])
                for index, group in enumerate(groups):
                    try:
                        provider.choose([arm_rows[arm][i] for i in group], contexts[arm])
                    except ProviderError as error:
                        failures.append(dict(repeat=repeat, arm=arm, prior=prior, group=index, error=str(error)))
                for name, value in provider.usage().items():
                    usage[name] += value
    result = dict(completed=not failures, planned_http_requests=planned, failed_requests=failures, usage=usage,
                  journals=str((directory / "journals").resolve()),
                  note=("Failed requests are not resent; analyze excludes their request groups from every arm."
                        if failures else None))
    write_json(directory / "run.json", result)
    return result


def _answers(directory, plan_record):
    """Validated answers and the request groups that cannot be used.

    Returns ({(repeat, arm, prior): {cell_id: probabilities}}, {group index: reason}).
    A group with a missing or failed result in any arm, prior or repeat is left
    out everywhere, so every arm is still compared on the same cells. A result
    that fails its checksum or belongs to another journal is an error.
    """
    settings, tables, unusable, records = plan_record["settings"], {}, {}, {}
    payloads = read_json(directory / "payloads.json")
    for repeat in range(settings["repeats"]):
        for arm in settings["arms"]:
            for prior in settings["prior_modes"]:
                journal, identity = _journal(directory, plan_record, repeat, arm, prior)
                marker = journal / "identity.json"
                if marker.exists() and read_json(marker) != {"identity": identity}:
                    raise ValueError(f"The journal for {arm}/{prior}/repeat {repeat} belongs to another probe or repeat.")
                for index, body in enumerate(payloads[arm][prior]):
                    request = digest(body)
                    path = journal / f"{request}.result.json"
                    if not path.exists():
                        unusable.setdefault(index, f"no result for {arm}/{prior}/repeat {repeat}")
                        continue
                    record = read_json(path)
                    checksum = record.pop("digest", None)
                    if checksum != digest(record) or record.get("request_hash") != request or record.get("identity") != identity:
                        raise ValueError(f"A journal result for {arm}/{prior}/repeat {repeat} was edited or copied from another journal.")
                    if record.get("status") != "valid":
                        unusable.setdefault(index, f"failed request for {arm}/{prior}/repeat {repeat}")
                        continue
                    records[(repeat, arm, prior, index)] = validate_response(record["response"], body)
    for (repeat, arm, prior, index), answers in records.items():
        if index in unusable:
            continue
        table = tables.setdefault((repeat, arm, prior), {})
        for position, cell_id in enumerate(plan_record["groups"][index]):
            table[cell_id] = answers[f"q{position}"]["probabilities"]
    if len(unusable) == len(plan_record["groups"]):
        raise ValueError("No request group has a valid result in every arm, prior mode and repeat; run the probe first.")
    return tables, unusable


def _mean_distribution(distributions):
    return {a: mean(d[a] for d in distributions) for a in distributions[0]}


def _modal(distribution):
    return max(sorted(distribution), key=lambda a: distribution[a])


def _bootstrap(clusters, replicates, rng):
    """Percentile interval of the mean, resampling request groups rather than cells.

    Questions sent in one request can share request-level variation, so the
    request group is the resampling unit. Fewer than two groups give no interval.
    """
    clusters = [values for values in clusters if values]
    if len(clusters) < 2:
        return None
    means = sorted(mean(value for cluster in rng.choices(clusters, k=len(clusters)) for value in cluster)
                   for _ in range(replicates))
    return [means[int(.025 * (replicates - 1))], means[int(.975 * (replicates - 1))]]


def _ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks, i = [0.0] * len(values), 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def spearman(x, y):
    """Rank correlation with average ranks for ties; None when either side is constant."""
    if len(x) < 3:
        return None
    rx, ry = _ranks(x), _ranks(y)
    mx, my = mean(rx), mean(ry)
    sx, sy = sum((a - mx) ** 2 for a in rx), sum((b - my) ** 2 for b in ry)
    if sx == 0 or sy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(rx, ry)) / (sx * sy) ** .5


def analyze(output_dir, *, bundle=None, bootstrap_replicates=2000, bootstrap_seed=0, output=None):
    """Compare arms on the frozen question bank; offline. Refuses to overwrite an analysis."""
    if type(bootstrap_replicates) is not int or not 100 <= bootstrap_replicates <= 100000:
        raise ValueError("bootstrap_replicates must be an integer from 100 to 100000.")
    directory, plan_record, bundle, _ = _load(output_dir, bundle)
    output = Path(output) if output else directory / "analysis.json"
    if output.exists():
        raise ValueError("Probe analysis exists; choose a new --output path.")
    settings = plan_record["settings"]
    answers, unusable = _answers(directory, plan_record)
    usable_groups = [group for index, group in enumerate(plan_record["groups"]) if index not in unusable]
    cell_ids = [cell_id for group in usable_groups for cell_id in group]
    repeats, lineage_of = range(settings["repeats"]), {c["cell_id"]: c["lineage"] for c in plan_record["cells"]}
    pairs_of_repeats = list(combinations(repeats, 2))
    rng = random.Random(int(digest([PROBE_VERSION, "bootstrap", bootstrap_seed]), 16))
    result = dict(probe_version=PROBE_VERSION, plan_fingerprint=plan_record["fingerprint"], model=settings["model"],
                  prompt_versions=plan_record["prompt_versions"], n_cells=len(cell_ids), n_cells_planned=len(plan_record["cells"]),
                  repeats=settings["repeats"], request_groups_used=len(usable_groups),
                  request_groups_excluded=[dict(group=index, cells=len(plan_record["groups"][index]), reason=reason)
                                           for index, reason in sorted(unusable.items())],
                  by_prior={})
    rows, _ = decision_rows(initialize(bundle, seed=settings["seed"]), bundle["guidebook"], bundle["context"],
                            seed=settings["seed"], step=0)
    row_of = {row["cell_id"]: row for row in rows}
    policy = bundle["context"].get("evidence_policy")
    shown = policy["provider_evidence"] if policy else "all"
    for prior in settings["prior_modes"]:
        def within(arm, cell_id):
            values = [js_divergence(answers[(i, arm, prior)][cell_id], answers[(j, arm, prior)][cell_id])
                      for i, j in pairs_of_repeats]
            return mean(values) if values else None

        def modal_agreement(arm, cell_id):
            """Share of repeat pairs with the same most-preferred action; within baseline, or arm against baseline."""
            if arm == "baseline":
                pairs = [(answers[(i, arm, prior)][cell_id], answers[(j, arm, prior)][cell_id]) for i, j in pairs_of_repeats]
            else:
                pairs = [(answers[(i, arm, prior)][cell_id], answers[(j, "baseline", prior)][cell_id])
                         for i in repeats for j in repeats]
            return mean(_modal(a) == _modal(b) for a, b in pairs) if pairs else None

        base_mean = {cid: _mean_distribution([answers[(r, "baseline", prior)][cid] for r in repeats]) for cid in cell_ids}
        noise = [v for v in (within("baseline", cid) for cid in cell_ids) if v is not None]
        section = dict(
            replicate_noise=dict(
                unit="Jensen-Shannon divergence, bits, between identical baseline requests in independent journals",
                mean=mean(noise) if noise else None, max=max(noise) if noise else None, n_cells=len(noise),
                modal_action_agreement=(mean(modal_agreement("baseline", cid) for cid in cell_ids) if pairs_of_repeats else None),
                note=None if noise else "One repeat only: replicate noise is unmeasured and arm contrasts cannot be judged against it."),
            arms={})
        for arm in settings["arms"][1:]:
            exchanged = plan_record["exchanges"].get(arm)
            included = [cid for cid in cell_ids if exchanged is None or exchanged.get(cid) is not None]
            between, by_lineage, excess_of = [], defaultdict(list), {}
            for cid in included:
                value = mean(js_divergence(answers[(i, arm, prior)][cid], answers[(j, "baseline", prior)][cid])
                             for i in repeats for j in repeats)
                between.append(value)
                by_lineage[lineage_of[cid]].append(value)
                inside = [v for v in (within("baseline", cid), within(arm, cid)) if v is not None]
                if inside:
                    excess_of[cid] = value - mean(inside)
            excess = list(excess_of.values())
            arm_mean = {cid: _mean_distribution([answers[(r, arm, prior)][cid] for r in repeats]) for cid in included}
            section["arms"][arm] = dict(
                description=ARMS[arm], n_cells=len(included), n_cells_unchanged_by_design=len(cell_ids) - len(included),
                js_between_arm_and_baseline=dict(mean=mean(between) if between else None, max=max(between) if between else None,
                    by_lineage={k: dict(n=len(v), mean=mean(v)) for k, v in sorted(by_lineage.items())}),
                js_excess_over_replicates=dict(
                    definition="per cell: mean JS(arm repeat, baseline repeat) minus the average of the mean JS between "
                               "baseline repeats and the mean JS between repeats of the arm",
                    mean=mean(excess) if excess else None,
                    bootstrap_ci95_over_request_groups=_bootstrap(
                        [[excess_of[cid] for cid in group if cid in excess_of] for group in usable_groups],
                        bootstrap_replicates, rng) if excess else None,
                    fraction_of_cells_above_zero=mean(v > 0 for v in excess) if excess else None),
                modal_action_agreement_with_baseline=(mean(modal_agreement(arm, cid) for cid in included) if included else None),
                event_mass=mean(1 - arm_mean[cid]["REST"] for cid in included) if included else None,
                baseline_event_mass_same_cells=mean(1 - base_mean[cid]["REST"] for cid in included) if included else None)
        # Does the stated preference follow the evidence the rulebook declares for that action?
        tracking = []
        pairs = sorted({(lineage_of[cid], action) for cid in cell_ids for action in base_mean[cid] if action != "REST"})
        for lineage, action in pairs:
            rule = declared_rule(bundle["context"], lineage, action)
            if rule is None:
                continue
            name, sign = rule["evidence"], (-1.0 if rule["direction"] == "opposes" else 1.0)
            zs, shares, priors = [], [], []
            for cid in cell_ids:
                spec = row_of[cid]["observation"]["t0_evidence"].get(name, {})
                if lineage_of[cid] != lineage or action not in base_mean[cid] or spec.get("z") is None:
                    continue
                event = 1 - base_mean[cid]["REST"]
                prior_event = 1 - row_of[cid]["prior"]["REST"]
                if event <= 0 or prior_event <= 0:
                    continue
                zs.append(sign * spec["z"])
                shares.append(base_mean[cid][action] / event)
                priors.append(row_of[cid]["prior"][action] / prior_event)
            if len(zs) >= 3:
                tracking.append(dict(lineage=lineage, action=action, evidence=name, kind=rule["kind"],
                                     direction=rule["direction"], n_cells=len(zs),
                                     evidence_shown_to_provider=shown == "all" or name in shown,
                                     spearman_jev_share_vs_z=spearman(zs, shares),
                                     spearman_guidebook_share_vs_z=spearman(zs, priors)))
        section["evidence_tracking"] = dict(
            definition="Across probed cells of a lineage: rank correlation between the declared evidence z, signed so that "
                       "positive means agreement with the rule's direction, and the share of event preference given to the "
                       "action, P(action)/(1-P(REST)), in the baseline arm. The guidebook column is the same quantity for "
                       "the local prior; it is expected to be positive when beta is nonzero, although menus, gates and "
                       "other multipliers also differ between cells.",
            rows=tracking,
            note="Small n; descriptive only. With the prior shown, agreement can come from the prior rather than the "
                 "evidence. When evidence_shown_to_provider is false, the model never saw that item.")
        result["by_prior"][prior] = section
    result["interpretation"] = ("Compare each arm's divergence from baseline with the divergence between identical repeats. "
                                "An excess near zero means the manipulation did not move the answers beyond replicate noise "
                                "on this question bank. Nothing here validates a biological prediction.")
    result["fingerprint"] = digest(result)
    write_json(output, result)
    return result
