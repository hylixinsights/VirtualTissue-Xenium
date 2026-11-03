"""Auditable seed × arm execution with separate uncertainty and seed envelopes.

No inferential unit is a simulated cell. Failed or truncated runs are retained
in the manifest and never silently enter aggregates. There is no automatic
provider retry or automatic continuation of an existing ensemble directory.
"""
from __future__ import annotations

from collections import defaultdict
import csv
from copy import deepcopy
import json
import math
from pathlib import Path
import random
import re
from statistics import median

from .inputs import digest, read_json, validate_prepared, write_json
from .readouts import extract_readouts, load_readout_spec
from .simulation import run


def quantile(values, q):
    """Linear interpolation between ordered observations (type 7)."""
    ordered = sorted(values)
    if not ordered:
        return None
    position = (len(ordered) - 1) * q
    lower, upper = math.floor(position), math.ceil(position)
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def summarize_values(values, *, bootstrap_replicates=2000, seed=0):
    """Percentile bootstrap CI for the median, separate from seed quantiles."""
    if type(bootstrap_replicates) is not int or bootstrap_replicates < 100:
        raise ValueError("Use at least 100 bootstrap replicates.")
    finite = [float(v) for v in values if v is not None]
    if any(not math.isfinite(v) for v in finite):
        raise ValueError("Readouts must be finite numbers or null.")
    n = len(finite)
    result = dict(n=n, n_missing=len(values)-n, median=median(finite) if n else None,
                  seed_q025=quantile(finite, .025), seed_q975=quantile(finite, .975),
                  median_ci95_low=None, median_ci95_high=None,
                  interval_status="insufficient_seeds" if n < 2 else "ok")
    if n >= 2:
        rng = random.Random(seed)
        medians = [median(rng.choices(finite, k=n)) for _ in range(bootstrap_replicates)]
        result.update(median_ci95_low=quantile(medians, .025), median_ci95_high=quantile(medians, .975))
    return result


def standardized_paired_effect(differences, *, bootstrap_replicates=2000, seed=0):
    """Cohen dz across paired simulation seeds, not cells or patients.

The denominator is the sample SD (ddof=1) of arm-minus-baseline differences.
Resampling whole differences preserves pairing. Undefined bootstrap ratios are
counted, never dropped to manufacture a conditional confidence interval.
"""
    if type(bootstrap_replicates) is not int or bootstrap_replicates < 100:
        raise ValueError("Use at least 100 bootstrap replicates.")
    finite = [float(v) for v in differences if v is not None]
    if any(not math.isfinite(v) for v in finite):
        raise ValueError("Paired differences must be finite numbers or null.")
    n = len(finite)
    result = dict(paired_mean_difference=math.fsum(finite)/n if n else None,
        paired_difference_sd=None, paired_dz=None, paired_dz_ci95_low=None, paired_dz_ci95_high=None,
        paired_dz_n=n, paired_dz_n_missing=len(differences)-n,
        paired_dz_status="insufficient_seeds", paired_dz_ci_status="insufficient_seeds",
        paired_dz_bootstrap_valid=0, paired_dz_bootstrap_undefined=0)
    if n < 2:
        return result

    def estimate(sample):
        mean = math.fsum(sample)/n
        # Scaling avoids overflow when squaring large finite differences.
        scale = max(abs(value-mean) for value in sample)
        sd = scale * math.sqrt(math.fsum(((value-mean)/scale)**2 for value in sample)/(n-1)) if scale else 0.0
        dz = mean/sd if sd else None
        return sd, dz if dz is not None and math.isfinite(dz) else None

    sd, dz = estimate(finite)
    result.update(paired_difference_sd=sd, paired_dz=dz)
    if dz is None:
        result.update(paired_dz_status="undefined_sd_or_ratio", paired_dz_ci_status="undefined_original_sd_or_ratio")
        return result
    rng = random.Random(seed)
    boot = [estimate(rng.choices(finite, k=n))[1] for _ in range(bootstrap_replicates)]
    valid = [value for value in boot if value is not None]
    result.update(paired_dz_status="ok", paired_dz_bootstrap_valid=len(valid),
                  paired_dz_bootstrap_undefined=bootstrap_replicates-len(valid))
    if len(valid) == bootstrap_replicates:
        result.update(paired_dz_ci95_low=quantile(valid, .025), paired_dz_ci95_high=quantile(valid, .975),
                      paired_dz_ci_status="ok")
    else:
        result["paired_dz_ci_status"] = "undefined_bootstrap_ratios"
    return result


def _key(row):
    return (round(row["time_hours"], 12), row["readout"], row["lineage"])


def aggregate_readouts(run_records, *, bootstrap_replicates=2000, bootstrap_seed=0, baseline=None):
    """Aggregate only complete matched seed sets, keeping absent values null.

Any missing/failed seed disqualifies its arm from aggregation. Paired contrasts
are computed only for the same seeds and exact common time points.
"""
    arm_runs = defaultdict(list)
    for record in run_records:
        arm_runs[record["arm"]].append(record)
    eligible = {arm: runs for arm, runs in arm_runs.items() if all(r["completed"] for r in runs)}
    summaries, contrasts, indexes = [], [], {}
    for arm, runs in sorted(eligible.items()):
        if len({r["seed"] for r in runs}) != len(runs):
            raise ValueError("Duplicate seed within an arm.")
        indexes[arm] = {r["seed"]: {_key(row): row for row in r["readouts"]} for r in runs}
        keys = set.union(*(set(index) for index in indexes[arm].values()))
        for key in sorted(keys):
            values = [index.get(key, {}).get("value") for index in indexes[arm].values()]
            stat_seed = int(digest([bootstrap_seed, arm, key]), 16)
            summaries.append(dict(arm=arm, time_hours=key[0], readout=key[1], lineage=key[2],
                **summarize_values(values, bootstrap_replicates=bootstrap_replicates, seed=stat_seed)))
    baseline = baseline if baseline is not None else next(iter(arm_runs), None)
    if baseline in indexes:
        reference = indexes[baseline]
        for arm, by_seed in sorted(indexes.items()):
            if arm == baseline:
                continue
            if set(by_seed) != set(reference):
                raise ValueError("Paired comparisons require identical seed sets.")
            common = set.intersection(*(set(index) for index in [*reference.values(), *by_seed.values()]))
            for key in sorted(common):
                diffs = []
                for seed in sorted(reference):
                    a, b = by_seed[seed][key]["value"], reference[seed][key]["value"]
                    diffs.append(None if a is None or b is None else a-b)
                stat_seed = int(digest([bootstrap_seed, baseline, arm, "paired", key]), 16)
                contrasts.append(dict(arm=arm, baseline=baseline, time_hours=key[0], readout=key[1],
                    lineage=key[2], contrast="arm_minus_baseline_paired_by_seed",
                    **summarize_values(diffs, bootstrap_replicates=bootstrap_replicates, seed=stat_seed),
                    **standardized_paired_effect(diffs, bootstrap_replicates=bootstrap_replicates, seed=stat_seed)))
    return dict(summary=summaries, paired_differences=contrasts,
                excluded_arms=sorted(set(arm_runs)-set(eligible)))


def _write_csv(path, rows):
    if not rows:
        Path(path).write_text("")
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _validate_arms(arms):
    from .measurement import comparison_signature
    if not isinstance(arms, dict) or not arms:
        raise ValueError("arms must be a nonempty mapping of unique names to prepared bundles or paths.")
    bundles = {}
    for name, source in arms.items():
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", name):
            raise ValueError("Arm names must be 1–80 letters/digits, underscores or hyphens, starting with a letter/digit.")
        bundle = deepcopy(source) if isinstance(source, dict) else read_json(source)
        bundles[name] = validate_prepared(bundle)
    reference = next(iter(bundles.values()))
    ref_hash = reference["guidebook"]["provenance"]["input_sha256"]
    ref_root = reference["guidebook"].get("control", {}).get("source_input_sha256", ref_hash)
    for name, bundle in bundles.items():
        if bundle["context"].get("kernel", "cell_agent") != reference["context"].get("kernel", "cell_agent"):
            raise ValueError("Cell-agent and lattice kernels have different units and cannot share an ensemble comparison.")
        if (bundle["context"].get("kernel") == "lattice" and
                bundle["context"]["lattice"]["capacity"] != reference["context"]["lattice"]["capacity"]):
            raise ValueError("Matched lattice arms require the same occupancy capacity per bin for comparable units.")
        if comparison_signature(bundle.get("measurement")) != comparison_signature(reference.get("measurement")):
            raise ValueError("Matched ensemble arms require the same measurement semantics (platform, units and preparation).")
        provenance = bundle["guidebook"]["provenance"]
        ref_provenance = reference["guidebook"]["provenance"]
        for key in ("ontology_version", "ontology_sha256", "kinetics_version", "kinetics_sha256"):
            if provenance.get(key) != ref_provenance.get(key):
                raise ValueError(f"Matched ensemble arms differ in {key}; declare parameter changes through the same versioned model.")
        own_hash = provenance["input_sha256"]
        root = bundle["guidebook"].get("control", {}).get("source_input_sha256", own_hash)
        if root != ref_root:
            raise ValueError(f"Arm {name!r} has unrelated measured input; use a separate ensemble for another ROI.")
        if bundle["context"]["organism"] != reference["context"]["organism"]:
            raise ValueError("Cannot pair arms from different organisms.")
        if bundle["guidebook"]["provenance"]["gene_sets_sha256"] != reference["guidebook"]["provenance"]["gene_sets_sha256"]:
            raise ValueError("Matched arms require the same gene signatures.")
    return bundles


def run_ensemble(arms, output_dir, *, seeds=(0, 1, 2, 3, 4), steps=None, horizon_hours=None,
                 mode="local", jev_weight=.5, max_population=2000, bootstrap_replicates=2000,
                 bootstrap_seed=0, readouts_path=None, provider_factory=None, arm_options=None):
    """Execute each arm using each seed and persist run and analysis provenance.

Exactly one of steps/horizon_hours is required. provider_factory is called once
per run with arm, bundle, output_dir, seed, steps and jev_weight keywords; it
must produce an independently journaled provider. No provider is created in
    local mode. arm_options may override mode and jev_weight by arm name so a
    local baseline can be compared to Jev arms. An existing output directory is
    never reused or overwritten.
"""
    bundles = _validate_arms(arms)
    seeds = list(seeds)
    if not seeds or any(type(seed) is not int for seed in seeds) or len(set(seeds)) != len(seeds):
        raise ValueError("Provide a nonempty set of unique integer seeds.")
    if type(bootstrap_seed) is not int:
        raise ValueError("bootstrap_seed must be an integer.")
    summarize_values([], bootstrap_replicates=bootstrap_replicates)
    if (steps is None) == (horizon_hours is None):
        raise ValueError("Specify exactly one of steps or horizon_hours.")
    if steps is not None and (type(steps) is not int or not 1 <= steps <= 10000):
        raise ValueError("steps must be an integer from 1 to 10000.")
    if horizon_hours is not None and (type(horizon_hours) not in (int, float) or not math.isfinite(horizon_hours) or horizon_hours <= 0):
        raise ValueError("horizon_hours must be finite and positive.")
    if mode not in ("local", "jev"):
        raise ValueError("mode must be local or jev.")
    if type(jev_weight) not in (int, float) or not math.isfinite(jev_weight) or not 0 <= jev_weight <= 1:
        raise ValueError("jev_weight must be finite in [0,1].")
    if type(max_population) is not int or max_population < max(len(b["cells"]) for b in bundles.values()):
        raise ValueError("max_population must cover every starting arm.")
    arm_options = {} if arm_options is None else arm_options
    if not isinstance(arm_options, dict) or set(arm_options) - set(bundles):
        raise ValueError("arm_options must map existing arm names to mode/jev_weight settings.")
    settings = {}
    for name in bundles:
        options = arm_options.get(name, {})
        if not isinstance(options, dict) or set(options) - {"mode", "jev_weight"}:
            raise ValueError("arm_options accepts only mode and jev_weight.")
        setting = dict(mode=mode, jev_weight=jev_weight)
        setting.update(options)
        if setting["mode"] not in ("local", "jev") or type(setting["jev_weight"]) not in (int, float) or not math.isfinite(setting["jev_weight"]) or not 0 <= setting["jev_weight"] <= 1:
            raise ValueError("Each arm needs local/jev mode and a finite jev_weight in [0,1].")
        settings[name] = setting
    needs_provider = any(setting["mode"] == "jev" for setting in settings.values())
    if needs_provider and provider_factory is None:
        raise ValueError("Jev arms require an explicit provider_factory.")
    if not needs_provider and provider_factory is not None:
        raise ValueError("A provider_factory is only permitted when an arm uses Jev.")
    arm_steps = {}
    for name, bundle in bundles.items():
        ratio = horizon_hours / bundle["context"]["time_step_hours"] if horizon_hours is not None else steps
        if not math.isclose(ratio, round(ratio), rel_tol=0, abs_tol=1e-9) or not 1 <= round(ratio) <= 10000:
            raise ValueError(f"Horizon must equal 1–10000 whole time steps for arm {name!r}.")
        arm_steps[name] = int(round(ratio))
    if steps is not None and len({b["context"]["time_step_hours"] for b in bundles.values()}) > 1:
        raise ValueError("Different time steps require horizon_hours to compare the same simulated duration.")
    spec, spec_hash = load_readout_spec(readouts_path)
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=False)
    configuration = dict(schema_version="3.0", ensemble_version="2", seeds=seeds,
        arms={name: dict(bundle_fingerprint=b["fingerprint"], steps=arm_steps[name], **settings[name],
                        time_step_hours=b["context"]["time_step_hours"], context=b["context"],
                        control=b["guidebook"].get("control")) for name, b in bundles.items()},
        mode=mode, jev_weight=jev_weight, max_population=max_population, bootstrap_replicates=bootstrap_replicates,
        bootstrap_seed=bootstrap_seed, readout_sha256=spec_hash,
        analysis_unit="simulation seed; not cell or biological replicate", baseline=next(iter(bundles)),
        time_alignment="exact common committed hours; no interpolation or extrapolation",
        uncertainty="percentile bootstrap confidence interval of median, separately from empirical 2.5–97.5% seed quantiles",
        standardized_paired_effect="Cohen dz = mean(arm-baseline paired seed differences)/sample SD of differences (ddof=1); percentile bootstrap of paired seeds; null for undefined denominator or any undefined bootstrap ratio; simulation-seed effect only",
        failure_policy="all runs recorded; any incomplete seed excludes its entire arm from aggregates")
    configuration["fingerprint"] = digest(configuration)
    write_json(directory / "config.json", configuration)
    write_json(directory / "readouts.json", spec)
    run_records, public_runs = [], []
    for arm, bundle in bundles.items():
        write_json(directory / "bundles" / f"{arm}.json", bundle)
        for seed in seeds:
            output = directory / "runs" / arm / f"seed_{seed}"
            record = dict(arm=arm, seed=seed, directory=str(output.relative_to(directory)), completed=False)
            try:
                setting = settings[arm]
                provider = provider_factory(arm=arm, bundle=bundle, output_dir=output, seed=seed,
                                            steps=arm_steps[arm], jev_weight=setting["jev_weight"]) if setting["mode"] == "jev" else None
                report = run(bundle, output, steps=arm_steps[arm], seed=seed, mode=setting["mode"], provider=provider,
                             jev_weight=setting["jev_weight"], max_population=max_population)
                record.update(completed=bool(report["completed"]), completed_steps=report["completed_steps"],
                              requested_steps=arm_steps[arm], stop_reason=report["stop_reason"],
                              error=report.get("error"), usage=report.get("usage", {}),
                              config_hash=report.get("config_hash"), bundle_fingerprint=bundle["fingerprint"])
                records = extract_readouts(output, context=bundle["context"], readouts_path=readouts_path)
                _write_csv(output / "readouts.csv", records)
                run_records.append(dict(record, readouts=records))
            except (ValueError, RuntimeError, OSError) as exc:
                record.update(completed=False, stop_reason="ensemble_run_error", error=str(exc))
                run_records.append(dict(record, readouts=[]))
            public_runs.append(record)
            write_json(directory / "runs.json", public_runs)
    aggregates = aggregate_readouts(run_records, bootstrap_replicates=bootstrap_replicates,
                                    bootstrap_seed=bootstrap_seed, baseline=configuration["baseline"])
    _write_csv(directory / "summary.csv", aggregates["summary"])
    _write_csv(directory / "paired_differences.csv", aggregates["paired_differences"])
    tidy = [dict(arm=r["arm"], seed=r["seed"], completed=r["completed"], **row)
            for r in run_records for row in r["readouts"]]
    _write_csv(directory / "readouts.csv", tidy)
    failures = [r for r in public_runs if not r["completed"]]
    warnings = (["Fewer than five seeds: engineering check only; empirical envelopes are unstable."] if len(seeds) < 5 else [])
    warnings += ["Intervals describe simulation-seed uncertainty, not patient variability or biological validity."]
    usage_keys = ("http_attempts", "cache_hits", "input_tokens", "output_tokens", "attempts_with_unknown_usage")
    by_arm_usage = {arm: {key: sum(r.get("usage", {}).get(key, 0) for r in public_runs if r["arm"] == arm)
                          for key in usage_keys} for arm in bundles}
    ledger = dict(by_arm=by_arm_usage, total={key: sum(row[key] for row in by_arm_usage.values()) for key in usage_keys},
                  currency_cost=None, note="Realized provider usage; no currency price is assumed. Unknown usage remains separately counted.")
    report = dict(schema_version="3.0", completed=not failures, fingerprint=configuration["fingerprint"],
        n_requested_runs=len(seeds)*len(bundles), n_completed_runs=len(public_runs)-len(failures),
        failures=failures, excluded_arms=aggregates["excluded_arms"], warnings=warnings,
        readout_sha256=spec_hash, runs=public_runs, cost_ledger=ledger,
        outputs=dict(summary="summary.csv", paired_differences="paired_differences.csv", readouts="readouts.csv"),
        biological_validation="not established")
    write_json(directory / "report.json", report)
    return report
