"""Explicit independent-input sensitivity designs and local simulation runner.

LHS is a screening design; it is never passed to the Sobol estimators. Sobol
indices use paired A/B/AB_i evaluations of an independent-input pick-freeze
design. The estimand is the response averaged over a fixed set of run seeds.
"""
from __future__ import annotations

from copy import deepcopy
import math
from pathlib import Path
import random
from statistics import mean, pvariance

from .ensemble import quantile, run_ensemble
from .inputs import digest, prepare, read_json, validate_prepared, write_json
from .readouts import extract_readouts, load_readout_spec


def _parameters(bounds):
    if not isinstance(bounds, dict) or not bounds:
        raise ValueError("bounds must map parameter names to [low, high] or a bound specification.")
    params = {}
    for name, value in bounds.items():
        if not isinstance(name, str) or not name:
            raise ValueError("Parameter names must be nonempty strings.")
        if isinstance(value, (list, tuple)) and len(value) == 2:
            spec = dict(low=value[0], high=value[1], distribution="uniform")
        elif isinstance(value, dict) and set(value) <= {"low", "high", "distribution"}:
            spec = dict(value)
            spec.setdefault("distribution", "uniform")
        else:
            raise ValueError("Each bound needs low, high and optional uniform/uniform_integer distribution.")
        low, high = spec.get("low"), spec.get("high")
        if any(type(v) not in (int, float) or not math.isfinite(v) for v in (low, high)) or low >= high:
            raise ValueError("Parameter bounds must be finite and strictly increasing.")
        if spec["distribution"] not in ("uniform", "uniform_integer"):
            raise ValueError("Only independent uniform and uniform_integer inputs are implemented.")
        if spec["distribution"] == "uniform_integer" and any(type(v) is not int for v in (low, high)):
            raise ValueError("uniform_integer requires integer inclusive bounds.")
        params[name] = spec
    return params


def make_design(bounds, *, n=200, method="sobol", seed=0):
    """Return named parameter points and hash; no simulation is run here.

sobol here names the target indices, not a quasi-random Sobol sequence: the
design uses independent pseudo-random A/B matrices with pick-freeze hybrids.
"""
    params = _parameters(bounds)
    if type(n) is not int or n < 2 or type(seed) is not int:
        raise ValueError("Use n >= 2 and an integer design seed.")
    if method not in ("sobol", "lhs"):
        raise ValueError("method must be sobol or lhs.")
    rng, names = random.Random(seed), list(params)
    def transform(name, u):
        spec = params[name]
        if spec["distribution"] == "uniform_integer":
            return spec["low"] + min(spec["high"]-spec["low"], int(u*(spec["high"]-spec["low"]+1)))
        return spec["low"] + u*(spec["high"]-spec["low"])
    points = []
    if method == "lhs":
        columns = {}
        for name in names:
            strata = list(range(n))
            rng.shuffle(strata)
            columns[name] = [transform(name, (i+rng.random())/n) for i in strata]
        points = [dict(id=f"LHS_{i:06d}", row=i, matrix="LHS",
                       parameters={name: columns[name][i] for name in names}) for i in range(n)]
    else:
        for i in range(n):
            a = {name: transform(name, rng.random()) for name in names}
            b = {name: transform(name, rng.random()) for name in names}
            points.extend([dict(id=f"A_{i:06d}", row=i, matrix="A", parameters=a),
                           dict(id=f"B_{i:06d}", row=i, matrix="B", parameters=b)])
            for j, name in enumerate(names):
                hybrid = dict(a)
                hybrid[name] = b[name]
                points.append(dict(id=f"AB{j}_{i:06d}", row=i, matrix=f"AB{j}", parameters=hybrid))
    result = dict(schema_version="3.0", design_version="1", method=method, n_base=n, seed=seed,
                  parameters=params, n_evaluations=len(points), points=points,
                  assumption="independent declared input distributions",
                  sampling="independent pseudo-random A/B with one-column swaps" if method == "sobol" else "Latin hypercube")
    result["fingerprint"] = digest(result)
    return result


def analyze_sobol(design, evaluations, *, bootstrap_replicates=1000, seed=0):
    """First/total order indices with paired-row bootstrap confidence intervals.

evaluations maps each point ID to one finite scalar. For stochastic simulations
this must be the mean over the SAME declared run seeds at every design point.
Indices outside [0,1] are retained as Monte-Carlo error, never clipped.
"""
    if design.get("method") != "sobol":
        raise ValueError("Sobol indices require the A/B/AB_i design, not LHS or arbitrary parameter draws.")
    if digest({k: v for k, v in design.items() if k != "fingerprint"}) != design.get("fingerprint"):
        raise ValueError("Sensitivity design changed after generation.")
    ids = [point["id"] for point in design["points"]]
    if len(set(ids)) != len(ids) or set(evaluations) != set(ids):
        raise ValueError("Exactly one complete evaluation is required for every design point.")
    if any(type(v) not in (int, float) or not math.isfinite(v) for v in evaluations.values()):
        raise ValueError("Missing, failed or nonfinite evaluations cannot enter Sobol analysis.")
    if type(bootstrap_replicates) is not int or bootstrap_replicates < 100 or type(seed) is not int:
        raise ValueError("Use at least 100 bootstrap replicates and an integer seed.")
    n, names = design["n_base"], list(design["parameters"])
    rows = {matrix: [None]*n for matrix in ("A", "B", *(f"AB{j}" for j in range(len(names))))}
    for point in design["points"]:
        rows[point["matrix"]][point["row"]] = float(evaluations[point["id"]])
    if any(value is None for values in rows.values() for value in values):
        raise ValueError("Incomplete pick-freeze design.")
    def estimate(indices, j):
        a, b, ab = ([rows[m][i] for i in indices] for m in ("A", "B", f"AB{j}"))
        center = mean(a+b)
        variance = pvariance(a+b)
        if variance == 0:
            return None, None
        first = mean((bv-center)*(abv-av) for av, bv, abv in zip(a, b, ab)) / variance
        total = .5*mean((av-abv)**2 for av, abv in zip(a, ab)) / variance
        return first, total
    rng = random.Random(seed)
    samples = [rng.choices(range(n), k=n) for _ in range(bootstrap_replicates)]
    results = {}
    for j, name in enumerate(names):
        first, total = estimate(range(n), j)
        reps = [estimate(indices, j) for indices in samples]
        f = [r[0] for r in reps if r[0] is not None]
        t = [r[1] for r in reps if r[1] is not None]
        results[name] = dict(S1=first, ST=total, S1_ci95=[quantile(f, .025), quantile(f, .975)],
            ST_ci95=[quantile(t, .025), quantile(t, .975)], bootstrap_valid=len(f),
            status="constant_response_undefined" if first is None else "ok")
    return dict(design_fingerprint=design["fingerprint"], n_base=n, indices=results,
        estimator="centered Saltelli first order and squared-difference total effect; paired base-row percentile bootstrap",
        references=["https://doi.org/10.1016/j.cpc.2009.09.018",
                    "https://salib.readthedocs.io/en/latest/_modules/SALib/analyze/sobol.html"],
        note="Indices depend on declared independent distributions; finite-sample estimates are not clipped to [0,1].")


def _override_context(context, values):
    context = deepcopy(context)
    for path, value in values.items():
        parts = path.split(".")
        allowed = ((len(parts) == 3 and parts[0] == "action_priors") or
                   (len(parts) in (2, 3) and parts[:2] == ["scoring", "beta"]) or
                   (len(parts) == 2 and parts[0] == "dynamics") or path == "neighbor_radius_um")
        if not allowed:
            raise ValueError(f"Unsupported local sensitivity parameter {path!r}; allowed: action_priors.LINEAGE.ACTION, scoring.beta[.ACTION], dynamics.KEY, neighbor_radius_um.")
        if value < 0 and parts[:2] != ["scoring", "beta"]:
            raise ValueError("Hazards, beta and spatial/dynamics magnitudes must be nonnegative in this sensitivity runner.")
        target = context
        for key in parts[:-1]:
            if key not in target:
                target[key] = {}
            if not isinstance(target[key], dict):
                raise ValueError(f"Cannot set {path}: {key} is not a mapping. Declare a per-action beta map first.")
            target = target[key]
        if parts[0] == "dynamics" and parts[-1] not in target:
            raise ValueError("Sensitivity can vary only an existing declared dynamics field.")
        target[parts[-1]] = value
    return context


def run_sensitivity(bundle, bounds, output_dir, *, horizon_hours, seeds=(0, 1, 2), n=200,
                    method="sobol", design_seed=0, readout="living_cells", lineage="malignant",
                    max_population=2000, bootstrap_replicates=1000):
    """Local F-10 runner. The response is the final readout MEAN across run seeds.

Hazard overrides are absolute per-hour values, not multipliers. Pooling alpha
is deliberately unavailable: it has no effect in local mode and needs a
separately costed Jev design. Every failed or undefined point blocks indices.
"""
    bundle = read_json(bundle) if isinstance(bundle, (str, Path)) else deepcopy(bundle)
    validate_prepared(bundle)
    if bundle["context"]["kinetics"]["mode"] != "hazard_per_hour":
        raise ValueError("Rate sensitivity requires hazard_per_hour kinetics.")
    if type(horizon_hours) not in (int, float) or not math.isfinite(horizon_hours) or horizon_hours <= 0:
        raise ValueError("horizon_hours must be finite and positive.")
    spec, _ = load_readout_spec()
    if readout not in {entry["name"] for entry in spec["readouts"]}:
        raise ValueError("Unknown readout; choose a name in analysis/readouts.yaml.")
    lattice_readouts = {"occupied_bins", "total_occupancy", "cumulative_occupancy_gained", "cumulative_occupancy_lost"}
    if bundle["context"].get("kernel") == "lattice" and readout not in lattice_readouts:
        raise ValueError("Local lattice sensitivity requires an occupancy/bin readout, not a cell or provider readout.")
    if bundle["context"].get("kernel", "cell_agent") != "lattice" and readout in lattice_readouts:
        raise ValueError("Occupancy/bin readouts apply only to the lattice kernel.")
    if readout == "jev_prior_js_bits":
        raise ValueError("Jev divergence is undefined in a local sensitivity sweep.")
    if lineage not in {"ALL", *(c["lineage"] for c in bundle["cells"])}:
        raise ValueError("Choose ALL or a lineage present in the source bundle.")
    seeds = list(seeds)
    if not seeds or len(set(seeds)) != len(seeds) or any(type(s) is not int for s in seeds):
        raise ValueError("Sensitivity needs unique integer simulation seeds.")
    design = make_design(bounds, n=n, method=method, seed=design_seed)
    for point in design["points"]:
        _override_context(bundle["context"], point["parameters"])
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/"design.json", design)
    write_json(output/"source_bundle.json", bundle)
    config = dict(source_fingerprint=bundle["fingerprint"], design_fingerprint=design["fingerprint"],
        seeds=seeds, horizon_hours=horizon_hours, readout=readout, lineage=lineage,
        response="mean final readout over the same simulation seeds at every parameter point",
        max_population=max_population, mode="local", hazard_override_semantics="absolute per-hour hazards before scoring/modifiers",
        pooling_alpha="not varied: no effect in local mode")
    write_json(output/"config.json", config)
    evaluations, records = {}, []
    for point in design["points"]:
        record = dict(id=point["id"], parameters=point["parameters"], completed=False)
        try:
            context = _override_context(bundle["context"], point["parameters"])
            extra = {}
            if bundle.get("schema_version") == "3.0":
                extra = dict(measurement=bundle.get("measurement"), annotation_provenance=bundle.get("annotation_provenance"),
                             measured_evidence=bundle.get("measured_evidence"))
            prepared = prepare(bundle["cells"], context, assayed_genes=context.get("assayed_genes"),
                gene_sets=bundle["guidebook"]["provenance"]["gene_sets"], **extra)
            folder = output/"points"/point["id"]
            result = run_ensemble({"response": prepared}, folder, seeds=seeds, horizon_hours=horizon_hours,
                                  bootstrap_replicates=100, max_population=max_population)
            if not result["completed"]:
                raise ValueError("One or more simulation seeds did not reach the requested horizon.")
            values = []
            for seed in seeds:
                rows = extract_readouts(folder/"runs/response"/f"seed_{seed}")
                found = [r for r in rows if r["readout"] == readout and r["lineage"] == lineage
                         and math.isclose(r["time_hours"], horizon_hours, rel_tol=0, abs_tol=1e-9)]
                if len(found) != 1 or found[0]["value"] is None:
                    raise ValueError("Requested final readout is missing or undefined.")
                values.append(found[0]["value"])
            evaluations[point["id"]] = mean(values)
            record.update(completed=True, seed_values=values, response=mean(values))
        except (ValueError, RuntimeError, OSError) as exc:
            record["error"] = str(exc)
        records.append(record)
        write_json(output/"evaluations.json", records)
    complete = all(record["completed"] for record in records)
    indices = analyze_sobol(design, evaluations, bootstrap_replicates=bootstrap_replicates,
                           seed=design_seed) if complete and method == "sobol" else None
    if indices is not None:
        write_json(output/"sobol.json", indices)
    report = dict(schema_version="3.0", completed=complete, design_fingerprint=design["fingerprint"],
        n_parameter_points=len(records), n_complete_points=len(evaluations), n_requested_runs=len(records)*len(seeds),
        method=method, sobol_indices_written=indices is not None,
        failures=[record for record in records if not record["completed"]],
        interpretation="Sensitivity of a local model under declared independent parameter distributions, not calibration or biological validation.",
        limitation="A finite seed mean approximates the expected stochastic response; convergence in design size and simulation replicates remains required.")
    write_json(output/"report.json", report)
    return report
