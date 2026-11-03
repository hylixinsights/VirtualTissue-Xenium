#!/usr/bin/env python3
"""Small, local-only v3 clock gate with an analytical birth/death benchmark.

Run: python tools/validate_clock_v3.py --output outputs/clock-validation

This checks a homogeneous malignant population with no expression modulation,
interaction effects, motion or crowding. It does NOT validate all lineage/state
trajectories, tissue structure, biology, or Jev clock semantics. A hidden-prior
Jev arm uses raw alpha=1 action preferences and is not constrained to preserve
these local guidebook hazards. Seeds are simulation replicates, not biological
replicates; the same seed labels are not a coupled cell-hour random stream.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from virtualtissue_xenium.v2.frame_store import iter_frames
from virtualtissue_xenium.v2.inputs import prepare, write_json
from virtualtissue_xenium.v2.simulation import decision_rows, initialize, run

BIRTH_PER_HOUR = .006
DEATH_PER_HOUR = .001


def branching_moments(n0, birth, death, dt, steps):
    """Exact moments for the kernel's independent one-event branching process.

    Each alive parent contributes zero, one or two alive descendants with
    probabilities p_death, p_rest, p_birth. Conditional offspring are independent
    under the synthetic local-only fixture; no living-population cap may bind.
    """
    total = birth + death
    any_event = -math.expm1(-total * dt)
    p_birth, p_death = any_event * birth / total, any_event * death / total
    offspring_mean = 1 + p_birth - p_death
    offspring_variance = p_birth + p_death - (p_birth - p_death) ** 2
    mean = n0 * offspring_mean ** steps
    if steps == 0:
        return float(n0), 0.0
    if math.isclose(offspring_mean, 1.0, abs_tol=1e-14):
        variance = n0 * offspring_variance * steps
    else:
        variance = n0 * offspring_variance * offspring_mean ** (steps - 1) * (offspring_mean ** steps - 1) / (offspring_mean - 1)
    return mean, variance


def make_bundle(n, dt):
    context = dict(schema_version="3.0", organism="human", tissue="synthetic malignant population",
        disease="engineering clock benchmark", assay="scrna", kernel="cell_agent", interaction_mode="well_mixed",
        sampling=dict(design="subsampled", description="Identical synthetic malignant agents; not biological data", missing_populations=[]),
        perturbations=[], recognition=False, time_step_hours=dt,
        kinetics=dict(mode="hazard_per_hour", max_step_probability=.25),
        dynamics=dict(movement_enabled=False, crowding_enabled=False),
        scoring=dict(method="depth_conditioned_auc", beta=0, min_transcripts=0),
        action_priors={"malignant": {"PROLIFERATE": BIRTH_PER_HOUR, "APOPTOSIS": DEATH_PER_HOUR, "SECRETE": 0}})
    cells = [dict(cell_id=f"mal_{i:04d}", cell_type="synthetic malignant", lineage="malignant", sample_id="clock",
                  expression={"ACTB": 100}) for i in range(n)]
    bundle = prepare(cells, context, assayed_genes=["ACTB"])
    rows, _ = decision_rows(initialize(bundle), bundle["guidebook"], bundle["context"], seed=0, step=0)
    expected = {"PROLIFERATE": BIRTH_PER_HOUR, "APOPTOSIS": DEATH_PER_HOUR}
    if any(r["effective_hazards_per_hour"] != expected for r in rows):
        raise RuntimeError("Clock fixture acquired extra hazards or modifiers; update its explicit analytical contract.")
    return bundle


def draw_figure(result, directory):
    # Set a writable cache before importing matplotlib, useful in restricted CI.
    os.environ.setdefault("MPLCONFIGDIR", str(directory / ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.7), gridspec_kw={"width_ratios": [1.35, 1]})
    colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7", "#E69F00", "#56B4E9"]
    p = result["parameters"]
    horizon, n0 = p["horizon_hours"], p["initial_cells"]
    times = [horizon * i / 100 for i in range(101)]
    left.plot(times, [n0 * math.exp((BIRTH_PER_HOUR - DEATH_PER_HOUR) * t) for t in times],
              color="#222222", linewidth=2, linestyle="--", label="Continuous birth–death mean")
    for index, arm in enumerate(result["arms"]):
        color = colors[index % len(colors)]
        curve = arm["curve"]
        x = [r["time_hours"] for r in curve]
        y = [r["observed_mean"] for r in curve]
        left.plot(x, y, color=color, linewidth=1.5, label=f"Δt = {arm['dt_hours']:g} h")
        left.fill_between(x, [r["discrete_expected_mean"] - r["mc_tolerance"] for r in curve],
                          [r["discrete_expected_mean"] + r["mc_tolerance"] for r in curve], color=color, alpha=.055)
        terminal = arm["terminal"]
        jitter = [(j - (len(arm["runs"]) - 1) / 2) * .035 for j in range(len(arm["runs"]))]
        right.scatter([index + delta for delta in jitter], [r["terminal_living"] for r in arm["runs"]],
                      color=color, alpha=.6, s=22, zorder=3)
        right.errorbar(index, terminal["discrete_expected_mean"], yerr=terminal["mc_tolerance"],
                       color=color, marker="_", markersize=12, capsize=4, linewidth=1.3, zorder=2)
        right.scatter(index, terminal["observed_mean"], color=color, marker="D", s=38, edgecolors="white", linewidths=.6, zorder=4)
    right.axhline(result["continuous_expected_terminal"], color="#222222", linestyle="--", linewidth=1.4)
    right.set_xticks(range(len(result["arms"])), [f"{a['dt_hours']:g}" for a in result["arms"]])
    right.set_xlabel("Time step Δt (hours)")
    right.set_ylabel(f"Living agents at {horizon:g} hours")
    right.set_title("Common horizon: per-seed outcomes")
    left.set_xlabel("Simulated time (hours)")
    left.set_ylabel("Living malignant agents")
    left.set_title("Local guidebook-only trajectories")
    left.legend(frameon=False, fontsize=8, loc="upper left")
    for axis in (left, right):
        axis.grid(axis="y", alpha=.15)
    fig.suptitle(f"v3 clock engineering gate · {n0} initial agents · {len(p['seeds'])} seeds", fontsize=12, x=.07, ha="left")
    fig.text(.07, .025, f"Bars/ribbons: exact discrete expectation ± {p['mc_sigma']:g} Monte Carlo standard errors of the seed mean.\n"
             "Dots: seeds; diamonds: means. Synthetic local-only check, not biological or Jev validation.",
             fontsize=8, color="#444444")
    fig.subplots_adjust(left=.07, right=.98, top=.83, bottom=.22, wspace=.27)
    fig.savefig(directory / "clock_validation.png", dpi=180, facecolor="white")
    fig.savefig(directory / "clock_validation.svg", facecolor="white")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cells", type=int, default=80)
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    parser.add_argument("--time-steps", type=float, nargs="+", default=[.5, 1, 2, 4])
    parser.add_argument("--horizon-hours", type=float, default=24)
    parser.add_argument("--mc-sigma", type=float, default=4)
    parser.add_argument("--max-relative-discretization-bias", type=float, default=.01)
    args = parser.parse_args(argv)
    if args.cells < 2 or len(args.seeds) < 2 or len(set(args.seeds)) != len(args.seeds):
        parser.error("Use at least two initial agents and two unique simulation seeds.")
    if (not math.isfinite(args.horizon_hours) or args.horizon_hours <= 0 or not math.isfinite(args.mc_sigma)
            or args.mc_sigma <= 0 or not math.isfinite(args.max_relative_discretization_bias)
            or not 0 < args.max_relative_discretization_bias < 1):
        parser.error("Horizon and sigma must be finite and positive; relative bias bound must lie in (0,1).")
    for dt in args.time_steps:
        if not math.isfinite(dt) or dt <= 0 or not math.isclose(args.horizon_hours / dt, round(args.horizon_hours / dt), abs_tol=1e-9):
            parser.error("Every finite positive time step must divide the same horizon exactly.")
    if len(set(args.time_steps)) != len(args.time_steps):
        parser.error("Time steps must be unique.")
    directory = args.output.resolve()
    if directory.exists() and any(directory.iterdir()):
        parser.error("Output directory must be empty; prior validation evidence is never overwritten.")
    directory.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    expected = args.cells * math.exp((BIRTH_PER_HOUR - DEATH_PER_HOUR) * args.horizon_hours)
    result = dict(schema_version="1.0", gate="v3_simple_birth_death_clock", mode="local", passed=False,
        parameters=dict(initial_cells=args.cells, seeds=args.seeds, time_steps_hours=args.time_steps,
                        horizon_hours=args.horizon_hours, mc_sigma=args.mc_sigma,
                        max_relative_discretization_bias=args.max_relative_discretization_bias,
                        birth_per_hour=BIRTH_PER_HOUR, death_per_hour=DEATH_PER_HOUR),
        continuous_expected_terminal=expected,
        gate_rule="At the common terminal horizon, every arm: |observed seed mean-exact discrete mean| <= sigma*sqrt(exact discrete variance/n_seeds), "
                  "and relative deterministic discrete/continuous bias <= declared bound; all runs complete without capacity events.",
        uncertainty="A conservative engineering Monte Carlo tolerance from the known branching process, not a biological confidence interval or an estimated seed envelope.",
        limitations=["Only homogeneous malignant birth/death under local guidebook rates is checked.",
                     "No expression modulation, crowding, motility, cell-cell effects or perturbation.",
                     "Does not validate all tissue/state curves or biological rates.",
                     "Jev alpha=1 raw preferences are not constrained to preserve local hazard clock semantics.",
                     "Seed labels repeat across time steps but RNG streams are cell/step, not coupled cell/hour.",
                     "Failure or success is an engineering gate; no tests on simulated cells as independent biological observations."],
        environment=dict(python=platform.python_version(), script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), arms=[])
    source_hashes = {}
    for code_dir in (ROOT / "src/virtualtissue_xenium/v2", ROOT / "src/virtualtissue_xenium/v3"):
        for source in sorted(code_dir.rglob("*")):
            if source.is_file() and source.suffix in {".py", ".json", ".yaml"}:
                source_hashes[str(source.relative_to(ROOT))] = hashlib.sha256(source.read_bytes()).hexdigest()
    result["environment"]["source_file_sha256"] = source_hashes
    try:
        result["environment"]["git_revision"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
        result["environment"]["working_tree_modified"] = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip())
    except (OSError, subprocess.CalledProcessError):
        result["environment"]["git_revision"] = None
    for dt in args.time_steps:
        steps = int(round(args.horizon_hours / dt))
        name = f"dt_{dt:g}".replace(".", "p")
        bundle = make_bundle(args.cells, dt)
        write_json(directory / "prepared" / f"{name}.json", bundle)
        arm = dict(dt_hours=dt, steps=steps, bundle_fingerprint=bundle["fingerprint"], runs=[], curve=[])
        observed = []
        for seed in args.seeds:
            destination = directory / "runs" / name / f"seed_{seed}"
            report = run(bundle, destination, steps=steps, seed=seed, max_population=max(2000, args.cells * 10))
            counts = [sum(c["status"] == "alive" for c in frame["cells"].values()) for frame in iter_frames(destination / "frames.jsonl")]
            observed.append(counts)
            arm["runs"].append(dict(seed=seed, directory=str(destination.relative_to(directory)), completed=report["completed"],
                                    terminal_living=report["living"], capacity_events=report["events"].get("population_cap_blocked", 0)))
        complete = all(r["completed"] and r["capacity_events"] == 0 for r in arm["runs"])
        for step in range(min(len(x) for x in observed)):
            mean, variance = branching_moments(args.cells, BIRTH_PER_HOUR, DEATH_PER_HOUR, dt, step)
            arm["curve"].append(dict(time_hours=step * dt, observed_mean=statistics.mean(x[step] for x in observed),
                                     discrete_expected_mean=mean, discrete_variance_per_seed=variance,
                                     mc_tolerance=args.mc_sigma * math.sqrt(variance / len(args.seeds))))
        terminal = arm["curve"][-1]
        terminal["continuous_expected_mean"] = expected
        terminal["absolute_discretization_bias"] = abs(terminal["discrete_expected_mean"] - expected)
        terminal["relative_discretization_bias"] = terminal["absolute_discretization_bias"] / expected
        terminal["total_tolerance_around_continuous_mean"] = terminal["absolute_discretization_bias"] + terminal["mc_tolerance"]
        arm["terminal"] = terminal
        arm["passed"] = bool(complete and math.isclose(terminal["time_hours"], args.horizon_hours, rel_tol=1e-12, abs_tol=1e-12)
            and abs(terminal["observed_mean"] - terminal["discrete_expected_mean"]) <= terminal["mc_tolerance"]
            and terminal["relative_discretization_bias"] <= args.max_relative_discretization_bias)
        result["arms"].append(arm)
    result["passed"] = all(a["passed"] for a in result["arms"])
    result["elapsed_seconds"] = time.monotonic() - started
    write_json(directory / "validation.json", result)
    draw_figure(result, directory)
    print(json.dumps(dict(passed=result["passed"], arms=len(result["arms"]), runs=sum(len(a["runs"]) for a in result["arms"]),
                          elapsed_seconds=round(result["elapsed_seconds"], 2), report=str(directory / "validation.json"),
                          figure=str(directory / "clock_validation.png")), indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
