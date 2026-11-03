"""Regenerate seed-summary figures from a completed ensemble's CSV and manifest.

Requires matplotlib (optional analysis dependency); does not rerun simulations.
Example: python analysis/figures.py --ensemble outputs/ensemble --output figures
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import textwrap


def create_figures(ensemble_dir, output_dir, *, readouts=None, lineages=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    source, output = Path(ensemble_dir), Path(output_dir)
    report = json.loads((source/"report.json").read_text())
    config = json.loads((source/"config.json").read_text())
    if not report.get("completed") or report.get("excluded_arms"):
        raise ValueError("Figures require a completed ensemble without excluded arms.")
    if report.get("fingerprint") != config.get("fingerprint"):
        raise ValueError("Ensemble report and configuration fingerprints differ.")
    with (source/"summary.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    available = sorted({row["readout"] for row in rows})
    chosen = list(readouts) if readouts else [name for name in ("living_cells", "occupied_bins", "total_occupancy",
                         "cumulative_kills", "malignant_activated_cd8_contact_fraction") if name in available]
    if not chosen or set(chosen)-set(available):
        raise ValueError("Choose readouts present in this ensemble summary.")
    if any(not re.fullmatch(r"[a-zA-Z0-9_]+", name) for name in chosen):
        raise ValueError("Unexpected readout name.")
    output.mkdir(parents=True, exist_ok=False)
    generated = []
    for readout in chosen:
        subset = [row for row in rows if row["readout"] == readout]
        groups = sorted({row["lineage"] for row in subset})
        if lineages:
            groups = [name for name in lineages if name in groups]
        if not groups:
            continue
        columns = min(2, len(groups))
        fig, axes = plt.subplots(math.ceil(len(groups)/columns), columns,
                                figsize=(5.2*columns, 3.5*math.ceil(len(groups)/columns)),
                                squeeze=False, layout="constrained")
        for axis, lineage in zip(axes.flat, groups):
            for arm in sorted({row["arm"] for row in subset}):
                points = sorted((row for row in subset if row["arm"] == arm and row["lineage"] == lineage),
                                key=lambda row: float(row["time_hours"]))
                if not points:
                    continue
                def values(key):
                    return [float(row[key]) if row.get(key) else float("nan") for row in points]
                times, center = values("time_hours"), values("median")
                line, = axis.plot(times, center, label=arm, linewidth=1.7)
                color = line.get_color()
                axis.fill_between(times, values("seed_q025"), values("seed_q975"), color=color, alpha=.10)
                axis.fill_between(times, values("median_ci95_low"), values("median_ci95_high"), color=color, alpha=.25)
            axis.set_title(lineage)
            axis.set_xlabel("Simulated time (hours)")
            axis.set_ylabel(textwrap.fill(readout.replace("_", " "), width=30))
            axis.grid(alpha=.2)
            axis.legend(frameon=False)
        for axis in list(axes.flat)[len(groups):]:
            axis.set_visible(False)
        fig.suptitle(readout.replace("_", " ").capitalize(), fontsize=14)
        caption = "Line: seed median · dark band: median 95% bootstrap CI"
        caption += ("\n" if columns == 1 else " · ") + "Light band: empirical 95% seed envelope"
        fig.supxlabel(caption, fontsize=8)
        for extension in ("svg", "png"):
            name = f"{readout}.{extension}"
            fig.savefig(output/name, dpi=180)
            generated.append(name)
        plt.close(fig)
    manifest = dict(ensemble=str(source.resolve()), ensemble_fingerprint=report["fingerprint"],
        readout_sha256=report["readout_sha256"], generated=generated,
        source_sha256={name: hashlib.sha256((source/name).read_bytes()).hexdigest()
                       for name in ("config.json", "report.json", "summary.csv")},
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        matplotlib_version=matplotlib.__version__, warnings=report.get("warnings", []),
        interpretation="Simulation seed summaries only; no biological replication or validation claim.")
    (output/"FIGURE_MANIFEST.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ensemble", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--readouts", nargs="+")
    parser.add_argument("--lineages", nargs="+")
    args = parser.parse_args(argv)
    result = create_figures(args.ensemble, args.output, readouts=args.readouts, lineages=args.lineages)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
