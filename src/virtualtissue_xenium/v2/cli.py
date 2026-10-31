"""Ordered collaborator workflow: vtx2 prepare -> inspect -> run -> compare."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

from .inputs import digest, load_cells, prepare, read_json, write_json
from .simulation import preflight, run


def cmd_prepare(args):
    context = read_json(args.context)
    if args.command == "prepare-xenium":
        from .datasets import load_xenium_cells
        cells, genes = load_xenium_cells(args.input, args.annotations, max_cells=args.max_cells,
                                       center=args.center, radius_um=args.radius_um)
    else:
        cells, genes = load_cells(args.input, cell_type_key=args.cell_type_key, lineage_key=args.lineage_key,
            sample_key=args.sample_key, layer=args.layer, coordinate_scale=args.coordinate_scale, max_cells=args.max_cells)
        if len(cells) > args.max_cells:
            raise ValueError("Subset the input first or explicitly increase --max-cells.")
    gene_sets = read_json(args.gene_sets) if args.gene_sets else None
    bundle = prepare(cells, context, assayed_genes=genes, gene_sets=gene_sets)
    report = preflight(bundle)
    sidecar = Path(args.output).with_suffix(".guidebook.json")
    if Path(args.output).exists() or sidecar.exists():
        raise ValueError("Prepared output or guidebook sidecar exists; choose a new bundle path.")
    write_json(args.output, bundle)
    write_json(sidecar, bundle["guidebook"])
    return dict(prepared=args.output, fingerprint=bundle["fingerprint"], preflight=report)


def cmd_run(args):
    from .provider import JevProvider
    bundle = read_json(args.prepared)
    facts = preflight(bundle, model=args.model, max_questions=args.max_questions, seed=args.seed)
    if args.mode == "jev":
        key = os.environ.get("TYPESAFE_API_KEY") or os.environ.get("JEV_API_KEY")
        if not key:
            raise ValueError("Set TYPESAFE_API_KEY or JEV_API_KEY before Jev mode.")
        identity = dict(bundle_fingerprint=bundle["fingerprint"], steps=args.steps, seed=args.seed,
                        mode=args.mode, jev_weight=args.jev_weight, max_population=args.max_population,
                        model=args.model, max_questions=args.max_questions)
        provider = JevProvider(key, Path(args.output)/"jev_journal", identity, model=args.model,
            max_attempts=args.max_http_attempts, workers=args.workers, max_questions=args.max_questions)
    else:
        provider = None
    report = run(bundle, args.output, steps=args.steps, seed=args.seed, mode=args.mode,
                 provider=provider, jev_weight=args.jev_weight, max_population=args.max_population, resume=args.resume)
    return dict(**report, replay=str(Path(args.output).resolve()/"replay.html"))


def compare(paths):
    """Require the same measurements; explicitly describe changed context/guidebook."""
    if len(paths) < 2:
        raise ValueError("Compare requires at least two run directories.")
    results = []
    for path in paths:
        folder = Path(path)
        config, report = read_json(folder/"config.json"), read_json(folder/"report.json")
        guidebook, frames = read_json(folder/"guidebook.json"), read_json(folder/"frames.json")
        results.append(dict(path=str(folder), config=config, report=report, guidebook=guidebook, frames=frames))
    ref = results[0]
    for result in results[1:]:
        for key in ("seed", "mode", "jev_weight", "model", "kernel", "max_population", "assumptions",
                    "provider_name", "provider_prompt_version", "max_questions"):
            if ref["config"][key] != result["config"][key]:
                raise ValueError(f"Unmatched run setting: {key}")
        if ref["guidebook"]["provenance"]["input_sha256"] != result["guidebook"]["provenance"]["input_sha256"]:
            raise ValueError("Runs have different measured input cells/expression/geometry.")
        if ref["guidebook"]["provenance"]["gene_sets_sha256"] != result["guidebook"]["provenance"]["gene_sets_sha256"]:
            raise ValueError("Runs have different gene sets.")
        for key in ("guidebook_version", "assayed_genes", "panel_source", "score_method", "symbol_normalization",
                    "default_action_weights", "recognition_licensed_kill_weight", "time_semantics"):
            if ref["guidebook"]["provenance"][key] != result["guidebook"]["provenance"][key]:
                raise ValueError(f"Unmatched guidebook scoring assumption: {key}")
        for key in ("time_step_hours", "interaction_mode", "neighbor_radius_um", "max_neighbors", "sampling", "organism", "assay", "tissue", "disease"):
            if ref["guidebook"]["context"][key] != result["guidebook"]["context"][key]:
                raise ValueError(f"Unmatched biological/context setting: {key}")
    shared = min(r["report"]["completed_steps"] for r in results)
    output = []
    for r in results:
        frame = next(f for f in r["frames"] if f["step"] == shared)
        cells = list(frame["cells"].values())
        output.append(dict(path=r["path"], shared_step=shared, time_hours=frame["time_hours"],
            completed=r["report"]["completed"], stop_reason=r["report"]["stop_reason"],
            living=sum(c["status"] == "alive" for c in cells), dead=sum(c["status"] == "dead" for c in cells),
            activated=sum(c["status"] == "alive" and c["activated"] for c in cells),
            exhausted=sum(c["status"] == "alive" and c["exhausted"] for c in cells),
            context=r["guidebook"]["context"],
            guidebook_fingerprint=digest(r["guidebook"])))
    return dict(shared_step=shared, runs=output,
                interpretation="Matched stochastic assumptions experiment; no causal treatment efficacy claim. Repeat across seeds.")


def cmd_dataset(args):
    from . import datasets
    if args.action == "list":
        return datasets.list_datasets()
    if not args.dataset or not args.output:
        raise ValueError("dataset download/verify requires DATASET_ID and OUTPUT_DIR.")
    if args.action == "download":
        return datasets.download_dataset(args.dataset, args.output, max_bytes=args.max_mib*2**20,
                                         include_optional=args.include_optional)
    return datasets.verify_dataset(args.dataset, args.output, include_optional=args.include_optional)


def cmd_annotation_template(args):
    from .datasets import write_annotation_template
    return dict(template=str(write_annotation_template(args.input, args.output, max_cells=args.max_cells,
                         center=args.center, radius_um=args.radius_um)))


def parser():
    p = argparse.ArgumentParser(prog="vtx2", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("prepare", "prepare-xenium"):
        s = sub.add_parser(name, help="Build a reviewed context-conditioned guidebook and prepared JSON.")
        s.add_argument("input"); s.add_argument("context"); s.add_argument("output")
        s.add_argument("--gene-sets"); s.add_argument("--max-cells", type=int, default=1000)
        if name == "prepare-xenium":
            s.add_argument("--annotations", required=True)
            s.add_argument("--center", nargs=2, type=float); s.add_argument("--radius-um", type=float)
        else:
            s.add_argument("--cell-type-key", default="cell_type"); s.add_argument("--lineage-key", default="lineage")
            s.add_argument("--sample-key"); s.add_argument("--layer"); s.add_argument("--coordinate-scale", type=float, default=1)
        s.set_defaults(func=cmd_prepare)
    s = sub.add_parser("inspect", help="Validate a prepared bundle; estimate first-round calls, no API use.")
    s.add_argument("prepared"); s.set_defaults(func=lambda a: preflight(read_json(a.prepared)))
    s = sub.add_parser("run", help="Synchronous local or Jev simulation; local is the default.")
    s.add_argument("prepared"); s.add_argument("output"); s.add_argument("--mode", choices=("local", "jev"), default="local")
    s.add_argument("--steps", type=int, default=12); s.add_argument("--seed", type=int, default=21)
    s.add_argument("--max-population", type=int, default=2000); s.add_argument("--jev-weight", type=float, default=.5)
    s.add_argument("--model", default="jev-1.14.0"); s.add_argument("--max-http-attempts", type=int, default=100)
    s.add_argument("--workers", type=int, default=2); s.add_argument("--max-questions", type=int, default=20)
    s.add_argument("--resume", action="store_true"); s.set_defaults(func=cmd_run)
    s = sub.add_parser("compare", help="Compare matched arms at their latest common completed step.")
    s.add_argument("runs", nargs="+"); s.set_defaults(func=lambda a: compare(a.runs))
    s = sub.add_parser("dataset", help="Verified public-data manifests and bounded downloads.")
    s.add_argument("action", choices=("list", "download", "verify")); s.add_argument("dataset", nargs="?")
    s.add_argument("output", nargs="?"); s.add_argument("--max-mib", type=int, default=100)
    s.add_argument("--include-optional", action="store_true"); s.set_defaults(func=cmd_dataset)
    s = sub.add_parser("annotation-template", help="Blank reviewed labels with exact Xenium IDs/vendor clusters.")
    s.add_argument("input"); s.add_argument("output"); s.add_argument("--max-cells", type=int, default=1000)
    s.add_argument("--center", nargs=2, type=float); s.add_argument("--radius-um", type=float)
    s.set_defaults(func=cmd_annotation_template)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        result = args.func(args)
        print(json.dumps(result, indent=2, default=str, allow_nan=False))
        if isinstance(result, dict) and result.get("completed") is False:
            return 2
        return 0
    except (ValueError, KeyError, OSError, RuntimeError) as e:
        print(f"vtx2: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
