"""Command-line interface: ``vtx <command>`` (or ``python -m virtualtissue_xenium <command>``).

The commands run the same engine functions as the notebooks, without widgets:

  vtx download-example DATA_DIR           Part 1 input: the public 10x sample (~503 MiB)
  vtx own-data XENIUM_OUTS DATA_DIR ...   Part 1 input: register your own Xenium run
  vtx annotation-template DATA_DIR OUT    Draft a cluster annotation map for your sample
  vtx check-annotations OUT [--data-dir]  Validate that map
  vtx prepare DATA_DIR OUTPUT_BASE ...    Part 1: AUCell, regions, export bundle ZIP(s)
  vtx check BUNDLE_ZIP                    Part 2: verify a bundle and preview the first round (no API calls)
  vtx run BUNDLE_ZIP RESULTS_DIR ...      Part 3: run Jev (or the free local preview)
  vtx explore RESULTS                     Part 4: summarise a results ZIP/folder, plot, open the replay

Run ``vtx <command> -h`` for every option.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import virtualtissue_xenium as vtx

SCENARIOS = ("baseline_unknown", "declared_recognition", "recognition_off")


def _print_json(obj) -> None:
    print(json.dumps(obj, indent=2, default=str))


# --------------------------------------------------------------------------- Part 1 inputs
def cmd_download_example(a):
    from .download import download_public_lung5k
    download_public_lung5k(a.data_dir)


def cmd_register_download(a):
    from .download import register_manual_download
    report = register_manual_download(a.data_dir)
    sys.exit(0 if report["ok"] else 1)


def cmd_verify(a):
    from .download import verify_data_folder
    report = verify_data_folder(a.data_dir, check_hashes=not a.quick)
    _print_json(report)
    sys.exit(0 if report["ok"] else 1)


def cmd_own_data(a):
    from .own_data import create_data_folder
    create_data_folder(a.xenium_outs, a.data_dir, a.sample_id, include_optional=a.include_optional,
                       source_page=a.source_page, license=a.license, attribution=a.attribution)


def cmd_annotation_template(a):
    from .own_data import write_annotation_template
    path = write_annotation_template(a.data_dir, a.out, top_n=a.top_n)
    print(f"Template written: {path}\nEdit label/compartment/runtime_label/markers for every cluster, "
          f"then: vtx check-annotations {path} --data-dir {a.data_dir}")


def cmd_check_annotations(a):
    from .own_data import validate_annotations
    problems = validate_annotations(a.annotations, a.data_dir)
    print("Valid." if not problems else "Problems:\n  " + "\n  ".join(problems))
    sys.exit(0 if not problems else 1)


# --------------------------------------------------------------------------- Part 1
def cmd_prepare(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from . import own_data

    prep = vtx.prepare
    from prepared_runtime import validate_bundle

    own_data.require_annotation_map(a.data_dir, a.annotations)
    run_dir = Path(a.output_base) / datetime.now(timezone.utc).strftime("run_%Y%m%d_%H%M%S_%f")
    if a.annotations:
        run_dir.mkdir(parents=True, exist_ok=True)
        own_data.use_cluster_annotations(prep, a.annotations, run_dir)
    print("Results:", run_dir)
    sample = prep.load_sample(a.data_dir, run_dir, verify_hashes=not a.no_verify_hashes, min_counts=20, min_genes=10)
    prep.save_json(run_dir / "environment.json", prep.environment())
    prep.annotate(sample, {})
    own_data.require_tumour_interface(sample)
    plt.close(prep.plot_marker_review(sample))
    prep.score_aucell(sample, n_up=a.aucell_top_ranks, min_overlap=5, batch_size=4096)
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.hist(sample["obs"].n_genes_by_counts, bins=100, color="#397f9c")
    ax.axvline(a.aucell_top_ranks, color="#bc582d", label="Minimum detected genes for AUCell")
    ax.set(xlabel="Detected genes per cell", ylabel="Number of cells"); ax.legend()
    fig.tight_layout(); fig.savefig(run_dir / "rna_depth_qc.png", dpi=150); plt.close(fig)
    prep.read_geometry(sample)
    table = prep.propose_regions(sample, n_options=a.n_regions, core_target=a.core_target,
                                 max_total=a.max_total, halo_um=a.halo)
    plt.close(prep.plot_overview(sample))
    plt.close(prep.plot_regions(sample))
    import pandas as pd
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(table[["region_id", "n_core", "n_halo", "n_total", "tumor_candidate", "fibroblast", "immune",
                     "stroma", "cd8_provisional", "macrophage_provisional", "x_um", "y_um"]].to_string(index=False))
    print("Maps: candidate_regions.png and whole_section_regions.png in", run_dir)
    regions = list(a.region or [])
    if a.center:
        prep.manual_region(sample, a.center[0], a.center[1], a.core_target, a.max_total, a.halo)
        regions.append("MANUAL")
    if a.list_only or not regions:
        print("No region exported (use --region R03, or --center X Y).")
        return
    known = {r["region_id"] for r in sample["regions"]}
    for region_id in regions:
        if region_id not in known:
            sys.exit(f"Unknown region {region_id}; available: {sorted(known)}")
        stamp = datetime.now(timezone.utc).strftime("%H%M%S_%f")
        bundle = prep.export_region(sample, region_id, tumor_mode=a.tumor_mode,
                                    destination=run_dir / f"simulation_input_{region_id}_{stamp}")
        report = validate_bundle(bundle, smoke_steps=1)
        prep.write_bundle_manifest(bundle)
        zip_path = prep.zip_bundle(bundle)
        print(f"{region_id}: {report['total_cells']} cells ({report['acting_cells']} acting, "
              f"{report['context_cells']} context) -> {zip_path}")


# --------------------------------------------------------------------------- Parts 2-4
def _settings(a):
    return dict(model=a.model, scenario=a.scenario, seed=a.seed, horizon_hours=a.hours,
                max_decisions=a.max_decisions, max_population=a.max_population)


def cmd_check(a):
    wf = vtx.workflow
    workdir = Path(a.workdir)
    bundle = wf.extract_bundle(a.bundle_zip, workdir / "extracted")
    settings = _settings(a)
    wf.effective_contract(bundle, settings)
    facts, first_request = wf.preflight(bundle, settings)
    from .experiments import experiment_notes
    for note in experiment_notes(settings, facts):
        print(note)
    out = workdir / "preflight"
    out.mkdir(parents=True, exist_ok=True)
    wf.write_json(out / "preflight.json", facts)
    wf.write_json(out / "first_request_preview.json", first_request)
    _print_json({k: facts[k] for k in ("total_cells", "acting_cells", "context_cells", "acting_composition",
                                       "initial_cell_questions", "initial_http_batches", "growth_headroom",
                                       "static_preflight_passed")})
    print("Verified bundle:", bundle["atlas"]["patch_id"], "| ZIP SHA-256:", bundle["zip_sha256"])
    print("Preflight files:", out)


def cmd_run(a):
    from .colab import get_jev_key
    wf = vtx.workflow
    workdir = Path(a.results_dir)
    bundle = wf.extract_bundle(a.bundle_zip, workdir / "extracted")
    settings = _settings(a)
    wf.effective_contract(bundle, settings)
    from .experiments import experiment_notes
    for note in experiment_notes(settings):
        print(note)
    mode = "jev" if a.mode == "jev" else "local_preview"
    run_dir = workdir / (a.run_name + ("_local_preview" if mode == "local_preview" else ""))
    key = ""
    if mode == "jev":
        facts, _ = wf.preflight(bundle, settings)
        print(f"First round: {facts['initial_cell_questions']} cell questions in {facts['initial_http_batches']} "
              f"HTTP requests. Allowance for this run: {a.max_http_attempts} requests.")
        if not a.yes:
            answer = input("Send Jev requests now? Type 'yes' to continue: ")
            if answer.strip().lower() != "yes":
                sys.exit("Cancelled. No request was sent.")
        key = get_jev_key(prompt=sys.stdin.isatty())
        if not key:
            sys.exit("No API key found (set TYPESAFE_API_KEY or JEV_API_KEY). No request was sent.")
    try:
        result = wf.run(bundle, settings, run_dir, api_key=key, max_http_attempts=a.max_http_attempts, mode=mode)
    finally:
        key = ""
    _print_json({k: result["report"][k] for k in ("completed", "stop_reason", "checkpoint_time_min", "usage", "provider_error")})
    print("Replay:", result["folder"] / "replay.html")
    print("Results ZIP:", result["zip_path"])


def cmd_explore(a):
    from . import results as R
    res = R.load_results(a.results)
    _print_json(R.summary(res))
    import pandas as pd
    with pd.option_context("display.width", 200, "display.max_columns", 40):
        print("\nLiving acting cells per type (first and last saved frame):")
        pop = R.population_by_type(res)
        print(pop.iloc[[0, -1]].to_string())
        print("\nChoices by cell type:")
        print(R.choice_counts(res).to_string())
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        R.plot_overview(res).savefig(a.plot, dpi=150)
        print("Overview figure:", a.plot)
    print("Replay:", res["replay"])
    if a.open:
        import webbrowser
        webbrowser.open(res["replay"].resolve().as_uri())


# --------------------------------------------------------------------------- parser
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="vtx", description="VirtualTissue-Xenium command line (Parts 1-4).")
    p.add_argument("--version", action="version", version=f"virtualtissue-xenium {vtx.__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("download-example", help="Download the public 10x Lung 5K sample (~503 MiB).")
    s.add_argument("data_dir"); s.set_defaults(func=cmd_download_example)
    s = sub.add_parser("register-download", help="Register files you extracted yourself from the full 10x ZIP.")
    s.add_argument("data_dir"); s.set_defaults(func=cmd_register_download)
    s = sub.add_parser("verify", help="Check a data folder against its DOWNLOAD_MANIFEST.json.")
    s.add_argument("data_dir"); s.add_argument("--quick", action="store_true", help="Sizes only, no SHA-256.")
    s.set_defaults(func=cmd_verify)
    s = sub.add_parser("own-data", help="Create a Part 1 data folder from your own Xenium output.")
    s.add_argument("xenium_outs"); s.add_argument("data_dir"); s.add_argument("--sample-id", required=True)
    s.add_argument("--include-optional", action="store_true"); s.add_argument("--source-page")
    s.add_argument("--license"); s.add_argument("--attribution"); s.set_defaults(func=cmd_own_data)
    s = sub.add_parser("annotation-template", help="Draft a cluster annotation map for your sample.")
    s.add_argument("data_dir"); s.add_argument("out"); s.add_argument("--top-n", type=int, default=15)
    s.set_defaults(func=cmd_annotation_template)
    s = sub.add_parser("check-annotations", help="Validate a cluster annotation map.")
    s.add_argument("annotations"); s.add_argument("--data-dir"); s.set_defaults(func=cmd_check_annotations)

    s = sub.add_parser("prepare", help="Part 1: AUCell, candidate regions and bundle export.")
    s.add_argument("data_dir", help="Folder with raw/, SOURCE.json and DOWNLOAD_MANIFEST.json.")
    s.add_argument("output_base", help="Results folder (outside DATA_DIR); a run_<timestamp> subfolder is created.")
    s.add_argument("--region", action="append", help="Region to export, e.g. R03 (repeat for several).")
    s.add_argument("--center", nargs=2, type=float, metavar=("X_UM", "Y_UM"), help="Also export a manual region centred here.")
    s.add_argument("--list-only", action="store_true", help="Propose regions and save maps; export nothing.")
    s.add_argument("--tumor-mode", default="exploratory_candidates", choices=["exploratory_candidates", "context_until_review"])
    s.add_argument("--annotations", help="Your cluster annotation map (required for samples other than the public example).")
    s.add_argument("--core-target", type=int, default=750); s.add_argument("--max-total", type=int, default=1000)
    s.add_argument("--halo", type=float, default=30); s.add_argument("--n-regions", type=int, default=8)
    s.add_argument("--aucell-top-ranks", type=int, default=100)
    s.add_argument("--no-verify-hashes", action="store_true", help="Skip SHA-256 checks of the raw files (faster).")
    s.set_defaults(func=cmd_prepare)

    def sim_options(s):
        s.add_argument("--scenario", default="baseline_unknown", choices=SCENARIOS)
        s.add_argument("--model", default="jev-1.14.0"); s.add_argument("--hours", type=float, default=24)
        s.add_argument("--seed", type=int, default=21); s.add_argument("--max-population", type=int, default=1000)
        s.add_argument("--max-decisions", type=int, default=10000)

    s = sub.add_parser("check", help="Part 2: verify a bundle ZIP and preview the first round (no API calls).")
    s.add_argument("bundle_zip"); s.add_argument("--workdir", default="vtx_check"); sim_options(s)
    s.set_defaults(func=cmd_check)
    s = sub.add_parser("run", help="Part 3: run the simulation (Jev or free local preview).")
    s.add_argument("bundle_zip"); s.add_argument("results_dir"); s.add_argument("--run-name", required=True)
    s.add_argument("--mode", default="local", choices=["local", "jev"],
                   help="'local' = programmed rules, no API (default); 'jev' = paid Jev Choice calls.")
    s.add_argument("--max-http-attempts", type=int, default=600)
    s.add_argument("--yes", action="store_true", help="Do not ask for confirmation before sending Jev requests.")
    sim_options(s); s.set_defaults(func=cmd_run)
    s = sub.add_parser("explore", help="Part 4: summarise a results ZIP or run folder.")
    s.add_argument("results"); s.add_argument("--plot", help="Save an overview figure (PNG).")
    s.add_argument("--open", action="store_true", help="Open replay.html in the default browser.")
    s.set_defaults(func=cmd_explore)
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "hours", None) is not None and float(args.hours).is_integer():
        args.hours = int(args.hours)
    from .own_data import PreparationScopeError
    try:
        args.func(args)
    except PreparationScopeError as exc:
        parser.exit(2, f"Preparation stopped: {exc}\n")


if __name__ == "__main__":
    main()
