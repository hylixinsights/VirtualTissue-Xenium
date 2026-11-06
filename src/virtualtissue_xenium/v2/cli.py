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
    if args.polygon_key and context.get("schema_version") != "3.0":
        raise ValueError("--polygon-key requires schema 3.0.")
    from .measurement import read_lineage_map, organism_key
    lineage_map = read_lineage_map(args.lineage_map) if args.lineage_map else None
    measurement = read_json(args.measurement) if args.measurement else None
    annotation_provenance = None
    if context.get("schema_version") == "3.0" and Path(args.input).suffix == ".json":
        source_object = read_json(args.input)
        if not args.measurement:
            measurement = source_object.get("measurement")
        annotation_provenance = source_object.get("annotation_provenance")
        if annotation_provenance is not None and not isinstance(annotation_provenance, dict):
            raise ValueError("annotation_provenance must be an object.")
    if args.lineage_map:
        import hashlib
        annotation_provenance = dict(annotation_provenance or {}, lineage_map_file_sha256=hashlib.sha256(Path(args.lineage_map).read_bytes()).hexdigest())
    if args.platform:
        if context.get("schema_version") != "3.0":
            raise ValueError("--platform requires a 3.0 context.")
        context["platform"] = args.platform
        context["assay"] = args.platform
    source_path = Path(args.input) / "SOURCE.json"
    if source_path.exists():
        source = read_json(source_path)
        if organism_key(source.get("organism")) != organism_key(context.get("organism")):
            raise ValueError("Dataset manifest organism contradicts the context.")
        if context.get("schema_version") == "3.0":
            measurement = dict(measurement or {}, organism=source["organism"], source_manifest_sha256=digest(source))
    native_platform = context.get("platform", context.get("assay"))
    if (context.get("schema_version") == "3.0" and native_platform in ("atera", "visium_hd", "stereo_seq")
            and Path(args.input).suffix not in (".json", ".h5ad")):
        if args.command == "prepare-xenium":
            raise ValueError("Use prepare --platform for non-Xenium native files.")
        from ..v3.adapters import load as load_platform
        cells, genes, measurement = load_platform(native_platform, args.input, args.annotations,
            dict(max_cells=args.max_cells, center=args.center, radius_um=args.radius_um),
            lineage_map=lineage_map, measurement=measurement)
    elif args.command == "prepare-xenium":
        from .datasets import load_xenium_cells
        cells, genes = load_xenium_cells(args.input, args.annotations, max_cells=args.max_cells,
                                       center=args.center, radius_um=args.radius_um, lineage_map=lineage_map)
    else:
        cells, genes = load_cells(args.input, cell_type_key=args.cell_type_key, lineage_key=args.lineage_key,
            sample_key=args.sample_key, layer=args.layer, coordinate_scale=args.coordinate_scale, max_cells=args.max_cells,
            lineage_map=lineage_map, polygon_key=args.polygon_key)
        if len(cells) > args.max_cells:
            raise ValueError("Subset the input first or explicitly increase --max-cells.")
    if args.polygons:
        if context.get("schema_version") != "3.0":
            raise ValueError("Polygon input requires schema 3.0.")
        import hashlib
        from ..v3.object_input import attach_polygons
        cells, geometry_join = attach_polygons(cells, read_json(args.polygons), coordinate_scale=args.polygon_scale)
        measurement = dict(measurement or {}, geometry_input=dict(geometry_join,
            source_file=Path(args.polygons).name, source_sha256=hashlib.sha256(Path(args.polygons).read_bytes()).hexdigest()))
    gene_sets = read_json(args.gene_sets) if args.gene_sets else None
    if isinstance(gene_sets, dict):
        # Keys starting with "_" are file-level notes (provenance, review status), not gene sets.
        gene_sets = {k: v for k, v in gene_sets.items() if not k.startswith("_")}
    bundle = prepare(cells, context, assayed_genes=genes, gene_sets=gene_sets, measurement=measurement,
                     lineage_map=lineage_map, annotation_provenance=annotation_provenance)
    report = preflight(bundle)
    sidecar = Path(args.output).with_suffix(".guidebook.json")
    if Path(args.output).exists() or sidecar.exists():
        raise ValueError("Prepared output or guidebook sidecar exists; choose a new bundle path.")
    write_json(args.output, bundle)
    write_json(sidecar, bundle["guidebook"])
    return dict(prepared=args.output, fingerprint=bundle["fingerprint"], preflight=report)


def cmd_cost(args):
    """Estimate Jev HTTP calls and input tokens from step-0 payloads; no API use."""
    from .provider import estimate_cost
    from .simulation import decision_rows, initialize
    bundle = read_json(args.prepared)
    preflight(bundle, model=args.model, max_questions=args.max_questions, prior_mode=args.jev_prior)
    rows, _ = decision_rows(initialize(bundle), bundle["guidebook"], bundle["context"], seed=args.seed, step=0)
    active = [r for r in rows if len(r["actions"]) > 1]
    return estimate_cost(active, bundle["context"], model=args.model, steps=args.steps, max_questions=args.max_questions,
                         prior_mode=args.jev_prior, max_http_attempts=args.max_http_attempts)


def cmd_run(args):
    from .provider import JevProvider
    bundle = read_json(args.prepared)
    facts = preflight(bundle, model=args.model, max_questions=args.max_questions, seed=args.seed, prior_mode=args.jev_prior)
    effective_weight = 1.0 if args.mode == "jev" and args.jev_prior == "hidden" else args.jev_weight
    if args.mode == "jev":
        key = os.environ.get("TYPESAFE_API_KEY") or os.environ.get("JEV_API_KEY")
        if not key:
            raise ValueError("Set TYPESAFE_API_KEY or JEV_API_KEY before Jev mode.")
        identity = dict(bundle_fingerprint=bundle["fingerprint"], steps=args.steps, seed=args.seed,
                        mode=args.mode, jev_weight=effective_weight, max_population=args.max_population,
                        model=args.model, max_questions=args.max_questions, jev_prior=args.jev_prior)
        provider = JevProvider(key, Path(args.output)/"jev_journal", identity, model=args.model,
            max_attempts=args.max_http_attempts, workers=args.workers, max_questions=args.max_questions,
            prior_mode=args.jev_prior, schema_version=bundle["schema_version"])
    else:
        provider = None
    report = run(bundle, args.output, steps=args.steps, seed=args.seed, mode=args.mode,
                 provider=provider, jev_weight=effective_weight, max_population=args.max_population, resume=args.resume)
    return dict(**report, replay=str(Path(args.output).resolve()/"replay.html"))


def compare(paths, *, allow_ontology_diff=False):
    """Compare recorded artifacts; an ontology override never licenses a new run."""
    import math
    from .measurement import comparison_signature
    from .readouts import iter_frames

    if len(paths) < 2:
        raise ValueError("Compare requires at least two run directories.")
    if len({Path(path).resolve() for path in paths}) != len(paths):
        raise ValueError("Compare requires distinct run directories; duplicate runs would double-count usage.")
    results = []
    for path in paths:
        folder = Path(path)
        config, report = read_json(folder/"config.json"), read_json(folder/"report.json")
        guidebook = read_json(folder/"guidebook.json")
        results.append(dict(path=str(folder), config=config, report=report, guidebook=guidebook))
    ref = results[0]
    ontology_differences = []
    for result in results[1:]:
        reference = ref["guidebook"]["provenance"]
        candidate = result["guidebook"]["provenance"]
        # Use the frozen declarations, never the ontology installed today.
        changed_hashes = []
        left_hashes, right_hashes = reference.get("ontology_hashes") or {}, candidate.get("ontology_hashes") or {}
        for key in sorted(set(left_hashes) | set(right_hashes)):
            if left_hashes.get(key) != right_hashes.get(key):
                changed_hashes.append(dict(asset=key, reference=left_hashes.get(key), candidate=right_hashes.get(key)))
        changed_metadata = []
        for key in ("ontology_sha256", "kinetics_sha256", "ontology_version"):
            if reference.get(key) != candidate.get(key):
                changed_metadata.append(dict(field=key, reference=reference.get(key), candidate=candidate.get(key)))
        left_version = ref["guidebook"]["context"].get("ontology_version")
        right_version = result["guidebook"]["context"].get("ontology_version")
        if left_version != right_version:
            changed_metadata.append(dict(field="context.ontology_version", reference=left_version, candidate=right_version))
        left_rows, right_rows = reference.get("lineage_rules"), candidate.get("lineage_rules")
        row_differences = []
        if isinstance(left_rows, dict) and isinstance(right_rows, dict):
            for name in sorted(set(left_rows) | set(right_rows)):
                if left_rows.get(name) != right_rows.get(name):
                    row_differences.append(dict(lineage=name, reference=left_rows.get(name), candidate=right_rows.get(name)))
        elif left_rows != right_rows:
            changed_metadata.append(dict(field="lineage_rules_snapshot", reference_available=isinstance(left_rows, dict),
                                         candidate_available=isinstance(right_rows, dict)))
        if changed_hashes or changed_metadata or row_differences:
            if not allow_ontology_diff:
                raise ValueError("Unmatched ontology hashes or recorded declarations; --allow-ontology-diff permits forensic artifact comparison only.")
            ontology_differences.append(dict(reference_path=ref["path"], candidate_path=result["path"],
                changed_asset_hashes=changed_hashes, changed_metadata=changed_metadata,
                recorded_lineage_row_differences=row_differences,
                lineage_rows_available=isinstance(left_rows, dict) and isinstance(right_rows, dict),
                unavailable_rows="Complete historical kinetics and other asset rows are not stored in run artifacts; only recorded lineage declarations can be compared. No rows are reconstructed from the currently installed ontology."))
        if comparison_signature(ref["config"].get("measurement")) != comparison_signature(result["config"].get("measurement")):
            raise ValueError("Unmatched measurement unit or acquisition semantics.")
        for key in ("seed", "mode", "jev_weight", "model", "kernel", "max_population", "assumptions",
                    "provider_name", "provider_prompt_version", "jev_prior_mode", "max_questions", "provider_sampling_parameters"):
            if ref["config"].get(key) != result["config"].get(key):
                raise ValueError(f"Unmatched run setting: {key}")
        control = result["guidebook"].get("control") or ref["guidebook"].get("control") or {}
        inputs = {reference["input_sha256"], candidate["input_sha256"]}
        if len(inputs) > 1 and control.get("source_input_sha256") not in inputs:
            raise ValueError("Runs have different measured input cells/expression/geometry and no declared control links them.")
        if reference["gene_sets_sha256"] != candidate["gene_sets_sha256"]:
            raise ValueError("Runs have different gene sets.")
        for key in ("guidebook_version", "assayed_genes", "panel_source", "score_method", "symbol_normalization",
                    "default_action_weights", "recognition_licensed_kill_weight", "time_semantics", "kinetics_mode"):
            if reference.get(key) != candidate.get(key):
                raise ValueError(f"Unmatched guidebook scoring assumption: {key}")
        for key in ("time_step_hours", "kinetics", "scoring", "interaction_mode", "neighbor_radius_um", "max_neighbors", "sampling", "organism", "assay", "tissue", "disease", "platform", "dynamics", "kernel", "lattice"):
            if key in ("tissue", "disease") and control.get("kind") == "jev_scrambled_context":
                continue
            if key == "interaction_mode" and control.get("kind") == "well_mixed_equivalent":
                continue
            if key == "scoring" and control.get("kind") == "beta_zero":
                continue
            if ref["guidebook"]["context"].get(key) != result["guidebook"]["context"].get(key):
                raise ValueError(f"Unmatched biological/context setting: {key}")
    shared = min(r["report"]["completed_steps"] for r in results)
    output, usage_rows = [], []
    shared_hours = None
    usage_keys = ("http_attempts", "cache_hits", "input_tokens", "output_tokens", "attempts_with_unknown_usage")
    for r in results:
        frame = None
        gained = lost = 0.0
        for current in iter_frames(r["path"]):
            if current["step"] > shared:
                break
            for event in current.get("events", []):
                if event["kind"] == "occupancy_gain":
                    gained += event["added_occupancy"]
                elif event["kind"] == "occupancy_loss":
                    lost += event["removed_occupancy"]
            if current["step"] == shared:
                frame = current
                break
        if frame is None:
            raise ValueError(f"Run {r['path']} is missing the common completed frame {shared}.")
        if shared_hours is None:
            shared_hours = frame["time_hours"]
        elif not math.isclose(shared_hours, frame["time_hours"], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("Common completed step has unmatched recorded time_hours.")
        cells = list(frame["cells"].values())
        context = r["guidebook"]["context"]
        row = dict(path=r["path"], shared_step=shared, time_hours=frame["time_hours"],
            completed=r["report"]["completed"], stop_reason=r["report"]["stop_reason"],
            context=context, control=r["guidebook"].get("control"),
            jev_prior_divergence=r["report"].get("jev_prior_divergence", {}).get("by_lineage"),
            jev_prior_divergence_scope="entire recorded run, which may exceed the shared endpoint",
            guidebook_fingerprint=digest(r["guidebook"]))
        if context.get("kernel") == "lattice":
            occupied = [c for c in cells if c["occupancy"] > 0]
            by_lineage = {}
            for cell in occupied:
                by_lineage[cell["lineage"]] = by_lineage.get(cell["lineage"], 0.0) + cell["occupancy"]
            row.update(kernel="lattice", occupied_bins=len(occupied), empty_bins=len(cells)-len(occupied),
                activated_occupied_bins=sum(c["activated"] for c in occupied),
                exhausted_occupied_bins=sum(c["exhausted"] for c in occupied),
                occupancy=dict(total=sum(by_lineage.values()), by_lineage=dict(sorted(by_lineage.items())),
                    capacity_per_bin=context["lattice"]["capacity"],
                    unit="abstract occupancy units; a full bin has the declared capacity, not a cell count"),
                cumulative_occupancy_gained=gained, cumulative_occupancy_lost=lost)
        else:
            row.update(living=sum(c["status"] == "alive" for c in cells), dead=sum(c["status"] == "dead" for c in cells),
                activated=sum(c["status"] == "alive" and c["activated"] for c in cells),
                exhausted=sum(c["status"] == "alive" and c["exhausted"] for c in cells))
        recorded = r["report"].get("usage") or {}
        usage = {key: recorded.get(key, 0 if r["config"].get("mode") == "local" else None) for key in usage_keys}
        if any(value is not None and (type(value) is not int or value < 0) for value in usage.values()):
            raise ValueError(f"Invalid recorded provider usage in {r['path']}.")
        usage_rows.append(dict(path=r["path"], **usage))
        row["usage"] = usage
        output.append(row)
    total = {key: sum(row[key] for row in usage_rows if row[key] is not None) for key in usage_keys}
    missing = {key: sum(row[key] is None for row in usage_rows) for key in usage_keys}
    ledger = dict(by_run=usage_rows, known_totals=total, runs_missing_fields=missing, currency_cost=None,
        token_totals_complete=not any(missing[key] for key in ("input_tokens", "output_tokens", "attempts_with_unknown_usage")) and total["attempts_with_unknown_usage"] == 0,
        scope="Actual recorded usage for entire runs, including attempts beyond the shared endpoint; no estimated tokens or currency prices. Missing fields and unknown-usage attempts are not counted as known tokens.")
    return dict(shared_step=shared, time_hours=shared_hours, runs=output, cost_ledger=ledger,
        allow_ontology_diff=allow_ontology_diff, ontology_differences=ontology_differences,
        interpretation=("Forensic artifact comparison with changed ontology declarations; differences confound the outcomes and are not a matched experiment. No causal treatment efficacy claim."
                        if ontology_differences else
                        "Matched stochastic assumptions experiment; no causal treatment efficacy claim. Repeat across seeds."))

def cmd_control(args):
    from .controls import CONTROLS, apply_control
    if args.kind == "list":
        return CONTROLS
    if not args.output:
        raise ValueError("control requires an OUTPUT path for the transformed bundle.")
    bundle = apply_control(read_json(args.prepared), args.kind, seed=args.seed,
                           donor_context=read_json(args.donor_context) if args.donor_context else None)
    if Path(args.output).exists():
        raise ValueError("Control output exists; choose a new bundle path.")
    write_json(args.output, bundle)
    return dict(prepared=args.output, fingerprint=bundle["fingerprint"], control=bundle["guidebook"]["control"],
                preflight=preflight(bundle))


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


def cmd_ensemble(args):
    """Paths inside an arm manifest are resolved relative to that manifest."""
    from .ensemble import run_ensemble
    from .provider import JevProvider
    manifest_path = Path(args.manifest).resolve()
    manifest = read_json(manifest_path)
    if not isinstance(manifest, dict) or not isinstance(manifest.get("arms"), dict):
        raise ValueError("Ensemble manifest requires an arms map: name -> prepared path or object.")
    arms, options, prior_modes = {}, {}, {}
    for name, spec in manifest["arms"].items():
        if isinstance(spec, str):
            spec = {"prepared": spec}
        if not isinstance(spec, dict) or "prepared" not in spec or set(spec) - {"prepared", "mode", "jev_weight", "jev_prior"}:
            raise ValueError("Each ensemble arm declares prepared, and optionally mode, jev_weight, jev_prior.")
        arms[name] = str((manifest_path.parent / spec["prepared"]).resolve())
        mode = spec.get("mode", "local")
        prior = spec.get("jev_prior", "shown")
        if prior not in ("shown", "hidden"):
            raise ValueError("jev_prior must be shown or hidden.")
        options[name] = {"mode": mode, "jev_weight": 1.0 if prior == "hidden" else spec.get("jev_weight", .5)}
        prior_modes[name] = prior
    key = os.environ.get("TYPESAFE_API_KEY") or os.environ.get("JEV_API_KEY")
    if any(o["mode"] == "jev" for o in options.values()) and not key:
        raise ValueError("Set TYPESAFE_API_KEY or JEV_API_KEY for declared Jev ensemble arms.")

    def provider_factory(*, arm, bundle, output_dir, seed, steps, jev_weight):
        identity = dict(bundle_fingerprint=bundle["fingerprint"], arm=arm, seed=seed, steps=steps,
                        mode="jev", jev_weight=jev_weight, max_population=args.max_population,
                        model=args.model, prior_mode=prior_modes[arm])
        return JevProvider(key, Path(output_dir)/"jev_journal", identity, model=args.model,
                           max_attempts=args.max_http_attempts, workers=args.workers,
                           prior_mode=prior_modes[arm], schema_version=bundle["schema_version"])
    return run_ensemble(arms, args.output, seeds=args.seeds, steps=args.steps,
                        horizon_hours=args.horizon_hours, max_population=args.max_population,
                        bootstrap_replicates=args.bootstrap_replicates, bootstrap_seed=args.bootstrap_seed,
                        arm_options=options, provider_factory=provider_factory if any(o["mode"] == "jev" for o in options.values()) else None,
                        readouts_path=args.readouts)



def cmd_sensitivity(args):
    from .sensitivity import run_sensitivity
    return run_sensitivity(args.prepared, read_json(args.bounds), args.output,
        horizon_hours=args.horizon_hours, seeds=args.seeds, n=args.n, method=args.method,
        design_seed=args.design_seed, readout=args.readout, lineage=args.lineage,
        max_population=args.max_population, bootstrap_replicates=args.bootstrap_replicates)


def cmd_roi_candidates(args):
    from ..v3.annotation import select_rois
    obj = read_json(args.input)
    result = select_rois(obj["cells"], cells_per_roi=args.cells_per_roi, grid_um=args.grid_um,
                         max_candidates=args.max_candidates, contact_radius_um=args.contact_radius_um,
                         core_max_immune_fraction=args.core_max_immune_fraction,
                         core_min_malignant_fraction=args.core_min_malignant_fraction)
    if Path(args.output).exists():
        raise ValueError("ROI candidate output exists; choose a new path.")
    write_json(args.output, result)
    return result


def cmd_annotation_agreement(args):
    from ..v3.annotation import write_annotation_comparison
    return write_annotation_comparison(args.first, args.second, args.output, id_scope=args.id_scope)


def cmd_platform_template(args):
    from ..v3.adapters import write_template
    measurement = read_json(args.measurement) if args.measurement else None
    if measurement and measurement.get("centroids_csv"):
        measurement["centroids_csv"] = str((Path(args.measurement).resolve().parent / measurement["centroids_csv"]).resolve())
    path = write_template(args.platform, args.input, args.output,
                          dict(max_cells=args.max_cells, center=args.center, radius_um=args.radius_um),
                          measurement=measurement)
    return {"template": str(path), "annotation_status": "blank; requires investigator review"}


def _roi_source(args):
    """Apply an explicit label-to-permission mapping before hashing a raw object."""
    source = read_json(args.input)
    if args.lineage_map:
        if "guidebook" in source or "fingerprint" in source:
            raise ValueError("A prepared source already freezes its lineage mapping; map the raw source instead.")
        from .measurement import apply_lineage_map, read_lineage_map
        mapping = read_lineage_map(args.lineage_map)
        source["cells"] = apply_lineage_map(source["cells"], mapping)
        source["annotation_provenance"] = dict(source.get("annotation_provenance") or {},
            lineage_map=mapping, lineage_map_sha256=digest(mapping))
    return source


def cmd_roi_catalog(args):
    from ..v3.roi_catalog import write_roi_catalog
    return write_roi_catalog(_roi_source(args), args.output, cells_per_roi=args.cells_per_roi,
                             grid_um=args.grid_um, max_candidates=args.max_candidates)


def cmd_select_roi(args):
    from copy import deepcopy
    from ..v3.roi_catalog import extract_roi_selection, validate_roi_selection
    source, selection = _roi_source(args), read_json(args.selection)
    selection = validate_roi_selection(source, selection)
    selected = extract_roi_selection(source, selection)
    is_prepared = "guidebook" in source and "fingerprint" in source
    if args.context:
        context = read_json(args.context)
    elif is_prepared:
        context = deepcopy(source["context"])
    else:
        raise ValueError("Selecting from raw JSON requires --context; only the chosen ROI will be prepared.")
    if context.get("schema_version") != "3.0":
        raise ValueError("User-selected ROI bundles require schema 3.0.")
    source_sampling = deepcopy(context.get("sampling", {}))
    missing = list(dict.fromkeys([*source_sampling.get("missing_populations", []), "outside the user-selected ROI"]))
    context["sampling"] = dict(design="subsampled", description="ROI explicitly chosen by the user from a recorded composition catalog.",
                               missing_populations=missing, source_sampling=source_sampling)
    if is_prepared:
        measurement = deepcopy(source["measurement"])
        gene_sets = source["guidebook"]["provenance"]["gene_sets"]
        genes = source["guidebook"]["provenance"]["assayed_genes"]
        annotation = deepcopy(source.get("annotation_provenance", {}))
    else:
        measurement = read_json(args.measurement) if args.measurement else deepcopy(source.get("measurement"))
        gene_sets = read_json(args.gene_sets) if args.gene_sets else None
        genes = source.get("assayed_genes")
        annotation = deepcopy(source.get("annotation_provenance", {}))
        annotation.setdefault("source", "user_supplied_object")
        annotation.setdefault("verification", "not_independently_verified")
    if measurement:
        measurement.pop("mean_transcripts_per_unit", None)
    annotation["roi_selection"] = selection
    bundle = prepare(selected, context, assayed_genes=genes, gene_sets=gene_sets,
                     measurement=measurement, annotation_provenance=annotation)
    preflight_result = preflight(bundle)
    sidecar = Path(args.output).with_suffix(".guidebook.json")
    if Path(args.output).exists() or sidecar.exists():
        raise ValueError("Selected ROI output exists; choose a new path.")
    write_json(args.output, bundle)
    write_json(sidecar, bundle["guidebook"])
    return dict(prepared=args.output, fingerprint=bundle["fingerprint"], selection=selection,
                preflight=preflight_result, scoring="Recomputed for the selected ROI; catalog coverage is not silently reused.")


def cmd_expression_audit(args):
    from ..v3.expression_audit import audit_expression
    if Path(args.output).exists():
        raise ValueError("Expression audit output exists; choose a new path.")
    result = audit_expression(read_json(args.prepared), seed=args.seed)
    write_json(args.output, result)
    return result


def cmd_replay(args):
    """Refresh presentation from saved snapshots; never rerun a decision round."""
    from .replay import write_replay, write_replay_stream

    directory = Path(args.run_directory)
    report = read_json(directory / "report.json")
    if not isinstance(report, dict) or report.get("schema_version") not in {"2.0", "3.0"}:
        raise ValueError("Replay refresh requires a saved schema-2.0 or schema-3.0 run report.")
    output = Path(args.output) if args.output else directory / "replay.html"
    if output.suffix.lower() != ".html":
        raise ValueError("Replay output must end in .html; scientific exports are not replay destinations.")
    if output.is_symlink():
        raise ValueError("Replay output must not be a symbolic link to another saved artifact.")
    if output.exists() and not output.is_file():
        raise ValueError("Replay output must be an HTML file path, not an existing directory.")
    protected = [path for path in directory.iterdir() if path.is_file() and path.suffix.lower() != ".html"]
    if output.exists() and any(output.samefile(path) for path in protected):
        raise ValueError("Replay output collides with a saved artifact; choose another .html path.")
    stream = directory / "frames.jsonl"
    array = directory / "frames.json"
    if report.get("schema_version") == "3.0":
        if not stream.is_file():
            raise ValueError("A v3 replay requires the saved frames.jsonl stream.")
        sidecars = output.with_name(output.stem + "_frames")
        if sidecars.is_symlink() or (sidecars.exists() and any(path.is_symlink() for path in sidecars.iterdir())):
            raise ValueError("Replay sidecars must not be symbolic links to other saved artifacts.")
        output.parent.mkdir(parents=True, exist_ok=True)
        write_replay_stream(output, stream, report)
    else:
        frames = read_json(array)
        output.parent.mkdir(parents=True, exist_ok=True)
        write_replay(output, frames, report)
    return dict(replay=str(output.resolve()), source=str(directory.resolve()),
                refreshed_from_saved_frames=True)


def parser():
    p = argparse.ArgumentParser(prog="vtx2", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("prepare", "prepare-xenium"):
        s = sub.add_parser(name, help="Build a reviewed context-conditioned guidebook and prepared JSON.")
        s.add_argument("input"); s.add_argument("context"); s.add_argument("output")
        s.add_argument("--gene-sets"); s.add_argument("--max-cells", type=int, default=1000)
        s.add_argument("--lineage-map"); s.add_argument("--measurement")
        s.add_argument("--polygons", help="Exact-ID polygon map or GeoJSON FeatureCollection; schema 3.0 only.")
        s.add_argument("--polygon-scale", type=float, default=1.0, help="Convert external polygon coordinates to micrometres.")
        s.set_defaults(polygon_key=None)
        s.add_argument("--platform", choices=("xenium", "scrna", "snrna", "atera", "visium_hd", "stereo_seq"))
        if name == "prepare-xenium":
            s.add_argument("--annotations", required=True)
            s.add_argument("--center", nargs=2, type=float); s.add_argument("--radius-um", type=float)
        else:
            s.add_argument("--cell-type-key", default="cell_type"); s.add_argument("--lineage-key", default="lineage")
            s.add_argument("--sample-key"); s.add_argument("--layer"); s.add_argument("--coordinate-scale", type=float, default=1)
            s.add_argument("--annotations")
            s.add_argument("--polygon-key", help="AnnData.uns/obs key containing polygons; uses coordinate-scale.")
            s.add_argument("--center", nargs=2, type=float); s.add_argument("--radius-um", type=float)
        s.set_defaults(func=cmd_prepare)
    s = sub.add_parser("inspect", help="Validate a prepared bundle; estimate first-round calls, no API use.")
    s.add_argument("prepared"); s.set_defaults(func=lambda a: preflight(read_json(a.prepared)))
    s = sub.add_parser("run", help="Synchronous local or Jev simulation; local is the default.")
    s.add_argument("prepared"); s.add_argument("output"); s.add_argument("--mode", choices=("local", "jev"), default="local")
    s.add_argument("--steps", type=int, default=12); s.add_argument("--seed", type=int, default=21)
    s.add_argument("--max-population", type=int, default=2000); s.add_argument("--jev-weight", type=float, default=.5)
    s.add_argument("--model", default="jev-1.14.0"); s.add_argument("--max-http-attempts", type=int, default=100)
    s.add_argument("--workers", type=int, default=2); s.add_argument("--max-questions", type=int, default=20)
    s.add_argument("--jev-prior", choices=("shown", "hidden"), default="shown",
                   help="hidden omits the guidebook prior from the prompt (anchoring control); journals never mix modes.")
    s.add_argument("--resume", action="store_true"); s.set_defaults(func=cmd_run)
    s = sub.add_parser("replay", help="Refresh a saved run's offline replay without rerunning the simulation or calling Jev.")
    s.add_argument("run_directory", help="Saved run folder containing report.json and frames.jsonl (v3) or frames.json (v2).")
    s.add_argument("--output", help="Replay HTML destination; defaults to RUN_DIRECTORY/replay.html.")
    s.set_defaults(func=cmd_replay)
    s = sub.add_parser("cost", help="Estimate Jev HTTP calls and input tokens for a run; no API use.")
    s.add_argument("prepared"); s.add_argument("--steps", type=int, default=12); s.add_argument("--seed", type=int, default=21)
    s.add_argument("--model", default="jev-1.14.0"); s.add_argument("--max-questions", type=int, default=20)
    s.add_argument("--max-http-attempts", type=int, default=100)
    s.add_argument("--jev-prior", choices=("shown", "hidden"), default="shown"); s.set_defaults(func=cmd_cost)
    s = sub.add_parser("compare", help="Compare matched arms at their latest common completed step.")
    s.add_argument("runs", nargs="+")
    s.add_argument("--allow-ontology-diff", action="store_true",
                   help="Forensic comparison of recorded ontology changes; does not permit run/resume with stale ontology.")
    s.set_defaults(func=lambda a: compare(a.runs, allow_ontology_diff=a.allow_ontology_diff))
    s = sub.add_parser("control", help="Write a null-model control bundle (shuffled expression, labels, coordinates, ...).")
    s.add_argument("prepared"); s.add_argument("output", nargs="?")
    s.add_argument("--kind", required=True, help="Control name, or 'list' to print the available controls.")
    s.add_argument("--seed", type=int, default=0); s.add_argument("--donor-context")
    s.set_defaults(func=cmd_control)
    s = sub.add_parser("dataset", help="Verified public-data manifests and bounded downloads.")
    s.add_argument("action", choices=("list", "download", "verify")); s.add_argument("dataset", nargs="?")
    s.add_argument("output", nargs="?"); s.add_argument("--max-mib", type=int, default=100)
    s.add_argument("--include-optional", action="store_true"); s.set_defaults(func=cmd_dataset)
    s = sub.add_parser("annotation-template", help="Blank reviewed labels with exact Xenium IDs/vendor clusters.")
    s.add_argument("input"); s.add_argument("output"); s.add_argument("--max-cells", type=int, default=1000)
    s.add_argument("--center", nargs=2, type=float); s.add_argument("--radius-um", type=float)
    s.set_defaults(func=cmd_annotation_template)
    s = sub.add_parser("ensemble", help="Run matched arms across seeds and aggregate predefined readouts.")
    s.add_argument("manifest"); s.add_argument("output")
    time = s.add_mutually_exclusive_group(required=True)
    time.add_argument("--steps", type=int); time.add_argument("--horizon-hours", type=float)
    s.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    s.add_argument("--bootstrap-replicates", type=int, default=2000)
    s.add_argument("--bootstrap-seed", type=int, default=0); s.add_argument("--readouts")
    s.add_argument("--max-population", type=int, default=2000)
    s.add_argument("--model", default="jev-1.14.0")
    s.add_argument("--max-http-attempts", type=int, default=100)
    s.add_argument("--workers", type=int, default=2)
    s.set_defaults(func=cmd_ensemble)
    s = sub.add_parser("sensitivity", help="Local parameter sweep with LHS or a valid pick-freeze Sobol design.")
    s.add_argument("prepared"); s.add_argument("bounds"); s.add_argument("output")
    s.add_argument("--horizon-hours", type=float, required=True)
    s.add_argument("--seeds", type=int, nargs="+", default=[0,1,2]); s.add_argument("--n", type=int, default=200)
    s.add_argument("--method", choices=("lhs","sobol"), default="sobol")
    s.add_argument("--design-seed", type=int, default=0); s.add_argument("--readout", default="living_cells")
    s.add_argument("--lineage", default="malignant"); s.add_argument("--max-population", type=int, default=2000)
    s.add_argument("--bootstrap-replicates", type=int, default=1000); s.set_defaults(func=cmd_sensitivity)
    s = sub.add_parser("roi-candidates", help="Select reproducible candidate ROIs from reviewed labels and measured coordinates.")
    s.add_argument("input"); s.add_argument("output"); s.add_argument("--cells-per-roi", type=int, default=1000)
    s.add_argument("--grid-um", type=float, default=500); s.add_argument("--max-candidates", type=int, default=50)
    s.add_argument("--contact-radius-um", type=float, default=30)
    s.add_argument("--core-max-immune-fraction", type=float, default=.05)
    s.add_argument("--core-min-malignant-fraction", type=float, default=.5); s.set_defaults(func=cmd_roi_candidates)
    s = sub.add_parser("annotation-agreement", help="Compare independent annotations by exact cell ID and report Cohen kappa.")
    s.add_argument("first"); s.add_argument("second"); s.add_argument("output")
    s.add_argument("--id-scope", choices=("exact","intersection"), default="exact"); s.set_defaults(func=cmd_annotation_agreement)
    s = sub.add_parser("jev-audit", help="Measure repeated-payload and batch-order variance from existing journals; no API calls.")
    s.add_argument("journals", nargs="+"); s.add_argument("--output", required=True)
    s.set_defaults(func=lambda a: __import__("virtualtissue_xenium.v2.provider_audit", fromlist=["audit_journals"]).audit_journals(a.journals, a.output))
    s = sub.add_parser("platform-template", help="Export blank annotation rows from native Atera, Visium HD or Stereo-seq geometry.")
    s.add_argument("platform", choices=("atera", "visium_hd", "stereo_seq"))
    s.add_argument("input"); s.add_argument("output"); s.add_argument("--measurement")
    s.add_argument("--max-cells", type=int, default=1000); s.add_argument("--center", nargs=2, type=float)
    s.add_argument("--radius-um", type=float); s.set_defaults(func=cmd_platform_template)
    s = sub.add_parser("roi-catalog", help="Present composition and measured geometry; the user chooses a region offline.")
    s.add_argument("input"); s.add_argument("output")
    s.add_argument("--lineage-map", help="Explicit CSV mapping existing cell_type labels to model lineages; use the same map for select-roi.")
    s.add_argument("--cells-per-roi", type=int, default=1000); s.add_argument("--grid-um", type=float, default=500)
    s.add_argument("--max-candidates", type=int, default=50); s.set_defaults(func=cmd_roi_catalog)
    s = sub.add_parser("select-roi", help="Verify an explicit catalog selection and prepare only the selected ROI.")
    s.add_argument("input"); s.add_argument("selection"); s.add_argument("output")
    s.add_argument("--context"); s.add_argument("--measurement"); s.add_argument("--gene-sets")
    s.add_argument("--lineage-map", help="Use the same label mapping used to generate the raw-object catalogue.")
    s.set_defaults(func=cmd_select_roi)
    s = sub.add_parser("expression-audit", help="Measure fixed-state expression effects versus label-only and shuffle controls; no API calls.")
    s.add_argument("prepared"); s.add_argument("output"); s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=cmd_expression_audit)
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
