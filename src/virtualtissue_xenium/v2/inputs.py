"""Strict adapters for annotated AnnData and a small, transparent JSON format."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import tempfile


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(value, indent=2, allow_nan=False).encode() + b"\n"
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def validate_cells(cells):
    if not isinstance(cells, list) or not cells:
        raise ValueError("Provide a nonempty list of annotated cells.")
    result, seen = [], set()
    for record in cells:
        c = dict(record)
        for key in ("cell_id", "cell_type", "lineage"):
            if not isinstance(c.get(key), str) or not c[key].strip():
                raise ValueError(f"Every cell needs a nonempty string {key}.")
        if c["cell_id"] in seen:
            raise ValueError("Cell IDs must be unique; prefix by sample when merging samples.")
        seen.add(c["cell_id"])
        if not isinstance(c.get("expression"), dict):
            raise ValueError("Each cell needs an expression dictionary, even if empty.")
        for gene, value in c["expression"].items():
            if not isinstance(gene, str) or not gene or type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError("Expression must map gene symbols to finite nonnegative counts.")
        if ("x" in c) != ("y" in c):
            raise ValueError("Both x and y coordinates are required together.")
        for key in ("x", "y"):
            if key in c and (type(c[key]) not in (int, float) or not math.isfinite(c[key])):
                raise ValueError("Coordinates must be finite numbers in micrometres.")
        c.setdefault("sample_id", "sample")
        for key in ("purity", "occupancy"):
            if key in c and (type(c[key]) not in (int, float) or not math.isfinite(c[key]) or c[key] < 0):
                raise ValueError(f"{key} must be a finite nonnegative number.")
        if "purity" in c and c["purity"] > 1:
            raise ValueError("purity must be between zero and one.")
        if not isinstance(c["sample_id"], str) or not c["sample_id"]:
            raise ValueError("sample_id must be a nonempty string.")
        initial = c.get("initial_state", {})
        if not isinstance(initial, dict) or set(initial) - {"activated", "exhausted"} or any(type(x) is not bool for x in initial.values()):
            raise ValueError("initial_state accepts only boolean activated/exhausted assumptions.")
        if "measured_scores" in c:
            scores = c["measured_scores"]
            if not isinstance(scores, dict) or any(
                    not isinstance(name, str) or not name
                    or (value is not None and (type(value) not in (int, float) or not math.isfinite(value)))
                    for name, value in scores.items()):
                raise ValueError("measured_scores must map score names to finite numbers, or null when not calculated.")
        result.append(c)
    return sorted(result, key=lambda c: c["cell_id"])


def from_anndata(adata, *, cell_type_key="cell_type", lineage_key="lineage", sample_key=None,
                 layer=None, spatial_key="spatial", coordinate_scale=1.0, max_cells=1000, lineage_map=None, polygon_key=None,
                 score_keys=None):
    """Accept an in-memory annotated object. Never infer lineage, malignancy, or XY.

    X (or the chosen layer) must hold nonnegative counts, with unique gene symbols
    in var_names. No densification of the full matrix is performed.

    ``score_keys`` names numeric ``obs`` columns holding per-cell scores computed
    outside the simulator. NaN means not calculated and is kept as null, never
    as zero. Their meaning comes from the measured-evidence manifest, not here.
    """
    score_keys = list(score_keys or [])
    if len(set(score_keys)) != len(score_keys):
        raise ValueError("score_keys must be unique obs column names.")
    for key in score_keys:
        if key not in adata.obs:
            raise ValueError(f"Missing measured-score column: {key}")
    if adata.n_obs > max_cells:
        raise ValueError(f"Select a region/subset first: {adata.n_obs} cells exceeds max_cells={max_cells}.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("AnnData cell IDs and gene symbols must be unique.")
    if sample_key is None and "sample_id" in adata.obs:
        sample_key = "sample_id"
    for key in (cell_type_key, *([] if lineage_map else [lineage_key]), *([sample_key] if sample_key else [])):
        if key not in adata.obs or adata.obs[key].isna().any():
            raise ValueError(f"Missing annotation column or values: {key}")
    if not math.isfinite(coordinate_scale) or coordinate_scale <= 0:
        raise ValueError("coordinate_scale must be positive (conversion to micrometres).")
    matrix = adata.layers[layer] if layer else adata.X
    if matrix is None:
        raise ValueError("The selected expression matrix is empty.")
    coords = adata.obsm.get(spatial_key)
    if coords is not None and (len(coords.shape) != 2 or coords.shape[1] < 2):
        raise ValueError("Spatial coordinates need at least two columns.")
    genes, cells = list(map(str, adata.var_names)), []
    for i, cid in enumerate(adata.obs_names):
        row = matrix[i:i+1]
        values = row.toarray().ravel() if hasattr(row, "toarray") else row.reshape(-1)
        # Validate ALL features, including zeros. NaN must never silently disappear.
        if any(not math.isfinite(float(v)) or float(v) < 0 for v in values):
            raise ValueError("Use nonnegative finite counts, not scaled/residual expression.")
        c = dict(cell_id=str(cid), cell_type=str(adata.obs.iloc[i][cell_type_key]),
                 lineage=(lineage_map.get(str(adata.obs.iloc[i][cell_type_key]), "") if lineage_map else str(adata.obs.iloc[i][lineage_key])),
                 expression={g: float(v) for g, v in zip(genes, values) if v != 0},
                 sample_id=str(adata.obs.iloc[i][sample_key]) if sample_key else "sample")
        if coords is not None:
            c.update(x=float(coords[i, 0]) * coordinate_scale, y=float(coords[i, 1]) * coordinate_scale)
        if score_keys:
            scores = {}
            for key in score_keys:
                value = adata.obs.iloc[i][key]
                if type(value).__name__ in ("bool", "bool_"):
                    raise ValueError(f"Measured-score column {key} must be numeric, not boolean.")
                if value is None or type(value).__name__ in ("NAType", "NaTType"):
                    scores[key] = None  # pandas missing value: not calculated
                    continue
                try:
                    value = float(value)
                except (TypeError, ValueError):
                    raise ValueError(f"Measured-score column {key} must be numeric.") from None
                if math.isinf(value):
                    raise ValueError(f"Measured-score column {key} contains an infinite value.")
                scores[key] = None if math.isnan(value) else value
            c["measured_scores"] = scores
        cells.append(c)
    if polygon_key:
        from ..v3.object_input import attach_anndata_polygons
        cells, _ = attach_anndata_polygons(cells, adata, polygon_key, coordinate_scale=coordinate_scale)
    return validate_cells(cells), genes


def load_cells(path, **kwargs):
    path = Path(path)
    if path.suffix == ".h5ad":
        # read_h5ad(backed=...) may eagerly load layers. Enforce the bound before opening AnnData.
        import h5py
        with h5py.File(path, "r") as handle:
            obs = handle["obs"]
            index = obs.attrs.get("_index", "_index")
            if isinstance(index, bytes):
                index = index.decode()
            if index not in obs:
                raise ValueError("Cannot verify AnnData observation count before loading; export a modern H5AD first.")
            n_obs = obs[index].shape[0]
            if n_obs > kwargs.get("max_cells", 1000):
                raise ValueError(f"Select a region/subset first: {n_obs} cells exceeds max_cells={kwargs.get('max_cells', 1000)}.")
        import anndata
        obj = anndata.read_h5ad(path, backed="r")
        try:
            return from_anndata(obj, **kwargs)
        finally:
            obj.file.close()
    obj = read_json(path)
    if not isinstance(obj, dict) or set(obj) - {"cells", "assayed_genes", "description", "measurement", "annotation_provenance", "context", "measured_evidence"}:
        raise ValueError("JSON input requires cells, with optional assayed_genes, description, measurement, annotation_provenance, context and measured_evidence.")
    cells = obj.get("cells")
    if kwargs.get("lineage_map"):
        from .measurement import apply_lineage_map
        cells = apply_lineage_map(cells, kwargs["lineage_map"])
    return validate_cells(cells), obj.get("assayed_genes")


def prepare(cells, context, *, assayed_genes=None, gene_sets=None, measurement=None, lineage_map=None, annotation_provenance=None,
            measured_evidence=None):
    from .context import validate_context
    from .guidebook import build_guidebook
    if lineage_map is not None:
        from .measurement import apply_lineage_map
        cells = apply_lineage_map(cells, lineage_map)
    cells = validate_cells(cells)
    context = dict(context)
    v3 = context.get("schema_version") == "3.0"
    symbol = (lambda g: g) if v3 else str.upper
    if assayed_genes is not None:
        if "assayed_genes" in context and {symbol(g) for g in context["assayed_genes"]} != {symbol(g) for g in assayed_genes}:
            raise ValueError("Input and context disagree about the assayed gene panel.")
        context["assayed_genes"] = sorted(assayed_genes)
    context = validate_context(context)
    geometry = None
    if v3:
        from .measurement import validate_measurement
        if any(any(k in cell for k in ("polygon", "polygon_um", "segmentation_polygon", "t0_geometry", "geometry")) for cell in cells):
            from ..v3.geometry import normalize_geometry
            cells, geometry = normalize_geometry(cells)
        genes = context.get("assayed_genes", sorted({g for c in cells for g in c["expression"]}))
        raw_measurement = measurement if measurement is not None else context.get("measurement")
        measurement = validate_measurement(raw_measurement, cells, genes, context)
        if measurement["unit"] == "binned_spot":
            from copy import deepcopy
            cells = deepcopy(cells)
            for cell in cells:
                if (measurement["deconvolution_method"] == "not_declared" or
                        cell.get("purity", 0) < measurement["purity_threshold"]):
                    cell.setdefault("annotation_before_purity_gate", {k: cell[k] for k in ("cell_type", "lineage")})
                    cell["cell_type"] = cell["lineage"] = "unknown"
                    cell["annotation_gate"] = "insufficient_or_undeclared_bin_purity"
    if measured_evidence is not None or any("measured_scores" in cell for cell in cells):
        if not v3:
            raise ValueError("Measured scores require schema 3.0.")
        from ..v5.evidence import validate_measured_evidence
        measured_evidence = validate_measured_evidence(measured_evidence, cells)
    guidebook = build_guidebook(cells, context, gene_sets=gene_sets, measured_evidence=measured_evidence)
    content = dict(schema_version=context["schema_version"], cells=cells, context=guidebook["context"], guidebook=guidebook)
    if measured_evidence is not None:
        content["measured_evidence"] = measured_evidence
        guidebook["provenance"]["measured_evidence_sha256"] = digest(measured_evidence)
    if v3:
        from .measurement import validate_measurement
        content["measurement"] = measurement
        guidebook["provenance"]["measurement_sha256"] = digest(content["measurement"])
        if geometry is not None:
            content["geometry"] = geometry
            guidebook["provenance"]["geometry_sha256"] = digest(geometry)
        annotation_provenance = dict(annotation_provenance or {})
        annotation_provenance.setdefault("source", "user_supplied_object")
        annotation_provenance.setdefault("verification", "not_independently_verified")
        if lineage_map is not None:
            annotation_provenance = dict(annotation_provenance or {}, lineage_map=lineage_map, lineage_map_sha256=digest(lineage_map))
        if annotation_provenance is not None:
            content["annotation_provenance"] = annotation_provenance
            guidebook["provenance"]["annotation_sha256"] = digest(annotation_provenance)
    return dict(content, fingerprint=digest(content))


def validate_prepared(bundle):
    if not isinstance(bundle, dict) or bundle.get("schema_version") not in ("2.0", "3.0"):
        raise ValueError("Expected a v2 or v3 prepared JSON bundle.")
    content = {k: v for k, v in bundle.items() if k != "fingerprint"}
    if digest(content) != bundle.get("fingerprint"):
        raise ValueError("Prepared bundle changed; rebuild it after reviewing context/guidebook edits.")
    cells = validate_cells(bundle["cells"])
    from .context import validate_context
    if validate_context(bundle["context"]) != bundle["context"]:
        raise ValueError("Context is not normalized; rebuild the prepared bundle.")
    if bundle["context"] != bundle["guidebook"]["context"]:
        raise ValueError("Bundle and guidebook contexts differ.")
    if bundle["schema_version"] != bundle["context"]["schema_version"]:
        raise ValueError("Bundle and context schema versions differ.")
    v3 = bundle["schema_version"] == "3.0"
    from .guidebook import BASE_WEIGHTS
    rules = None
    if v3:
        from ..v3.ontology import lineage_rules, load_ontology
        from .measurement import validate_measurement
        rules = lineage_rules(bundle["context"]["ontology_version"])
        expected_hashes = load_ontology(bundle["context"]["ontology_version"])["hashes"]
        if bundle["guidebook"]["provenance"].get("ontology_hashes") != expected_hashes:
            raise ValueError("Prepared ontology hashes differ from the installed runtime; restore the pinned data or rebuild the bundle.")
        measurement = bundle.get("measurement")
        if not isinstance(measurement, dict) or validate_measurement(measurement, cells, bundle["guidebook"]["provenance"]["assayed_genes"], bundle["context"]) != measurement:
            raise ValueError("A normalized measurement record is required for v3.")
        if digest(measurement) != bundle["guidebook"]["provenance"].get("measurement_sha256"):
            raise ValueError("Measurement provenance does not match the guidebook.")
        if "measured_evidence" in bundle or any("measured_scores" in cell for cell in cells):
            from ..v5.evidence import validate_measured_evidence
            manifest = bundle.get("measured_evidence")
            if (manifest is None or validate_measured_evidence(manifest, cells) != manifest
                    or digest(manifest) != bundle["guidebook"]["provenance"].get("measured_evidence_sha256")):
                raise ValueError("Measured-evidence manifest is missing, not normalized or differs from the guidebook provenance.")
        policy = bundle["context"].get("evidence_policy")
        if (policy is None) != ("evidence_policy_sha256" not in bundle["guidebook"]["provenance"]) or (
                policy is not None and digest(policy) != bundle["guidebook"]["provenance"]["evidence_policy_sha256"]):
            raise ValueError("Evidence policy differs from the one recorded in the guidebook; rebuild the bundle.")
        if "geometry" in bundle:
            from ..v3.geometry import normalize_geometry
            normalized, geometry = normalize_geometry(cells)
            if normalized != cells or geometry != bundle["geometry"] or digest(geometry) != bundle["guidebook"]["provenance"].get("geometry_sha256"):
                raise ValueError("Polygon geometry or its provenance differs from the normalized input.")
        if measurement["unit"] == "binned_spot":
            for cell in cells:
                if (measurement["deconvolution_method"] == "not_declared" or cell.get("purity", 0) < measurement["purity_threshold"]):
                    if cell["lineage"] != "unknown":
                        raise ValueError("A mixed bin without declared adequate purity must remain unknown.")
    entries = bundle["guidebook"]["cells"]
    if set(entries) != {c["cell_id"] for c in cells}:
        raise ValueError("Guidebook cells do not match the input cell IDs.")
    for c in cells:
        entry = entries[c["cell_id"]]
        if entry.get("lineage") != c["lineage"] or entry.get("cell_type") != c["cell_type"] or c["lineage"] not in (rules if v3 else BASE_WEIGHTS):
            raise ValueError("Guidebook identity/lineage differs from input annotations.")
        allowed = set(rules[c["lineage"]]["allowed_actions"] if v3 else BASE_WEIGHTS[c["lineage"]])
        if v3 and not bundle["context"]["recognition"]:
            allowed.discard("KILL")
        if not v3 and c["lineage"] == "T_CD8" and bundle["context"]["recognition"]:
            allowed.add("KILL")
        actions = entry.get("actions", {})
        if not isinstance(actions, dict) or "REST" not in actions or set(actions) - allowed:
            raise ValueError("Guidebook includes a lineage-incompatible action or omits REST.")
        for spec in actions.values():
            weight = spec.get("weight")
            if type(weight) not in (int, float) or not math.isfinite(weight) or weight <= 0 or not isinstance(spec.get("description"), str):
                raise ValueError("Guidebook action requires a positive finite weight and description.")
    return bundle
