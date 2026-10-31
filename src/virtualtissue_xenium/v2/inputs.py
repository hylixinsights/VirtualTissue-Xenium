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
        if not isinstance(c["sample_id"], str) or not c["sample_id"]:
            raise ValueError("sample_id must be a nonempty string.")
        initial = c.get("initial_state", {})
        if not isinstance(initial, dict) or set(initial) - {"activated", "exhausted"} or any(type(x) is not bool for x in initial.values()):
            raise ValueError("initial_state accepts only boolean activated/exhausted assumptions.")
        result.append(c)
    return sorted(result, key=lambda c: c["cell_id"])


def from_anndata(adata, *, cell_type_key="cell_type", lineage_key="lineage", sample_key=None,
                 layer=None, spatial_key="spatial", coordinate_scale=1.0, max_cells=1000):
    """Accept an in-memory annotated object. Never infer lineage, malignancy, or XY.

    X (or the chosen layer) must hold nonnegative counts, with unique gene symbols
    in var_names. No densification of the full matrix is performed.
    """
    if adata.n_obs > max_cells:
        raise ValueError(f"Select a region/subset first: {adata.n_obs} cells exceeds max_cells={max_cells}.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("AnnData cell IDs and gene symbols must be unique.")
    if sample_key is None and "sample_id" in adata.obs:
        sample_key = "sample_id"
    for key in (cell_type_key, lineage_key, *([sample_key] if sample_key else [])):
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
                 lineage=str(adata.obs.iloc[i][lineage_key]),
                 expression={g: float(v) for g, v in zip(genes, values) if v != 0},
                 sample_id=str(adata.obs.iloc[i][sample_key]) if sample_key else "sample")
        if coords is not None:
            c.update(x=float(coords[i, 0]) * coordinate_scale, y=float(coords[i, 1]) * coordinate_scale)
        cells.append(c)
    return validate_cells(cells), genes


def load_cells(path, **kwargs):
    path = Path(path)
    if path.suffix == ".h5ad":
        import anndata
        obj = anndata.read_h5ad(path, backed="r")
        try:
            return from_anndata(obj, **kwargs)
        finally:
            obj.file.close()
    obj = read_json(path)
    if not isinstance(obj, dict) or set(obj) - {"cells", "assayed_genes", "description"}:
        raise ValueError("JSON input requires cells and optional assayed_genes/description.")
    return validate_cells(obj.get("cells")), obj.get("assayed_genes")


def prepare(cells, context, *, assayed_genes=None, gene_sets=None):
    from .context import validate_context
    from .guidebook import build_guidebook
    cells = validate_cells(cells)
    context = dict(context)
    if assayed_genes is not None:
        if "assayed_genes" in context and {g.upper() for g in context["assayed_genes"]} != {g.upper() for g in assayed_genes}:
            raise ValueError("Input and context disagree about the assayed gene panel.")
        context["assayed_genes"] = sorted(assayed_genes)
    context = validate_context(context)
    guidebook = build_guidebook(cells, context, gene_sets=gene_sets)
    content = dict(schema_version="2.0", cells=cells, context=guidebook["context"], guidebook=guidebook)
    return dict(content, fingerprint=digest(content))


def validate_prepared(bundle):
    if not isinstance(bundle, dict) or bundle.get("schema_version") != "2.0":
        raise ValueError("Expected a v2 prepared JSON bundle.")
    content = {k: v for k, v in bundle.items() if k != "fingerprint"}
    if digest(content) != bundle.get("fingerprint"):
        raise ValueError("Prepared bundle changed; rebuild it after reviewing context/guidebook edits.")
    cells = validate_cells(bundle["cells"])
    from .context import validate_context
    if validate_context(bundle["context"]) != bundle["context"]:
        raise ValueError("Context is not normalized; rebuild the prepared bundle.")
    if bundle["context"] != bundle["guidebook"]["context"]:
        raise ValueError("Bundle and guidebook contexts differ.")
    from .guidebook import BASE_WEIGHTS
    entries = bundle["guidebook"]["cells"]
    if set(entries) != {c["cell_id"] for c in cells}:
        raise ValueError("Guidebook cells do not match the input cell IDs.")
    for c in cells:
        entry = entries[c["cell_id"]]
        if entry.get("lineage") != c["lineage"] or entry.get("cell_type") != c["cell_type"] or c["lineage"] not in BASE_WEIGHTS:
            raise ValueError("Guidebook identity/lineage differs from input annotations.")
        allowed = set(BASE_WEIGHTS[c["lineage"]])
        if c["lineage"] == "T_CD8" and bundle["context"]["recognition"]:
            allowed.add("KILL")
        actions = entry.get("actions", {})
        if not isinstance(actions, dict) or "REST" not in actions or set(actions) - allowed:
            raise ValueError("Guidebook includes a lineage-incompatible action or omits REST.")
        for spec in actions.values():
            weight = spec.get("weight")
            if type(weight) not in (int, float) or not math.isfinite(weight) or weight <= 0 or not isinstance(spec.get("description"), str):
                raise ValueError("Guidebook action requires a positive finite weight and description.")
    return bundle
