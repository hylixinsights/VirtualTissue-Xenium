#!/usr/bin/env python3
"""Export existing annotated H5AD or raw JSON to a spatial ROI input, offline.

Example:
  python tools/export_object_v3.py object.h5ad object.json \
      --cell-type-key existing_labels --lineage-map lineages.csv --layer counts

This exports every input cell, including low-count cells, without fitting scores,
inferring annotations, selecting an ROI or calling Jev. H5AD var_names preserve
the entire assayed panel, including genes with zero observed counts. Raw JSON
must declare its complete assayed_genes list; a sparse feature union is not a
complete panel. Use a count matrix: finite nonnegative values are checked, but
normalization history cannot be established from positive numbers alone.
Supplied polygon QC never certifies a true biological boundary.

The default limit is 1000 cells. H5AD observation count is checked before AnnData
or expression layers are loaded. Larger objects require an explicit --max-cells;
raising it can consume substantial memory because H5AD layers may be loaded by
AnnData and the output stores all sparse per-cell dictionaries. No truncation or
automatic subsampling occurs. JSON itself must be parsed before its count can be
checked. Select a subset in the source application if the full object is too big.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from virtualtissue_xenium.v2.inputs import digest, load_cells, read_json, validate_cells
from virtualtissue_xenium.v2.measurement import read_lineage_map
from virtualtissue_xenium.v3.geometry import normalize_geometry
from virtualtissue_xenium.v3.object_input import attach_polygons

EXPORT_VERSION = "raw-object-export-1"


def file_sha256(path):
    """Hash input bytes with bounded memory, independently of JSON normalization."""
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def _file_record(path):
    path = Path(path)
    return {"name": path.name, "sha256": file_sha256(path)}


def _positive_scale(value, name):
    if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive.")


def _write_new_json(path, value):
    """Exclusive creation refuses existing files, directories and symlinks."""
    data = json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    # Serialize before opening so invalid metadata never leaves a partial export.
    with path.open("x", encoding="utf-8") as handle:
        try:
            handle.write(data)
        except BaseException:
            path.unlink(missing_ok=True)
            raise


def export_object(input_path, output_path, *, max_cells=1000, cell_type_key="cell_type",
                  lineage_key="lineage", sample_key=None, layer=None, spatial_key="spatial",
                  coordinate_scale=1.0, polygon_key=None, lineage_map=None, polygons=None,
                  polygon_scale=1.0, measurement=None):
    """Return and save a raw ROI input; optional map/measurement arguments are paths.

    Measurement declarations are retained rather than biologically validated.
    Their full schema is checked when the user-selected ROI is prepared with its
    explicit context. Polygon-derived centroids are allowed and recorded; absent
    positions with no valid polygon are refused for this spatial export.
    """
    source, output = Path(input_path), Path(output_path)
    if output.exists() or output.is_symlink():
        raise ValueError("Output already exists; choose a new JSON path.")
    if not source.is_file():
        raise ValueError("Input must be an existing H5AD or raw JSON file.")
    if type(max_cells) is not int or max_cells < 1:
        raise ValueError("max_cells must be a positive integer; no automatic subsampling is performed.")
    _positive_scale(coordinate_scale, "coordinate_scale")
    _positive_scale(polygon_scale, "polygon_scale")
    if polygon_key and polygons:
        raise ValueError("Choose polygon_key or an external polygons file, not both.")
    if polygon_scale != 1 and not polygons:
        raise ValueError("polygon_scale applies only to an external polygons file.")
    kind = source.suffix.lower()
    if kind not in {".h5ad", ".json"}:
        raise ValueError("Use an .h5ad or raw .json input file.")
    if kind == ".h5ad" and source.suffix != ".h5ad":
        raise ValueError("The H5AD filename must use the .h5ad extension.")
    if kind == ".json" and (polygon_key or layer or sample_key or coordinate_scale != 1
                            or cell_type_key != "cell_type" or lineage_key != "lineage"
                            or spatial_key != "spatial"):
        raise ValueError("Annotation/coordinate column options and polygon_key apply only to H5AD; raw JSON uses its declared fields in micrometres.")
    mapping = read_lineage_map(lineage_map) if lineage_map else None
    raw = read_json(source) if kind == ".json" else {}
    if not isinstance(raw, dict):
        raise ValueError("Raw JSON must be an object with cells and assayed_genes.")
    if kind == ".json" and isinstance(raw.get("cells"), list) and len(raw["cells"]) > max_cells:
        raise ValueError(f"Select a region/subset first: {len(raw['cells'])} cells exceeds max_cells={max_cells}.")
    cells, panel = load_cells(source, max_cells=max_cells, cell_type_key=cell_type_key,
        lineage_key=lineage_key, sample_key=sample_key, layer=layer, spatial_key=spatial_key,
        coordinate_scale=coordinate_scale, lineage_map=mapping, polygon_key=polygon_key)
    if len(cells) > max_cells:
        raise ValueError(f"Input cells exceed max_cells={max_cells}; no cells were exported.")
    if (not isinstance(panel, list) or not panel or any(not isinstance(g, str) or not g.strip() or g != g.strip() for g in panel)
            or len(set(panel)) != len(panel)):
        raise ValueError("Declare a complete, nonempty, unique assayed_genes panel; it cannot be reconstructed from nonzero counts.")
    observed = {gene for cell in cells for gene in cell["expression"]}
    if observed - set(panel):
        raise ValueError("The declared assayed_genes panel omits observed expression features.")
    polygon_join = None
    if polygons:
        cells, polygon_join = attach_polygons(cells, read_json(polygons), coordinate_scale=polygon_scale)
    cells, geometry_qc = normalize_geometry(cells, source=f"raw_object_export:{source.name}")
    cells = validate_cells(cells)
    if any("x" not in cell or "y" not in cell for cell in cells):
        raise ValueError("Spatial ROI export requires finite x/y for every cell or a valid supplied polygon with a derived centroid.")

    annotation = deepcopy(raw.get("annotation_provenance") or {})
    if not isinstance(annotation, dict):
        raise ValueError("annotation_provenance must be a JSON object.")
    annotation.setdefault("source", "user_supplied_existing_labels")
    annotation.setdefault("verification", "not_independently_verified")
    audit = dict(version=EXPORT_VERSION, input_file=_file_record(source), input_format=kind[1:],
        cells_exported=len(cells), panel_features=len(panel), panel_source="h5ad_var_names" if kind == ".h5ad" else "declared_raw_assayed_genes",
        max_cells=max_cells, subset_policy="all_input_cells_no_selection_or_truncation",
        scoring_status="not_computed", annotation_status="existing_declared_labels_not_verified_by_exporter",
        expression_semantics="user_declared_counts; finite_nonnegative_values_checked; normalization_history_not_verified",
        coordinate_units="micrometres", geometry_qc=geometry_qc,
        parameters=dict(cell_type_key=cell_type_key, lineage_key=lineage_key,
            sample_key=sample_key, sample_fallback="sample_id column when present, otherwise literal sample",
            layer=layer, spatial_key=spatial_key, coordinate_scale_to_um=coordinate_scale,
            polygon_key=polygon_key, external_polygon_scale_to_um=polygon_scale))
    if "object_export" in annotation:
        audit["previous_export"] = annotation["object_export"]
    annotation["object_export"] = audit
    if mapping is not None:
        audit["lineage_map_file"] = _file_record(lineage_map)
        audit["lineage_map"] = deepcopy(mapping)
        audit["lineage_map_sha256"] = digest(mapping)
        annotation["lineage_map"] = deepcopy(mapping)
        annotation["lineage_map_sha256"] = digest(mapping)
    if polygons:
        audit["polygons_file"] = _file_record(polygons)
        audit["polygon_join"] = polygon_join
    result = dict(cells=cells, assayed_genes=panel,
                  description=raw.get("description", "Raw spatial object with existing labels; exported without scoring or ROI selection."),
                  annotation_provenance=annotation)
    if "context" in raw:
        result["context"] = deepcopy(raw["context"])
    if measurement is not None:
        declared_measurement = read_json(measurement)
        audit["measurement_file"] = _file_record(measurement)
        if "measurement" in raw:
            audit["source_measurement_before_override"] = deepcopy(raw["measurement"])
    else:
        declared_measurement = deepcopy(raw.get("measurement"))
    if declared_measurement is not None:
        if not isinstance(declared_measurement, dict):
            raise ValueError("measurement must be a JSON object; full validation occurs during preparation.")
        result["measurement"] = declared_measurement
        audit["measurement_sha256"] = digest(declared_measurement)
    audit["measurement_status"] = "declared_not_normalized" if declared_measurement is not None else "not_declared"
    _write_new_json(output, result)
    return result


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path)
    p.add_argument("output", type=Path, help="New raw JSON path; existing outputs are refused.")
    p.add_argument("--max-cells", type=int, default=1000, help="Explicit input limit, not a subsampling target; larger values can consume substantial memory.")
    p.add_argument("--cell-type-key", default="cell_type")
    p.add_argument("--lineage-key", default="lineage")
    p.add_argument("--sample-key")
    p.add_argument("--layer", help="H5AD count layer; omit to use X. Finite nonnegative values required; normalization history is not inferred.")
    p.add_argument("--spatial-key", default="spatial")
    p.add_argument("--coordinate-scale", type=float, default=1.0, help="H5AD coordinate and polygon-key conversion to micrometres.")
    p.add_argument("--polygon-key", help="H5AD uns/obs polygon key; no segmentation is inferred.")
    p.add_argument("--lineage-map", type=Path, help="Explicit CSV with cell_type,lineage columns.")
    p.add_argument("--polygons", type=Path, help="Exact-cell-ID polygon map or GeoJSON FeatureCollection.")
    p.add_argument("--polygon-scale", type=float, default=1.0, help="External polygon conversion to micrometres.")
    p.add_argument("--measurement", type=Path, help="JSON measurement declaration preserved for later preparation.")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if args.max_cells > 1000:
        print("Explicit large-object export: all cells up to --max-cells are retained; memory use can be substantial, including H5AD layers. No automatic subset is selected.", file=sys.stderr)
    values = vars(args).copy()
    source, output = values.pop("input"), values.pop("output")
    try:
        result = export_object(source, output, **values)
    except ImportError as exc:
        p.error(f"H5AD export requires the project's xenium optional dependencies (including anndata and h5py): {exc}")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        p.error(str(exc))
    audit = result["annotation_provenance"]["object_export"]
    print(json.dumps(dict(output=str(output.resolve()), cells=len(result["cells"]), panel_features=len(result["assayed_genes"]),
        geometry_qc=audit["geometry_qc"], scoring="not_computed", roi="not_selected"), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
