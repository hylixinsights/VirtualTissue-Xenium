#!/usr/bin/env python3
"""Build a raw ROI object for one LUAD3 TLS candidate region. Offline; no API calls.

  python examples/v5/luad_tls/build_inputs.py --data ../data/Xenium --candidate 5 \
      --output ../runs/luad_B5/object.json

Reads ``LUAD3_TLS_GO20_UCell.h5ad`` and ``cell_boundaries.csv`` from --data and
writes one JSON object holding, for every cell of the chosen geodesic candidate:

- expression from the declared layer (default ``purified_raw``, SPLIT-corrected);
- the existing ``first_type`` label, unchanged (lineages come from a lineage map
  at preparation time, never from this script);
- ``obsm['spatial']`` centroids in micrometres and the Xenium boundary polygon,
  joined by exact cell ID;
- the 20 saved GO Biological Process UCell scores as ``measured_scores``; a
  score that was not calculated stays null and is never written as zero;
- a ``measured_evidence`` manifest built from the object's own ``uns`` tables.

It verifies the H5AD checksum against ``h5ad_validation.json`` and the region
size against ``TLS_number_of_cells.csv``. It selects no cells by expression,
fits nothing and changes no label. The data are not packaged with the software.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

H5AD = "LUAD3_TLS_GO20_UCell.h5ad"
CAVEATS = {
    "GOBP_MITOTIC_CELL_CYCLE": "Proliferation proxy from a broad GO set; scored genes include signalling regulators as well as cycling markers.",
    "GOBP_RESPONSE_TO_OXYGEN_LEVELS": "Response to oxygen levels is broader than hypoxia.",
    "GOBP_TYPE_II_INTERFERON_PRODUCTION": "Interferon-gamma production, not an interferon-response signature.",
    "GOBP_APOPTOTIC_SIGNALING_PATHWAY": "Contains pro- and anti-apoptotic regulators; the score has no single direction.",
}
GENERAL_CAVEAT = ("Expression-signature score on a 302-gene targeted panel; not an enrichment statistic or a direct measure of "
                  "pathway activity. GO sets can include regulators of both signs. Compare cells within a process, not across processes.")


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def region_counts(path, candidate):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row["Candidato"] == f"B{candidate}":
                return dict(core=int(row["Core B"]), surrounding=int(row["Entorno"]), total=int(row["Total"]))
    raise ValueError(f"Candidate B{candidate} is not listed in {Path(path).name}.")


def manifest(adata, source):
    """Describe each saved score from the object's own gene-set and scoring tables."""
    table, scoring = adata.uns["go20_gene_sets"], adata.uns["go20_scoring"]
    provenance = json.loads(text(scoring["source_provenance_json"]))
    parameters = dict(max_rank=int(scoring["max_rank"]), ties_method=text(scoring["ties_method"]),
                      missing_genes=text(scoring["missing_genes"]))
    scores = {}
    for _, row in table.iterrows():
        genes = [g for g in text(row["genes_available"]).split(";") if g]
        if len(genes) != int(row["n_genes_available"]):
            raise ValueError(f"Gene list and count disagree for {row['geneset']}.")
        scores[text(row["score_column"])] = dict(
            method=f"{text(scoring['method'])} {text(scoring['version'])}",
            description=f"{text(row['process'])} ({text(row['go_id'])}); UCell rank score on a 0-1 scale.",
            parameters=parameters, expression_layer=text(row["expression_layer"]),
            gene_set=text(row["geneset"]), gene_set_id=text(row["go_id"]), gene_set_version=text(row["geneset_version"]),
            source_url=text(row["source_url"]), genes_scored=genes, n_genes_set=int(row["n_genes_geneset"]),
            scored_population=text(scoring["population"]),
            missing_value="null means not calculated; zero is a calculated zero score",
            **({"caveat": CAVEATS[text(row["geneset"])]} if text(row["geneset"]) in CAVEATS else {}))
    return dict(manifest_version="v5-measured-evidence-1", scores=scores, source=source, caveat=GENERAL_CAVEAT,
                gene_set_collection=dict(url=provenance["source_url"], sha256=provenance["source_sha256"]),
                relation_to_expression="Computed from the same t0 expression layer as the cells' counts; correlated with any "
                                       "gene-set evidence the simulator computes from those counts.")


def polygons(path, wanted):
    """{cell_id: closed ring} for the wanted IDs, in file order; long CSV, exact integer IDs."""
    rings = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if header != ["cell_id", "vertex_x", "vertex_y"]:
            raise ValueError("cell_boundaries.csv must have cell_id, vertex_x, vertex_y columns.")
        for cell_id, x, y in reader:
            if cell_id in wanted:
                rings.setdefault(cell_id, []).append([float(x), float(y)])
    return rings


def build(data, candidate, *, layer="purified_raw", label_key="first_type", core_only=False, verify_hash=True,
          with_scores=True):
    import anndata
    import numpy as np

    data = Path(data)
    path = data / H5AD
    validation = json.loads((data / "h5ad_validation.json").read_text(encoding="utf8"))
    checksum = file_sha256(path) if verify_hash else None
    if verify_hash and checksum != validation["sha256"]:
        raise ValueError("H5AD checksum differs from h5ad_validation.json; restore the validated file.")
    expected = region_counts(data / "TLS_number_of_cells.csv", candidate)
    adata = anndata.read_h5ad(path, backed="r")
    try:
        obs = adata.obs
        if not (obs.index.astype(str) == obs["cell_id"].astype(str)).all() or not obs.index.is_unique:
            raise ValueError("obs_names must equal the unique cell_id column.")
        if layer not in adata.layers:
            raise ValueError(f"Layer {layer!r} is absent from the object.")
        in_region = (obs["geodesic_candidate_id"] == candidate).to_numpy()
        core = (obs["is_candidate_core"] == 1).to_numpy()
        if int(in_region.sum()) != expected["total"] or int((in_region & core).sum()) != expected["core"]:
            raise ValueError("Region size disagrees with TLS_number_of_cells.csv.")
        rows = np.flatnonzero(in_region & core) if core_only else np.flatnonzero(in_region)
        columns = [text(c) for c in adata.uns["go20_scoring"]["score_columns"]]
        source = dict(file=H5AD, sha256=checksum or validation["sha256"], sha256_verified_now=bool(verify_hash),
                      n_cells=int(adata.n_obs), n_genes=int(adata.n_vars),
                      original_object=text(adata.uns["jev_export"]["source_filename"]))
        evidence = manifest(adata, source)
        if sorted(evidence["scores"]) != sorted(columns):
            raise ValueError("Score columns and gene-set table disagree.")
        genes = [str(g) for g in adata.var_names]
        ids = [str(i) for i in obs.index[rows]]
        matrix = adata.layers[layer][rows]
        matrix = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
        if not np.isfinite(matrix).all() or (matrix < 0).any():
            raise ValueError("The expression layer must hold finite nonnegative values.")
        xy = np.asarray(adata.obsm["spatial"])[rows]
        units = text(adata.uns["spatial_coordinate_metadata"]["units"])
        if units != "micrometers":
            raise ValueError("obsm['spatial'] must be declared in micrometres.")
        if obs[label_key].iloc[rows].isna().any():
            raise ValueError(f"Column {label_key!r} has missing labels in this region.")
        labels = obs[label_key].astype(str).to_numpy()[rows]
        scored = obs["go20_scored"].to_numpy()[rows]
        values = obs[columns].to_numpy(dtype=float)[rows]
        zone = core[rows]
    finally:
        adata.file.close()
    rings = polygons(data / "cell_boundaries.csv", set(ids))
    cells, zones = [], {}
    for i, cell_id in enumerate(ids):
        scores = {}
        for name, value in zip(columns, values[i]):
            if math.isinf(value):
                raise ValueError("Saved scores must be finite or NaN.")
            scores[name] = None if math.isnan(value) else float(value)
        if bool(scored[i]) != any(v is not None for v in scores.values()):
            raise ValueError("go20_scored disagrees with the saved score columns.")
        cell = dict(cell_id=cell_id, cell_type=str(labels[i]), lineage="unmapped", sample_id="LUAD3",
                    x=float(xy[i, 0]), y=float(xy[i, 1]),
                    expression={g: float(v) for g, v in zip(genes, matrix[i]) if v != 0})
        if with_scores:
            cell["measured_scores"] = scores
        if cell_id in rings:
            cell["polygon"] = rings[cell_id]
        cells.append(cell)
        zones[cell_id] = "core" if zone[i] else "surrounding"
    region = dict(candidate=f"B{candidate}", definition="geodesic TLS candidate, core plus surroundings",
                  cells=len(cells), expected=expected, core_only=core_only,
                  core_cells=sum(v == "core" for v in zones.values()),
                  cells_with_polygon=len(rings), polygon_join="exact cell_id",
                  scored_cells=int(scored.sum()), measured_scores_included=with_scores)
    result = dict(
        cells=cells, assayed_genes=genes,
        description=f"LUAD3 geodesic TLS candidate B{candidate}; existing labels, SPLIT-corrected expression, saved GO20 UCell scores.",
        measurement=dict(platform="xenium", unit="segmented_cell", panel_type="targeted",
                         segmentation_method="Xenium cell_boundaries.csv as supplied; not independently verified",
                         expression_layer=f"{layer} (SPLIT-corrected expression; non-integer values)",
                         coordinate_source="obsm['spatial'] = [x_centroid, y_centroid], micrometres"),
        annotation_provenance=dict(source="existing obs column", annotation_column=label_key,
                                   verification="not_independently_verified", source_object=source, region=region,
                                   note="A second annotation column (cell_type_fine) disagrees with first_type for many cells; "
                                        "labels are used as declared and mapped to lineages by an explicit CSV."))
    if with_scores:
        result["measured_evidence"] = evidence
    return result, zones, region


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True, type=Path, help="Folder holding the LUAD3 H5AD, boundaries and validation files.")
    p.add_argument("--candidate", required=True, type=int, help="Geodesic candidate number, 1-28 (for example 5 for B5).")
    p.add_argument("--output", required=True, type=Path, help="New object JSON path; existing files are refused.")
    p.add_argument("--layer", default="purified_raw")
    p.add_argument("--label-key", default="first_type")
    p.add_argument("--core-only", action="store_true", help="Keep only core cells of the candidate.")
    p.add_argument("--no-scores", action="store_true",
                   help="Leave the saved GO scores out: the object then prepares without an evidence policy, as in v4.")
    p.add_argument("--skip-hash", action="store_true", help="Do not recompute the H5AD checksum (it is still recorded as unverified).")
    args = p.parse_args(argv)
    zones_path = args.output.with_suffix(".zones.json")
    if args.output.exists() or zones_path.exists():
        p.error("Output exists; choose a new path.")
    try:
        obj, zones, region = build(args.data, args.candidate, layer=args.layer, label_key=args.label_key,
                                   core_only=args.core_only, verify_hash=not args.skip_hash,
                                   with_scores=not args.no_scores)
    except ImportError as exc:
        p.error(f"This script needs the xenium optional dependencies (anndata, numpy, pandas): {exc}")
    except (OSError, ValueError, KeyError) as exc:
        p.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for path, value in ((args.output, obj), (zones_path, zones)):
        with path.open("x", encoding="utf8") as handle:
            json.dump(value, handle, allow_nan=False)
            handle.write("\n")
    print(json.dumps(dict(output=str(args.output.resolve()), zones=str(zones_path.resolve()), region=region,
                          measured_scores=len(obj.get("measured_evidence", {}).get("scores", {})),
                          panel_genes=len(obj["assayed_genes"])), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
