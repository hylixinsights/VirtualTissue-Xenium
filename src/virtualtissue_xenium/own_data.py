"""Use your own Xenium run in Part 1.

Part 1 (``xenium_prepare``) was written and validated for the public example, so it
expects two things that a fresh Xenium output folder does not have yet:

1. A *data folder* with ``raw/`` plus ``SOURCE.json`` and ``DOWNLOAD_MANIFEST.json``
   (file sizes and SHA-256). Part 1 refuses to run if a listed file is missing or
   changed. :func:`create_data_folder` or :func:`write_manifest` create these files.
2. A *cluster annotation map* that turns the 10x graph-based clusters of **your**
   sample into cell types. The packaged map (``engine/resources/cluster_annotations.json``)
   describes the 32 clusters of the public sample only. :func:`write_annotation_template`
   drafts a map for your clusters (with their top differentially expressed genes);
   you fill in the labels, check it with :func:`validate_annotations` and activate it
   with :func:`use_cluster_annotations`.

None of these helpers changes the engine code. ``use_cluster_annotations`` only points
the engine's resource folder (``xenium_prepare.RESOURCES``) to a copy in which the
cluster map is replaced.

pandas/pyarrow/h5py are imported inside the functions that need them.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import tarfile
from datetime import datetime, timezone
from pathlib import Path

from . import ENGINE_DIR

PUBLIC_SAMPLE_ID = "Xenium_Prime_Human_Lung_Cancer_FFPE"

REQUIRED_FILES = {
    "cell_feature_matrix.h5": "Sparse cell-feature count matrix (10x HDF5); only 'Gene Expression' features are used.",
    "cells.parquet": "Cell IDs, XY centroids (micrometres), transcript counts, cell and nucleus areas.",
    "cell_boundaries.parquet": "Simplified cell polygons (cell_id, vertex_x, vertex_y).",
    "nucleus_boundaries.parquet": "Simplified nucleus polygons; copied into the bundle for the selected cells.",
    "analysis.tar.gz": "Xenium Onboard Analysis secondary analysis; graph-based clusters and differential expression are used.",
}
OPTIONAL_FILES = {
    "experiment.xenium": "Run manifest (provenance only).",
    "gene_panel.json": "Panel definition (provenance only).",
    "metrics_summary.csv": "Run-level QC metrics (provenance only).",
    "analysis_summary.html": "Vendor web summary (provenance only).",
    "cells.zarr.zip": "Full segmentation masks (provenance only; Part 1 uses the parquet polygons).",
}

#: Values accepted by the engine (see xenium_prepare.COLORS and runtime_template.json).
COMPARTMENTS = ["tumor_candidate", "fibroblast", "immune", "stroma", "epithelial_other", "unknown"]
RUNTIME_LABELS = ["Tumor", "CAF", "Fibroblast", "Endothelial", "Macrophage", "T_CD8_cytotoxic",
                  "T_CD4", "B_cell", "Plasma_cell"]
RUNTIME_LABEL_NOTES = {
    "Tumor": "Epithelial tumour candidate. Can REMAIN_QUIESCENT, ENTER_CYCLE, UPREGULATE_PDL1, SECRETE_VEGFA, metabolic shifts.",
    "CAF": "Cancer-associated fibroblast. DEPOSIT_MATRIX, SECRETE_TGFB, SECRETE_VEGFA at a tumour interface.",
    "Fibroblast": "Fibroblast. DEPOSIT_MATRIX, SECRETE_TGFB at a tumour interface.",
    "Endothelial": "Endothelial cell. STAY_QUIESCENT or SPROUT (only if angiogenesis is enabled).",
    "Macrophage": "Generic macrophage. APPROACH_CORPSE, EFFEROCYTOSE, SECRETE_CXCL9, SECRETE_VEGFA_TGFB, PATROL.",
    "T_CD8_cytotoxic": "Cytotoxic T cell. Migration, APPROACH_TUMOR, HOLD_CONTACT, DEGRANULATE, PROLIFERATE, REST.",
    "T_CD4": "Helper T cell. HELP_CD8 at a stimulated CD8 contact, MIGRATE_UP_CXCL9.",
    "B_cell": "Physical body only: no action menu is implemented.",
    "Plasma_cell": "Physical body only: no action menu is implemented.",
    None: "null = the cells stay as physical context (geometry, contacts) and never decide.",
}
CLUSTERS_MEMBER = "analysis/clustering/gene_expression_graphclust/clusters.csv"
DIFFEXP_MEMBER = "analysis/diffexp/gene_expression_graphclust/differential_expression.csv"


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# 1. Data folder
# ---------------------------------------------------------------------------
def check_xenium_files(raw_dir) -> dict:
    """Open the required files and check the structure Part 1 relies on.

    Returns a small summary (cells, genes, clusters). Raises ``ValueError`` with a
    readable message when something is missing.
    """
    import h5py
    import pandas as pd
    import pyarrow.parquet as pq

    raw_dir = Path(raw_dir)
    missing = [name for name in REQUIRED_FILES if not (raw_dir / name).is_file()]
    if missing:
        raise ValueError(f"Missing required Xenium files in {raw_dir}: {', '.join(missing)}. "
                         "If your output has an extracted 'analysis/' folder instead of analysis.tar.gz, "
                         "create it with:  tar -czf analysis.tar.gz analysis")
    with h5py.File(raw_dir / "cell_feature_matrix.h5", "r") as f:
        if "matrix" not in f:
            raise ValueError("cell_feature_matrix.h5 has no 'matrix' group (10x HDF5 format expected).")
        types = f["matrix/features/feature_type"].asstr()[:]
        n_cells = int(f["matrix/shape"][1])
    n_genes = int((types == "Gene Expression").sum())
    if n_genes == 0:
        raise ValueError("No 'Gene Expression' features found in cell_feature_matrix.h5.")
    needed = {"cell_id", "x_centroid", "y_centroid", "transcript_counts", "cell_area", "nucleus_area"}
    columns = set(pq.read_schema(raw_dir / "cells.parquet").names)
    if needed - columns:
        raise ValueError(f"cells.parquet lacks columns: {sorted(needed - columns)} "
                         "(Xenium Onboard Analysis 2.0 or later is expected; the example used 3.0).")
    for name in ("cell_boundaries.parquet", "nucleus_boundaries.parquet"):
        cols = set(pq.read_schema(raw_dir / name).names)
        if {"cell_id", "vertex_x", "vertex_y"} - cols:
            raise ValueError(f"{name} must contain cell_id, vertex_x and vertex_y.")
    with tarfile.open(raw_dir / "analysis.tar.gz", "r:gz") as archive:
        names = set(archive.getnames())
        for member in (CLUSTERS_MEMBER, DIFFEXP_MEMBER):
            if member not in names:
                raise ValueError(f"analysis.tar.gz lacks {member}.")
        clusters = pd.read_csv(archive.extractfile(CLUSTERS_MEMBER))
    return dict(cells=n_cells, gene_expression_features=n_genes, other_features=int(len(types) - n_genes),
                graphclust_clusters=int(clusters.Cluster.nunique()))


def write_manifest(data_dir, sample_id: str, *, source_page: str | None = None, license: str | None = None,
                   attribution: str | None = None, notes: str | None = None) -> Path:
    """Write ``DOWNLOAD_MANIFEST.json`` and ``SOURCE.json`` for files already in ``<data_dir>/raw``.

    Every regular file in ``raw/`` is listed with its size and SHA-256 (Part 1
    re-checks them on every run, so do not edit files afterwards).
    """
    data_dir = Path(data_dir)
    raw = data_dir / "raw"
    summary = check_xenium_files(raw)
    roles = {**REQUIRED_FILES, **OPTIONAL_FILES}
    now = datetime.now(timezone.utc).isoformat()
    files = []
    for path in sorted(p for p in raw.iterdir() if p.is_file() and not p.name.startswith(".")):
        files.append(dict(path=f"raw/{path.name}", role=roles.get(path.name, "Additional file (not used by Part 1)."),
                          bytes=path.stat().st_size, sha256=_sha256(path), source_archive_member=path.name,
                          registered_utc=now))
    manifest = dict(sample_id=sample_id, source_page=source_page, files=files, complete=True,
                    created_by="virtualtissue_xenium.own_data.write_manifest")
    source = dict(sample_id=sample_id, source_page=source_page, license_as_stated=license, attribution=attribution,
                  notes=notes, registered_utc=now, xenium_summary=summary,
                  scope="Original Xenium output files; no normalization, annotation or scoring performed.")
    (data_dir / "DOWNLOAD_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (data_dir / "SOURCE.json").write_text(json.dumps(source, indent=2) + "\n")
    return data_dir / "DOWNLOAD_MANIFEST.json"


def create_data_folder(xenium_output_dir, data_dir, sample_id: str, *, include_optional: bool = False,
                       source_page: str | None = None, license: str | None = None,
                       attribution: str | None = None, progress=print) -> Path:
    """Copy the files Part 1 needs from a Xenium output folder and register them.

    ``xenium_output_dir`` is the instrument/XOA output folder (often called ``outs``).
    Nothing in it is modified. ``data_dir`` must be a new or empty folder.
    """
    src, data_dir = Path(xenium_output_dir), Path(data_dir)
    raw = data_dir / "raw"
    if raw.exists() and any(raw.iterdir()):
        raise FileExistsError(f"{raw} is not empty; choose a new data folder.")
    raw.mkdir(parents=True, exist_ok=True)
    wanted = list(REQUIRED_FILES) + (list(OPTIONAL_FILES) if include_optional else [])
    for name in wanted:
        if (src / name).is_file():
            progress(f"  copying {name}")
            shutil.copy2(src / name, raw / name)
        elif name in REQUIRED_FILES:
            raise FileNotFoundError(f"{name} not found in {src}.")
    write_manifest(data_dir, sample_id, source_page=source_page, license=license, attribution=attribution)
    progress(f"Data folder ready: {data_dir}")
    return data_dir


# ---------------------------------------------------------------------------
# 2. Cluster annotation map
# ---------------------------------------------------------------------------
#: Canonical markers summarised per cluster (same groups as Part 1's marker heat map).
MARKER_GROUPS = {
    "Epithelium": ["EPCAM", "CDH1", "MUC1", "MAL2", "FOLR1", "KRT8", "KRT18", "NKX2-1"],
    "Fibroblast": ["PDGFRA", "TCF21", "FAP", "COL5A1", "DPT", "COL1A1"],
    "T / NK": ["CD3E", "CD3G", "CD4", "CD8A", "KLRD1", "PRF1", "FOXP3"],
    "Myeloid": ["CSF1R", "CD68", "MRC1", "MARCO", "CD163"],
    "B / plasma": ["MS4A1", "CD79A", "MZB1", "JCHAIN"],
    "Endothelial": ["CLDN5", "CDH5", "ROBO4", "PECAM1", "VWF"],
}


def cluster_marker_detection(data_dir, markers=None):
    """Fraction of cells of each 10x cluster with ≥ 1 count of each marker gene.

    Returns a DataFrame (rows: clusters, columns: markers present in the panel).
    Use it to decide each cluster's label and runtime label.
    """
    import h5py
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    raw = Path(data_dir) / "raw"
    groups = MARKER_GROUPS if markers is None else markers
    wanted = [g for genes in groups.values() for g in genes] if isinstance(groups, dict) else list(groups)
    with h5py.File(raw / "cell_feature_matrix.h5", "r") as f:
        m = f["matrix"]
        barcodes = m["barcodes"].asstr()[:]
        names = m["features/name"].asstr()[:]
        types = m["features/feature_type"].asstr()[:]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:]))
    index = {n: i for i, (n, t) in enumerate(zip(names, types)) if t == "Gene Expression"}
    present = [g for g in dict.fromkeys(wanted) if g in index]
    detected = (X.tocsr()[[index[g] for g in present]] > 0).T.tocsr()  # cells × markers
    with tarfile.open(raw / "analysis.tar.gz", "r:gz") as archive:
        clusters = pd.read_csv(archive.extractfile(CLUSTERS_MEMBER)).set_index("Barcode").Cluster
    labels = clusters.reindex(barcodes).fillna(-1).astype(int).to_numpy()
    rows = {}
    for c in sorted(set(labels) - {-1}):
        mask = labels == c
        rows[int(c)] = np.asarray(detected[mask].mean(axis=0)).ravel()
    return pd.DataFrame.from_dict(rows, orient="index", columns=present).rename_axis("cluster")


def write_annotation_template(data_dir, out_json, *, top_n: int = 15, marker_summary: bool = True) -> Path:
    """Draft a cluster annotation map for your sample.

    Each record lists the cluster, its size, its top differentially expressed genes
    (same filter as Part 1: mean counts > 0.05, adjusted p < 0.01, sorted by log2
    fold change) and, with ``marker_summary``, the fraction of its cells in which each
    canonical marker is detected. ``label`` is a placeholder: set ``label``,
    ``compartment``, ``runtime_label`` and ``markers`` for every cluster, then run
    :func:`validate_annotations`. A CSV with the marker table is written next to the JSON.
    """
    import pandas as pd

    raw = Path(data_dir) / "raw"
    with tarfile.open(raw / "analysis.tar.gz", "r:gz") as archive:
        clusters = pd.read_csv(archive.extractfile(CLUSTERS_MEMBER))
        de = pd.read_csv(archive.extractfile(DIFFEXP_MEMBER))
    counts = clusters.Cluster.value_counts()
    detection = cluster_marker_detection(data_dir) if marker_summary else None
    records = []
    for c in sorted(counts.index):
        prefix = f"Cluster {c} "
        top = []
        if prefix + "Mean Counts" in de:
            eligible = de[(de[prefix + "Mean Counts"] > .05) & (de[prefix + "Adjusted p value"] < .01)]
            top = eligible.sort_values(prefix + "Log2 fold change", ascending=False)["Feature Name"].head(top_n).tolist()
        record = dict(cluster=int(c), n_cells=int(counts[c]), top_de_genes=top)
        if detection is not None and int(c) in detection.index:
            row = detection.loc[int(c)]
            record["marker_detection_fraction"] = {g: round(float(v), 3) for g, v in row.items() if v >= .05}
        record.update(label=f"Cluster_{int(c)}_TO_REVIEW", compartment="unknown", runtime_label=None,
                      markers=[], status="needs_review",
                      source="Template from 10x graph-based clustering; fill in after marker review.")
        records.append(record)
    if detection is not None:
        detection.round(3).to_csv(Path(out_json).with_suffix(".markers.csv"))
    out = Path(out_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(records, indent=1) + "\n")
    return out


def validate_annotations(annotations, data_dir=None) -> list[str]:
    """Return a list of problems (empty list = valid).

    ``annotations`` is a path or the loaded list. With ``data_dir``, every 10x cluster
    of the sample must have exactly one record.
    """
    records = json.loads(Path(annotations).read_text()) if not isinstance(annotations, list) else annotations
    problems = []
    seen = set()
    for i, r in enumerate(records):
        where = f"record {i} (cluster {r.get('cluster')})"
        if not isinstance(r.get("cluster"), int):
            problems.append(f"{where}: 'cluster' must be an integer.")
            continue
        if r["cluster"] in seen:
            problems.append(f"{where}: duplicated cluster.")
        seen.add(r["cluster"])
        if not isinstance(r.get("label"), str) or not r["label"]:
            problems.append(f"{where}: 'label' must be a non-empty string.")
        elif r["label"].endswith("_TO_REVIEW"):
            problems.append(f"{where}: label is still the template placeholder.")
        if r.get("compartment") not in COMPARTMENTS:
            problems.append(f"{where}: compartment must be one of {COMPARTMENTS}.")
        if r.get("runtime_label") not in RUNTIME_LABELS + [None]:
            problems.append(f"{where}: runtime_label must be one of {RUNTIME_LABELS} or null.")
        if not isinstance(r.get("markers"), list):
            problems.append(f"{where}: 'markers' must be a list of gene symbols (it may be empty).")
    if data_dir is not None:
        import pandas as pd
        with tarfile.open(Path(data_dir) / "raw" / "analysis.tar.gz", "r:gz") as archive:
            sample_clusters = set(pd.read_csv(archive.extractfile(CLUSTERS_MEMBER)).Cluster.astype(int))
        if sample_clusters - seen:
            problems.append(f"clusters without a record: {sorted(sample_clusters - seen)}")
        if seen - sample_clusters:
            problems.append(f"records for clusters absent from the sample: {sorted(seen - sample_clusters)}")
    return problems


def use_cluster_annotations(prep, annotations_json, workdir) -> Path:
    """Make Part 1 use your cluster map instead of the packaged one.

    ``prep`` is the imported ``xenium_prepare`` module (``vtx.prepare``). A copy of
    the engine resources is written to ``<workdir>/resources_custom`` with
    ``cluster_annotations.json`` replaced, and ``prep.RESOURCES`` is pointed to it
    for the current Python session. Call this *before* ``prep.annotate``.
    """
    problems = validate_annotations(annotations_json)
    if problems:
        raise ValueError("Invalid cluster annotation map:\n  " + "\n  ".join(problems))
    target = Path(workdir) / "resources_custom"
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(ENGINE_DIR / "resources", target)
    shutil.copyfile(annotations_json, target / "cluster_annotations.json")
    prep.RESOURCES = target
    return target


def use_packaged_annotations(prep) -> Path:
    """Undo :func:`use_cluster_annotations` (back to the public-sample map)."""
    prep.RESOURCES = ENGINE_DIR / "resources"
    return prep.RESOURCES


def require_annotation_map(data_dir, annotations_json: str | None) -> None:
    """Stop early if a non-example sample would be annotated with the example's map.

    The packaged map only fits the public sample. On another sample with 32 or fewer
    clusters it would silently assign wrong cell types, so a custom map is required.
    """
    source = json.loads((Path(data_dir) / "DOWNLOAD_MANIFEST.json").read_text())
    if source.get("sample_id") != PUBLIC_SAMPLE_ID and not annotations_json:
        raise ValueError(
            f"This data folder is sample {source.get('sample_id')!r}, not the public example. "
            "Provide your own cluster annotation map (CLUSTER_ANNOTATIONS_JSON). "
            "Create a draft with virtualtissue_xenium.own_data.write_annotation_template(DATA_DIR, 'my_clusters.json'). "
            "See docs/06_your_own_xenium_data.md.")
