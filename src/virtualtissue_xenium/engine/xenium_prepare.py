"""One-sample Xenium preparation: immutable inputs, AUCell, visual ROI selection.

No live inference and no implicit reuse of LUAD3 sample calibration.
"""
from __future__ import annotations

import copy
import gzip
import hashlib
import importlib.metadata
import json
import math
import tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.stats import rankdata, spearmanr
import shapely
from shapely.geometry import mapping

ROOT = Path(__file__).resolve().parent
RESOURCES = ROOT / "resources"
VERSION = "lung5k-preparation-1.1.0"
COLORS = {"tumor_candidate": "#ed8758", "fibroblast": "#d2bd58", "immune": "#51b7ca",
          "stroma": "#9d8ada", "epithelial_other": "#8baf83", "unknown": "#949ca8"}
GROUP_NAMES = {"tumor_candidate": "Epithelial / tumor candidate", "fibroblast": "Fibroblasts / CAF-like",
               "immune": "Immune", "stroma": "Other stroma / vessels",
               "epithelial_other": "Other epithelium", "unknown": "Unresolved"}


def file_hash(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def fingerprint(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def save_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(obj, indent=2, sort_keys=True, allow_nan=False).encode() + b"\n"
    if path.suffix == ".gz":
        raw = gzip.compress(raw, mtime=0)
    tmp = path.with_name(path.name + ".partial")
    tmp.write_bytes(raw)
    tmp.replace(path)


def load_json(path):
    path = Path(path)
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes())


def environment():
    import platform
    names = ["numpy", "scipy", "pandas", "h5py", "pyarrow", "shapely", "anndata", "decoupler", "numba", "matplotlib"]
    return {"python": platform.python_version(), "platform": platform.platform(),
            "packages": {n: importlib.metadata.version(n) for n in names}}


def load_sample(data_dir, output_dir, verify_hashes=True, min_counts=20, min_genes=10):
    """Read the vendor sparse matrix; join every table by external string ID."""
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    if output_dir.resolve() == data_dir.resolve() or data_dir.resolve() in output_dir.resolve().parents:
        raise ValueError("OUTPUT_DIR must be outside the immutable input folder.")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_json(data_dir / "DOWNLOAD_MANIFEST.json")
    if not manifest.get("complete"):
        raise ValueError("Download manifest is incomplete.")
    required = {"cell_feature_matrix.h5", "cells.parquet", "cell_boundaries.parquet", "analysis.tar.gz"}
    inventory = {}
    for item in manifest["files"]:
        path = data_dir / item["path"]
        if not path.is_file() or path.stat().st_size != item["bytes"]:
            raise ValueError(f"Missing or changed input: {path.name}")
        if verify_hashes and file_hash(path) != item["sha256"]:
            raise ValueError(f"Input checksum mismatch: {path.name}")
        inventory[path.name] = item["sha256"]
    if not required <= inventory.keys():
        raise ValueError("Required original files are absent from the manifest.")
    raw = data_dir / "raw"
    with h5py.File(raw / "cell_feature_matrix.h5", "r") as f:
        m = f["matrix"]
        barcodes = m["barcodes"].asstr()[:]
        names = m["features/name"].asstr()[:]
        gene_ids = m["features/id"].asstr()[:]
        types = m["features/feature_type"].asstr()[:]
        # A CSC feature-by-cell matrix transposes to cell-by-feature CSR.
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                          shape=tuple(m["shape"][:])).T.tocsr()
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("Duplicate matrix cell IDs.")
    keep = types == "Gene Expression"
    X = X[:, keep].astype(np.float32)
    X.eliminate_zeros()
    genes = np.asarray(names[keep], dtype=str)
    if len(set(genes)) != len(genes):
        raise ValueError("Duplicate gene symbols require an explicit mapping; none was applied.")
    if not np.isfinite(X.data).all() or (X.data < 0).any():
        raise ValueError("Counts must be finite and nonnegative.")
    obs = pq.read_table(raw / "cells.parquet").to_pandas().set_index("cell_id", drop=False)
    if not obs.index.is_unique or set(obs.index) != set(barcodes):
        raise ValueError("Cell metadata do not join one-to-one to expression.")
    obs = obs.loc[barcodes].copy()
    obs.index.name = "external_cell_id"
    obs["n_genes_by_counts"] = np.diff(X.indptr)
    obs["counts_gene_expression"] = np.asarray(X.sum(axis=1)).ravel()
    if not np.allclose(obs.counts_gene_expression, obs.transcript_counts):
        raise ValueError("Gene-expression row sums differ from vendor transcript counts.")
    if not np.isfinite(obs[["x_centroid", "y_centroid"]].to_numpy()).all():
        raise ValueError("Nonfinite centroid coordinates.")
    obs["qc_pass"] = (obs.transcript_counts >= min_counts) & (obs.n_genes_by_counts >= min_genes) & (obs.cell_area > 0)
    # Retain failed-QC objects as physical context; never delete their geometry.
    obs["qc_reason"] = np.where(obs.qc_pass, "pass", "low_counts_or_genes_or_invalid_area")
    with tarfile.open(raw / "analysis.tar.gz", "r:gz") as archive:
        member = archive.extractfile("analysis/clustering/gene_expression_graphclust/clusters.csv")
        clusters = pd.read_csv(member, dtype={"Barcode": str}).set_index("Barcode")
        de = pd.read_csv(archive.extractfile("analysis/diffexp/gene_expression_graphclust/differential_expression.csv"))
    if not clusters.index.is_unique or len(set(clusters.index) - set(obs.index)):
        raise ValueError("Invalid clustering identifiers.")
    obs["vendor_cluster"] = clusters.Cluster.reindex(obs.index).fillna(-1).astype(int)
    # IDs are stable under input row reordering and retain an explicit bijection.
    id_map = {s: i + 1 for i, s in enumerate(sorted(barcodes))}
    obs["internal_id"] = [id_map[s] for s in obs.index]
    sample = dict(X=X, obs=obs, genes=genes, gene_ids=gene_ids[keep], data_dir=data_dir,
                  output_dir=output_dir, inventory=inventory, de=de, sample_id=manifest["sample_id"])
    save_json(output_dir / "input_audit.json", dict(sample_id=sample["sample_id"], version=VERSION,
              cells=len(obs), genes=len(genes), qc_pass=int(obs.qc_pass.sum()),
              qc_parameters=dict(min_counts=min_counts, min_genes=min_genes),
              count_layer="original Xenium Q20 cell-assigned gene counts; no SPLIT correction",
              units="micrometres", excluded_feature_types=dict(Counter(types[~keep])),
              input_sha256=inventory, all_input_hashes_verified=verify_hashes))
    print(f"Loaded {len(obs):,} cells × {len(genes):,} genes; {obs.qc_pass.sum():,} pass initial RNA QC.")
    return sample


def annotate(sample, overrides=None):
    """Editable, sample-specific cluster interpretation; not a trained classifier."""
    # Rerunning annotations invalidates references and selections derived from them.
    for key in ["scores", "percentiles", "depth_percentiles", "calibration", "score_annotation_fingerprint", "regions"]:
        sample.pop(key, None)
    records = load_json(RESOURCES / "cluster_annotations.json")
    known = {r["cluster"]: r for r in records}
    if set(sample["obs"].vendor_cluster.unique()) - (set(known) | {-1}):
        raise ValueError("Unknown vendor clusters: this annotation map belongs to one frozen sample.")
    for key, edit in (overrides or {}).items():
        key = int(key)
        if key not in known or set(edit) - {"label", "compartment", "runtime_label", "status"}:
            raise ValueError("Invalid annotation override.")
        known[key] = {**known[key], **edit, "status": edit.get("status", "user_edited")}
        if known[key]["compartment"] not in COLORS:
            raise ValueError("Unknown compartment.")
    obs = sample["obs"]
    for col in ("label", "compartment", "runtime_label", "status"):
        name = {"label": "cell_type_provisional", "status": "annotation_status"}.get(col, col)
        obs[name] = [known.get(int(c), {}).get(col) for c in obs.vendor_cluster]
    obs["compartment"] = obs.compartment.fillna("unknown")
    obs["cell_type_provisional"] = obs.cell_type_provisional.fillna("Unassigned")
    obs["annotation_status"] = obs.annotation_status.fillna("unassigned")
    # Keep marker detection transparent for individual cells, including conflicts.
    gene_index = {g: j for j, g in enumerate(sample["genes"])}
    groups = {
        "T": ["CD3E", "CD3G", "CD2", "CD247", "ZAP70"],
        "epithelial": ["EPCAM", "CDH1", "MUC1", "MAL2", "FOLR1"],
        "fibroblast": ["PDGFRA", "TCF21", "COL5A1", "COL5A2", "DPT", "FAP"],
        "myeloid": ["CSF1R", "CD68", "MRC1", "MARCO", "LILRB2"],
    }
    for name, markers in groups.items():
        indices = [gene_index[g] for g in markers if g in gene_index]
        obs[f"markers_{name}_detected"] = np.asarray((sample["X"][:, indices] > 0).sum(axis=1)).ravel()
    obs["mixed_marker_flag"] = (obs.markers_epithelial_detected >= 2) & ((obs.markers_T_detected >= 2) | (obs.markers_myeloid_detected >= 2))
    table = []
    for c, row in sorted(known.items()):
        mask = obs.vendor_cluster == c
        prefix = f"Cluster {c} "
        de = sample["de"]
        eligible = de[(de[prefix + "Mean Counts"] > .05) & (de[prefix + "Adjusted p value"] < .01)]
        top = eligible.sort_values(prefix + "Log2 fold change", ascending=False)["Feature Name"].head(12).tolist()
        table.append({**row, "n_cells": int(mask.sum()), "markers": ", ".join(row["markers"]), "top_genes": ", ".join(top),
                      "mixed_marker_fraction": round(float(obs.loc[mask, "mixed_marker_flag"].mean()), 4)})
    table = pd.DataFrame(table)
    table.to_csv(sample["output_dir"] / "annotation_review.csv", index=False)
    sample["annotation_table"] = table
    sample["annotation_fingerprint"] = fingerprint(known)
    print("Provisional compartments:", obs.compartment.value_counts().to_dict())
    return table


def score_aucell(sample, n_up=100, min_overlap=5, batch_size=4096):
    """Official decoupler AUCell, panel-only universe, no zero-ranked scores."""
    import decoupler as dc
    from anndata import AnnData
    if dc.__version__ != "2.3.0":
        raise ValueError("This notebook was tested with decoupler==2.3.0.")
    import numba
    numba.set_num_threads(min(2, numba.get_num_threads()))
    definitions = load_json(RESOURCES / "go20_msigdb_2026.1.Hs.json")
    genes = set(sample["genes"])
    rows, networks = [], []
    for name, definition in sorted(definitions.items()):
        full = set(definition["geneSymbols"])
        overlap = sorted(genes & full)
        rows.append(dict(signature=name, go_id=definition["exactSource"], full_genes=len(full), panel_genes=len(overlap),
                         panel_coverage=len(overlap) / len(full), available=len(overlap) >= min_overlap,
                         measured_genes=";".join(overlap)))
        if len(overlap) >= min_overlap:
            networks.extend(dict(source=name, target=g) for g in overlap)
    coverage = pd.DataFrame(rows).set_index("signature")
    coverage.to_csv(sample["output_dir"] / "signature_coverage.csv")
    if not 2 <= int(n_up) < len(genes):
        raise ValueError("Invalid AUCell rank cutoff.")
    obs = sample["obs"]
    valid = obs.qc_pass.to_numpy() & (obs.n_genes_by_counts.to_numpy() >= n_up)
    names = sorted(coverage.index[coverage.available])
    scores = np.full((len(obs), len(names)), np.nan, dtype=np.float32)
    network = pd.DataFrame(networks)
    use = np.flatnonzero(valid)
    for start in range(0, len(use), batch_size):
        idx = use[start:start + batch_size]
        batch = AnnData(sample["X"][idx], obs=pd.DataFrame(index=obs.index[idx]),
                        var=pd.DataFrame(index=sample["genes"]))
        dc.mt.aucell(batch, network, tmin=min_overlap, n_up=int(n_up), empty=False,
                     bsize=batch_size, verbose=False)
        result = batch.obsm["score_aucell"]
        scores[idx] = result.reindex(index=obs.index[idx], columns=names).to_numpy(dtype=np.float32)
        if start == 0 or (start // batch_size) % 10 == 0:
            print(f"AUCell: {min(start + len(idx), len(use)):,}/{len(use):,} eligible cells", flush=True)
    if not np.isfinite(scores[valid]).all() or (scores[valid] < -1e-6).any() or (scores[valid] > 1 + 1e-6).any():
        raise ValueError("Invalid AUCell output.")
    percentiles = np.full_like(scores, np.nan)
    depth_percentiles = np.full_like(scores, np.nan)
    obs["aucell_available"] = valid
    obs["score_depth_decile"] = -1
    obs["depth_reference_n"] = 0
    calibration = {}
    for label in sorted(obs.cell_type_provisional.unique()):
        indices = np.flatnonzero(valid & (obs.cell_type_provisional.to_numpy() == label))
        if len(indices) < 20:
            continue
        depth = obs.n_genes_by_counts.to_numpy()[indices]
        edges = np.quantile(depth, np.linspace(0, 1, 11))
        bins = np.searchsorted(edges[1:-1], depth, side="right")
        percentiles[indices] = rankdata(scores[indices], method="average", axis=0) / len(indices)
        for b in range(10):
            group = indices[bins == b]
            if len(group) >= 30:
                depth_percentiles[group] = rankdata(scores[group], method="average", axis=0) / len(group)
                # A constant signature has no identifiable relative ordering.
                constant = np.ptp(scores[group], axis=0) == 0
                depth_percentiles[np.ix_(group, np.flatnonzero(constant))] = np.nan
            obs.iloc[group, obs.columns.get_loc("score_depth_decile")] = b
            obs.iloc[group, obs.columns.get_loc("depth_reference_n")] = len(group)
        calibration[label] = dict(n_reference=len(indices), depth_edges=edges.tolist(),
                                  scores={s: dict(quantiles=np.quantile(scores[indices, j], np.linspace(0, 1, 21)).tolist()) for j, s in enumerate(names)})
    # Missing depth-conditioned references remain unavailable; no silent fallback.
    sample.update(scores=scores, percentiles=percentiles, depth_percentiles=depth_percentiles,
                  signature_names=names, coverage=coverage, n_up=int(n_up), calibration=calibration)
    sample["score_annotation_fingerprint"] = sample["annotation_fingerprint"]
    audit = dict(method="AUCell", implementation="decoupler.mt.aucell", implementation_version=dc.__version__,
                 n_up=int(n_up), min_panel_overlap=min_overlap, scored_cells=int(valid.sum()),
                 unscored_cells=int((~valid).sum()), universe="all 5001 assayed Gene Expression features, including undetected genes",
                 tie_policy="decoupler 2.3.0 seeded feature permutation (seed=0), followed by ordinal ranks; deterministic, not R AUCell random per-cell ties",
                 zero_policy="Cells with fewer than n_up detected genes retain missing scores; zero-count genes cannot enter a scored cell's top ranks.",
                 normalisation="Original counts; within-cell ranks are invariant to library-size scaling and log1p. No cross-gene scaling.",
                 reference="all eligible cells of this sample, by provisional type; additionally by within-type detection-depth decile (minimum n=30); constant signatures have no relative percentile",
                 annotation_fingerprint=sample["annotation_fingerprint"],
                 caveat="Relative expression-signature scores, not measured pathway activation; depth conditioning is descriptive, not causal correction.",
                 gene_set_definition_sha256=file_hash(RESOURCES / "go20_msigdb_2026.1.Hs.json"))
    save_json(sample["output_dir"] / "aucell_provenance.json", audit)
    save_json(sample["output_dir"] / "score_calibration.json", {"provenance": audit, "types": calibration})
    obs.to_parquet(sample["output_dir"] / "whole_sample_annotations.parquet")
    np.savez_compressed(sample["output_dir"] / "whole_sample_aucell.npz", external_cell_ids=obs.index.to_numpy(dtype=str),
                        signature_names=np.asarray(names), aucell=scores, within_type=percentiles, within_type_depth=depth_percentiles)
    return coverage


def read_geometry(sample):
    """Keep original outlines and repaired polygonal components separately."""
    tab = pq.read_table(sample["data_dir"] / "raw/cell_boundaries.parquet").to_pandas()
    ids = tab.cell_id.to_numpy()
    starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])
    ends = np.r_[starts[1:], len(ids)]
    names = ids[starts]
    if len(set(names)) != len(names):
        raise ValueError("Non-contiguous or multipart vendor boundaries need an explicit importer.")
    coords = np.ascontiguousarray(tab[["vertex_x", "vertex_y"]].to_numpy(dtype=float))
    if not np.isfinite(coords).all() or (ends - starts < 4).any():
        raise ValueError("Invalid polygon vertices.")
    rings = np.r_[starts, len(ids)]
    polygons = shapely.from_ragged_array(shapely.GeometryType.POLYGON, coords, (rings, np.arange(len(names) + 1)))
    order = pd.Index(names).get_indexer(sample["obs"].index)
    if (order < 0).any():
        raise ValueError("Some measured cells have no boundary; refusing to invent geometry.")
    original = polygons[order]
    valid = shapely.is_valid(original)
    repaired = original.copy()
    repaired[~valid] = shapely.make_valid(original[~valid])
    usable = ~shapely.is_empty(repaired) & np.isin(shapely.get_type_id(repaired), [3, 6]) & (shapely.area(repaired) > 0)
    sample.update(polygons=repaired, original_polygons=original, polygon_repaired=~valid,
                  geometry_usable=usable, tree=shapely.STRtree(repaired))
    sample["obs"]["geometry_usable"] = usable
    sample["obs"]["polygon_repaired"] = ~valid
    save_json(sample["output_dir"] / "geometry_audit.json", dict(cells=len(original), repaired=int((~valid).sum()),
              usable=int(usable.sum()), unhandled=int((~usable).sum()), source="vendor simplified cell polygons, in micrometres",
              masks_available="raw/cells.zarr.zip; not reconstructed or substituted by polygons"))
    print(f"Geometry: {usable.sum():,} usable original/repaired outlines; {(~valid).sum():,} repairs recorded.")


def make_region(sample, center, core_target=750, max_total=1000, halo_um=30):
    """Contiguous nearest-neighbour disc; no cherry-picking of individual labels."""
    if not 100 <= core_target <= max_total <= 1000 or not np.isfinite(halo_um) or halo_um < 30:
        raise ValueError("Require 100 <= core_target <= max_total <= 1000 and finite halo >= 30 um.")
    xy = sample["obs"][["x_centroid", "y_centroid"]].to_numpy()
    distances = np.linalg.norm(xy - np.asarray(center), axis=1)
    nearest = np.argsort(distances, kind="stable")[:max_total]
    # Include all centroids up to the selected radius, including failed-QC objects.
    target = int(core_target)
    while target >= 100:
        radius = float(distances[nearest[target - 1]])
        core = np.flatnonzero(distances <= radius + 1e-9)
        if not sample["geometry_usable"][core].all():
            return None
        envelope = shapely.buffer(shapely.union_all(sample["polygons"][core]), halo_um)
        selected = np.unique(np.r_[core, sample["tree"].query(envelope, predicate="intersects")]).astype(int)
        if len(selected) <= max_total and sample["geometry_usable"][selected].all():
            return dict(center_um=list(map(float, center)), radius_um=radius, core=core,
                        selected=selected, halo=np.setdiff1d(selected, core), halo_um=halo_um)
        target -= max(10, int((len(selected) - max_total) * .8))
    return None


def region_statistics(sample, region):
    obs = sample["obs"].iloc[region["core"]]
    counts = obs.compartment.value_counts().to_dict()
    n = len(obs)
    fractions = {g: counts.get(g, 0) / n for g in COLORS}
    targets = np.asarray([.35, .20, .30, .10, .03, .02])
    f = np.asarray([fractions[g] for g in COLORS])
    # Selection uses t0 composition, never simulated killing or final tumour count.
    quality = 1 - float(np.abs(f - targets).sum()) / 2
    good = fractions["tumor_candidate"] >= .15 and fractions["fibroblast"] >= .08 and fractions["immune"] >= .15
    usable = obs.qc_pass & ~obs.mixed_marker_flag
    score_available = obs.get("aucell_available", pd.Series(False, index=obs.index))
    return dict(n_core=n, n_halo=len(region["halo"]), n_total=len(region["selected"]),
                radius_um=round(region["radius_um"], 1), balance_score=round(quality, 4), balance_pass=good,
                qc_fraction=round(float(obs.qc_pass.mean()), 3), mixed_fraction=round(float(obs.mixed_marker_flag.mean()), 3),
                aucell_fraction=round(float(score_available.mean()), 3),
                cd8_provisional=int((obs.runtime_label == "T_CD8_cytotoxic").sum()),
                macrophage_provisional=int((obs.runtime_label == "Macrophage").sum()),
                supported_qc_core=int((obs.runtime_label.notna() & obs.qc_pass & ~obs.mixed_marker_flag).sum()),
                **{g: counts.get(g, 0) for g in COLORS})


def propose_regions(sample, n_options=8, core_target=750, max_total=1000, halo_um=30, grid_um=180):
    """Find balanced discs, then retain spatially different alternatives."""
    obs = sample["obs"]
    xy = obs[["x_centroid", "y_centroid"]].to_numpy()
    tumor_xy = xy[(obs.compartment == "tumor_candidate").to_numpy() & obs.qc_pass.to_numpy()]
    if not len(tumor_xy):
        raise ValueError("No tumour-candidate epithelium; review annotation before proposing tumour interfaces.")
    grid = np.unique(np.floor(tumor_xy / grid_um).astype(int), axis=0)
    centers = (grid + .5) * grid_um
    kd = cKDTree(xy)
    _, idx = kd.query(centers, k=min(core_target, len(obs)))
    compartments = obs.compartment.to_numpy()
    quick = []
    for c, inds in zip(centers, idx):
        cnt = Counter(compartments[inds])
        f = np.array([cnt.get(g, 0) / len(inds) for g in COLORS])
        if f[0] >= .10 and f[1] >= .05 and f[2] >= .10:
            quick.append((1 - np.abs(f - [.35, .20, .30, .10, .03, .02]).sum()/2, c))
    choices = []
    for _, center in sorted(quick, key=lambda x: -x[0]):
        if any(np.linalg.norm(center - np.asarray(r["center_um"])) < .8 * r["radius_um"] for r in choices):
            continue
        r = make_region(sample, center, core_target, max_total, halo_um)
        if r is None:
            continue
        stats = region_statistics(sample, r)
        if not stats["balance_pass"] or stats["qc_fraction"] < .8:
            continue
        # Prevent near-duplicates by exact cell overlap (Jaccard).
        if any(len(np.intersect1d(r["core"], x["core"])) / len(np.union1d(r["core"], x["core"])) > .3 for x in choices):
            continue
        r.update(stats)
        choices.append(r)
        if len(choices) >= n_options:
            break
    if not choices:
        raise ValueError("No region meets the declared balance/QC criteria. Review annotations or use a manual centre; criteria were not silently relaxed.")
    choices.sort(key=lambda r: -r["balance_score"])
    for i, r in enumerate(choices, 1):
        r["region_id"] = f"R{i:02d}"
    sample["regions"] = choices
    table = pd.DataFrame([{k: v for k, v in r.items() if k not in ("core", "selected", "halo", "center_um")} |
                         dict(x_um=r["center_um"][0], y_um=r["center_um"][1]) for r in choices])
    table.to_csv(sample["output_dir"] / "candidate_regions.csv", index=False)
    save_json(sample["output_dir"] / "region_selection_rules.json", dict(max_total_cells=max_total, core_target=core_target,
              halo_um=halo_um, min_tumor_candidate_fraction=.15, min_fibroblast_fraction=.08, min_immune_fraction=.15,
              maximum_core_jaccard=.3, selection_basis="Only measured t0 geometry, provisional composition and RNA QC; no simulated outcomes.",
              annotation_fingerprint=sample["annotation_fingerprint"]))
    return table


def plot_regions(sample, region_ids=None, save=True):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Patch, Circle
    choices = sample["regions"]
    if region_ids is not None:
        choices = [r for r in choices if r["region_id"] in region_ids]
    n = len(choices)
    columns = min(n, 3)
    fig, axes = plt.subplots(math.ceil(n / columns), columns,
                             figsize=(6 if n == 1 else 16, 5.6 * math.ceil(n / columns)), squeeze=False)
    obs = sample["obs"]
    for ax, r in zip(axes.flat, choices):
        vertices, colors = [], []
        core_set = set(r["core"])
        for i in r["selected"]:
            geom = sample["polygons"][i]
            parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
            for part in parts:
                vertices.append(np.asarray(part.exterior.coords))
                colors.append(COLORS[obs.iloc[i].compartment] if i in core_set else "#e2e5e9")
        ax.add_collection(PolyCollection(vertices, facecolors=colors, edgecolors="#ffffff", linewidths=.12))
        ax.autoscale_view();ax.set_aspect("equal");ax.invert_yaxis()
        ax.add_patch(Circle(r["center_um"], r["radius_um"], fill=False, edgecolor="#26344a", linestyle="--", lw=1))
        cx, cy = r["center_um"]; radius = r["radius_um"] + r["halo_um"] + 20
        ax.set_xlim(cx-radius, cx+radius); ax.set_ylim(cy+radius, cy-radius)
        ax.set_title(f"{r['region_id']} · {r['n_core']} core + {r['n_halo']} context = {r['n_total']}\n"
                     f"Tumor candidates {r['tumor_candidate']} · Fibroblasts {r['fibroblast']} · Immune {r['immune']}", fontsize=10)
        ax.set_xlabel("X (µm)");ax.set_ylabel("Y (µm)")
    for ax in list(axes.flat)[n:]:
        ax.set_visible(False)
    fig.legend(handles=[Patch(facecolor=c, label=GROUP_NAMES[g]) for g, c in COLORS.items()], loc="lower center", ncol=3)
    fig.suptitle("Choose a measured region · provisional cell identities · grey = surrounding context", fontsize=14)
    fig.tight_layout(rect=[0,.055,1,.96])
    if save:
        fig.savefig(sample["output_dir"] / "candidate_regions.png", dpi=170, bbox_inches="tight")
    return fig


def plot_overview(sample, save=True):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    obs = sample["obs"]
    fig, ax = plt.subplots(figsize=(13, 9))
    for group, color in COLORS.items():
        part = obs[obs.compartment == group]
        ax.scatter(part.x_centroid, part.y_centroid, s=.4, color=color, alpha=.65, label=GROUP_NAMES[group], rasterized=True)
    for r in sample.get("regions", []):
        ax.add_patch(Circle(r["center_um"], r["radius_um"], fill=False, color="#17243a", lw=1.4))
        ax.text(*r["center_um"], r["region_id"], fontsize=10, weight="bold", bbox=dict(facecolor="white", alpha=.9, edgecolor="none"))
    ax.set_aspect("equal");ax.invert_yaxis();ax.legend(markerscale=8, loc="upper right")
    ax.set(xlabel="X (µm)", ylabel="Y (µm)", title="Whole Xenium section and candidate regions")
    if save:
        fig.savefig(sample["output_dir"] / "whole_section_regions.png", dpi=170, bbox_inches="tight")
    return fig


def new_sample_contract(sample):
    """Recalculate the structural proxy; do not transfer the LUAD3 calibration."""
    contract = load_json(RESOURCES / "runtime_template.json")
    obs = sample["obs"]
    ratios = (obs.nucleus_area / obs.cell_area).to_numpy(dtype=float)
    valid = np.isfinite(ratios) & (ratios > 0) & (ratios <= 1) & obs.qc_pass.to_numpy()
    if valid.sum() < 30:
        raise ValueError("Insufficient valid nuclear areas for structural initialization.")
    global_q10 = float(np.quantile(ratios[valid], .1))
    stats, floors = {}, {}
    supported = contract["mechanics"]["anchoring"]["class_of_label"]
    for label in supported:
        if label in {"occluded", "Tumor_proliferating"}:
            continue
        mask = valid & (obs.runtime_label.to_numpy() == label)
        values = ratios[mask]
        enough = len(values) >= 30
        q10 = float(np.quantile(values, .1)) if enough else global_q10
        floors[label] = float(np.clip(q10, .1, .8))
        stats[label] = dict(n=len(values), q10=q10, reference="same_sample_type" if enough else "same_sample_global_fallback")
    constraint = contract["mechanics"]["compression_constraint"]
    constraint["area_floor_fraction_by_label"] = floors
    constraint["fallback_fraction"] = float(np.clip(global_q10, .1, .8))
    contract["calibration"] = dict(sample_id=sample["sample_id"], reference="whole new sample, measured nuclear/cell area ratio",
        source_sha256=sample["inventory"]["cells.parquet"], annotation_fingerprint=sample["annotation_fingerprint"],
        n_valid=int(valid.sum()), ratio_statistics=stats, global_q10=global_q10,
        scope="2D structural proxy; neither physical stiffness nor kinetic rates are identified by one static section.")
    contract["status"] = "New-sample structural proxy; exploratory dynamics and provisional identities, no validated tumour trajectory."
    return contract


def export_region(sample, region_id, tumor_mode="exploratory_candidates", destination=None):
    """Export raw counts, XY, polygons, AUCell and an explicitly adapted atlas."""
    import anndata as ad
    import shutil
    if sample.get("score_annotation_fingerprint") != sample.get("annotation_fingerprint"):
        raise ValueError("Annotations changed or AUCell has not been calculated. Recompute scores and regions before export.")
    if tumor_mode not in {"exploratory_candidates", "context_until_review"}:
        raise ValueError("Choose exploratory_candidates or context_until_review.")
    region = next(r for r in sample["regions"] if r["region_id"] == region_id)
    if len(region["selected"]) > 1000:
        raise ValueError("Region exceeds the total cell limit, including context.")
    selected = region["selected"]
    obs = sample["obs"].iloc[selected].copy()
    core = set(region["core"])
    contract = new_sample_contract(sample)
    supported = set(contract["mechanics"]["anchoring"]["class_of_label"]) - {"occluded", "Tumor_proliferating"}
    path = Path(destination) if destination else sample["output_dir"] / f"simulation_input_{region_id}_{tumor_mode}"
    if path.exists():
        raise FileExistsError(f"Preserve the existing export; choose a new destination: {path}")
    if sample["data_dir"].resolve() == path.resolve() or sample["data_dir"].resolve() in path.resolve().parents:
        raise ValueError("Export must be outside the input folder.")
    path.mkdir(parents=True)
    cells, roles, reasons, labels = [], [], [], []
    for i in selected:
        row = sample["obs"].iloc[i]
        label = row.runtime_label
        reasons_i = []
        if not row.qc_pass: reasons_i.append("RNA_QC_failed")
        if row.mixed_marker_flag: reasons_i.append("mixed_lineage_markers")
        if label not in supported: reasons_i.append("no_supported_runtime_identity")
        if row.compartment == "tumor_candidate" and tumor_mode == "context_until_review":
            reasons_i.append("tumour_identity_not_confirmed")
        enabled = not reasons_i
        role = ("roi" if enabled else "occluded") if i in core else "halo"
        # Unknown labels remain physical context; do not assign a false lineage.
        identity = dict(candidates=[dict(label=str(label), lineage=str(row.compartment), source="provisional_vendor_cluster_interpretation")],
                        agreement="provisional_single_source") if enabled else None
        geom = sample["polygons"][i]
        signatures = {}
        for j, name in enumerate(sample["signature_names"]):
            raw = sample["scores"][i, j]
            pct = sample["depth_percentiles"][i, j]
            # The runtime receives evidence only with a sufficient reference group.
            if np.isfinite(raw) and np.isfinite(pct):
                signatures[name] = dict(method="AUCell", raw=float(raw),
                    percentile_within_type=float(sample["percentiles"][i,j]), percentile_within_type_depth=float(pct),
                    coverage_flag="targeted_panel_relative_expression", panel_coverage=float(sample["coverage"].loc[name,"panel_coverage"]))
        expression = sample["X"].getrow(i)
        cell = dict(cell_id=int(row.internal_id), external_cell_id=str(row.cell_id), role=role,
                    geometry=mapping(geom), area_um2=float(geom.area), x=float(row.x_centroid), y=float(row.y_centroid),
                    polygon_repaired=bool(row.polygon_repaired), identity=identity,
                    annotation=dict(provisional_type=row.cell_type_provisional, compartment=row.compartment,
                                    vendor_cluster=int(row.vendor_cluster), status=row.annotation_status,
                                    context_reasons=reasons_i, tumor_mode=tumor_mode),
                    genes_detected=sorted(sample["genes"][expression.indices].tolist()),
                    panel_genes=sample["genes"].tolist(), signatures=signatures, contacts=[], neighbour_ids_30um=[],
                    qc=dict(scored=bool(row.aucell_available), qc_pass=bool(row.qc_pass),
                            n_genes_by_counts=int(row.n_genes_by_counts), transcript_counts=int(row.transcript_counts),
                            mixed_marker_flag=bool(row.mixed_marker_flag), depth_reference_n=int(row.depth_reference_n)),
                    region=dict(region_id=region_id, geometric_core=bool(i in core)),
                    observation_provenance="Original Xenium t0 counts; AUCell (decoupler 2.3.0), all 5001 assayed genes, rank cutoff %s. Percentiles conditional on provisional type and RNA detection-depth decile, reference entire sample. Missing evidence is not zero. No observed malignancy or functional competence." % sample["n_up"])
        cells.append(cell);roles.append(role);reasons.append(";".join(reasons_i));labels.append(label if enabled else "occluded")
    if not any(c["role"] == "roi" for c in cells):
        raise ValueError("No supported acting cells remain; review the chosen region and annotations.")
    shapes = sample["polygons"][selected]
    tree = shapely.STRtree(shapes)
    xy = obs[["x_centroid","y_centroid"]].to_numpy()
    kd = cKDTree(xy)
    for j, cell in enumerate(cells):
        cell["neighbour_ids_30um"] = [cells[k]["cell_id"] for k in kd.query_ball_point(xy[j],30) if k != j]
        for k in tree.query(shapely.buffer(shapes[j],1.01)):
            if j == k: continue
            gap = float(shapely.distance(shapes[j],shapes[k]))
            if gap > 1: continue
            interface = shapely.area(shapely.intersection(shapely.buffer(shapes[j],.75,quad_segs=2),
                                                         shapely.buffer(shapes[k],.75,quad_segs=2))) / 1.5
            cell["contacts"].append(dict(cell_id=cells[k]["cell_id"],gap_um=gap,interface_um=float(interface),
                                         outside_patch=cells[k]["role"] != "roi",occluded=cells[k]["identity"] is None))
        cell["contacts"].sort(key=lambda e:e["cell_id"])
    atlas = dict(schema_version="1.0.0", preparation_schema="xenium-aucell-input-1.0.0",
                 patch_id=f"{sample['sample_id']}_{region_id}",cells=cells,
                 geometry=dict(center_um=region["center_um"],radius_um=region["radius_um"],halo_um=region["halo_um"],
                               roi_bounds_um=list(shapely.union_all(sample["polygons"][region["core"]]).bounds),
                               coordinate_source="Original Xenium centroids and simplified polygons, micrometres",dimension=2),
                 summary=dict(n_total=len(cells),n_geometric_core=len(core),n_roi=roles.count("roi"),
                              n_context=len(cells)-roles.count("roi"),roles=dict(Counter(roles)),
                              composition=dict(Counter(c["identity"]["candidates"][0]["label"] for c in cells if c["role"]=="roi")),
                              provisional_core_composition=region_statistics(sample,region)),
                 provenance=dict(sample_id=sample["sample_id"],source_url=load_json(sample["data_dir"]/'SOURCE.json').get('source_page'),
                    input_sha256=sample["inventory"],annotation_fingerprint=sample["annotation_fingerprint"],
                    preparation_version=VERSION,preparation_source_sha256=file_hash(Path(__file__)),score_method="AUCell",tumor_mode=tumor_mode,
                    expression_layer="original_raw_counts",contact_definition="distance between vendor simplified cell polygons <= 1 micrometre",
                    contact_interpretation="geometric proximity, not demonstrated receptor binding",
                    runtime_contract_fingerprint=fingerprint(contract),
                    source_hashes_verified=load_json(sample["output_dir"]/'input_audit.json')['all_input_hashes_verified'],
                    reference_scope="New 10x sample; not the legacy LUAD3 reference"),
                 signatures={name:dict(method="AUCell",n_up=sample["n_up"],universe_size=len(sample["genes"])) for name in sample["signature_names"]})
    atlas["fingerprint"] = fingerprint(atlas)
    save_json(path/'atlas.json.gz',atlas)
    save_json(path/'runtime_contract.json',contract)
    obs['simulation_role']=roles;obs['simulation_label']=labels;obs['context_reasons']=reasons
    obs['in_geometric_core']=[i in core for i in selected]
    obs['tumor_mode']=tumor_mode
    # H5AD keeps missing numeric observations as NaN, not fabricated zeros.
    for col in obs.select_dtypes(include=['object','string']).columns:
        obs[col]=obs[col].fillna('unavailable').astype(str)
    obs.index=obs.index.astype(str)
    var=pd.DataFrame(dict(gene_id=sample['gene_ids']),index=pd.Index(sample['genes'],name='gene_symbol'))
    data=ad.AnnData(X=sample['X'][selected].copy(),obs=obs,var=var)
    data.obsm['spatial']=xy.astype(float)
    for key,values in [('AUCell',sample['scores']),('AUCell_pct_within_type',sample['percentiles']),('AUCell_pct_within_type_depth',sample['depth_percentiles'])]:
        data.obsm[key]=pd.DataFrame(values[selected],index=obs.index,columns=sample['signature_names'])
    data.uns['preparation']=dict(version=VERSION,sample_id=sample['sample_id'],region_id=region_id,
        coordinate_units='micrometres',X_semantics='Original unnormalized Xenium gene counts; no SPLIT correction',
        annotation='Provisional cluster interpretation; not malignant-cell validation',tumor_mode=tumor_mode,
        signatures=sample['signature_names'],aucell_n_up=sample['n_up'],atlas_fingerprint=atlas['fingerprint'])
    data.write_h5ad(path/'region.h5ad',compression='gzip')
    obs.to_csv(path/'cells.csv',index=False)
    pd.DataFrame(dict(cell_id=obs.internal_id,external_cell_id=obs.cell_id,
        original_wkb=shapely.to_wkb(sample['original_polygons'][selected]),
        simulation_wkb=shapely.to_wkb(shapes),repaired=sample['polygon_repaired'][selected])).to_parquet(path/'geometry.parquet',index=False)
    for kind in ['cell','nucleus']:
        tab=pq.read_table(sample['data_dir']/f'raw/{kind}_boundaries.parquet',filters=[('cell_id','in',obs.cell_id.tolist())]).to_pandas()
        tab.to_parquet(path/f'{kind}_boundaries.parquet',index=False)
    sample['coverage'].to_csv(path/'signature_coverage.csv')
    sample['annotation_table'].to_csv(path/'annotation_review.csv',index=False)
    for name in ['score_calibration.json','aucell_provenance.json','input_audit.json','geometry_audit.json','region_selection_rules.json']:
        shutil.copy2(sample['output_dir']/name,path/name)
    for name in ['SOURCE.json','DOWNLOAD_MANIFEST.json']:
        shutil.copy2(sample['data_dir']/name,path/name)
    shutil.copy2(RESOURCES/'runtime_origin.json',path/'runtime_origin.json')
    shutil.copy2(ROOT/'prepared_runtime.py',path/'prepared_runtime.py')
    shutil.copytree(ROOT/'xenium_runtime',path/'xenium_runtime',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    if (ROOT/'requirements-colab.txt').exists(): shutil.copy2(ROOT/'requirements-colab.txt',path/'requirements-colab.txt')
    (path/'START_HERE.txt').write_text(
        'Prepared Xenium region '+region_id+'\n\n'
        'region.h5ad: raw counts in X; XY in obsm/spatial (micrometres); AUCell and reference percentiles in obsm.\n'
        'cell_id is the original Xenium string; internal_id is a reversible numeric simulation identifier.\n'
        'Original cell/nucleus boundaries are retained; geometry.parquet also records original and repaired WKB.\n'
        'atlas.json.gz and runtime_contract.json are inputs for the bundled prepared_runtime loader.\n'
        'The older Colab_Prototype notebook hard-codes LUAD3 and UCell. Do not replace its files or bypass its reference checks.\n'
        'Instead, add this export folder to sys.path, then from prepared_runtime import load_inputs, initialize.\n'
        'initialize(folder, scenario="static") constructs the initial state; scenario="S0" uses exploratory local rules.\n'
        'No HTTP provider or paid calls are included. New sample structural references replace LUAD3 statistics.\n'
        'This bundle contains provisional identities. Tumour candidate handling: '+tumor_mode+'.\n'
        'Malignancy, antigen recognition, receptor protein competence and kinetic rates are not established by these data.\n'
        'The local policy is a programmed test, not a trained biological agent.\n',encoding='utf8')
    write_bundle_manifest(path)
    print(f'Exported {path}: {len(cells)} total, {roles.count("roi")} acting, {len(cells)-roles.count("roi")} context.')
    return path


def write_bundle_manifest(path):
    path=Path(path)
    files={p.relative_to(path).as_posix():dict(sha256=file_hash(p),bytes=p.stat().st_size)
           for p in sorted(path.rglob('*')) if p.is_file() and p.name!='BUNDLE_MANIFEST.json' and '__pycache__' not in p.parts}
    save_json(path/'BUNDLE_MANIFEST.json',dict(format='xenium-aucell-bundle-1',files=files))


def zip_bundle(path):
    import zipfile
    path=Path(path)
    write_bundle_manifest(path)
    dest=path.with_suffix('.zip')
    with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(path.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:
                z.write(p,arcname=(Path(path.name)/p.relative_to(path)).as_posix())
    return dest


def plot_marker_review(sample):
    """Visual audit of measured markers; no external classifier or imputation."""
    import matplotlib.pyplot as plt
    groups={
        'Epithelium':['EPCAM','CDH1','MUC1','MAL2','FOLR1'],
        'Fibroblast':['PDGFRA','TCF21','FAP','COL5A1','DPT'],
        'T / NK':['CD3E','CD3G','CD4','CD8A','KLRD1','PRF1'],
        'Myeloid':['CSF1R','CD68','MRC1','MARCO'],
        'B / plasma':['MS4A1','CD79A','MZB1'],
        'Endothelial':['CLDN5','CDH5','ROBO4'],
    }
    genes=[g for v in groups.values() for g in v if g in sample['genes']]
    columns=pd.Index(sample['genes']).get_indexer(genes)
    detected=(sample['X'][:,columns]>0).astype(float)
    rows=[];labels=[]
    for c in sorted(sample['obs'].vendor_cluster.unique()):
        mask=(sample['obs'].vendor_cluster==c).to_numpy()
        rows.append(np.asarray(detected[mask].mean(axis=0)).ravel())
        labels.append(f"{c:02d} · {sample['obs'].loc[mask,'cell_type_provisional'].iloc[0]}")
    fig,ax=plt.subplots(figsize=(13,12))
    im=ax.imshow(np.asarray(rows),aspect='auto',vmin=0,vmax=1,cmap='Blues')
    ax.set_xticks(range(len(genes)),genes,rotation=90);ax.set_yticks(range(len(rows)),labels,fontsize=9)
    ax.set_title('Marker detection by vendor cluster · provisional interpretation')
    fig.colorbar(im,ax=ax,label='Fraction of cells with at least one detected transcript')
    fig.tight_layout();fig.savefig(sample['output_dir']/'annotation_marker_review.png',dpi=160)
    return fig


def manual_region(sample,x_um,y_um,core_target=750,max_total=1000,halo_um=30):
    if not np.isfinite([x_um,y_um]).all(): raise ValueError('Use finite coordinates in micrometres.')
    region=make_region(sample,[float(x_um),float(y_um)],core_target,max_total,halo_um)
    if region is None: raise ValueError('No usable region fits the total cell limit here; choose another centre.')
    region.update(region_statistics(sample,region))
    region['region_id']='MANUAL'
    sample['regions']=[r for r in sample['regions'] if r['region_id']!='MANUAL']+[region]
    return region
