"""Pinned public data downloads and strict annotated Xenium input adapter.

Only small individual files are downloaded. Annotation is an investigator input,
never inferred from vendor cluster numbers. Raw matrices stay sparse on disk;
only requested cell columns of the 10x genes-by-cells CSC matrix are read.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import tarfile
import urllib.request


def list_datasets():
    """Return independent copies of the packaged, checksum-pinned source records."""
    return json.loads(Path(__file__).with_name("dataset_manifests.json").read_text())


def get_dataset(dataset_id):
    for record in list_datasets():
        if record["id"] == dataset_id:
            return deepcopy(record)
    raise ValueError(f"Unknown dataset {dataset_id!r}; choose {[x['id'] for x in list_datasets()]}")


def _sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def verify_dataset(dataset_id, data_dir, *, include_optional=False):
    """Check size AND SHA-256 for every requested file. No network is used."""
    dataset = get_dataset(dataset_id)
    root = Path(data_dir) / "raw"
    results = []
    for item in dataset["files"]:
        if item.get("optional") and not include_optional:
            continue
        path = root / item["path"]
        ok = path.is_file() and path.stat().st_size == item["bytes"] and _sha256(path) == item["sha256"]
        results.append(dict(path=item["path"], ok=ok))
    return dict(dataset_id=dataset_id, ok=all(x["ok"] for x in results), files=results)


def download_dataset(dataset_id, data_dir, *, max_bytes=100 * 2**20, include_optional=False):
    """Download size-bounded, pinned files atomically; reuse only matching files.

    max_bytes bounds the total selected dataset, including cached files. Every
    response is bounded independently even if the server omits Content-Length.
    Failed/oversized downloads never become final files or complete receipts.
    """
    dataset = get_dataset(dataset_id)
    files = [x for x in dataset["files"] if include_optional or not x.get("optional")]
    if type(max_bytes) is not int or max_bytes < sum(x["bytes"] for x in files):
        raise ValueError("Selected files exceed max_bytes; increase it explicitly after inspecting the manifest.")
    root = Path(data_dir)
    (root / "raw").mkdir(parents=True, exist_ok=True)
    # Never leave a stale success receipt if a changed/corrupt cached file fails.
    (root / "DOWNLOAD_RECEIPT.json").unlink(missing_ok=True)
    for item in files:
        path = root / "raw" / item["path"]
        if path.is_file() and path.stat().st_size == item["bytes"] and _sha256(path) == item["sha256"]:
            continue
        temporary = path.with_name(path.name + ".partial")
        try:
            request = urllib.request.Request(item["url"], headers={"User-Agent": "virtualtissue-xenium/2"})
            h, n = hashlib.sha256(), 0
            with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as stream:
                while True:
                    block = response.read(min(2**20, item["bytes"] + 1 - n))
                    if not block:
                        break
                    n += len(block)
                    if n > item["bytes"]:
                        raise ValueError(f"Download exceeds declared size: {item['path']}")
                    h.update(block)
                    stream.write(block)
            if n != item["bytes"] or h.hexdigest() != item["sha256"]:
                raise ValueError(f"Download checksum/size mismatch: {item['path']}")
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)
    from .inputs import write_json
    write_json(root / "SOURCE.json", dataset)
    write_json(root / "DOWNLOAD_RECEIPT.json", verify_dataset(dataset_id, root, include_optional=include_optional))
    return root


def _raw_dir(data_dir):
    path = Path(data_dir)
    return path / "raw" if (path / "raw").is_dir() else path


def _cells_table(data_dir):
    import pandas as pd
    root = _raw_dir(data_dir)
    if (root / "cells.parquet").is_file():
        obs = pd.read_parquet(root / "cells.parquet")
    elif (root / "cells.csv.gz").is_file():
        obs = pd.read_csv(root / "cells.csv.gz", dtype={"cell_id": str})
    elif (root / "cells.csv").is_file():
        obs = pd.read_csv(root / "cells.csv", dtype={"cell_id": str})
    else:
        raise ValueError("Provide cells.parquet, cells.csv.gz or cells.csv next to cell_feature_matrix.h5.")
    required = {"cell_id", "x_centroid", "y_centroid"}
    if not required <= set(obs) or obs[list(required)].isna().any().any():
        raise ValueError("Cell table requires nonmissing cell_id, x_centroid and y_centroid.")
    obs["cell_id"] = obs["cell_id"].astype(str)
    if obs.cell_id.duplicated().any() or obs.cell_id.str.strip().eq("").any():
        raise ValueError("Cell table IDs must be nonempty and unique.")
    import numpy as np
    if not np.isfinite(obs[["x_centroid", "y_centroid"]].to_numpy(dtype=float)).all():
        raise ValueError("All Xenium coordinates must be finite micrometres.")
    return obs.set_index("cell_id")


def _decode(values):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in values]


def _selected_obs(data_dir, *, max_cells=1000, center=None, radius_um=None):
    import h5py
    import numpy as np
    if type(max_cells) is not int or max_cells <= 0:
        raise ValueError("max_cells must be a positive integer.")
    obs = _cells_table(data_dir)
    with h5py.File(_raw_dir(data_dir) / "cell_feature_matrix.h5", "r") as f:
        barcodes = _decode(f["matrix/barcodes"][:])
    if len(set(barcodes)) != len(barcodes) or set(barcodes) != set(obs.index):
        raise ValueError("Matrix and cell table must contain exactly the same unique cell IDs.")
    obs = obs.loc[barcodes].copy()
    obs["matrix_column"] = np.arange(len(barcodes))
    if center is None:
        if radius_um is not None:
            raise ValueError("radius_um requires an explicit center in micrometres.")
        if len(obs) > max_cells:
            raise ValueError(f"{len(obs)} cells exceeds max_cells={max_cells}; select an explicit center/ROI first.")
    else:
        if len(center) != 2 or any(not math.isfinite(float(x)) for x in center):
            raise ValueError("center must be two finite coordinates in micrometres.")
        distance = ((obs.x_centroid - float(center[0]))**2 + (obs.y_centroid - float(center[1]))**2)
        if radius_um is not None:
            if not math.isfinite(float(radius_um)) or radius_um <= 0:
                raise ValueError("radius_um must be finite and positive.")
            obs = obs.loc[distance <= radius_um**2]
        # Stable tie-break by ID; every ROI selection is reproducible.
        obs = obs.assign(_distance=distance.loc[obs.index]).sort_values(["_distance", "matrix_column"], kind="stable")
        obs = obs.iloc[:max_cells].drop(columns="_distance")
        if obs.empty:
            raise ValueError("The selected region contains no cells.")
    return obs, set(barcodes)


def write_annotation_template(data_dir, out, *, max_cells=1000, center=None, radius_um=None):
    """Write exact selected cell IDs, optional vendor clusters, and blank reviewed labels.

    A graph-clustering file may be supplied in analysis.tar.gz. It is read in
    place, never extracted. Blank cell_type/lineage entries intentionally fail
    load_xenium_cells until a reviewer fills them in.
    """
    import pandas as pd
    obs, _ = _selected_obs(data_dir, max_cells=max_cells, center=center, radius_um=radius_um)
    table = pd.DataFrame({"cell_id": obs.index, "cluster_id": "", "cell_type": "", "lineage": ""})
    archive = _raw_dir(data_dir) / "analysis.tar.gz"
    member = "analysis/clustering/gene_expression_graphclust/clusters.csv"
    if archive.is_file():
        with tarfile.open(archive, "r:gz") as tf:
            info = tf.getmember(member)
            if not info.isfile() or info.size > 32 * 2**20:
                raise ValueError("Unexpected clustering member type or size.")
            with tf.extractfile(info) as stream:
                clusters = pd.read_csv(stream, dtype={"Barcode": str})
        if set(clusters.columns) != {"Barcode", "Cluster"} or clusters.Barcode.duplicated().any():
            raise ValueError("Vendor graph clustering has an unexpected schema or duplicate IDs.")
        cluster_map = clusters.set_index("Barcode").Cluster
        # Vendor analysis may omit low-quality cells present in the raw matrix.
        # Missing cluster IDs stay blank; measured cells are never silently lost.
        table["cluster_id"] = table.cell_id.map(cluster_map).fillna("")
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    return path


def load_xenium_cells(data_dir, annotations_path, max_cells=1000, center=None, radius_um=None):
    """Return (validated cell records, assayed Gene Expression symbols).

    annotations_path is a CSV with cell_id, cell_type, lineage. It may cover the
    whole sample or the explicitly selected ROI, but unknown IDs, duplicate IDs,
    missing selected labels and blank reviewed labels are errors. Coordinates
    are measured centroids in micrometres. Cell boundaries are not inferred.
    """
    import h5py
    import numpy as np
    import pandas as pd
    from .inputs import validate_cells
    obs, all_ids = _selected_obs(data_dir, max_cells=max_cells, center=center, radius_um=radius_um)
    annotations = pd.read_csv(annotations_path, dtype=str, keep_default_na=False)
    required = {"cell_id", "cell_type", "lineage"}
    if not required <= set(annotations):
        raise ValueError("Annotation CSV requires cell_id, cell_type and lineage.")
    if annotations.cell_id.duplicated().any() or not set(annotations.cell_id) <= all_ids:
        raise ValueError("Annotation IDs must be unique and present in the measured matrix.")
    annotations = annotations.set_index("cell_id")
    if not set(obs.index) <= set(annotations.index):
        raise ValueError("Every selected cell needs reviewed annotation; rows are missing.")
    annotations = annotations.loc[obs.index]
    if annotations[["cell_type", "lineage"]].apply(lambda s: s.str.strip().eq("")).any().any():
        raise ValueError("Fill every selected cell_type and lineage; vendor cluster numbers are not annotations.")
    records = []
    with h5py.File(_raw_dir(data_dir) / "cell_feature_matrix.h5", "r") as file:
        matrix = file["matrix"]
        names = _decode(matrix["features/name"][:])
        kinds = _decode(matrix["features/feature_type"][:])
        gene_indices = {i: name for i, (name, kind) in enumerate(zip(names, kinds)) if kind == "Gene Expression"}
        genes = list(gene_indices.values())
        if not genes or len(set(genes)) != len(genes) or any(not g.strip() for g in genes):
            raise ValueError("Gene Expression feature symbols must be nonempty and unique.")
        shape = tuple(int(x) for x in matrix["shape"][:])
        indptr = matrix["indptr"][:]
        if shape != (len(names), len(all_ids)) or len(kinds) != len(names) or len(indptr) != len(all_ids) + 1:
            raise ValueError("Expected standard 10x genes-by-cells CSC shape and pointers.")
        if len(matrix["data"]) != len(matrix["indices"]) or indptr[0] != 0 or indptr[-1] != len(matrix["data"]) or (np.diff(indptr.astype('int64')) < 0).any():
            raise ValueError("Malformed 10x sparse matrix pointers.")
        for cid, row in obs.iterrows():
            column = int(row.matrix_column)
            start, end = int(indptr[column]), int(indptr[column+1])
            indices, values = matrix["indices"][start:end], matrix["data"][start:end]
            if len(set(indices)) != len(indices) or any(int(i) < 0 or int(i) >= len(names) for i in indices):
                raise ValueError("Malformed sparse column indices.")
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError("Count values must be finite and nonnegative.")
            expression = {gene_indices[int(i)]: float(v) for i, v in zip(indices, values) if int(i) in gene_indices and v != 0}
            label = annotations.loc[cid]
            records.append(dict(cell_id=str(cid), cell_type=label.cell_type.strip(), lineage=label.lineage.strip(),
                                x=float(row.x_centroid), y=float(row.y_centroid), expression=expression,
                                sample_id=Path(data_dir).name))
    return validate_cells(records), genes


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list")
    for name in ("download", "verify"):
        p = sub.add_parser(name)
        p.add_argument("dataset_id"); p.add_argument("data_dir")
        p.add_argument("--include-optional", action="store_true")
        if name == "download": p.add_argument("--max-mib", type=int, default=100)
    args = parser.parse_args(argv)
    if args.command == "list":
        print(json.dumps(list_datasets(), indent=2))
    elif args.command == "verify":
        report = verify_dataset(args.dataset_id, args.data_dir, include_optional=args.include_optional)
        print(json.dumps(report, indent=2))
        return 0 if report["ok"] else 1
    else:
        path = download_dataset(args.dataset_id, args.data_dir, include_optional=args.include_optional,
                                max_bytes=args.max_mib * 2**20)
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
