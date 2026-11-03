"""Deterministic ROI proposals and reviewer agreement on investigator labels.

These tools never assign a lineage, establish malignancy, or claim that supplied
labels have been independently reviewed. Positions are measured micrometres.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import heapq
import json
import math
from pathlib import Path

from .ontology import lineage_rules

IMMUNE_LINEAGES = frozenset({"T_CD8", "T_CD4", "T_reg", "B", "NK", "myeloid", "microglia",
                            "plasma", "mast", "neutrophil", "dendritic"})


def _positive(value, field):
    if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{field} must be a positive finite number.")
    return float(value)


def _checked_cells(cells):
    if not isinstance(cells, list) or not cells:
        raise ValueError("ROI selection needs a nonempty list of already annotated cells.")
    seen, grouped, rules = set(), defaultdict(list), lineage_rules()
    for cell in cells:
        if not isinstance(cell, dict):
            raise ValueError("Each annotated cell must be an object.")
        for key in ("cell_id", "sample_id", "cell_type", "lineage"):
            if not isinstance(cell.get(key), str) or not cell[key].strip():
                raise ValueError(f"Each annotated cell requires a nonempty {key}.")
        if cell["cell_id"] in seen:
            raise ValueError("cell_id must be globally unique; qualify IDs by sample before joining.")
        seen.add(cell["cell_id"])
        if cell["lineage"] not in rules:
            raise ValueError(f"Unsupported lineage: {cell['lineage']!r}.")
        for coordinate in ("x", "y"):
            value = cell.get(coordinate)
            if type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError("ROI selection requires finite measured x/y coordinates in micrometres.")
        grouped[cell["sample_id"]].append(cell)
    for sample, members in grouped.items():
        lineages = {c["lineage"] for c in members}
        if not {"malignant", "T_CD8"} <= lineages:
            raise ValueError(f"Sample {sample!r} needs reviewed malignant and T_CD8 labels; unknown labels cannot identify the melanoma ROIs.")
    return grouped


def _contact_count(cells, radius):
    """Count distinct CD8–malignant pairs at centroid distance <= radius."""
    grid = defaultdict(list)
    for c in cells:
        if c["lineage"] == "malignant":
            grid[(math.floor(c["x"] / radius), math.floor(c["y"] / radius))].append(c)
    total, r2 = 0, radius * radius
    for c in cells:
        if c["lineage"] != "T_CD8":
            continue
        bx, by = math.floor(c["x"] / radius), math.floor(c["y"] / radius)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                total += sum((c["x"] - m["x"]) ** 2 + (c["y"] - m["y"]) ** 2 <= r2
                             for m in grid[(bx + dx, by + dy)])
    return total


def select_rois(cells, *, cells_per_roi=1000, grid_um=500.0, max_candidates=50,
                contact_radius_um=30.0, core_max_immune_fraction=.05,
                core_min_malignant_fraction=.5, require_disjoint=True):
    """Propose invasive-front/core ROIs from supplied labels, per sample.

    Candidate centres are occupied square-grid nodes (origin 0,0), rounded
    deterministically half upwards. If more than max_candidates exist, evenly
    spaced indices of the lexicographic node list are retained. Each ROI is the
    nearest cells_per_roi cells; ties use cell_id. A missing eligible core/front
    is represented as null with a reason, never silently relaxed.
    """
    if type(cells_per_roi) is not int or cells_per_roi < 1:
        raise ValueError("cells_per_roi must be a positive integer.")
    if type(max_candidates) is not int or not 1 <= max_candidates <= 10000:
        raise ValueError("max_candidates must be an integer from 1 to 10000.")
    if type(require_disjoint) is not bool:
        raise ValueError("require_disjoint must be boolean.")
    grid_um = _positive(grid_um, "grid_um")
    contact_radius_um = _positive(contact_radius_um, "contact_radius_um")
    for name, value in (("core_max_immune_fraction", core_max_immune_fraction),
                        ("core_min_malignant_fraction", core_min_malignant_fraction)):
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{name} must be between zero and one.")
    grouped = _checked_cells(cells)
    samples = {}
    for sample, members in sorted(grouped.items()):
        if len(members) < cells_per_roi:
            raise ValueError(f"Sample {sample!r} contains fewer than {cells_per_roi} cells; explicitly choose a smaller ROI.")
        nodes = sorted({(math.floor(c["x"] / grid_um + .5) * grid_um,
                         math.floor(c["y"] / grid_um + .5) * grid_um) for c in members})
        if len(nodes) > max_candidates:
            indices = [0] if max_candidates == 1 else [i * (len(nodes) - 1) // (max_candidates - 1) for i in range(max_candidates)]
            nodes = [nodes[i] for i in indices]
        candidates = []
        for x, y in nodes:
            selected = heapq.nsmallest(cells_per_roi, members,
                key=lambda c: ((c["x"] - x) ** 2 + (c["y"] - y) ** 2, c["cell_id"]))
            counts = Counter(c["lineage"] for c in selected)
            candidates.append(dict(center_um=[x, y], cell_ids=sorted(c["cell_id"] for c in selected),
                n_cells=len(selected), cd8_malignant_contacts=_contact_count(selected, contact_radius_um),
                immune_fraction=sum(n for lineage, n in counts.items() if lineage in IMMUNE_LINEAGES) / len(selected),
                malignant_fraction=counts["malignant"] / len(selected), unknown_fraction=counts["unknown"] / len(selected),
                lineage_counts=dict(sorted(counts.items())),
                radius_to_farthest_cell_um=max(math.hypot(c["x"] - x, c["y"] - y) for c in selected)))
        ranked = sorted(candidates, key=lambda c: (-c["cd8_malignant_contacts"], c["center_um"]))
        front = ranked[0] if ranked and ranked[0]["cd8_malignant_contacts"] > 0 else None
        front_ids = set(front["cell_ids"]) if front else set()
        core_candidates = [c for c in candidates if c["immune_fraction"] < core_max_immune_fraction
                           and c["malignant_fraction"] >= core_min_malignant_fraction
                           and (not require_disjoint or not front_ids.intersection(c["cell_ids"]))]
        core = min(core_candidates, key=lambda c: (-c["malignant_fraction"], c["immune_fraction"], c["center_um"])) if core_candidates else None
        samples[sample] = dict(invasive_front=front, tumour_core=core, candidates=candidates,
            complete=front is not None and core is not None,
            missing_reasons=(["no_CD8_malignant_contacts_at_declared_radius"] if front is None else []) +
                            (["no_core_satisfies_declared_fractions_and_overlap_rule"] if core is None else []))
    return dict(schema_version="3.0", purpose="ROI proposals from investigator-supplied annotations",
        annotation_review="Not performed by this tool; malignant identity and independent review are investigator responsibilities.",
        parameters=dict(cells_per_roi=cells_per_roi, grid_um=grid_um, grid_origin_um=[0., 0.],
            max_candidates=max_candidates, contact_radius_um=contact_radius_um,
            core_max_immune_fraction=core_max_immune_fraction,
            core_min_malignant_fraction=core_min_malignant_fraction, require_disjoint=require_disjoint,
            candidate_rule="occupied_grid_nodes_evenly_spaced_lexicographic_indices",
            core_rule="immune_fraction strictly below threshold; malignant_fraction at least threshold"),
        input_sha256=hashlib.sha256(json.dumps(sorted(cells, key=lambda c: c["cell_id"]), sort_keys=True,
                                              separators=(",", ":"), allow_nan=False).encode()).hexdigest(),
        samples=samples, complete=all(result["complete"] for result in samples.values()))


def _read_reviewer(path, label_column):
    path = Path(path)
    rows = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not {"cell_id", label_column} <= set(reader.fieldnames or []):
            raise ValueError(f"Reviewer CSV must contain cell_id and {label_column}.")
        for row in reader:
            cid, label = (row["cell_id"] or "").strip(), (row[label_column] or "").strip()
            if not cid or not label:
                raise ValueError("Reviewer IDs and labels must be filled; use unknown explicitly for unresolved cells.")
            if cid in rows:
                raise ValueError("Duplicate cell_id in reviewer CSV.")
            rows[cid] = label
    if not rows:
        raise ValueError("Reviewer CSV contains no annotations.")
    return rows, dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def _agreement(first, second):
    n = len(first)
    if not n:
        return dict(n=0, observed_agreement=None, expected_agreement=None, kappa=None,
                    kappa_status="undefined_no_common_cells")
    observed = sum(a == b for a, b in zip(first, second)) / n
    a, b = Counter(first), Counter(second)
    expected = sum(a[label] * b[label] for label in a.keys() | b.keys()) / (n * n)
    denominator = 1 - expected
    return dict(n=n, observed_agreement=observed, expected_agreement=expected,
                kappa=(observed - expected) / denominator if denominator else None,
                kappa_status="defined" if denominator else "undefined_zero_expected_disagreement")


def compare_reviewers(first_csv, second_csv, *, label_column="lineage", id_scope="exact"):
    """Cohen kappa on exact shared IDs, or an explicitly requested intersection."""
    if id_scope not in ("exact", "intersection"):
        raise ValueError("id_scope must be exact or intersection.")
    first, first_source = _read_reviewer(first_csv, label_column)
    second, second_source = _read_reviewer(second_csv, label_column)
    only_first, only_second = sorted(first.keys() - second.keys()), sorted(second.keys() - first.keys())
    if id_scope == "exact" and (only_first or only_second):
        raise ValueError("Reviewer cell IDs differ; resolve the join or explicitly choose id_scope='intersection'.")
    ids = sorted(first.keys() & second.keys())
    a, b = [first[cid] for cid in ids], [second[cid] for cid in ids]
    labels = sorted(set(a) | set(b))
    confusion = {label: dict(Counter(y for x, y in zip(a, b) if x == label)) for label in labels}
    return dict(schema_version="3.0", label_column=label_column, id_scope=id_scope,
        sources=[first_source, second_source], common_cell_ids=ids, excluded_only_first=only_first,
        excluded_only_second=only_second, overall=_agreement(a, b), confusion=confusion,
        by_lineage={label: _agreement([x == label for x in a], [x == label for x in b]) for label in labels},
        disagreements=[dict(cell_id=cid, reviewer_1=first[cid], reviewer_2=second[cid])
                       for cid in ids if first[cid] != second[cid]],
        interpretation="Agreement between supplied labels; not proof of biological truth or independent review. No simulated cells are used.")


def write_annotation_comparison(first_csv, second_csv, output, **kwargs):
    result = compare_reviewers(first_csv, second_csv, **kwargs)
    output = Path(output)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Annotation comparison output must be a new or empty directory.")
    output.mkdir(parents=True, exist_ok=True)
    (output / "agreement.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    with (output / "disagreements.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cell_id", "reviewer_1", "reviewer_2"])
        writer.writeheader()
        writer.writerows(result["disagreements"])
    return result
