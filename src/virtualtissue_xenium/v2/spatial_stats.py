"""Structural summaries with an explicit observation window and held-out ROIs.

Only rectangular windows in micrometres are supported. We never infer a window
from the observed cell bounding box or silently clip cells that move outside it.
Cross-K uses translation edge correction; L=sqrt(K/pi). These exploratory
stationary-process statistics do not establish biological validity by themselves.
"""
from __future__ import annotations

from collections import Counter
import math

from .ensemble import quantile


def _window(window):
    if not isinstance(window, (list, tuple)) or len(window) != 4:
        raise ValueError("window must be [xmin, ymin, xmax, ymax] in micrometres.")
    if any(type(v) not in (int, float) or not math.isfinite(v) for v in window):
        raise ValueError("Window bounds must be finite numbers.")
    xmin, ymin, xmax, ymax = window
    if xmax <= xmin or ymax <= ymin:
        raise ValueError("Window must have positive area.")
    return tuple(window)


def cross_l(cells, radii_um, *, window, sample_id, source_lineage="T_CD8", target_lineage="malignant"):
    """Translation-corrected cross-K/L for distinct lineages in one sample.

K = area²/(n_source*n_target) * sum(1{distance<=r}/overlap_area).
The observation window must contain all selected living cells, and radii are
restricted to half the shorter side to avoid unstable long-distance estimates.
"""
    xmin, ymin, xmax, ymax = _window(window)
    radii = list(radii_um)
    if not radii or any(type(r) not in (int, float) or not math.isfinite(r) or r < 0 for r in radii) or radii != sorted(set(radii)):
        raise ValueError("Radii must be unique increasing nonnegative numbers.")
    width, height = xmax-xmin, ymax-ymin
    if max(radii) > min(width, height)/2:
        raise ValueError("Maximum radius must not exceed half the shorter window side.")
    if source_lineage == target_lineage:
        raise ValueError("This cross-type estimator requires distinct source and target lineages.")
    if not isinstance(sample_id, str) or not sample_id:
        raise ValueError("An explicit sample_id is required.")
    selected = [c for c in (cells.values() if isinstance(cells, dict) else cells)
                if c.get("status", "alive") == "alive" and c["sample_id"] == sample_id]
    for cell in selected:
        if any(type(cell.get(k)) not in (int, float) or not math.isfinite(cell[k]) for k in ("x", "y")):
            raise ValueError("Spatial statistics require finite measured centroids.")
        if not (xmin <= cell["x"] <= xmax and ymin <= cell["y"] <= ymax):
            raise ValueError("Living cell lies outside the declared observation window; do not silently change ROI boundaries.")
    source = [c for c in selected if c["lineage"] == source_lineage]
    target = [c for c in selected if c["lineage"] == target_lineage]
    area = width*height
    pairs = []
    for a in source:
        for b in target:
            dx, dy = abs(a["x"]-b["x"]), abs(a["y"]-b["y"])
            distance = math.hypot(dx, dy)
            if distance <= max(radii):
                overlap = (width-dx)*(height-dy)
                if overlap <= 0:
                    raise ValueError("Undefined translation edge weight at the window boundary.")
                pairs.append((distance, 1/overlap))
    curve = []
    for radius in radii:
        k = area**2/(len(source)*len(target))*sum(w for d, w in pairs if d <= radius) if source and target else None
        l = math.sqrt(k/math.pi) if k is not None else None
        curve.append(dict(radius_um=radius, K_um2=k, L_um=l, L_minus_r_um=l-radius if l is not None else None))
    return dict(sample_id=sample_id, window_um=list(window), source_lineage=source_lineage,
        target_lineage=target_lineage, n_source=len(source), n_target=len(target), curve=curve,
        status="ok" if source and target else "missing_lineage", edge_correction="translation_rectangle",
        assumption="stationary marked spatial point process within the declared rectangular observation window",
        references=["https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kcross.html",
                    "https://search.r-project.org/CRAN/refmans/spatstat.explore/html/edge.Trans.html"])


def cluster_sizes(cells, *, radius_um, sample_id, lineage="malignant"):
    """Connected-component sizes of a declared same-lineage radius graph."""
    if type(radius_um) not in (int, float) or not math.isfinite(radius_um) or radius_um <= 0:
        raise ValueError("Cluster radius must be finite and positive.")
    selected = [c for c in (cells.values() if isinstance(cells, dict) else cells)
                if c.get("status", "alive") == "alive" and c["sample_id"] == sample_id and c["lineage"] == lineage]
    if any(any(type(c.get(k)) not in (int, float) or not math.isfinite(c[k]) for k in ("x", "y")) for c in selected):
        raise ValueError("Clusters require finite spatial centroids.")
    parents = list(range(len(selected)))
    def find(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i
    for i, a in enumerate(selected):
        for j in range(i):
            b = selected[j]
            if math.hypot(a["x"]-b["x"], a["y"]-b["y"]) <= radius_um:
                parents[find(i)] = find(j)
    sizes = sorted(Counter(find(i) for i in range(len(selected))).values())
    return dict(sample_id=sample_id, lineage=lineage, radius_um=radius_um, sizes=sizes,
                n_cells=len(selected), n_clusters=len(sizes),
                definition="connected components, including singleton cells; graph is uncapped and does not cross samples")


def held_out_envelope(observed, held_out):
    """Pointwise empirical held-out ROI band, not a global significance test."""
    if not isinstance(held_out, list) or len(held_out) < 2:
        raise ValueError("At least two independent held-out ROI summaries are required.")
    if observed.get("status") != "ok" or any(r.get("status") != "ok" for r in held_out):
        raise ValueError("Undefined spatial curves cannot establish a reference envelope.")
    for ref in held_out:
        for key in ("source_lineage", "target_lineage", "edge_correction"):
            if observed[key] != ref[key]:
                raise ValueError(f"Spatial summaries differ in {key}.")
        if [r["radius_um"] for r in observed["curve"]] != [r["radius_um"] for r in ref["curve"]]:
            raise ValueError("All spatial curves need identical declared radii.")
    rows = []
    for index, point in enumerate(observed["curve"]):
        values = [ref["curve"][index]["L_minus_r_um"] for ref in held_out]
        low, high = quantile(values, .025), quantile(values, .975)
        rows.append(dict(radius_um=point["radius_um"], observed_L_minus_r_um=point["L_minus_r_um"],
                         reference_q025=low, reference_q975=high,
                         inside_pointwise_envelope=low <= point["L_minus_r_um"] <= high))
    return dict(n_reference_rois=len(held_out), curve=rows,
        interpretation="pointwise empirical held-out ROI range; not simultaneous coverage, a p-value, or clinical validation",
        warning="Few reference ROIs make tail quantiles unstable; preregister sampling and report all ROIs.")
