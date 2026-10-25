"""t=0 fidelity of the simulated geometry against the measured Xenium polygons."""
from collections import defaultdict

import numpy as np

from .geometry import power_geometry

TUMOUR = {"Tumor", "Tumor_proliferating"}


def _quantiles(values):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return None
    q = np.quantile(values, [.25, .5, .75])
    return dict(n=int(len(values)), q25=round(float(q[0]), 4), median=round(float(q[1]), 4), q75=round(float(q[2]), 4))


def _f1(measured, simulated):
    tp = len(measured & simulated)
    precision = tp / len(simulated) if simulated else None
    recall = tp / len(measured) if measured else None
    f1 = 2 * tp / (len(measured) + len(simulated)) if measured or simulated else None
    return dict(measured=len(measured), simulated=len(simulated), shared=tp,
                precision=None if precision is None else round(precision, 4),
                recall=None if recall is None else round(recall, 4), f1=None if f1 is None else round(f1, 4))


def fidelity_report(atlas, tissue, geometry=None, idx=None, label="fitted"):
    """IoU per type and contact-graph F1 at t=0 (pairs involving at least one ROI cell)."""
    if geometry is None:
        idx, geometry = tissue.geometry()
    local = {int(g): k for k, g in enumerate(idx)}
    pos = {int(cid): g for g, cid in enumerate(tissue.ids.tolist())}
    iou = defaultdict(list)
    for g, k in local.items():
        if tissue.origins[g] != "roi":
            continue
        measured = tissue.polygons[g]
        simulated = geometry.polygon(k)
        union = measured.union(simulated).area
        value = measured.intersection(simulated).area / union if union > 0 else 0.
        iou[tissue.labels[g]].append(value)
        iou["all_roi"].append(value)
    measured_pairs = set()
    for cell in atlas["cells"]:
        if cell["role"] != "roi":
            continue
        for edge in cell["contacts"]:
            if edge["gap_um"] <= 1 and edge["cell_id"] in pos:
                measured_pairs.add(tuple(sorted((cell["cell_id"], edge["cell_id"]))))
    i, j, _ = geometry.contacts(tissue.g["contact_max_gap_um"])
    simulated_pairs = set()
    for a, b in zip(idx[i].tolist(), idx[j].tolist()):
        if tissue.origins[a] == "roi" or tissue.origins[b] == "roi":
            simulated_pairs.add(tuple(sorted((int(tissue.ids[a]), int(tissue.ids[b])))))
    lab = {int(cid): tissue.labels[g] for g, cid in enumerate(tissue.ids.tolist())}

    def cd8_tumour(pairs):
        return {p for p in pairs if {lab[p[0]], lab[p[1]]} & {"T_CD8_cytotoxic"} and {lab[p[0]], lab[p[1]]} & TUMOUR}

    def roi_only(pairs):
        return {p for p in pairs if tissue.origins[pos[p[0]]] == "roi" and tissue.origins[pos[p[1]]] == "roi"}

    area_ratio = geometry.area[[local[g] for g in local if tissue.origins[g] == "roi"]] / \
        tissue.measured_area[[g for g in local if tissue.origins[g] == "roi"]]
    return dict(variant=label,
                iou_by_type={k: _quantiles(v) for k, v in sorted(iou.items())},
                contacts_any_roi=_f1(measured_pairs, simulated_pairs),
                contacts_roi_roi=_f1(roi_only(measured_pairs), roi_only(simulated_pairs)),
                contacts_cd8_tumour=_f1(cd8_tumour(measured_pairs), cd8_tumour(simulated_pairs)),
                area_ratio_simulated_over_measured=_quantiles(area_ratio),
                definitions=dict(measured_contact="polygon gap <= 1 um (H04)",
                                 simulated_contact=f"power-adjacent and gap <= {tissue.g['contact_max_gap_um']} um (shared interface = gap 0)",
                                 iou="area(measured ∩ simulated) / area(measured ∪ simulated), ROI cells"))


def unfitted_report(atlas, tissue):
    """Same metrics with equal-area radii (kappa = 1/pi): shows what the fit contributes."""
    idx = np.flatnonzero(tissue.present)
    geometry = power_geometry(tissue.xy[idx], np.sqrt(tissue.measured_area[idx] / np.pi), tissue.g["dummy_sites"])
    return fidelity_report(atlas, tissue, geometry, idx, label="equal_area_radius_unfitted")
