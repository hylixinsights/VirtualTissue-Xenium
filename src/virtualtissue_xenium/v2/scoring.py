"""Depth-conditioned, within-type signature evidence (lane B of the v3 plan).

Method ``depth_conditioned_auc``:

1. Per cell, an AUCell-style recovery AUC over the assayed panel: genes are
   ranked by count (descending; ties share their average rank; undetected genes
   are unranked), and the AUC of the signature's recovery curve within the top
   ``max_rank`` positions is normalised by its maximum. Deterministic: no random
   tie breaking.
2. Per sample, cells are binned by total counts into ``depth_bins`` quantile bins.
3. Each cell's AUC becomes a mid-rank percentile within its reference group:
   (cell_type x depth bin) when that group has at least ``min_reference_cells``
   scored cells, otherwise its cell_type (flagged ``cell_type_no_depth``),
   otherwise no score (``insufficient_reference``; multiplier exactly 1).
4. z = probit(percentile), clipped to [-3, 3].
5. The action multiplier is m = exp(beta * z) / mean_ref(exp(beta * z)): it
   redistributes the hazard among cells of the same reference group and
   preserves the group's mean hazard, so the kinetics table keeps its meaning as
   a population-average rate. beta = 0 is the label-only model.

This is a relative within-type score, not AUCell as published (which ranks the
whole transcriptome with random tie breaks), not protein activity and not a fate
probability. Panels, depth and composition still limit what it can resolve.
"""
from __future__ import annotations

from collections import defaultdict
import math
from statistics import NormalDist

SCORING_METHODS = ("rank_legacy", "depth_conditioned_auc")
SCORING_DEFAULTS = dict(method="depth_conditioned_auc", auc_top_fraction=0.05, auc_min_rank=50,
                        depth_bins=10, min_reference_cells=30, min_overlap=5, beta=1.0, min_transcripts=0)
Z_CLIP = 3.0
_NORMAL = NormalDist()


def max_rank_for(panel_size, scoring):
    """Top-of-ranking window: the larger of a panel fraction and a floor, at most the panel."""
    return max(1, min(panel_size, max(math.ceil(scoring["auc_top_fraction"] * panel_size), scoring["auc_min_rank"])))


def recovery_auc(expression, panel, genes, max_rank):
    """Normalised recovery AUC of ``genes`` within the top ``max_rank`` panel genes of one cell."""
    by_value = defaultdict(list)
    for gene, value in expression.items():
        if value > 0 and gene in panel:
            by_value[value].append(gene)
    ranks, offset = {}, 0
    for value in sorted(by_value, reverse=True):
        tied = by_value[value]
        rank = offset + (len(tied) + 1) / 2
        for gene in tied:
            ranks[gene] = rank
        offset += len(tied)
    area = sum(max(0.0, max_rank - ranks[g] + 1) for g in genes if g in ranks)
    best = sum(max_rank - i + 1 for i in range(1, min(len(genes), max_rank) + 1))
    return min(1.0, area / best) if best > 0 else 0.0


def depth_bins(cells, n_bins):
    """Quantile bin of total counts within each sample; ties broken by cell_id (deterministic)."""
    by_sample = defaultdict(list)
    for cell in cells:
        by_sample[cell.get("sample_id", "sample")].append(cell)
    bins = {}
    for members in by_sample.values():
        ordered = sorted(members, key=lambda c: (sum(c["expression"].values()), c["cell_id"]))
        n = len(ordered)
        for i, cell in enumerate(ordered):
            bins[cell["cell_id"]] = min(n_bins - 1, i * n_bins // n)
    return bins


def _percentile(value, reference):
    less = sum(v < value for v in reference)
    equal = sum(v == value for v in reference)
    return (less + 0.5 * equal) / len(reference)


def beta_for(scoring, action):
    beta = scoring["beta"]
    return float(beta.get(action, 0.0)) if isinstance(beta, dict) else float(beta)


def score_population(cells, panel, signatures, scoring):
    """Return {cell_id: {signature: evidence}} and per-signature reference z lists.

    Cells listed in ``scoring['_qc_fail']`` are neither scored nor used as references.
    """
    max_rank = max_rank_for(len(panel), scoring)
    qc_fail = scoring.get("_qc_fail", set())
    scored = [c for c in cells if c["cell_id"] not in qc_fail]
    bins = depth_bins(scored, scoring["depth_bins"])
    evidence = {c["cell_id"]: {} for c in cells}
    references = {}
    for name, genes in signatures.items():
        overlap = sorted(set(genes) & panel)
        enough_panel = len(overlap) >= scoring["min_overlap"]
        auc = {c["cell_id"]: recovery_auc(c["expression"], panel, overlap, max_rank)
               for c in scored} if enough_panel else {}
        groups = defaultdict(list)
        for c in scored:
            if enough_panel:
                groups[(c["cell_type"], bins[c["cell_id"]])].append(auc[c["cell_id"]])
                groups[(c["cell_type"], None)].append(auc[c["cell_id"]])
        z_by_ref = defaultdict(list)
        for c in cells:
            cid = c["cell_id"]
            detected = [g for g in overlap if c["expression"].get(g, 0) > 0]
            record = dict(score=None, auc=None, percentile=None, z=None, reference=None, reference_size=0,
                          depth_bin=bins.get(cid), assayed_genes=overlap, detected_genes=detected,
                          signature_size=len(genes), panel_coverage=len(overlap) / len(genes),
                          detected_fraction_of_assayed=len(detected) / len(overlap) if overlap else None,
                          max_rank=max_rank)
            if cid in qc_fail:
                record["status"] = "qc_fail_below_min_transcripts"
            elif not enough_panel:
                record["status"] = "insufficient_panel_overlap"
            else:
                record["auc"] = auc[cid]
                key = (c["cell_type"], bins[cid])
                if len(groups[key]) < scoring["min_reference_cells"]:
                    key = (c["cell_type"], None)
                if len(groups[key]) < scoring["min_reference_cells"]:
                    record["status"] = "insufficient_reference"
                    record["reference_size"] = len(groups[key])
                else:
                    p = _percentile(auc[cid], groups[key])
                    z = max(-Z_CLIP, min(Z_CLIP, _NORMAL.inv_cdf(min(max(p, 1e-6), 1 - 1e-6))))
                    record.update(score=p, percentile=p, z=z, reference_size=len(groups[key]),
                                  reference="cell_type_x_depth_bin" if key[1] is not None else "cell_type_no_depth",
                                  status="available")
                    ref_id = f"{key[0]}|{'all' if key[1] is None else key[1]}"
                    record["reference_id"] = ref_id
                    z_by_ref[ref_id].append(z)
            evidence[cid][name] = record
        references[name] = z_by_ref
    return evidence, references


def expression_multiplier(record, references_for_signature, beta):
    """exp(beta z) normalised to mean 1 over the cell's reference group; 1 when unscored."""
    if record.get("z") is None or beta == 0:
        return 1.0, 1.0
    zs = references_for_signature[record["reference_id"]]
    normaliser = sum(math.exp(beta * z) for z in zs) / len(zs)
    return math.exp(beta * record["z"]) / normaliser, normaliser
