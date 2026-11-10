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

Opt-in ``tiered_panel_auc`` retains sample/type/depth references when supported,
uses a whole sample/type reference when depth strata are too small, then uses
explicitly weak within-cell panel evidence for rare types. This last tier has a
bounded multiplier without group normalization: it can change a rare type's
mean hazard and does not claim calibration from other reference cells.
"""
from __future__ import annotations

from collections import defaultdict
import math
import json
from statistics import NormalDist

SCORING_METHODS = ("rank_legacy", "depth_conditioned_auc", "tiered_panel_auc")
SCORING_DEFAULTS = dict(method="depth_conditioned_auc", auc_top_fraction=0.05, auc_min_rank=50,
                        depth_bins=10, min_reference_cells=30, min_overlap=5, beta=1.0, min_transcripts=0)
Z_CLIP = 3.0
TIERED_PANEL_VERSION = "tiered-panel-auc-1"
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
    if scoring.get("method") == "tiered_panel_auc":
        return _score_population_tiered(cells, panel, signatures, scoring)
    if "reference_scope" in scoring:
        return _score_population_v3(cells, panel, signatures, scoring)
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


def reference_partition(eligible, bins, scoring):
    """{cell_id: (sample, cell_type, depth bin or None)}: one non-overlapping reference per cell.

    With explicit cell_type fallback, the entire type falls back when any depth
    stratum is undersized. Shared by the internal AUC and imported measured
    scores, so both kinds of evidence are standardised against the same cells.
    """
    groups = defaultdict(list)
    for c in eligible:
        sample = c.get("sample_id", "sample") if scoring["reference_scope"] == "sample" else "__pooled__"
        groups[(sample, c["cell_type"])].append(c)
    partition = {}
    for group, members in groups.items():
        strata = defaultdict(list)
        for c in members:
            strata[bins[c["cell_id"]]].append(c)
        fallback = (scoring["reference_fallback"] == "cell_type" and
                    any(len(m) < scoring["min_reference_cells"] for m in strata.values()))
        for c in members:
            partition[c["cell_id"]] = (*group, None if fallback else bins[c["cell_id"]])
    return partition


def _score_population_v3(cells, panel, signatures, scoring):
    """Partition references by sample and type; never silently borrow another sample.

    With explicit cell_type fallback, the entire type falls back when any depth
    stratum is undersized. Percentiles and normalisation then use the same
    non-overlapping partition. The default is null for undersized strata.
    """
    max_rank = max_rank_for(len(panel), scoring)
    qc_fail = scoring.get("_qc_fail", set())
    eligible = [c for c in cells if c["cell_id"] not in qc_fail]
    bins = depth_bins(eligible, scoring["depth_bins"])
    partition = reference_partition(eligible, bins, scoring)
    evidence = {c["cell_id"]: {} for c in cells}
    references = {}
    for name, genes in signatures.items():
        overlap = sorted(set(genes) & panel)
        coverage = len(overlap) / len(genes)
        enough = len(overlap) >= scoring["min_overlap"]
        auc = {c["cell_id"]: recovery_auc(c["expression"], panel, overlap, max_rank)
               for c in eligible} if enough else {}
        reference_values = defaultdict(list)
        for cid, value in auc.items():
            reference_values[partition[cid]].append(value)
        zs = defaultdict(list)
        for c in cells:
            cid = c["cell_id"]
            detected = [g for g in overlap if c["expression"].get(g, 0) > 0]
            record = dict(score=None, auc=None, percentile=None, z=None, reference=None,
                          reference_size=0, depth_bin=bins.get(cid), assayed_genes=overlap,
                          detected_genes=detected, signature_size=len(genes), panel_coverage=coverage,
                          coverage_warning="less_than_half_signature_assayed" if coverage < .5 else None,
                          detected_fraction_of_assayed=len(detected)/len(overlap) if overlap else None,
                          max_rank=max_rank, informative=False)
            if cid in qc_fail:
                record["status"] = "qc_fail_below_min_transcripts"
            elif not enough:
                record["status"] = "insufficient_panel_overlap"
            else:
                key = partition[cid]
                values = reference_values[key]
                record.update(auc=auc[cid], reference_size=len(values))
                if len(values) < scoring["min_reference_cells"]:
                    record["status"] = "insufficient_reference"
                else:
                    percentile = _percentile(auc[cid], values)
                    z = max(-Z_CLIP, min(Z_CLIP, _NORMAL.inv_cdf(min(max(percentile, 1e-6), 1-1e-6))))
                    reference_id = json.dumps(key, ensure_ascii=False, separators=(",", ":"))
                    record.update(score=percentile, percentile=percentile, z=z, status="available",
                                  reference="cell_type_x_depth_bin" if key[-1] is not None else "cell_type_no_depth",
                                  reference_id=reference_id, reference_sample=key[0],
                                  informative=min(values) != max(values))
                    zs[reference_id].append(z)
            evidence[cid][name] = record
        references[name] = zs
    return evidence, references


def panel_null_moments(expression, panel, n_signature, max_rank):
    """Exact AUC moments for uniform size-k subsets of this cell's panel.

The finite population consists of observed rank-recovery weights (including
zero weights for undetected genes). This is an explicit random-gene-set model,
not a null distribution learned from biological reference cells. No normal
approximation, significance probability or biological activity is inferred.
"""
    panel = set(panel)
    size = len(panel)
    if type(n_signature) is not int or not 1 <= n_signature <= size:
        raise ValueError("Signature overlap must contain 1 through panel-size genes.")
    if type(max_rank) is not int or not 1 <= max_rank <= size:
        raise ValueError("max_rank must be an integer from 1 through panel size.")
    groups = defaultdict(list)
    for gene in panel:
        value = expression.get(gene, 0)
        if value > 0:
            groups[value].append(gene)
    weights, offset = [], 0
    for value in sorted(groups, reverse=True):
        count = len(groups[value])
        average_rank = offset + (count + 1)/2
        weights.extend([max(0.0, max_rank-average_rank+1)] * count)
        offset += count
    weights.extend([0.0] * (size-len(weights)))
    mean_weight = math.fsum(weights)/size
    variance_weight = math.fsum((weight-mean_weight)**2 for weight in weights)/size
    best = sum(max_rank-i+1 for i in range(1, min(n_signature, max_rank)+1))
    mean_auc = n_signature*mean_weight/best
    # Var(sum of k draws without replacement) = k(N-k)/(N-1) * population Var.
    variance_auc = (n_signature*(size-n_signature)/(size-1)*variance_weight/best**2) if size > 1 else 0.0
    return dict(mean_auc=mean_auc, sd_auc=math.sqrt(variance_auc),
                panel_size=size, signature_overlap=n_signature,
                assumption="uniform same-size gene subsets of this measured panel, conditional on this cell's ranks/ties/undetected genes; not biological calibration")


def _score_population_tiered(cells, panel, signatures, scoring):
    """Opt-in depth -> type -> cell-panel evidence, with no cross-type pooling.

Empirical groups form a non-overlapping partition: if any depth stratum for a
sample/type is too small, that whole type uses its sample-specific reference.
Rare types then use an explicitly weak panel-based heuristic. Unlike empirical
tiers, their multiplier does not preserve the group's mean hazard.
"""
    if scoring.get("reference_scope") != "sample":
        raise ValueError("tiered_panel_auc requires explicit sample-scoped v3 references.")
    bound = scoring.get("panel_max_multiplier", 2.0)
    if type(bound) not in (int, float) or not math.isfinite(bound) or not 1 <= bound <= 3:
        raise ValueError("panel_max_multiplier must lie in [1,3].")
    settings = dict(scoring, reference_fallback="cell_type")
    evidence, references = _score_population_v3(cells, panel, signatures, settings)
    cells_by_id = {c["cell_id"]: c for c in cells}
    for cid, signatures_for_cell in evidence.items():
        for record in signatures_for_cell.values():
            record.update(scoring_version=TIERED_PANEL_VERSION, calibration="empirical_within_sample_type",
                          multiplier_normalization="reference_group_mean_one", panel_max_multiplier=float(bound))
            if record["status"] != "insufficient_reference":
                if record["status"] != "available":
                    record.update(calibration="unavailable", multiplier_normalization="neutral")
                continue
            # Do not present a gene-set reference as n biological reference cells.
            record.update(candidate_cell_reference_size=record["reference_size"], reference_size=0,
                          reference="panel_gene_set_null", reference_sample=cells_by_id[cid].get("sample_id", "sample"),
                          calibration="weak_panel_calibration", multiplier_normalization="bounded_without_cell_normalization",
                          percentile=None, score=None, informative=False)
            moments = panel_null_moments(cells_by_id[cid]["expression"], panel,
                                         len(record["assayed_genes"]), record["max_rank"])
            record["panel_null"] = moments
            if not record["detected_genes"]:
                record["status"] = "no_detected_signature_genes_neutral"
            elif moments["sd_auc"] == 0:
                record["status"] = "constant_panel_null_neutral"
            else:
                raw_z = (record["auc"]-moments["mean_auc"])/moments["sd_auc"]
                z = max(-Z_CLIP, min(Z_CLIP, raw_z))
                record.update(status="available_panel_heuristic", score=record["auc"], z=z,
                              panel_null_z_raw=raw_z, informative=z != 0,
                              interpretation="Relative panel-rank evidence under a random-gene-set assumption; not measured activity, cell-reference calibration or a fate probability.")
    return evidence, references


def expression_multiplier(record, references_for_signature, beta):
    """Empirical mean-one effect or bounded panel heuristic; 1 when unscored."""
    if record.get("z") is None or beta == 0:
        return 1.0, 1.0
    if record.get("multiplier_normalization") == "bounded_without_cell_normalization":
        limit = math.log(record["panel_max_multiplier"])
        return math.exp(max(-limit, min(limit, beta*record["z"]))), 1.0
    zs = references_for_signature[record["reference_id"]]
    normaliser = sum(math.exp(beta * z) for z in zs) / len(zs)
    return math.exp(beta * record["z"]) / normaliser, normaliser
