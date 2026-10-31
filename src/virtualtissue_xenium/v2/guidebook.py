"""Expression evidence and explicitly assumed action preferences, without an LLM.

Scores use mean ranks among positively detected assayed genes, with zero for an
assayed but undetected gene. This is a small descriptive signature score, not
AUCell, protein activity, transition frequency, or a calibrated fate probability.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import hashlib
import json
import math

from .context import ACTIONS, LINEAGES, modifier_matches, validate_context

DEFAULT_GENE_SETS = {
    "cell_cycle": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "CDK1", "CCNB1", "UBE2C"],
    "cytotoxic": ["NKG7", "PRF1", "GZMB", "GNLY", "CTSW", "GZMH"],
    "activation": ["CD69", "IL2RA", "IFNG", "TNFRSF9", "FOS", "JUN"],
    "exhaustion_associated": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"],
    "apoptosis_associated": ["BAX", "BAK1", "FAS", "CASP3", "CASP8", "BBC3", "PMAIP1"],
    "inflammation": ["IL1B", "TNF", "IL6", "CXCL8", "CCL2", "NFKBIA"],
    "antigen_presentation": ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "B2M", "HLA-A"],
}
DESCRIPTIONS = {
    "REST": "Maintain the current simulated state for one step.",
    "PROLIFERATE": "Attempt one daughter cell; the runtime enforces the living population cap, without physical space or mechanics.",
    "APOPTOSIS": "Commit this cell to simulated apoptosis under an assumed action model.",
    "ACTIVATE": "Increase a simulated activation state; this does not measure protein activity.",
    "EXHAUST": "Increase a simulated exhaustion state under an explicit heuristic.",
    "RECOVER": "Decrease a simulated exhaustion state under an explicit heuristic.",
    "SECRETE": "Emit an abstract local simulated inflammatory signal.",
    "CLEAR_CORPSE": "Attempt removal of a neighboring simulated corpse.",
    "KILL": "Attempt killing of a recognized live malignant graph neighbor; local gates apply, and proximity does not establish contact.",
}
# These weights describe illustrative preferences per decision, not rates.
BASE_WEIGHTS = {
    "unknown": {"REST": 1.0},
    "T_CD8": {"REST": 1.0, "PROLIFERATE": .01, "ACTIVATE": .15, "EXHAUST": .02, "RECOVER": .02, "SECRETE": .08, "APOPTOSIS": .005},
    "T_CD4": {"REST": 1.0, "PROLIFERATE": .01, "ACTIVATE": .15, "EXHAUST": .02, "RECOVER": .02, "SECRETE": .08, "APOPTOSIS": .005},
    "B": {"REST": 1.0, "PROLIFERATE": .01, "ACTIVATE": .1, "SECRETE": .04, "APOPTOSIS": .005},
    "myeloid": {"REST": 1.0, "PROLIFERATE": .005, "ACTIVATE": .1, "SECRETE": .1, "CLEAR_CORPSE": .1, "APOPTOSIS": .005},
    "microglia": {"REST": 1.0, "PROLIFERATE": .005, "ACTIVATE": .1, "SECRETE": .1, "CLEAR_CORPSE": .1, "APOPTOSIS": .005},
    "fibroblast": {"REST": 1.0, "PROLIFERATE": .02, "SECRETE": .08, "APOPTOSIS": .005},
    "endothelial": {"REST": 1.0, "PROLIFERATE": .01, "SECRETE": .04, "APOPTOSIS": .005},
    "epithelial": {"REST": 1.0, "PROLIFERATE": .02, "SECRETE": .04, "APOPTOSIS": .005},
    "malignant": {"REST": 1.0, "PROLIFERATE": .08, "SECRETE": .08, "APOPTOSIS": .005},
    "neuron": {"REST": 1.0, "ACTIVATE": .04, "APOPTOSIS": .005},
    "astrocyte": {"REST": 1.0, "PROLIFERATE": .005, "ACTIVATE": .06, "SECRETE": .04, "APOPTOSIS": .005},
    "oligodendrocyte": {"REST": 1.0, "APOPTOSIS": .005},
}
ACTION_SIGNATURE = {
    "PROLIFERATE": "cell_cycle", "ACTIVATE": "activation", "EXHAUST": "exhaustion_associated",
    "APOPTOSIS": "apoptosis_associated", "SECRETE": "inflammation", "KILL": "cytotoxic",
}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _validate_cells(cells):
    if not isinstance(cells, list) or not cells:
        raise ValueError("cells must be a nonempty list of annotated cells.")
    checked, ids, type_lineages = [], set(), {}
    for cell in cells:
        if not isinstance(cell, dict):
            raise ValueError("Every cell must be an object.")
        for field in ("cell_id", "cell_type", "lineage"):
            if not isinstance(cell.get(field), str) or not cell[field].strip():
                raise ValueError(f"Every cell needs a nonempty {field}.")
        if cell["cell_id"] in ids:
            raise ValueError("cell_id values must be unique.")
        ids.add(cell["cell_id"])
        if cell["lineage"] not in LINEAGES:
            raise ValueError(f"Unsupported lineage {cell['lineage']!r}; provide an explicit supported mapping.")
        previous = type_lineages.setdefault(cell["cell_type"], cell["lineage"])
        if previous != cell["lineage"]:
            raise ValueError("A cell_type cannot map to multiple lineages; refine the annotation or mapping.")
        expression = cell.get("expression")
        if not isinstance(expression, dict):
            raise ValueError("Every cell needs an expression gene-to-nonnegative-count map.")
        normalized = {}
        for gene, value in expression.items():
            if not isinstance(gene, str) or not gene.strip():
                raise ValueError("Expression gene names must be nonempty strings.")
            symbol = gene.strip().upper()
            if symbol in normalized:
                raise ValueError("Expression genes must be unique after case normalization.")
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError("Expression values must be finite nonnegative numbers, not booleans.")
            normalized[symbol] = float(value)
        checked.append(dict(cell, expression=normalized))
    return checked


def _rank_evidence(expression, panel, gene_sets):
    # Average ranks for ties; no arbitrary gene-name tie breaking.
    values = defaultdict(list)
    for gene, value in expression.items():
        if value > 0:
            values[value].append(gene)
    ranks, offset = {}, 0
    n_positive = sum(map(len, values.values()))
    for value in sorted(values):
        genes = values[value]
        rank = (offset + (len(genes) + 1) / 2) / n_positive
        ranks.update({gene: rank for gene in genes})
        offset += len(genes)
    result = {}
    for name, genes in gene_sets.items():
        overlap = sorted(set(genes) & panel)
        detected = [gene for gene in overlap if expression.get(gene, 0) > 0]
        score = sum(ranks.get(gene, 0) for gene in overlap) / len(overlap) if len(overlap) >= 3 else None
        result[name] = {
            "score": score, "assayed_genes": overlap, "detected_genes": detected,
            "signature_size": len(genes), "panel_coverage": len(overlap) / len(genes),
            "detected_fraction_of_assayed": len(detected) / len(overlap) if overlap else None,
            "status": "available" if score is not None else "insufficient_panel_overlap",
        }
    return result


def build_guidebook(cells: list[dict], context: dict, gene_sets: dict | None = None) -> dict:
    """Build auditable per-cell action weights and per-type summaries.

    A diagnosis or intervention is passed to Jev as context, but changes these
    weights only through an explicit ``context_modifiers`` rule. Missing panel
    support remains null and cannot be interpreted as biological inactivity.
    """
    c = validate_context(context)
    checked = _validate_cells(cells)
    if gene_sets is None:
        if c["organism"].casefold() not in ("human", "homo sapiens"):
            raise ValueError("Default signatures use human symbols. Supply reviewed gene_sets for another organism.")
        gene_sets = DEFAULT_GENE_SETS
    if not isinstance(gene_sets, dict) or not gene_sets:
        raise ValueError("gene_sets must be a nonempty name-to-gene-list map.")
    signatures = {}
    for name, genes in gene_sets.items():
        if not isinstance(name, str) or not name.strip() or not isinstance(genes, list) or not genes:
            raise ValueError("Each gene set needs a nonempty name and gene list.")
        if any(not isinstance(g, str) or not g.strip() for g in genes):
            raise ValueError("Gene-set symbols must be nonempty strings.")
        signatures[name] = sorted({gene.strip().upper() for gene in genes})
    observed = {g for cell in checked for g in cell["expression"]}
    explicit_panel = "assayed_genes" in c
    panel = set(c.get("assayed_genes", observed))
    if not panel:
        raise ValueError("Declare assayed_genes when the expression matrix has no features.")
    if explicit_panel and observed - panel:
        raise ValueError("Expression contains genes outside context.assayed_genes.")
    records, grouped = {}, defaultdict(list)
    for cell in sorted(checked, key=lambda x: x["cell_id"]):
        lineage = cell["lineage"]
        evidence = _rank_evidence(cell["expression"], panel, signatures)
        weights = deepcopy(BASE_WEIGHTS[lineage])
        if lineage == "T_CD8" and c["recognition"]:
            weights["KILL"] = .1
        # Overrides cannot create a lineage-incompatible action. This prevents
        # a typo or context modifier from making mature neurons divide.
        overrides = c["action_priors"].get(lineage, {})
        if set(overrides) - set(weights):
            raise ValueError(f"Action override is incompatible with lineage {lineage} or recognition licensing.")
        weights.update(overrides)
        prior_weights = dict(weights)
        audit = []
        for action, weight in list(weights.items()):
            signature = ACTION_SIGNATURE.get(action)
            score = evidence.get(signature, {}).get("score")
            if score is not None:
                multiplier = 1 + score
                weights[action] = weight * multiplier
                audit.append({"action": action, "source": "expression_rank_heuristic", "signature": signature,
                              "multiplier": multiplier, "assumption": "weight *= 1 + score; uncalibrated"})
        for index, rule in enumerate(c["context_modifiers"]):
            applies = lineage in rule["lineages"] and modifier_matches(rule, c)
            if applies and rule["action"] in weights:
                weights[rule["action"]] *= rule["multiplier"]
            audit.append({"source": "declared_context_modifier", "rule_index": index,
                          "applied": applies and rule["action"] in weights,
                          "action": rule["action"], "multiplier": rule["multiplier"], "reason": rule["reason"]})
        if any(not math.isfinite(w) for w in weights.values()):
            raise ValueError("Action weights overflowed; reduce prior or modifier magnitudes.")
        total = sum(weights.values())
        if not math.isfinite(total) or total <= 0 or weights["REST"] <= 0:
            raise ValueError("Finite action weights and a positive REST fallback are required.")
        actions = {a: {"weight": w, "prior_weight": prior_weights[a], "probability": w / total, "description": DESCRIPTIONS[a],
                       "provenance": "uncalibrated_assumed_preference_per_decision"}
                   for a, w in sorted(weights.items()) if w > 0}
        expression_total = sum(cell["expression"].values())
        if not math.isfinite(expression_total):
            raise ValueError("Expression total overflowed; supply finite counts within numeric range.")
        record = {"cell_type": cell["cell_type"], "lineage": lineage, "actions": actions,
                  "evidence": evidence, "weight_adjustments": audit,
                  "detected_genes": sum(v > 0 for v in cell["expression"].values()),
                  "expression_total": expression_total,
                  "local_gates_required": ["KILL requires a recognized, live malignant neighbor and effector competence",
                                           "CLEAR_CORPSE requires a neighboring corpse",
                                           "PROLIFERATE requires remaining living population capacity; there is no mechanics or space exclusion"]}
        records[cell["cell_id"]] = record
        grouped[cell["cell_type"]].append(record)
    by_type = {}
    for name, group in sorted(grouped.items()):
        by_type[name] = {
            "n_cells": len(group), "lineage": group[0]["lineage"],
            "actions": sorted({a for item in group for a in item["actions"]}),
            "signatures": {signature: {
                "mean_score": (sum(values) / len(values)) if values else None,
                "available_cells": len(values), "missing_cells": len(group) - len(values)}
                for signature in signatures
                for values in [[item["evidence"][signature]["score"] for item in group
                                if item["evidence"][signature]["score"] is not None]]},
        }
    return {
        "schema_version": "2.0", "context": c, "cells": records, "by_cell_type": by_type,
        "provenance": {
            "guidebook_version": "v2-guidebook-1", "context_sha256": _digest(c),
            "default_action_weights": deepcopy(BASE_WEIGHTS),
            "recognition_licensed_kill_weight": .1,
            "input_sha256": _digest(checked), "gene_sets_sha256": _digest(signatures),
            "gene_sets": signatures, "assayed_genes": sorted(panel),
            "panel_source": "context.assayed_genes" if explicit_panel else "input_feature_union_incomplete_if_sparse",
            "score_method": "mean_positive_expression_rank_with_undetected_assayed_genes_zero; minimum_overlap=3",
            "symbol_normalization": "uppercase; no orthology inference",
            "probability_semantics": "Uncalibrated heuristic action preferences conditional on the offered menu; not learned fate probabilities or rates.",
            "time_semantics": "One preference distribution per decision; changing the decision interval changes the model and requires calibration.",
            "limitations": [
                "Expression is fixed baseline evidence; dynamic states are simulated, not remeasured.",
                "Exhaustion-associated and apoptosis-associated transcripts do not establish functional exhaustion or impending death.",
                "No disease or treatment name guarantees proliferation, killing, exhaustion, or therapeutic response.",
                "Enrichment and missing populations affect interactions and composition; no whole-tissue reconstruction is inferred.",
                "Default weights and expression/context multipliers are assumptions requiring experimental calibration.",
            ],
        },
    }
