"""Strict, portable experiment declarations for the general v2 workflow.

The context is a declaration made by the investigator. Disease and treatment
names never silently turn on a mechanism or specify a desired outcome.
"""
from __future__ import annotations

from copy import deepcopy
import math

KINETICS_MODES = ("legacy_per_decision", "hazard_per_hour")

LINEAGES = (
    "T_CD8", "T_CD4", "B", "myeloid", "fibroblast", "endothelial",
    "epithelial", "malignant", "neuron", "astrocyte", "oligodendrocyte",
    "microglia", "unknown",
)
ACTIONS = (
    "REST", "PROLIFERATE", "APOPTOSIS", "ACTIVATE", "EXHAUST", "RECOVER",
    "SECRETE", "CLEAR_CORPSE", "KILL",
)


def _text(value, field):
    if not isinstance(value, str) or not value.strip() or len(value) > 2000:
        raise ValueError(f"{field} must be nonempty text (at most 2000 characters).")
    return value.strip()


def _number(value, field, *, positive=False):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{field} must be a finite number.")
    if value < 0 or (positive and value == 0):
        raise ValueError(f"{field} must be {'positive' if positive else 'nonnegative'}.")
    return float(value)


def _text_list(value, field):
    if not isinstance(value, list):
        raise ValueError(f"{field} must be a list (use [] when none are declared).")
    result = [_text(x, field) for x in value]
    if len(result) != len(set(result)):
        raise ValueError(f"{field} contains duplicate entries.")
    return result


def _finite_json(value, field="context"):
    """Reject JSON-unrepresentable state, including nonfinite hidden metadata."""
    if value is None or isinstance(value, (str, bool)):
        return
    if type(value) in (int, float):
        if not math.isfinite(value):
            raise ValueError(f"{field} contains a nonfinite value.")
        return
    if isinstance(value, list):
        for item in value:
            _finite_json(item, field)
        return
    if isinstance(value, dict) and all(isinstance(k, str) for k in value):
        for key, item in value.items():
            _finite_json(item, f"{field}.{key}")
        return
    raise ValueError(f"{field} must contain only JSON values and string keys.")


def validate_context(context: dict) -> dict:
    """Return an independent normalized declaration, or fail before an API call.

    Optional ``kinetics`` selects how action weights are interpreted:
    ``legacy_per_decision`` (v2 default: weights are preferences per decision,
    independent of ``time_step_hours``) or ``hazard_per_hour`` (weights are
    hazards in events per hour; the kernel converts them with
    ``p = 1 - exp(-k * dt)`` and refuses a step whose summed probability
    exceeds ``max_step_probability``, default 0.25).

    Optional ``action_priors`` overrides assumed action weights by lineage.
    Optional ``context_modifiers`` are explicit, auditable multiplicative
    assumptions, with optional exact-match ``when`` clauses. They are never
    inferred from a drug or diagnosis name.
    """
    if not isinstance(context, dict):
        raise ValueError("Context must be a JSON object.")
    _finite_json(context)
    c = deepcopy(context)
    if c.get("schema_version") != "2.0":
        raise ValueError("context.schema_version must be '2.0'.")
    for field in ("tissue", "organism", "disease"):
        c[field] = _text(c.get(field), field)
    if c.get("assay") not in ("xenium", "scrna", "snrna"):
        raise ValueError("assay must be xenium, scrna, or snrna.")
    sampling = c.get("sampling")
    if not isinstance(sampling, dict):
        raise ValueError("sampling must declare design, description, and missing_populations.")
    if sampling.get("design") not in ("whole_tissue", "enriched", "subsampled"):
        raise ValueError("sampling.design must be whole_tissue, enriched, or subsampled.")
    sampling["description"] = _text(sampling.get("description"), "sampling.description")
    sampling["missing_populations"] = _text_list(
        sampling.get("missing_populations"), "sampling.missing_populations")
    if not isinstance(c.get("perturbations"), list):
        raise ValueError("perturbations must be a list (use [] for a control).")
    names = []
    for perturbation in c["perturbations"]:
        if not isinstance(perturbation, dict):
            raise ValueError("Every perturbation must be a JSON object.")
        for field in ("type", "name"):
            perturbation[field] = _text(perturbation.get(field), f"perturbation.{field}")
        perturbation["targets"] = _text_list(perturbation.get("targets"), "perturbation.targets")
        perturbation["mechanism_assumptions"] = _text_list(
            perturbation.get("mechanism_assumptions"), "perturbation.mechanism_assumptions")
        perturbation.setdefault("status", "hypothetical")
        if perturbation["status"] not in ("observed", "hypothetical"):
            raise ValueError("perturbation.status must be observed or hypothetical.")
        for field in ("dose", "duration"):
            if field in perturbation:
                perturbation[field] = _text(perturbation[field], f"perturbation.{field}")
        if "evidence" in perturbation:
            perturbation["evidence"] = _text_list(perturbation["evidence"], "perturbation.evidence")
        names.append(perturbation["name"])
    if len(names) != len(set(names)):
        raise ValueError("Perturbation names must be unique.")
    if type(c.get("recognition")) is not bool:
        raise ValueError("recognition must explicitly be true or false; it is an assumption.")
    c["time_step_hours"] = _number(c.get("time_step_hours"), "time_step_hours", positive=True)
    kinetics = c.setdefault("kinetics", {"mode": "legacy_per_decision"})
    if not isinstance(kinetics, dict) or kinetics.get("mode") not in KINETICS_MODES:
        raise ValueError(f"kinetics.mode must be one of {KINETICS_MODES}.")
    if set(kinetics) - {"mode", "max_step_probability"}:
        raise ValueError("kinetics accepts only mode and max_step_probability.")
    if kinetics["mode"] == "hazard_per_hour":
        kinetics["max_step_probability"] = _number(kinetics.get("max_step_probability", 0.25),
                                                   "kinetics.max_step_probability", positive=True)
        if kinetics["max_step_probability"] > 1:
            raise ValueError("kinetics.max_step_probability must be at most 1.")
    elif "max_step_probability" in kinetics:
        raise ValueError("kinetics.max_step_probability applies only to hazard_per_hour.")
    c.setdefault("interaction_mode", "spatial" if c["assay"] == "xenium" else "well_mixed")
    if c["interaction_mode"] not in ("spatial", "well_mixed"):
        raise ValueError("interaction_mode must be spatial or well_mixed.")
    c["neighbor_radius_um"] = _number(c.get("neighbor_radius_um", 30), "neighbor_radius_um", positive=True)
    c.setdefault("max_neighbors", 20)
    if type(c["max_neighbors"]) is not int or not 1 <= c["max_neighbors"] <= 100:
        raise ValueError("max_neighbors must be an integer from 1 to 100.")
    if "assayed_genes" in c:
        genes = _text_list(c["assayed_genes"], "assayed_genes")
        if not genes or len({g.upper() for g in genes}) != len(genes):
            raise ValueError("assayed_genes must be nonempty and unique after case normalization.")
        c["assayed_genes"] = sorted(g.upper() for g in genes)
    priors = c.setdefault("action_priors", {})
    if not isinstance(priors, dict):
        raise ValueError("action_priors must map lineages to action weights.")
    for lineage, weights in priors.items():
        if lineage not in LINEAGES or not isinstance(weights, dict):
            raise ValueError("action_priors keys must be supported lineages with action maps.")
        for action, value in weights.items():
            if action not in ACTIONS:
                raise ValueError(f"Unknown action prior: {action}.")
            weights[action] = _number(value, f"action_priors.{lineage}.{action}")
        if "REST" in weights and weights["REST"] <= 0:
            raise ValueError("REST must retain a positive weight.")
    modifiers = c.setdefault("context_modifiers", [])
    if not isinstance(modifiers, list):
        raise ValueError("context_modifiers must be a list.")
    for rule in modifiers:
        if not isinstance(rule, dict) or rule.get("action") not in ACTIONS:
            raise ValueError("Every context modifier must name a supported action.")
        rule["lineages"] = _text_list(rule.get("lineages"), "modifier.lineages")
        if not rule["lineages"] or set(rule["lineages"]) - set(LINEAGES):
            raise ValueError("Modifier lineages must be a nonempty list of supported lineages.")
        rule["multiplier"] = _number(rule.get("multiplier"), "modifier.multiplier")
        rule["reason"] = _text(rule.get("reason"), "modifier.reason")
        if rule["action"] == "REST" and rule["multiplier"] == 0:
            raise ValueError("A context modifier cannot remove the REST fallback.")
        when = rule.setdefault("when", {})
        if not isinstance(when, dict) or set(when) - {"tissue", "disease", "perturbation"}:
            raise ValueError("Modifier when supports only tissue, disease, and perturbation.")
        for key, value in when.items():
            when[key] = _text(value, f"modifier.when.{key}")
        if "perturbation" in when and when["perturbation"] not in names:
            # An absent named intervention is a legitimate inactive rule for a
            # matched control; its inactive status is retained in the guidebook.
            pass
    return c


def modifier_matches(rule: dict, context: dict) -> bool:
    """Evaluate exact declarations, never diagnosis synonyms or free text."""
    when = rule.get("when", {})
    for key in ("tissue", "disease"):
        if key in when and context[key] != when[key]:
            return False
    return ("perturbation" not in when or
            when["perturbation"] in {p["name"] for p in context["perturbations"]})
