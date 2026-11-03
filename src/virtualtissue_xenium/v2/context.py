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


def _int(value, field, lo, hi):
    if type(value) is not int or not lo <= value <= hi:
        raise ValueError(f"{field} must be an integer from {lo} to {hi}.")
    return value


def _scoring(scoring, *, v3=False, actions=ACTIONS):
    """Validate the expression-evidence method; fill documented defaults."""
    from .scoring import SCORING_DEFAULTS, SCORING_METHODS
    if not isinstance(scoring, dict) or scoring.get("method") not in SCORING_METHODS:
        raise ValueError(f"scoring.method must be one of {SCORING_METHODS}.")
    if scoring["method"] == "rank_legacy":
        if set(scoring) - {"method"}:
            raise ValueError("scoring.method rank_legacy takes no other settings.")
        return dict(scoring)
    defaults = dict(SCORING_DEFAULTS)
    if v3:
        defaults.update(auc_min_rank=1, reference_scope="sample", reference_fallback="null")
    if scoring["method"] == "tiered_panel_auc":
        if not v3:
            raise ValueError("tiered_panel_auc requires schema 3.0.")
        defaults["panel_max_multiplier"] = 2.0
        defaults["reference_fallback"] = "cell_type"
    unknown = set(scoring) - set(defaults)
    if unknown:
        raise ValueError(f"Unknown scoring settings: {sorted(unknown)}.")
    s = dict(defaults, **scoring)
    if v3:
        if s["reference_scope"] not in ("sample", "pooled"):
            raise ValueError("scoring.reference_scope must be sample or pooled.")
        if s["reference_fallback"] not in ("null", "cell_type"):
            raise ValueError("scoring.reference_fallback must be null or cell_type.")
    if s["method"] == "tiered_panel_auc":
        if s["reference_fallback"] != "cell_type":
            raise ValueError("tiered_panel_auc requires its explicit cell_type/panel hierarchy; use depth_conditioned_auc for strict null fallback.")
        if s["reference_scope"] != "sample":
            raise ValueError("tiered_panel_auc keeps references within each sample.")
        s["panel_max_multiplier"] = _number(s["panel_max_multiplier"], "scoring.panel_max_multiplier", positive=True)
        if not 1 <= s["panel_max_multiplier"] <= 3:
            raise ValueError("scoring.panel_max_multiplier must be between 1 and 3.")
    s["auc_top_fraction"] = _number(s["auc_top_fraction"], "scoring.auc_top_fraction", positive=True)
    if s["auc_top_fraction"] > 1:
        raise ValueError("scoring.auc_top_fraction must be at most 1.")
    s["auc_min_rank"] = _int(s["auc_min_rank"], "scoring.auc_min_rank", 1, 100000)
    s["depth_bins"] = _int(s["depth_bins"], "scoring.depth_bins", 1, 50)
    s["min_reference_cells"] = _int(s["min_reference_cells"], "scoring.min_reference_cells", 2, 100000)
    s["min_overlap"] = _int(s["min_overlap"], "scoring.min_overlap", 1, 1000)
    s["min_transcripts"] = _number(s["min_transcripts"], "scoring.min_transcripts")
    beta = s["beta"]
    if isinstance(beta, dict):
        if set(beta) - set(actions) or "REST" in beta:
            raise ValueError("scoring.beta keys must be non-REST actions.")
        if any(type(v) not in (int, float) or not math.isfinite(v) for v in beta.values()):
            raise ValueError("scoring.beta values must be finite numbers (negative allowed).")
        # Actions not listed get beta = 0 (no expression effect).
        s["beta"] = {a: float(v) for a, v in sorted(beta.items())}
    elif type(beta) in (int, float) and math.isfinite(beta):
        s["beta"] = float(beta)
    else:
        raise ValueError("scoring.beta must be a finite number or an action-to-number map.")
    for value in ([s["beta"]] if isinstance(s["beta"], float) else s["beta"].values()):
        if abs(value) > 3:
            raise ValueError("scoring.beta magnitudes above 3 imply >400-fold hazard ratios; declare them only with evidence.")
    return s


def validate_context(context: dict) -> dict:
    """Return an independent normalized declaration, or fail before an API call.

    Optional ``kinetics`` selects how action weights are interpreted:
    ``legacy_per_decision`` (v2 default: weights are preferences per decision,
    independent of ``time_step_hours``) or ``hazard_per_hour`` (weights are
    hazards in events per hour; the kernel converts them with
    competing risks: P(event) = 1 - exp(-sum(k) * dt), with shares k/sum(k).
    The guard bounds sum(k) * dt by ``max_step_probability``, default 0.25).

    Optional ``scoring`` selects the expression-evidence method:
    ``rank_legacy`` (2.0.0a1 default) or ``depth_conditioned_auc`` (see
    scoring.py; settings auc_top_fraction, auc_min_rank, depth_bins,
    min_reference_cells, min_overlap, beta, min_transcripts).

    Optional ``action_priors`` overrides assumed action weights by lineage.
    Optional ``context_modifiers`` are explicit, auditable multiplicative
    assumptions, with optional exact-match ``when`` clauses. They are never
    inferred from a drug or diagnosis name.
    """
    if not isinstance(context, dict):
        raise ValueError("Context must be a JSON object.")
    _finite_json(context)
    c = deepcopy(context)
    if c.get("schema_version") not in ("2.0", "3.0"):
        raise ValueError("context.schema_version must be '2.0' or '3.0'.")
    v3 = c["schema_version"] == "3.0"
    lineages, actions = LINEAGES, ACTIONS
    if v3:
        from ..v3.ontology import lineage_rules, load_platform, load_profile, load_ontology
        from .dynamics import validate_dynamics
        c.setdefault("ontology_version", "3.0")
        rules = lineage_rules(c["ontology_version"])
        lineages = tuple(rules)
        actions = tuple(sorted({a for r in rules.values() for a in r["allowed_actions"]}))
        c.setdefault("platform", c.get("assay"))
        platform = load_platform(c["platform"])
        if c["platform"] != c.get("assay"):
            raise ValueError("platform and assay must identify the same measurement platform.")
        dynamics = c.get("dynamics", {})
        if dynamics is None:
            dynamics = {}
        if isinstance(dynamics, dict):
            # This is a fallback for absent lineage radii, never a bin size.
            dynamics = {"default_radius_um": platform["mean_radius_um"], **dynamics}
        c["dynamics"] = validate_dynamics(dynamics)
        c.setdefault("kernel", "cell_agent")
        if c["kernel"] not in ("cell_agent", "lattice"):
            raise ValueError("kernel must be cell_agent or lattice.")
        if c["kernel"] == "lattice":
            from .lattice import validate_lattice
            c["lattice"] = validate_lattice(c.get("lattice", {}))
        if "tissue_profile" in c:
            load_profile(c["tissue_profile"])
        if "geometry" in c:
            from ..v3.geometry import validate_geometry_config
            c["geometry"] = validate_geometry_config(c["geometry"])
    for field in ("tissue", "organism", "disease"):
        c[field] = _text(c.get(field), field)
    assays = ("xenium", "scrna", "snrna", "atera", "visium_hd", "stereo_seq") if v3 else ("xenium", "scrna", "snrna")
    if c.get("assay") not in assays:
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
    kinetics = c.setdefault("kinetics", {"mode": "hazard_per_hour" if v3 else "legacy_per_decision"})
    if v3 and (not isinstance(kinetics, dict) or kinetics.get("mode") != "hazard_per_hour"):
        raise ValueError("v3 requires hazard_per_hour; use schema 2.0 for legacy kinetics.")
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
    c.setdefault("interaction_mode", "well_mixed" if c["assay"] in ("scrna", "snrna") else "spatial")
    if c["interaction_mode"] not in ("spatial", "well_mixed"):
        raise ValueError("interaction_mode must be spatial or well_mixed.")
    if v3 and c["kernel"] == "lattice" and c["interaction_mode"] != "spatial":
        raise ValueError("The lattice kernel requires measured spatial bin coordinates.")
    c["neighbor_radius_um"] = _number(c.get("neighbor_radius_um", platform["neighbor_radius_um"] if v3 else 30),
                                       "neighbor_radius_um", positive=True)
    c.setdefault("max_neighbors", 20)
    if type(c["max_neighbors"]) is not int or not 1 <= c["max_neighbors"] <= 100:
        raise ValueError("max_neighbors must be an integer from 1 to 100.")
    if v3 and c["assay"] in ("scrna", "snrna") and c["interaction_mode"] != "well_mixed":
        raise ValueError("Dissociated scRNA/snRNA requires well_mixed interactions.")
    if "assayed_genes" in c:
        genes = _text_list(c["assayed_genes"], "assayed_genes")
        normalize_symbol = (lambda g: g) if v3 else str.upper
        if not genes or len({normalize_symbol(g) for g in genes}) != len(genes):
            raise ValueError("assayed_genes must be nonempty and unique after case normalization.")
        c["assayed_genes"] = sorted(normalize_symbol(g) for g in genes)
    scoring_default = {"method": "rank_legacy"}
    if "scoring" in c and not isinstance(c["scoring"], dict):
        raise ValueError("scoring must be an object.")
    if v3:
        scoring_default = {"method": "depth_conditioned_auc", "min_transcripts": platform.get("min_transcripts", 20),
                           "auc_top_fraction": platform.get("auc_top_fraction", .05),
                           "auc_min_rank": platform.get("auc_min_rank", 1),
                           "beta": load_ontology(c["ontology_version"])["kinetics"]["beta_by_action"]}
        if c.get("scoring", {}).get("method", "depth_conditioned_auc") not in ("depth_conditioned_auc", "tiered_panel_auc"):
            raise ValueError("v3 requires depth_conditioned_auc or tiered_panel_auc scoring.")
    c["scoring"] = _scoring(dict(scoring_default, **c.get("scoring", {})), v3=v3, actions=actions)
    priors = c.setdefault("action_priors", {})
    if not isinstance(priors, dict):
        raise ValueError("action_priors must map lineages to action weights.")
    for lineage, weights in priors.items():
        if lineage not in lineages or not isinstance(weights, dict):
            raise ValueError("action_priors keys must be supported lineages with action maps.")
        for action, value in weights.items():
            if action not in actions:
                raise ValueError(f"Unknown action prior: {action}.")
            weights[action] = _number(value, f"action_priors.{lineage}.{action}")
        if "REST" in weights and weights["REST"] <= 0:
            raise ValueError("REST must retain a positive weight.")
    modifiers = c.setdefault("context_modifiers", [])
    if not isinstance(modifiers, list):
        raise ValueError("context_modifiers must be a list.")
    for rule in modifiers:
        if not isinstance(rule, dict) or rule.get("action") not in actions:
            raise ValueError("Every context modifier must name a supported action.")
        rule["lineages"] = _text_list(rule.get("lineages"), "modifier.lineages")
        if not rule["lineages"] or set(rule["lineages"]) - set(lineages):
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
