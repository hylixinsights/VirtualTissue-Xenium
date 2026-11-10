"""Declared evidence-to-action rules and externally measured per-cell scores.

The v4 guidebook connects six actions to six packaged gene sets through one
global table. On a targeted panel those gene sets can be unassayed, which
leaves every expression multiplier at exactly 1. An ``evidence_policy`` states,
per lineage and action, which evidence may modulate the assumed hazard:

- ``gene_set``: a gene set scored by the simulator on the assayed panel;
- ``measured_score``: a per-cell score computed outside the simulator and
  imported with a provenance manifest (for example saved UCell columns).

What does not change: the action menu, every gate, the kinetic and beta tables
and the multiplier law ``exp(beta*z) / mean_reference(exp(beta*z))``. A rule
cannot add an action, and at most one rule may address a lineage/action pair,
so correlated evidence from the same counts is never multiplied twice.

Measured scores are standardised like the internal recovery AUC: mid-rank
percentile within the same sample, declared type and depth stratum, then
probit, clipped. A missing score stays missing; it is never read as zero.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import json
import math

POLICY_VERSION = "v5-evidence-policy-1"
MANIFEST_VERSION = "v5-measured-evidence-1"
GUIDEBOOK_VERSION = "v5-guidebook-1"
KINDS = ("gene_set", "measured_score")
DIRECTIONS = ("supports", "opposes")
DEFAULT_REVIEW = "requires_biological_review"
# The v4 table, restated as rules. A policy holding exactly these rules for
# every lineage reproduces the v4 hazards; see tests/test_v5_evidence_policy.py.
V4_RULES = (("PROLIFERATE", "cell_cycle"), ("ACTIVATE", "activation"), ("EXHAUST", "exhaustion_associated"),
            ("APOPTOSIS", "apoptosis_associated"), ("SECRETE", "inflammation"), ("KILL", "cytotoxic"))


def _text(value, field, limit=2000):
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise ValueError(f"{field} must be nonempty text (at most {limit} characters).")
    return value.strip()


def v4_equivalent_policy(name="v4_equivalent"):
    """The implicit v4 mapping written out as an explicit policy."""
    return dict(policy_version=POLICY_VERSION, name=name, review_status=DEFAULT_REVIEW,
                rules=[dict(action=action, lineages="*", evidence=evidence, kind="gene_set", direction="supports",
                            rationale="v4 global action-to-signature table, restated without change.")
                       for action, evidence in V4_RULES])


def validate_policy(policy, lineage_rules):
    """Return a normalised policy or fail; never infers a rule from a name."""
    if not isinstance(policy, dict):
        raise ValueError("evidence_policy must be an object.")
    unknown = set(policy) - {"policy_version", "name", "description", "review_status", "rules", "provider_evidence"}
    if unknown:
        raise ValueError(f"Unknown evidence_policy settings: {sorted(unknown)}.")
    if policy.get("policy_version") != POLICY_VERSION:
        raise ValueError(f"evidence_policy.policy_version must be {POLICY_VERSION!r}.")
    result = dict(policy_version=POLICY_VERSION, name=_text(policy.get("name"), "evidence_policy.name", 200),
                  review_status=_text(policy.get("review_status", DEFAULT_REVIEW), "evidence_policy.review_status", 500))
    if "description" in policy:
        result["description"] = _text(policy["description"], "evidence_policy.description")
    rules = policy.get("rules")
    if not isinstance(rules, list):
        raise ValueError("evidence_policy.rules must be a list (use [] for a label-only rulebook).")
    normalised, claimed = [], {}
    for index, rule in enumerate(rules):
        where = f"evidence_policy.rules[{index}]"
        if not isinstance(rule, dict):
            raise ValueError(f"{where} must be an object.")
        unknown = set(rule) - {"action", "lineages", "evidence", "kind", "direction", "rationale", "review_status", "reference"}
        if unknown:
            raise ValueError(f"{where} has unknown settings: {sorted(unknown)}.")
        action = rule.get("action")
        if not isinstance(action, str) or action == "REST":
            raise ValueError(f"{where}.action must be a non-REST action.")
        permitted = sorted(name for name, row in lineage_rules.items() if action in row["allowed_actions"])
        if not permitted:
            raise ValueError(f"{where}.action {action!r} is not permitted for any lineage.")
        lineages = rule.get("lineages")
        if lineages == "*":
            targets = permitted
        else:
            if (not isinstance(lineages, list) or not lineages or len(set(lineages)) != len(lineages)
                    or any(not isinstance(name, str) for name in lineages)):
                raise ValueError(f"{where}.lineages must be '*' or a nonempty list of unique lineage names.")
            unsupported = sorted(set(lineages) - set(lineage_rules))
            if unsupported:
                raise ValueError(f"{where}.lineages names unsupported lineages: {unsupported}.")
            forbidden = sorted(set(lineages) - set(permitted))
            if forbidden:
                raise ValueError(f"{where}: {action} is not a permitted action of {forbidden}; a rule cannot add an action.")
            targets = lineages = sorted(lineages)
        if rule.get("kind") not in KINDS:
            raise ValueError(f"{where}.kind must be one of {KINDS}.")
        if rule.get("direction", "supports") not in DIRECTIONS:
            raise ValueError(f"{where}.direction must be one of {DIRECTIONS}.")
        for lineage in targets:
            if (lineage, action) in claimed:
                raise ValueError(f"{where} and rules[{claimed[(lineage, action)]}] both address {lineage}/{action}; "
                                 "declare one evidence source per lineage and action.")
            claimed[(lineage, action)] = index
        item = dict(action=action, lineages=lineages, evidence=_text(rule.get("evidence"), f"{where}.evidence", 200),
                    kind=rule["kind"], direction=rule.get("direction", "supports"),
                    rationale=_text(rule.get("rationale"), f"{where}.rationale"),
                    review_status=_text(rule.get("review_status", DEFAULT_REVIEW), f"{where}.review_status", 500))
        if "reference" in rule:
            item["reference"] = _text(rule["reference"], f"{where}.reference")
        normalised.append(item)
    result["rules"] = normalised
    shown = policy.get("provider_evidence", "all")
    if shown != "all":
        if (not isinstance(shown, list) or len(set(shown)) != len(shown)
                or any(not isinstance(name, str) or not name.strip() for name in shown)):
            raise ValueError("evidence_policy.provider_evidence must be 'all' or a list of unique evidence names.")
        shown = sorted(shown)
    result["provider_evidence"] = shown
    return result


def rule_index(policy, lineage_rules):
    """{(lineage, action): rule} after expanding '*' to the lineages that permit the action."""
    index = {}
    for rule in policy["rules"]:
        lineages = (sorted(name for name, row in lineage_rules.items() if rule["action"] in row["allowed_actions"])
                    if rule["lineages"] == "*" else rule["lineages"])
        for lineage in lineages:
            index[(lineage, rule["action"])] = rule
    return index


def declared_rule(context, lineage, action):
    """The rule (evidence, kind, direction) that may modulate this lineage/action, or None."""
    policy = context.get("evidence_policy")
    if policy is None:
        from ..v2.guidebook import ACTION_SIGNATURE
        name = ACTION_SIGNATURE.get(action)
        return dict(evidence=name, kind="gene_set", direction="supports") if name else None
    for rule in policy["rules"]:
        if rule["action"] == action and (rule["lineages"] == "*" or lineage in rule["lineages"]):
            return rule
    return None


def declared_evidence(context, lineage, action):
    """(evidence name, kind) that may modulate this lineage/action, or (None, None)."""
    rule = declared_rule(context, lineage, action)
    return (rule["evidence"], rule["kind"]) if rule else (None, None)


def validate_measured_evidence(manifest, cells):
    """Normalise the provenance manifest and require it to cover every imported score."""
    used = sorted({name for cell in cells for name in cell.get("measured_scores", {})})
    if manifest is None:
        if used:
            raise ValueError("Cells carry measured_scores without a measured_evidence manifest; "
                             "declare the method and definition of every score.")
        return None
    if not isinstance(manifest, dict) or not isinstance(manifest.get("scores"), dict) or not manifest["scores"]:
        raise ValueError("measured_evidence requires a nonempty scores map: name -> method and description.")
    if manifest.get("manifest_version", MANIFEST_VERSION) != MANIFEST_VERSION:
        raise ValueError(f"measured_evidence.manifest_version must be {MANIFEST_VERSION!r}.")
    result = deepcopy(manifest)
    result["manifest_version"] = MANIFEST_VERSION
    if "caveat" in result:
        result["caveat"] = _text(result["caveat"], "measured_evidence.caveat")
    for name, spec in result["scores"].items():
        where = f"measured_evidence.scores[{name!r}]"
        if not isinstance(name, str) or not name.strip() or name != name.strip():
            raise ValueError("Measured score names must be nonempty strings without surrounding spaces.")
        if not isinstance(spec, dict):
            raise ValueError(f"{where} must be an object.")
        spec["method"] = _text(spec.get("method"), f"{where}.method", 500)
        spec["description"] = _text(spec.get("description"), f"{where}.description")
        genes = spec.get("genes_scored")
        if genes is not None and (not isinstance(genes, list) or not genes or len(set(genes)) != len(genes)
                                  or any(not isinstance(g, str) or not g.strip() for g in genes)):
            raise ValueError(f"{where}.genes_scored must be a nonempty list of unique gene symbols when given.")
        size = spec.get("n_genes_set")
        if size is not None and (type(size) is not int or size < 1 or (genes is not None and size < len(genes))):
            raise ValueError(f"{where}.n_genes_set must be a positive integer covering genes_scored.")
    undeclared = sorted(set(used) - set(result["scores"]))
    if undeclared:
        raise ValueError(f"Measured scores without a manifest entry: {undeclared}.")
    json.dumps(result, allow_nan=False)  # provenance must be finite JSON
    return result


def check_references(policy, signatures, manifest):
    """Every rule and provider-evidence name must identify declared evidence of its kind."""
    measured = set((manifest or {}).get("scores", {}))
    both = sorted(set(signatures) & measured)
    if both:
        raise ValueError(f"Gene-set and measured-score names collide: {both}.")
    for index, rule in enumerate(policy["rules"]):
        available = signatures if rule["kind"] == "gene_set" else measured
        if rule["evidence"] not in available:
            raise ValueError(f"evidence_policy.rules[{index}] references {rule['kind']} {rule['evidence']!r}, "
                             "which is not declared for this bundle.")
    if policy["provider_evidence"] != "all":
        unknown = sorted(set(policy["provider_evidence"]) - set(signatures) - measured)
        if unknown:
            raise ValueError(f"evidence_policy.provider_evidence names undeclared evidence: {unknown}.")


def score_measured(cells, manifest, scoring, qc_fail):
    """Standardise imported scores within the same references as the internal AUC.

    Returns ({cell_id: {name: record}}, {name: {reference_id: [z, ...]}}), the
    shapes ``scoring.expression_multiplier`` already consumes. An imported score
    has no panel null model, so a type too small for an empirical reference is
    ``insufficient_reference`` and stays neutral under every scoring method.
    """
    from ..v2.scoring import Z_CLIP, _NORMAL, _percentile, depth_bins, reference_partition
    eligible = [c for c in cells if c["cell_id"] not in qc_fail]
    bins = depth_bins(eligible, scoring["depth_bins"])
    settings = dict(scoring, reference_fallback="cell_type") if scoring["method"] == "tiered_panel_auc" else scoring
    evidence, references = {c["cell_id"]: {} for c in cells}, {}
    for name in sorted(manifest["scores"]):
        spec = manifest["scores"][name]
        genes = spec.get("genes_scored")
        # References are the cells that have this score: an unscored cell neither
        # serves as a reference nor decides whether a depth stratum is large enough.
        scored = [c for c in eligible if c.get("measured_scores", {}).get(name) is not None]
        partition = reference_partition(scored, bins, settings)
        groups = defaultdict(list)
        for c in scored:
            groups[partition[c["cell_id"]]].append(c["measured_scores"][name])
        zs = defaultdict(list)
        for c in cells:
            cid = c["cell_id"]
            value = c.get("measured_scores", {}).get(name)
            record = dict(kind="measured_score", score=None, measured_value=value, percentile=None, z=None,
                          reference=None, reference_size=0, depth_bin=bins.get(cid), method=spec["method"],
                          signature_size=spec.get("n_genes_set"), n_genes_scored=len(genes) if genes else None,
                          panel_coverage=(len(genes) / spec["n_genes_set"] if genes and spec.get("n_genes_set") else None),
                          informative=False, calibration="unavailable", multiplier_normalization="neutral")
            if cid in qc_fail:
                record["status"] = "qc_fail_below_min_transcripts"
            elif value is None:
                record["status"] = "not_measured"
            else:
                key = partition[cid]
                values = groups[key]
                record["reference_size"] = len(values)
                if len(values) < scoring["min_reference_cells"]:
                    record["status"] = "insufficient_reference"
                else:
                    percentile = _percentile(value, values)
                    z = max(-Z_CLIP, min(Z_CLIP, _NORMAL.inv_cdf(min(max(percentile, 1e-6), 1 - 1e-6))))
                    reference_id = json.dumps(key, ensure_ascii=False, separators=(",", ":"))
                    record.update(score=percentile, percentile=percentile, z=z, status="available",
                                  reference="cell_type_x_depth_bin" if key[-1] is not None else "cell_type_no_depth",
                                  reference_id=reference_id, reference_sample=key[0],
                                  informative=min(values) != max(values),
                                  calibration="empirical_within_sample_type",
                                  multiplier_normalization="reference_group_mean_one")
                    zs[reference_id].append(z)
            evidence[cid][name] = record
        references[name] = zs
    return evidence, references


def action_evidence(bundle):
    """Per action: how many cells that are offered it have evidence able to modulate it.

    Read from the guidebook's recorded adjustments, so it describes what the
    rulebook actually applied. ``label_only`` lists actions whose hazard is the
    same for every cell of a lineage on this panel.
    """
    context, entries = bundle["context"], bundle["guidebook"]["cells"]
    offered, usable, modulated, sources, statuses = Counter(), Counter(), Counter(), defaultdict(set), defaultdict(Counter)
    for entry in entries.values():
        adjusted = {a["action"]: a for a in entry.get("weight_adjustments", []) if "signature" in a}
        for action in entry["actions"]:
            if action == "REST":
                continue
            offered[action] += 1
            name, kind = declared_evidence(context, entry["lineage"], action)
            if name is None:
                statuses[action]["no_rule_declared"] += 1
                continue
            sources[action].add((name, kind))
            adjustment = adjusted.get(action)
            if adjustment is None:
                statuses[action]["evidence_absent_from_bundle"] += 1
                continue
            statuses[action][adjustment["status"]] += 1
            if adjustment.get("z") is not None:
                usable[action] += 1
            if adjustment.get("multiplier", 1.0) != 1.0:
                modulated[action] += 1
    by_action = {action: dict(cells_offered=offered[action], cells_with_usable_evidence=usable[action],
                              cells_with_changed_hazard=modulated[action],
                              evidence=[dict(name=name, kind=kind) for name, kind in sorted(sources[action])],
                              status_counts=dict(sorted(statuses[action].items())))
                 for action in sorted(offered)}
    label_only = sorted(action for action, row in by_action.items() if row["cells_with_usable_evidence"] == 0)
    warnings = []
    for action in label_only:
        row = by_action[action]
        if row["evidence"]:
            reasons = ", ".join(f"{k}: {v}" for k, v in row["status_counts"].items())
            warnings.append(f"{action}: declared evidence {[e['name'] for e in row['evidence']]} is unusable for all "
                            f"{row['cells_offered']} cells offered this action ({reasons}); its hazard is label-only.")
    return dict(policy=(context["evidence_policy"]["name"] if context.get("evidence_policy") else "v4_default_table"),
                by_action=by_action, label_only_actions=label_only, warnings=warnings,
                note="Usable means a standardised value exists; beta = 0 or a neutral value can still leave the hazard unchanged. "
                     "Actions without a declared rule are label-only by design.")


def decision_evidence_association(run_directory):
    """Do the cells that take an action carry more of its declared evidence than the cells offered it?

    The multiplier preserves each reference group's mean hazard, so evidence is
    not expected to change how many events happen; it changes which cells have
    them. For every lineage/action with declared evidence this compares the mean
    t0 evidence z of the decisions that selected the action with the mean z of
    all decisions in which the action was offered. A label-only rulebook gives a
    difference near zero. Reads saved files only.

    This is a check that sampled decisions follow the declared evidence, in
    local or Jev mode. It is not evidence that the rule is biologically right.
    """
    from pathlib import Path
    from statistics import mean, stdev
    from ..v2.inputs import read_json
    from ..v2.simulation import read_decisions
    directory = Path(run_directory)
    guidebook = read_json(directory / "guidebook.json")
    context, entries = guidebook["context"], guidebook["cells"]
    if context.get("schema_version") != "3.0":
        raise ValueError("Evidence association requires a schema 3.0 run.")
    offered, selected, z_cache = defaultdict(list), defaultdict(list), {}
    for decision in read_decisions(directory, expand=False):
        observation = decision["observation"]
        lineage = observation["lineage"]
        for action in decision["prior"]:
            if action == "REST":
                continue
            name, kind = declared_evidence(context, lineage, action)
            if name is None:
                continue
            reference = observation.get("t0_evidence_ref")
            if reference is not None:
                key = (reference["guidebook_id"], name)
                if key not in z_cache:
                    z_cache[key] = (entries[reference["guidebook_id"]].get("evidence", {}).get(name) or {}).get("z")
                z = z_cache[key]
            else:
                z = (observation.get("t0_evidence", {}).get(name) or {}).get("z")
            if z is None:
                continue
            offered[(lineage, action, name, kind)].append(z)
            if decision["action"] == action:
                selected[(lineage, action, name, kind)].append(z)
    rows = []
    for key in sorted(offered):
        lineage, action, name, kind = key
        chosen, pool = selected.get(key, []), offered[key]
        rows.append(dict(lineage=lineage, action=action, evidence=name, kind=kind,
                         decisions_offered=len(pool), decisions_selected=len(chosen), mean_z_offered=mean(pool),
                         mean_z_selected=mean(chosen) if chosen else None,
                         difference=(mean(chosen) - mean(pool)) if chosen else None,
                         standard_error_selected=(stdev(chosen) / len(chosen) ** .5) if len(chosen) > 1 else None))
    report = read_json(directory / "report.json")
    return dict(run=str(directory.resolve()), mode=report.get("mode"), seed=report.get("seed"),
                completed_steps=report.get("completed_steps"),
                policy=(context["evidence_policy"]["name"] if context.get("evidence_policy") else "v4_default_table"),
                control=(guidebook.get("control") or {}).get("kind"), rows=rows,
                unit="t0 evidence z of the acting cell's ancestral guidebook entry; one value per decision, so a cell "
                     "counts once per step in which the action was offered",
                interpretation="A positive difference means selected decisions carried more of the declared evidence than "
                               "the decisions offered the action. Expect about zero under beta_zero. Under an expression "
                               "shuffle the difference stays positive: z is the evidence in this run's guidebook, which is "
                               "the shuffled evidence that drove its hazards; a shuffle breaks the link between evidence and "
                               "position, not between evidence and decision. Decisions of one cell across steps are not "
                               "independent; the standard error is descriptive.")
