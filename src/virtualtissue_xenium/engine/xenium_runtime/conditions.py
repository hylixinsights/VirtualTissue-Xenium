"""Three-valued, exhaustive evaluation of the reference rule vocabulary.

Reference rule eligibility is kept separate from an executable action. Missing
measurements, subtype disagreement and unavailable handlers remain visible.
"""
from collections import Counter

ALIASES = {
    "Tumor": {"Tumor", "Tumor_proliferating"},
    "Tumor_proliferating": {"Tumor", "Tumor_proliferating"},
    "T_CD8_cytotoxic": {"T_CD8_cytotoxic", "T_CD8"},
    "T_CD4": {"T_CD4"},
    "B_cell": {"B_cell"},
    "Plasma_cell": {"Plasma_cell", "B_cell"},
    "Macrophage": {"Macrophage", "Myeloid"},
    "CAF": {"CAF", "Fibroblast", "Stromal"},
    "Fibroblast": {"CAF", "Fibroblast", "Stromal"},
    "Endothelial": {"Endothelial"},
}


def identity_gate(cell, accepted, policy="intersection"):
    if not cell or cell.get("role") != "roi":
        return False
    candidates = cell["identity"]["candidates"]
    primary = candidates[0]["label"]
    if primary not in accepted:
        return False
    if policy == "primary_annotation":
        return True
    # A broad secondary label may support the primary capability, but a primary
    # B cell is never promoted to a plasma cell merely by sharing a lineage.
    permitted = set().union(*(ALIASES.get(label, {label}) for label in accepted))
    return all(c["label"] in permitted and c.get("lineage") is not None for c in candidates)


def compare(actual, op, expected):
    if op not in {"eq", "neq", "in", "contains", "contains_any", "lte", "gte"}:
        raise ValueError(f"Unsupported condition operator: {op}")
    if actual is None or expected is None:
        return "unknown"
    operations = {
        "eq": lambda: actual == expected,
        "neq": lambda: actual != expected,
        "in": lambda: actual in expected,
        "contains": lambda: expected in actual,
        "contains_any": lambda: bool(set(actual) & set(expected)),
        "lte": lambda: actual <= expected,
        "gte": lambda: actual >= expected,
    }
    if op not in operations:
        raise ValueError(f"Unsupported condition operator: {op}")
    return "true" if operations[op]() else "false"


def aggregate(results):
    if "false" in results:
        return "false"
    return "unknown" if "unknown" in results else "true"


def value_for(path, cell, target, contact, state, scenario, by_id):
    own = state.get(cell["cell_id"], {})
    other = state.get(target["cell_id"], {}) if target else {}
    if path == "cell.viability_state":
        return own.get("viability", "alive")
    if path == "target.viability_state":
        return other.get("viability", "alive") if target else None
    if path == "contact.polygon_gap_um":
        return contact.get("gap_um") if contact else None
    if path == "cell.genes_detected":
        return cell["genes_detected"]
    if path == "target.genes_detected":
        return target.get("genes_detected", []) if target else None
    if path == "cell.identity.agreement":
        return cell["identity"]["agreement"]
    if path == "target.identity.canonical_type":
        return target["identity"]["candidates"][0]["label"] if target and target.get("identity") else None
    if path == "recognition.cognate_pMHC_I":
        return scenario.get("recognition")
    if path == "cell.effector_competence":
        return scenario.get("effector_competence")
    if path == "cell.degranulation_cooldown_complete":
        return own.get("cooldown_complete", True)
    if path == "cell.processing_capacity_available":
        return own.get("processing_capacity_available") if scenario.get("clearance_competence") else None
    if path == "neighbourhood_30um.cxcl13_source_present":
        # RNA observation can flag a candidate source, never a secreted gradient.
        neighbours = cell.get("neighbour_ids_30um", [])
        return any("CXCL13" in by_id[i].get("genes_detected", []) for i in neighbours if i in by_id)
    if path == "geometry.displaceable_partner_available":
        return False  # fixed-geometry runtime has no displacement implementation
    if path.startswith("perturbation.") or path == "identity_scenario.treg_declared":
        return False
    return None


def evaluate(rule, cell, target, contact, state, scenario, by_id, contract):
    checks = []

    def add(path, result, actual=None, expected=None, evidence="runtime"):
        checks.append(dict(path=path, result=result, actual=actual, expected=expected, evidence=evidence))

    policy = scenario["identity_policy"]
    add("identity.capability", "true" if identity_gate(cell, rule["eligible_identities"], policy) else "false", evidence=policy)
    if not rule.get("enabled_by_default", True):
        add("rule.enabled", "false", evidence="disabled in this runtime")
    add("actor.available", "true" if state.get(cell["cell_id"], {}).get("viability", "alive") == "alive" else "false")
    if target is not None:
        add("target.in_active_domain", "true" if target.get("role") == "roi" else "false", evidence="frozen boundary")
    for condition in rule["required_conditions"].get("all", []):
        path, op, expected = condition["path"], condition["op"], condition["value"]
        actual = value_for(path, cell, target, contact, state, scenario, by_id)
        if path == "target.identity.canonical_type" and target is not None:
            result = "true" if identity_gate(target, expected, policy) else "false"
        elif isinstance(expected, str) and expected.startswith("scenario."):
            expected = scenario.get(expected.split(".", 1)[1])
            result = compare(actual, op, expected)
        else:
            result = compare(actual, op, expected)
        if path.endswith("genes_detected"):
            named = expected if isinstance(expected, list) else [expected]
            # Preserve non-assayed versus assayed-negative distinctions.
            subject = target if path.startswith("target") else cell
            assayed = set(subject.get("panel_genes", [])) if subject else set()
            if result != "true" and any(g not in assayed for g in named):
                result = "unknown"
            actual = {g: "detected" if g in (actual or []) else "not_detected" if g in assayed else "not_assayed" for g in named}
        add(path, result, actual, expected, condition["measurability"])
    if rule["required_conditions"].get("any"):
        raise ValueError("Reference rule contains an unimplemented any-group; do not silently ignore it")
    reference_result = aggregate([c["result"] for c in checks])
    own = state.get(cell["cell_id"], {})
    add("runtime.actor_unoccupied", "true" if own.get("program") is None else "false", evidence="scheduler")
    add("runtime.decision_interval_due", "true" if own.get("decision_due", True) else "false", evidence="scheduler")
    for gate in contract["extra_gates"].get(rule["rule_id"], []):
        actual = state.get(cell["cell_id"], {}).get("pulses_remaining", 0) > 0 if gate == "runtime.effector_store_available" else None
        add(gate, compare(actual, "eq", True), actual, True, "SIM or missing implementation")
    implemented = rule["rule_id"] in contract["handler_rules"]
    add("runtime.handler", "true" if implemented else "false", implemented, True)
    return dict(rule_id=rule["rule_id"], cell_id=cell["cell_id"], target_id=target["cell_id"] if target else None,
                reference_result=reference_result, result=aggregate([c["result"] for c in checks]), checks=checks)


def cell_report(cell, by_id, rules, contract, state=None, scenario=None):
    state = state or {}
    scenario = scenario or contract["scenarios"]["baseline"]
    reports = []
    for rule in rules:
        if not identity_gate(cell, rule["eligible_identities"], "primary_annotation"):
            reports.append(dict(rule_id=rule["rule_id"], result="false", reference_result="false", checks=[dict(path="identity.primary_type", result="false", evidence="observed annotation")], targets=[]))
            continue
        needs_contact = any(c["path"].startswith(("contact.", "target.")) for c in rule["required_conditions"].get("all", []))
        contacts = cell.get("contacts", []) if needs_contact else [None]
        candidates = [evaluate(rule, cell, by_id.get(e["cell_id"]) if e else None, e, state, scenario, by_id, contract) for e in contacts]
        best = min(candidates, key=lambda r: ({"true": 0, "unknown": 1, "false": 2}[r["result"]], sum(c["result"] == "false" for c in r["checks"]), sum(c["result"] == "unknown" for c in r["checks"]))) if candidates else dict(result="false", reference_result="false", checks=[dict(path="contact.exists", result="false", evidence="measured geometry")])
        reports.append(dict(rule_id=rule["rule_id"], result=best["result"], reference_result=best["reference_result"], checks=best["checks"], targets=candidates))
    return reports


def eligibility_summary(atlas, rules, contract):
    by_id = {c["cell_id"]: c for c in atlas["cells"]}
    summary = {r["rule_id"]: dict(true=0, false=0, unknown=0, common_blocks={}) for r in rules}
    blocks = {r["rule_id"]: Counter() for r in rules}
    for cell in atlas["cells"]:
        if cell["role"] != "roi":
            continue
        for report in cell_report(cell, by_id, rules, contract):
            rid = report["rule_id"]
            summary[rid][report["result"]] += 1
            blocks[rid].update(c["path"] for c in report["checks"] if c["result"] != "true")
    for rid in summary:
        summary[rid]["common_blocks"] = dict(blocks[rid].most_common(6))
    return summary
