"""Fixed-state expression counterfactuals: functional influence, not validation.

No transitions are sampled and no provider is called. Labels, geometry, clocks,
neighbours and target ordering are identical in every arm. Only the prepared
expression mapping or its local beta multiplier changes.
"""
from collections import Counter, defaultdict
from copy import deepcopy
from statistics import mean

from ..v2.controls import apply_control
from ..v2.guidebook import ACTION_SIGNATURE
from ..v2.inputs import digest, prepare, validate_prepared
from ..v2.simulation import decision_rows, initialize, js_divergence


def _coverage(bundle):
    grouped = defaultdict(list)
    for entry in bundle["guidebook"]["cells"].values():
        grouped[(entry["lineage"], entry["cell_type"])].append(entry)
    result = []
    for (lineage, cell_type), entries in sorted(grouped.items()):
        signatures = sorted({name for entry in entries for name in entry["evidence"]})
        result.append(dict(lineage=lineage, cell_type=cell_type, n_cells=len(entries),
            n_cells_with_informative_evidence=sum(any(e.get("informative", False) for e in entry["evidence"].values()) for entry in entries),
            n_qc_rest_only=sum(entry.get("qc") != "pass" for entry in entries),
            signatures={name: dict(
                status_counts=dict(sorted(Counter(entry["evidence"][name]["status"] for entry in entries).items())),
                reference_counts=dict(sorted(Counter(entry["evidence"][name].get("reference") or "unavailable" for entry in entries).items())),
                n_with_z=sum(entry["evidence"][name].get("z") is not None for entry in entries),
                n_informative=sum(entry["evidence"][name].get("informative", False) for entry in entries),
                n_without_score=sum(entry["evidence"][name].get("z") is None for entry in entries)) for name in signatures}))
    return result


def _contrast(reference, candidate, bundle, other):
    if set(reference) != set(candidate):
        raise ValueError("Counterfactuals must score the same living cells.")
    rows = []
    for cid, original in reference.items():
        altered = candidate[cid]
        # Expression can change QC eligibility, but cannot change geometry/state.
        for key in ("targets", "corpses", "target_rule"):
            if original[key] != altered[key]:
                raise ValueError("Fixed-state audit acquired different target information.")
        for key in ("cell_type", "lineage", "state", "local", "step", "time_hours"):
            if original["observation"][key] != altered["observation"][key]:
                raise ValueError("Fixed-state audit acquired different labels, state or neighbours.")
        probabilities = {}
        for action in sorted(set(original["prior"]) | set(altered["prior"])):
            before, after = original["prior"].get(action, 0.0), altered["prior"].get(action, 0.0)
            hazard_before = original.get("effective_hazards_per_hour", {}).get(action)
            hazard_after = altered.get("effective_hazards_per_hour", {}).get(action)
            changed = abs(after-before) > 1e-12
            direct = hazard_before != hazard_after
            probabilities[action] = dict(reference=before, counterfactual=after, delta=after-before,
                reference_hazard_per_hour=hazard_before, counterfactual_hazard_per_hour=hazard_after,
                declared_signature=ACTION_SIGNATURE.get(action),
                local_path=("expression_hazard_or_qc" if direct else "competing_risk_redistribution") if changed else "unchanged")
        evidence_changed = digest(original["observation"]["t0_evidence"]) != digest(altered["observation"]["t0_evidence"])
        def provider_evidence(row):
            # Match the compact v3 payload's z, status and signature definitions.
            return {name: {key: spec.get(key) for key in ("z", "status", "reference", "assayed_genes", "signature_size", "panel_coverage")}
                    for name, spec in row["observation"]["t0_evidence"].items()}
        provider_evidence_changed = digest(provider_evidence(original)) != digest(provider_evidence(altered))
        maximum = max(abs(row["delta"]) for row in probabilities.values())
        rows.append(dict(cell_id=cid, cell_type=original["observation"]["cell_type"],
            lineage=original["observation"]["lineage"],
            reference_qc=bundle["guidebook"]["cells"][cid]["qc"], counterfactual_qc=other["guidebook"]["cells"][cid]["qc"],
            action_probabilities=probabilities, max_abs_probability_change=maximum,
            total_variation=.5*sum(abs(row["delta"]) for row in probabilities.values()),
            js_bits=js_divergence(original["prior"], altered["prior"]),
            probability_changed=maximum > 1e-12,
            offered_menu_changed=set(original["actions"]) != set(altered["actions"]),
            expression_evidence_changed=evidence_changed,
            baseline_has_multiple_actions=len(original["actions"]) > 1,
            jev_input_routes=dict(shown_prior_changes=maximum > 1e-12,
                                  direct_expression_evidence_changes=provider_evidence_changed,
                                  opaque_context_hash_changes=digest(bundle["context"]) != digest(other["context"]),
                                  jev_response_effect_measured=False)))
    by_type = defaultdict(list)
    for row in rows:
        by_type[(row["lineage"], row["cell_type"])].append(row)
    def summary(items):
        return dict(n_cells=len(items), n_probability_changed=sum(r["probability_changed"] for r in items),
            n_evidence_changed=sum(r["expression_evidence_changed"] for r in items),
            n_menu_changed=sum(r["offered_menu_changed"] for r in items),
            n_baseline_multiple_actions=sum(r["baseline_has_multiple_actions"] for r in items),
            mean_js_bits=mean(r["js_bits"] for r in items) if items else None,
            mean_total_variation=mean(r["total_variation"] for r in items) if items else None,
            max_abs_probability_change=max((r["max_abs_probability_change"] for r in items), default=None))
    return dict(summary=summary(rows), cells=rows,
        by_type=[dict(lineage=lineage, cell_type=cell_type, **summary(items))
                 for (lineage, cell_type), items in sorted(by_type.items())])


def audit_expression(bundle, seed=0):
    """Return baseline-vs-beta0/within-type/global-shuffle decision contrasts.

The declared seed initializes the one frozen t0 snapshot, target ordering and
the two deterministic permutations. A singleton within-type permutation may
be identical; that is reported rather than forcing an artificial difference.
Beta zero disables local expression modulation, not the evidence sent to Jev.
"""
    validate_prepared(bundle)
    if bundle["schema_version"] != "3.0":
        raise ValueError("Expression audit requires a schema 3.0 prepared bundle.")
    if type(seed) is not int:
        raise ValueError("Audit seed must be an integer.")
    context = deepcopy(bundle["context"])
    context["scoring"]["beta"] = 0.0
    beta_zero = prepare(bundle["cells"], context,
        assayed_genes=bundle["guidebook"]["provenance"]["assayed_genes"],
        gene_sets=bundle["guidebook"]["provenance"]["gene_sets"],
        measurement=bundle.get("measurement"), annotation_provenance=bundle.get("annotation_provenance"))
    arms = dict(baseline=bundle, beta_zero=beta_zero,
        shuffle_expression_within_type=apply_control(bundle, "shuffle_expression_within_type", seed=seed),
        shuffle_expression_global=apply_control(bundle, "shuffle_expression_global", seed=seed))
    snapshot = initialize(bundle, seed=seed)
    indexes, graphs = {}, {}
    for name, prepared in arms.items():
        decisions, graph = decision_rows(snapshot, prepared["guidebook"], prepared["context"], seed=seed, step=0)
        indexes[name] = {row["cell_id"]: row for row in decisions}
        graphs[name] = graph
    if any(graph != graphs["baseline"] for graph in graphs.values()):
        raise ValueError("Expression counterfactual unexpectedly changed the frozen graph.")
    result = dict(schema_version="3.0", audit_version="expression-fixed-state-1", seed=seed,
        source_fingerprint=bundle["fingerprint"], scoring_method=bundle["context"]["scoring"]["method"],
        step=0, time_hours=0.0, frozen_snapshot_sha256=digest(snapshot), frozen_graph_sha256=digest(graphs["baseline"]),
        arms={name: dict(fingerprint=prepared["fingerprint"], coverage_by_type=_coverage(prepared)) for name, prepared in arms.items()},
        comparisons={name: _contrast(indexes["baseline"], index, bundle, arms[name])
                     for name, index in indexes.items() if name != "baseline"},
        jev_response_effect_measured=False,
        interpretation="Functional influence of expression on declared local probabilities at the same t0 state; not biological validation, calibrated activity, temporal prediction or evidence that Jev responds to expression.",
        limitations=[
            "Beta-zero keeps t0 expression evidence available to Jev; it removes the local hazard multiplier only.",
            "Changing beta also changes the full-context hash included in a Jev payload. A clean provider attribution experiment must control this metadata difference separately.",
            "A changed prior affects the shown-prior route; changed evidence can affect either prior mode. Actual Jev response attribution needs matched independent provider probes.",
            "Within-type shuffling of a singleton is necessarily unchanged. Global shuffling can change QC, depth and cross-type biology; menu changes are reported.",
            "Panel-based rare-type evidence is an explicitly bounded modelling heuristic without a biological reference population.",
            "The audit does not optimize beta, signatures or rates to force a difference. Gates and unsupported signatures can legitimately leave probabilities unchanged."])
    result["fingerprint"] = digest(result)
    return result
