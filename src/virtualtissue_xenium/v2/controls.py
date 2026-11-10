"""Null-model controls: one declared transformation of a prepared bundle.

Each control rebuilds the guidebook from transformed cells (or a transformed
context) with the same gene sets and scoring rules, and records what was
changed so ``vtx2 compare`` can match the arm to its source run. Overlapping
seed envelopes alone do not establish equivalence or absence of information;
coverage, effect-size uncertainty and the declared comparison matter.
"""
from __future__ import annotations

from copy import deepcopy
import random

from .inputs import digest, prepare, validate_prepared

CONTROLS = {
    "jev_scrambled_context": "Replace descriptive tissue/disease/perturbation text from a declared donor context; numeric mechanisms unchanged.",
    "shuffle_expression_within_type": "Permute expression vectors among cells of the same cell_type; labels, lineage and coordinates stay.",
    "shuffle_expression_global": "Permute expression vectors among all cells; labels, lineage and coordinates stay.",
    "label_permutation": "Permute (cell_type, lineage) pairs among cells; expression and coordinates stay.",
    "coordinate_shuffle": "Permute measured (x, y) among cells of the same sample; labels and expression stay.",
    "no_modifier_treatment": "Keep the declared perturbations but remove every context_modifier (drug name alone).",
    "well_mixed_equivalent": "Same cells and context with interaction_mode = well_mixed (coordinates ignored).",
    "beta_zero": "Same cells and context with scoring.beta = 0: label-only kinetics, the reference for expression effects.",
}


# A cell's counts and the scores imported from those counts move together, so an
# expression shuffle cannot leave measured evidence attached to the original cell.
_T0_MEASUREMENTS = ("expression", "measured_scores")


def _rng(bundle, kind, seed):
    return random.Random(int(digest([bundle["fingerprint"], "control", kind, seed]), 16))


def _permute_field(cells, groups, field_names, rng):
    """Permute the named fields among cells within each group, in place."""
    for members in groups.values():
        values = [{f: c.get(f) for f in field_names} for c in members]
        rng.shuffle(values)
        for cell, value in zip(members, values):
            for f in field_names:
                if value[f] is None:
                    cell.pop(f, None)
                else:
                    cell[f] = value[f]


def apply_control(bundle: dict, kind: str, *, seed: int = 0, donor_context=None) -> dict:
    """Return a new prepared bundle for the named control.

    The result carries ``guidebook.control`` with the kind, seed, the source
    bundle fingerprint and the source input digest, and is otherwise a normal
    prepared bundle for ``vtx2 run``.
    """
    validate_prepared(bundle)
    if kind not in CONTROLS:
        raise ValueError(f"Unknown control {kind!r}; choose one of {sorted(CONTROLS)}.")
    if type(seed) is not int:
        raise ValueError("Control seed must be an integer.")
    cells = deepcopy(bundle["cells"])
    context = deepcopy(bundle["context"])
    rng = _rng(bundle, kind, seed)
    if kind == "jev_scrambled_context":
        if not isinstance(donor_context, dict):
            raise ValueError("jev_scrambled_context requires a donor_context from another declared arm.")
        from .context import validate_context, modifier_matches
        donor = validate_context(donor_context)
        if donor["organism"].casefold() != context["organism"].casefold():
            raise ValueError("Scrambled context must preserve organism.")
        # Freeze whether each numerical modifier was active before replacing text.
        context["context_modifiers"] = [dict(rule, when={}) for rule in context["context_modifiers"] if modifier_matches(rule, context)]
        for field in ("tissue", "disease", "perturbations"):
            context[field] = deepcopy(donor[field])
        context.pop("sample_context", None)
        if "sample_context" in donor:
            context["sample_context"] = deepcopy(donor["sample_context"])
    elif kind == "shuffle_expression_within_type":
        groups = {}
        for c in cells:
            groups.setdefault((c["sample_id"], c["cell_type"]) if bundle["schema_version"] == "3.0" else c["cell_type"], []).append(c)
        _permute_field(cells, groups, _T0_MEASUREMENTS, rng)
    elif kind == "shuffle_expression_global":
        _permute_field(cells, {"all": cells}, _T0_MEASUREMENTS, rng)
    elif kind == "label_permutation":
        groups = {"all": cells}
        if bundle["schema_version"] == "3.0":
            groups = {}
            measurement = bundle["measurement"]
            for c in cells:
                if measurement["unit"] == "binned_spot" and (measurement["deconvolution_method"] == "not_declared" or
                        c.get("purity", 0) < measurement["purity_threshold"]):
                    continue
                groups.setdefault(c["sample_id"], []).append(c)
        _permute_field(cells, groups, ("cell_type", "lineage"), rng)
    elif kind == "coordinate_shuffle":
        if any("x" not in c for c in cells):
            raise ValueError("coordinate_shuffle requires measured x/y for every cell.")
        groups = {}
        for c in cells:
            groups.setdefault(c["sample_id"], []).append(c)
        _permute_field(cells, groups, ("x", "y"), rng)
    elif kind == "no_modifier_treatment":
        context["context_modifiers"] = []
    elif kind == "well_mixed_equivalent":
        context["interaction_mode"] = "well_mixed"
    elif kind == "beta_zero":
        if context["scoring"]["method"] not in ("depth_conditioned_auc", "tiered_panel_auc"):
            raise ValueError("beta_zero needs depth_conditioned_auc or tiered_panel_auc in the source bundle.")
        context["scoring"] = dict(context["scoring"], beta=0.0)
    provenance = bundle["guidebook"]["provenance"]
    gene_sets = provenance["gene_sets"]
    assayed = context.get("assayed_genes")
    extra = {}
    if bundle["schema_version"] == "3.0":
        extra = dict(measurement=deepcopy(bundle["measurement"]), annotation_provenance=bundle.get("annotation_provenance"),
                     measured_evidence=deepcopy(bundle.get("measured_evidence")))
    result = prepare(cells, context, assayed_genes=assayed, gene_sets=gene_sets, **extra)
    record = dict(kind=kind, description=CONTROLS[kind], seed=seed,
                  source_fingerprint=bundle["fingerprint"], source_input_sha256=provenance["input_sha256"],
                  source_context_sha256=provenance["context_sha256"],
                  interpretation="Null-model contrast. Overlapping seed envelopes do not prove equivalence or absence of information; assess coverage and paired effect uncertainty.")
    if kind == "jev_scrambled_context":
        record["donor_context_sha256"] = digest(donor)
    if kind == "label_permutation" and bundle["schema_version"] == "3.0":
        record["permutation_scope"] = "within sample; binned units failing purity remain fixed unknown"
    if kind in ("shuffle_expression_within_type", "shuffle_expression_global") and bundle.get("measured_evidence") is not None:
        record["measured_scores"] = "permuted together with the expression vector they were computed from"
    result["guidebook"]["control"] = record
    content = {k: v for k, v in result.items() if k != "fingerprint"}
    return dict(content, fingerprint=digest(content))
