"""Null-model controls: one declared transformation of a prepared bundle.

Each control rebuilds the guidebook from transformed cells (or a transformed
context) with the same gene sets and scoring rules, and records what was
changed so ``vtx2 compare`` can match the arm to its source run. A control
arm that does not move the readouts outside the seed envelope of the real arm
is evidence that the corresponding input carries no information in this model.
"""
from __future__ import annotations

from copy import deepcopy
import random

from .inputs import digest, prepare, validate_prepared

CONTROLS = {
    "shuffle_expression_within_type": "Permute expression vectors among cells of the same cell_type; labels, lineage and coordinates stay.",
    "shuffle_expression_global": "Permute expression vectors among all cells; labels, lineage and coordinates stay.",
    "label_permutation": "Permute (cell_type, lineage) pairs among cells; expression and coordinates stay.",
    "coordinate_shuffle": "Permute measured (x, y) among cells of the same sample; labels and expression stay.",
    "no_modifier_treatment": "Keep the declared perturbations but remove every context_modifier (drug name alone).",
    "well_mixed_equivalent": "Same cells and context with interaction_mode = well_mixed (coordinates ignored).",
}


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


def apply_control(bundle: dict, kind: str, *, seed: int = 0) -> dict:
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
    if kind == "shuffle_expression_within_type":
        groups = {}
        for c in cells:
            groups.setdefault(c["cell_type"], []).append(c)
        _permute_field(cells, groups, ("expression",), rng)
    elif kind == "shuffle_expression_global":
        _permute_field(cells, {"all": cells}, ("expression",), rng)
    elif kind == "label_permutation":
        _permute_field(cells, {"all": cells}, ("cell_type", "lineage"), rng)
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
    provenance = bundle["guidebook"]["provenance"]
    gene_sets = provenance["gene_sets"]
    assayed = context.get("assayed_genes")
    result = prepare(cells, context, assayed_genes=assayed, gene_sets=gene_sets)
    record = dict(kind=kind, description=CONTROLS[kind], seed=seed,
                  source_fingerprint=bundle["fingerprint"], source_input_sha256=provenance["input_sha256"],
                  source_context_sha256=provenance["context_sha256"],
                  interpretation="Null model. A difference from the source arm within the seed envelope "
                                 "means the transformed input carries no information in this model.")
    result["guidebook"]["control"] = record
    content = {k: v for k, v in result.items() if k != "fingerprint"}
    return dict(content, fingerprint=digest(content))
