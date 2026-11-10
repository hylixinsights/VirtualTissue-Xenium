"""Collect investigator-declared sample context without inferring mechanisms."""
from __future__ import annotations

from copy import deepcopy

from ..v2.context import validate_context
from ..v3.ontology import lineage_rules


def _answer(input_fn, question, default=None, *, optional=False):
    suffix = f" [{default}]" if default else (" [optional]" if optional else "")
    try:
        value = input_fn(question + suffix + ": ").strip()
    except EOFError as exc:
        raise ValueError("context-init needs interactive input; no context was written.") from exc
    value = value or default
    if not value and not optional:
        raise ValueError(f"An answer is required: {question}")
    return value


def _restrictions(text, rules):
    """Parse only explicitly supported lineage/action pairs, before asking reasons."""
    if not text:
        return []
    pairs = []
    for entry in text.split(","):
        parts = [part.strip() for part in entry.split(":")]
        if len(parts) != 2 or not all(parts):
            raise ValueError("Restrictions must be comma-separated LINEAGE:ACTION pairs.")
        lineage, action = parts[0], parts[1].upper()
        if lineage not in rules:
            raise ValueError(f"Unsupported restriction lineage: {lineage!r}.")
        if action == "REST":
            raise ValueError("A sample restriction cannot disable the REST fallback.")
        if action not in rules[lineage]["allowed_actions"]:
            raise ValueError(f"{action!r} is not a permitted action for lineage {lineage!r}.")
        if (lineage, action) in pairs:
            raise ValueError(f"Duplicate restriction: {lineage}:{action}.")
        pairs.append((lineage, action))
    return pairs


def collect_context(template, *, input_fn=None, output_fn=None):
    """Return a validated independent context; blank answers preserve template values.

    Free text is sample metadata. Only separately confirmed, reasoned restrictions
    become zero-valued context modifiers. No action or numerical rate is added.
    """
    input_fn = input if input_fn is None else input_fn
    output_fn = print if output_fn is None else output_fn
    context = validate_context(template)
    if context["schema_version"] != "3.0":
        raise ValueError("context-init requires a schema-3.0 template.")
    context = deepcopy(context)
    output_fn("Sample context questionnaire. Press Enter to keep each template value.")
    output_fn("Descriptions do not change action permissions. Additional restrictions below are explicit assumptions.")
    output_fn("Mature neurons cannot divide in the ontology; brain glia, progenitors and tumour cells can differ.")
    output_fn("SECRETE emits an abstract inflammatory or suppressive signal. Mucus secretion is not a modelled action.")
    output_fn("Treatment declarations, recognition and numerical settings come from the template; review them before using it.")
    output_fn("Template perturbations: " + (", ".join(p["name"] for p in context["perturbations"]) or "none declared"))
    output_fn("Use 'unknown' for unavailable sample information; an unknown value does not imply healthy tissue.")
    context["organism"] = _answer(input_fn, "Which organism is this sample from?", context["organism"])
    context["tissue"] = _answer(input_fn, "Which tissue or organ was sampled?", context["tissue"])
    context["disease"] = _answer(input_fn, "What disease or pathology is present (use 'none' for healthy tissue)?", context["disease"])
    context["sampling"]["description"] = _answer(
        input_fn, "How was the sample or region selected, and what material does it represent?",
        context["sampling"]["description"])
    sample = dict(context.get("sample_context", {}))
    questions = (
        ("anatomical_region", "Which anatomical region or compartment was sampled?"),
        ("developmental_stage", "What is the age or developmental stage of the sample?"),
        ("cell_state_notes", "Which cell subtypes, injury, regeneration or tumour states affect interpretation?"),
        ("secretion_notes", "Which secretory cell populations or products are known (for example, mucus-producing cells)?"),
    )
    for field, question in questions:
        value = _answer(input_fn, question, sample.get(field), optional=True)
        if value:
            sample[field] = value
    if sample:
        context["sample_context"] = sample
    rules = lineage_rules(context["ontology_version"])
    output_fn("Supported lineages: " + ", ".join(sorted(rules)))
    output_fn("Supported actions by lineage: " + "; ".join(
        f"{lineage}: {', '.join(action for action in row['allowed_actions'] if action != 'REST')}"
        for lineage, row in sorted(rules.items()) if row["allowed_actions"] != ["REST"]))
    requested = _answer(input_fn,
        "Which additional actions should be disabled for this sample? Enter LINEAGE:ACTION pairs separated by commas, or leave blank",
        optional=True)
    for lineage, action in _restrictions(requested, rules):
        reason = _answer(input_fn, f"Why should {lineage}:{action} be disabled in this sample?")
        context["context_modifiers"].append(dict(action=action, lineages=[lineage], multiplier=0, reason=reason))
    return validate_context(context)
