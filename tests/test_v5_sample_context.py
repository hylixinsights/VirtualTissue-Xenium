"""Sample declarations reach both decision paths without inventing biology."""
from copy import deepcopy

import pytest

from virtualtissue_xenium.v2.context import validate_context
from virtualtissue_xenium.v2.controls import apply_control
from virtualtissue_xenium.v2.inputs import prepare
from virtualtissue_xenium.v2.provider import payload
from virtualtissue_xenium.v2.simulation import decision_rows, initialize
from virtualtissue_xenium.v5.probe import _arm_context
from v5_helpers import PANEL, bundle, cells, context


DETAILS = dict(anatomical_region="cerebral cortex", developmental_stage="adult",
               cell_state_notes="Mature neurons with resident glia; no tumour declared.",
               secretion_notes="No mucus-producing population declared.")


def rows_of(prepared):
    return decision_rows(initialize(prepared), prepared["guidebook"], prepared["context"], seed=1, step=0)[0]


def test_sample_details_are_validated_and_copied_without_inserting_defaults():
    plain = validate_context(context())
    assert "sample_context" not in plain
    raw = context(sample_context=deepcopy(DETAILS))
    checked = validate_context(raw)
    assert checked["sample_context"] == DETAILS
    raw["sample_context"]["anatomical_region"] = "changed"
    assert checked["sample_context"] == DETAILS


@pytest.mark.parametrize("value", [None, [], "brain", {"unrecognised": "value"},
                                  {"anatomical_region": " "}, {"developmental_stage": False},
                                  {"secretion_notes": "x" * 2001}])
def test_malformed_sample_details_fail_before_any_request(value):
    with pytest.raises(ValueError, match="sample_context"):
        validate_context(context(sample_context=value))


@pytest.mark.parametrize("with_policy", [False, True])
def test_sample_details_reach_the_rulebook_and_provider_with_or_without_policy(with_policy):
    plain = bundle() if with_policy else prepare(cells(scores=False), context(), assayed_genes=PANEL)
    detailed_context = dict(plain["context"], sample_context=deepcopy(DETAILS))
    prepared = prepare(plain["cells"], detailed_context, assayed_genes=PANEL,
                       measured_evidence=plain.get("measured_evidence"))
    assert prepared["guidebook"]["context"]["sample_context"] == DETAILS
    assert prepared["fingerprint"] != plain["fingerprint"]
    # Descriptive text by itself never changes the local prior or action menu.
    assert prepared["guidebook"]["cells"] == plain["guidebook"]["cells"]
    body = payload(rows_of(prepared), prepared["context"], "jev-1.14.0", "hidden")
    assert body["state"]["context"]["sample_context"] == DETAILS
    assert "not mucus" in body["state"]["sample_context_semantics"]
    assert "sample_context_semantics" not in payload(rows_of(plain), plain["context"], "jev-1.14.0")["state"]


def test_brain_context_keeps_lineage_specific_proliferation_permissions():
    source = cells(scores=False)
    examples = []
    for lineage in ("neuron", "oligodendrocyte", "astrocyte", "microglia", "malignant"):
        for cell in source[:6]:
            examples.append(dict(cell, cell_id=f"{lineage}_{cell['cell_id']}", lineage=lineage, cell_type=lineage))
    prepared = prepare(examples, context(tissue="brain", sample_context=DETAILS), assayed_genes=PANEL)
    actions = {record["lineage"]: record["actions"] for record in prepared["guidebook"]["cells"].values()}
    assert all("PROLIFERATE" not in actions[name] for name in ("neuron", "oligodendrocyte"))
    assert all("PROLIFERATE" in actions[name] for name in ("astrocyte", "microglia", "malignant"))
    assert all("MUCUS" not in action for record in actions.values() for action in record)


def test_declared_restriction_survives_high_expression_controls_and_provider_menu():
    prepared = bundle(sample_context=DETAILS, context_modifiers=[dict(
        action="PROLIFERATE", lineages=["malignant"], multiplier=0,
        reason="Explicit sample-specific restriction for this engineering test.")])
    # Even the cells with the highest imported cycling score cannot bypass a hard restriction.
    assert any(item["evidence"]["S_CYCLE"]["z"] > 1 for item in prepared["guidebook"]["cells"].values())
    for arm in (prepared, apply_control(prepared, "beta_zero"),
                apply_control(prepared, "shuffle_expression_within_type", seed=4)):
        for row in rows_of(arm):
            if row["observation"]["lineage"] == "malignant":
                assert "PROLIFERATE" not in row["actions"]
        body = payload(rows_of(arm), arm["context"], "jev-1.14.0", "hidden")
        for question in body["questions"].values():
            if question["instructions"]["observation"]["lineage"] == "malignant":
                assert "PROLIFERATE" not in question["criteria"]


@pytest.mark.parametrize("donor_details", [None, dict(anatomical_region="airway", secretion_notes="Goblet cells present.")])
def test_context_exchange_replaces_or_removes_original_sample_details(donor_details):
    prepared = bundle(tissue="brain", sample_context=DETAILS)
    donor = context(tissue="lung", disease="control")
    if donor_details is not None:
        donor["sample_context"] = donor_details
    probe_context = _arm_context("context_exchanged", prepared["context"], donor)
    control_context = apply_control(prepared, "jev_scrambled_context", donor_context=donor)["context"]
    for changed in (probe_context, control_context):
        assert changed.get("sample_context") == donor_details
        assert changed["tissue"] == "lung"
    assert prepared["context"]["sample_context"] == DETAILS


def test_probe_can_exchange_only_sample_details():
    source = validate_context(context(sample_context=DETAILS))
    donor = dict(source, sample_context=dict(DETAILS, developmental_stage="developing"))
    changed = _arm_context("context_exchanged", source, donor)
    assert changed["sample_context"]["developmental_stage"] == "developing"
