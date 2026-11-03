"""Compact evidence and context-text controls retain audited mechanisms."""
from copy import deepcopy
import pytest
from virtualtissue_xenium.v2.provider import payload, JevProvider, ProviderError
from virtualtissue_xenium.v2.inputs import canonical, prepare, read_json
from virtualtissue_xenium.v2.controls import apply_control
from virtualtissue_xenium.v2.simulation import run
from test_v3_contracts import context, cells


def test_payload_shares_definitions_and_preserves_null_evidence():
    evidence = {f"signature{i}": dict(z=None if i == 0 else .5, status="insufficient_reference" if i == 0 else "available",
                assayed_genes=[f"GENE{j}" for j in range(30)], detected_genes=[f"GENE{j}" for j in range(20)],
                signature_size=50, panel_coverage=.6) for i in range(6)}
    rows = [dict(cell_id=str(i), actions={"REST":"Rest", "ACTIVATE":"Activate"},
                 prior={"REST":.9,"ACTIVATE":.1}, observation={"cell_id":str(i),"t0_evidence":evidence}) for i in range(20)]
    old = payload(rows, dict(context(), schema_version="2.0"), "jev-1.14.0")
    compact = payload(rows, context(), "jev-1.14.0", prior_mode="hidden")
    assert len(canonical(compact)) < len(canonical(old)) / 2
    question = compact["questions"]["q0"]["instructions"]
    assert question["observation"]["t0_evidence"]["signature0"] is None
    assert "guidebook_prior" not in question
    assert len(compact["state"]["signature_definitions"]) == 6
    assert rows[0]["observation"]["t0_evidence"] == evidence
    assert compact["state"]["prompt_version"] != old["state"]["prompt_version"]


def test_provider_schema_is_bound_to_journal(tmp_path):
    JevProvider("", tmp_path, {}, transport=lambda _: {}, schema_version="3.0")
    with pytest.raises(ProviderError, match="another run"):
        JevProvider("", tmp_path, {}, transport=lambda _: {}, schema_version="2.0")


def test_scrambled_text_preserves_local_numeric_mechanisms(tmp_path):
    c = context(context_modifiers=[dict(action="PROLIFERATE", lineages=["malignant"],
                multiplier=2, reason="declared test", when={"tissue":"skin"})])
    original = prepare(cells(), c)
    with pytest.raises(ValueError, match="donor_context"):
        apply_control(original, "jev_scrambled_context")
    donor = context(tissue="brain", disease="donor disease")
    changed = apply_control(original, "jev_scrambled_context", donor_context=donor)
    assert changed["context"]["tissue"] == "brain"
    assert changed["guidebook"]["control"]["donor_context_sha256"]
    assert original["guidebook"]["cells"] == changed["guidebook"]["cells"]
    run(original, tmp_path/"original", steps=4, seed=7)
    run(changed, tmp_path/"scrambled", steps=4, seed=7)
    assert read_json(tmp_path/"original/frames.json") == read_json(tmp_path/"scrambled/frames.json")


def test_treatment_name_without_modifiers_is_an_exact_v3_local_null(tmp_path):
    population = cells()
    for i, cell in enumerate(population):
        if i % 2:
            cell.update(cell_type="CD8", lineage="T_CD8")
    c = context()
    treated = deepcopy(c)
    treated["perturbations"] = [dict(type="drug", name="anti-PD1", targets=["PDCD1"],
        mechanism_assumptions=["Illustrative declared assumption"], status="hypothetical")]
    treated["context_modifiers"] = [dict(action="EXHAUST", lineages=["T_CD8"], multiplier=.5,
        reason="Test only", when={"perturbation":"anti-PD1"})]
    baseline = prepare(population, c)
    null = apply_control(prepare(population, treated), "no_modifier_treatment")
    run(baseline, tmp_path/"baseline", steps=16, seed=21)
    run(null, tmp_path/"null", steps=16, seed=21)
    assert read_json(tmp_path/"baseline/frames.json") == read_json(tmp_path/"null/frames.json")
