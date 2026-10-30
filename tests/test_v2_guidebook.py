"""Biological boundaries and evidence/provenance guarantees for the v2 guidebook."""
from copy import deepcopy
import json

import pytest

from virtualtissue_xenium.v2.context import validate_context
from virtualtissue_xenium.v2.guidebook import build_guidebook


def context():
    return {
        "schema_version": "2.0", "tissue": "skin", "organism": "human", "disease": "melanoma",
        "assay": "xenium", "sampling": {
            "design": "subsampled", "description": "Illustrative selected region, not a biological validation.",
            "missing_populations": ["cells outside the selected region"],
        },
        "perturbations": [], "recognition": False, "time_step_hours": 1,
        "assayed_genes": ["MKI67", "TOP2A", "PCNA", "ACTB", "PDCD1"],
    }


def cells(lineage="malignant"):
    return [{"cell_id": "c1", "cell_type": "reviewed_annotation", "lineage": lineage,
             "expression": {"MKI67": 5, "TOP2A": 4, "PCNA": 3, "ACTB": 10}, "x": 0, "y": 0}]


def test_treatment_name_is_context_and_does_not_force_a_response():
    control = context()
    treated = deepcopy(control)
    treated["perturbations"] = [{"type": "treatment", "name": "anti-PD1", "targets": ["PDCD1"],
                                 "mechanism_assumptions": ["Receptor occupancy is not measured."]}]
    a = build_guidebook(cells("T_CD8"), control)
    b = build_guidebook(cells("T_CD8"), treated)
    assert b["context"]["perturbations"][0]["status"] == "hypothetical"
    assert a["cells"]["c1"]["actions"] == b["cells"]["c1"]["actions"]
    assert "KILL" not in b["cells"]["c1"]["actions"]
    assert a["provenance"]["context_sha256"] != b["provenance"]["context_sha256"]


def test_explicit_modifier_changes_only_the_declared_action_and_lineage():
    c = context()
    c["context_modifiers"] = [{"action": "ACTIVATE", "lineages": ["T_CD8"], "multiplier": 2,
                              "reason": "Sensitivity analysis, not estimated treatment efficacy.",
                              "when": {"perturbation": "anti-PD1"}}]
    inactive = build_guidebook(cells("T_CD8"), c)
    c["perturbations"] = [{"type": "treatment", "name": "anti-PD1", "targets": ["PDCD1"],
                           "mechanism_assumptions": []}]
    active = build_guidebook(cells("T_CD8"), c)
    a, b = inactive["cells"]["c1"], active["cells"]["c1"]
    assert b["actions"]["ACTIVATE"]["weight"] == 2 * a["actions"]["ACTIVATE"]["weight"]
    assert b["actions"]["REST"]["weight"] == a["actions"]["REST"]["weight"]
    assert b["weight_adjustments"][-1]["applied"] is True
    assert a["weight_adjustments"][-1]["applied"] is False


def test_missing_panel_evidence_stays_null_and_zero_is_not_missing():
    c = context()
    obj = cells()
    obj[0]["expression"] = {"ACTB": 10}
    g = build_guidebook(obj, c)
    evidence = g["cells"]["c1"]["evidence"]
    assert evidence["cell_cycle"]["score"] == 0
    assert evidence["exhaustion_associated"]["score"] is None
    assert evidence["exhaustion_associated"]["panel_coverage"] == pytest.approx(1/6)
    assert evidence["cytotoxic"]["score"] is None
    assert evidence["cytotoxic"]["panel_coverage"] == 0


def test_rank_evidence_is_scale_invariant_and_inputs_are_unchanged():
    c, obj = context(), cells()
    original = deepcopy((c, obj))
    g = build_guidebook(obj, c)
    scaled = deepcopy(obj)
    scaled[0]["expression"] = {k: v * 100 for k, v in scaled[0]["expression"].items()}
    assert g["cells"]["c1"]["evidence"] == build_guidebook(scaled, c)["cells"]["c1"]["evidence"]
    assert (c, obj) == original
    assert sum(a["probability"] for a in g["cells"]["c1"]["actions"].values()) == pytest.approx(1)
    json.dumps(g, allow_nan=False)


def test_lineage_and_recognition_are_hard_boundaries():
    c = context()
    assert set(build_guidebook(cells("unknown"), c)["cells"]["c1"]["actions"]) == {"REST"}
    assert "PROLIFERATE" not in build_guidebook(cells("neuron"), c)["cells"]["c1"]["actions"]
    c["recognition"] = True
    assert "KILL" in build_guidebook(cells("T_CD8"), c)["cells"]["c1"]["actions"]
    assert "KILL" not in build_guidebook(cells("T_CD4"), c)["cells"]["c1"]["actions"]
    c["action_priors"] = {"neuron": {"PROLIFERATE": 1}}
    with pytest.raises(ValueError, match="incompatible"):
        build_guidebook(cells("neuron"), c)


def test_type_annotation_cannot_hide_a_lineage_conflict():
    obj = cells()
    obj.append(dict(obj[0], cell_id="c2", lineage="T_CD8"))
    with pytest.raises(ValueError, match="multiple lineages"):
        build_guidebook(obj, context())


@pytest.mark.parametrize("field,value", [
    ("recognition", None), ("recognition", "false"), ("time_step_hours", 0),
    ("time_step_hours", float("nan")), ("time_step_hours", True),
    ("max_neighbors", 0), ("max_neighbors", 101), ("interaction_mode", "infer_contacts"),
])
def test_ambiguous_or_invalid_context_is_rejected(field, value):
    c = context()
    c[field] = value
    with pytest.raises(ValueError):
        validate_context(c)


def test_panel_contradiction_nonfinite_and_duplicate_input_rejected():
    obj = cells()
    obj[0]["expression"]["UNDECLARED"] = 1
    with pytest.raises(ValueError, match="outside"):
        build_guidebook(obj, context())
    obj = cells()
    obj[0]["expression"]["ACTB"] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        build_guidebook(obj, context())
    obj = cells()
    obj[0]["expression"]["Actb"] = 1
    with pytest.raises(ValueError, match="unique"):
        build_guidebook(obj, context())
    with pytest.raises(ValueError, match="unique"):
        build_guidebook(cells() + cells(), context())


def test_species_requires_reviewed_gene_sets():
    c = context()
    c["organism"] = "mouse"
    with pytest.raises(ValueError, match="another organism"):
        build_guidebook(cells(), c)
    g = build_guidebook(cells(), c, {"custom_cycle": ["Mki67", "Top2a", "Pcna"]})
    assert g["cells"]["c1"]["evidence"]["custom_cycle"]["score"] is not None


def test_observed_perturbation_metadata_is_validated_and_descriptive():
    c = context()
    c["perturbations"] = [{"type": "treatment", "name": "anti-PD1", "status": "observed",
                           "targets": ["PDCD1"], "mechanism_assumptions": [],
                           "dose": "Dose not available", "duration": "Baseline specimen before second dose",
                           "evidence": ["https://example.org/source-record"]}]
    normalized = validate_context(c)
    assert normalized["perturbations"][0] == c["perturbations"][0]
    assert build_guidebook(cells(), c)["cells"]["c1"]["actions"] == build_guidebook(cells(), context())["cells"]["c1"]["actions"]
    c["perturbations"][0]["status"] = "guaranteed_response"
    with pytest.raises(ValueError, match="status"):
        validate_context(c)
    c["perturbations"][0]["status"] = "observed"
    c["perturbations"][0]["dose"] = 5
    with pytest.raises(ValueError, match="dose"):
        validate_context(c)
