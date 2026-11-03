"""Data-driven permissions and species contracts, beyond the original 13 lineages."""
import hashlib
import json

import pytest

from virtualtissue_xenium.v2.guidebook import BASE_WEIGHTS, DEFAULT_GENE_SETS, build_guidebook
from virtualtissue_xenium.v3.ontology import (
    DATA_ROOT, lineage_rule, load_ontology, load_platform, load_signatures, validate_ontology,
)


def context(**overrides):
    c = dict(schema_version="3.0", tissue="skin", organism="human", disease="melanoma",
             assay="xenium", time_step_hours=.1, recognition=True, perturbations=[],
             sampling=dict(design="subsampled", description="ontology fixture", missing_populations=[]))
    c.update(overrides)
    return c


def cell(lineage="NK", expression=None):
    return dict(cell_id=lineage, cell_type=lineage, lineage=lineage, sample_id="s1",
                x=0., y=0., expression=expression or {"MKI67": 40., "ACTB": 40.})


def test_new_lineages_are_data_driven_and_unknown_stays_rest_only():
    ontology = load_ontology()
    new = {"NK", "T_reg", "plasma", "mast", "pericyte", "melanocyte", "keratinocyte", "neutrophil", "dendritic"}
    assert new <= set(ontology["lineages"])
    book = build_guidebook([cell(name) for name in ontology["lineages"]], context())
    assert "KILL" in book["cells"]["NK"]["actions"]
    assert "MOVE" in book["cells"]["NK"]["actions"]
    assert "KILL" not in book["cells"]["T_reg"]["actions"]
    assert book["cells"]["T_reg"]["lineage_rules"]["secretion_kind"] == "suppressive"
    assert "MOVE" not in book["cells"]["neuron"]["actions"]
    assert "PROLIFERATE" not in book["cells"]["plasma"]["actions"]
    assert set(book["cells"]["unknown"]["actions"]) == {"REST"}
    assert book["provenance"]["ontology_hashes"]["lineages_sha256"]
    assert book["schema_version"] == "3.0"


def test_permissions_cannot_be_added_by_an_override_or_bad_data():
    with pytest.raises(ValueError, match="incompatible"):
        build_guidebook([cell("T_reg")], context(action_priors={"T_reg": {"KILL": .1}}))
    with pytest.raises(ValueError, match="allowed_actions"):
        ontology = load_ontology()
        ontology["lineages"]["NK"].pop("allowed_actions")
        validate_ontology(ontology, ontology["kinetics"])
    with pytest.raises(ValueError, match="Kinetics actions"):
        ontology = load_ontology()
        ontology["kinetics"]["hazards"]["neuron"]["PROLIFERATE"] = {"hazard_per_hour": .1}
        validate_ontology(ontology, ontology["kinetics"])


def test_recognition_movement_and_profile_are_explicit():
    book = build_guidebook([cell("NK")], context(recognition=False, dynamics={"movement_enabled": False}))
    assert not ({"KILL", "MOVE"} & set(book["cells"]["NK"]["actions"]))
    book = build_guidebook([cell("neuron")], context(tissue_profile="skin_melanoma"))
    assert book["provenance"]["unexpected_lineages"] == ["neuron"]
    assert load_platform("xenium")["min_transcripts"] == 20
    with pytest.raises(ValueError, match="profile organism"):
        build_guidebook([cell("neuron")], context(tissue_profile="brain_cortex"))


def test_mouse_symbols_are_explicit_and_not_uppercased():
    sets, provenance = load_signatures("mouse")
    assert "Mki67" in sets["cell_cycle"] and "MKI67" not in sets["cell_cycle"]
    assert provenance["ortholog_table_sha256"]
    book = build_guidebook([cell("microglia", {"Mki67": 40., "Actb": 40.})], context(organism="mouse"))
    assert book["provenance"]["assayed_genes"] == ["Actb", "Mki67"]
    assert "Mki67" in book["provenance"]["gene_sets"]["cell_cycle"]
    with pytest.raises(ValueError, match="outside context.assayed_genes"):
        build_guidebook([cell(expression={"MKI67": 40.})], context(organism="mouse", assayed_genes=["Mki67"]))


def test_explicit_mgi_pairs_cover_mapped_membership_and_exclude_ambiguity():
    human, _ = load_signatures("human")
    mouse, _ = load_signatures("mouse")
    table = json.loads((DATA_ROOT / "signatures/human_mouse_orthologs_v3.json").read_text())
    mapping = {row["human_symbol"]: row["mouse_symbol"] for row in table["pairs"]}
    assert len(mapping) == len(table["pairs"])
    assert all(row["mouse_mgi_id"].startswith("MGI:") for row in table["pairs"])
    for name, genes in human.items():
        if name != "MHC_II":
            assert set(mouse[name]) == {mapping[g] for g in genes if g in mapping}
    assert any(row["human_symbol"] == "GNLY" for row in table["excluded"])


def test_data_loads_are_independent_and_legacy_export_is_exact():
    rule = lineage_rule("NK")
    rule["allowed_actions"].append("INVENTED")
    assert "INVENTED" not in lineage_rule("NK")["allowed_actions"]
    snapshot = json.loads((DATA_ROOT / "ontology/lineages_v2_legacy.json").read_text())
    assert snapshot["weights"] == BASE_WEIGHTS
    assert snapshot["signatures"] == DEFAULT_GENE_SETS
    ontology = load_ontology()
    assert ontology["hashes"]["lineages_sha256"] == hashlib.sha256((DATA_ROOT / "ontology/lineages_v3.json").read_bytes()).hexdigest()
    from virtualtissue_xenium.v2 import datasets
    assert json.loads((DATA_ROOT / "datasets/manifests.json").read_text()) == datasets.list_datasets()
    assert ontology["hashes"]["dataset_registry_sha256"] == hashlib.sha256((DATA_ROOT / "datasets/manifests.json").read_bytes()).hexdigest()
    assert "sum(k_a*m_a)" in ontology["kinetics"]["step_probability"]
    assert all(row["confidence"] == "assumed" for rows in ontology["kinetics"]["hazards"].values() for row in rows.values())


def test_unassayed_signature_is_warning_not_inactive_biology():
    book = build_guidebook([cell()], context())
    assert "MHC_II" in book["provenance"]["coverage_warnings"]
    assert book["cells"]["NK"]["evidence"]["MHC_II"]["score"] is None
