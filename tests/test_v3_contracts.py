"""Scientific boundary cases for the v3 measurement and reference contracts."""
from copy import deepcopy
import math
import pytest

from virtualtissue_xenium.v2.context import validate_context
from virtualtissue_xenium.v2.inputs import prepare, validate_prepared, digest
from virtualtissue_xenium.v2.scoring import score_population, expression_multiplier
from virtualtissue_xenium.v2.measurement import read_lineage_map


def context(**overrides):
    result = dict(schema_version="3.0", tissue="skin", organism="human", disease="control",
                  assay="xenium", sampling=dict(design="subsampled", description="test", missing_populations=[]),
                  perturbations=[], recognition=False, time_step_hours=.5)
    result.update(overrides)
    return result


def cells(n=20, sample="a"):
    return [dict(cell_id=f"{sample}{i}", cell_type="tumour", lineage="malignant", sample_id=sample,
                 expression={"MKI67": i+1., "TOP2A": 2., "PCNA": 3., "MCM2": 4., "MCM5": 5., "ACTB": 40.-i},
                 x=float(i)*12, y=0.) for i in range(n)]


def test_reference_minimum_is_not_met_by_borrowing_other_samples():
    scoring = validate_context(context(scoring={"depth_bins": 1}))["scoring"]
    population = cells()+cells(sample="b")
    panel = set(population[0]["expression"])
    evidence, _ = score_population(population, panel, {"cycle": sorted(panel)}, scoring)
    assert {e["cycle"]["status"] for e in evidence.values()} == {"insufficient_reference"}
    pooled, _ = score_population(population, panel, {"cycle": sorted(panel)}, dict(scoring, reference_scope="pooled"))
    assert {e["cycle"]["status"] for e in pooled.values()} == {"available"}


def test_explicit_type_fallback_uses_one_partition_for_percentiles_and_normalisation():
    population = cells(20)
    scoring = validate_context(context(scoring={"depth_bins": 4, "min_reference_cells": 10,
                                                "reference_fallback": "cell_type", "auc_min_rank": 6}))["scoring"]
    signatures = {"cycle": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5"]}
    evidence, references = score_population(population, set(population[0]["expression"]), signatures, scoring)
    values = [e["cycle"] for e in evidence.values()]
    assert {e["reference"] for e in values} == {"cell_type_no_depth"}
    assert len({e["reference_id"] for e in values}) == 1
    assert {e["reference_size"] for e in values} == {20}
    assert math.isclose(sum(expression_multiplier(e, references["cycle"], 1)[0] for e in values)/20, 1)


def test_v3_measurement_roundtrip_and_false_declarations_rejected():
    bundle = prepare(cells(), context(), assayed_genes=list(cells()[0]["expression"]))
    assert bundle["schema_version"] == "3.0"
    assert bundle["measurement"]["n_features"] == 6
    assert validate_prepared(bundle) is bundle
    with pytest.raises(ValueError, match="organism"):
        prepare(cells(), context(), measurement={"organism": "mouse"})
    with pytest.raises(ValueError, match="n_features"):
        prepare(cells(), context(), measurement={"n_features": 5000})
    corrupt = deepcopy(bundle)
    corrupt["measurement"]["n_features"] = 5000
    corrupt["fingerprint"] = digest({k:v for k,v in corrupt.items() if k != "fingerprint"})
    with pytest.raises(ValueError, match="n_features"):
        validate_prepared(corrupt)


def test_v3_symbols_exact_and_lineage_map_fingerprinted(tmp_path):
    population = cells(1)
    population[0]["expression"] = {"Mki67": 30.}
    bundle = prepare(population, context(organism="mouse"), assayed_genes=["Mki67"],
                     gene_sets={"cell_cycle": ["Mki67"]}, lineage_map={"tumour": "malignant"})
    assert bundle["context"]["assayed_genes"] == ["Mki67"]
    assert bundle["guidebook"]["provenance"]["gene_sets"]["cell_cycle"] == ["Mki67"]
    assert bundle["annotation_provenance"]["lineage_map_sha256"] == digest({"tumour":"malignant"})
    path = tmp_path/"map.csv"
    path.write_text("cell_type,lineage\ntumour,malignant\ntumour,epithelial\n")
    with pytest.raises(ValueError, match="unique"):
        read_lineage_map(path)


def test_binned_measurement_cannot_be_misrepresented_as_a_cell():
    with pytest.raises(ValueError, match="lattice"):
        prepare(cells(), context(), measurement={"unit": "binned_spot", "unit_size_um":8})
    with pytest.raises(ValueError, match="well_mixed"):
        validate_context(context(assay="scrna", interaction_mode="spatial"))


def test_prepared_bundle_cannot_silently_use_a_changed_runtime_ontology(monkeypatch):
    from virtualtissue_xenium.v3 import ontology
    bundle = prepare(cells(), context())
    original = ontology.load_ontology
    def changed(version="3.0"):
        record = original(version)
        record["hashes"]["lineages_sha256"] = "changed-runtime-data"
        return record
    monkeypatch.setattr(ontology, "load_ontology", changed)
    with pytest.raises(ValueError, match="ontology hashes"):
        validate_prepared(bundle)
