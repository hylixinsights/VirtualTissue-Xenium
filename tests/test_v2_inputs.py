"""Actual AnnData alignment, missingness and sampling boundaries for v2 input."""
from copy import deepcopy

import pytest

ad = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sp = pytest.importorskip("scipy.sparse")

from virtualtissue_xenium.v2.inputs import from_anndata, load_cells, prepare, validate_prepared
from virtualtissue_xenium.v2.simulation import initialize, interaction_graph


def context():
    return dict(schema_version="2.0", tissue="skin", organism="human", disease="melanoma",
                assay="xenium", sampling=dict(design="subsampled", description="Input test fixture",
                missing_populations=[]), perturbations=[], recognition=False, time_step_hours=1)


def annotated():
    obj = ad.AnnData(
        X=sp.csr_matrix([[8., 0., 2.], [1., 7., 0.]]),
        obs=pd.DataFrame({"cell_type": ["reviewed tumor", "reviewed CD8"],
                          "lineage": ["malignant", "T_CD8"],
                          "sample_id": ["section_A", "section_B"]}, index=["cell_z", "cell_a"]),
        var=pd.DataFrame(index=["ACTB", "MKI67", "TOP2A"]),
    )
    obj.obsm["spatial"] = np.array([[0., 0.], [1., 0.]])
    return obj


def test_anndata_alignment_survives_sorted_cell_export_and_backed_read(tmp_path):
    obj = annotated()
    cells, genes = from_anndata(obj)
    assert genes == ["ACTB", "MKI67", "TOP2A"]
    assert [cell["cell_id"] for cell in cells] == ["cell_a", "cell_z"]
    assert cells[0]["expression"] == {"ACTB": 1., "MKI67": 7.}
    assert cells[0]["lineage"] == "T_CD8" and cells[0]["x"] == 1.
    assert cells[1]["expression"] == {"ACTB": 8., "TOP2A": 2.}
    path = tmp_path / "annotated.h5ad"
    obj.write_h5ad(path)
    assert load_cells(path) == (cells, genes)


@pytest.mark.parametrize("duplicate", ["cells", "genes"])
def test_duplicate_identifiers_rejected_before_alignment(duplicate):
    obj = annotated()
    if duplicate == "cells":
        obj.obs_names = ["same", "same"]
    else:
        obj.var_names = ["same", "same", "different"]
    with pytest.raises(ValueError, match="unique"):
        from_anndata(obj)


def test_missing_annotations_and_unknown_guidebook_cell_ids_are_rejected():
    obj = annotated()
    obj.obs.loc["cell_a", "lineage"] = None
    with pytest.raises(ValueError, match="Missing annotation"):
        from_anndata(obj)
    cells, genes = from_anndata(annotated())
    bundle = prepare(cells, context(), assayed_genes=genes)
    entry = bundle["guidebook"]["cells"].pop("cell_a")
    bundle["guidebook"]["cells"]["unknown_id"] = entry
    from virtualtissue_xenium.v2.inputs import digest
    bundle["fingerprint"] = digest({k: v for k, v in bundle.items() if k != "fingerprint"})
    with pytest.raises(ValueError, match="cell IDs"):
        validate_prepared(bundle)


@pytest.mark.parametrize("value", [-1., float("nan"), float("inf")])
def test_scaled_or_nonfinite_expression_is_not_silently_dropped(value):
    obj = annotated()
    matrix = obj.X.toarray()
    matrix[0, 1] = value
    obj.X = matrix
    with pytest.raises(ValueError, match="nonnegative finite"):
        from_anndata(obj)


def test_explicit_counts_layer_and_coordinate_conversion():
    obj = annotated()
    obj.layers["counts"] = obj.X.copy()
    obj.X = np.full(obj.shape, -1.)  # scaled X is invalid, but counts remain valid.
    cells, _ = from_anndata(obj, layer="counts", coordinate_scale=2)
    assert cells[0]["expression"] == {"ACTB": 1., "MKI67": 7.}
    assert cells[0]["x"] == 2.


@pytest.mark.parametrize("mode", ["spatial", "well_mixed"])
def test_canonical_sample_column_prevents_false_cross_section_edges(mode):
    cells, genes = from_anndata(annotated())
    assert {cell["sample_id"] for cell in cells} == {"section_A", "section_B"}
    c = context()
    c["interaction_mode"] = mode
    bundle = prepare(cells, c, assayed_genes=genes)
    graph = interaction_graph(initialize(bundle), bundle["context"])
    assert graph == {"cell_a": [], "cell_z": []}


def test_case_normalized_panel_agreement_and_declared_panel_conflicts():
    cells, genes = from_anndata(annotated())
    c = context()
    c["assayed_genes"] = [gene.lower() for gene in genes]
    bundle = prepare(cells, c, assayed_genes=genes)
    assert bundle["context"]["assayed_genes"] == sorted(genes)
    c["assayed_genes"] = ["ACTB"]
    with pytest.raises(ValueError, match="disagree"):
        prepare(cells, c, assayed_genes=genes)


def test_fingerprint_detects_guidebook_change():
    cells, genes = from_anndata(annotated())
    bundle = prepare(cells, context(), assayed_genes=genes)
    modified = deepcopy(bundle)
    modified["guidebook"]["cells"]["cell_z"]["actions"]["REST"]["weight"] *= 2
    with pytest.raises(ValueError, match="changed"):
        validate_prepared(modified)


def test_full_5000_feature_panel_retains_assayed_zeros_and_rowwise_reads(monkeypatch):
    # A controlled 5K-size feature panel: synthetic features plus known human
    # signature symbols. This checks the adapter, not a biological dataset.
    genes = ["ACTB", "MKI67", "TOP2A", "PCNA"] + [f"TEST_GENE_{i}" for i in range(4996)]
    matrix = sp.csr_matrix(([10., 5.], ([0, 1], [0, 1])), shape=(2, 5000))
    obj = ad.AnnData(matrix, obs=annotated().obs.copy(), var=pd.DataFrame(index=genes))
    obj.obsm["spatial"] = annotated().obsm["spatial"]
    calls, original = [], sp.csr_matrix.toarray
    def observed_toarray(self, *args, **kwargs):
        calls.append(self.shape)
        return original(self, *args, **kwargs)
    monkeypatch.setattr(sp.csr_matrix, "toarray", observed_toarray)
    cells, panel = from_anndata(obj)
    assert calls == [(1, 5000), (1, 5000)]
    assert len(panel) == 5000
    assert all(len(cell["expression"]) == 1 for cell in cells)
    bundle = prepare(cells, context(), assayed_genes=panel)
    evidence = bundle["guidebook"]["cells"]["cell_a"]["evidence"]["cell_cycle"]
    assert evidence["status"] == "available"
    assert evidence["assayed_genes"] == ["MKI67", "PCNA", "TOP2A"]
    assert evidence["detected_genes"] == ["MKI67"]
    assert evidence["score"] == pytest.approx(1/3)
    assert len(bundle["guidebook"]["provenance"]["assayed_genes"]) == 5000
