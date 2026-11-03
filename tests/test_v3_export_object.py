"""Raw object export retains supplied evidence before ROI choice or scoring."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from virtualtissue_xenium.v2.inputs import digest, load_cells, read_json, write_json
from virtualtissue_xenium.v3.roi_catalog import build_roi_catalog

SPEC = importlib.util.spec_from_file_location("export_object_v3", Path(__file__).resolve().parents[1] / "tools/export_object_v3.py")
exporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exporter)


def h5ad_fixture(tmp_path, *, with_lineage=True, with_polygons=True, with_coordinates=True):
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    obs = {"existing_label": ["immune", "tumour"], "sample": ["patient_section", "patient_section"]}
    if with_lineage:
        obs["existing_lineage"] = ["T_CD8", "malignant"]
    obj = ad.AnnData(np.array([[200., 0., 0.], [300., 0., 0.]]),
        obs=pd.DataFrame(obs, index=["c1", "c2"]), var=pd.DataFrame(index=["A", "B", "always_zero"]))
    obj.layers["counts"] = np.array([[1., 0., 0.], [0., 0., 0.]])
    if with_coordinates:
        obj.obsm["measured_xy"] = np.array([[10., 10.], [30., 10.]])
    if with_polygons:
        obj.uns["provided_outlines"] = {
            "c1": np.array([[0., 0.], [20., 0.], [20., 20.], [0., 20.]]),
            "c2": np.array([[20., 0.], [40., 20.], [40., 0.], [20., 20.]]),
        }
    path = tmp_path / "supplied.h5ad"
    obj.write_h5ad(path)
    return path


def h5ad_options(**extra):
    return dict(cell_type_key="existing_label", lineage_key="existing_lineage", sample_key="sample",
                layer="counts", spatial_key="measured_xy", **extra)


def raw_input():
    return dict(cells=[dict(cell_id="x", cell_type="existing label", lineage="myeloid", sample_id="s",
                           expression={"A": 1.0}, x=3., y=4.)], assayed_genes=["A", "unobserved"],
                measurement={"unit": "nucleus", "segmentation_method": "supplied imperfect polygons"},
                annotation_provenance={"source": "existing automated pipeline", "pipeline_version": "1"},
                description="Original object description")


def test_real_h5ad_preserves_low_counts_complete_panel_labels_and_polygon_qc(tmp_path, monkeypatch):
    source = h5ad_fixture(tmp_path)
    measurement = {"unit": "nucleus", "segmentation_method": "user supplied", "organism": "human"}
    write_json(tmp_path / "measurement.json", measurement)
    import virtualtissue_xenium.v2.guidebook as guidebook
    monkeypatch.setattr(guidebook, "build_guidebook", lambda *a, **k: pytest.fail("Raw export must not fit a guidebook"))
    result = exporter.export_object(source, tmp_path / "raw.json", polygon_key="provided_outlines",
        coordinate_scale=.5, measurement=tmp_path / "measurement.json", **h5ad_options())
    assert read_json(tmp_path / "raw.json") == result
    assert result["assayed_genes"] == ["A", "B", "always_zero"]
    assert result["measurement"] == measurement
    assert [(c["cell_type"], c["lineage"], c["sample_id"]) for c in result["cells"]] == [
        ("immune", "T_CD8", "patient_section"), ("tumour", "malignant", "patient_section")]
    assert [c["expression"] for c in result["cells"]] == [{"A": 1.}, {}]
    assert result["cells"][0]["x"] == 5
    assert result["cells"][0]["t0_geometry"]["polygon_um"][1] == [10., 0.]
    assert result["cells"][1]["t0_geometry"]["status"] == "invalid_polygon_centroid_fallback"
    assert result["cells"][1]["t0_geometry"]["polygon_um"] is None
    assert result["annotation_provenance"]["verification"] == "not_independently_verified"
    audit = result["annotation_provenance"]["object_export"]
    assert audit["input_file"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert audit["measurement_sha256"] == digest(measurement)
    assert audit["geometry_qc"]["valid_polygons"] == audit["geometry_qc"]["invalid_polygon_centroid_fallback"] == 1
    assert audit["subset_policy"] == "all_input_cells_no_selection_or_truncation"
    assert "normalization_history_not_verified" in audit["expression_semantics"]
    assert "guidebook" not in result and "fingerprint" not in result and "roi_selection" not in result
    assert load_cells(tmp_path / "raw.json")[1] == result["assayed_genes"]
    before = deepcopy(result)
    catalogue = build_roi_catalog(result, cells_per_roi=1, grid_um=10)
    assert catalogue["candidates"] and result == before


def test_explicit_mapping_without_lineage_column_and_h5ad_bound_before_anndata(tmp_path, monkeypatch):
    source = h5ad_fixture(tmp_path, with_lineage=False, with_polygons=False)
    mapping = tmp_path / "mapping.csv"
    mapping.write_text("cell_type,lineage\nimmune,T_CD8\ntumour,malignant\n")
    result = exporter.export_object(source, tmp_path / "raw.json", lineage_map=mapping, **h5ad_options())
    assert [c["cell_type"] for c in result["cells"]] == ["immune", "tumour"]
    assert [c["lineage"] for c in result["cells"]] == ["T_CD8", "malignant"]
    audit = result["annotation_provenance"]["object_export"]
    assert audit["lineage_map_sha256"] == digest({"immune": "T_CD8", "tumour": "malignant"})
    assert audit["lineage_map_file"]["sha256"] == hashlib.sha256(mapping.read_bytes()).hexdigest()
    import anndata
    monkeypatch.setattr(anndata, "read_h5ad", lambda *a, **k: pytest.fail("Bound must precede AnnData/layer loading"))
    with pytest.raises(ValueError, match="exceeds max_cells=1"):
        exporter.export_object(source, tmp_path / "too_big.json", max_cells=1, lineage_map=mapping, **h5ad_options())
    assert not (tmp_path / "too_big.json").exists()


def test_external_polygons_scale_is_separate_and_does_not_change_source(tmp_path):
    source = h5ad_fixture(tmp_path, with_polygons=False)
    shapes = {"c1": [[0, 0], [10, 0], [10, 10], [0, 10]]}
    polygons = tmp_path / "polygons.json"
    write_json(polygons, shapes)
    result = exporter.export_object(source, tmp_path / "raw.json", coordinate_scale=.5,
        polygons=polygons, polygon_scale=2, **h5ad_options())
    assert result["cells"][0]["x"] == 5
    assert result["cells"][0]["t0_geometry"]["polygon_um"][1] == [20., 0.]
    audit = result["annotation_provenance"]["object_export"]
    assert audit["polygon_join"]["cells_without_polygon"] == 1
    assert audit["polygons_file"]["sha256"] == hashlib.sha256(polygons.read_bytes()).hexdigest()
    assert read_json(polygons) == shapes


def test_raw_json_preserves_declared_metadata_and_rejects_overwrite_before_loading(tmp_path, monkeypatch):
    raw = raw_input()
    raw["context"] = {"assay": "xenium", "organism": "human"}
    source, output = tmp_path / "input.json", tmp_path / "output.json"
    write_json(source, raw)
    result = exporter.export_object(source, output)
    assert result["cells"] == raw["cells"]
    assert result["measurement"] == raw["measurement"]
    assert result["context"] == raw["context"]
    assert result["description"] == raw["description"]
    assert result["annotation_provenance"]["source"] == "existing automated pipeline"
    assert result["annotation_provenance"]["pipeline_version"] == "1"
    assert read_json(source) == raw
    before = output.read_bytes()
    monkeypatch.setattr(exporter, "load_cells", lambda *a, **k: pytest.fail("Existing output must fail before loading"))
    with pytest.raises(ValueError, match="already exists"):
        exporter.export_object(source, output)
    assert output.read_bytes() == before


@pytest.mark.parametrize("panel", [None, [], ["unobserved"], ["A", "A"]])
def test_raw_json_does_not_invent_an_assayed_panel(tmp_path, panel):
    raw = raw_input()
    if panel is None:
        del raw["assayed_genes"]
    else:
        raw["assayed_genes"] = panel
    write_json(tmp_path / "input.json", raw)
    with pytest.raises(ValueError, match="assayed_genes"):
        exporter.export_object(tmp_path / "input.json", tmp_path / "output.json")
    assert not (tmp_path / "output.json").exists()


def test_spatial_export_derives_only_a_valid_polygon_centroid_or_fails(tmp_path):
    raw = raw_input()
    raw["cells"][0].pop("x")
    raw["cells"][0].pop("y")
    write_json(tmp_path / "missing.json", raw)
    with pytest.raises(ValueError, match="requires finite x/y"):
        exporter.export_object(tmp_path / "missing.json", tmp_path / "bad_output.json")
    raw["cells"][0]["polygon"] = [[0, 0], [10, 0], [10, 10], [0, 10]]
    write_json(tmp_path / "polygon.json", raw)
    result = exporter.export_object(tmp_path / "polygon.json", tmp_path / "good_output.json")
    assert result["cells"][0]["x"] == result["cells"][0]["y"] == 5.
    assert result["cells"][0]["t0_geometry"]["anchor_source"] == "polygon_centroid"
    assert result["annotation_provenance"]["object_export"]["geometry_qc"]["derived_centroids"] == 1


def test_cli_reports_no_scoring_or_roi_and_warns_for_explicit_large_limit(tmp_path, capsys):
    write_json(tmp_path / "input.json", raw_input())
    assert exporter.main([str(tmp_path / "input.json"), str(tmp_path / "output.json"), "--max-cells", "2000"]) == 0
    output = capsys.readouterr()
    report = json.loads(output.out)
    assert report["cells"] == 1 and report["scoring"] == "not_computed" and report["roi"] == "not_selected"
    assert "memory use can be substantial" in output.err
    with pytest.raises(SystemExit) as error:
        exporter.main([str(tmp_path / "input.json"), str(tmp_path / "output.json")])
    assert error.value.code == 2


def test_geometry_source_conflicts_and_invalid_scales_refused(tmp_path):
    write_json(tmp_path / "input.json", raw_input())
    for options in ({"max_cells": 0}, {"coordinate_scale": float("nan")}, {"polygon_scale": 2},
                    {"polygons": "unused.json", "polygon_key": "outlines"}, {"layer": "counts"}):
        with pytest.raises(ValueError):
            exporter.export_object(tmp_path / "input.json", tmp_path / "output.json", **options)
    assert not (tmp_path / "output.json").exists()
