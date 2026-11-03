"""ROI choices bind to exact data and remain human choices, including unknowns."""
from copy import deepcopy
import json
import re
import shutil
import subprocess

import pytest

from virtualtissue_xenium.v2.inputs import digest, prepare
from virtualtissue_xenium.v3.roi_catalog import (
    build_roi_catalog, extract_roi_selection, validate_roi_selection, write_roi_catalog,
)


def bundle(*, all_unknown=False, polygon=True, label="supplied type"):
    context = dict(schema_version="3.0", tissue="supplied tissue", disease="unspecified", organism="human",
        assay="xenium", time_step_hours=.1, recognition=False, perturbations=[],
        sampling=dict(design="subsampled", description="User-supplied labels", missing_populations=[]),
        scoring=dict(depth_bins=2, min_reference_cells=2, min_transcripts=20, min_overlap=1))
    cells = [dict(cell_id=f"{sample}-{i}", sample_id=sample, cell_type=label if i < 3 else "rare" if i == 3 else "unknown",
                  lineage="unknown" if all_unknown or i == 4 else "fibroblast",
                  x=float(i)*20, y=0. if sample == "a" else 50.,
                  expression={"ACTB": 5. if i == 4 else 30.+i, "MKI67": float(i)})
             for sample in ("a", "b") for i in range(5)]
    if polygon:
        cells[0]["polygon"] = [[-4, -4], [4, -4], [0, 4]] if polygon is True else [[0, 0], [1, 1]]
    return prepare(cells, context, assayed_genes=["ACTB", "MKI67"], gene_sets={"cell_cycle": ["MKI67"]},
                   annotation_provenance={"origin": "user_supplied", "description": "No independent review claimed"})


def test_catalogue_is_deterministic_does_not_select_and_accepts_unknown_labels():
    source = bundle(all_unknown=True)
    before = deepcopy(source)
    first = build_roi_catalog(source, cells_per_roi=3, grid_um=25, max_candidates=2)
    assert first == build_roi_catalog(source, cells_per_roi=3, grid_um=25, max_candidates=2)
    assert source == before
    assert first["selected_candidate"] is None
    assert len(first["candidates"]) == 4
    assert first["samples"]["a"]["summary"]["unknown_fraction"] == 1
    for candidate in first["candidates"]:
        assert len(candidate["selection"]["cell_ids"]) == 3
        assert all(cid.startswith(candidate["sample_id"]) for cid in candidate["selection"]["cell_ids"])
        assert candidate["summary"]["by_lineage"] == {"unknown": {"count": 3, "fraction": 1.}}
        assert candidate["selection"]["window"]["is_statistical_observation_window"] is False


def test_coverage_qc_geometry_and_small_sample_are_explicit():
    source = bundle()
    catalog = build_roi_catalog(source, cells_per_roi=1000, grid_um=100, max_candidates=1)
    candidate = catalog["candidates"][0]
    summary = candidate["summary"]
    assert candidate["undersized"]
    assert summary["n_units"] == 5
    assert summary["qc_pass"] == 4
    assert summary["qc_fail"] == 1
    assert summary["unknown_fraction"] == .2
    assert summary["polygons"] == 1 and summary["centroid_fallback"] == 4
    assert summary["transcript_total"] == sum(sum(c["expression"].values()) for c in source["cells"][:5])
    assert summary["source_score_coverage"]["cell_cycle"]["statuses"]["qc_fail_below_min_transcripts"] == 1
    rare = next(row for row in summary["potential_references"] if row["cell_type"] == "rare")
    assert rare["n_qc_pass"] == 1
    assert rare["cells_in_strata_meeting_minimum"] == 0
    assert "subsetting" in summary["reference_note"]
    assert catalog["cells"][0]["geometry"] == "supplied_polygon"
    assert catalog["cells"][1]["geometry"] == "centroid_fallback"


def test_selection_roundtrip_returns_copies_without_changing_source():
    source = bundle()
    selection = build_roi_catalog(source, cells_per_roi=3, grid_um=25)["candidates"][1]["selection"]
    checked = validate_roi_selection(source, selection)
    assert checked == selection and checked is not selection
    cells = extract_roi_selection(source, selection)
    assert [c["cell_id"] for c in cells] == selection["cell_ids"]
    cells[0]["expression"]["ACTB"] = 999
    assert all(c["expression"]["ACTB"] != 999 for c in source["cells"])


@pytest.mark.parametrize("field,value", [
    ("data_sha256", "wrong"), ("source_fingerprint", "wrong"),
    ("cell_ids", ["a-0", "a-0"]), ("cell_ids", ["a-0", "a-1", "missing"]),
    ("sample_id", "b"), ("selection_declared_by", "automatic"),
    ("window", {"kind": "rectangle", "bounds_um": [0, 0, 1, 1]}),
])
def test_selection_rejects_changed_source_membership_and_window(field, value):
    source = bundle()
    selection = deepcopy(build_roi_catalog(source, cells_per_roi=3, grid_um=25)["candidates"][0]["selection"])
    selection[field] = value
    with pytest.raises(ValueError):
        validate_roi_selection(source, selection)


def test_bad_polygon_fails_instead_of_silently_drawing_a_centroid():
    source = bundle()
    source["cells"][0]["polygon"] = [[0, 0], [1, 1]]
    source["fingerprint"] = digest({k: v for k, v in source.items() if k != "fingerprint"})
    with pytest.raises(ValueError, match="[Pp]olygon|geometry"):
        build_roi_catalog(source)


def test_invalid_supplied_polygon_is_visible_as_a_centroid_fallback():
    source = bundle(polygon="invalid")
    catalog = build_roi_catalog(source)
    assert catalog["cells"][0]["geometry"] == "invalid_polygon_centroid_fallback"
    assert catalog["cells"][0]["polygon"] is None
    assert catalog["samples"]["a"]["summary"]["invalid_polygons"] == 1
    assert catalog["samples"]["a"]["summary"]["polygons"] == 0


def test_raw_catalog_does_not_prepare_or_compute_expression_scores(monkeypatch):
    source = bundle()
    raw = {"cells": deepcopy(source["cells"]), "assayed_genes": ["ACTB", "MKI67"],
           "description": "Supplied data with no histology or independent review"}
    original = deepcopy(raw)
    from virtualtissue_xenium.v2 import guidebook
    monkeypatch.setattr(guidebook, "build_guidebook", lambda *a, **k: pytest.fail("Raw catalogue must not score the whole section"))
    catalog = build_roi_catalog(raw, cells_per_roi=3, grid_um=25, max_candidates=2)
    assert catalog["source_kind"] == "raw_json"
    assert catalog["source_fingerprint"] == digest(raw)
    assert catalog["score_coverage_status"] == "not_computed"
    assert catalog["candidates"][0]["summary"]["source_score_coverage"] == {}
    assert "diagnostic" in catalog["candidates"][0]["summary"]["qc_interpretation"]
    selection = catalog["candidates"][1]["selection"]
    selected = extract_roi_selection(raw, selection)
    assert [c["cell_id"] for c in selected] == selection["cell_ids"]
    assert raw == original
    raw["description"] = "Changed source description"
    with pytest.raises(ValueError, match="fingerprint"):
        validate_roi_selection(raw, selection)


def test_weak_panel_evidence_is_available_but_not_counted_as_empirical_reference():
    source = bundle()
    context = deepcopy(source["context"])
    context["scoring"]["method"] = "tiered_panel_auc"
    context["scoring"]["reference_fallback"] = "cell_type"
    tiered = prepare(source["cells"], context, assayed_genes=["ACTB", "MKI67"],
                     gene_sets={"cell_cycle": ["MKI67"]})
    catalog = build_roi_catalog(tiered, cells_per_roi=5, max_candidates=1)
    coverage = catalog["candidates"][0]["summary"]["source_score_coverage"]["cell_cycle"]
    assert coverage["weak_panel_available"] > 0
    assert coverage["available"] == coverage["empirical_available"] + coverage["weak_panel_available"]
    assert coverage["source_reference_size_min"] == 0
    assert catalog["candidates"][0]["selection"]["diagnostics"]["source_score_coverage"]["cell_cycle"] == coverage


def test_raw_invalid_nonfinite_polygon_has_stable_hash_and_visible_fallback():
    source = bundle(polygon=False)
    raw = {"cells": deepcopy(source["cells"]), "assayed_genes": ["ACTB", "MKI67"]}
    raw["cells"][0]["polygon"] = [[0, 0], [float("nan"), 1], [1, 0]]
    catalog = build_roi_catalog(raw, cells_per_roi=3)
    json.dumps(catalog, allow_nan=False)
    assert catalog["cells"][0]["geometry"] == "invalid_polygon_centroid_fallback"
    selection = catalog["candidates"][0]["selection"]
    assert validate_roi_selection(raw, selection) == selection
    assert extract_roi_selection(raw, selection)[0]["t0_geometry"]["invalid_reason"] == "nonfinite_or_invalid_xy_vertices"
    raw["description"] = float("nan")
    with pytest.raises(ValueError, match="metadata"):
        build_roi_catalog(raw)


def test_html_is_standalone_escapes_labels_and_never_overwrites(tmp_path):
    label = "</script><script>throw Error('injected')</script>"
    result = write_roi_catalog(bundle(label=label), tmp_path/"catalog.html", cells_per_roi=3, grid_um=25)
    html = (tmp_path/"catalog.html").read_text()
    assert result["selected_candidate"] is None
    assert label not in html
    assert "\\u003c/script>" in html
    assert "fetch(" not in html and '<script src=' not in html
    embedded = re.search(r'<script id="catalog-data" type="application/json">(.*?)</script>', html, re.S).group(1)
    assert json.loads(embedded) == json.loads((tmp_path/"catalog.catalog.json").read_text())
    assert json.loads(embedded)["cells"][0]["cell_type"] == label
    with pytest.raises(ValueError, match="exists"):
        write_roi_catalog(bundle(), tmp_path/"catalog.html")


def test_actual_html_player_requires_choice_and_exports_valid_descriptor(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for the actual offline player's DOM test")
    source = bundle()
    write_roi_catalog(source, tmp_path/"catalog.html", cells_per_roi=3, grid_um=25, max_candidates=2)
    html = (tmp_path/"catalog.html").read_text()
    data = re.search(r'<script id="catalog-data" type="application/json">(.*?)</script>', html, re.S).group(1)
    script = re.search(r'<script>\s*(.*?)</script>', html, re.S).group(1)
    harness = r'''
const vm=require('node:vm');const fs=require('node:fs');const assert=require('node:assert/strict');
const data=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));
class Element{constructor(tag='div'){this.tag=tag;this.children=[];this.style={};this.dataset={};this.value='';this.checked=false;this.disabled=false;this.classList={toggle:()=>{}}}append(...items){this.children.push(...items);if(this.tag==='select'&&!this.value&&items[0])this.value=items[0].value}replaceChildren(...items){this.children=[];this.append(...items)}click(){if(this.onclick)this.onclick()}getBoundingClientRect(){return {width:800,height:530,left:0,top:0}}getContext(){return new Proxy({}, {get:(target,key)=>target[key]||(()=>{})})}}
const elements={};for(const id of ['sample','download','map','catalog-data','legend','metrics','composition','coverage','references','provenance','detail-title','evidence-status','selection-status','choices','overview','zoom','polygons','clear','pointer'])elements[id]=new Element(id==='sample'?'select':'div');elements['catalog-data'].textContent=JSON.stringify(data);elements.polygons.checked=true;
let exported=null;const document={getElementById:id=>elements[id],createElement:tag=>new Element(tag),createTextNode:text=>({textContent:text})};
const sandbox={document,window:{devicePixelRatio:1,addEventListener:()=>{}},Blob:class{constructor(parts){exported=parts.join('')}},URL:{createObjectURL:()=>"blob:local",revokeObjectURL:()=>{}},setTimeout:fn=>fn(),console};
vm.runInNewContext(fs.readFileSync(process.argv[3],'utf8'),sandbox);
assert.equal(elements.download.disabled,true);elements.download.onclick();assert.equal(exported,null);
assert.equal(elements.choices.children.length,2);elements.choices.children[1].onclick();assert.equal(elements.download.disabled,false);elements.download.onclick();
const selection=JSON.parse(exported);assert.deepEqual(selection,data.candidates[1].selection);
elements.zoom.checked=true;elements.zoom.onchange();elements.polygons.checked=false;elements.polygons.onchange();
elements.sample.value='b';elements.sample.onchange();assert.equal(elements.download.disabled,true);
elements.choices.children[0].onclick();elements.clear.onclick();assert.equal(elements.download.disabled,true);
fs.writeFileSync(process.argv[4],JSON.stringify(selection));
'''
    (tmp_path/"catalog-data.json").write_text(data)
    (tmp_path/"player.js").write_text(script)
    (tmp_path/"harness.cjs").write_text(harness)
    subprocess.run([node, str(tmp_path/"harness.cjs"), str(tmp_path/"catalog-data.json"),
                    str(tmp_path/"player.js"), str(tmp_path/"selection.json")], check=True, capture_output=True, text=True)
    assert validate_roi_selection(source, tmp_path/"selection.json")["cell_ids"]
