"""Supplied polygons remain immutable evidence; modelled shapes stay explicit."""
from copy import deepcopy
import json
import math
from pathlib import Path
import re
import shutil
import subprocess

import pytest

from virtualtissue_xenium.v3.geometry import (normalize_geometry, projected_polygon, distance_between,
    validate_geometry_records, validate_geometry_config)
from virtualtissue_xenium.v2.inputs import prepare, validate_prepared
from virtualtissue_xenium.v2.controls import apply_control
from virtualtissue_xenium.v2.simulation import initialize, interaction_graph, decision_rows, commit, run
from virtualtissue_xenium.v2.replay import write_replay_stream


def square(x, y=0, radius=5):
    return [[x-radius,y-radius],[x+radius,y-radius],[x+radius,y+radius],[x-radius,y+radius]]


def cell(cid="a", x=0, polygon=None, lineage="malignant", **kw):
    result=dict(cell_id=cid,cell_type=lineage,lineage=lineage,expression={"ACTB":100},sample_id="s",x=x,y=0)
    if polygon is not None:
        result["polygon"]=polygon
    result.update(kw)
    return result


def context(metric="boundary", **extra):
    result=dict(schema_version="3.0",organism="human",tissue="skin",disease="synthetic geometry test",assay="xenium",
        sampling=dict(design="subsampled",description="Provided polygon engineering fixture",missing_populations=[]),
        recognition=True,perturbations=[],time_step_hours=.5,neighbor_radius_um=6,geometry=dict(interaction_distance=metric))
    result.update(extra)
    return result


def test_polygon_vertices_area_and_source_are_preserved_without_mutation():
    cells=[cell(polygon=square(0))]
    original=deepcopy(cells)
    result,metadata=normalize_geometry(cells,source="object.uns.segmentation")
    record=result[0]["t0_geometry"]
    assert cells==original and "polygon" not in result[0]
    assert record["polygon_um"]==square(0)+[square(0)[0]]
    assert record["area_um2"]==100 and record["centroid_um"]==[0,0]
    assert json.loads(record["raw_input_json"])==original[0]["polygon"]
    assert metadata["valid_polygons"]==1 and metadata["derived_centroids"]==0
    assert normalize_geometry(result)==(result,metadata)


def test_valid_polygon_can_supply_a_missing_centroid_but_does_not_replace_an_existing_one():
    c=cell(polygon=square(10))
    c.pop("x");c.pop("y")
    result,meta=normalize_geometry([c])
    assert (result[0]["x"],result[0]["y"])==(10,0) and meta["derived_centroids"]==1
    result,meta=normalize_geometry([cell(x=100,polygon=square(10))])
    assert result[0]["x"]==100 and result[0]["t0_geometry"]["status"]=="valid"
    assert meta["qc_warnings"]["provided_centroid_outside_polygon"]==1


@pytest.mark.parametrize("polygon,reason", [
    ([[0,0],[2,2],[0,2],[2,0]],"self_intersection"),
    ([[0,0],[float('nan'),2],[2,0]],"nonfinite_or_invalid_xy_vertices"),
    ([[0,0],[1,0],[2,0]],"zero_or_nonfinite_area"),
    ({"type":"Polygon","coordinates":[square(0),square(0,radius=1)]},"holes_or_multiple_rings_not_supported"),
    ({"type":"MultiPolygon","coordinates":[[square(0)]]},"unsupported_geometry_type")])
def test_invalid_polygons_are_retained_as_inert_source_without_fabricated_contours(polygon,reason):
    result,meta=normalize_geometry([cell(polygon=polygon)])
    record=result[0]["t0_geometry"]
    assert record["status"]=="invalid_polygon_centroid_fallback" and record["invalid_reason"]==reason
    assert record["polygon_um"] is None and projected_polygon(result[0]) is None
    assert meta["invalid_polygon_centroid_fallback"]==1
    json.dumps(result,allow_nan=False)
    assert normalize_geometry(result)[0]==result


def test_declared_scale_applies_to_raw_vertices_and_qc_tampering_is_rejected():
    result,_=normalize_geometry([cell(polygon={"type":"Polygon","coordinates":[square(0)]})],coordinate_units="pixels",coordinate_scale=2)
    assert result[0]["t0_geometry"]["area_um2"]==400
    altered=deepcopy(result);altered[0]["t0_geometry"]["area_um2"]=100
    with pytest.raises(ValueError,match="area/centroid"):
        validate_geometry_records(altered)
    with pytest.raises(ValueError,match="interaction_distance"):
        validate_geometry_config({"interaction_distance":"inferred_contact"})


def test_boundary_graph_finds_large_polygons_that_centroid_grid_would_miss():
    cells=[cell("a",0,square(0,radius=10)),cell("b",25,square(25,radius=10))]
    boundary=prepare(cells,context(),assayed_genes=["ACTB"])
    centroid=prepare(cells,context("centroid"),assayed_genes=["ACTB"])
    assert interaction_graph(initialize(boundary),boundary["context"])=={"a":["b"],"b":["a"]}
    assert interaction_graph(initialize(centroid),centroid["context"])=={"a":[],"b":[]}
    assert distance_between(*initialize(boundary).values())==(5,"boundary")
    boundary["cells"][1]["sample_id"]="other"
    assert interaction_graph(initialize(boundary),boundary["context"])=={"a":[],"b":[]}


def test_overlap_is_zero_and_invalid_pair_uses_declared_centroid_fallback():
    normalized,_=normalize_geometry([cell("a",0,square(0,radius=10)),cell("b",0,square(0,radius=1))])
    assert distance_between(*normalized)==(0,"boundary")
    normalized,_=normalize_geometry([cell("a",0,square(0,radius=10)),cell("b",25,[[0,0],[1,0],[2,0]])])
    assert distance_between(*normalized)==(25,"centroid_fallback")


def test_targets_are_sorted_by_boundary_gap_not_centroid_when_declared():
    cells=[cell("cd8",0,square(0),lineage="T_CD8",initial_state={"activated":True}),
           cell("near_center",15,square(15,radius=1)),cell("near_outline",30,square(30,radius=24))]
    b=prepare(cells,context(neighbor_radius_um=10),assayed_genes=["ACTB"])
    rows,_=decision_rows(initialize(b),b["guidebook"],b["context"],seed=9,step=0)
    cd8=next(r for r in rows if r["cell_id"]=="cd8")
    assert cd8["targets"]==["near_outline","near_center"]
    assert cd8["target_rule"]=="nearest_boundary_with_centroid_fallback_seeded_tie"


def test_motion_and_daughters_translate_shape_without_altering_t0_source():
    b=prepare([cell("cd8",0,square(0),lineage="T_CD8",initial_state={"activated":True})],context(),assayed_genes=["ACTB"])
    state=initialize(b);original=deepcopy(state["cd8"]["t0_geometry"])
    state["cd8"].update(moving_until=1,heading_radians=0,activated_until=100)
    moved,events=commit(state,[],max_population=10,seed=3,step=0,context=b["context"])
    assert moved["cd8"]["t0_geometry"]==original
    assert projected_polygon(moved["cd8"])[0][0]==original["polygon_um"][0][0]+3
    divided,events=commit(state,[dict(cell_id="cd8",action="PROLIFERATE",targets=[],corpses=[])],max_population=10,seed=3,step=0,context=b["context"])
    daughter=next(c for cid,c in divided.items() if cid!="cd8")
    assert daughter["t0_geometry"]==original and daughter["geometry_model"]=="synthetic_inherited_rigid_translation"
    assert next(e for e in events if e["kind"]=="division")["daughter_shape"]=="synthetic_inherited_rigid_translation"


def test_coordinate_shuffle_preserves_source_shapes_and_translates_runtime_outlines():
    b=prepare([cell("a",0,square(0)),cell("b",50,square(50))],context(),assayed_genes=["ACTB"])
    controlled=apply_control(b,"coordinate_shuffle",seed=2)
    validate_prepared(controlled)
    original={c["cell_id"]:c for c in b["cells"]}
    for c in controlled["cells"]:
        old=original[c["cell_id"]]
        assert c["t0_geometry"]==old["t0_geometry"]
        assert projected_polygon(c)[0][0]-projected_polygon(old)[0][0]==c["x"]-old["x"]


def test_sidecar_geometry_is_shared_and_rehydrates_losslessly(tmp_path):
    normal,_=normalize_geometry([cell("cell",0,square(0))])
    frames=[]
    for i in range(3):
        c=dict(normal[0],x=i,status="alive",geometry_model="provided_t0_rigid_translation")
        frames.append(dict(step=i,time_hours=i*.5,cells={"cell":c},events=[]))
    source=tmp_path/"frames.jsonl";source.write_text(''.join(json.dumps(f)+'\n' for f in frames))
    write_replay_stream(tmp_path/"replay.html",source,dict(schema_version="3.0"))
    script=(tmp_path/"replay_frames/chunk_00000.js").read_text()
    args=json.loads('['+script.split('(',1)[1].rsplit(');',1)[0]+']')
    _,loaded,shapes=args
    assert len(shapes)==1 and script.count('"raw_input_json"')==1
    for f in loaded:
        for c in f["cells"].values():
            c["t0_geometry"]=shapes[c.pop("t0_geometry_ref")]
    assert loaded==frames
    metadata=json.loads((tmp_path/"replay_frames/manifest.json").read_text())["replay"]
    assert metadata["has_polygons"] and metadata["bounds"]==[-5,7,-5,5]


def test_polygon_player_rehydrates_draws_and_selects_supplied_outlines(tmp_path):
    node=shutil.which("node")
    if node is None: pytest.skip("Node is optional for player DOM tests.")
    normalized,_=normalize_geometry([cell("cell",0,square(0))])
    frames=[]
    for i in range(151):
        c=dict(normalized[0],status="alive",x=i,y=-i,geometry_model="provided_t0_rigid_translation")
        frames.append(dict(step=i,time_hours=i*.5,cells={"cell":c},events=[] if i==0 else [dict(kind="move",cell_id="cell",from_xy=[i-1,1-i],to_xy=[i,-i])]))
    source=tmp_path/"frames.jsonl";source.write_text(''.join(json.dumps(f)+'\n' for f in frames))
    report=dict(schema_version="3.0",kernel="cell_agent",preflight=dict(interaction_mode="spatial"),mode="local",completed=True,completed_steps=150,seed=1,events={},usage={})
    write_replay_stream(tmp_path/"replay.html",source,report)
    completed=subprocess.run([node,str(Path(__file__).with_name('test_v3_replay_loader.cjs')),str(tmp_path),'provided'],capture_output=True,text=True,timeout=15)
    assert completed.returncode==0,completed.stderr
    assert json.loads(completed.stdout)["polygonVertexCounts"]["5"]>0
