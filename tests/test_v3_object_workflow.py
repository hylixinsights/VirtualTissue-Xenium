"""User-provided geometry and existing labels travel through the actual CLI."""
from copy import deepcopy
import json

import pytest

from virtualtissue_xenium.v2.cli import main
from virtualtissue_xenium.v2.inputs import prepare, validate_prepared, from_anndata, read_json, write_json
from virtualtissue_xenium.v3.object_input import attach_polygons
from virtualtissue_xenium.v3.roi_catalog import build_roi_catalog
from test_v3_contracts import context, cells


def object_context():
    return context(time_step_hours=.25,
                   scoring=dict(method="tiered_panel_auc"),
                   geometry=dict(interaction_distance="boundary", invalid_polygon_policy="centroid_fallback"))


def population():
    values = cells(8)
    for cell in values:
        x,y=cell['x'],cell['y']
        cell['polygon']=[[x-4,y-3],[x+4,y-3],[x+4,y+3],[x-4,y+3]]
    return values


def test_raw_object_to_explicit_selection_to_prepared_cli(tmp_path):
    raw={'cells':population(), 'assayed_genes':list(population()[0]['expression']),
         'measurement': {'unit': 'nucleus', 'segmentation_method': 'supplied imperfect outlines'},
         'annotation_provenance': {'source': 'existing automated object labels'}}
    raw_before=deepcopy(raw)
    write_json(tmp_path/'object.json',raw)
    ctx = object_context()
    ctx['sampling'] = dict(design='enriched', description='immune enriched', missing_populations=['epithelial'])
    write_json(tmp_path/'context.json',ctx)
    assert main(['roi-catalog',str(tmp_path/'object.json'),str(tmp_path/'regions.html'),
                 '--cells-per-roi','4','--grid-um','20'])==0
    catalog=build_roi_catalog(raw,cells_per_roi=4,grid_um=20)
    assert raw==raw_before
    candidate=catalog['candidates'][0]
    selection=candidate['selection'] if 'selection' in candidate else candidate
    write_json(tmp_path/'selection.json',selection)
    assert main(['select-roi',str(tmp_path/'object.json'),str(tmp_path/'selection.json'),str(tmp_path/'selected.json'),
                 '--context',str(tmp_path/'context.json')])==0
    selected=read_json(tmp_path/'selected.json')
    validate_prepared(selected)
    assert len(selected['cells'])==4
    assert {c['cell_type'] for c in selected['cells']}=={'tumour'}
    assert selected['geometry']['valid_polygons']==4
    assert selected['annotation_provenance']['verification']=='not_independently_verified'
    assert selected['annotation_provenance']['source']=='existing automated object labels'
    assert selected['measurement']['unit']=='nucleus'
    assert selected['context']['sampling']['source_sampling']==ctx['sampling']
    assert 'epithelial' in selected['context']['sampling']['missing_populations']
    assert main(['expression-audit',str(tmp_path/'selected.json'),str(tmp_path/'audit.json')])==0
    assert (tmp_path/'audit.json').is_file()


def test_external_geojson_exact_join_and_unit_scale():
    raw={'type':'FeatureCollection','features':[
        {'type':'Feature','properties':{'cell_id':'a0'},'geometry':{'type':'Polygon','coordinates':[[[0,0],[20,0],[20,20],[0,20],[0,0]]]}}]}
    joined,meta=attach_polygons(cells(2),raw,coordinate_scale=.5)
    assert joined[0]['polygon']['coordinates'][0][1]==[10,0]
    assert meta['matched_cells']==1 and meta['cells_without_polygon']==1
    assert 'polygon' not in cells(2)[0]
    with pytest.raises(ValueError,match='No polygon IDs'):
        attach_polygons(cells(1),{'wrong':[[0,0],[1,0],[0,1]]})


def test_anndata_polygons_use_same_explicit_coordinate_conversion():
    ad=pytest.importorskip('anndata')
    np=pytest.importorskip('numpy')
    pd=pytest.importorskip('pandas')
    obj=ad.AnnData(np.array([[10.,20.]]),obs=pd.DataFrame({'cell_type':['existing'], 'lineage':['myeloid']},index=['c']),
                   var=pd.DataFrame(index=['A','B']))
    obj.obsm['spatial']=np.array([[10.,10.]])
    obj.uns['outlines']={'c':np.array([[0.,0.],[20.,0.],[20.,20.],[0.,20.]])}
    imported,genes=from_anndata(obj,polygon_key='outlines',coordinate_scale=.5)
    assert imported[0]['cell_type']=='existing'
    assert imported[0]['x']==5 and imported[0]['polygon'][1]==[10.,0.]
    assert genes==['A','B']


def test_bad_polygon_is_a_flag_not_a_fake_contour():
    pop=population()
    pop[0]['polygon']=[[0,0],[float('nan'),1],[1,0]]
    bundle=prepare(pop,object_context())
    validate_prepared(bundle)
    assert bundle['geometry']['invalid_polygon_centroid_fallback']==1
    assert bundle['cells'][0]['t0_geometry']['polygon_um'] is None
    json.dumps(bundle,allow_nan=False)


def test_geometry_alias_feature_scaling_and_duplicate_source():
    pop = cells(1)
    feature = {'type':'Feature', 'properties':{}, 'geometry':{'type':'Polygon',
               'coordinates': [[[0,0],[10,0],[0,10],[0,0]]]}}
    joined, _ = attach_polygons(pop, {'a0': feature}, coordinate_scale=.5)
    joined[0]['segmentation_polygon'] = joined[0].pop('polygon')
    bundle = prepare(joined, object_context())
    validate_prepared(bundle)
    assert bundle['cells'][0]['t0_geometry']['area_um2'] == 12.5
    for alias in ('geometry', 'segmentation_polygon'):
        with pytest.raises(ValueError, match='already has geometry'):
            attach_polygons([dict(pop[0], **{alias:feature})], {'a0': feature})


def test_raw_prepare_preserves_declared_metadata(tmp_path):
    raw = {'cells':population(), 'assayed_genes':list(population()[0]['expression']),
           'measurement':{'unit':'nucleus'}, 'annotation_provenance':{'source':'existing labels'}}
    write_json(tmp_path/'raw.json',raw)
    write_json(tmp_path/'context.json',object_context())
    assert main(['prepare', str(tmp_path/'raw.json'), str(tmp_path/'context.json'), str(tmp_path/'bundle.json')])==0
    bundle=read_json(tmp_path/'bundle.json')
    validate_prepared(bundle)
    assert bundle['measurement']['unit']=='nucleus'
    assert bundle['annotation_provenance']['source']=='existing labels'


def test_catalog_explicit_lineage_mapping_keeps_existing_labels(tmp_path):
    raw = {'cells':population()}
    for cell in raw['cells']:
        del cell['lineage']
    write_json(tmp_path/'raw.json',raw)
    write_json(tmp_path/'context.json',object_context())
    mapping=tmp_path/'map.csv'
    mapping.write_text('cell_type,lineage\ntumour,malignant\n')
    assert main(['roi-catalog',str(tmp_path/'raw.json'),str(tmp_path/'regions.html'),
                 '--cells-per-roi','4','--grid-um','20','--lineage-map',str(mapping)])==0
    catalog=read_json(tmp_path/'regions.catalog.json')
    write_json(tmp_path/'selection.json',catalog['candidates'][0]['selection'])
    assert main(['select-roi',str(tmp_path/'raw.json'),str(tmp_path/'selection.json'),str(tmp_path/'bundle.json'),
                 '--context',str(tmp_path/'context.json'),'--lineage-map',str(mapping)])==0
    bundle=read_json(tmp_path/'bundle.json')
    assert {c['cell_type'] for c in bundle['cells']}=={'tumour'}
    assert bundle['annotation_provenance']['lineage_map']=={'tumour':'malignant'}
    assert all('lineage' not in c for c in read_json(tmp_path/'raw.json')['cells'])
