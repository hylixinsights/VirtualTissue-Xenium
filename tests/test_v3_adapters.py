"""Native-format fixtures catch ID/coordinate/count and purity errors."""
import json
import numpy as np
import pandas as pd
import h5py
import pytest
from virtualtissue_xenium.v3.adapters import atera, visium_hd, stereo_seq
from virtualtissue_xenium.v2.inputs import prepare
from test_v3_contracts import context


def annotations(path, text="cell_id,cell_type,lineage,purity\na,tumour,malignant,0.9\nb,tumour,malignant,0.4\n"):
    path.write_text(text)
    return path


def test_visium_xy_scale_exact_matrix_join_and_purity(tmp_path):
    spatial=tmp_path/'spatial'; spatial.mkdir()
    (spatial/'scalefactors_json.json').write_text(json.dumps(dict(microns_per_pixel=.5,bin_size_um=8)))
    pd.DataFrame(dict(barcode=['a','b'],pxl_col_in_fullres=[0.,16.],pxl_row_in_fullres=[0.,0.])).to_parquet(spatial/'tissue_positions.parquet')
    with h5py.File(tmp_path/'filtered_feature_bc_matrix.h5','w') as f:
        m=f.create_group('matrix'); ft=m.create_group('features')
        m['barcodes']=np.array(['b','a'],dtype='S'); ft['name']=np.array(['MKI67','Control'],dtype='S')
        ft['feature_type']=np.array(['Gene Expression','Negative Control Probe'],dtype='S')
        m['shape']=[2,2]; m['indptr']=[0,2,3]; m['indices']=[0,1,0];m['data']=[30,999,20]
    ann=annotations(tmp_path/'annotation.csv')
    cells,genes,measurement=visium_hd.load(tmp_path,ann,measurement={'deconvolution_method':'declared test'})
    assert genes==['MKI67']
    assert [(c['cell_id'],c['x'],c['expression']) for c in cells]==[('a',0.,{'MKI67':20.}),('b',8.,{'MKI67':30.})]
    bundle=prepare(cells,context(assay='visium_hd',kernel='lattice'),assayed_genes=genes,measurement=measurement)
    assert [c['lineage'] for c in bundle['cells']]==['malignant','unknown']
    assert bundle['measurement']['is_single_cell_resolution'] is False
    (spatial/'scalefactors_json.json').write_text('{"bin_size_um":8}')
    with pytest.raises(ValueError,match='microns_per_pixel'):
        visium_hd.load(tmp_path,ann)


def test_atera_reads_csr_without_using_vendor_annotations(tmp_path):
    import zarr
    from zarr.storage import ZipStore
    with ZipStore(str(tmp_path/'cell_feature_matrix.zarr.zip'),mode='w') as store:
        v3 = int(zarr.__version__.split('.')[0]) >= 3
        g=zarr.open_group(store=store,mode='w',**({'zarr_format':2} if v3 else {}))
        for name,values in {'obs/barcode':['b','a'], 'obs/centroid_column':[10.,0.], 'obs/centroid_row':[0.,0.],
                            'var/feature_name':['MKI67','Control'], 'var/feature_type/categories':['Gene Expression','Control'],
                            'var/feature_type/codes':[0,1], 'X/indptr':[0,2,3], 'X/indices':[0,1,0], 'X/data':[30,999,20]}.items():
            (g.create_array if v3 else g.create_dataset)(name,data=np.asarray(values))
        g['X'].attrs.update({'encoding-type':'csr_matrix','shape':[2,2]})
    cells,genes,m=atera.load(tmp_path,annotations(tmp_path/'ann.csv'))
    assert genes==['MKI67'] and cells[0]['expression']=={'MKI67':20.}
    assert cells[1]['x']==10 and m['platform']=='atera'
    with pytest.raises(ValueError,match='reviewed annotations'):
        atera.load(tmp_path,None)


def test_gem_sums_counts_and_requires_explicit_scale_and_cell_centroids(tmp_path):
    gem=tmp_path/'data.gem'
    gem.write_text('#BinType=Cell\ngeneID\tx\ty\tMIDCount\tCellID\nMki67\t2\t4\t2\ta\nMki67\t3\t4\t3\ta\nActb\t5\t4\t4\tb\n')
    ann=annotations(tmp_path/'ann.csv','cell_id,cell_type,lineage,x,y\na,tumour,malignant,1,2\nb,tumour,malignant,2.5,2\n')
    with pytest.raises(ValueError,match='microns_per_coordinate'):
        stereo_seq.load(gem,ann,measurement={'unit':'segmented_cell'})
    cells,genes,m=stereo_seq.load(gem,ann,measurement={'unit':'segmented_cell','microns_per_coordinate':.5,'assayed_genes':['Actb','Mki67']})
    assert cells[0]['expression']=={'Mki67':5.} and cells[0]['x']==1.
    assert genes==['Actb','Mki67']
    ann=annotations(tmp_path/'no_xy.csv','cell_id,cell_type,lineage\na,tumour,malignant\nb,tumour,malignant\n')
    with pytest.raises(ValueError,match='centroids'):
        stereo_seq.load(gem,ann,measurement={'unit':'segmented_cell','microns_per_coordinate':.5,'assayed_genes':['Actb','Mki67']})
