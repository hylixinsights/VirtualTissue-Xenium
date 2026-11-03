"""Space Ranger HD binned matrix adapter; no invented segmentation or purity."""
from pathlib import Path
import json
import math
from .common import annotation_rows, finish, select_indices, sparse_row, strings, template_for_records

SOURCE = 'https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs'


def template(path, output, roi=None, *, measurement=None):
    """Read bin barcodes and calibrated positions without reading matrix counts."""
    import h5py
    import pandas as pd
    root = Path(path)
    settings = dict(measurement or {})
    if settings.get('unit', 'binned_spot') != 'binned_spot':
        raise ValueError('This Visium HD adapter exports fixed-bin templates, not inferred cells.')
    scales = json.loads((root/'spatial/scalefactors_json.json').read_text())
    factor = settings.get('microns_per_pixel', scales.get('microns_per_pixel'))
    size = settings.get('unit_size_um', scales.get('bin_size_um'))
    if type(factor) not in (int, float) or not math.isfinite(factor) or factor <= 0:
        raise ValueError('Visium HD requires explicit calibrated microns_per_pixel.')
    if type(size) not in (int, float) or not math.isfinite(size) or size <= 0:
        raise ValueError('Declare Visium HD bin_size_um / measurement.unit_size_um.')
    positions = pd.read_parquet(root/'spatial/tissue_positions.parquet',
        columns=['barcode', 'pxl_col_in_fullres', 'pxl_row_in_fullres'])
    if positions.barcode.duplicated().any():
        raise ValueError('Invalid/duplicate Visium HD tissue positions.')
    positions = positions.set_index('barcode')
    with h5py.File(root/'filtered_feature_bc_matrix.h5', 'r') as handle:
        ids = strings(handle['matrix/barcodes'][:])
    if not set(ids) <= set(positions.index):
        raise ValueError('Missing coordinates for matrix barcodes.')
    xy = [(float(positions.loc[cid, 'pxl_col_in_fullres']) * factor,
           float(positions.loc[cid, 'pxl_row_in_fullres']) * factor) for cid in ids]
    selected = select_indices(ids, xy, **(roi or {}))
    return template_for_records([dict(cell_id=ids[i], x=xy[i][0], y=xy[i][1]) for i in selected],
                                output, unit='binned_spot')


def load(path, annotations, roi=None, *, lineage_map=None, measurement=None):
    import h5py
    import pandas as pd
    root=Path(path)
    settings=dict(measurement or {})
    scales=json.loads((root/'spatial/scalefactors_json.json').read_text())
    factor=settings.get('microns_per_pixel', scales.get('microns_per_pixel'))
    if type(factor) not in (float,int) or factor <= 0:
        raise ValueError('Visium HD requires explicit calibrated microns_per_pixel; pixel positions are not micrometres.')
    size=settings.get('unit_size_um', scales.get('bin_size_um'))
    if type(size) not in (float,int) or size <= 0:
        raise ValueError('Declare Visium HD bin_size_um / measurement.unit_size_um.')
    obs=pd.read_parquet(root/'spatial/tissue_positions.parquet')
    needed={'barcode','pxl_col_in_fullres','pxl_row_in_fullres'}
    if not needed <= set(obs) or obs.barcode.duplicated().any():
        raise ValueError('Invalid/duplicate Visium HD tissue positions.')
    obs=obs.set_index('barcode')
    with h5py.File(root/'filtered_feature_bc_matrix.h5','r') as f:
        m=f['matrix']; ids=strings(m['barcodes'][:]); names=strings(m['features/name'][:])
        kinds=strings(m['features/feature_type'][:])
        if not set(ids) <= set(obs.index):
            raise ValueError('Missing coordinates for matrix barcodes.')
        if len(kinds)!=len(names):
            raise ValueError('Feature names/types mismatch.')
        gi={i:g for i,(g,k) in enumerate(zip(names,kinds)) if k=='Gene Expression'}
        genes=list(gi.values())
        if not genes or len(set(genes))!=len(genes):
            raise ValueError('Gene symbols must be present and unique.')
        xy=[(float(obs.loc[cid,'pxl_col_in_fullres'])*factor,float(obs.loc[cid,'pxl_row_in_fullres'])*factor) for cid in ids]
        selected=select_indices(ids,xy,**(roi or {})); labels=annotation_rows(annotations,set(ids),lineage_map)
        pointers=m['indptr'][:]
        if tuple(m['shape'][:])!=(len(names),len(ids)) or len(pointers)!=len(ids)+1 or pointers[0]!=0 or pointers[-1]!=len(m['data']) or len(m['indices'])!=len(m['data']) or any(b<a for a,b in zip(pointers,pointers[1:])):
            raise ValueError('Malformed genes-by-barcodes CSC matrix.')
        rows=[]
        for i in selected:
            a,b=int(pointers[i]),int(pointers[i+1]); rows.append(sparse_row(m['indices'][a:b],m['data'][a:b],gi,len(names)))
    settings.update(platform='visium_hd',unit='binned_spot',unit_size_um=float(size), panel_type='whole_transcriptome',
                    segmentation_method='none; fixed Space Ranger bins',is_single_cell_resolution=False,
                    microns_per_pixel=float(factor),adapter='visium_hd-binned-1',format_source=SOURCE)
    return finish(ids,xy,selected,rows,genes,labels,sample_id=root.name,metadata=settings)
