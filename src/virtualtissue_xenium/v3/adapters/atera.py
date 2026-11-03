"""Atera v1 AnnData Zarr CSR adapter, tested on format fixtures only."""
from pathlib import Path
from .common import annotation_rows, finish, select_indices, sparse_row, strings, template_for_records

SOURCE='https://www.10xgenomics.com/support/software/atera-onboard-analysis/latest/analysis/output-zarr-cfmatrix'


def _values(node):
    if hasattr(node, 'shape'):
        return node[:]
    categories=node['categories'][:]; codes=node['codes'][:]
    if any(int(c)<0 or int(c)>=len(categories) for c in codes):
        raise ValueError('Missing or invalid categorical codes in Atera metadata.')
    return [categories[int(c)] for c in codes]


def _geometry(root):
    ids = strings(_values(root['obs/barcode']))
    columns = _values(root['obs/centroid_column'])
    rows = _values(root['obs/centroid_row'])
    if len(columns) != len(rows):
        raise ValueError('Atera centroid row/column lengths differ.')
    return ids, list(zip(columns, rows))


def template(path, output, roi=None, *, measurement=None):
    """Export measured cell centroids; X/feature arrays are never accessed."""
    import zarr
    from zarr.storage import ZipStore
    if measurement is not None and measurement.get('unit', 'segmented_cell') != 'segmented_cell':
        raise ValueError('Atera native cell templates require unit=segmented_cell.')
    path = Path(path)
    file = path/'cell_feature_matrix.zarr.zip' if path.is_dir() else path
    store = ZipStore(str(file), mode='r')
    try:
        ids, xy = _geometry(zarr.open_group(store=store, mode='r'))
        selected = select_indices(ids, xy, **(roi or {}))
    finally:
        store.close()
    return template_for_records([dict(cell_id=ids[i], x=float(xy[i][0]), y=float(xy[i][1]))
                                 for i in selected], output, unit='segmented_cell')


def load(path, annotations, roi=None, *, lineage_map=None, measurement=None):
    import zarr
    from zarr.storage import ZipStore
    path=Path(path)
    file=path/'cell_feature_matrix.zarr.zip' if path.is_dir() else path
    store=ZipStore(str(file),mode='r')
    try:
        root=zarr.open_group(store=store,mode='r')
        ids, xy = _geometry(root)
        names=strings(_values(root['var/feature_name'])); kinds=strings(_values(root['var/feature_type']))
        if len(names)!=len(kinds):
            raise ValueError('Atera feature names/types mismatch.')
        settings=dict(measurement or {})
        gene_type=settings.get('gene_feature_type','Gene Expression')
        gi={i:g for i,(g,k) in enumerate(zip(names,kinds)) if k==gene_type}
        genes=list(gi.values())
        if not genes or len(set(genes))!=len(genes):
            raise ValueError('No unique gene features. Declare actual measurement.gene_feature_type from Atera categories.')
        selected=select_indices(ids,xy,**(roi or {})); labels=annotation_rows(annotations,set(ids),lineage_map)
        x=root['X']; pointers=x['indptr'][:]
        if x.attrs.get('encoding-type')!='csr_matrix' or tuple(x.attrs.get('shape',()))!=(len(ids),len(names)) or len(pointers)!=len(ids)+1 or pointers[0]!=0 or pointers[-1]!=x['data'].shape[0] or x['indices'].shape[0]!=x['data'].shape[0] or any(b<a for a,b in zip(pointers,pointers[1:])):
            raise ValueError('Expected validated cells-by-features CSR AnnData Zarr; CSC is not accepted.')
        rows=[]
        for i in selected:
            a,b=int(pointers[i]),int(pointers[i+1]); rows.append(sparse_row(x['indices'][a:b],x['data'][a:b],gi,len(names)))
    finally:
        store.close()
    settings.update(platform='atera',unit='segmented_cell',panel_type='whole_transcriptome',is_single_cell_resolution=True,
                    segmentation_method=settings.get('segmentation_method','Atera onboard; per-cell method not imported'),
                    adapter='atera-csr-zarr-1',format_source=SOURCE)
    return finish(ids,xy,selected,rows,genes,labels,sample_id=path.stem,metadata=settings)
