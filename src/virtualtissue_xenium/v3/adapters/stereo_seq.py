"""Streaming GEM loader; requires explicit coordinate scale and measured bin/cell semantics."""
from collections import defaultdict
from pathlib import Path
import csv
import gzip
import math
from .common import annotation_rows, finish, select_indices, template_for_records, bounded_coordinate_records

SOURCE='https://cdn-newfile.stomics.tech/Stereo-seq_Analysis_Workflow_File_Format_A4.2.pdf'


def _geometry_rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith('#')), delimiter='\t')
        yield reader.fieldnames or []
        yield from reader


def template(path, output, roi=None, *, measurement=None):
    """Stream bin geometry, or use explicit segmented-cell centroids in µm.

    The gene universe and counts are irrelevant to a blank annotation template.
    For segmented GEM, a transcript-weighted centroid is never substituted for a
    measured segmentation centroid: measurement.centroids_csv is required.
    """
    path, settings = Path(path), dict(measurement or {})
    if path.suffix not in ('.gem', '.gz', '.tsv'):
        raise ValueError('Templates read GEM/TSV(.gz); export GEF with the vendor tools first.')
    factor = settings.get('microns_per_coordinate')
    if type(factor) not in (int, float) or not math.isfinite(factor) or factor <= 0:
        raise ValueError('Declare measurement.microns_per_coordinate; coordinate scale is never inferred.')
    unit = settings.get('unit')
    if unit not in ('segmented_cell', 'binned_spot'):
        raise ValueError('Declare Stereo-seq measurement.unit as segmented_cell or binned_spot.')
    if unit == 'binned_spot':
        def positions():
            rows = _geometry_rows(path)
            header = next(rows)
            if not {'x', 'y'} <= set(header):
                raise ValueError('GEM bin templates require measured x/y columns.')
            for row in rows:
                try:
                    x, y = float(row['x']) * factor, float(row['y']) * factor
                except (ValueError, TypeError) as exc:
                    raise ValueError('GEM bin coordinates must be numeric.') from exc
                yield f"{row['x']}:{row['y']}", x, y
        records = bounded_coordinate_records(positions(), roi)
    else:
        centroids_path = settings.get('centroids_csv')
        if not isinstance(centroids_path, str) or not centroids_path:
            raise ValueError('Segmented GEM templates require measurement.centroids_csv with measured cell_id,x,y in micrometres; transcript centroids are not inferred.')
        ids, xy = [], []
        with Path(centroids_path).open(newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle)
            if not {'cell_id', 'x', 'y'} <= set(reader.fieldnames or []):
                raise ValueError('Measured centroid CSV requires cell_id,x,y in micrometres.')
            for row in reader:
                if not row.get('cell_id'):
                    raise ValueError('Measured centroid IDs must be nonempty.')
                try:
                    xy.append((float(row['x']), float(row['y'])))
                except (ValueError, TypeError) as exc:
                    raise ValueError('Measured centroids must be numeric micrometre coordinates.') from exc
                ids.append(row['cell_id'])
        selected = select_indices(ids, xy, **(roi or {}))
        requested, observed = {ids[i] for i in selected}, set()
        rows = _geometry_rows(path)
        if 'CellID' not in next(rows):
            raise ValueError('Segmented GEM templates require measured CellID.')
        for row in rows:
            if row['CellID'] in requested and row['CellID'] not in ('0', '-1', ''):
                observed.add(row['CellID'])
        if requested != observed:
            raise ValueError('Selected measured centroid IDs absent from GEM: ' + str(sorted(requested-observed)[:5]))
        records = [dict(cell_id=ids[i], x=float(xy[i][0]), y=float(xy[i][1])) for i in selected]
    return template_for_records(records, output, unit=unit)


def load(path, annotations, roi=None, *, lineage_map=None, measurement=None):
    path=Path(path); settings=dict(measurement or {})
    factor=settings.get('microns_per_coordinate')
    if type(factor) not in (int,float) or not math.isfinite(factor) or factor<=0:
        raise ValueError('Declare measurement.microns_per_coordinate from the Stereo-seq chip; no implicit 0.5 conversion.')
    unit=settings.get('unit')
    if unit not in ('segmented_cell','binned_spot'):
        raise ValueError('Declare Stereo-seq measurement.unit as segmented_cell or binned_spot.')
    declared_genes=settings.get('assayed_genes')
    if declared_genes is None and settings.get('feature_universe') != 'observed_GEM_union':
        raise ValueError('GEM omits undetected genes. Supply measurement.assayed_genes or explicitly declare feature_universe=observed_GEM_union for exploratory scoring.')
    if declared_genes is not None and (not isinstance(declared_genes,list) or not declared_genes or
            any(not isinstance(g,str) or not g.strip() for g in declared_genes) or len(set(declared_genes))!=len(declared_genes)):
        raise ValueError('measurement.assayed_genes must contain unique nonempty exact gene symbols.')
    if path.suffix not in ('.gem','.gz','.tsv'):
        raise ValueError('This adapter reads GEM/TSV(.gz). Export GEF with the vendor tools before import.')
    # Bound imported units by the annotation IDs; expression never materializes the whole GEM.
    with Path(annotations).open(newline='',encoding='utf-8-sig') as f:
        label_rows=list(csv.DictReader(f))
    if len(label_rows)>int((roi or {}).get('max_cells',1000)):
        raise ValueError('Subset the Stereo-seq annotation table to the intended bounded ROI first.')
    requested={r['cell_id'] for r in label_rows}
    expression=defaultdict(lambda:defaultdict(float)); positions={}; sums=defaultdict(lambda:[0.,0.,0.]); all_genes=set(); seen=set()
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8') as f:
        reader=csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t')
        if not {'geneID','x','y','MIDCount'} <= set(reader.fieldnames or []):
            raise ValueError('GEM requires geneID, x, y, MIDCount.')
        if unit=='segmented_cell' and 'CellID' not in reader.fieldnames:
            raise ValueError('Segmented GEM requires measured CellID.')
        for row in reader:
            gene=row.get('geneName') or row['geneID']
            if not gene:
                raise ValueError('GEM gene symbols must be nonempty.')
            all_genes.add(gene)
            x,y=float(row['x']),float(row['y']); count=float(row['MIDCount'])
            if any(not math.isfinite(v) for v in (x,y,count)) or count<0:
                raise ValueError('GEM coordinates/counts must be finite and counts nonnegative.')
            cid=row['CellID'] if unit=='segmented_cell' else f"{row['x']}:{row['y']}"
            if unit=='segmented_cell' and cid in ('0','-1',''):
                continue
            if cid not in requested:
                continue
            seen.add(cid); expression[cid][gene]+=count
            if unit=='binned_spot':
                positions[cid]=(x*factor,y*factor)
            elif count>0:
                acc=sums[cid];acc[0]+=x*count;acc[1]+=y*count;acc[2]+=count
    if requested-seen:
        raise ValueError('Annotation IDs absent from the GEM: '+str(sorted(requested-seen)[:5]))
    labels=annotation_rows(annotations,seen,lineage_map)
    if unit=='segmented_cell':
        # A transcript-weighted centroid is not a measured segmentation centroid.
        for cid,row in labels.items():
            if not row.get('x') or not row.get('y'):
                raise ValueError('Segmented GEM needs reviewed annotation x/y centroids in micrometres; transcript centroids are not substituted.')
            positions[cid]=(float(row['x']),float(row['y']))
    ids=sorted(seen);xy=[positions[cid] for cid in ids];selected=select_indices(ids,xy,**(roi or {}))
    settings.update(platform='stereo_seq',panel_type='whole_transcriptome',is_single_cell_resolution=unit=='segmented_cell',
                    segmentation_method=settings.get('segmentation_method','declared GEM CellID' if unit=='segmented_cell' else 'none; declared GEM bins'),
                    adapter='stereo-gem-1',format_source=SOURCE)
    if declared_genes is not None:
        if all_genes-set(declared_genes):
            raise ValueError('GEM features are outside the declared gene universe; resolve gene IDs/symbols explicitly.')
        genes=sorted(declared_genes)
        settings['feature_universe']='declared_assayed_genes'
    else:
        genes=sorted(all_genes)
        settings['coverage_warning']='Observed GEM feature union is not the full assayed transcriptome; cross-platform AUC comparisons require a declared complete universe.'
    settings['gene_identifier_policy']='Exact geneName if present, otherwise exact geneID; no ID-to-symbol inference.'
    return finish(ids,xy,selected,[dict(expression[ids[i]]) for i in selected],genes,labels,sample_id=path.stem,metadata=settings)
