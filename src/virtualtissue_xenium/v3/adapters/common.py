"""Shared bounded ROI selection, reviewed labels and sparse count validation."""
from pathlib import Path
import csv
import hashlib
import heapq
import math
from ...v2.inputs import validate_cells


def strings(values):
    return [v.decode() if isinstance(v, bytes) else str(v) for v in values]


def annotation_rows(path, ids, lineage_map=None):
    if path is None:
        raise ValueError("Provide reviewed annotations; clusters are not lineages.")
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        needed = {'cell_id', 'cell_type'} | (set() if lineage_map else {'lineage'})
        if not needed <= set(reader.fieldnames or []):
            raise ValueError(f"Annotations require {sorted(needed)}.")
        result = {}
        for row in reader:
            cid = row['cell_id']
            if cid in result or cid not in ids:
                raise ValueError("Annotation IDs must be unique and present in the measured matrix.")
            if lineage_map:
                row['lineage'] = lineage_map.get(row['cell_type'], '')
            if not row['cell_type'].strip() or not row['lineage'].strip():
                raise ValueError("Annotations must explicitly name cell_type and lineage (unknown is allowed).")
            result[cid] = row
    return result


def select_indices(ids, xy, *, center=None, radius_um=None, max_cells=1000):
    if type(max_cells) is not int or max_cells < 1:
        raise ValueError("max_cells must be a positive integer.")
    if len(ids) != len(xy) or len(set(ids)) != len(ids):
        raise ValueError("Coordinate rows and unique measured IDs must match exactly.")
    if any(len(pos) != 2 or any(not math.isfinite(float(v)) for v in pos) for pos in xy):
        raise ValueError("Coordinates must be finite measured XY in micrometres.")
    if center is None and (radius_um is not None or len(ids) > max_cells):
        raise ValueError("An explicit center is required for bounded ROI selection.")
    if center is not None and (len(center) != 2 or any(not math.isfinite(float(v)) for v in center)):
        raise ValueError("center must be two finite micrometre coordinates.")
    if radius_um is not None and (not math.isfinite(radius_um) or radius_um <= 0):
        raise ValueError("radius_um must be finite and positive.")
    distances = [((float(p[0])-center[0])**2+(float(p[1])-center[1])**2 if center else 0) for p in xy]
    chosen = sorted((i for i in range(len(ids)) if radius_um is None or distances[i] <= radius_um**2),
                    key=lambda i: (distances[i], ids[i]))[:max_cells]
    if not chosen:
        raise ValueError("The selected ROI contains no units.")
    return sorted(chosen, key=lambda i:ids[i])


def finish(ids, xy, selected, expressions, genes, annotations, *, sample_id, metadata):
    records = []
    for i, expression in zip(selected, expressions):
        cid = ids[i]
        if cid not in annotations:
            raise ValueError(f"Every selected unit needs an annotation: missing {cid}.")
        label = annotations[cid]
        cell = dict(cell_id=cid, cell_type=label['cell_type'].strip(), lineage=label['lineage'].strip(),
                    sample_id=sample_id, x=float(xy[i][0]), y=float(xy[i][1]), expression=expression)
        for field in ('purity', 'occupancy'):
            if label.get(field, '').strip():
                cell[field] = float(label[field])
        records.append(cell)
    return validate_cells(records), genes, metadata


def sparse_row(indices, values, gene_indices, total_features):
    if len(indices) != len(values) or len(set(map(int,indices))) != len(indices):
        raise ValueError("Malformed sparse row: duplicate indices or unequal lengths.")
    if any(int(i) < 0 or int(i) >= total_features for i in indices):
        raise ValueError("Sparse feature index out of bounds.")
    if any(not math.isfinite(float(v)) or v < 0 for v in values):
        raise ValueError("Counts must be finite and nonnegative.")
    return {gene_indices[int(i)]:float(v) for i,v in zip(indices,values) if int(i) in gene_indices and v != 0}


def template_for_records(records, output, *, unit):
    path=Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', newline='', encoding='utf-8') as f:
        writer=csv.DictWriter(f, fieldnames=['cell_id','x','y','measurement_unit','cell_type','lineage','purity'])
        writer.writeheader()
        writer.writerows(dict({k: c[k] for k in ('cell_id','x','y')}, measurement_unit=unit) for c in records)
    return path


class _FarthestFirst:
    """Heap root is the worst candidate, including deterministic cell-ID ties."""
    def __init__(self, key, record):
        self.key, self.record = key, record

    def __lt__(self, other):
        return self.key > other.key


def bounded_coordinate_records(records, roi=None):
    """Retain O(max_cells) positions from a repeated-ID coordinate stream.

    Used for GEM bins: their IDs derive directly from the coordinate strings,
    so a repeated ID has the same coordinates and needs no count accumulation.
    """
    settings = dict(roi or {})
    unknown = set(settings) - {'center', 'radius_um', 'max_cells'}
    if unknown:
        raise ValueError(f'Unknown ROI options: {sorted(unknown)}.')
    center, radius, limit = settings.get('center'), settings.get('radius_um'), settings.get('max_cells', 1000)
    if type(limit) is not int or limit < 1:
        raise ValueError('max_cells must be a positive integer.')
    if center is not None:
        if len(center) != 2 or any(not math.isfinite(float(v)) for v in center):
            raise ValueError('center must be two finite micrometre coordinates.')
        center = tuple(map(float, center))
    if radius is not None:
        if center is None:
            raise ValueError('An explicit center is required for bounded ROI selection.')
        if type(radius) not in (int, float) or not math.isfinite(radius) or radius <= 0:
            raise ValueError('radius_um must be finite and positive.')
    heap, selected = [], {}
    for cid, x, y in records:
        if not isinstance(cid, str) or not cid or not math.isfinite(x) or not math.isfinite(y):
            raise ValueError('Template coordinates and measured IDs must be finite and nonempty.')
        if cid in selected:
            old = selected[cid]
            if old['x'] != x or old['y'] != y:
                raise ValueError('Repeated measured ID has conflicting coordinates.')
            continue
        distance = (x-center[0])**2 + (y-center[1])**2 if center else 0.
        if radius is not None and distance > radius**2:
            continue
        record = dict(cell_id=cid, x=x, y=y)
        candidate = _FarthestFirst((distance, cid), record)
        if len(heap) < limit:
            heapq.heappush(heap, candidate)
            selected[cid] = record
        elif center is None:
            raise ValueError('An explicit center is required for bounded ROI selection.')
        elif candidate.key < heap[0].key:
            removed = heapq.heapreplace(heap, candidate)
            selected.pop(removed.record['cell_id'])
            selected[cid] = record
    if not selected:
        raise ValueError('The selected ROI contains no units.')
    return [selected[cid] for cid in sorted(selected)]
