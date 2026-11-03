"""Preserve supplied segmentation and use it without claiming ground-truth shape.

Valid single-ring polygons remain immutable at t0. Runtime shapes are rigid
translations of those measurements; inherited daughter shapes are explicitly
synthetic. Invalid or unsupported polygons are recorded and use centroids.
No contour repair, inferred segmentation or histological identity is performed.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
import hashlib
import json
import math
from numbers import Real

EPSILON = 1e-9
DEFAULT_CONFIG = dict(interaction_distance="centroid", invalid_polygon_policy="centroid_fallback")
RAW_KEYS = ("polygon", "polygon_um", "segmentation_polygon", "geometry")


def validate_geometry_config(raw=None):
    if raw is None:
        raw = {}
    if not isinstance(raw, dict) or set(raw) - set(DEFAULT_CONFIG):
        raise ValueError(f"geometry accepts only {sorted(DEFAULT_CONFIG)}.")
    result = dict(DEFAULT_CONFIG, **raw)
    if result["interaction_distance"] not in {"centroid", "boundary"}:
        raise ValueError("geometry.interaction_distance must be centroid or boundary.")
    if result["invalid_polygon_policy"] != "centroid_fallback":
        raise ValueError("Invalid polygons currently require explicit centroid_fallback policy.")
    return result


def _json_input(value):
    if hasattr(value, "__geo_interface__"):
        return _json_input(value.__geo_interface__)
    if hasattr(value, "tolist"):
        return _json_input(value.tolist())
    if isinstance(value, dict):
        return {str(k): _json_input(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_input(v) for v in value]
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    return {"unsupported_input_type": type(value).__name__}


def _number(value):
    return isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)


def _bbox(points):
    return (min(p[0] for p in points), max(p[0] for p in points), min(p[1] for p in points), max(p[1] for p in points))


def bbox_gap(a, b):
    return math.hypot(max(0.0, a[0]-b[1], b[0]-a[1]), max(0.0, a[2]-b[3], b[2]-a[3]))


def _cross(a, b, c):
    return (b[0]-a[0]) * (c[1]-a[1]) - (b[1]-a[1]) * (c[0]-a[0])


def _on_segment(a, b, p):
    tolerance = EPSILON * max(1.0, math.hypot(b[0]-a[0], b[1]-a[1]))
    return (abs(_cross(a, b, p)) <= tolerance and min(a[0], b[0])-EPSILON <= p[0] <= max(a[0], b[0])+EPSILON
            and min(a[1], b[1])-EPSILON <= p[1] <= max(a[1], b[1])+EPSILON)


def _segments_intersect(a, b, c, d):
    if bbox_gap(_bbox((a,b)), _bbox((c,d))) > EPSILON:
        return False
    ab_c, ab_d, cd_a, cd_b = _cross(a,b,c), _cross(a,b,d), _cross(c,d,a), _cross(c,d,b)
    if ((ab_c > 0 > ab_d or ab_d > 0 > ab_c) and (cd_a > 0 > cd_b or cd_b > 0 > cd_a)):
        return True
    return any((_on_segment(a,b,c), _on_segment(a,b,d), _on_segment(c,d,a), _on_segment(c,d,b)))


def _self_intersects(ring):
    # Sweep segment bounding boxes instead of comparing every distant edge.
    edges = [(min(a[0],b[0]), max(a[0],b[0]), i, a, b) for i,(a,b) in enumerate(zip(ring,ring[1:]))]
    active = []
    count = len(edges)
    for edge in sorted(edges):
        lo, hi, i, a, b = edge
        active = [old for old in active if old[1]+EPSILON >= lo]
        for _, _, j, c, d in active:
            if abs(i-j) == 1 or {i,j} == {0,count-1}:
                continue
            if _segments_intersect(a,b,c,d):
                return True
        active.append(edge)
    return False


def point_in_polygon(point, ring):
    inside = False
    for a,b in zip(ring,ring[1:]):
        if _on_segment(a,b,point):
            return True
        if (a[1] > point[1]) != (b[1] > point[1]):
            crossing_x = a[0] + (point[1]-a[1]) * (b[0]-a[0]) / (b[1]-a[1])
            if point[0] < crossing_x:
                inside = not inside
    return inside


def _analyse(raw, scale):
    if isinstance(raw, dict) and raw.get("type") == "Feature":
        raw = raw.get("geometry")
    if isinstance(raw, dict):
        if raw.get("type") != "Polygon":
            return None, "unsupported_geometry_type"
        coordinates = raw.get("coordinates")
        if not isinstance(coordinates, list) or len(coordinates) != 1:
            return None, "holes_or_multiple_rings_not_supported"
        raw = coordinates[0]
    if not isinstance(raw, list) or len(raw) < 3:
        return None, "fewer_than_three_vertices"
    if any(not isinstance(p,list) or len(p) != 2 or not all(_number(v) for v in p) for p in raw):
        return None, "nonfinite_or_invalid_xy_vertices"
    ring = [[float(p[0])*scale, float(p[1])*scale] for p in raw]
    if any(not all(math.isfinite(v) for v in p) for p in ring):
        return None, "scaled_vertices_overflow"
    if ring[-1] != ring[0]:
        ring.append(list(ring[0]))
    if len({tuple(p) for p in ring[:-1]}) < 3:
        return None, "fewer_than_three_distinct_vertices"
    if any(a == b for a,b in zip(ring,ring[1:])):
        return None, "duplicate_adjacent_vertices"
    if _self_intersects(ring):
        return None, "self_intersection"
    origin = ring[0]
    local = [[p[0]-origin[0],p[1]-origin[1]] for p in ring]
    crosses = [a[0]*b[1]-b[0]*a[1] for a,b in zip(local,local[1:])]
    twice_area = sum(crosses)
    if not math.isfinite(twice_area) or abs(twice_area) <= EPSILON:
        return None, "zero_or_nonfinite_area"
    centroid = [origin[k] + sum((a[k]+b[k])*cross for a,b,cross in zip(local,local[1:],crosses)) / (3*twice_area) for k in (0,1)]
    if not all(math.isfinite(v) for v in centroid):
        return None, "nonfinite_centroid"
    return dict(polygon_um=ring, area_um2=abs(twice_area)/2, centroid_um=centroid,
                bbox_um=list(_bbox(ring)), vertex_count=len(ring)-1), None


def normalize_geometry(cells, *, source="provided_object", coordinate_units="micrometres", coordinate_scale=1.0):
    """Return copied cells and QC metadata, preserving raw input as inert JSON text.

    ``coordinate_scale`` acts only on raw polygon vertices: input x/y must already
    be in micrometres. Existing normalized t0 records retain their original anchor,
    allowing coordinate-shuffle controls to translate shapes by changing x/y only.
    When both x/y are absent, a valid polygon supplies its declared derived centroid.
    """
    if not _number(coordinate_scale) or coordinate_scale <= 0:
        raise ValueError("Polygon coordinate_scale must be finite and positive.")
    if not isinstance(source,str) or not source or not isinstance(coordinate_units,str) or not coordinate_units:
        raise ValueError("Geometry source and coordinate_units must be explicit nonempty text.")
    result = deepcopy(cells)
    for cell in result:
        keys = [key for key in RAW_KEYS if cell.get(key) is not None]
        if not keys:
            if "t0_geometry" in cell:
                validate_geometry_records([cell])
            continue
        raw = _json_input(cell[keys[0]]) if len(keys) == 1 else {key:_json_input(cell[key]) for key in keys}
        raw_text = json.dumps(raw,sort_keys=True,separators=(",",":"),ensure_ascii=True,allow_nan=True)
        record = dict(version="provided-polygon-1", source=source, source_cell_id=cell["cell_id"], input_key=keys[0] if len(keys)==1 else "ambiguous",
                      raw_input_json=raw_text, input_sha256=hashlib.sha256(raw_text.encode()).hexdigest(),
                      coordinate_units=coordinate_units, coordinate_scale=float(coordinate_scale), warnings=[])
        details, failure = _analyse(raw, coordinate_scale) if len(keys)==1 else (None,"ambiguous_polygon_fields")
        if ("x" in cell) != ("y" in cell):
            raise ValueError("Provide both centroid x/y or neither when deriving from a polygon.")
        if "x" not in cell and details is not None:
            cell["x"],cell["y"] = details["centroid_um"]
            anchor_source = "polygon_centroid"
        else:
            anchor_source = "provided_centroid" if "x" in cell else "unavailable"
        anchor = [cell["x"],cell["y"]] if "x" in cell else None
        if anchor is not None and not all(_number(v) for v in anchor):
            raise ValueError("Provided centroid coordinates must be finite; invalid centroids are not silently replaced.")
        record.update(anchor_um=anchor, anchor_source=anchor_source)
        if details is None:
            record.update(status="invalid_polygon_centroid_fallback", polygon_um=None, area_um2=None, centroid_um=None,
                          bbox_um=None, vertex_count=0, invalid_reason=failure)
        else:
            record.update(details,status="valid",invalid_reason=None)
            record["centroid_offset_um"] = math.dist(anchor, details["centroid_um"])
            if not point_in_polygon(anchor, details["polygon_um"]):
                record["warnings"].append("provided_centroid_outside_polygon")
            if record["vertex_count"] > 200:
                record["warnings"].append("many_vertices_boundary_distance_may_be_expensive")
        for key in keys:
            cell.pop(key)
        cell["t0_geometry"] = record
    return result, summarize_geometry(result)


def validate_geometry_records(cells):
    for cell in cells:
        record = cell.get("t0_geometry")
        if record is None:
            continue
        if not isinstance(record,dict) or record.get("version") != "provided-polygon-1":
            raise ValueError("Unknown t0 geometry record; normalize supplied polygons first.")
        raw_text = record.get("raw_input_json")
        if not isinstance(raw_text,str) or hashlib.sha256(raw_text.encode()).hexdigest() != record.get("input_sha256"):
            raise ValueError("The preserved polygon input hash does not match its raw record.")
        scale = record.get("coordinate_scale")
        if not _number(scale) or scale <= 0:
            raise ValueError("Invalid preserved polygon coordinate scale.")
        details,failure = _analyse(json.loads(raw_text),scale) if record.get("input_key") != "ambiguous" else (None,"ambiguous_polygon_fields")
        expected_status = "valid" if details is not None else "invalid_polygon_centroid_fallback"
        if record.get("status") != expected_status or record.get("invalid_reason") != failure:
            raise ValueError("Preserved polygon QC status does not match its source geometry.")
        if details is not None:
            if any(record.get(k) != v for k,v in details.items()):
                raise ValueError("Normalized polygon vertices/area/centroid differ from preserved source geometry.")
            anchor = record.get("anchor_um")
            if not isinstance(anchor,list) or len(anchor)!=2 or not all(_number(v) for v in anchor):
                raise ValueError("Valid polygons require an immutable finite t0 anchor.")
        elif record.get("polygon_um") is not None:
            raise ValueError("Invalid polygons cannot supply a fabricated usable contour.")
    return cells


def summarize_geometry(cells):
    records = [c["t0_geometry"] for c in cells if "t0_geometry" in c]
    valid = [r for r in records if r.get("status") == "valid"]
    return dict(version="provided-polygon-1", cells=len(cells), supplied_polygons=len(records), valid_polygons=len(valid),
                invalid_polygon_centroid_fallback=len(records)-len(valid), centroid_only=len(cells)-len(records),
                derived_centroids=sum(r.get("anchor_source")=="polygon_centroid" for r in records),
                qc_warnings=dict(Counter(w for r in records for w in r.get("warnings",[]))),
                total_vertices=sum(r.get("vertex_count",0) for r in valid),
                sources=sorted({r.get("source","") for r in records}),
                source_hashes={c["cell_id"]:c["t0_geometry"]["input_sha256"] for c in cells if "t0_geometry" in c},
                interpretation="Supplied segmentation retained with QC, not a validated true contour; shape does not establish cell identity.",
                runtime_model="rigid translation of t0 shape; daughters inherit a synthetic copy, without shape mechanics",
                distance_cost="bounding-box candidate sweep; exact polygon pairs worst-case O(vertices_a*vertices_b); centroid fallback if either polygon is invalid or absent")


def projected_polygon(cell):
    geometry = cell.get("t0_geometry")
    if not geometry or geometry.get("status") != "valid" or "x" not in cell or "y" not in cell:
        return None
    anchor = geometry["anchor_um"]
    dx,dy = cell["x"]-anchor[0],cell["y"]-anchor[1]
    return [[p[0]+dx,p[1]+dy] for p in geometry["polygon_um"]]


def _point_segment_distance(p,a,b):
    dx,dy = b[0]-a[0],b[1]-a[1]
    norm = dx*dx+dy*dy
    fraction = max(0.0,min(1.0,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/norm)) if norm else 0.0
    return math.hypot(p[0]-a[0]-fraction*dx,p[1]-a[1]-fraction*dy)


def polygon_gap(a,b):
    """Minimum boundary gap in µm; intersections and containment have gap zero."""
    if point_in_polygon(a[0],b) or point_in_polygon(b[0],a):
        return 0.0
    best = math.dist(a[0],b[0])
    for p,q in zip(a,a[1:]):
        box_a = _bbox((p,q))
        for r,s in zip(b,b[1:]):
            if bbox_gap(box_a,_bbox((r,s))) > best:
                continue
            if _segments_intersect(p,q,r,s):
                return 0.0
            best = min(best,_point_segment_distance(p,r,s),_point_segment_distance(q,r,s),
                       _point_segment_distance(r,p,q),_point_segment_distance(s,p,q))
    return best


def distance_between(a,b, *, metric="boundary"):
    if metric not in {"boundary","centroid"}:
        raise ValueError("Distance metric must be centroid or boundary.")
    if metric == "boundary":
        first,second = projected_polygon(a),projected_polygon(b)
        if first is not None and second is not None:
            return polygon_gap(first,second),"boundary"
    return math.hypot(a["x"]-b["x"],a["y"]-b["y"]),"centroid_fallback" if metric=="boundary" else "centroid"


def boundary_candidates(cells,ids,radius):
    """Yield within-range pairs; bounding boxes include anchors for safe fallback."""
    polygons = {cid:projected_polygon(cells[cid]) for cid in ids}
    boxes = {cid:_bbox((polygons[cid] or [])+[[cells[cid]["x"],cells[cid]["y"]]]) for cid in ids}
    by_sample = {}
    for cid in ids:
        by_sample.setdefault(cells[cid]["sample_id"],[]).append(cid)
    for members in by_sample.values():
        ordered = sorted(members,key=lambda cid:(boxes[cid][0],cid))
        for i,cid in enumerate(ordered):
            a = boxes[cid]
            for j in range(i+1,len(ordered)):
                nid = ordered[j]
                b = boxes[nid]
                if b[0] > a[1]+radius:
                    break
                if bbox_gap(a,b) > radius:
                    continue
                if polygons[cid] is not None and polygons[nid] is not None:
                    distance,kind = polygon_gap(polygons[cid],polygons[nid]),"boundary"
                else:
                    distance,kind = distance_between(cells[cid],cells[nid])
                if distance <= radius:
                    yield distance,cid,nid,kind
