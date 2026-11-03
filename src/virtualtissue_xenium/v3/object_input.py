"""Explicit polygon joins for user objects; no segmentation or identity inference."""
from copy import deepcopy
import json
import math


def _plain(value):
    if hasattr(value, "__geo_interface__"):
        return _plain(value.__geo_interface__)
    if hasattr(value, "tolist"):
        return _plain(value.tolist())
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def polygon_map(raw):
    """Read an ID-to-ring map or GeoJSON FeatureCollection with cell_id properties."""
    raw = _plain(raw)
    if not isinstance(raw, dict):
        raise ValueError("Polygons must map exact cell IDs to rings or be a GeoJSON FeatureCollection.")
    if raw.get("type") == "FeatureCollection":
        result = {}
        features = raw.get("features")
        if not isinstance(features, list):
            raise ValueError("GeoJSON features must be a list.")
        for feature in features:
            if not isinstance(feature, dict) or feature.get("type") != "Feature":
                raise ValueError("GeoJSON must contain Feature records.")
            cid = feature.get("properties", {}).get("cell_id", feature.get("id"))
            if not isinstance(cid, str) or not cid or cid in result:
                raise ValueError("Polygon features require unique exact string cell_id values.")
            result[cid] = feature.get("geometry")
        return result
    return raw


def _scale_polygon(value, scale):
    value = _plain(value)
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except ValueError:
            return value  # Geometry QC records the invalid source without inventing vertices.
    if scale == 1:
        return value
    def coordinates(items):
        if not isinstance(items, list):
            return items
        if len(items) >= 2 and all(type(v) in (int, float) for v in items[:2]):
            return [items[0] * scale, items[1] * scale, *items[2:]]
        return [coordinates(item) for item in items]
    if isinstance(value, dict) and "coordinates" in value:
        return {**value, "coordinates": coordinates(value["coordinates"])}
    if isinstance(value, dict) and value.get("type") == "Feature":
        return {**value, "geometry": _scale_polygon(value.get("geometry"), scale)}
    return coordinates(value)


def attach_polygons(cells, raw, *, coordinate_scale=1.0):
    """Join only matching IDs and scale vertices to µm; unmatched IDs are reported."""
    if type(coordinate_scale) not in (int, float) or not math.isfinite(coordinate_scale) or coordinate_scale <= 0:
        raise ValueError("Polygon coordinate_scale must be finite and positive.")
    mapping = polygon_map(raw)
    result = deepcopy(cells)
    ids = {c["cell_id"] for c in result}
    matched = ids & set(mapping)
    if not matched:
        raise ValueError("No polygon IDs match the selected input cells.")
    for cell in result:
        if cell["cell_id"] in mapping:
            if any(key in cell for key in ("polygon", "polygon_um", "segmentation_polygon", "geometry", "t0_geometry")):
                raise ValueError("A cell already has geometry; choose one explicit polygon source.")
            cell["polygon"] = _scale_polygon(mapping[cell["cell_id"]], coordinate_scale)
    return result, dict(join="exact_cell_id", coordinate_scale_to_um=coordinate_scale,
                        matched_cells=len(matched), cells_without_polygon=len(ids-matched),
                        unused_source_polygon_ids=sorted(set(mapping)-ids))


def attach_anndata_polygons(cells, adata, key, *, coordinate_scale=1.0):
    """Explicit key: uns[key] ID map/GeoJSON, or obs[key] per-cell rings/JSON."""
    if key in adata.uns:
        raw = adata.uns[key]
    elif key in adata.obs:
        raw = {str(cid): _plain(value) for cid, value in zip(adata.obs_names, adata.obs[key])}
    else:
        raise ValueError(f"Polygon key {key!r} is absent from AnnData.uns and AnnData.obs.")
    return attach_polygons(cells, raw, coordinate_scale=coordinate_scale)
