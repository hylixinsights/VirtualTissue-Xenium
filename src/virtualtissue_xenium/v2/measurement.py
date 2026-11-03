"""Versioned measured-input contract; declarations never infer cell identity."""
from copy import deepcopy
import math


def organism_key(value):
    aliases = {"human": "human", "homo sapiens": "human", "mouse": "mouse", "mus musculus": "mouse"}
    return aliases.get(str(value).strip().casefold(), str(value).strip().casefold())


def comparison_signature(measurement):
    """Measurement semantics that cannot be changed in a matched-arm comparison."""
    if measurement is None:
        return None
    keys = ("platform", "organism", "unit", "unit_size_um", "panel_type", "n_features",
            "is_single_cell_resolution", "segmentation_method", "deconvolution_method", "purity_threshold")
    return {k: measurement.get(k) for k in keys}


def validate_measurement(raw, cells, genes, context):
    m = deepcopy({} if raw is None else raw)
    if not isinstance(m, dict):
        raise ValueError("measurement must be an object.")
    platform = context["platform"]
    m.setdefault("platform", platform)
    m.setdefault("organism", context["organism"])
    if m["platform"] != platform or organism_key(m["organism"]) != organism_key(context["organism"]):
        raise ValueError("Measurement platform/organism contradicts the context.")
    if platform in ("visium_hd", "stereo_seq") and "unit" not in m:
        raise ValueError("Declare measurement.unit for Visium HD/Stereo-seq; bins must not be assumed to be cells.")
    m.setdefault("unit", "nucleus" if platform == "snrna" else "segmented_cell")
    if m["unit"] not in ("segmented_cell", "binned_spot", "nucleus"):
        raise ValueError("measurement.unit must be segmented_cell, nucleus, or binned_spot.")
    m.setdefault("unit_size_um", None)
    if m["unit_size_um"] is not None and (type(m["unit_size_um"]) not in (int, float) or
            not math.isfinite(m["unit_size_um"]) or m["unit_size_um"] <= 0):
        raise ValueError("measurement.unit_size_um must be null or positive micrometres.")
    m.setdefault("panel_type", "targeted" if platform == "xenium" else "whole_transcriptome")
    if m["panel_type"] not in ("targeted", "whole_transcriptome"):
        raise ValueError("measurement.panel_type must be targeted or whole_transcriptome.")
    m.setdefault("segmentation_method", "not_declared")
    if not isinstance(m["segmentation_method"], str) or not m["segmentation_method"].strip():
        raise ValueError("Declare segmentation_method, or use not_declared.")
    single = m["unit"] in ("segmented_cell", "nucleus")
    m.setdefault("is_single_cell_resolution", single)
    if type(m["is_single_cell_resolution"]) is not bool or m["is_single_cell_resolution"] != single:
        raise ValueError("measurement resolution must agree with the unit; binned spots are not cells.")
    if not single and (context["kernel"] != "lattice" or m["unit_size_um"] is None):
        raise ValueError("Binned spots require the lattice kernel and an explicit unit_size_um.")
    if not single:
        m.setdefault("purity_threshold", .7)
        if type(m["purity_threshold"]) not in (int, float) or not .7 <= m["purity_threshold"] <= 1:
            raise ValueError("Binned-unit purity_threshold must be from 0.7 to 1.")
        m.setdefault("deconvolution_method", "not_declared")
        if not isinstance(m["deconvolution_method"], str) or not m["deconvolution_method"].strip():
            raise ValueError("Declare deconvolution_method or use not_declared.")
    if single and context["kernel"] != "cell_agent":
        raise ValueError("Cell/nucleus units require the cell_agent kernel.")
    count = len(genes)
    mean = sum(sum(c["expression"].values()) for c in cells) / len(cells)
    for key, measured in (("n_features", count), ("mean_transcripts_per_unit", mean)):
        if key in m and (type(m[key]) not in (int, float) or not math.isfinite(m[key]) or
                         not math.isclose(m[key], measured, rel_tol=1e-10, abs_tol=1e-10)):
            raise ValueError(f"measurement.{key} disagrees with the selected input.")
        m[key] = measured
    m["coordinate_units"] = "micrometres" if context["interaction_mode"] == "spatial" else "not_spatial"
    if context["interaction_mode"] == "spatial" and any("x" not in c or "y" not in c for c in cells):
        raise ValueError("Spatial input requires measured x/y in micrometres for every unit.")
    return m


def read_lineage_map(path):
    import csv
    from pathlib import Path
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not {"cell_type", "lineage"} <= set(reader.fieldnames or []):
            raise ValueError("lineage_map.csv requires cell_type and lineage columns.")
        result = {}
        for row in reader:
            cell_type, lineage = row["cell_type"].strip(), row["lineage"].strip()
            if not cell_type or not lineage or cell_type in result:
                raise ValueError("Lineage mapping requires unique, nonempty cell types and lineages.")
            result[cell_type] = lineage
    if not result:
        raise ValueError("The lineage mapping is empty.")
    return result


def apply_lineage_map(cells, mapping):
    result = deepcopy(cells)
    for c in result:
        if c.get("cell_type") not in mapping:
            raise ValueError(f"Missing lineage mapping for {c.get('cell_type')!r}.")
        c["lineage"] = mapping[c["cell_type"]]
    return result
