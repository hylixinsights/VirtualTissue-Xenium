"""Generate a clearly synthetic, annotated polygon object; select no ROI.

Run from the repository with PYTHONPATH=src. The 96 artificial cells exercise
geometry QC and both common-type and rare-type scoring. They are not a fitted
tissue model, a measured Xenium sample, or biological validation data.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import random

from virtualtissue_xenium.v2.inputs import write_json
from virtualtissue_xenium.v3.ontology import lineage_rules, load_signatures


GENERATOR_VERSION = "synthetic-polygon-object-1"
DEFAULT_SEED = 20261101
PATCHES = (
    ((15.0, 15.0), dict(malignant=24, T_CD8=4, myeloid=2, fibroblast=2)),
    ((155.0, 30.0), dict(malignant=12, T_CD8=10, myeloid=6, T_reg=2, NK=2)),
    ((75.0, 125.0), dict(malignant=4, T_CD8=4, myeloid=4, fibroblast=6,
                         endothelial=6, T_reg=2, NK=2, B=4)),
)
PROGRAMS = {
    "malignant": ("cell_cycle", "hypoxia"),
    "T_CD8": ("cytotoxic", "activation", "exhaustion_associated"),
    "myeloid": ("inflammation", "antigen_presentation"),
    "fibroblast": ("TGFb_response",),
    "endothelial": ("hypoxia",),
    "T_reg": ("activation", "TGFb_response"),
    "NK": ("cytotoxic", "activation"),
    "B": ("MHC_II", "antigen_presentation"),
}


def create_object(seed=DEFAULT_SEED, *, include_invalid_polygon=False):
    """Return a raw object with deterministic synthetic provenance and context."""
    rng = random.Random(seed)
    signatures, signature_source = load_signatures("human")
    rules = lineage_rules()
    if set(PROGRAMS) - set(rules):
        raise ValueError("The demo requires supported v3 lineage labels.")
    background = [f"SYNTHETIC_BACKGROUND_{i:02d}" for i in range(48)]
    biological_panel = sorted({gene for genes in signatures.values() for gene in genes})
    panel = sorted(set(biological_panel) | set(background))
    cells = []
    for patch_index, ((origin_x, origin_y), composition) in enumerate(PATCHES):
        labels = [lineage for lineage, count in composition.items() for _ in range(count)]
        rng.shuffle(labels)
        for local_index, lineage in enumerate(labels):
            index = len(cells)
            x = round(origin_x + (local_index % 8) * 15 + rng.uniform(-1.5, 1.5), 4)
            y = round(origin_y + (local_index // 8) * 15 + rng.uniform(-1.5, 1.5), 4)
            depth = (1, 2, 5)[index % 3]
            expression = {gene: rng.randint(1, 4) * depth for gene in rng.sample(background, 24)}
            expression.update({gene: rng.randint(1, 3) * depth for gene in rng.sample(biological_panel, 12)})
            # These are assigned fixture patterns, not inferred biological states.
            # Alternate strong/weak patterns and total counts within each label.
            for program_index, program in enumerate(PROGRAMS[lineage]):
                amplitude = (2, 7, 18, 3)[(local_index + program_index) % 4]
                for gene in signatures[program][:8]:
                    expression[gene] = (amplitude + rng.randint(0, 3)) * depth
            if index % 7 == 0:
                for gene in signatures["apoptosis_associated"]:
                    expression[gene] = rng.randint(4, 10) * depth
            radius = rng.uniform(4.5, 6.3)
            phase = rng.uniform(0, math.tau)
            ring = []
            for vertex in range(8):
                angle = phase + math.tau * vertex / 8
                scale = rng.uniform(.85, 1.15)
                ring.append([round(x + radius * scale * math.cos(angle), 4),
                             round(y + radius * scale * math.sin(angle), 4)])
            ring.append(list(ring[0]))
            cells.append(dict(cell_id=f"SYNTHETIC_{index:03d}", sample_id="SYNTHETIC_ONLY",
                              cell_type=lineage, lineage=lineage, x=x, y=y,
                              expression=expression, polygon=ring,
                              synthetic_patch=f"artificial_patch_{patch_index + 1}"))

    # One declared low-count object exercises QC; its label and shape remain data.
    low_count = cells[-1]
    low_count["expression"] = {background[0]: 2, background[1]: 1}
    low_count["synthetic_qc_case"] = "deliberate_total_3_below_min_transcripts_20"
    invalid_ids = []
    if include_invalid_polygon:
        cell = cells[0]
        x, y = cell["x"], cell["y"]
        cell["polygon"] = [[x-5, y-5], [x+5, y+5], [x-5, y+5], [x+5, y-5], [x-5, y-5]]
        cell["synthetic_qc_case"] = "deliberate_self_intersection_expected_centroid_fallback"
        invalid_ids.append(cell["cell_id"])

    context = json.loads(Path(__file__).with_name("context_user_object.json").read_text())
    context.update(tissue="SYNTHETIC polygon layout; no measured tissue",
                   disease="SYNTHETIC engineering fixture; no disease claim")
    context["sampling"] = dict(design="subsampled",
        description="96 generated objects in three arbitrary patches; user must choose an ROI from its composition.",
        missing_populations=["No biological source population exists for this synthetic fixture"])
    return dict(
        description="SYNTHETIC ONLY: arbitrary counts, labels, coordinates and polygons for engineering checks. No ROI selected; no biological or longitudinal validation.",
        cells=cells, assayed_genes=panel, context=context,
        measurement=dict(platform="xenium", organism="human", unit="segmented_cell",
                         panel_type="targeted", coordinate_units="micrometres",
                         segmentation_method="synthetic_closed_eight_vertex_polygons_not_an_assay",
                         is_single_cell_resolution=True, synthetic=True,
                         acquisition="No instrument acquisition; Xenium contract used for engineering only"),
        annotation_provenance=dict(
            source="synthetic_generator", verification="not_independently_verified",
            synthetic=True, generator=GENERATOR_VERSION, seed=seed,
            label_source="Assigned existing v3 lineage names, not inferred from expression or shape",
            expression_source="Arbitrary integer counts and assigned signature patterns; not fitted to biological data",
            signature_file_sha256=signature_source["file_sha256"],
            background_symbols="SYNTHETIC_BACKGROUND_* are artificial panel features",
            expected_invalid_polygon_cells=invalid_ids,
            expected_low_count_cells=[low_count["cell_id"]],
            expected_scoring="Whole-type empirical references for malignant; weak panel heuristic for rare types; one QC failure",
            roi_selection="none; requires explicit user selection",
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New raw JSON object path.")
    parser.add_argument("--context-output", type=Path, help="Optional new context JSON for explicit CLI --context.")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--include-invalid-polygon", action="store_true",
                        help="Add one declared self-intersecting ring to test centroid fallback.")
    args = parser.parse_args()
    destinations = [args.output] + ([args.context_output] if args.context_output else [])
    if len({path.resolve() for path in destinations}) != len(destinations):
        parser.error("Object and context output must be different paths.")
    if any(path.exists() for path in destinations):
        parser.error("Choose new output paths; this generator does not overwrite existing files.")
    obj = create_object(args.seed, include_invalid_polygon=args.include_invalid_polygon)
    write_json(args.output, obj)
    if args.context_output:
        write_json(args.context_output, obj["context"])
    print(json.dumps(dict(output=str(args.output), synthetic=True, cells=len(obj["cells"]),
                          composition=dict(Counter(c["lineage"] for c in obj["cells"])),
                          invalid_polygons=len(obj["annotation_provenance"]["expected_invalid_polygon_cells"]),
                          selected_roi=None), sort_keys=True))


if __name__ == "__main__":
    main()
