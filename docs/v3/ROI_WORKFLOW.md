# Choose a region using supplied annotations and geometry

The ROI catalogue is an offline inspection and selection tool. It accepts a raw
JSON object with `cells` and optional `assayed_genes`/`annotation_provenance`, or
a spatial schema 3.0 prepared bundle, including `unknown` labels.
It does not require H&E, independent reviewers, or a malignant/CD8 population.
Record where labels came from in `annotation_provenance`; provenance is a
declaration, not a claim that the labels were independently verified. The tool
does not infer malignancy or change labels to create a desired composition.

```bash
vtx2 roi-catalog full_sample.json regions.html --cells-per-roi 1000 --grid-um 500
```

Open `regions.html` directly in a browser. Choose a sample, inspect the map and
candidate compositions, then select a candidate by its list entry or map centre.
No candidate is selected initially. The download button becomes available only
after a selection. The downloaded `*.selection.json` contains exact IDs, sample,
geometry bounds, generator parameters and source hashes. The catalogue also
writes `regions.catalog.json` for reproducible inspection. No network service or
Jev request is involved, and the source bundle is never overwritten.

```bash
vtx2 select-roi full_sample.json roi_example.selection.json selected_region.json --context context.json
```

Selection is validated against the original source fingerprint and data hash,
then exact candidate membership, centre and bounds. Using a selection with a
different source or editing its membership is rejected. Preparing selected cells
again may change scoring references: inspect the resulting prepared bundle
before simulating. The original data remain intact. `--context` is required when
selecting from raw input; a prepared source already has a context. Raw inspection
does not build a guidebook or calculate scores for the entire section: only the
selected cells are prepared for simulation. Cells need IDs, supplied type/lineage,
sample IDs, expression count dictionaries and measured coordinates (or valid
polygons from which the geometry helper derives an explicitly recorded centroid).
This route does not provide a metadata-only H5AD browser; export the supported
JSON object or use a bounded prepared input.

Existing `cell_type` labels are retained. If the object does not already contain
the model's `lineage` permission categories, supply an explicit
`cell_type,lineage` CSV using `--lineage-map` on **both** `roi-catalog` and
`select-roi`. This maps existing labels to allowed model actions; it does not
perform a new annotation or require manual review of individual cells. The map
is recorded and hashed with the source. Unsupported or uncertain categories can
be explicitly mapped to `unknown`. A prepared source already freezes this map.

Selection preserves declared measurement units and annotation provenance. It
also records the original sampling description under `sampling.source_sampling`
and retains its declared missing populations alongside the exclusion of units
outside the chosen ROI. Scores and transcript summaries are recomputed for the
selected cells. JSON `context` metadata is retained in the source fingerprint;
the explicit `--context` file determines the simulation settings for raw input.

## What candidates and metrics mean

Candidates use occupied nodes of a grid anchored at `(0,0)` in micrometres,
rounded half upwards. Per sample, nodes are sorted spatially and at most 50
evenly spaced indices are retained by default. Each centre selects its nearest
1,000 units, breaking distance ties by exact cell ID. The parameters are
configurable. Samples with fewer units are retained and visibly marked as
undersized. Candidates can overlap, including identical memberships; they are
alternatives for the user, not independent replicates or a biological ranking.

The map displays supplied cell polygons when available. Missing polygons use
measured centroids, identified explicitly. Polygon coordinates must use the same
micrometre coordinate system as centroids. Candidate membership uses centroid
distance; polygons do not implicitly include or exclude additional cells.
Normalized `t0_geometry` records are supported, including the recorded invalid-
polygon centroid fallback. Geometry QC and invalid-contour counts remain visible;
the catalogue does not repair contours or replace their recorded provenance.
The selection's bounding rectangle is a display bound, explicitly marked
`is_statistical_observation_window: false`; it must not be treated automatically
as the observation window for spatial-statistical inference.

The catalogue reports lineage counts/fractions, unknown fraction, QC pass/fail,
total and median transcript counts, detected genes, and per-signature source
availability and informative coverage. An available score is not necessarily
informative: a constant reference distribution has no relative variation.
Source score summaries describe the supplied prepared bundle. Empirical
availability and `available_panel_heuristic` are counted separately. The latter
is explicitly weak panel calibration, not empirical support from additional
cells or a validated biological effect. These diagnostics and the supplied
composition are also retained in the downloaded selection descriptor. Counts of
QC-passing units by type and within-ROI depth stratum indicate potential local
reference support; they do not promise a score for every reference strategy.
The catalogue does not recompute or silently change the scoring method.

For raw input, signature coverage is explicitly `not_computed`. Transcript and
reference counts use visible descriptive defaults (20 transcripts, 10 depth bins,
30 reference units), not a fitted score or a certified QC decision. They are not
copied into the simulation context. The actual context and platform determine
QC and scoring when the selected ROI is prepared. The raw source fingerprint
hashes its metadata together with normalized cells/geometry; invalid polygon
numbers remain preserved inside inert raw-input text rather than nonfinite JSON
numbers. A second hash covers only the normalized cells used in selection. Both
are checked on extraction. Nonfinite metadata outside geometry is rejected.

Unknown-only catalogues are supported and exportable, but `unknown` remains
subject to the model's explicit permissions. Inspecting or selecting a region
does not convert uncertain labels into biological evidence. The earlier
front/core `roi-candidates` workflow remains a separate, opt-in geometric rule;
this catalogue does not call it or auto-select its highest-contact region.

## Python interface

```python
from virtualtissue_xenium.v3.roi_catalog import (
    build_roi_catalog, write_roi_catalog, validate_roi_selection,
    extract_roi_selection,
)
catalog = build_roi_catalog(bundle, cells_per_roi=1000, grid_um=500, max_candidates=50)
paths = write_roi_catalog(bundle, "regions.html", cells_per_roi=1000)
descriptor = validate_roi_selection(bundle, "chosen.selection.json")
selected_cells = extract_roi_selection(bundle, descriptor)
```

The extraction function returns independent copies in ID order. It does not
prepare a bundle or mutate the source. All source data and the HTML catalogue
remain local; the standalone HTML embeds geometry, supplied labels and summary
metrics, so share it only where those data are intended to be shared.
