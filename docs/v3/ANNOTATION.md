# Existing labels and optional annotation utilities

The current workflow uses the labels already present in the supplied object.
It does not require new automatic annotation, H&E or independent manual review.
Record the originating pipeline/column, any explicit cell-type-to-lineage map
and its provenance. A map translates existing categories into supported model
permissions; it does not certify malignancy or create evidence for an identity.
Blank/missing labels cannot be invented. `unknown` remains an explicit category.

Use the [ROI catalogue](ROI_WORKFLOW.md) to inspect supplied polygons or centroid
fallback, composition, QC and expression support, then choose a region yourself.
The catalogue has no selected candidate initially and also accepts unknown-only
objects. Unknown-only historical smoke inputs remain useful engineering inputs,
but cannot demonstrate lineage-specific biological interactions.

The [confirmed scientific scope](SCIENTIFIC_SCOPE.md) supersedes the mandatory
annotation/review and automatic front/core selection in the historical plan.
The utilities below remain available **only when explicitly chosen** for a study.
Their requirements are not prerequisites for the catalogue or simulator.

## Optional historical melanoma ROI rule

Python API: `virtualtissue_xenium.v3.annotation.select_rois(cells, ...)`.
Each input cell needs a unique `cell_id`, `sample_id`, `cell_type`, supported
`lineage`, and measured `x`/`y` coordinates in micrometres. Selection is performed
separately per sample. Every sample must contain supplied `T_CD8` and `malignant`
labels. Preserve their declared origin. These labels need not be newly reviewed,
but the software does not establish their biological validity. H&E or reviewer
records may be included when available; neither is required by this utility.

When this optional function is called, its default rule is:

1. Map measured centroids to occupied nodes of a 500 µm grid anchored at `(0,0)`.
   Rounding is explicitly half upwards. Sort nodes by `(x,y)`; if there are more
   than 50, retain 50 evenly spaced indices of this ordering. This declared rule
   is deterministic and avoids searching unoccupied parts of a bounding box.
2. At each retained centre, select the nearest 1,000 cells in that sample.
   Distance ties use `cell_id`. A sample with fewer cells is rejected unless the
   investigator explicitly requests a smaller ROI.
3. Count distinct CD8–malignant centroid pairs within 30 µm. Select the centre
   with the highest count as the proposed invasive front; ties use `(x,y)`.
   If there are no contacts, return no front and a reason.
4. Select a proposed core with immune fraction **strictly below 5%** and
   malignant fraction **at least 50%**. By default its cells must not overlap the
   selected front. Rank eligible cores by descending malignant fraction, then
   ascending immune fraction and centre. If none qualifies, return no core and
   a reason; thresholds are never silently relaxed.

Parameters, all evaluated candidates, selected IDs, radius to the farthest
selected cell, unknown fractions and input hash are returned. A non-overlapping
core is an additional conservative default relative to the initial plan and may
be disabled explicitly. Centroid proximity is a geometric contact proxy, not
evidence of antigen recognition or a histological invasion boundary. The rule
must be declared before interpreting its outcomes. This automatic geometric
proposal is not the current catalogue's user-selection workflow, and the word
“front” does not establish a histological invasion boundary.

```python
from virtualtissue_xenium.v3.annotation import select_rois
proposals = select_rois(annotated_cells, cells_per_roi=1000)
assert proposals["complete"]  # Otherwise inspect samples[*].missing_reasons.
```

This tool selects ROIs for the melanoma design. It does not invent matched brain
cortical regions or substitute lung data for melanoma.

## Optional independent-review agreement

Use this tool only if two annotation files already exist for a comparison you
want to perform. The project does not require recruiting reviewers or generating
these files to use existing object labels.

`compare_reviewers(first_csv, second_csv, label_column="lineage", id_scope="exact")`
requires `cell_id` and the chosen label column. Blank labels and duplicate IDs
are rejected; unresolved cells must be labelled `unknown` explicitly. By default
the two CSVs must contain exactly the same IDs. An explicit
`id_scope="intersection"` restricts analysis to common IDs and lists excluded
IDs on both sides.

The result includes the confusion table, observed/expected agreement, Cohen's
κ, one-versus-rest κ per lineage and all disagreements. A zero expected-
disagreement denominator produces `null`, including two identical constant
label vectors; no common cells also produces `null`. `unknown` is retained as
an explicit category rather than silently discarded. There are no cell-level
significance tests or claims of biological truth.

`write_annotation_comparison(first_csv, second_csv, output_directory, ...)`
writes `agreement.json` and `disagreements.csv`, refusing to overwrite an occupied
directory. Input file hashes are recorded. If labels are subsequently adjudicated,
retain both original files and the adjudication record. Agreement measures
consistency between files, not biological truth or proof of independent review.
An unchanged supplied label can remain in the exploratory model with its origin
declared; sensitivity to uncertainty must be assessed when a paper claim depends
on that identity.
