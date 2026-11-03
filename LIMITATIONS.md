# Interpretation and remaining evidence

The current scope accepts the **existing labels and segmentation supplied in the
user's object**. H&E and independent manual annotation review are not requirements
for using the software or its exploratory paper workflow. Label origin and any
available verification are recorded as provenance. The user chooses an ROI from
its composition and geometry. The [scientific scope](docs/v3/SCIENTIFIC_SCOPE.md)
supersedes the mandatory annotation and automatic front/core requirements in the
historical plan.

Three claims remain hypotheses: expression changes simulated dynamics beyond
labels; Jev preferences respond to context beyond prior anchoring; and selected
simulated structural statistics are compatible with reserved real regions.
Software tests alone establish none of these empirical findings. A seven-day
simulated trajectory is not a demonstrated prediction of the tissue seven days
later, and a declared treatment modifier is not evidence of treatment efficacy.

## Inputs and geometry

- Labels, including `malignant`, determine model permissions but are not certified
  identities. An explicit lineage map translates supplied categories; it does not
  establish their biological correctness. Unknown labels remain conservative.
  Label perturbations can measure sensitivity without estimating the true error
  rate of the original annotation.
- Supplied polygons represent segmented units, not verified cell membranes.
  Geometry QC checks computational properties and records invalid/missing
  contours and centroid fallback. It cannot establish correct transcript
  assignment or detect every segmentation error. Runtime contours are rigid
  translations; daughter contours are explicitly synthetic inherited copies.
- User-selected candidate ROIs may overlap and are not independent replicates.
  Their displayed bounding rectangles are not statistical observation windows.
  The current cross-L estimator requires a declared rectangular window; matching
  a candidate selection to a valid observation design remains study work.

## Expression, rates and decisions

- Expression remains fixed at t0. Empirical sample/type/depth references describe
  relative evidence. Their mean-normalised multipliers preserve the baseline
  reference mean, not future events or absolute expression shifts between samples.
  The strict scoring method leaves unsupported references null. The optional
  `tiered_panel_auc` policy uses the whole type within a sample when necessary,
  then a bounded panel-based heuristic for rare types. That last tier is weak
  calibration: it neither guarantees depth-confound removal nor preserves the
  type's mean rate. Report empirical, weak and unavailable evidence separately.
- Available evidence can affect only permitted, currently relevant actions.
  A changed probability need not change a particular sampled action. A
  probability audit demonstrates a software mechanism; outcome contrasts across
  seeds and permutations are still needed to evaluate an effect on dynamics.
  Missing panel genes or absent detected signature expression are not fabricated
  into evidence. In Jev mode, beta-zero removes the prior's expression multiplier
  but does not alone remove expression from the provider payload.
- Default hazards, speeds, state durations, crowding and occupancy rules remain
  assumptions. Verified reference identities do not calibrate numerical rates to
  a supplied tissue. Sensitivity analysis can identify dependence on assumptions;
  it does not turn them into measured parameters. Population doubling, cell-cycle
  duration and turnover are not interchangeable rate measurements.
- The local model allows one competing-risk event per unit per step. The clock
  test covers restricted homogeneous growth; integrated interactions/states need
  their own step-size checks. Jev-hidden forces alpha=1 and its returned event
  mass need not preserve the local hazard clock. Provider replicate variability
  and dependence on the simulation step must be measured separately.
- Movement, signals, Treg suppression, crowding and repulsion are simplified
  mechanisms. Repulsion does not guarantee nonoverlap. Molecular diffusion,
  pharmacokinetics, inferred antigen recognition and evolving expression are not
  implemented.

## Study and software boundaries

- Seed/bootstrap uncertainty concerns the simulation. Permutation variability,
  provider replicates and biological sampling are different sources of variation.
  Multiple regions from one section do not establish patient-level generalisation.
- Structural compatibility needs reference regions reserved before parameter
  adjustment, declared statistics/windows, static-input comparisons and controls.
  Without those data the structural claim is unassessed. Agreement with a static
  spatial reference does not validate a future trajectory.
- Lattice outputs are occupied bins and abstract occupancy, not cell counts.
  New native-platform adapters have format tests, but public-data and comparable-
  unit cross-platform validation remain incomplete. Direct GEF input and mask
  aggregation are unimplemented. Observed-only GEM gene unions do not establish
  complete-transcriptome score comparability.
- Raw-JSON ROI browsing avoids whole-section scoring; the catalogue still loads
  the supplied object/geometry in memory. It is not a metadata-only H5AD browser.
  Replay loads at most two 50-frame chunks; metadata, identifiers and complete
  audit files still grow. Copy its sidecar folder with the HTML.
- Synthetic demos and historical unknown-only smoke tests do not constitute
  expression, Jev-context or structural findings on a real supplied object. No
  paid Jev study results were generated by this implementation.

## Consequences for claims

Listing a limitation does not resolve it. Test dependence on labels, QC,
segmentation geometry, scoring tier, rates and step size where a claim depends
on them. Predeclare alternatives and report all tested outcomes; do not tune
beta, signatures or ROI selection merely to obtain a positive contrast.

If evidence is unavailable, omit the corresponding empirical claim or mark it
unassessed. If a result is sensitive to an assumption, restrict the claim to the
tested setting. If controls show little or no effect, report that finding with
its uncertainty; overlapping seed bands alone do not prove equivalence or absence
of information. A reproducible exploratory methods/resource article can retain
negative findings without claiming biological prediction or validated treatment
effects. See [paper/README.md](paper/README.md) for the claim-to-evidence contract.
