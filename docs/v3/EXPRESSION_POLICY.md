# Expression evidence for common and rare types

The optional schema-3.0 method `tiered_panel_auc` lets supported expression influence compatible local actions even when a declared cell type has too few reference cells. It uses the existing supplied labels, counts and segmentation. Those inputs are declarations of uncertain quality; this method does not make their biological identities or boundaries independently verified.

The existing `depth_conditioned_auc` method and schema-2.0 behavior are preserved. Select the new method explicitly in a new context and rebuild preparation; do not edit a prepared fingerprint in place.

```json
{
  "scoring": {
    "method": "tiered_panel_auc",
    "reference_scope": "sample",
    "reference_fallback": "cell_type",
    "depth_bins": 10,
    "min_reference_cells": 30,
    "min_overlap": 5,
    "panel_max_multiplier": 2
  }
}
```

Unspecified values retain the schema-validated platform/kinetics defaults; these defaults are not biologically calibrated. This example does not change signatures, rate values or beta to force a positive result. The optional `panel_max_multiplier` is a modelling assumption constrained to `[1, 3]`; the default is 2. A value of 1 disables the panel fallback's local multiplier while retaining its evidence for inspection. Beta zero disables local expression multipliers in all tiers.

For this method, `reference_fallback` defaults to `cell_type`. An explicit `null` fallback is rejected because it contradicts the tier hierarchy; use `depth_conditioned_auc` with `reference_fallback: "null"` for the strict policy. When adapting an already normalized strict context, change both declarations and rebuild preparation.

## Reference hierarchy

Every signature first needs the declared minimum overlap with the assayed panel; every cell must pass the transcript-count QC gate. Missing genes are not replaced with invented measurements.

1. **Sample × declared type × depth:** use empirical AUC percentiles when every depth stratum for that sample/type reaches the reference minimum.
2. **Sample × declared type:** if any stratum is too small but the type has enough eligible cells overall, that entire type uses its sample-specific reference without depth conditioning. The whole-type fallback creates one consistent, non-overlapping reference partition; cells are not partly normalized against overlapping groups.
3. **Cell's measured panel:** when the type still has too few cells, compare its signature's rank-recovery AUC to an explicit random-gene-set model of that cell's panel. This is weak panel evidence, not a percentile of a reference-cell population.

The new method always keeps cell references within a sample and never pools distinct types or parent lineages. Tier 2 loses depth matching and is identified as `cell_type_no_depth`. Tier 3 has `reference_size: 0`, `candidate_cell_reference_size` identifying the undersized cell group, `reference: panel_gene_set_null`, and `calibration: weak_panel_calibration`. Its genes are a mathematical reference, not extra cells.

## Panel calculation and bounds

For each measured gene, let `w` be its contribution to the top-rank recovery curve. Positive counts have deterministic average ranks for ties; undetected genes contribute zero. Let `N` be the number of assayed genes, `k` the signature overlap and `B` the largest possible recovery sum for `k` genes in the declared rank window. Raw signature AUC is `sum(signature w) / B`.

Under uniformly selected size-`k` gene subsets **without replacement** from this panel, the exact moments are:

```text
E[AUC]   = k × mean(w) / B
Var[AUC] = k × (N − k) / (N − 1) × population_variance(w) / B²
z_raw    = (observed_AUC − E[AUC]) / sqrt(Var[AUC])
z        = clip(z_raw, −3, 3)
```

The implementation verifies these moments against exhaustive gene-subset enumeration in a small tied/zero-count fixture. It does not infer a normal distribution, p-value, validated pathway activity or protein abundance. Uniform random gene sets are an explicit reference assumption; targeted-panel composition and systematic differences between genes can violate its biological relevance. This is not external calibration or complete correction of depth/batch effects.

When no signature gene is detected, panel variance is zero, QC fails or overlap is insufficient, the panel tier stays unscored and its local multiplier is exactly 1. Non-detection is not treated as evidence of pathway inactivity. With detected genes and nonzero reference variance, positive or negative `z` describes relative rank enrichment or depletion under the stated assumption.

The panel tier applies:

```text
m = exp(clip(beta × z, −log(panel_max_multiplier), +log(panel_max_multiplier)))
```

Thus its default multiplier is within `[0.5, 2]`. It is not divided by the one rare cell's own value, which would cancel every effect. **Its group mean is not preserved**, so the assumed baseline rate becomes a reference level rather than the exact rare-type population average. This changed interpretation is recorded explicitly.

Empirical tiers continue using `exp(beta × z) / mean_reference(exp(beta × z))`, preserving their reference-group mean. The panel bound does not constrain empirical-tier multipliers. Gates, permitted actions, crowding and competing risks still apply after expression scoring. A score cannot license a forbidden action or override an `unknown` lineage's REST-only menu.

## Fixed-state influence audit

`virtualtissue_xenium.v3.expression_audit.audit_expression(bundle, seed=0)` returns a JSON-compatible report without sampling transitions or calling Jev. It rebuilds beta-zero, within-type expression shuffle and global expression shuffle arms, then evaluates their local action probabilities on the **same initialized snapshot**, graph, target ordering, labels, positions and clocks.

The report includes coverage/status by type and signature; changed probabilities and hazards per cell/action; Jensen–Shannon divergence and total variation; changes in QC or offered actions; and hashes of the frozen state and graph. It distinguishes a directly changed action hazard from redistribution of probabilities under competing risks, including REST. A singleton's within-type shuffle is necessarily unchanged. Global shuffling can move QC/depth characteristics across types and is not a biologically realistic alternate tissue.

This audit answers “does the measured expression affect this model's decisions in this state?” It does not demonstrate that those decisions are correct biology, guarantee a sampled action changes, or establish effects over a trajectory. A legitimate result is zero influence because beta is zero, the signature is unsupported, ranks lack information, an action is gated or every permitted action is REST. Rates, signatures and beta are never optimized by this audit.

To support a claim about **simulated dynamics**, follow the fixed-state audit with matched multi-seed trajectories and predeclared endpoint/curve readouts, as specified in [the scientific protocol](SCIENTIFIC_SCOPE.md). Keep permutation realizations separate from simulation seeds. A difference in one initial probability is functional evidence about the implementation, not sufficient evidence for a trajectory-level claim. If the trajectory effect is uncertain, report the narrower probability-level result.

The frozen snapshot includes the supplied geometry. When boundary distance is enabled, polygons can change neighbours, target availability and local evidence, but these remain fixed within each expression contrast. Compare boundary-versus-centroid geometry separately to avoid attributing a changed action opportunity to expression. Current contact/distance readouts in `readouts.yaml` remain **centroid-based and independent of the capped action graph**, even when that graph uses polygon boundaries. They do not measure membrane contact or directly validate the boundary graph.

## Local and Jev attribution are distinct

Beta zero removes local expression modulation but keeps t0 evidence in the Jev observation. Shown-prior Jev can receive expression through both its numeric prior and direct evidence; hidden-prior Jev receives direct evidence and a constrained action menu. The audit reports which input routes change, matching the compact payload's score, status, reference basis and signature definitions. It does **not** claim that Jev's response changes.

Changing beta also changes the full-context hash included in a provider payload. A matched provider experiment must control that opaque metadata difference as well as the prior, direct evidence, ordering and state. Evidence-only, prior-only and combined interventions need independent recorded responses before attributing any Jev effect. No new provider calls occur as part of the local expression audit.

The implementation version recorded in evidence is `tiered-panel-auc-1`. Report empirical-tier and weak-panel-tier coverage separately in any figures; an available panel score must not be counted as 30 biological reference cells or as independently validated activity. The new method and audit have synthetic implementation checks. No real-data effect study or live Jev response experiment for these new features has been performed in this implementation phase.
