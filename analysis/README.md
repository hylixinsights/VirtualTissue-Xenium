# Analysis contract for v3

`readouts.yaml` is JSON-compatible YAML. It is byte-identical to the installed
`virtualtissue_xenium/v2/readouts_v3.json` asset. The ensemble records the SHA-256
and rejects definitions that differ from its implemented readout version. This
is a versioned analysis definition, not an external preregistration or evidence
that biological validation has been completed.

## Seed × arm ensembles

`vtx2 compare RUN_A RUN_B` compares recorded runs at the latest common completed
step and rejects different recorded ontology hashes or lineage declarations.
`--allow-ontology-diff` permits a forensic comparison of existing artifacts and
lists changed asset hashes and the differing lineage rows recorded in each
guidebook. Historical kinetics/other asset rows are unavailable when they were
not saved; the command does not reconstruct them from today's installation.
Such a contrast is explicitly marked as confounded by ontology changes. The
override does not relax preparation, execution or resume validation, nor other
measurement/execution compatibility checks. Lattice comparisons report abstract
occupancy and occupied bins instead of cell population counts.

Compare's `cost_ledger` sums actual reported provider usage over the **entire**
runs, including any work after the common comparison endpoint. It reports known
token subtotals, missing usage fields and attempts with unknown usage separately;
there is no assumed token price. Repeated paths are rejected to avoid counting
the same run twice. Provider divergence from each report also describes the
entire run and is labelled with that scope.

The Python entry point is `virtualtissue_xenium.v2.ensemble.run_ensemble`:

```python
run_ensemble(
    {"control": "control.json", "shuffled": "shuffled.json"},
    "outputs/ensemble", seeds=list(range(10)), horizon_hours=24,
    bootstrap_replicates=2000, bootstrap_seed=0,
)
```

The output directory must be new. Prepared bundles are copied into it, and the
configuration records their fingerprints, contexts, controls, run settings and
readout hash. Each seed/arm gets its own run directory and journal. Existing
ensembles are never overwritten or automatically retried. Different ROIs belong
in separate ensembles; controls must identify the same original measured input.

Use either `steps` or `horizon_hours`, not both. Different time steps require a
shared horizon divisible by each time step. Paired results use exact common
committed hours; there is no interpolation, horizon extrapolation, or implicit
conversion of an early stopped run into a completed one.

`report.json` identifies incomplete runs and excluded arms. A single incomplete
seed excludes that entire arm from aggregation; its available readouts remain
in `readouts.csv` with `completed=false`. Runs that stop when no living cells
remain are currently recorded as incomplete if the kernel did not produce the
requested horizon; the analysis does not manufacture additional frames.

## What the intervals mean

- `summary.csv`: median per arm, time, readout and lineage, using one value per
  simulation seed. `median_ci95_low/high` are percentile bootstrap confidence
  limits for the median, resampling seeds with replacement.
- `seed_q025/q975`: empirical seed quantiles with linear interpolation. These
  are seed envelopes, distinct from confidence limits for the median.
- `paired_differences.csv`: arm minus the first arm, calculated within the same
  seed and common time point, followed by the same summary operations. Additional
  `paired_dz*` columns report the standardized paired effect described below;
  raw medians and their intervals are unchanged.
- `readouts.csv`: all individual seed-level values, including explicit nulls,
  denominator counts and reasons for unavailability.

With fewer than two finite seeds, median confidence limits are undefined. With
fewer than five requested seeds, the report warns that the run is an engineering
check. Use at least ten seeds for the headline figures proposed in the plan.
Seeds quantify simulation variability; they do not substitute for patients,
sections, independently annotated ROIs or experimental replicates. There are no
tests treating simulated cells as independent biological replicates.

The standardized contrast is Cohen's `dz = mean(differences) / SD(differences)`,
using sample SD with `ddof=1`; `paired_mean_difference` and
`paired_difference_sd` expose both quantities. This is the paired-difference
standardizer described in [Lakens (2013), equation 6](https://www.frontiersin.org/journals/psychology/articles/10.3389/fpsyg.2013.00863/full).
Here each observation is a matched **simulation seed**, not a patient. It is not
a standardized median, independent-group effect, bias-corrected effect or
clinical effect size. Pairing is preserved by resampling whole differences,
using only seeds with both finite values. Compatible lattice arms may use this
calculation, but mixed cell/lattice units and different occupancy capacities
are rejected before execution.

`paired_dz_ci95_low/high` are percentile bootstrap limits, distinct from both
the raw median confidence interval and seed envelope. Fewer than two finite
pairs or zero original SD yields null `dz` and interval. A zero-SD bootstrap
sample makes that bootstrap ratio undefined, even if the original `dz` exists.
The implemented conservative policy then leaves the interval null and reports
`paired_dz_ci_status=undefined_bootstrap_ratios` and the valid/undefined replicate
counts; it does not silently drop samples or replace the denominator. This is
especially common with few seeds or discrete differences. Finite intervals are
reported only when all bootstrap ratios are defined. No p-values or generic
small/medium/large labels are inferred from these simulation contrasts.

Nulls remain null: absence of a provider decision is not zero JS divergence,
absence of living cells is not a zero activated fraction, and well-mixed runs
do not have measured spatial distance readouts. `n_missing` makes the number of
undefined seed values visible in every summary.

## Readout semantics

Population and state fractions use committed snapshots. Cumulative events use
committed events, excluding blocked proposals. Kills, clearances, divisions and
apoptoses are attributed to the actor lineage; an additional killed-target
readout attributes deaths to the target lineage. JS at snapshot `s+1` describes
the provider-scored decisions in round `s`; it is undefined at the initial
snapshot and for unscored/local decisions.

Spatial contact means centroid distance within the declared radius. It uses
same-sample coordinates independently of the graph's neighbor cap, and does not
claim membrane contact. The nearest-neighbor mean includes living CD8 cells
with at least one living malignant cell in the same sample; `n` records this
denominator. Spatial metrics are undefined for well-mixed runs.

## Jev arms

`provider_factory` is an explicit Python extension point; it receives `arm`,
`bundle`, `output_dir`, `seed`, `steps` and `jev_weight`. It must return a provider
with its own durable journal. `arm_options` may override `mode` and `jev_weight`
by arm, permitting a local baseline with Jev arms. The factory can select shown
or hidden prior by the declared arm name. Provider model, prompt version and
prior mode remain recorded by each run's configuration. No factory is called
for a local arm, and the ensemble performs no automatic retries.

Before interpreting figures, complete annotation, rate/signature review, the
local null controls, numerical checks on the integrated kernel, and the Jev
replicate experiments. An implemented readout or control is not a positive
experimental finding. The plan's expression/context/spatial claims remain
conditional on those results.

## Lattice units

Cell-agent and lattice arms cannot share an ensemble comparison. Lattice runs
emit `occupied_bins`, `total_occupancy`, `cumulative_occupancy_gained` and
`cumulative_occupancy_lost`, plus provider divergence when available. Occupancy
uses abstract units; a full bin has the declared `lattice.capacity`. The sum
is not normalized by capacity and is not a count of cells. Lattice
analysis does not emit cell counts, cell state fractions, cell event counts or
cell-centroid proximity metrics under misleading names.

## Sensitivity

`virtualtissue_xenium.v2.sensitivity.run_sensitivity` executes a local parameter
sweep with the same simulation seeds at every parameter point. Its response is
the mean final readout across those seeds. Hazard overrides through
`action_priors.LINEAGE.ACTION` are absolute hazards per hour, before the normal
score and modifier operations; all other rates remain at their declared base
values. Other supported paths are `scoring.beta[.ACTION]`, `dynamics.KEY` and
`neighbor_radius_um`. Pooling alpha is rejected in local sweeps because it has
no effect there; a Jev alpha sweep needs a separate, budgeted experiment.

Bounds are explicit independent uniform distributions. Use `[low, high]` for
continuous values, or `{ "low": 1, "high": 20, "distribution": "uniform_integer" }`
for inclusive integer bounds. `method="lhs"` produces a Latin-hypercube screen
without claiming Sobol indices. `method="sobol"` constructs independent
pseudorandom A/B matrices and one-column-swapped hybrids: N × (D + 2) parameter
points, each repeated across the requested simulation seeds. The name denotes
the target indices; this implementation does not use a Sobol low-discrepancy
sequence.

The first-order estimate uses centered outputs and the total-effect estimate
uses squared differences, following the [Saltelli et al. estimator](https://doi.org/10.1016/j.cpc.2009.09.018)
and [SALib's documented formulas](https://salib.readthedocs.io/en/latest/_modules/SALib/analyze/sobol.html).
Confidence limits resample paired base rows. Failed/missing points block index
calculation; a constant response gives undefined indices, and estimates outside
[0, 1] remain visible rather than being clipped. Increasing parameter design
size and simulation replicates is required to assess convergence. A finite
seed mean approximates the stochastic model's expected response and retains
Monte-Carlo noise.

## Spatial reference checks

`virtualtissue_xenium.v2.spatial_stats` provides `cross_l`, `cluster_sizes` and
`held_out_envelope`. Cross-L requires one sample, two distinct lineages and an
explicit rectangular observation window in micrometres. It applies translation
edge weights based on rectangle-overlap area and reports K, L and L−r. See the
primary [spatstat Kcross documentation](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kcross.html)
and [translation correction definition](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/edge.Trans.html).
The window is never inferred from cell extrema. Out-of-window living cells
cause an error; radii cannot exceed half the shorter window side. Irregular
tissue masks and inhomogeneous intensity correction are not implemented.

Clusters are uncapped, same-sample connected components at a declared radius,
including singleton cells. Held-out envelopes use empirical ROI quantiles at
each radius, with matching lineages and correction. They are pointwise
descriptive comparisons, not simultaneous confidence bands or p-values. ROI
independence, the observation mask and annotation quality remain part of the
experimental design, not properties established by these utilities.

## Rebuilding figures

`python analysis/figures.py --ensemble outputs/ensemble --output figures/run1`
reads a completed ensemble without rerunning it and writes PNG/SVG plots plus a
figure manifest containing source hashes, generator hash, version and warnings.
Install the optional plotting dependency `matplotlib` in the analysis
environment. The output directory must be new; failed/excluded-arm ensembles
are refused. `--readouts` and `--lineages` select the plotted results. The dark
band represents bootstrap uncertainty in the median; the light band represents
the empirical seed envelope. Undefined values remain gaps, not zeroes.
