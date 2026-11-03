# Frozen scientific scope — November 1, 2026

Explicit decision by the author on this date. This document supersedes the
requirements for H&E, mandatory independent review and automatic front/core
selection in the historical plan. The originals remain preserved as historical
records.

## Three claims to evaluate

1. **Individual expression changes simulated dynamics beyond cell labels.**
2. **Jev preferences respond to context and do not merely reproduce the prior.**
3. **Simulated structures remain compatible with specified statistics from real regions held out for evaluation.**

These are scientific objectives subject to results, not conclusions already
demonstrated. Freezing these questions does not require positive findings. The
paper will present detected effects, absence of detectable effects and inconclusive
results. **A proven prediction of the tissue seven days later will never be
claimed.** Simulated time, assumed parameters and initial structure are not
observations of the future.

**Evidence status:** the new features have passed synthetic implementation tests.
The three claims have not yet been evaluated on real data, and the live Jev
experiments in this protocol have not yet been performed. Synthetic examples,
constructed payloads and lists of limitations do not complete those steps.

## Accepted inputs and user workflow

The project is an exploratory simulator distributed through GitHub, intended for
objects supplied by other researchers. It uses expression, existing labels,
segmentation polygons and declared coordinates. It does not require H&E or
organizing independent annotation review. Labels may come from automated
pipelines; their origin and verification status must accompany the results.

Polygons are **observed segmented units**, not true membranes. Original contours
are retained; geometric measurements check computational consistency without
certifying biological quality. Geometric problems must be shown to the user.
Invalid contours are not silently repaired. Where allowed, the declared centroid
is used with a warning; without a usable position, the object cannot enter a
spatial simulation.

The user sees ROI candidates, a map, composition, cell counts, the fraction of
unknown labels, expression information and QC. **The user chooses the ROI.** No
region is selected automatically because it has more tumour contact or produces
a favorable result. The choice, IDs and source are recorded.

Malignant, CD8 and other class labels are input declarations that determine model
permissions. The software does not certify malignancy or antigen recognition.
Missing labels do not authorize inventing them; `unknown` remains conservative.
Independent review, when available, is additional information.

## Expression as an explicit mechanism

Expression must have a traceable computational path to weights, conditional menus
and decision probabilities. Its effect depends on the available evidence and the
actions permitted for the lineage. No signal will be fabricated for a panel
without coverage, a cell without counts or a signature without useful variation.

In addition to the previous strict method, the user-object profile offers
`tiered_panel_auc`: a sample/type/depth reference when sufficient; a within-sample
type reference when strata are small; and, for rare types, within-panel rank
evidence against an analytical reference of random gene sets. The last tier is
explicitly a **weak calibration heuristic** with a bounded multiplier. It is not
a population percentile, does not guarantee removal of depth bias and does not
preserve the population mean of rates. Controls use the same rule. For this
method, `reference_fallback` is `cell_type`; an explicit `null` declaration is
rejected. To retain the strict policy of no score when references are insufficient,
use `depth_conditioned_auc` with `reference_fallback: "null"`.

Reliable labels and large groups are not invented by pooling across lineages.
Signatures with insufficient coverage remain unavailable; expanding gene sets or
changing thresholds requires justification independent of the result.

An audit at identical states measures available evidence, expression-sensitive
actions, multipliers, probabilities with/without expression and shuffle controls.
A changed probability does not require a random realization to choose a different
action, nor does it prove improved biological fidelity. In Jev mode, available
evidence is sent with its basis/calibration; the model is not instructed to
artificially produce a positive contrast. β=0 removes modulation of the local
prior but does not by itself remove the direct evidence received by Jev.

## Minimum evidence for the paper

| Claim | Required evidence | What it does not establish |
|---|---|---|
| Expression changes dynamics | Probability audit at identical states, followed by effects on trajectory outcomes; β=0, within-type and global expression shuffle controls; multiple seeds and permutations with separately identified uncertainty | Better biological prediction, protein activity or absence of information when bands overlap |
| Jev responds to context | Identical states/menus, shown/hidden prior, swapped context, repeated responses; JS between conditions compared with JS between identical-payload replicates | A calibrated physical clock or an effective treatment |
| Structural compatibility | Reference regions held out before tuning; predefined statistics/window; comparison with a static reference and controls | The true future, reconstructed real cells or generalization across patients |

Exploratory mode allows any ROI to be chosen. For the paper's evaluation, held-out
regions cannot be reused to tune scores or rates or to select the main result.
If no held-out reference is available, the third claim remains unevaluated.
Structural analysis is not required in every software user's session.

## Protocol for completing the three evaluations

### 1. Freeze the design before the main results

Create a study manifest with IDs/hashes of the objects and chosen ROIs,
sample/section/individual of origin when known, the existing label map, polygon
provenance and QC/fallback policies. Do not fill in missing identities by
assumption. Reserve regions for structural evaluation before tuning parameters;
their expression and outcomes cannot be used to calibrate scoring in the
development arms. Prepare the chosen ROI again using the declared reference
scope, recording the scores actually used.

In the same manifest, fix the signatures, scoring method, β per action, panel
fallback bound, rates, geometric rules, mixing α, Jev version/prompt, time step,
horizon, population cap and failure policy. Declare primary outcomes, units,
lineages, final time point or curve summary, and the **smallest relevant
difference** for each contrast. Choose these differences based on the question
and outcome units before seeing separation between arms; there is no universal
JS, β or effect-size threshold that validates the study.

Record simulator seeds, permutation seeds and independent provider replicates
separately. Ten seeds for main results, ten realizations of each evaluated shuffle
and three independent journals per Jev condition are operational starting points
inherited/adapted from the plan, not guarantees of precision or power. Set a
precision target and an expansion rule before confirmatory analysis. Do not count
simulated cells, reused response pairs or ROIs from the same individual as
independent individuals.

Clearly separate development results from held-out evaluation. Adjustments made
after observing a result are exploratory and require a new held-out evaluation
to support a confirmatory conclusion. If the budget cannot support the planned
design, narrow or defer the claim and document the change.

### 2. Demonstrate expression's influence on the mechanism and trajectories

Run `expression-audit` first, without Jev calls. For each type/signature, report
panel coverage, unavailable scores, empirical references and use of the weak
heuristic; for each action, report gates, hazard and probability before and after
controls. The current audit compares the same initial state, not all future
states. An available score that does not change a compatible action must remain
visible as such.

Then run local ensembles with original data, β=0, expression shuffled within
type/sample, and expression shuffled globally. Use the same simulation seeds for
each contrast and multiple recorded permutations. Global shuffling can change QC
and mix type profiles; its effect therefore does not isolate only the spatial
location of expression. Singleton types do not permit an informative within-type
permutation. When the hypothesis depends on them, include coordinate controls,
well-mixed conditions, recognition disabled, or treatment text without a modifier.

Compare the trajectory outcomes defined in the manifest, in their original units,
with simulation and permutation uncertainty identified. Show `dz` only as a
complementary statistic for paired differences between seeds; degenerate
denominators or bootstraps remain undefined. Repeat under predeclared β values
and panel bounds, including β=0, and under scoring/QC alternatives justified
during development. Do not increase β or change genes to obtain a positive
contrast.

A change only in initial probabilities supports the narrow conclusion of
**functional influence on local decisions**. The dynamics claim requires an
effect on trajectory outcomes. An interval beyond the predefined relevant
difference can support that influence in the evaluated scenario; overlap with
the relevance bounds makes the result inconclusive. Practical equivalence
requires sufficient precision within the previously defined bounds. Overlap of
seed envelopes alone does not demonstrate absence of information; see the
principle of predefined bounds in [Lakens, 2017](https://journals.sagepub.com/doi/full/10.1177/1948550617697177).

### 3. Separate context, expression, prior and Jev variability

Construct a recorded set of states/questions, keeping labels, geometry, dynamic
state, menus, batch order and batch composition identical in the main contrasts.
Collect responses in new, independent journals; reusing the cache reproduces a
response rather than creating a replicate. Audit repetition and order using ten
repetitions of an identical payload and ten permuted orders, separately for shown
and hidden priors, as an operational provider pilot.

Compare the original context and a declared swapped context at the same states,
with the prior shown and hidden. To attribute expression's influence on Jev, add
separate ablations of direct evidence, the numerical prior and both. β=0 alone
does not remove direct evidence. Also control differences in the opaque context
hash so that the comparison does not change metadata beyond the intervention of
interest. These response collections and analyses still need to be performed;
the local audit does not replace them.

Calculate JS between responses from different conditions and JS between
replicates of the same payload, always over corresponding menus. Compare the
distributions and contrasts of these quantities, accounting for repetition by
question/journal. Probability variance does not have the same units as JS and
must not be used as a direct numerical threshold. Multiple pairs constructed
from the same responses are dependent, not additional replicates. A different
response to context demonstrates sensitivity to the presented context, not that
the context or response is biologically correct.

Study α on a predeclared grid, for example the original grid
`0, 0.25, 0.5, 0.75, 1`, in shown-prior mode. α=0 is the local reference and does
not require new calls to run local dynamics; hidden forces α=1 and does not
constitute an α sweep. At frozen states, mixing can be examined by reusing
recorded responses; this does not reproduce trajectories that would visit
different states under another α. Jev trajectories require independent collection
and their own study of time-step dependence. The local clock does not calibrate
them.

### 4. Evaluate structural compatibility in held-out regions

Before evaluation, fix which real regions represent the compartment, their
observation windows, edge rule, radii and statistics: for example, cross-L and
cluster-size distributions. Document the relationship between held-out and
simulated regions; choosing only favorable regions invalidates the evaluation.
ROI, section and individual are different levels of replication.

Compare the initial state, a static reference that retains t0, and simulated
trajectories at the registered horizon, including appropriate geometric controls.
Compatibility also achieved by the static reference does not demonstrate that
dynamics are necessary or more faithful. Do not redefine the window from the
extreme positions of simulated cells; fix how cells leaving it will be treated.

The current tools implement rectangular windows, translation edge correction and
pointwise empirical envelopes. Irregular masks, inhomogeneous-intensity
corrections and simultaneous curve inference require additional analysis before
they can be claimed. Current proximity/CD8–malignant readouts use **centroids**,
regardless of the neighbor cap and whether the action graph uses polygon-boundary
distance; they are not membrane-contact measurements or a direct validation of
that graph.

A descriptive comparison must be presented as descriptive. For a formal decision
about an entire curve or many outcomes, predefine a procedure that addresses
multiplicity and dependence; do not declare success because some curves fall
within pointwise bands. The distinction between pointwise bands and global
assessment is developed by [Myllymäki et al., 2017](https://academic.oup.com/jrsssb/article/79/2/381/7041937).
Without suitable held-out regions, claim 3 remains unevaluated. Even with a
favorable result, the conclusion is compatibility of the selected statistics,
not observation of the future tissue.

### 5. Complete robustness, cost and the reproducible package

Check time-step dependence in the complete scenario over the same horizon before
interpreting trajectories. Identical seed labels at different time steps do not
automatically constitute a random stream coupled by cell-hour. Do not disable the
hazard guard to make a coarse time step finish.

Declare parameter ranges and distributions for sensitivity analysis. LHS is for
exploration; the implemented Sobol indices require A/B/ABᵢ and
`N × (D + 2) × number of seeds` runs, with independent distributions. α does not
change the local model and requires a separate Jev study. Check convergence in
design size and number of seeds; sensitivity analysis does not turn an assumed
rate into a measured one.

The current ensemble separates seeds and arms but does not automatically estimate
a hierarchy of journals, permutations and individuals. That analysis must be
specified and implemented/verified for the chosen design. Do not flatten these
dimensions as if all were independent seeds. The current rule excludes an arm
from aggregation if any seed is incomplete; preserve failures and costs as well.
The journal-audit tool reads existing records and does not automatically perform
collection or all the JS contrasts described above.

Deliver input/ROI manifests, frozen contracts, coverage and failure tables,
results by seed/permutation/journal, observed costs, analysis/figure scripts and
a final table for each claim: evaluated with support, effect below the relevant
bound, inconclusive or unevaluated. Interpretation must follow the results; a
complete list of limitations does not replace the tests.

## Robustness possible without external biological ground truth

Record the fraction of missing/invalid polygons and compare boundary distance
with centroid distance; vary the radius and edge handling. These exercises
measure dependence on geometric assumptions. They do not detect every
segmentation error or correct transcripts assigned to the wrong cell.

Report the fraction of null scores, small references and panel-heuristic use;
repeat analyses under predeclared scoring alternatives. Examine sensitivity to
QC and uncertain labels without optimizing β to produce separation. Use
controlled input perturbations as robustness analysis, not as a validated
estimate of the real error distribution.

Testing the local clock does not validate the Jev clock. While Jev preferences
are applied directly per step, describe trajectories in simulated units and test
time-step dependence; do not equate simulated duration with proven temporal
prediction. Uncalibrated rates and mechanisms remain assumed.

## Limitations that must accompany the results

Segmentation and labels may be wrong; without H&E/review, these identities have
no independent verification. The panel is incomplete, scores are approximations,
rare types use limited information, expression remains fixed at t0, and subsequent
dynamics are generated by assumed rules. Moving polygons undergo rigid
translations; daughter shapes are declared approximations, not observed future
segmentations. Spatial and expression errors may be correlated through
segmentation.

Listing these limitations bounds interpretation. The corresponding tests are
still needed to support the three claims. A methods/resource paper may report
negative results and exploratory use without promising a true representation of
the tissue or guaranteeing editorial acceptance.
