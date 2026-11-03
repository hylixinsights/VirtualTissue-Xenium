# Paper scope and reproducible evidence

The [confirmed scope](../docs/v3/SCIENTIFIC_SCOPE.md) uses the supplied object's
existing labels and segmentation. It requires neither H&E nor a new manual
annotation campaign. ROI selection belongs to the user and is recorded by IDs
and source hashes. The three claims below are hypotheses to evaluate, not results
already obtained. No seven-day simulation will be described as a demonstrated
prediction of the real tissue seven days later.

## Evidence and decisions

| Hypothesis | Evidence to collect | Interpretation if evidence is insufficient |
|---|---|---|
| Expression changes simulated dynamics beyond labels | Fixed-state probability audit; local beta-zero, within-type and global expression controls; outcome contrasts over multiple seeds and independently declared permutations; empirical/weak/unavailable score coverage | Report only the measured mechanism or bounded outcome contrast. A probability change alone is not evidence of improved biological fidelity. Do not tune beta until a contrast appears. |
| Jev responds to context beyond prior anchoring | Identical measured states and menus; shown/hidden prior and swapped descriptive context; repeated identical-payload journals and batch-order audit; compare condition differences with replicate noise; record usage | Report anchoring, variability or inconclusive context effects. Beta-zero is not an expression-free Jev arm because the provider still receives evidence. Local clock validation does not validate Jev's time scale. |
| Selected structural statistics remain compatible with reserved real regions | Reserve reference regions before tuning; specify statistics and valid observation windows; compare trajectories with static input, reference distributions and controls; report all tested regions | Without reserved regions the claim is unassessed. Static spatial compatibility does not establish future prediction or generalisation across patients. |

For each study, freeze dataset/sample identifiers, existing label provenance and
explicit lineage mapping, ROI selection, signatures/orthologs, scoring policy,
QC, geometry handling, rates, context, horizon/step, readouts and analysis before
examining the main contrasts. Reserve structural reference regions at that stage.
The exploratory catalogue may be browsed freely; regions used to tune a model
cannot later serve as untouched validation data.

Use at least ten simulation seeds for the headline figures retained from the
plan, and assess whether estimates stabilise. Repeat null permutations rather
than treating one shuffled object as the complete null distribution. Provider
replicates need distinct journals and matched states; simulation seeds,
permutations, provider noise and biological samples must remain separately
identified. These are design commitments, not claims that any experiment has
already been run or sufficient statistical power has been established.

Inspect supplied geometry/QC and empirical versus weak panel evidence. Test
dependencies on label uncertainty, depth/QC policy, centroid versus boundary
geometry, rates, movement/crowding and simulation step when relevant. Such tests
measure robustness under specified perturbations; they do not recover true
segmentation or annotation error rates. Keep unsupported rates `assumed`.
Bibliographic identities, Ki67 scores and sensitivity analyses do not by
themselves identify an hourly tissue-specific division hazard.

Limitations must affect interpretation, not merely appear in a paragraph. If a
result depends strongly on an assumption, restrict the claim to that setting;
if a test is unavailable, omit the unsupported claim or report it unassessed.
Negative and inconclusive findings remain part of a methods/resource article.
Neither unchanged seed envelopes nor synthetic software fixtures establish
biological equivalence, absence of information or treatment efficacy.

## Rebuilding figures and retaining provenance

Figures are regenerated from complete ensemble outputs:

```bash
make figures ENSEMBLE=/path/to/completed/ensemble FIGURES=/path/to/new/figure-directory
```

The renderer writes PNG/SVG figures and source manifests. It distinguishes the
bootstrap confidence interval of the median from the empirical seed envelope.
Keep prepared inputs, selections, parameter/signature/ontology versions, hashes,
controls, failed runs, provider identities/journals, cost ledger and analysis
settings. The renderer covers ensemble readouts and the clock figure; final
spatial/sensitivity and study-specific panels still need their corresponding
experiments and rendering code.

See [analysis definitions](../analysis/README.md), [input/ROI workflow](../docs/v3/ROI_WORKFLOW.md),
[expression policy](../docs/v3/EXPRESSION_POLICY.md) and [limitations](../LIMITATIONS.md).
The existing synthetic demo and clock plots are engineering evidence only.
