# Jev feedback and the 0.1.3 changes

The 0.1.2 reports describe almost exclusively quiescent tumour decisions and a
failure when Jev's chosen action trails the reported maximum by 0.01. The rounding
failure is reproducible offline. A prompt effect is a plausible explanation for
the proliferation reports, but the available breast results and model provenance
do not establish which wording change caused it. Recovery of tumour growth with
the new prompt has **not been validated with live Jev**.

## Declared tissue in the default prompt

`vtx-choice-3` uses the tissue explicitly recorded in `SOURCE.json`, for example
human breast carcinoma. It restores the concise t0 AUCell semantics used before
0.1.2 and removes the sentence about frozen scores from the prompt. Scores still
remain fixed during the simulation; daughters inherit them. This is a wording
change, not dynamic transcriptomics. Tissue context can introduce literature
priors and is an explicit modelling assumption, not a measured cell-level state.

The global recognition flag and excluded-action diagnostics remain out of the
request. CD8 questions remain in separate batches, preserving the 0.1.2 control
isolation. Menus, argmax choice handling, clocks, wake rules and biological rates
are unchanged. No probability threshold was lowered to force growth.

```bash
vtx own-data /path/to/outs data/Breast --sample-id Breast5K \
  --tissue "human breast carcinoma"
vtx annotation-template data/Breast Breast_clusters.json
# Review and replace every placeholder annotation before preparing regions.
vtx check-annotations Breast_clusters.json --data-dir data/Breast
```

For an already prepared bundle without tissue metadata, use the same override in
preflight and the new run. In the notebook, set `TISSUE` (blank uses `SOURCE.json`).

```bash
vtx check Breast_R01.zip --tissue "human breast carcinoma" --workdir check_Breast
vtx run Breast_R01.zip results --run-name Breast_R01_013 \
  --tissue "human breast carcinoma" --mode local
```

No tissue is inferred from filenames. Undeclared data use an explicitly labelled
generic human tumour-interface context. The resolved tissue, its source and the
full prompt state are recorded in `preflight.json`, `run_config.json` and the
episode. Changing this context requires a new run name; `vtx compare` rejects
comparisons with different prompt contexts.

## Rounding without invalidating input bundles

The provider tolerates a chosen probability up to `0.01 + 1e-9` below the maximum
**before normalization**. Larger discrepancies remain errors. The chosen option
is tied with the maximum in a compatibility copy and this copy is normalized for
the unchanged kernel. This also handles `{A: 0.50, B: 0.49}` with choice B: simply
changing both validators to 0.01 would still reject its normalized gap.

The exact returned values are retained in `reported_probabilities`, including
their original sum, and the answer records `probability_repairs`. The journal
also retains the whitelisted raw response. Tables, replay probabilities and cycle
diagnostics use the reported values, not the adjusted compatibility copy. A
normalizable but rejected answer is retained for audit with credentials redacted;
unserializable/nonfinite responses receive a failed journal entry without that
response. Failed or interrupted requests are never automatically resent.

All bundled runtime files remain byte-identical. Existing 0.1.2 input bundles
remain usable, but saved runs must be resumed with their original software.
Use a new run name for this code; do not edit fingerprints or reuse old journals.

## Interpreting quiet cells

`vtx check`, `vtx run` and notebook preflight now report:

* Acting cells, cells without any AUCell signatures, eligible t0 questions,
  selected t0 questions and questions deferred by the decision cap, per type.
* Active CD8 cells touching live acting tumour cells, or with one within the
  sensing radius (centroid distance), and CD8 cells without a t0 question.

A cell with no eligible question cannot choose a movement or effector program.
Passive mechanical displacement is still possible. Nearby tumour alone is not
a wake trigger in this kernel. Missing AUCell is not the scheduler trigger;
inspect RNA depth, QC and gene-set overlap separately. `--aucell-top-ranks 50`
is available as a preparation sensitivity analysis, not a guarantee of activity
or a calibration for smaller gene panels.

`run_report.json`, `vtx explore` and both results notebooks expose
`tumour_cycle_diagnostics`: reported P(ENTER_CYCLE) min/median/mean/max and histogram,
choices by mitotic quartile, and the fraction of **binary** cycle decisions
within 0.1 of 0.5. These are decision-weighted summaries, not independent-cell
estimates. Nonbinary menus have a separate margin against the best alternative.
Local preview diagnostics describe its rule-based answers, not Jev evidence.

## Prompt validation before release

The free snapshot tests pin the 0.1.2 comparator and current prompt. The ablation
tool holds cells, menus, model and question groups fixed. The semantics-only
variant changes only semantics, including keeping the version key fixed.
The six variants are `legacy_v1`, `v2`, `v2_legacy_setting`,
`v2_legacy_semantics`, `v2_with_excluded` and the actual shipped `v3`.

```bash
python tools/jev_prompt_ablation.py R03.zip ablation_dry \
  --reference R03_jev04_t0_tumour_reference.json
# Paid experiment, only when explicitly authorized by the operator:
python tools/jev_prompt_ablation.py R03.zip ablation_paid \
  --reference R03_jev04_t0_tumour_reference.json --send --max-requests 24
# Smaller candidate canary: 20 cells, two repeats, normally four requests:
python tools/jev_prompt_ablation.py R03.zip canary_v3 \
  --variants v3 --send --max-requests 4
```

The reference JSON and R03 bundle are external evidence, not included in the
repository. Dry runs write payloads and a fingerprinted plan without contacting
Jev. Actual counts depend on payload sizes; a send exceeding the cap is refused
before the first request. Both byte and question limits are enforced. Requests
are journaled before sending. Existing journals refuse a new send, including
after interruption, and the first failure stops the experiment.

Before declaring the proliferation regression resolved or tagging a release,
record the authorized live ablation/canary and review these predefined criteria:

1. P(ENTER_CYCLE) increases across mitotic quartiles; the top-quartile mean is
   within 0.10 of the saved legacy reference.
2. Repeated legacy questions reproduce the saved reference within their repeat
   spread. Otherwise investigate model drift before attributing changes to text.
3. Non-CD8 requests remain identical between recognition arms. Test breast tissue
   separately; an R03 result cannot establish recovery for Breast 5K.

## Scope retained

This patch does not implement probability sampling, nearby-tumour CD8 wake rules,
region ranking by CD8 contact, or non-tumour tissue profiles. These would change
the model and require separate validation. `own-data` still registers extracted
folders; URL/ZIP/GEO import is not implemented here. The team's
`fetch_xenium_outs.py` was not present in `Revision_v1.2` and was not incorporated.
