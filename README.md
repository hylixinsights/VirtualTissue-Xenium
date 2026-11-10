# VirtualTissue-Xenium v5

A research simulator using supplied expression, existing cell labels and segmentation polygons. **5.0.0** is the current primary software release. The simulation mechanisms, kinetic assumptions, cell rendering and dark replay are those of v4, unchanged. V5 improves which measured evidence reaches the per-cell rulebook, what the Jev model is told about that evidence, and how the quality of its decisions is measured. See [what v5 adds](#what-v5-adds). See the [v4 release notes](docs/v4/RELEASE_NOTES.md).

**Scientific status:** default kinetic rates are **DRAFT assumptions, not calibrated biological rates**. Verified citations provide context, not validation of the numerical values. Passing software tests does not validate biological predictions. See [rate review](docs/v3/RATE_REVIEW.md) and [limitations](LIMITATIONS.md).

## Start here

On a Mac, install Python **3.12 or 3.13**, Node.js **22 or later**, and Git. Clone the repository in Terminal, then run the test:

```bash
git clone https://github.com/hylixinsights/VirtualTissue-Xenium.git
cd VirtualTissue-Xenium
./START_TEST.command
```

You can also double-click **START_TEST.command** inside the cloned folder. Use a Git clone for this workflow: the acceptance check needs Git history and will not work from GitHub's Download ZIP. The first run may install test dependencies in `.venv-test`; keep that folder with the checkout. Python and Node.js themselves are not installed by this script.

The test opens one **RESULT.html** page with PASS/FAIL, an interactive synthetic ROI catalogue and simulation replay. Results are saved beside the checkout in `VirtualTissue-Xenium-test-results`. No Jev calls or GitHub publishing occur. The automated simulation uses a fixed synthetic fixture; the catalogue starts without a user selection.

Linux users with the same Python/Node.js prerequisites and `zsh` can run `zsh START_TEST.command --no-open`, then open the result page manually.

Existing prepared `V3_teste` folders retain their own top-level **START_TEST.command** and pinned package. Use the current repository checkout for v5.

A passing result applies to the exact recorded source snapshot. Keep the scientific scope: expression influence, contextual Jev preferences and selected structural statistics must still be evaluated on appropriate real inputs. Existing labels are accepted; H&E and a new annotation review are not required. The model does not establish a proven prediction of tissue seven days later.

**Commands and compatibility:** use `vtx5` for the current workflow. `vtx4` and `vtx2` remain aliases for the same workflow; the legacy physical engine keeps `vtx`. V5 retains simulation schema **3.0**, the existing ontologies and the `docs/v3` and `examples/v3` paths. Existing contexts without a v5 evidence policy or the new optional sample-context fields produce the same frames, decisions, configuration and Jev requests as v4.0.0. The test suite verifies exact bytes against the pinned immutable v4 commit under the same interpreter and platform; see [verification details](docs/v5/VALIDATION.md).

## What v5 adds

- **Evidence that reaches the rulebook.** On a targeted panel the packaged marker sets can be unassayed, which leaves every expression multiplier at exactly 1. `vtx5 inspect` now reports, per action, how many cells have usable evidence and warns when an action is label-only. An optional, declared [evidence policy](docs/v5/EVIDENCE_POLICY.md) states which evidence may modulate which already-permitted action for which lineage, and can use per-cell scores computed outside the simulator, imported with their provenance.
- **English sample-context questions.** `vtx5 context-init sample_context.json --template examples/v3/context_user_object.json` asks about tissue, disease, anatomical region, developmental stage, relevant cell states and secretory populations. These declarations reach the rulebook and Jev; only separately declared restrictions change permitted actions. See [sample context](docs/v5/SAMPLE_CONTEXT.md).
- **A clearer prompt for Jev**, used only with an evidence policy: each value comes with its reference size and detected genes, definitions are sent once per request, and measured evidence is distinguished from simulated state.
- **Tools that measure Jev's decisions.** `vtx5 jev-probe` asks a frozen set of t0 questions under controlled manipulations and compares the answers with identical repeats. Every Jev run reports how much event probability the kernel declared and how much was applied.
- **A smaller, lossless decision log.** Frozen t0 evidence is stored once per run; `vtx5 decisions` restores the v4 layout exactly.
- **An opt-in pooling rule**, `--jev-pooling event_share`, which keeps the kernel's event probability and lets Jev choose only among events. The default is the v4 rule.

Details and limits are in the [5.0.0 release notes](docs/v5/RELEASE_NOTES.md) and in [Jev decisions](docs/v5/JEV_DECISIONS.md). None of this validates a biological prediction; the scientific status above is unchanged. Saved v3 trajectories can refresh their replay without a new simulation. Resume still requires matching configuration and fingerprints; retain pinned checkouts for incompatible historical runs.

## Choose an ROI from your object

For a raw JSON object containing cells, exact labels, counts, coordinates and
optional polygons:

```bash
vtx5 roi-catalog object.json regions.html --cells-per-roi 1000
# Open regions.html, inspect composition and explicitly download a selection.
vtx5 select-roi object.json selection.json selected.json --context context.json
vtx5 expression-audit selected.json expression_audit.json
vtx5 cost selected.json --steps 48
```

Create your context with `vtx5 context-init context.json --template examples/v3/context_user_object.json`
or edit that template directly. Review its treatment declarations, recognition and numerical assumptions.
It uses `tiered_panel_auc`: sufficient
cell references when available, then a bounded weak within-cell panel heuristic
for rare types. Missing panel support stays unavailable. No rule forces every
cell or every dataset to show a positive effect. Read the
[scoring policy](docs/v3/EXPRESSION_POLICY.md) before interpreting it.

AnnData imports may use `--polygon-key` for a named `uns`/`obs` field, or an
external exact-ID polygon/GeoJSON file through `--polygons`. See the
[ROI workflow](docs/v3/ROI_WORKFLOW.md). Existing labels need an explicit mapping
to supported action lineages; no new cell-type annotation is inferred.
For a raw object without `lineage`, use the same `--lineage-map labels.csv` on
both catalogue and selection commands. This retains `cell_type` and records the
mapping. Try the [synthetic polygon example](examples/v3/README_OBJECT_DEMO.md)
to inspect the workflow before importing real data.

## Run the synthetic example

Python 3.12 or 3.13, from this checkout:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e '.[xenium,test,spatial,analysis]'
python examples/v3/make_demo.py --output outputs/v3/demo
vtx5 inspect outputs/v3/demo/control.json
vtx5 ensemble outputs/v3/demo/arms.json outputs/v3/demo/ensemble \
  --horizon-hours 24 --seeds 0 1 2 3 4
make figures ENSEMBLE=outputs/v3/demo/ensemble FIGURES=outputs/v3/demo/figures
```

This creates 270 synthetic cells with heterogeneous expression, a β=0 arm and a within-type expression-shuffle arm. Readouts distinguish bootstrap confidence intervals of the median from the distribution across seeds. A completed run includes decisions, events, reports, frames and `replay.html`; the ensemble adds long-form CSVs, paired differences and provenance.

## Refresh the replay of an existing run

The original dark replay layout is restored with tissue indicators, the malignant-cell chart, event feed and cell inspector. The player retains the bounded frame cache introduced in v3.

Cells also use the original renderer: rounded cell bodies clipped against neighboring bodies, sizes in tissue coordinates, visible nuclei, context transparency and event glows. **Original cells** uses the original transition from supplied t0 segmentation to illustrative simulated silhouettes. Choose **Provided outlines** to inspect exact saved polygons, or **Centroids** for position markers. These display choices do not change the recorded coordinates, interaction graph, segmentation evidence or model dynamics.

Update your checkout and package in the environment used for the run:

```bash
git pull --ff-only
pip install -e .
vtx5 replay outputs/v3/ensemble_990_1350_48h/runs/recognition/seed_4
```

Open that run's `replay.html`. This command reads the saved frames and report and replaces only the replay HTML and its presentation sidecars. It does not rerun the simulation, call Jev or change decisions, readouts or scientific exports. Keep `replay.html` beside its `replay_frames` folder when moving it. Repeat the command for other arms/seeds as needed; `--output another/replay.html` writes a separate presentation copy. Schema-3.0 runs use `frames.jsonl`; schema-2.0 runs use `frames.json`.

## Use measured data

Supply existing `cell_type`/`lineage` annotations, exact gene symbols and measured coordinates in micrometres, recording their origin. For raw Xenium:

```bash
vtx5 prepare-xenium DATA_DIR context_v3.json prepared.json \
  --annotations existing_labels.csv --lineage-map lineage_map.csv \
  --center 1200 1600 --max-cells 1000
vtx5 inspect prepared.json
vtx5 run prepared.json outputs/sample_local --steps 48 --seed 21
```

The centre above is only an example. Prefer the composition catalogue for user selection. `roi-candidates` and `annotation-agreement` remain optional tools for studies using the earlier front/core or independent-review design. The bundled historical unknown-only annotations do not supply usable cell identities.

- [Annotation workflow](docs/v3/ANNOTATION.md)
- [Ontology, rate and signature provenance](docs/v3/ONTOLOGY.md)
- [Native platforms and measured-unit contracts](docs/v3/PLATFORMS.md)
- [Readouts, ensembles, sensitivity and spatial statistics](analysis/README.md)
- [Legacy own-data workflow](docs/v2/YOUR_DATA.md) and [Xenium dataset guide](docs/v2/DATASETS.md)

For scRNA/snRNA without measured XY, use well-mixed interactions. Unsegmented Visium/Stereo bins require the separate lattice kernel and declared purity; occupancy is not a cell count. New native adapters have format tests and still require validation on public real datasets.

## Jev experiments

```bash
vtx5 cost prepared.json --steps 48
vtx5 run prepared.json outputs/sample_jev --mode jev --jev-prior hidden \
  --steps 48 --seed 21 --max-http-attempts 500
```

Set `TYPESAFE_API_KEY` or `JEV_API_KEY` for live use and choose a budget from the cost estimate. The example budget is illustrative, not guaranteed sufficient. Hidden prior forces α=1. Compare shown/hidden/scrambled-context arms with independent journal replicates; `jev-audit` summarizes replicate and order variance without sending requests. Provider preferences are not calibrated cell-fate rates. Journals are configuration-bound; failed or unresolved attempts are never automatically resent.

## Validation and reproducibility

```bash
make test
make clock
python tools/v4_equivalence.py compare --same-platform-v4
vtx5 --help
```

The equivalence check requires full Git history; run `git fetch --unshallow` first in a shallow clone. The fixed macOS-recorded hash file remains a historical, environment-specific artifact. CI exposed a mismatch when comparing Linux outputs to that snapshot; the same-platform check retains strict equality without changing numerical logic.

The CI checks the current replay and the retained v2/v3 simulation contracts, exports coverage and produces the restricted local clock-validation figure. `make figures` regenerates PNG/SVG outputs with a hash manifest. No clinical or treatment-efficacy results are claimed by these software checks.

The immutable physical engine remains under `src/virtualtissue_xenium/engine` with the legacy `vtx` command. Its bundle format differs from `vtx5`/`vtx4`/`vtx2`. The general kernel retains the v3 fix that resolves competing death proposals without double-counting; incompatible old run fingerprints require their original checkout for resumption.

MIT software license. Dataset and signature source terms/attribution remain attached to their manifests. Cite the software using `CITATION.cff` and cite the actual data, biological resources and provider model used.
