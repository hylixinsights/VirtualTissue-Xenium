# VirtualTissue-Xenium v3

A research simulator using supplied expression, existing cell labels and segmentation polygons. Version **3.0.0a4** is an alpha prerelease for testing.

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

For the prepared `V3_teste` folder, use its top-level **START_TEST.command**. You do not need to read the internal documents to run the test.

A passing result applies to the exact recorded source snapshot. Keep the scientific scope: expression influence, contextual Jev preferences and selected structural statistics must still be evaluated on appropriate real inputs. Existing labels are accepted; H&E and a new annotation review are not required. The model does not establish a proven prediction of tissue seven days later.

**Compatibility:** create new prepared inputs and runs with v3. Changed kernel, prompt and ontology fingerprints prevent resuming older runs with this version; retain their original pinned checkout. The legacy `vtx` interface remains available; v3 uses `vtx2`.

## Choose an ROI from your object

For a raw JSON object containing cells, exact labels, counts, coordinates and
optional polygons:

```bash
vtx2 roi-catalog object.json regions.html --cells-per-roi 1000
# Open regions.html, inspect composition and explicitly download a selection.
vtx2 select-roi object.json selection.json selected.json --context context.json
vtx2 expression-audit selected.json expression_audit.json
vtx2 cost selected.json --steps 48
```

Use `examples/v3/context_user_object.json` as an explicit starting profile and
replace its tissue/organism declarations. It uses `tiered_panel_auc`: sufficient
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
vtx2 inspect outputs/v3/demo/control.json
vtx2 ensemble outputs/v3/demo/arms.json outputs/v3/demo/ensemble \
  --horizon-hours 24 --seeds 0 1 2 3 4
make figures ENSEMBLE=outputs/v3/demo/ensemble FIGURES=outputs/v3/demo/figures
```

This creates 270 synthetic cells with heterogeneous expression, a β=0 arm and a within-type expression-shuffle arm. Readouts distinguish bootstrap confidence intervals of the median from the distribution across seeds. A completed run includes decisions, events, reports, frames and `replay.html`; the ensemble adds long-form CSVs, paired differences and provenance.

## Refresh the replay of an existing run

The original dark replay layout is restored with tissue indicators, the malignant-cell chart, event feed and cell inspector. The v3 player retains its bounded frame cache. After updating your checkout and package in the environment used for the run:

```bash
git pull --ff-only
pip install -e .
vtx2 replay outputs/v3/ensemble_990_1350_48h/runs/recognition/seed_4
```

Open that run's `replay.html`. This command reads the saved frames and report and replaces only the replay HTML and its presentation sidecars. It does not rerun the simulation, call Jev or change decisions, readouts or scientific exports. Keep `replay.html` beside its `replay_frames` folder when moving it. Repeat the command for other arms/seeds as needed; `--output another/replay.html` writes a separate presentation copy. V3 uses `frames.jsonl`; v2 uses `frames.json`.

## Use measured data

Supply existing `cell_type`/`lineage` annotations, exact gene symbols and measured coordinates in micrometres, recording their origin. For raw Xenium:

```bash
vtx2 prepare-xenium DATA_DIR context_v3.json prepared.json \
  --annotations existing_labels.csv --lineage-map lineage_map.csv \
  --center 1200 1600 --max-cells 1000
vtx2 inspect prepared.json
vtx2 run prepared.json outputs/sample_local --steps 48 --seed 21
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
vtx2 cost prepared.json --steps 48
vtx2 run prepared.json outputs/sample_jev --mode jev --jev-prior hidden \
  --steps 48 --seed 21 --max-http-attempts 500
```

Set `TYPESAFE_API_KEY` or `JEV_API_KEY` for live use and choose a budget from the cost estimate. The example budget is illustrative, not guaranteed sufficient. Hidden prior forces α=1. Compare shown/hidden/scrambled-context arms with independent journal replicates; `jev-audit` summarizes replicate and order variance without sending requests. Provider preferences are not calibrated cell-fate rates. Journals are configuration-bound; failed or unresolved attempts are never automatically resent.

## Validation and reproducibility

```bash
make test
make clock
vtx2 --help
```

The CI runs v2/v3 tests, exports coverage and produces the restricted local clock-validation figure. `make figures` regenerates PNG/SVG outputs with a hash manifest. No clinical or treatment-efficacy results are claimed by these software checks.

The immutable physical engine remains under `src/virtualtissue_xenium/engine` with the legacy `vtx` command. Its bundle format differs from `vtx2`. The general v2/v3 kernel now resolves competing death proposals without double-counting; old run fingerprints require their original checkout for resumption.

MIT software license. Dataset and signature source terms/attribution remain attached to their manifests. Cite the software using `CITATION.cff` and cite the actual data, biological resources and provider model used.
