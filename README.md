# VirtualTissue-Xenium v2

**Build an inspectable virtual-tissue experiment from an annotated cell object, expression data and a declared biological context.**

VirtualTissue-Xenium is part of Helder Nakaya's [Virtual Tissue initiative](https://virtualtissue.org/). This repository explores tissues as communities of individual cell agents whose decisions are constrained by explicit biological rules.

Version **2.0.0a1** is a working research preview. Each cell has an expression-informed guidebook, lineage-compatible actions and simulated states. Jev evaluates permitted actions; a seeded sampler chooses actions; the kernel commits a complete round from a shared snapshot. Disease, treatment, mutation, tissue and sampling context are explicit inputs.

**Start with the [ordered workflow and task plan](docs/v2/WORKFLOW.md).** The [v0.1.3 code review](docs/v2/CODE_REVIEW_V0_1_3.md) explains the changes. The [validation plan](docs/v2/VALIDATION_PLAN.md) separates software checks from biological validation.

**V2 is the default workflow on `main`.** New users can follow the [step-by-step guide for their own data](docs/v2/YOUR_DATA.md), including annotation requirements, a complete context example, preparation, simulation and results. No development branch selection is needed.

| Your starting point | Start here |
|---|---|
| Your annotated Xenium, scRNA-seq or single-nucleus object | [Run your own data](docs/v2/YOUR_DATA.md) |
| Raw Xenium output with reviewed annotations | [Xenium files, annotation templates and ROI selection](docs/v2/DATASETS.md#2-create-reviewed-annotations-and-select-a-region) |
| Try a small example without an API key | [Synthetic quick start below](#first-run-no-dataset-or-api-key-needed) |
| Prefer a notebook | [V2 workflow notebook](notebooks/V2_Context_Guidebook_Workflow.ipynb) · [Open in Colab](https://colab.research.google.com/github/hylixinsights/VirtualTissue-Xenium/blob/main/notebooks/V2_Context_Guidebook_Workflow.ipynb) |

The guidebook currently uses **transparent heuristic preferences**, rather than probabilities learned from observed cell fates. Jev scores are also not calibrated cell-fate probabilities. Expression stays fixed at baseline; the simulator evolves discrete cell states, population and abstract local signals. Anti-PD1 exposure does not guarantee killing; a chronic diagnosis does not establish exhaustion. Those hypotheses need independent perturbation data.

## What v2 implements

| Input or operation | Implemented behavior |
|---|---|
| Annotated `.h5ad` / in-memory AnnData | Counts, reviewed cell type, explicit lineage, sample IDs and optional measured XY; sparse input stays sparse during selection |
| Xenium output | Exact cell-ID joins across counts, cells and reviewed annotation CSV; explicit bounded region selection |
| Context | Organism, tissue, disease, assay, enrichment/subsampling, missing populations, perturbation assumptions, recognition and timestep |
| Guidebook | Panel-aware signature evidence, per-cell action weights, per-type summaries, explicit contextual modifiers and provenance |
| Cell decisions | Reproducible sampling from guidebook preferences, optionally pooled with Jev preferences; local gates apply first |
| Interactions | Centroid-radius neighborhoods or explicitly assumed well-mixed neighborhoods; activation, exhaustion/recovery, abstract secretion, killing, clearance and division |
| Parallel execution | Concurrent bounded Jev requests; all questions use the same pre-update snapshot; deterministic conflict handling |
| Audit | Durable request intents/results, no automatic retry of unknown attempts, budgets, configuration-bound cache, last-committed-state snapshot, run lock, decision logs and offline replay |
| Comparisons | Same measurement input, scoring rules and simulation settings; latest common completed step; context differences exported |

The spatial kernel is a centroid graph. It does not model mechanics, molecular diffusion, pharmacokinetics or emergent transcriptomes. Living-cell capacity is a computational limit; deaths free capacity. Daughter positions are synthetic offsets. This general kernel is separate from the preserved, lung-specific v0.1.3 physical engine.

## First run: no dataset or API key needed

Python 3.12 or 3.13:

```bash
git clone --branch main https://github.com/hylixinsights/VirtualTissue-Xenium.git
cd VirtualTissue-Xenium
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'

vtx2 prepare examples/v2/synthetic_cells.json \
  examples/v2/context_melanoma_control.json outputs/v2/control.json
vtx2 inspect outputs/v2/control.json
vtx2 run outputs/v2/control.json outputs/v2/control --steps 12 --seed 21
```

Open `outputs/v2/control/replay.html`. This 24-cell example is **synthetic software-test data**.

Run an assumed treatment arm with the same seed:

```bash
vtx2 prepare examples/v2/synthetic_cells.json \
  examples/v2/context_melanoma_anti_pd1.json outputs/v2/anti_pd1.json
vtx2 run outputs/v2/anti_pd1.json outputs/v2/anti_pd1 --steps 12 --seed 21
vtx2 compare outputs/v2/control outputs/v2/anti_pd1
```

The treatment example declares an illustrative activation-weight modifier and leaves recognition off. It is a sensitivity experiment, not an efficacy estimate. Inspect every assumption and repeat across seeds before interpreting differences.

## Use your annotated object

Complete the installation above, then follow the [own-data guide](docs/v2/YOUR_DATA.md) to create a reviewed `annotated_region.h5ad` and a tissue-specific `context.json`. The guide includes a complete context template and the required lineage mapping. Do not reuse the synthetic example's tissue or treatment assumptions for a different experiment.

```bash
pip install -e '.[xenium]'
vtx2 prepare annotated_region.h5ad context.json outputs/v2/sample.json \
  --cell-type-key cell_type --lineage-key lineage --sample-key sample_id --layer counts
vtx2 inspect outputs/v2/sample.json
vtx2 run outputs/v2/sample.json outputs/v2/sample_local --steps 12
```

`obs` must contain reviewed cell types and a mapping to supported lineages. `var_names` must contain unique gene symbols. Select a bounded region/subset first; the default is at most 1,000 starting cells. Coordinates in `obsm['spatial']` must be micrometres, or supply `--coordinate-scale`. Separate samples never share graph edges. An existing canonical `obs['sample_id']` is automatically retained; otherwise supply `--sample-key` for multisample objects.

For scRNA-seq or single-nucleus data without measured XY, set `assay` to `scrna` or `snrna` and `interaction_mode` to `well_mixed`. The replay labels its layout as schematic. Enriched libraries describe the sampled compartment; absent populations are not reconstructed. Mouse or other species require reviewed gene sets; no orthology is inferred.

See the [context/guidebook contract and Python API](docs/v2/CONTEXT_AND_GUIDEBOOK.md), the [workflow](docs/v2/WORKFLOW.md), and `notebooks/V2_Context_Guidebook_Workflow.ipynb`.

Existing testers can update a clean checkout with `git switch main`, `git pull --ff-only origin main`, and `pip install -e '.[xenium]'` in their virtual environment. Users with local changes should preserve those changes before switching branches.

## Public skin and brain examples

The [dataset guide](docs/v2/DATASETS.md) includes verified direct downloads, checksum manifests, annotation tasks and limitations:

- **Human skin primary dermal melanoma, Xenium Prime 5K:** one untreated technical demonstration, 112,551 segmented cells and 5,006 panel genes. Download the required individual files instead of the multi-GB image bundle.
- **Mouse brain WT versus TgCRND8:** an independent disease-context engineering benchmark; paired manifests for the 5.7-month samples. This is not repeated measurement of the same cells or a treatment-response cohort.
- A published systemic-sclerosis skin cohort is documented as a future biological-validation candidate with its access and design constraints.

```bash
vtx2 dataset list
vtx2 dataset download melanoma_5k data/v2/melanoma_5k --include-optional
vtx2 annotation-template data/v2/melanoma_5k annotations.csv \
  --center 2500 2500 --max-cells 1000
# Review the actual ROI and fill cell_type/lineage; vendor clusters are not cell types.
vtx2 prepare-xenium data/v2/melanoma_5k context.json outputs/v2/melanoma.json \
  --annotations annotations.csv --center 2500 2500 --max-cells 1000
```

The example center is a coordinate-selection example, not a biologically chosen niche. Review the tissue map and use your own selection. Annotation templates retain vendor-QC omissions explicitly.

## Run Jev

The API transmits biological context and summarized evidence. Select the cells and metadata you intend to send. Keep API keys in environment variables or notebook secrets.

```bash
# Estimate calls and input tokens first; no API use. A 1,000-cell bundle needs
# about 50 HTTP calls per step, so 12 steps exceed the default budget of 100.
vtx2 cost outputs/v2/sample.json --steps 12 --max-http-attempts 100

# Set TYPESAFE_API_KEY or JEV_API_KEY in your environment, then pass a budget
# that covers the horizon reported by vtx2 cost.
vtx2 run outputs/v2/sample.json outputs/v2/sample_jev --mode jev \
  --steps 12 --seed 21 --max-http-attempts 700 --workers 2
```

A run stops before the first round its budget cannot cover and records `provider_stopped_before_round_commit`; it never falls back to local decisions. `--jev-prior hidden` omits the guidebook prior from the prompt so that anchoring can be measured; the report records the Jensen–Shannon divergence between Jev and the prior per lineage.

The model is pinned to `jev-1.14.0` by default. Use `--jev-weight 0.5` for the explicit normalized pooling rule `prior^(1-alpha) × Jev^alpha`; `1` uses Jev preferences, while local mode uses the guidebook alone. Pooling reinforces preferences and is not a Bayesian evidence update. All choices are sampled with stable per-cell/per-step seeds; Jev's returned argmax is retained for audit.

Inspect the preflight and `vtx2 cost` before spending API credits. A request cap may stop before the requested horizon. Resume with the **identical command plus `--resume`**: successful responses are replayed from t0 without new calls; a failed or unresolved attempt is never resent automatically. Do not delete journals to bypass that safeguard. A stale `.running.lock` requires checking that no process remains before manual removal.

## Tests and evidence

```bash
pip install -e '.[xenium,test]'   # test_v2_datasets needs h5py/pandas
pytest tests/test_v2_*.py -q
```

Null-model controls are one command away and `vtx2 compare` matches them to their source arm:

```bash
vtx2 control outputs/v2/control.json outputs/v2/control_shuffled.json \
  --kind shuffle_expression_within_type --seed 1
vtx2 run outputs/v2/control_shuffled.json outputs/v2/control_shuffled --steps 12 --seed 21
vtx2 compare outputs/v2/control outputs/v2/control_shuffled
vtx2 control outputs/v2/control.json --kind list
```

If a shuffled-expression arm stays inside the seed envelope of the real arm, per-cell expression carries no information in this model for that dataset; report that rather than the label "expression-informed". The 24-cell synthetic fixture has identical expression within each lineage, so it cannot show an expression effect by construction.

Recorded evidence is in [docs/v2](docs/v2/WORKFLOW.md#validation-evidence): verified real downloads; 1,000-cell adapter/simulation checks for melanoma and each brain sample; a small live Jev protocol test; and offline regression tests. Unknown real-data annotations deliberately permit REST only. These checks establish software behavior and input compatibility, not biological prediction.

## v0.1.3 compatibility

The original `vtx` commands, notebooks and hashed engine remain available. Install the legacy environment with `pip install -e '.[legacy,notebooks]'`. Follow the immutable [v0.1.3 README](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/README.md) and [legacy guides](docs/README.md). The v2 JSON input and output formats are distinct; old physical-engine bundles are not silently migrated.

MIT software license. Public datasets retain their own licenses and attribution requirements. Cite the software, data provider/study and the pinned Jev model. Contributions should include assumptions, supporting evidence and validation status.
