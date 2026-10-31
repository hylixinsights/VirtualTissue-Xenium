# Run v2 with your annotated data

Use this guide when you already have reviewed cell annotations in an AnnData `.h5ad` object. The order is **install → check annotations/counts → select a bounded region → declare context → prepare/review the guidebook → inspect → local run → Jev run → compare matched arms**. The `vtx2` workflow is the general v2 research preview; `vtx` and the older notebooks use the legacy lung-specific workflow.

For raw Xenium output without an annotated object, follow [Public datasets and the raw-Xenium adapter](DATASETS.md): it explains the exact-ID annotation CSV, `annotation-template`, region selection and `prepare-xenium` path. Vendor cluster numbers alone are not biological annotations.

## 1. Install the current main branch

Use Python **3.12 or 3.13**. From a terminal:

```bash
git clone --branch main https://github.com/hylixinsights/VirtualTissue-Xenium.git
cd VirtualTissue-Xenium
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[xenium]'
vtx2 --help
```

If your Python executable is named `python3.13`, use that when creating the environment. The `xenium` extra includes the AnnData input dependencies. A Jev key is unnecessary for preparation, inspection and local mode. Use the main branch version displayed in `pyproject.toml` when reporting a test, and record `git rev-parse HEAD` so collaborators can reproduce the exact code.

## 2. Check the annotated object

| Location | Required content |
|---|---|
| `adata.obs_names` | Unique cell IDs; preserve their alignment with expression and coordinates. Prefix by sample when merging objects. |
| `adata.obs['cell_type']` | Reviewed biological labels, with no missing entries. |
| `adata.obs['lineage']` | Explicit supported action lineage for each cell, with no missing entries. |
| `adata.obs['sample_id']` | Donor/section identifier separating independent tissues. Use distinct IDs for sections whose coordinates have independent origins. |
| `adata.var_names` | Unique gene symbols from the full assayed panel, including features with zero detected counts in the selected region. Ensembl IDs need a reviewed upstream symbol conversion. |
| `adata.layers['counts']` | Finite nonnegative, unscaled expression counts. If counts are in `adata.X`, omit `--layer counts` below. Do not pass centered/scaled residuals. |
| `adata.obsm['spatial']` | For spatial mode, measured XY coordinates with one row per cell. Default units are micrometres. |

The adapter automatically uses `obs['sample_id']` when present; the commands below name it explicitly. If your column has another name, use `--sample-key your_column`. Without either, all cells receive one sample ID, so ensure that input is a single tissue section. Nearby cells from independent sections must not become neighbors.

Choose lineage values from this exact list:

```text
T_CD8, T_CD4, B, myeloid, fibroblast, endothelial, epithelial, malignant,
neuron, astrocyte, oligodendrocyte, microglia, unknown
```

Map your reviewed labels to these roles explicitly. For example, a reviewer may approve `"CD8 T cell" → "T_CD8"` and `"fibroblast" → "fibroblast"`. A tumor-derived epithelial cell is not automatically malignant. Use `malignant` only with supporting annotation evidence. An unresolved population should have `lineage='unknown'`, which permits only `REST`. All cells carrying the same `cell_type` must map to the same lineage; refine ambiguous labels first. Mature neurons and oligodendrocytes cannot proliferate in the current action contract. See the full [action table and guidebook semantics](CONTEXT_AND_GUIDEBOOK.md).

For Seurat objects, first export a reviewed `.h5ad` using your established Seurat/AnnData conversion tool. V2 does not read `.rds` directly or supply a Seurat exporter. After conversion, check counts, gene symbols, cell-ID alignment, annotations, sample IDs and spatial units in the resulting AnnData; do not assume the converter placed counts in `X` or preserved your coordinate field.

## 3. Select a small, meaningful input

Start with a few hundred cells and at most **1,000 cells** for the first test. The adapter rejects larger inputs by default. Select a pathology-reviewed region or a predefined cell subset upstream, retaining relevant neighboring populations and all gene features. Randomly discarding cells changes the interaction network and must be declared as sampling.

For example, after defining a `reviewed_roi` boolean column upstream:

```python
import anndata as ad

adata = ad.read_h5ad("my_annotated_tissue.h5ad")
# Replace this sample value and ROI column with your reviewed selection.
mask = (adata.obs["sample_id"] == "donor1_section1") & adata.obs["reviewed_roi"].eq(True)
region = adata[mask].copy()
assert 0 < region.n_obs <= 1000, "Select a smaller reviewed region first."
region.write_h5ad("annotated_region.h5ad")
```

Do not restrict `var_names` to highly variable genes for this workflow. Panel overlap determines whether a signature is available; a filtered feature list cannot distinguish an unassayed gene from a discarded measurement. Keep a record of segmentation/count QC and the cells removed upstream.

## 4. Save a complete context JSON

Save the following as `context.json` and replace the biological declarations with your experiment. This is a human melanoma **control example**, not a diagnosis inferred by the software:

```json
{
  "schema_version": "2.0",
  "tissue": "skin melanoma",
  "organism": "human",
  "disease": "melanoma",
  "assay": "xenium",
  "sampling": {
    "design": "subsampled",
    "description": "Pathology-reviewed region from donor1_section1; replace with the actual selection and enrichment procedure.",
    "missing_populations": ["cells outside the selected region", "cells removed by upstream QC"]
  },
  "perturbations": [],
  "recognition": false,
  "time_step_hours": 1,
  "interaction_mode": "spatial",
  "neighbor_radius_um": 30,
  "max_neighbors": 20,
  "action_priors": {},
  "context_modifiers": []
}
```

Use `sampling.design: "enriched"` for an enriched library, or `"whole_tissue"` only when that describes your sampling design. Declare missing/depleted populations explicitly; enrichment is not corrected into whole-tissue cell frequencies. Default signatures support human gene symbols. For mouse or another species, supply reviewed gene sets with `--gene-sets reviewed_species_signatures.json`; capitalization does not establish orthology.

For dissociated single-cell data, use `assay: "scrna"`; for single nuclei, use `"snrna"`. When measured spatial positions are unavailable, set `interaction_mode: "well_mixed"` explicitly and describe sorting/nuclei selection in `sampling`. This samples hypothetical neighbors within each sample and does not reconstruct tissue anatomy. Do not put UMAP coordinates in `obsm['spatial']` to imitate physical distance. Nuclear expression does not establish whole-cell protein activity.

For a hypothetical intervention, add a perturbation with `type`, `name`, `status: "hypothetical"`, `targets` and `mechanism_assumptions`. Use `status: "observed"` for treatment actually applied to the source sample; distinguish this history from a proposed additional intervention. Optional `dose` and `duration` are descriptive text, and `evidence` is a list of references. These fields do not implement drug kinetics.

The [anti-PD1 example](../../examples/v2/context_melanoma_anti_pd1.json) includes an explicitly assumed activation-weight modifier. A drug name can also be supplied as prompt-only context without a modifier. Neither anti-PD1 nor a chronic-disease label guarantees effector activity or exhaustion. `recognition: true` is a separate investigator assumption: killing additionally needs an activated, nonexhausted CD8 and a live malignant graph neighbor. Initial activated/exhausted states default to false; a reviewed Python transformation can supply explicit `initial_state` flags as described in the [contract](CONTEXT_AND_GUIDEBOOK.md#6-native-anndata-and-initial-states).

## 5. Prepare, review and inspect

Run from the repository root, with `annotated_region.h5ad` and `context.json` at the paths shown:

```bash
vtx2 prepare annotated_region.h5ad context.json outputs/my_tissue/prepared.json \
  --cell-type-key cell_type --lineage-key lineage --sample-key sample_id \
  --layer counts --max-cells 1000
vtx2 inspect outputs/my_tissue/prepared.json
```

If measured coordinates are in another unit, pass `--coordinate-scale` to convert to micrometres; for example, `--coordinate-scale 1000` converts millimetres. The CLI expects `obsm['spatial']`; move a different measured XY array to that key upstream, retaining cell order. Omit `--layer counts` only when `X` contains the intended counts.

Preparation writes `prepared.json` and the readable sidecar `prepared.guidebook.json`. Review the sidecar before an API run:

- `context`: the normalized experiment declaration and explicit assumptions.
- `by_cell_type`: cell counts, compatible menus, signature means and unavailable evidence.
- `cells`: per-cell baseline evidence, panel coverage, action weights and applied modifiers.
- `provenance`: full panel, gene sets, scoring method, probability semantics and input/context hashes.

`inspect` validates the bundle and reports the first-round question/request estimate, lineage composition and isolated cells. Many unknown labels or isolated cells may indicate unresolved annotation or an unsuitable region/radius. Zero overlap stays missing, not biologically inactive. These descriptive ranks are not AUCell or learned transition probabilities. If you change assumptions, edit the source context or reviewed input and prepare a new bundle at a new path; do not edit the fingerprinted bundle or sidecar as if it updated the simulation input.

## 6. Run a free local check, then Jev

```bash
vtx2 run outputs/my_tissue/prepared.json outputs/my_tissue/local \
  --mode local --steps 3 --seed 21 --max-population 2000
```

Open `outputs/my_tissue/local/replay.html` in a browser. Check `report.json` for `completed`, `completed_steps`, `stop_reason`, population counts and events. A local run samples the guidebook's assumed preferences and uses no API. If all lineages are unknown, the expected result is REST-only behavior; this tests compatibility, not biological dynamics.

Set your own credential in the terminal environment, then run a small bounded Jev experiment:

```bash
export TYPESAFE_API_KEY="replace-with-your-own-key"
vtx2 run outputs/my_tissue/prepared.json outputs/my_tissue/jev \
  --mode jev --steps 3 --seed 21 --max-population 2000 \
  --model jev-1.14.0 --jev-weight 0.5 \
  --max-http-attempts 100 --max-questions 20 --workers 2
unset TYPESAFE_API_KEY
```

Use your normal credential manager/environment setup if available; keep real credentials out of notebooks, source files and shared command transcripts. `JEV_API_KEY` is also supported. The API receives context and cell-level summarized evidence; use data you are authorized to send to that service.

The request limit covers HTTP attempts, not cells or completed simulated hours. Inspection estimates the first round; later rounds depend on living cells and their menus. If the remaining budget cannot cover a full round, the run stops before committing that round and records the reason. It does not silently switch to local decisions. Read the report even if a replay was produced. Reuse a directory only with `--resume` and identical settings; successful cached responses are reused, while failed or unresolved attempts are not automatically resent. Preserve journals when diagnosing an interrupted run.

## 7. Inspect decisions and compare matched contexts

Each run directory contains:

| File | Inspect for |
|---|---|
| `report.json` | Completion/stop reason, completed horizon, populations, events and usage. |
| `config.json` | Exact settings, input fingerprint, model and explicit kernel assumptions. |
| `guidebook.json` | The actual guidebook used by the run. |
| `decisions.jsonl` | Allowed/excluded actions, observations, guidebook priors, Jev preferences and sampled action. |
| `frames.json`, `checkpoint.json` | Saved simulated states and last fully committed round. |
| `replay.html` | Offline visualization for inspecting the saved run. |
| `jev_journal/` (Jev mode) | Durable request intents/results and usage for audit/recovery. |

To test a new context, prepare the **same measured object** with a separately reviewed context file and a new output path, then run it with the same mode, seed, step count, model, batch setting, graph and population cap. For example, after saving `context_treatment.json`:

```bash
vtx2 prepare annotated_region.h5ad context_treatment.json outputs/my_tissue/treatment.json \
  --cell-type-key cell_type --lineage-key lineage --sample-key sample_id \
  --layer counts --max-cells 1000
vtx2 inspect outputs/my_tissue/treatment.json
vtx2 run outputs/my_tissue/treatment.json outputs/my_tissue/treatment_local \
  --mode local --steps 3 --seed 21 --max-population 2000
vtx2 compare outputs/my_tissue/local outputs/my_tissue/treatment_local
```

Local prompt-only arms have identical numerical preferences unless explicit rules change. To compare Jev context effects, run both arms in Jev mode with identical settings and compare those directories. `compare` requires matched measurements, panel/scoring, tissue/disease, sampling and key runtime settings; it compares the latest common completed step. It allows reviewed perturbation/recognition/weight changes and returns each context so you can see what changed. A biological disease-versus-healthy dataset contrast has different measurements and should not be presented as this matched perturbation experiment.

Repeat planned arms across seeds and donors, retaining prompt-only, explicit-modifier and local baselines. The current engine uses fixed baseline expression, centroid/well-mixed graphs, abstract signal feedback and discrete actions. It does not predict a new transcriptome, calibrated time-dependent fate probabilities, pharmacology or physical tissue mechanics. See [Validation plan](VALIDATION_PLAN.md) and [ordered development tasks](WORKFLOW.md) before interpreting simulated differences as treatment or disease predictions.
