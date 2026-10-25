# Your own Xenium data

To simulate your own sample, run **Part 1** on your Xenium output. You need two extra things that the public example ships with: a **data folder with a manifest**, and a **cluster annotation map** that says which cell type each of *your* 10x clusters is. Parts 2–4 do not change.

**Check suitability first.** This tutorial covers **human tumour tissue**, using human gene sets and the cell types listed in §5. Automatic proposals need tumour-candidate epithelium, fibroblasts and immune cells in the same region. The current notebook and CLI cannot reach manual selection if no automatic proposal succeeds. A healthy tissue, a different species, or a tumour section without the required mixture needs work beyond this tutorial (§6).

You do not need a Jev key to prepare your data. You do need someone who can interpret the clusters using markers, tissue context and any available pathology: the template does **not** classify cells automatically, and passing validation does not establish biological correctness.

```text
your Xenium "outs" folder ──(1) create_data_folder──► DATA_DIR (raw/ + manifests)
                          ──(2) write_annotation_template──► MySample_clusters.json (you fill it in)
                          ──(3) Part 1 with CLUSTER_ANNOTATIONS_JSON──► bundle ZIP ──► Parts 2–4
```

---

## 1. What your data must look like

| Requirement | Why |
|---|---|
| Output of **Xenium Onboard Analysis ≥ 2.0** (validated on 3.0.0, Prime 5K) | Part 1 reads the Parquet tables introduced by XOA 2 |
| `cell_feature_matrix.h5` with features of type `Gene Expression` | Only these are scored; controls and codewords are excluded |
| `cells.parquet` with `cell_id`, `x_centroid`, `y_centroid`, `transcript_counts`, `cell_area`, `nucleus_area` | XY in µm; areas set compression floors |
| `transcript_counts` equal to the sum of Gene Expression counts per cell | Checked by Part 1 to catch mismatched files (true for unmodified XOA output) |
| `cell_boundaries.parquet` and `nucleus_boundaries.parquet` (`cell_id`, `vertex_x`, `vertex_y`), one contiguous polygon per cell | Geometry and contacts. Cells without a boundary are refused |
| `analysis.tar.gz` containing `analysis/clustering/gene_expression_graphclust/clusters.csv` and `analysis/diffexp/gene_expression_graphclust/differential_expression.csv` | The annotation map is keyed by these clusters |
| Unique human gene symbols matching the panel | The packaged gene sets are MSigDB `2026.1.Hs`; no species conversion is performed. Duplicates would need an explicit mapping |
| **Tumour tissue** for automatic region proposals | Proposals are centred on tumour-candidate epithelium (§6) |

Use all five required files from **the same sample and XOA run**, with these names:

```text
your_xenium_outs/
├── cell_feature_matrix.h5
├── cells.parquet
├── cell_boundaries.parquet
├── nucleus_boundaries.parquet
└── analysis.tar.gz
```

If you have an extracted `analysis/` folder instead of `analysis.tar.gz`, run `tar -czf analysis.tar.gz analysis` in a terminal inside the outs folder. In a Colab code cell, use `!tar -czf "/full/path/to/outs/analysis.tar.gz" -C "/full/path/to/outs" analysis`, replacing both paths. A matrix alone, a Xenium Explorer subset, or an annotated `.h5ad` is not sufficient input for this preparation workflow.

## 2. Step 1 · Create the data folder

### 2.1 Choose where you will work

**In Colab:** upload the five files above to a folder in **your own Google Drive**, for example `MyDrive/MyXenium/MySample_outs/`. Wait for the uploads to finish. Open `Part1_Prepare_Xenium_Region.ipynb`, save a copy in Drive, and run **only step 1** until it reports *Ready*. If it asks for a restart, restart and rerun step 1. Add a cell immediately after setup with **+ Code**, then run:

```python
from pathlib import Path
from virtualtissue_xenium import colab as vtx_colab, download, own_data

vtx_colab.mount_drive()  # Authorize access to your Drive before using any Drive paths.

# Edit these paths for YOUR sample. Keep outputs outside DATA_DIR.
XENIUM_OUTS = Path("/content/drive/MyDrive/MyXenium/MySample_outs")
DATA_DIR = Path("/content/drive/MyDrive/VirtualTissue/MySample")
OUTPUT_BASE = Path("/content/drive/MyDrive/VirtualTissue/MySample_Prepared")
CLUSTER_ANNOTATIONS_JSON = Path("/content/drive/MyDrive/VirtualTissue/MySample_clusters.json")
SAMPLE_ID = "MySample_Region1"  # Use a different ID and paths for each sample.
```

**In local Jupyter:** first follow the [local installation guide](08_local_installation.md), open Part 1 and run setup. Add a cell immediately after setup with the block below **instead of** the Colab block. Replace `/absolute/path/to/...` with real paths on your computer. No Google account or Drive is needed.

```python
from pathlib import Path
from virtualtissue_xenium import download, own_data

XENIUM_OUTS = Path("/absolute/path/to/MySample_outs")
DATA_DIR = Path("/absolute/path/to/VirtualTissue/MySample")
OUTPUT_BASE = Path("/absolute/path/to/VirtualTissue/MySample_Prepared")
CLUSTER_ANNOTATIONS_JSON = Path("/absolute/path/to/VirtualTissue/MySample_clusters.json")
SAMPLE_ID = "MySample_Region1"
```

Use one path block, then continue through the following cells in the **same notebook session**. If you reconnect in a new runtime, repeat setup and the path block first. Terminal users can use the complete CLI sequence in §4.2.

### 2.2 Create or verify the prepared data folder

Add another code cell and run:

```python
import json

manifest_path = DATA_DIR / "DOWNLOAD_MANIFEST.json"
if manifest_path.is_file():
    # Resume preparation without overwriting existing raw files.
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("sample_id") != SAMPLE_ID:
        raise ValueError("DATA_DIR belongs to another sample. Choose a new DATA_DIR and SAMPLE_ID.")
    check = download.verify_data_folder(DATA_DIR)
    if not check["ok"]:
        raise ValueError("Input folder problems: " + "; ".join(check["problems"]))
else:
    own_data.check_xenium_files(XENIUM_OUTS)
    own_data.create_data_folder(
        XENIUM_OUTS, DATA_DIR, sample_id=SAMPLE_ID,
        include_optional=False,
        source_page=None, license=None, attribution="My lab, unpublished")

print("Data folder ready:", DATA_DIR)
```

This copies the five files into `DATA_DIR/raw/`, checks their structure, and writes `DOWNLOAD_MANIFEST.json` (size and SHA-256 of each file) and `SOURCE.json`. Your original outs folder is preserved. Allow disk space for this additional copy. The second run verifies the existing folder; it does not switch it to a different sample. If copying was interrupted and `raw/` is partly filled, choose a new `DATA_DIR` and rerun this step.

If you already arranged the files in `DATA_DIR/raw/` yourself, register them once with `own_data.write_manifest(DATA_DIR, sample_id=SAMPLE_ID)` instead of calling `create_data_folder`. Optional provenance files such as `analysis_summary.html` can be kept beside the original outs files for review.

> Do not edit files in `raw/` afterwards: Part 1 re-verifies them on every run. If you must replace a file, call `write_manifest` again.

## 3. Step 2 · Write the cluster annotation map

### 3.1 Draft it

```python
CLUSTER_ANNOTATIONS_JSON.parent.mkdir(parents=True, exist_ok=True)
if CLUSTER_ANNOTATIONS_JSON.exists():
    print("Keeping your existing cluster map:", CLUSTER_ANNOTATIONS_JSON)
else:
    own_data.write_annotation_template(DATA_DIR, CLUSTER_ANNOTATIONS_JSON)
    print("Template created:", CLUSTER_ANNOTATIONS_JSON)
```

You get `MySample_clusters.json`, containing a **JSON list with one record per cluster**, and `MySample_clusters.markers.csv`, a table of the fraction of each cluster's cells in which each canonical marker is detected. Rerunning this cell preserves your edits. A record initially looks like this:

```json
{
  "cluster": 11,
  "n_cells": 10263,
  "top_de_genes": ["SLC6A3", "NELL1", "HHLA2", "TMPRSS11E", "ABCC3", "…"],
  "marker_detection_fraction": {"EPCAM": 0.994, "CDH1": 0.978, "MUC1": 0.919, "NKX2-1": 0.976, "CD68": 0.115},
  "label": "Cluster_11_TO_REVIEW",
  "compartment": "unknown",
  "runtime_label": null,
  "markers": [],
  "status": "needs_review",
  "source": "Template from 10x graph-based clustering; fill in after marker review."
}
```

The 10x `analysis_summary.html` in your outs folder (UMAP, cluster plots) is also useful at this point.

### 3.2 Fill in every record

1. Open `MySample_clusters.markers.csv` in a spreadsheet viewer and review it with the top differential-expression genes in the JSON and the original `analysis_summary.html`.
2. **In Colab:** use the Files panel (folder icon), refresh it, and navigate to `drive/MyDrive/VirtualTissue/`. Double-click `MySample_clusters.json` to open the text editor. **Locally:** open the JSON in a plain-text editor such as VS Code or your Jupyter file browser.
3. Edit the fields below in **every record**, retaining the actual `cluster` numbers and the evidence fields generated for your sample. Save the file (`Ctrl+S` / `Cmd+S`) before validation. If editing a downloaded copy, upload it back to exactly the path in `CLUSTER_ANNOTATIONS_JSON`.
4. Keep the outer `[` and `]`, double quotes, commas between records, and JSON `null` (not Python `None`). Do not use a word processor.

| Field | Required | Allowed values | Meaning |
|---|---|---|---|
| `label` | yes | any text, e.g. `Epithelial_tumor_candidate`, `CAF`, `T_CD8_like` | What the cluster is (shown in tables and plots) |
| `compartment` | yes | `tumor_candidate`, `fibroblast`, `immune`, `stroma`, `epithelial_other`, `unknown` | Used for region proposals and colours. **Region proposals need `tumor_candidate`, `fibroblast` and `immune` cells** |
| `runtime_label` | yes (may be `null`) | `Tumor`, `CAF`, `Fibroblast`, `Endothelial`, `Macrophage`, `T_CD8_cytotoxic`, `T_CD4`, `B_cell`, `Plasma_cell`, or `null` | The simulated cell type. `null` = context (occupies space, never decides) |
| `markers` | yes (may be `[]`) | gene symbols | Markers that support your interpretation (documentation) |
| `status` | no | text | e.g. `reviewed_by_pathologist`, `tentative` |
| `source` | no | text | How you decided |

Use `null` for mixed, doublet-like, low-quality or unsupported clusters (e.g. NK cells, neutrophils, mast cells, pericytes, smooth muscle). Being conservative is better than letting a cluster act with the wrong menu. The packaged map of the public sample, [`engine/resources/cluster_annotations.json`](../src/virtualtissue_xenium/engine/resources/cluster_annotations.json), is a worked example: 16 of its 32 clusters are `null`.

**Example of a filled map — format only.** These three cluster numbers and assignments are illustrative. Your file must contain **all of your own clusters**, with assignments supported by your evidence; do not paste this list over the generated template. Epithelial markers alone do not prove malignancy.

```json
[
  {
    "cluster": 1,
    "label": "Epithelial_tumor_candidate",
    "compartment": "tumor_candidate",
    "runtime_label": "Tumor",
    "markers": ["EPCAM", "KRT8", "KRT18"],
    "status": "tentative",
    "source": "Example only: epithelial markers plus independent tissue review; malignancy remains provisional."
  },
  {
    "cluster": 2,
    "label": "Fibroblast_like",
    "compartment": "fibroblast",
    "runtime_label": "Fibroblast",
    "markers": ["COL1A1", "DCN", "PDGFRA"],
    "status": "tentative",
    "source": "Example only: fibroblast markers reviewed against other lineages."
  },
  {
    "cluster": 3,
    "label": "T_CD8_like",
    "compartment": "immune",
    "runtime_label": "T_CD8_cytotoxic",
    "markers": ["CD3E", "CD8A"],
    "status": "tentative",
    "source": "Example only: T-cell and CD8 markers; functional competence is not inferred from RNA."
  }
]
```

For an unresolved cluster, replace its placeholder label with a descriptive label such as `Unresolved`, set `runtime_label` to `null`, and record your uncertainty. Assign its compartment only when supported; otherwise use `unknown`. Do not relabel cells just to make region proposals pass.

### 3.3 Validate it

```python
problems = own_data.validate_annotations(CLUSTER_ANNOTATIONS_JSON, data_dir=DATA_DIR)
if problems:
    raise ValueError("Fix the cluster map before continuing:\n  " + "\n  ".join(problems))
print("Valid: every cluster has an allowed annotation.")
print("Copy these values into Part 1 step 2:")
print("DATA_DIR =", DATA_DIR)
print("OUTPUT_BASE =", OUTPUT_BASE)
print("DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING = False")
print("CLUSTER_ANNOTATIONS_JSON =", CLUSTER_ANNOTATIONS_JSON)
```

The validator checks that every cluster has exactly one record, that no placeholder labels remain, and that compartments and runtime labels are allowed. Fix all reported issues, save the JSON, and rerun validation. It checks file structure and coverage, **not whether your cell identities are biologically correct**.

## 4. Step 3 · Run Part 1 with your map

### 4.1 Continue in Colab or local Jupyter

In `Part1_Prepare_Xenium_Region.ipynb`, **replace the form values in step 2** with the four values printed above. The form has its own assignments, so editing the extra cells alone does not change it. For the Colab paths used here:

```text
DATA_DIR                  = /content/drive/MyDrive/VirtualTissue/MySample
OUTPUT_BASE               = /content/drive/MyDrive/VirtualTissue/MySample_Prepared
DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING = False   (unchecked)
CLUSTER_ANNOTATIONS_JSON  = /content/drive/MyDrive/VirtualTissue/MySample_clusters.json
```

In Colab, enter paths directly in the form fields. In local Jupyter, edit the code assignments and keep paths inside quotes, for example `DATA_DIR = "/absolute/path/to/VirtualTissue/MySample"`; use the Python boolean `False` for the download option. The printed lines above show values to enter, not a Python cell to paste unchanged.

Then run the rest of the notebook as usual. Step 3 prints `Custom cluster map: …/resources_custom`: your map is copied into the run folder and the engine reads it from there (`xenium_prepare.RESOURCES` is pointed to that copy). The packaged files are not changed. The notebook **refuses** to annotate a sample that is not the public example without a map. Otherwise the example's map would silently mislabel your clusters.

The Part 1 guide ([01](01_part1_prepare_region.md)) applies unchanged. Check the marker heat map in step 4. If a label is wrong, fix the file at `CLUSTER_ANNOTATIONS_JSON` (or use `ANNOTATION_OVERRIDES` for a quick test) and rerun from step 3.

In step 8, choose an ID that appears in **your** region gallery; `R03` from the public example may not exist in your sample. After export, copy the exact ZIP path printed by Part 1. Open Parts 2–4, set `BUNDLE_ZIP` to that path, choose a new `RUN_NAME` such as `MySample_R01_local01`, and first try `RUN_MODE = "Local preview (no API)"`. For a paid run, follow the [Jev key guide](07_jev_api_key.md), select `Jev`, choose your limits and start. For local Jupyter, also apply the [local result/backup paths](08_local_installation.md#3-the-notebooks-in-local-jupyter).

### 4.2 Complete command-line route

After [installing the package](08_local_installation.md), run the following from the repository folder. Replace the input path and sample name; these examples use new local folders:

```bash
vtx own-data /absolute/path/to/MySample_outs data/MySample --sample-id MySample_Region1
vtx annotation-template data/MySample MySample_clusters.json
```

**Stop here to edit and save `MySample_clusters.json` as described in §3.2.** Then validate it:

```bash
vtx check-annotations MySample_clusters.json --data-dir data/MySample
```

Continue only after validation says `Valid.` Propose regions with:

```bash
vtx prepare data/MySample outputs/MySample_prepared --annotations MySample_clusters.json --list-only
```

Inspect the maps in the printed output folder. Choose a listed region, for example `R01`, then:

```bash
vtx prepare data/MySample outputs/MySample_prepared --annotations MySample_clusters.json --region R01
BUNDLE=$(ls -t outputs/MySample_prepared/run_*/simulation_input_R01_*.zip | head -1)
vtx check "$BUNDLE"
vtx run "$BUNDLE" outputs/runs --run-name MySample_R01 --mode local --hours 24
vtx explore outputs/runs/MySample_R01_local_preview_results.zip --open
```

To run with Jev, set your own `TYPESAFE_API_KEY` as described in [07](07_jev_api_key.md), then:

```bash
vtx run "$BUNDLE" outputs/runs --run-name MySample_R01_jev01 --mode jev --hours 24 --max-http-attempts 100
```

The CLI prints a first-round estimate and asks for `yes` before sending requests. The cap must cover the first round; you can later increase only `--max-http-attempts` and resume with the same input, name and settings. Keep the exact input ZIP, cluster map and result folder. Recreating a template overwrites its destination, so skip `annotation-template` when continuing work on an edited map.

---

## 5. Runtime labels and what each type can do

The menu is built by the kernel for each cell at each decision. An option is offered only if all its gates pass; failed gates are sent to Jev under `excluded`.

| Runtime label | Options (gates in brackets) |
|---|---|
| `Tumor` | `REMAIN_QUIESCENT` · `ENTER_CYCLE` [cycles enabled, not cycling, oxygen, glucose, undamaged, space] · `UPREGULATE_PDL1` [IFN-γ signal] · `SECRETE_VEGFA` [hypoxia] · `GLYCOLYTIC_SHIFT` / `OXIDATIVE_SHIFT` [metabolism menu, off by default] |
| `T_CD8_cytotoxic` | `WAIT` · `REST` · `HOLD_POSITION` · `SCAN_RANDOM` · `MIGRATE_UP_CXCL9` [gradient] · `APPROACH_TUMOR:<id>` [live tumour within 60 µm, not in contact] · `HOLD_CONTACT:<id>` · `DISENGAGE:<id>` · `DEGRANULATE:<id>` [in contact, **recognition declared**, effector competence declared, effector pulses left, cooldown, not terminally exhausted, glucose] · `PROLIFERATE` [recent stimulation, CD4 help received, space] |
| `T_CD4` | `WAIT` · `HELP_CD8:<id>` [in contact with a stimulated CD8, competence declared] · `MIGRATE_UP_CXCL9` [gradient] |
| `Macrophage` | `WAIT` · `APPROACH_CORPSE:<id>` [corpse within 120 µm] · `EFFEROCYTOSE:<id>` [contacted corpse, **clearance competence declared**, capacity 4] · `SECRETE_CXCL9` [IFN-γ] · `SECRETE_VEGFA_TGFB` [lactate or hypoxia] · `PATROL` |
| `CAF` | `WAIT` · `DEPOSIT_MATRIX` · `SECRETE_TGFB` · `SECRETE_VEGFA` [all at a tumour interface] |
| `Fibroblast` | `WAIT` · `DEPOSIT_MATRIX` · `SECRETE_TGFB` [at a tumour interface] |
| `Endothelial` | `STAY_QUIESCENT` · `SPROUT` [angiogenesis enabled, off by default] |
| `B_cell`, `Plasma_cell` | no menu: present as bodies, never asked |

Cells with fewer than two permitted options are not asked.

## 6. Known limitations

* **Tumour-centred proposals.** `propose_regions` places candidate centres over tumour-candidate epithelium and requires ≥ 15 % tumour candidates, ≥ 8 % fibroblasts and ≥ 15 % immune cells. A non-tumour tissue produces *"No tumour-candidate epithelium"* or *"No region meets the declared balance/QC criteria"*, and a manual centre cannot be used without at least one proposal. Supporting other tissues needs a change in the engine, not just new data.
* **Nine simulated cell types.** Other populations can only be context. Adding a type requires new menus and contract entries in the engine.
* **Panel size.** The default `AUCELL_TOP_RANKS = 100` is ≈ 2 % of a 5,000-gene panel. For a ~300–500-gene panel use a proportionally smaller value (≈ 5 % of the panel, e.g. 15–25), and expect some of the 20 signatures to be unavailable (each needs ≥ 5 panel genes; see `signature_coverage.csv`). A few provenance strings inside the engine name "all 5001 assayed genes". They are descriptive text only, and the actual panel size is recorded in `input_audit.json` and the atlas `signatures` block.
* **Other platforms** (CosMx, MERSCOPE, Visium HD) are not supported by Part 1. Their outputs would have to be converted to the Xenium file formats above, including a graph-based clustering in `analysis.tar.gz`.
* Annotations stay **provisional**. Record how you made them in `status` and `source`. They travel with the bundle in `annotation_review.csv`.
