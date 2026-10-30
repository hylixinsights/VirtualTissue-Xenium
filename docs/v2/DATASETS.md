# Public data for the v2 workflow

Start with **human skin melanoma for a rich expression panel**, then use an **age-matched mouse brain genotype contrast** to test transfer to a different tissue. Downloads and adapters are implemented; reviewed biological cell annotations and experimental calibration remain investigator tasks. An annotation template deliberately leaves cell types and lineages blank.

## 1. First dataset: public Xenium Prime 5K dermal melanoma

[10x Genomics source and license](https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-skin): one FFPE primary dermal melanoma donor, 112,551 segmented cells, and a development panel with 5,006 Gene Expression targets. The source describes necrosis in the section. This is a public technical example without a paired treatment arm. It supports initialization and hypothesis generation, rather than estimating anti-PD1 efficacy. Data license: CC BY 4.0; credit 10x Genomics and identify processing/annotation changes.

The source flags **ADGRG1, NMBR, OSMR, OSTF1, and RGS8** for nonspecific binding in this development panel. Exclude these from evidence modules or explicitly flag them during review. The loader retains measured values and their feature identities so this decision stays auditable.

| Selected individual file | Bytes | Purpose |
|---|---:|---|
| `cell_feature_matrix.h5` | 44,041,809 | Sparse counts; Gene Expression features selected explicitly |
| `cells.parquet` | 2,289,204 | Measured cell IDs and centroid coordinates |
| `gene_panel.json` | 2,933,339 | Assay design/probe metadata |
| `metrics_summary.csv` | 2,189 | Vendor run QC |
| `experiment.xenium` | 2,323 | Experiment metadata |
| `analysis.tar.gz` (optional) | 30,851,022 | Vendor graph clusters and differential expression for annotation review |

The essential download is **49,268,864 bytes (47.0 MiB)**; with clustering it is **80,119,886 bytes (76.4 MiB)**. Image data are unnecessary for the input adapter. These direct file URLs, byte sizes, and independently measured SHA-256 values are pinned in the [source manifest](../../examples/v2/datasets/melanoma_5k.json). The full output archive is 7,556,661,756 bytes and is intentionally avoided.

```bash
# Run after installing this repository; no API key is needed.
python -m virtualtissue_xenium.v2.datasets list
python -m virtualtissue_xenium.v2.datasets download melanoma_5k data/v2/melanoma_5k --include-optional
python -m virtualtissue_xenium.v2.datasets verify melanoma_5k data/v2/melanoma_5k --include-optional
```

Downloads have a 100 MiB default ceiling, verify every file's size and SHA-256, and reuse matching files. They write `SOURCE.json` and `DOWNLOAD_RECEIPT.json`. Checksums were computed from files retrieved over HTTPS on 2026-10-30; they are reproducibility pins, not a separate vendor attestation. Raw data live under ignored `data/` and are excluded from Git.

## 2. Create reviewed annotations and select a region

1. Open the sample in Xenium Explorer and inspect morphology, panel coverage, quality, graph-cluster markers, and cell types. A graph-cluster number does not establish identity or malignancy. Retain uncertain cells as `lineage=unknown`.
2. Choose a biologically meaningful region, record its center in micrometres, and include neighboring stromal/immune cells. Region selection is an explicit sampling choice; a cropped region does not represent the entire section.
3. Generate an exact-ID annotation template. For example, the following center is a reproducible **software demonstration**, positioned at the first vendor cell; replace it after morphology review:

```python
from virtualtissue_xenium.v2.datasets import write_annotation_template

write_annotation_template(
    "data/v2/melanoma_5k", "outputs/melanoma_annotations.csv",
    max_cells=1000, center=(174.199570, 2299.126221),
)
```

4. Fill every selected `cell_type` and `lineage` in the CSV, using reviewed reference mapping and markers. Vendor graph clusters, if downloaded, appear only in the `cluster_id` column. Some raw cells have no vendor cluster after vendor QC; those entries remain blank and cells are preserved for the investigator's QC decision. An annotation CSV may cover the region or entire measured sample. Unknown IDs, duplicate IDs, and missing selected labels are errors.
5. Load the measured sparse counts and coordinates using the same region selection:

```python
from virtualtissue_xenium.v2.datasets import load_xenium_cells
from virtualtissue_xenium.v2.inputs import prepare, read_json, write_json

cells, genes = load_xenium_cells(
    "data/v2/melanoma_5k", "outputs/melanoma_annotations.csv",
    max_cells=1000, center=(174.199570, 2299.126221),
)
context = read_json("my_reviewed_context.json")
# Declare the actual selected region and sampling limitations in this context.
bundle = prepare(cells, context, assayed_genes=genes)
write_json("outputs/melanoma_prepared.json", bundle)
```

`max_cells` is a computational ceiling. Without an explicit center the loader refuses a larger section; with a center it takes the nearest cells up to that limit, optionally within `radius_um`. Ties follow the vendor matrix order. It never fabricates geometry, cell labels, or excluded populations. Interactions in the portable engine use centroids, not the original polygon mechanics. Keep independent donors together in donor-level train/validation splits; do not split adjacent cells between training and testing.

## 3. Second dataset: age-matched WT and TgCRND8 mouse brain

The [official disease-model dataset](https://www.10xgenomics.com/datasets/xenium-in-situ-analysis-of-alzheimers-disease-mouse-model-brain-coronal-sections-from-one-hemisphere-over-a-time-course-1-standard) contains six FFPE male mouse brain hemispheres: WT and APP-overexpressing TgCRND8 at three ages, with a 248-gene brain base panel plus 99 custom glial/plaque targets (one base target was removed in final validation). Data are CC BY 4.0. Amyloid IF was performed after Xenium. This is a disease-model resource generated on a prototype instrument and reanalyzed with XOA 1.4.0.

Use the **5.7-month WT/Tg pair first**. It avoids the unequal oldest ages in the full series. Each genotype/age has one animal, so the pair is a technical and biological plausibility benchmark with limited replication. The cells are measured in independent animals, not tracked over time; age differences cannot identify cell-level transition rates. Mouse gene symbols differ from human symbols and the compact human demonstration modules must be replaced with reviewed mouse modules.

| ID | Count matrix | Cell table | Metrics | Total bytes |
|---|---:|---:|---:|---:|
| `brain_tg_5_7m` | 4,634,237 | 3,055,560 | 850 | 7,690,647 |
| `brain_wt_5_7m` | 5,362,551 | 3,045,252 | 857 | 8,408,660 |

```bash
python -m virtualtissue_xenium.v2.datasets download brain_tg_5_7m data/v2/brain_tg_5_7m
python -m virtualtissue_xenium.v2.datasets download brain_wt_5_7m data/v2/brain_wt_5_7m
```

Use the same annotation/template adapter. Review cell types against a mouse brain reference, map them explicitly to supported lineages, keep neurons nonproliferative, and inspect panel coverage per action before running. Compare initial measured states and spatial niches first. Later simulated differences are model-dependent hypotheses; they are not observed longitudinal disease progression.

## 4. A published human chronic-skin-disease cohort for later validation

[Multimodal analyses of early, untreated systemic sclerosis skin identify a proinflammatory vascular niche of macrophage-fibroblast signaling](https://doi.org/10.1172/jci.insight.198954) provides Xenium Prime 5K skin data in [GEO GSE312932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312932). For example, [GSM9357164, SSc10](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9357164) lists a 2.5 MB H5 count matrix, cell tables, and boundaries independently of the large images. This is a better candidate for donor-level chronic-disease validation after the technical workflow works. The paper also analyzes nintedanib in 3D skin-like tissues; do not label the human Xenium biopsies as paired nintedanib treatment data.

The accession and availability were verified, but this cohort is **not yet packaged or checksum-pinned here**. Obtain sample metadata, reviewed cell annotations, donor/condition mapping, and reuse terms from the original source before integrating. A GEO download being public does not automatically assign the 10x public-demo CC BY license to independent investigators' data.

An additional attractive resource is [Hidden immune memory niches in inflammatory skin diseases](https://pmc.ncbi.nlm.nih.gov/articles/PMC13041989/), with Xenium skin and post-treatment analyses. At review time its [author repository](https://github.com/haniffalab/spatial_skin_atlas) still labeled the data portal “TBC”; it is therefore a follow-up candidate, not a claimed ready-to-download example.

## 5. What was actually verified in this revision

All essential melanoma files, its clustering archive, and the 5.7-month WT/Tg brain matrices and cell tables were downloaded and checksum-verified. The adapter loaded 1,000 measured cells from each dataset, excluded control feature types, and kept exact ID/coordinate/count alignment. For these **engineering smoke tests only**, cells were explicitly labeled unreviewed with lineage `unknown`; no biological annotations were invented and no biological response validation is claimed. Results and file hashes are in [DATA_DOWNLOAD_VALIDATION.json](DATA_DOWNLOAD_VALIDATION.json).

Biological acceptance requires reviewed labels/modules, donor-separated benchmarking, predefined outcome metrics, comparison against expression-only and fixed-policy baselines, held-out perturbation data, and uncertainty assessment. Neither these initial snapshots nor a Jev answer alone establishes a calibrated transition probability.

## 6. Collaborator annotation-review task table

Complete these review tasks before interpreting a biological simulation. The downloaded vendor graph clusters only help organize the review.

| Order | Task | Suggested owner | Deliverable and acceptance check | Dependency |
|---|---|---|---|---|
| 1 | Verify sample, assay, panel, disease, treatment, and donor metadata | Data owner | Source accession, sample/donor IDs, and a documented measured-versus-assumed context; no treatment label inferred from a paper title | Download receipt |
| 2 | Inspect transcript depth, segmentation, controls, tissue morphology, and excluded probes | Spatial analyst + pathologist | QC report and explicit exclusion rules; do not assign biological meaning to an empty expression profile | 1 |
| 3 | Review a species- and tissue-matched reference and assign broad cell identities | Cell annotation analyst | Reviewed annotation table/object with markers and reference provenance; uncertain identities remain `unknown` | 1–2 |
| 4 | Review malignant identity and fine immune/glial states | Tissue-domain specialist | Evidence supporting each refined identity; morphology/CNV/mutation or orthogonal evidence where appropriate; avoid malignancy or functional exhaustion claims based on one transcript | 3 |
| 5 | Map every reviewed cell type to a supported action lineage | Tissue-domain specialist + simulation analyst | Explicit `cell_type`→`lineage` mapping with one lineage per label; mature neurons cannot gain proliferation through a generic context rule | 3–4 |
| 6 | Review signature coverage and species symbols | Molecular analyst | Versioned gene sets and panel overlap table; no implicit human-to-mouse orthology conversion; unavailable modules remain missing | 2–5 |
| 7 | Select a pathology-reviewed region and declare sampling limitations | Spatial analyst + pathologist | Recorded center/radius or annotated-object subset, preserved neighbors, and explicit missing-population/enrichment assumptions | 2–5 |
| 8 | Review action priors, context modifiers, recognition, controls, and outcome metrics | Experimental biologist + simulation analyst | Signed-off configuration and matched controls; expected outcomes stay hypotheses and are not forced into the prompt | 5–7 |
| 9 | Prepare, inspect preflight, run a small local preview, then a bounded Jev experiment | Simulation analyst | Prepared fingerprint, guidebook, preflight, journals, and complete stop/status report; zero unresolved input errors | 8 |
| 10 | Evaluate against held-out donor/perturbation measurements | Independent collaborator | Predefined metrics, baseline comparisons, sensitivity analysis, uncertainty, and a written claim limited to the available evidence | 9 |

QC filtering can be performed upstream in a reviewed annotated AnnData object and passed through `from_anndata`; retain its filter report. The direct raw-Xenium adapter preserves selected measured cells and requires exact matrix/cell-table ID agreement. It does not perform morphology-based QC or quietly remove cells missing vendor clusters.

The [real-data smoke-test report](REAL_DATA_SMOKE_TEST.json) records complete guidebook/preflight/local-run results for reproducible engineering regions. Its explicit `unknown` annotations are a compatibility test only and must be replaced before tasks 3–10 can support biological conclusions.
