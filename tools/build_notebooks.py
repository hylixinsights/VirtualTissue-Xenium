"""Generate the three notebooks in ../notebooks from readable source.

    python tools/build_notebooks.py

Edit this file, not the .ipynb files, then rebuild. The analysis cells reuse the
validated cells of the original notebooks (Xenium_Lung5K_Prepare_AUCell.ipynb v1.1
and Xenium_Lung5K_Run_Jev_NoDrive.ipynb v1.1); only the setup, path and
selection cells are new, so that the code can come from GitHub or Google Drive and
inputs/results can live in Google Drive.
"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notebooks"
REPO_URL = "https://github.com/hylixinsights/VirtualTissue-Xenium"
DOCS = REPO_URL + "/blob/main/docs/"


class Notebook:
    def __init__(self, name: str):
        self.name, self.cells = name, []

    def md(self, text: str) -> None:
        self.cells.append(dict(cell_type="markdown", metadata={}, source=textwrap.dedent(text).strip() + "\n"))

    def code(self, text: str, hidden: bool = False) -> None:
        meta = {"cellView": "form"} if hidden else {}
        self.cells.append(dict(cell_type="code", metadata=meta, execution_count=None, outputs=[],
                               source=textwrap.dedent(text).strip() + "\n"))

    def write(self) -> Path:
        stem = Path(self.name).stem.lower().replace("-", "_")
        for i, cell in enumerate(self.cells):
            cell["id"] = f"{stem[:24]}-{i:02d}"
            cell["metadata"]["id"] = cell["id"]
            if cell["cell_type"] == "code":
                compile(cell["source"], cell["id"], "exec")
        nb = dict(nbformat=4, nbformat_minor=5, cells=self.cells, metadata=dict(
            colab=dict(name=self.name, provenance=[]),
            kernelspec=dict(display_name="Python 3", language="python", name="python3"),
            language_info=dict(name="python")))
        OUT.mkdir(exist_ok=True)
        path = OUT / self.name
        path.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
        return path


# ---------------------------------------------------------------------------
# Shared setup cell
# ---------------------------------------------------------------------------
def setup_code(bootstrap_module: str | None, extra: str = "", install_default: bool = True) -> str:
    return '''
    #@title Get the VirtualTissue-Xenium code and install the pinned dependencies
    # ✏️ EDIT only if needed — where the code comes from:
    #   github        clone the public repository (default)
    #   google_drive  copy the VirtualTissue-Xenium folder you placed in Google Drive
    #   installed     the package is already installed with pip (local Jupyter)
    CODE_SOURCE = "github"  #@param ["github", "google_drive", "installed"]
    GITHUB_URL = "''' + REPO_URL + '''.git"  #@param {type:"string"}
    GITHUB_BRANCH = "main"  #@param {type:"string"}
    DRIVE_CODE_DIR = "/content/drive/MyDrive/VirtualTissue-Xenium"  #@param {type:"string"}
    ''' + (('INSTALL_DEPENDENCIES = ' + str(install_default) + '  #@param {type:"boolean"}') if bootstrap_module else '') + '''

    import importlib, os, shutil, subprocess, sys
    from pathlib import Path
    print("Running Python:", sys.version.split()[0], "| executable:", sys.executable)
    IN_COLAB = "google.colab" in sys.modules or "COLAB_RELEASE_TAG" in os.environ
    # Temporary Colab disk. Keep the code here, not in Drive.
    CODE_DIR = Path("/content/VirtualTissue-Xenium")
    if not IN_COLAB and CODE_SOURCE != "installed":
        # Local Jupyter: use the repository that contains this notebook.
        here = Path.cwd().resolve()
        CODE_DIR = here.parent if (here.parent / "src" / "virtualtissue_xenium").is_dir() else here
        CODE_SOURCE = "local_repository"
    if CODE_SOURCE == "github":
        if not (CODE_DIR / "src" / "virtualtissue_xenium").is_dir():
            subprocess.check_call(["git", "clone", "--quiet", "--depth", "1", "--branch", GITHUB_BRANCH,
                                   GITHUB_URL, str(CODE_DIR)])
    elif CODE_SOURCE == "google_drive":
        from google.colab import drive
        drive.mount("/content/drive")
        if not (Path(DRIVE_CODE_DIR) / "src" / "virtualtissue_xenium").is_dir():
            raise FileNotFoundError("DRIVE_CODE_DIR must be the VirtualTissue-Xenium folder "
                                    "(the one that contains src/, notebooks/ and README.md).")
        shutil.rmtree(CODE_DIR, ignore_errors=True)
        shutil.copytree(DRIVE_CODE_DIR, CODE_DIR,
                        ignore=shutil.ignore_patterns(".git", "__pycache__", "*.pyc", ".ipynb_checkpoints"))
    if CODE_SOURCE != "installed" and str(CODE_DIR / "src") not in sys.path:
        sys.path.insert(0, str(CODE_DIR / "src"))
    import virtualtissue_xenium as vtx
    from virtualtissue_xenium import colab as vtx_colab
    print("VirtualTissue-Xenium", vtx.__version__, "| engine:", vtx.ENGINE_DIR)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/vtx_matplotlib")
    ''' + extra + '''
    ''' + ('''
    import ''' + bootstrap_module + '''
    importlib.reload(''' + bootstrap_module + ''')
    SETUP = ''' + bootstrap_module + '''.setup_environment(vtx.ENGINE_DIR, install_packages=INSTALL_DEPENDENCIES)
    SETUP_READY = SETUP["ready"]
    ''' if bootstrap_module else '''
    # This notebook only reads saved results: pandas and matplotlib from the runtime are enough.
    SETUP_READY = True
    ''')


def build_part1() -> Path:
    nb = Notebook("Part1_Prepare_Xenium_Region.ipynb")
    nb.md('''
    # Part 1 · From a 10x Xenium sample to a simulation-ready region

    **VirtualTissue-Xenium** turns a Xenium section into a tissue of AI agents: every cell starts from its
    measured state (segmentation, XY position, gene expression → 20 AUCell scores) and, in Parts 2–4, asks
    **Jev** which permitted action to take next. This notebook is **Part 1**:

    | Part | Notebook | What happens |
    |---|---|---|
    | **1** | **this one** | Xenium output → QC → provisional cell types → AUCell → candidate regions → **bundle ZIP** |
    | 2 | `Part2-4_Run_Jev_Simulation.ipynb` | load and verify the bundle ZIP (yours or from Part 1) |
    | 3 | same notebook | choose scenario, duration and API limits; run Jev |
    | 4 | same notebook + `Part4_Explore_Saved_Results.ipynb` | interactive replay, figures, tables |

    **Example data:** *Xenium Prime Human Lung Cancer FFPE, experiment 2* from
    [10x Genomics](https://www.10xgenomics.com/datasets/xenium-human-lung-cancer-post-xenium-technote)
    (CC BY 4.0): 278,328 cells, 5,001 genes (Xenium Prime 5K panel). Only 10 files (~503 MiB) of the
    22.6 GB archive are needed; step 2 can download them for you.

    ### Before you start (Colab)
    1. **Runtime → Change runtime type → Python 3, CPU.** No GPU. The standard runtime is enough
       (peak memory ≈ 3–4 GB; AUCell on the whole section takes a few minutes).
    2. Run the cells **in order**. Every value you may change is a form field marked **✏️**.
    3. The four ✏️ places in this notebook: **code source** (step 1, only if not using GitHub),
       **DATA_DIR / OUTPUT_BASE** (step 2), optional **region-size parameters** (step 3), and
       **REGION_TO_EXPORT** (step 8 — this is where you choose `R03`, for example).

    **No API key is needed for Part 1.** Save a copy of this notebook in your own Drive to keep your settings.

    **Using your own sample?** This workflow targets human tumour tissue with tumour-candidate epithelium,
    fibroblasts and immune cells. It requires the five original Xenium files and a reviewed cluster map;
    it does not classify cells automatically. Run step 1, then **pause before step 2** and follow the
    [own-data guide](''' + DOCS + '''06_your_own_xenium_data.md). It mounts Drive, defines your paths, creates
    the data folder, and shows how to edit and validate the cluster map. Copy the four values it prints into
    step 2 below. If no automatic proposal meets the tumour-interface criteria, this notebook cannot proceed
    to manual region selection.

    **Local Jupyter:** input/output form paths do not change automatically. Before step 2, replace `DATA_DIR`
    and `OUTPUT_BASE` with local paths (for example `../data/Xenium_Public_Lung5K` and `../outputs/prepared`
    when this notebook runs from the repository's `notebooks/` folder). For your own sample, also set the
    local cluster-map path and uncheck the public download.

    Full guide: [docs/01_part1_prepare_region.md](''' + DOCS + '''01_part1_prepare_region.md) ·
    Colab paths and secrets: [docs/05_google_colab.md](''' + DOCS + '''05_google_colab.md)

    **Interpretation:** cluster labels are provisional, marker-based interpretations of the 10x clusters, not a
    validated classifier. “Tumour-candidate epithelium” does not establish malignancy. Antigen recognition,
    protein competence and kinetic rates are not measured by this workflow.
    ''')
    nb.md('''
    ## 1 · Get the code and install dependencies
    Installs the exact package versions the engine was validated with (Python 3.12 or 3.13).
    If this cell prints **ONE SESSION RESTART IS REQUIRED**, choose **Runtime → Restart session** and run this
    cell again (your Drive files are not affected). Continue only when it reports *Ready*.

    ✏️ Leave `CODE_SOURCE = github` once the repository is public. Before that, upload the whole
    `VirtualTissue-Xenium` folder to Google Drive and choose `google_drive`.
    ''')
    nb.code(setup_code("colab_setup_prepare", extra='''os.environ["NUMBA_NUM_THREADS"] = "2"'''
                       ) + '''
    if SETUP_READY:
        prep = vtx.prepare
        from virtualtissue_xenium import download as vtx_download, own_data as vtx_own
        print("Ready. Preparation engine:", prep.VERSION)
    else:
        prep = None
    ''', hidden=True)
    nb.md('''
    ## 2 · ✏️ Where the data are and where results go
    **DATA_DIR** is the folder that contains `raw/`, `SOURCE.json` and `DOWNLOAD_MANIFEST.json`
    (not the `raw/` folder itself). **OUTPUT_BASE** must be outside DATA_DIR; each run creates a new
    `run_<timestamp>` subfolder there. Paths starting with `/content/drive/MyDrive/` are in your Google Drive
    (Drive is mounted automatically and Colab asks for permission once).

    * **Public example:** keep `DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING` checked. If DATA_DIR is empty, the 10 original
      files (~503 MiB) are downloaded from 10x into it and verified by SHA-256 (a few minutes, once).
    * **Your own Xenium run:** complete the
      [own-data guide](''' + DOCS + '''06_your_own_xenium_data.md) first, including mounting Drive before
      accessing its files and validating your edited cluster map. In this form, replace `DATA_DIR`,
      `OUTPUT_BASE` and `CLUSTER_ANNOTATIONS_JSON` with your paths, and **uncheck
      `DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING`**. The form's assignments replace variables from earlier cells,
      so make these changes here too. Keep the public example and each new sample in separate data folders.
    ''')
    nb.code('''
    #@title Data and output folders
    if not globals().get("SETUP_READY", False):
        raise RuntimeError("Complete step 1 first. If a restart was requested, select Runtime > Restart session and rerun step 1.")
    # ✏️ EDIT: folder that contains raw/, SOURCE.json and DOWNLOAD_MANIFEST.json
    DATA_DIR = "/content/drive/MyDrive/VirtualTissue/Xenium_Public_Lung5K"  #@param {type:"string"}
    # ✏️ EDIT: where Part 1 writes its results (must be outside DATA_DIR)
    OUTPUT_BASE = "/content/drive/MyDrive/VirtualTissue/Xenium_Lung5K_Prepared"  #@param {type:"string"}
    DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING = True  #@param {type:"boolean"}
    # ✏️ EDIT only for your own sample: your cluster annotation map (JSON). Empty = public example map.
    CLUSTER_ANNOTATIONS_JSON = ""  #@param {type:"string"}

    if (DATA_DIR + OUTPUT_BASE + CLUSTER_ANNOTATIONS_JSON).count("/content/drive/") and IN_COLAB:
        vtx_colab.mount_drive()
    if not (Path(DATA_DIR) / "DOWNLOAD_MANIFEST.json").is_file():
        if DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING:
            print("DATA_DIR has no DOWNLOAD_MANIFEST.json: downloading the public example from 10x Genomics...")
            vtx_download.download_public_lung5k(DATA_DIR)
        else:
            raise FileNotFoundError("DOWNLOAD_MANIFEST.json was not found. Set DATA_DIR to the folder containing raw/ and "
                                    "DOWNLOAD_MANIFEST.json, not to the raw/ subfolder.")
    check = vtx_download.verify_data_folder(DATA_DIR, check_hashes=False)
    if not check["ok"]:
        raise ValueError("Input folder problems: " + "; ".join(check["problems"]))
    vtx_own.require_annotation_map(DATA_DIR, CLUSTER_ANNOTATIONS_JSON)
    print(f"Input: {DATA_DIR} ({check['files']} files present; SHA-256 is checked in step 3)")
    ''')
    nb.md('''
    ## 3 · Region size and AUCell parameters, then load the sample
    The defaults reproduce the example: regions of at most **1,000 cells including context**, a target of 750
    central cells plus a 30 µm context halo, 8 proposals, AUCell top-100 ranks. The program reduces the central
    area when needed to keep the total ≤ 1,000. All cells inside the selected area are kept, including
    unclassified and low-RNA cells.

    This cell checks every input file against its SHA-256 (set `VERIFY_INPUT_HASHES` off to save ~1 min).
    ''')
    nb.code('''
    #@title ✏️ Parameters, integrity check, expression, XY and clusters
    from datetime import datetime, timezone
    from IPython.display import display, Markdown, clear_output
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    MAX_TOTAL_CELLS = 1000 #@param {type:"integer"}
    CORE_TARGET = 750 #@param {type:"integer"}
    HALO_UM = 30 #@param {type:"number"}
    N_REGIONS = 8 #@param {type:"integer"}
    AUCELL_TOP_RANKS = 100 #@param {type:"integer"}
    VERIFY_INPUT_HASHES = True #@param {type:"boolean"}
    if not 100 <= CORE_TARGET <= MAX_TOTAL_CELLS <= 1000:
        raise ValueError('Use 100 <= CORE_TARGET <= MAX_TOTAL_CELLS <= 1000.')
    if not 1 <= N_REGIONS <= 12: raise ValueError('Choose between 1 and 12 region options.')
    RUN_DIR=Path(OUTPUT_BASE)/datetime.now(timezone.utc).strftime('run_%Y%m%d_%H%M%S_%f')
    if CLUSTER_ANNOTATIONS_JSON:
        RUN_DIR.mkdir(parents=True, exist_ok=True)
        print('Custom cluster map:', vtx_own.use_cluster_annotations(prep, CLUSTER_ANNOTATIONS_JSON, RUN_DIR))
    else:
        vtx_own.use_packaged_annotations(prep)
    print('Input:',DATA_DIR)
    print('Results:',RUN_DIR)
    sample=prep.load_sample(DATA_DIR,RUN_DIR,verify_hashes=VERIFY_INPUT_HASHES,min_counts=20,min_genes=10)
    prep.save_json(RUN_DIR/'environment.json',prep.environment())
    display(pd.DataFrame({'value':[len(sample['obs']),len(sample['genes']),int(sample['obs'].qc_pass.sum())]},
                         index=['Original cells','Gene Expression features','Passed initial RNA QC']))
    ''')
    nb.md('''
    ## 4 · Review provisional cell annotations
    The table lists original differential-expression results and markers supporting each interpretation. The plot
    shows the fraction of cells with each marker detected. A sparse panel will not detect every marker in every cell.

    Mixed groups and cell types without a supported simulation program remain physical context, as do cells with
    conflicting lineage markers. Macrophages use a generic provisional class: **SPP1 is absent from this panel**, so
    that subtype is not inferred. Review the detailed immune composition: “immune” does not necessarily mean CD8.

    ✏️ To revise an annotation, edit `ANNOTATION_OVERRIDES` and rerun this cell **and every following cell**.
    Example: `{14: {"label":"Unresolved", "compartment":"unknown", "runtime_label":None}}`.
    ''')
    nb.code('''
    #@title Editable annotations and marker review
    ANNOTATION_OVERRIDES = {}
    annotation_table=prep.annotate(sample,ANNOTATION_OVERRIDES)
    display(annotation_table[['cluster','n_cells','label','compartment','runtime_label','markers','top_genes','mixed_marker_fraction']])
    fig=prep.plot_marker_review(sample)
    display(fig);plt.close(fig)
    ''')
    nb.md('''
    ## 5 · Calculate 20 AUCell scores across the entire sample
    Uses [decoupler.mt.aucell 2.3.0](https://decoupler.scverse.org/en/latest/api/generated/decoupler.mt.aucell.html)
    and 20 human GO:BP definitions from MSigDB 2026.1, on **original counts** (control probes and codewords excluded).
    The ranking universe is the whole assayed panel; the default cutoff is **100 ranks (~2% of the panel)**. Cells
    with fewer than 100 detected genes keep missing scores, so zero-count genes never enter a scored cell's top ranks.

    Outputs: raw AUCell, percentiles within provisional type, and percentiles within type and detection-depth decile
    (the latter is what the simulated cells see). **Percentiles across different types do not measure absolute
    activity differences**; these are relative RNA-enrichment scores, not demonstrated protein activity.
    ''')
    nb.code('''
    #@title Calculate AUCell, reference distributions, and quality checks
    from virtualtissue_xenium.own_data import require_tumour_interface, preflight_regions
    require_tumour_interface(sample)
    # Fail before scoring the whole section if no usable region exists.
    preflight_regions(prep,sample,n_options=N_REGIONS,core_target=CORE_TARGET,
                      max_total=MAX_TOTAL_CELLS,halo_um=HALO_UM)
    coverage=prep.score_aucell(sample,n_up=AUCELL_TOP_RANKS,min_overlap=5,batch_size=4096)
    display(coverage[['full_genes','panel_genes','panel_coverage','available']])
    print('Cells with AUCell:',int(sample['obs'].aucell_available.sum()))
    print('Cells without AUCell:',int((~sample['obs'].aucell_available).sum()))
    fig,ax=plt.subplots(figsize=(8,3))
    ax.hist(sample['obs'].n_genes_by_counts,bins=100,color='#397f9c')
    ax.axvline(AUCELL_TOP_RANKS,color='#bc582d',label='Minimum detected genes for AUCell')
    ax.set(xlabel='Detected genes per cell',ylabel='Number of cells');ax.legend()
    fig.tight_layout();fig.savefig(RUN_DIR/'rna_depth_qc.png',dpi=150)
    display(fig);plt.close(fig)
    ''')
    nb.md('''
    ## 6 · Find candidate regions with mixed populations
    Each region is a spatial disc, not a hand-picked list of cell types. Automatic proposals require at least 15%
    tumour-candidate epithelium, 8% fibroblasts and 15% immune cells in the central area; ranking favours about
    35%/20%/30% plus ~10% other stroma. Options overlap little (cell-set Jaccard ≤ 0.30) and are named
    **R01, R02, …** in order of balance score. The grey halo is surrounding context and counts toward the
    1,000-cell limit. Geometry comes from the original simplified 10x polygons.
    ''')
    nb.code('''
    #@title Read segmentation and propose regions
    # Geometry was checked before AUCell; refresh score-dependent summaries now.
    region_table=prep.propose_regions(sample,n_options=N_REGIONS,core_target=CORE_TARGET,
                                     max_total=MAX_TOTAL_CELLS,halo_um=HALO_UM)
    display(region_table[['region_id','n_core','n_halo','n_total','tumor_candidate','fibroblast','immune','stroma','cd8_provisional','macrophage_provisional','supported_qc_core','aucell_fraction','x_um','y_um']])
    fig=prep.plot_overview(sample);display(fig);plt.close(fig)
    ''')
    nb.code('''
    #@title Region gallery — inspect and compare
    fig=prep.plot_regions(sample);display(fig);plt.close(fig)
    ''')
    nb.md('''
    ## 7 · Browse a region (optional)
    Use the dropdown to enlarge any option and see its cell-type composition. You can also enter a centre in
    micrometres (read it from the whole-section map) and click **Preview manual region**; it becomes the
    region `MANUAL`. A manual choice still respects the cell limit and warns if it fails the composition criteria.

    This browser only *shows* regions. **The region that is exported is the one you type in step 8.**
    ''')
    nb.code('''
    #@title Browse regions
    import ipywidgets as widgets
    try:
        from google.colab import output
        output.enable_custom_widget_manager()
    except ImportError:
        pass
    selector=widgets.Dropdown(options=[r['region_id'] for r in sample['regions']],value=sample['regions'][0]['region_id'],description='Region:')
    preview=widgets.Output()
    def show_selected(change=None):
        with preview:
            clear_output(wait=True)
            r=next(r for r in sample['regions'] if r['region_id']==selector.value)
            if not r['balance_pass']:
                display(Markdown('**This manual region does not meet all three minimum composition criteria.**'))
            display(Markdown(f"**{r['region_id']}: {r['n_total']} cells = {r['n_core']} central + {r['n_halo']} surrounding context**"))
            fig=prep.plot_regions(sample,[selector.value],save=False);display(fig);plt.close(fig)
            core_obs=sample['obs'].iloc[r['core']]
            detail=core_obs.groupby(['cell_type_provisional','compartment'],observed=True).agg(
                cells=('internal_id','size'),with_AUCell=('aucell_available','sum'),marker_conflicts=('mixed_marker_flag','sum'))
            display(detail)
    selector.observe(show_selected,names='value')
    x_input=widgets.FloatText(value=float(region_table.iloc[0].x_um),description='X (µm):')
    y_input=widgets.FloatText(value=float(region_table.iloc[0].y_um),description='Y (µm):')
    manual_button=widgets.Button(description='Preview manual region',layout=widgets.Layout(width='230px'))
    manual_status=widgets.Output()
    def choose_manual(_):
        with manual_status:
            clear_output(wait=True)
            try:
                r=prep.manual_region(sample,x_input.value,y_input.value,CORE_TARGET,MAX_TOTAL_CELLS,HALO_UM)
                selector.options=[r['region_id'] for r in sample['regions']]
                selector.value='MANUAL';show_selected()
            except Exception as exc:
                print(str(exc))
    manual_button.on_click(choose_manual)
    display(widgets.VBox([selector,widgets.HBox([x_input,y_input]),manual_button,manual_status,preview]))
    show_selected()
    ''')
    nb.md('''
    ## 8 · ✏️ Choose the region and export the bundle ZIP
    **Set `REGION_TO_EXPORT`** to the region you want to simulate — e.g. `R03` — or `MANUAL` after previewing a
    manual centre in step 7. Rerun only this cell to export another region (each export gets its own folder).

    **TUMOR_MODE** · `exploratory_candidates`: tumour-candidate epithelium acts as `Tumor` (an explicit,
    recorded assumption). `context_until_review`: those cells stay as context. No cell is automatically labelled
    proliferating tumour.

    The export contains `region.h5ad` (raw counts, XY, AUCell), original cell/nucleus boundaries, `atlas.json.gz`,
    `runtime_contract.json`, provenance and checksums — see
    [docs/02_part2_input_bundle.md](''' + DOCS + '''02_part2_input_bundle.md). A static 10-minute step checks format,
    IDs, geometry and numerical initialisation (not biology). **The ZIP is what you give to Part 2.**
    ''')
    nb.code('''
    #@title ✏️ Export the chosen region
    REGION_TO_EXPORT = "R01"  #@param ["R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "MANUAL"] {allow-input: true}
    TUMOR_MODE = "exploratory_candidates"  #@param ["exploratory_candidates", "context_until_review"]
    from prepared_runtime import validate_bundle
    available=[r['region_id'] for r in sample['regions']]
    selected_id=REGION_TO_EXPORT.strip().upper()
    if selected_id not in available:
        raise ValueError(f'{selected_id} is not available. Choose one of {available} '
                         '(MANUAL exists only after "Preview manual region" in step 7).')
    export_stamp=datetime.now(timezone.utc).strftime('%H%M%S_%f')
    bundle=prep.export_region(sample,selected_id,tumor_mode=TUMOR_MODE,
                              destination=RUN_DIR/f'simulation_input_{selected_id}_{export_stamp}')
    print('Verifying initialization and one static step...')
    report=validate_bundle(bundle,smoke_steps=1)
    prep.write_bundle_manifest(bundle)
    zip_path=prep.zip_bundle(bundle)
    display(pd.DataFrame([{k:v for k,v in report.items() if k not in ['fidelity','scope']}]))
    print('Prepared folder:',bundle)
    print('ZIP bundle:',zip_path)
    BUNDLE_ZIP_FOR_PART2=str(zip_path)
    display(Markdown(f"**Next:** open `Part2-4_Run_Jev_Simulation.ipynb` and set `BUNDLE_ZIP` to  \\n`{zip_path}`  \\n"
                     f"(or simply to the folder `{RUN_DIR}` — the newest bundle inside is used)."))
    ''')
    nb.code('''
    #@title Optional: also download the ZIP to your computer
    DOWNLOAD_ZIP = False #@param {type:"boolean"}
    if DOWNLOAD_ZIP and IN_COLAB:
        from google.colab import files
        files.download(str(zip_path))
    else:
        print('The ZIP is saved at:',zip_path)
    ''')
    nb.md('''
    ## Next step
    Open **`Part2-4_Run_Jev_Simulation.ipynb`**. It loads the ZIP, verifies every checksum, shows the measured t0
    tissue, estimates the first round of Jev requests (without calling the API), runs the simulation and opens
    the interactive replay.

    Data: 10x Genomics, experiment 2, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Gene sets: MSigDB
    C5:GO:BP 2026.1.Hs (20 processes). [AUCell method](https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html).
    ''')
    return nb.write()


def build_part2_4() -> Path:
    nb = Notebook("Part2-4_Run_Jev_Simulation.ipynb")
    nb.md('''
    # Parts 2–4 · Run a Jev virtual-tissue simulation

    | Part | In this notebook | ✏️ What you set |
    |---|---|---|
    | **2 · Input** | load the bundle ZIP, verify it, look at the measured tissue | `INPUT_SOURCE`, `BUNDLE_ZIP`, results folder |
    | **3 · Run** | scenario, duration, seed, API limits, preflight estimate, Jev key, run | `SCENARIO`, `HORIZON_HOURS`, `MAX_HTTP_ATTEMPTS`, … , `START_RUN` |
    | **4 · Results** | interactive replay, summary figures, download | — |

    **You need** a prepared bundle ZIP (`simulation_input_<REGION>_<time>.zip`) from Part 1, and for real runs a
    **Jev (TypeSafe) API key** stored in Colab **Secrets** as `TYPESAFE_API_KEY` (or `JEV_API_KEY`).
    Without a key you can still run the free **local preview** (programmed rules, no API).

    **Colab runtime:** Python 3 · CPU. Run the cells in order. Guides:
    [Part 2](''' + DOCS + '''02_part2_input_bundle.md) ·
    [Part 3](''' + DOCS + '''03_part3_run_simulation.md) ·
    [Part 4](''' + DOCS + '''04_part4_results.md) ·
    [Jev key](''' + DOCS + '''07_jev_api_key.md)

    The experiment is an **uncalibrated mechanism simulation**, not a patient prediction. Measured t0 RNA scores
    inform local choices; fields, kinetic rates, recognition and effector competence are assumptions. There is no
    automatic fallback from a failed Jev call to the local policy.
    ''')
    nb.md('''
    ## 1 · Get the code and install dependencies
    If a restart is requested, choose **Runtime → Restart session** and rerun this cell. No Jev request is sent here.
    ''')
    nb.code(setup_code("colab_setup") + '''
    if not SETUP_READY:
        raise RuntimeError("Restart the session and rerun this setup cell before continuing.")
    import simulation_workflow as workflow
    print("Setup complete. Simulation engine:", workflow.VERSION, "| No Jev requests have been sent.")
    ''', hidden=True)
    nb.md('''
    # PART 2 · Provide and verify the input
    The input is **one bundle ZIP** produced by Part 1 (`export_region` → `zip_bundle`). Keep it compressed. It
    contains the region's cells (≤ 1,000 including context), their polygons, contacts, provisional identities and
    t0 AUCell evidence, plus the runtime contract and checksums. The simulator verifies every checksum and refuses
    any bundle whose embedded runtime differs from this engine; Python files inside the ZIP are never executed.
    Using your own data? Run Part 1 on your Xenium output first —
    see [docs/06_your_own_xenium_data.md](''' + DOCS + '''06_your_own_xenium_data.md).

    ✏️ **INPUT_SOURCE**
    * `google_drive` — `BUNDLE_ZIP` is a path in Drive. It may be the ZIP itself or the Part 1 `run_…` folder
      (then the newest `simulation_input_*.zip` inside is used). **To resume, choose the exact original ZIP**
      so a later export is not selected by mistake.
    * `upload` — a file chooser appears; pick the ZIP on your computer.
    * `url` — download the ZIP from `BUNDLE_URL` (for example a GitHub release asset).

    ✏️ **SAVE_RESULTS_TO_DRIVE** keeps a copy of the run in `DRIVE_RESULTS_DIR`, refreshed after every simulated
    hour and at the end. If the runtime is recycled, rerun with the **same input ZIP, RUN_NAME, settings and
    software environment**: the latest successful copy is restored and the responses it contains are reused.
    Paid answers since that backup may be lost and requested again. Before the first successful backup there
    may be no run to restore; failed backups can extend the gap. Saving this notebook does not save session files.

    **Local Jupyter:** keep `INPUT_SOURCE = google_drive` (here it means a file path), replace `BUNDLE_ZIP`
    with your local ZIP or output folder, and **uncheck `SAVE_RESULTS_TO_DRIVE`** or change `DRIVE_RESULTS_DIR`
    to your local backup folder. The `/content/drive/...` defaults are not converted automatically.
    The `upload` chooser requires Colab. Working files are stored in `vtx_session/`.
    ''')
    nb.code('''
    #@title ✏️ Input bundle and results folder
    INPUT_SOURCE = "google_drive"  #@param ["google_drive", "upload", "url"]
    BUNDLE_ZIP = "/content/drive/MyDrive/VirtualTissue/Xenium_Lung5K_Prepared"  #@param {type:"string"}
    BUNDLE_URL = ""  #@param {type:"string"}
    SAVE_RESULTS_TO_DRIVE = True  #@param {type:"boolean"}
    DRIVE_RESULTS_DIR = "/content/drive/MyDrive/VirtualTissue/Simulation_Results"  #@param {type:"string"}

    import urllib.request
    # Working folders on the temporary session disk (fast, local). Drive keeps the copies.
    SESSION_ROOT = Path("/content/Xenium_Jev") if IN_COLAB else Path("vtx_session").resolve()
    INPUT_DIR = SESSION_ROOT / "inputs"
    RESULT_ROOT = SESSION_ROOT / "results"
    for folder in (INPUT_DIR, RESULT_ROOT):
        folder.mkdir(parents=True, exist_ok=True)
    if IN_COLAB and (SAVE_RESULTS_TO_DRIVE or INPUT_SOURCE == "google_drive"):
        vtx_colab.mount_drive()
    if INPUT_SOURCE == "google_drive":
        source = Path(BUNDLE_ZIP)
        if source.is_dir():
            candidates = sorted(source.rglob("simulation_input_*.zip"), key=lambda p: p.stat().st_mtime)
            if not candidates:
                raise FileNotFoundError(f"No simulation_input_*.zip under {source}.")
            source = candidates[-1]
            print("Newest bundle in that folder:", source)
        if not source.is_file() or source.suffix != ".zip":
            raise FileNotFoundError(f"BUNDLE_ZIP must be a .zip file or a folder that contains one: {source}")
        INPUT_ZIP = INPUT_DIR / source.name
        if not INPUT_ZIP.is_file() or INPUT_ZIP.stat().st_size != source.stat().st_size:
            shutil.copyfile(source, INPUT_ZIP)
    elif INPUT_SOURCE == "upload":
        from google.colab import files
        previous_directory = Path.cwd()
        try:
            os.chdir(INPUT_DIR)
            uploaded = files.upload()
        finally:
            os.chdir(previous_directory)
        if len(uploaded) != 1:
            raise ValueError("Select exactly one prepared simulation ZIP.")
        name, content = next(iter(uploaded.items()))
        if not name.lower().endswith(".zip"):
            raise ValueError("Select the simulation_input ZIP file.")
        INPUT_ZIP = INPUT_DIR / Path(name).name
        INPUT_ZIP.write_bytes(content)
        del content, uploaded
    elif INPUT_SOURCE == "url":
        if not BUNDLE_URL.startswith("https://"):
            raise ValueError("BUNDLE_URL must be an https:// link to a simulation_input ZIP.")
        INPUT_ZIP = INPUT_DIR / Path(BUNDLE_URL.split("?")[0]).name
        if not INPUT_ZIP.is_file():
            urllib.request.urlretrieve(BUNDLE_URL, INPUT_ZIP)
    else:
        raise ValueError("Unknown INPUT_SOURCE.")
    BACKUP_ROOT = Path(DRIVE_RESULTS_DIR) if SAVE_RESULTS_TO_DRIVE else None
    BUNDLE = workflow.extract_bundle(INPUT_ZIP, SESSION_ROOT / "extracted")
    print("Verified:", BUNDLE["atlas"]["patch_id"])
    print("Input SHA-256:", BUNDLE["zip_sha256"])
    print("Working copy of the input:", INPUT_ZIP)
    print("Results (session disk):", RESULT_ROOT)
    print("Results copy in Drive:", BACKUP_ROOT or "disabled — download the results ZIP in Part 4")
    ''')
    nb.md('''
    ### What is inside the bundle
    `atlas.json.gz` holds one record per cell: `role` (`roi` = acting cell, `occluded` = context inside the region,
    `halo` = surrounding context), polygon `geometry`, `contacts` (polygon distance ≤ 1 µm), `identity`
    (provisional label) and `signatures` (20 AUCell scores with within-type/depth percentiles).
    `region.h5ad` holds the raw counts and XY for the same cells. Full schema:
    [docs/02_part2_input_bundle.md](''' + DOCS + '''02_part2_input_bundle.md).
    ''')
    nb.code('''
    #@title Inspect the bundle and the measured t0 tissue (no API calls)
    import json
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch, Patch
    from shapely.geometry import shape
    from IPython.display import display
    atlas = BUNDLE["atlas"]
    roles = Counter(c["role"] for c in atlas["cells"])
    print(f"{atlas['patch_id']}: {len(atlas['cells'])} cells · roles {dict(roles)}")
    print("Tumour handling:", atlas["provenance"].get("tumor_mode"), "| score method:", atlas["provenance"].get("score_method"))
    print("Signatures (t0 AUCell):", ", ".join(sorted(atlas["signatures"])))
    files = json.loads((BUNDLE["folder"] / "BUNDLE_MANIFEST.json").read_text())["files"]
    display(pd.DataFrame([dict(file=k, kilobytes=round(v["bytes"] / 1024, 1)) for k, v in files.items()
                          if not k.startswith("xenium_runtime/")]).set_index("file"))
    acting = Counter((c.get("identity") or {}).get("candidates", [{}])[0].get("label") for c in atlas["cells"] if c["role"] == "roi")
    display(pd.DataFrame(sorted(acting.items()), columns=["acting cell type", "cells"]).set_index("acting cell type"))
    labels = sorted(acting)
    colors = {label: plt.get_cmap("tab20")(i / max(1,len(labels)-1)) for i,label in enumerate(labels)}
    fig, ax = plt.subplots(figsize=(10,9))
    for cell in atlas["cells"]:
        active = cell["role"] == "roi"
        identity = cell.get("identity") or {}
        label = (identity.get("candidates") or [{"label":"context"}])[0]["label"]
        geometry = shape(cell["geometry"])
        parts = [geometry] if geometry.geom_type == "Polygon" else geometry.geoms
        for polygon in parts:
            from shapely.geometry.polygon import orient
            polygon = orient(polygon, sign=1.0)
            rings = [polygon.exterior, *polygon.interiors]
            vertices, codes = [], []
            for ring in rings:
                points = list(ring.coords)
                vertices.extend(points)
                codes.extend([MplPath.MOVETO] + [MplPath.LINETO]*(len(points)-2) + [MplPath.CLOSEPOLY])
            ax.add_patch(PathPatch(MplPath(vertices,codes), facecolor=colors.get(label,"#999999"),
                                  alpha=.85 if active else .22, edgecolor="white", linewidth=.15))
    ax.autoscale_view(); ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set(xlabel="X (µm)", ylabel="Y (µm)", title=f"{atlas['patch_id']} · measured segmentation · provisional actor labels")
    ax.legend(handles=[Patch(color=colors[label],label=label) for label in labels] + [Patch(color="#cccccc", label="Context")],
              loc="upper left", bbox_to_anchor=(1.01,1), frameon=False)
    plt.tight_layout(); plt.show()
    MEASURED_T0_FIGURE = fig
    ''')
    nb.md('''
    # PART 3 · Configure and run the simulation

    **To test CD8-mediated tumour killing, explicitly select `declared_recognition`.** Compare with
    `recognition_off` using the same region, seed, model, duration and population ceiling. The default
    `baseline_unknown` does not enable immune killing, even with a longer horizon. Hypoxic necrosis is separate.

    **Scenario**
    * `baseline_unknown` — recognition and effector/clearance competence stay unknown; cytotoxic killing is not licensed.
    * `declared_recognition` — assume antigen recognition and effector/clearance competence. Enables contact-dependent
      actions; it does not guarantee killing or infer recognition from RNA.
    * `recognition_off` — same declared competences with recognition explicitly false (negative control).

    **Limits** (the run stops at whichever comes first, and still saves a valid replay up to that point)
    * `HORIZON_HOURS` — simulated time (multiple of 10 min, max 7 days). Decisions are checkpointed every 10 min.
    * `MAX_HTTP_ATTEMPTS` — every POST to Jev counts, including failures, across resumptions. Each POST carries
      up to 20 cell questions (≤ 40 KB). **This is your cost control.** It is not a dollar cap.
    * `MAX_CELL_DECISIONS` — total cell questions; if a round would exceed it, a seeded subset is asked.
    * `MAX_POPULATION` — body ceiling including dead/removed bodies (default 1,000). The starting region already
      uses most of it (e.g. R03: 989 → only 11 new cells). Raise it to study growth.

    `RUN_NAME` names the results folder; for your own sample, use a name such as `MySample_R01_jev01`.
    Keep the name and every setting unchanged to **resume**; use a new name
    for a different scenario, duration or seed. Rough scale for R03: the first round asks 280 questions in 26 requests;
    a 24 h validation run with mocked answers needed 122 requests. Real Jev choices may need more, so start with a
    small cap and resume with a larger one (saved answers are reused).
    ''')
    nb.code('''
    #@title ✏️ Experiment settings
    RUN_NAME = "R03_jev_run01" #@param {type:"string"}
    RUN_MODE = "Jev" #@param ["Jev", "Local preview (no API)"]
    SCENARIO = "baseline_unknown" #@param ["baseline_unknown", "declared_recognition", "recognition_off"]
    JEV_MODEL = "jev-1.14.0" #@param ["jev-1.14.0", "jev-1.13.0"]
    HORIZON_HOURS = 24 #@param {type:"number"}
    SEED = 21 #@param {type:"integer"}
    MAX_POPULATION = 1000 #@param {type:"integer"}
    MAX_HTTP_ATTEMPTS = 600 #@param {type:"integer"}
    MAX_CELL_DECISIONS = 10000 #@param {type:"integer"}
    import re
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", RUN_NAME):
        raise ValueError("RUN_NAME must use letters, numbers, underscores, or hyphens.")
    if type(MAX_HTTP_ATTEMPTS) is not int or MAX_HTTP_ATTEMPTS < 1:
        raise ValueError("The HTTP limit must be a positive integer.")
    # Whole hours are stored as integers so that the run configuration (and resuming) is stable.
    horizon_hours = int(HORIZON_HOURS) if float(HORIZON_HOURS).is_integer() else HORIZON_HOURS
    SETTINGS = dict(model=JEV_MODEL, scenario=SCENARIO, seed=SEED, horizon_hours=horizon_hours,
                    max_decisions=MAX_CELL_DECISIONS, max_population=MAX_POPULATION)
    MODE = "jev" if RUN_MODE == "Jev" else "local_preview"
    RUN_DIR = RESULT_ROOT / (RUN_NAME + ("_local_preview" if MODE == "local_preview" else ""))
    workflow.effective_contract(BUNDLE, SETTINGS)
    print("Output folder:", RUN_DIR)
    print("Run mode:", RUN_MODE)
    ''')
    nb.md('''
    ### Preflight: check before spending API credits
    Makes **zero API calls**: one static mechanical step, the number of questions and HTTP requests of the
    **first** round, and a sample request saved as JSON so you can see exactly what Jev receives. Later rounds
    depend on the trajectory, so this is not a forecast of the total.

    Each question contains the cell's permitted menu, its t0 raw AUCell values and percentiles (whole-sample
    type + detection-depth reference; missing scores stay missing), local contacts and simulated state.
    ''')
    nb.code('''
    #@title Preflight estimate and request preview (no API calls)
    PREFLIGHT, FIRST_REQUEST = workflow.preflight(BUNDLE, SETTINGS)
    from virtualtissue_xenium.experiments import experiment_notes
    for note in experiment_notes(SETTINGS, PREFLIGHT):
        print(note)
    display(pd.DataFrame([dict(cell_type=k, actors=v) for k,v in PREFLIGHT["acting_composition"].items()]))
    print(json.dumps({k:PREFLIGHT[k] for k in ("total_cells", "acting_cells", "context_cells", "initial_cell_questions", "initial_http_batches", "growth_headroom", "static_preflight_passed")}, indent=2))
    if MAX_HTTP_ATTEMPTS < PREFLIGHT["initial_http_batches"] and MODE == "jev":
        print("The allowance cannot cover the first round. Increase it before starting; otherwise no POST is sent.")
    PREVIEW_DIR = RESULT_ROOT / "preflight" / RUN_NAME
    PREVIEW_DIR.mkdir(parents=True,exist_ok=True)
    workflow.write_json(PREVIEW_DIR / "preflight.json", PREFLIGHT)
    workflow.write_json(PREVIEW_DIR / "first_request_preview.json", FIRST_REQUEST)
    MEASURED_T0_FIGURE.savefig(PREVIEW_DIR / "measured_t0.png", dpi=160, bbox_inches="tight")
    print("Request preview:", PREVIEW_DIR / "first_request_preview.json")
    if FIRST_REQUEST:
        q = next(iter(FIRST_REQUEST["questions"].values()))
        print("Example question — permitted options:", list(q["criteria"]))
    ''')
    nb.md('''
    ### Add your Jev key and start
    1. Click the **key icon (Secrets)** in Colab's left sidebar → **Add new secret** → Name `TYPESAFE_API_KEY`
       (or `JEV_API_KEY`) → paste the key as Value → switch **Notebook access** on.
       Locally, set the environment variable `TYPESAFE_API_KEY` instead. If no key is found, a hidden prompt appears.
       **Never paste a key into a code cell.** The key is not written to the journal, replay or results.
    2. ✏️ Check **START_RUN** and run the cell. Checking it authorises spending up to `MAX_HTTP_ATTEMPTS` requests.

    Every Jev batch is journaled before sending. There are **no automatic retries**. A whole round is validated before
    any decision is applied; API errors, malformed answers or model drift stop the run, which is still saved.
    Checkpoints are written every simulated hour. To resume after an interruption, rerun the notebook from the top
    with the same input ZIP, `RUN_NAME`, settings and software environment. The engine replays from t0 using
    the responses in the available journal and sends only missing batches. After runtime loss, responses absent
    from the last successful Drive copy may be requested and paid for again. An interrupted request whose
    entry survives in the journal has unknown billing status and is never resent automatically.
    ''')
    nb.code('''
    #@title Optional: unlock after a terminated session
    # Leave false during a normal run. Use only after stopping all other sessions writing this run.
    CONFIRM_PREVIOUS_SESSION_STOPPED = False #@param {type:"boolean"}
    if CONFIRM_PREVIOUS_SESSION_STOPPED:
        (RUN_DIR / ".running.lock").unlink(missing_ok=True)
        print("Stale lock removed. Existing journal and checkpoints were preserved.")
    ''')
    nb.code('''
    #@title ✏️ Start or resume the simulation
    START_RUN = False #@param {type:"boolean"}
    def start_selected_run():
        if not START_RUN:
            print("Ready. Review the scenario and limits, then check START_RUN and rerun this cell.")
            return None
        key = ""
        try:
            if MODE == "jev":
                key = vtx_colab.get_jev_key()
                if not key:
                    raise ValueError("No API key provided. No Jev request was sent.")
            if BACKUP_ROOT and vtx_colab.restore_run(BACKUP_ROOT, RUN_DIR):
                print("Restored this run from", BACKUP_ROOT / RUN_DIR.name, "- responses in this backup will be reused.")
            progress = vtx_colab.progress_with_backup(RUN_DIR, BACKUP_ROOT)
            return workflow.run(BUNDLE, SETTINGS, RUN_DIR, api_key=key,
                                max_http_attempts=MAX_HTTP_ATTEMPTS, mode=MODE, progress=progress)
        finally:
            key = ""
    RESULT = start_selected_run()
    if RESULT:
        print(json.dumps(RESULT["report"], indent=2))
    ''')
    nb.md('''
    # PART 4 · Explore and keep the results

    The replay supports a time slider, play/pause, cell types, current actions, t0 AUCell percentiles, simulated
    field overlays (IFN-γ, CXCL9, O₂, TGF-β) and a cell inspector with each cell's decision history and Jev
    probabilities. The t0 view uses measured boundaries; later shapes and smooth motion render the recorded
    simulated states. Full guide: [docs/04_part4_results.md](''' + DOCS + '''04_part4_results.md).

    The results ZIP contains `replay.html` (opens offline in any browser), `episode.json.gz`, the Jev journal,
    `metrics.csv` / `events.csv` / `decisions.csv`, the input provenance, the effective contract and a run report.
    ''')
    nb.code('''
    #@title Interactive replay
    report_path = RUN_DIR / "run_report.json"
    if report_path.exists():
        report = json.loads(report_path.read_text())
        print("Completed:", report["completed"], "| Stop reason:", report["stop_reason"])
        print("Saved simulation hours:", report["checkpoint_time_min"] / 60)
        print("API usage:", report["usage"])
        if report.get("provider_error"):
            print("Provider:", report["provider_error"])
        vtx_colab.show_replay(RUN_DIR / "replay.html")
    else:
        print("Run the simulation cell first. No saved result is available for this RUN_NAME.")
    ''')
    nb.code('''
    #@title Summary figures and tables
    from virtualtissue_xenium import results as vtx_results
    if report_path.exists():
        RES = vtx_results.load_results(RUN_DIR)
        display(pd.Series(vtx_results.summary(RES), name="value").to_frame())
        fig = vtx_results.plot_overview(RES)
        fig.savefig(RUN_DIR / "overview.png", dpi=150)
        plt.show()
        display(vtx_results.choice_counts(RES))
    ''')
    nb.code('''
    #@title Keep the results: Drive copy and/or download
    DOWNLOAD_RESULTS = True #@param {type:"boolean"}
    result_zip = RUN_DIR.with_name(RUN_DIR.name + "_results.zip")
    if BACKUP_ROOT and RUN_DIR.exists():
        print("Copy in Google Drive:", vtx_colab.backup_run(RUN_DIR, BACKUP_ROOT))
    if DOWNLOAD_RESULTS and result_zip.is_file() and IN_COLAB:
        from google.colab import files
        files.download(str(result_zip))
        print("Download requested. Check your browser's Downloads folder before closing this session.")
    elif not result_zip.is_file():
        print("No result ZIP exists yet. Start the simulation first, then rerun this cell.")
    else:
        print("Results ZIP:", result_zip)
    ''')
    nb.md('''
    ## Method and references
    Jev chooses among the engine's licensed actions through the TypeSafe Choice API
    ([System One API](https://docs.typesafe.ai/api/system-one), [Choice schema](https://docs.typesafe.ai/primitives/choice),
    [pinned models](https://docs.typesafe.ai/models)). The selected action is used directly; returned probabilities rank
    behaviours and are not reaction rates. The numerical engine controls movement, contacts, clocks, fields,
    cell-cycle progression and outcomes. Parameter changes are saved in `effective_contract.json` beside the original
    contract fingerprint. To revisit saved results later without re-running anything, use
    `Part4_Explore_Saved_Results.ipynb`.
    ''')
    return nb.write()


def build_part4() -> Path:
    nb = Notebook("Part4_Explore_Saved_Results.ipynb")
    nb.md('''
    # Part 4 · Explore saved simulation results

    Opens one or more results (`<RUN_NAME>_results.zip` or the run folder) **without calling the API and without
    re-running the simulation**: interactive replay, headline numbers, cell populations over time, events,
    decisions per cell type, the decision history of a single cell, and a comparison table across runs.

    Runs on any Colab runtime (CPU) or locally. Guide: [docs/04_part4_results.md](''' + DOCS + '''04_part4_results.md).
    Tip: `replay.html` inside the ZIP also opens directly in Chrome/Firefox/Safari, with no notebook at all.
    ''')
    nb.code(setup_code(None, install_default=False) + '''
    if not SETUP_READY:
        raise RuntimeError("Restart the session and rerun this setup cell before continuing.")
    from virtualtissue_xenium import results as vtx_results
    print("Ready.")
    ''', hidden=True)
    nb.md('''
    ## ✏️ Which result to open
    `RESULTS_PATH` can be a results ZIP or a run folder, in Drive (`/content/drive/MyDrive/...`) or on the session
    disk. With `upload`, choose the ZIP from your computer.

    **Local Jupyter:** leave `RESULTS_SOURCE = google_drive` and replace `RESULTS_PATH` with a local path;
    this option reads local files without mounting Drive. For example, use
    `vtx_session/results/R03_jev_run01_results.zip` when the simulation ran from the same notebook folder.
    A local preview has the suffix `_local_preview_results.zip`. The `upload` chooser requires Colab.
    ''')
    nb.code('''
    #@title ✏️ Result to open
    RESULTS_SOURCE = "google_drive"  #@param ["google_drive", "upload"]
    RESULTS_PATH = "/content/drive/MyDrive/VirtualTissue/Simulation_Results/R03_jev_run01_results.zip"  #@param {type:"string"}
    import json
    import pandas as pd
    import matplotlib.pyplot as plt
    from IPython.display import display
    WORK = Path("/content/vtx_explore") if IN_COLAB else Path("vtx_explore").resolve()
    WORK.mkdir(parents=True, exist_ok=True)
    results_path = Path(RESULTS_PATH)
    if RESULTS_SOURCE == "upload":
        from google.colab import files
        uploaded = files.upload()
        name, content = next(iter(uploaded.items()))
        results_path = WORK / Path(name).name
        results_path.write_bytes(content)
    elif RESULTS_PATH.startswith("/content/drive/") and IN_COLAB:
        vtx_colab.mount_drive()
    extract_to = WORK / results_path.stem
    shutil.rmtree(extract_to, ignore_errors=True)
    RES = vtx_results.load_results(results_path, extract_to=extract_to)
    display(pd.Series(vtx_results.summary(RES), name="value").to_frame())
    ''')
    nb.code('''
    #@title Interactive replay
    vtx_colab.show_replay(RES["replay"])
    ''')
    nb.code('''
    #@title Overview: populations, decisions, events, choices
    fig = vtx_results.plot_overview(RES)
    plt.show()
    display(vtx_results.population_by_type(RES).iloc[[0, -1]])
    display(vtx_results.choice_counts(RES))
    display(vtx_results.event_counts(RES, bin_hours=6))
    ''')
    nb.md('''
    ## ✏️ Follow one cell
    Find a cell id in the replay's inspector (click a cell) or in the decisions table, then set `CELL_ID`.
    Each row is one decision: the chosen action, Jev's probability for every permitted option, and why the
    cell was asked (`wake_reasons`).
    ''')
    nb.code('''
    #@title ✏️ Decision history of one cell
    DECISIONS = vtx_results.decision_table(RES)
    CELL_ID = 0  #@param {type:"integer"}
    cell = CELL_ID or (int(DECISIONS.cell_id.iloc[0]) if not DECISIONS.empty else 0)
    print("Cell", cell, "(type:", vtx_results.labels_by_body(RES).get(cell), ")")
    pd.set_option("display.max_colwidth", 200)
    display(vtx_results.cell_history(RES, cell))
    ''')
    nb.md('''
    ## ✏️ Compare several runs (optional)
    List result ZIPs or folders, separated by commas. Use the same region, seed, provider,
    model and settings under `declared_recognition` and `recognition_off`, or identical
    scenarios for a repeat check. Results are compared at the latest common saved time.
    Local preview and Jev runs are different providers and cannot form a matched pair.
    ''')
    nb.code('''
    #@title ✏️ Comparison table
    COMPARE_PATHS = ""  #@param {type:"string"}
    paths = [p.strip() for p in COMPARE_PATHS.replace(",", "\\n").splitlines() if p.strip()]
    if paths:
        try:
            display(vtx_results.compare_runs(paths))
        except ValueError as error:
            print("Comparison unavailable:", error)
    else:
        print("Add paths to COMPARE_PATHS (comma-separated) to compare runs.")
    ''')
    return nb.write()


if __name__ == "__main__":
    for path in (build_part1(), build_part2_4(), build_part4()):
        print("wrote", path.relative_to(ROOT))
