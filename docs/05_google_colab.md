# Google Colab, step by step

This guide covers everything that is specific to Colab: where the code comes from, which runtime to choose, every path you may need to edit, where the Jev key goes, how Google Drive is used, and what happens when Colab disconnects.

---

## 1. Open a notebook

Use any of these:

* **From GitHub (once the repository is public):** go to [colab.research.google.com](https://colab.research.google.com) and choose **File → Open notebook → GitHub**. Type `hylixinsights/VirtualTissue-Xenium` and click a notebook under `notebooks/`. Or open these links directly:
  - `https://colab.research.google.com/github/hylixinsights/VirtualTissue-Xenium/blob/main/notebooks/Part1_Prepare_Xenium_Region.ipynb`
  - `https://colab.research.google.com/github/hylixinsights/VirtualTissue-Xenium/blob/main/notebooks/Part2-4_Run_Jev_Simulation.ipynb`
  - `https://colab.research.google.com/github/hylixinsights/VirtualTissue-Xenium/blob/main/notebooks/Part4_Explore_Saved_Results.ipynb`
* **From your computer:** choose **File → Upload notebook** and select the `.ipynb` file.
* **From Google Drive:** double-click the `.ipynb` in Drive, then choose **Open with → Google Colaboratory**.

A notebook opened from GitHub is a read-only copy. Use **File → Save a copy in Drive** if you want to keep your edits and outputs.

For the public example, start with Part 1 and keep its default data paths. **For your own sample, run Part 1 setup and then pause before step 2:** follow the [own-data guide](06_your_own_xenium_data.md) to mount Drive, prepare your files and review your cluster map. The GitHub setup cell itself does not mount Drive. A Jev key is needed only for a paid run in the simulation notebook.

## 2. Choose the runtime

**Runtime → Change runtime type → Runtime type: Python 3 · Hardware accelerator: CPU.** No GPU is used.

* The first cell prints the Python version that is **actually** running. The engine supports **3.12 and 3.13**. If Colab's default ever moves to a newer Python, choose an earlier **Runtime version** in the same dialog. The setup cell stops with a clear message if the version is unsupported.
* Standard RAM is enough. Part 1 peaks at about 4 GB on the public example; Parts 2–4 need less than 2 GB.

## 3. The first cell: code and dependencies

Every notebook starts with the same hidden cell (*Get the VirtualTissue-Xenium code…*). Double-click its title to see the code.

| ✏️ Field | Default | When to change it |
|---|---|---|
| `CODE_SOURCE` | `github` | `google_drive` if the repository is not public yet, or if you use a modified copy. `installed` only for local Jupyter with `pip install` |
| `GITHUB_URL`, `GITHUB_BRANCH` | the official repository, `main` | To use a fork or another branch/tag |
| `DRIVE_CODE_DIR` | `/content/drive/MyDrive/VirtualTissue-Xenium` | Where you uploaded the repository folder: the folder that contains `src/`, `notebooks/` and `README.md` |
| `INSTALL_DEPENDENCIES` | `True` | Leave it on in Colab |

What it does:

1. It clones the repository, or copies it from Drive, to **`/content/VirtualTissue-Xenium`** on the temporary session disk. Code is never run from Drive.
2. It installs the exact versions from `src/virtualtissue_xenium/engine/requirements-colab.txt`, which takes about 1–2 minutes.
3. If an already-loaded library changed, it prints **ONE SESSION RESTART IS REQUIRED**. Choose **Runtime → Restart session** (not *Disconnect and delete*) and run the first cell again. It must end with *Ready* / *Setup complete*.

**Uploading the code to Drive (for `google_drive`):** download or clone the repository and upload the whole `VirtualTissue-Xenium` folder to `MyDrive/`. Uploading through drive.google.com works. Keep the folder structure intact.

## 4. Google Drive layout and every path

The notebooks use one folder in your Drive, `MyDrive/VirtualTissue/`. Colab asks once for permission to mount Drive.

```text
MyDrive/
├── VirtualTissue-Xenium/                       code — only for CODE_SOURCE = google_drive
└── VirtualTissue/
    ├── Xenium_Public_Lung5K/                   ◀ Part 1 DATA_DIR (filled automatically by the download)
    │   ├── SOURCE.json
    │   ├── DOWNLOAD_MANIFEST.json
    │   └── raw/                                10 original 10x files, ≈ 503 MiB
    ├── Xenium_Lung5K_Prepared/                 ◀ Part 1 OUTPUT_BASE
    │   └── run_20261025_203601_120520/         one folder per Part 1 run (UTC time)
    │       ├── candidate_regions.png, whole_section_regions.png, …
    │       ├── simulation_input_R03_203812_123456/
    │       └── simulation_input_R03_203812_123456.zip   ◀ Part 2 BUNDLE_ZIP
    └── Simulation_Results/                     ◀ Part 2–4 DRIVE_RESULTS_DIR
        ├── R03_jev_run01/                      copy of the run folder (refreshed every simulated hour)
        └── R03_jev_run01_results.zip           ◀ Part 4 RESULTS_PATH
```

### All the fields you may edit

| Notebook | ✏️ Field | Default | Meaning |
|---|---|---|---|
| Part 1 | `DATA_DIR` | `/content/drive/MyDrive/VirtualTissue/Xenium_Public_Lung5K` | Folder with `raw/` + manifests. **Not** `raw/` itself |
| Part 1 | `OUTPUT_BASE` | `/content/drive/MyDrive/VirtualTissue/Xenium_Lung5K_Prepared` | Must be outside `DATA_DIR` |
| Part 1 | `DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING` | ✔ | Download the example into an empty `DATA_DIR` |
| Part 1 | `CLUSTER_ANNOTATIONS_JSON` | empty | Your cluster map (own data only) |
| Part 1 | `MAX_TOTAL_CELLS`, `CORE_TARGET`, `HALO_UM`, `N_REGIONS`, `AUCELL_TOP_RANKS`, `VERIFY_INPUT_HASHES` | 1000, 750, 30, 8, 100, ✔ | Region size and scoring ([01](01_part1_prepare_region.md#parameters-step-3)) |
| Part 1 | `ANNOTATION_OVERRIDES` | `{}` | Session-only relabelling of clusters |
| Part 1 | **`REGION_TO_EXPORT`** | `R01` | **The region to simulate: `R01`…`R08` or `MANUAL`.** Type `R03` here |
| Part 1 | `TUMOR_MODE` | `exploratory_candidates` | Or `context_until_review` |
| Part 2 | `INPUT_SOURCE` | `google_drive` | Or `upload`, `url` |
| Part 2 | **`BUNDLE_ZIP`** | `/content/drive/MyDrive/VirtualTissue/Xenium_Lung5K_Prepared` | A bundle ZIP, or a folder (newest bundle inside is used) |
| Part 2 | `BUNDLE_URL` | empty | For `INPUT_SOURCE = url` |
| Part 2 | `SAVE_RESULTS_TO_DRIVE`, `DRIVE_RESULTS_DIR` | ✔, `/content/drive/MyDrive/VirtualTissue/Simulation_Results` | Keep a resumable copy in Drive |
| Part 3 | `RUN_NAME`, `RUN_MODE`, `SCENARIO`, `JEV_MODEL`, `HORIZON_HOURS`, `SEED`, `MAX_POPULATION`, **`MAX_HTTP_ATTEMPTS`**, `MAX_CELL_DECISIONS` | see [03](03_part3_run_simulation.md#1-experiment-settings) | The experiment and its limits |
| Part 3 | **`START_RUN`** | ☐ | Must be checked to send any request |
| Part 3 | `CONFIRM_PREVIOUS_SESSION_STOPPED` | ☐ | Only to remove a stale lock |
| Part 4 | `DOWNLOAD_RESULTS` | ✔ | Download the results ZIP at the end |
| Explorer | `RESULTS_SOURCE`, **`RESULTS_PATH`**, `CELL_ID`, `COMPARE_PATHS` | — | Which results to open |

You can use any Drive folders you like. Just change the defaults consistently. Paths are case-sensitive. Files on the session disk, such as `/content/Xenium_Jev/`, disappear when the runtime is recycled. Files in a successfully mounted Drive persist.

**Shared drives:** paths look like `/content/drive/Shareddrives/<drive name>/…`. **Shortcuts** to folders shared with you are not always followed by Colab, so copy the data into your own Drive instead.

## 5. The Jev API key (Colab Secrets)

1. Open your copy of **`Part2-4_Run_Jev_Simulation.ipynb` first**, then click the **key icon 🔑** in the left sidebar (**Secrets**).
2. Click **+ Add new secret**. Enter Name `TYPESAFE_API_KEY` (the name `JEV_API_KEY` is also accepted) and paste the key as the Value.
3. Switch **Notebook access** on for this notebook. Each new notebook needs its own access toggle.

The notebook reads the secret only inside the run cell and never prints or saves it. Never paste a key into a code cell or a form field: notebooks are often shared, and their outputs are saved in Drive. More in [07_jev_api_key.md](07_jev_api_key.md).

## 6. Running safely

* In a fresh copy of the simulation notebook, **Run all** (`Runtime → Run all`) stops short of sending requests because `START_RUN` is unchecked by default. In a saved copy, review the checkbox first: if you left it checked, **Run all** can start or resume a paid run. Check the box and run that cell when you are ready.
* Cells with a form (`#@title`) hide their code. Change values in the form on the right, then run the cell.
* Long runs: keep the browser tab open. Free Colab sessions can disconnect or be recycled. `SAVE_RESULTS_TO_DRIVE` preserves periodic copies, but responses recorded since the last successful backup can be lost (§7). Check the output for backup-failure messages.
* `files.upload()` is slow for files over ~100 MB. Prefer Drive for inputs.

## 7. When Colab disconnects

| What you see | What to do |
|---|---|
| "Runtime disconnected", **same** runtime reconnects (files under `/content` still there) | Rerun the setup, input and settings cells, then the run cell with the same `RUN_NAME`: the run resumes from the journal |
| A **new** runtime (files under `/content` gone) | Rerun the notebook with the exact original bundle ZIP, `RUN_NAME` and settings. The run cell restores the latest successful copy from `DRIVE_RESULTS_DIR/<RUN_NAME>` and reuses the responses it contains. Check for the `Restored this run from …` message |
| "This run is locked by another or interrupted session" | Make sure no other tab is running it, check `CONFIRM_PREVIOUS_SESSION_STOPPED` in the unlock cell, run it, then rerun the run cell |
| "This run folder belongs to another configuration" | You changed a setting, or Colab changed the Python or package versions. Use a new `RUN_NAME` |

The Drive copy is taken after each completed simulated hour and at the end. **Only responses present in that copy can be reused in a new runtime.** If the runtime disappears between successful hourly backups, the responses since the last copy can be lost, and resuming may request and pay for those batches again. Before the first successful backup there may be no saved run. If a backup fails, the gap can exceed one simulated hour. A browser disconnect that leaves the runtime and its files intact does not by itself lose the journal.

Use the original ZIP path, not a folder that may now select a different export. Keep the same Python/package versions and CPU architecture; otherwise the engine refuses to resume. A failed or unresolved request is not resent if its journal entry was preserved. An entry lost with the runtime cannot provide that protection. The recovered request count excludes lost entries, so check your TypeSafe account usage after a runtime loss if spending is a concern.

## 8. Saving outputs

* Part 1 writes directly to `OUTPUT_BASE` in Drive.
* Parts 2–4 work on the fast session disk and copy to `DRIVE_RESULTS_DIR`. The last cell refreshes the copy and downloads the results ZIP (`DOWNLOAD_RESULTS`).
* **Saving the notebook does not save files under `/content`.** Only Drive copies and downloads persist.
* Open `replay.html` after downloading it. Drive's own preview does not execute it.
