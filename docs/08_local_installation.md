# Local installation and the `vtx` command line

Everything that runs in Colab also runs on macOS, Linux or Windows (WSL recommended), on a laptop or a server.

## 1. Install

Requirements: **Python 3.12 or 3.13**, `git`, about 2 GB of disk for the environment, plus ≈ 0.5 GB for the example data. Part 1 peaks at ≈ 4 GB of RAM.

```bash
git clone https://github.com/hylixinsights/VirtualTissue-Xenium.git
cd VirtualTissue-Xenium
python3.13 -m venv .venv
source .venv/bin/activate              # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -e ".[notebooks,test]"     # exact validated versions + JupyterLab + pytest
vtx --version
```

The dependencies are **pinned exactly** (`numpy==2.3.5`, `decoupler==2.3.0`, …) because the engine was validated with them and resumable runs check the environment. Use a dedicated virtual environment. With `conda`, create an environment with `python=3.13` and use `pip` inside it.

`pip install -r requirements.txt` installs the same pins without the package itself. Without installing it, put the code on the path with `sys.path.insert(0, "src")`.

## 2. The command line

```text
vtx download-example DATA_DIR                 download + verify the public example (~503 MiB)
vtx register-download DATA_DIR                register files you extracted from the full 10x ZIP
vtx verify DATA_DIR [--quick]                 re-check a data folder against its manifest
vtx own-data OUTS DATA_DIR --sample-id ID     create a data folder from your Xenium output
vtx annotation-template DATA_DIR OUT.json     draft your cluster annotation map (+ OUT.markers.csv)
vtx check-annotations OUT.json --data-dir D   validate it
vtx prepare DATA_DIR OUTPUT_BASE [--region R03 …] [--center X Y] [--list-only] [--annotations OUT.json]
vtx check BUNDLE.zip [--scenario … --hours …]           Part 2: verify + preflight (no API)
vtx run BUNDLE.zip RESULTS_DIR --run-name NAME [--mode local|jev] [--scenario …] [--hours 24]
        [--seed 21] [--max-http-attempts 600] [--max-decisions 10000] [--max-population 1000]
        [--model jev-1.14.0] [--yes]
vtx explore RESULTS.zip|RUN_DIR [--plot overview.png] [--open]
```

`python -m virtualtissue_xenium …` is equivalent to `vtx …`. Every command has `-h`.

A full session:

```bash
vtx download-example data/Xenium_Public_Lung5K
vtx prepare data/Xenium_Public_Lung5K outputs/prepared --list-only          # look at the PNG maps
vtx prepare data/Xenium_Public_Lung5K outputs/prepared --region R03
BUNDLE=$(ls -t outputs/prepared/run_*/simulation_input_R03_*.zip | head -1)
vtx check "$BUNDLE"
vtx run "$BUNDLE" outputs/runs --run-name R03_local --mode local --hours 24
export TYPESAFE_API_KEY=...
vtx run "$BUNDLE" outputs/runs --run-name R03_jev01 --mode jev --hours 24 --max-http-attempts 100   # pilot
vtx run "$BUNDLE" outputs/runs --run-name R03_jev01 --mode jev --hours 24 --max-http-attempts 600   # resume, same answers reused
vtx explore outputs/runs/R03_jev01_results.zip --plot outputs/R03_jev01.png --open
```

## 3. The notebooks in local Jupyter

```bash
jupyter lab notebooks/
```

The first cell detects that it is not running in Colab and uses the repository that contains the notebook, whatever `CODE_SOURCE` says. Choose `CODE_SOURCE = installed` if you installed the package elsewhere. Set `INSTALL_DEPENDENCIES = False` if your environment already has the pins (recommended locally). Otherwise it runs `pip install -r …/requirements-colab.txt` in the current environment.

Paths to change:

**Change these before running the data/input cells.** Local execution does not replace the `/content/drive/...` values shown in the forms. The relative paths below assume you started Jupyter with `jupyter lab notebooks/` from the repository folder and opened a notebook there; absolute paths also work.

| Notebook | Field | Local example |
|---|---|---|
| Part 1 | `DATA_DIR`, `OUTPUT_BASE` | `../data/Xenium_Public_Lung5K`, `../outputs/prepared` |
| Part 1, own data | `CLUSTER_ANNOTATIONS_JSON`, `DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING` | Your reviewed map's local path; download unchecked. Follow the [own-data setup](06_your_own_xenium_data.md) first |
| Part 2–4 | `INPUT_SOURCE = google_drive` (it simply means "a path") and `BUNDLE_ZIP` | `../outputs/prepared` |
| Part 2–4 | `SAVE_RESULTS_TO_DRIVE` | ☐, or ✔ with `DRIVE_RESULTS_DIR` = any backup folder (e.g. a synced Dropbox/OneDrive folder) |
| Explorer | `RESULTS_PATH` | `vtx_session/results/R03_jev_run01_results.zip` |

Leave the explorer's `RESULTS_SOURCE = google_drive` for a local file path too. Despite the field names, these path options work without Google Drive in local Jupyter. Do not choose `upload` locally: that chooser requires Colab. For a local preview, the results ZIP has the suffix `_local_preview_results.zip`.

Working files go to `vtx_session/` next to the notebook. The key comes from the `TYPESAFE_API_KEY` environment variable of the shell that started Jupyter, or from the hidden prompt.

## 4. Python API

```python
import virtualtissue_xenium as vtx
vtx.prepare     # xenium_prepare        – Part 1
vtx.workflow    # simulation_workflow   – Parts 2–4
vtx.runtime     # prepared_runtime      – load_inputs, initialize, validate_bundle
from virtualtissue_xenium import download, own_data, colab, results
```

Importing `virtualtissue_xenium` inserts `src/virtualtissue_xenium/engine` at the front of `sys.path`, so the engine modules can also be imported by their original names (`import simulation_workflow`). Do **not** also put an extracted bundle folder on `sys.path`, because it contains its own copy of `prepared_runtime.py` and `xenium_runtime/`. `vtx.check_engine_origin()` raises if an engine module was loaded from somewhere else.

## 5. Running the tests

```bash
pytest                                                    # unit + engine integrity tests (< 1 min)
VTX_TEST_BUNDLE=path/to/simulation_input_R03_….zip pytest # + local run / explore on a real bundle
python tests/validate_r03.py --output /tmp/r03 --zip path/to/simulation_input_R03_….zip
        # original R03 integration test: mocked Jev, budget stops, resume, failures, unsafe ZIPs
python tests/run_notebooks_mocked.py --output /tmp/nb --data-dir data/Xenium_Public_Lung5K
        # executes all three notebooks with a mocked Colab and mocked Jev (needs the example data;
        # use --bundle BUNDLE.zip instead of --data-dir to skip Part 1). PIP_NO_INDEX=1 is set,
        # so the pinned packages must already be installed.
```

| Test file | What it protects |
|---|---|
| `test_engine_integrity.py` | Engine files byte-identical to `ORIGINAL_SOURCES.json`; pyproject pins = engine requirements; one shared engine folder; stdlib-only package import |
| `test_provider.py` (original) | Journal, budget atomicity, resume without new calls, response validation, key redaction |
| `test_preparation.py` (original) | AUCell direction and rank invariance, tie handling across batches, region geometry, missing-score policy, stale-annotation refusal |
| `test_colab_setup.py` | Python 3.12/3.13 check, restart detection (both bootstraps), notebooks valid and English, region choice is an explicit parameter |
| `test_helpers.py` | Byte-range download (local server), checksum refusal, own-data folder + template + validation + custom map read by the engine, Drive backup/restore, key lookup order, CLI run/explore |
