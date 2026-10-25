# Troubleshooting

Find the message you see (or part of it) in the left column.

## Setup (first cell)

| Message | Cause | Fix |
|---|---|---|
| `This session is running Python 3.x … supports Python 3.12 and 3.13` | Unsupported interpreter | **Runtime → Disconnect and delete runtime**, then reconnect with Python 3 / CPU. If Colab's default Python changed, pick an earlier *Runtime version* in *Change runtime type* |
| `ONE SESSION RESTART IS REQUIRED` | A library that was already loaded got a different version | **Runtime → Restart session**, then run the first cell again. It happens at most once per runtime |
| `Dependency verification failed …` | A pinned package failed to install or import | Read the error at the end. Usually a temporary PyPI or network problem: run the cell again. Locally, use a fresh virtual environment |
| `fatal: repository … not found` / git clone fails | Repository not public yet, wrong `GITHUB_URL`, or no network | Use `CODE_SOURCE = google_drive` with the repository uploaded to Drive |
| `DRIVE_CODE_DIR must be the VirtualTissue-Xenium folder …` | `DRIVE_CODE_DIR` points to the wrong level | Point to the folder that contains `src/`, `notebooks/`, `README.md` |
| `Module 'xenium_runtime' was imported from …, not from …/engine` | Another copy (e.g. an extracted bundle on `sys.path`) was imported first | Restart the session; import `virtualtissue_xenium` before anything else; do not add bundle folders to `sys.path` |
| `Complete step 1 first …` | A later cell was run before setup reported Ready | Run the first cell until it says *Ready* |

## Part 1

| Message | Cause | Fix |
|---|---|---|
| `DOWNLOAD_MANIFEST.json was not found …` | `DATA_DIR` is wrong (often points to `raw/`) or the download is off | Point to the folder containing `raw/`; or check `DOWNLOAD_PUBLIC_EXAMPLE_IF_MISSING` |
| `Download failed after 5 attempts` / `does not support byte-range requests` | Network problem or a changed 10x URL | Run again (finished files are kept), or download the full ZIP and use `register_manual_download` ([01](01_part1_prepare_region.md#option-b--download-the-full-archive-yourself)) |
| `Checksum mismatch for …; the partial file was removed` | Corrupted transfer, or 10x replaced the file | Run again; if it persists, the published file changed. Report it in an issue |
| `Missing or changed input: …` / `Input checksum mismatch: …` | A raw file is missing, truncated, or was modified | Re-download or re-copy it. For own data, call `write_manifest` again only if you changed the file on purpose |
| `This data folder is sample '…', not the public example. Provide your own cluster annotation map` | Own data without `CLUSTER_ANNOTATIONS_JSON` | Follow [06](06_your_own_xenium_data.md) §3 |
| `Invalid cluster annotation map: …` | Placeholder labels, unknown compartment or runtime label, missing cluster | Fix the listed records; validate with `vtx check-annotations` |
| `NameError: name 'DATA_DIR' is not defined` in an added own-data cell | The path-definition cell was skipped or the runtime restarted | Rerun setup, then the Colab or local path block in [06](06_your_own_xenium_data.md#2-step-1--create-the-data-folder), before creation/template/validation cells |
| Your own outs folder or cluster JSON cannot be found | Drive not mounted, upload incomplete, or a relative path points elsewhere | Run the own-data path block, confirm the upload is complete, and use the same `CLUSTER_ANNOTATIONS_JSON` path for creation, editing, validation and Part 1 step 2 |
| `…/raw is not empty; choose a new data folder` | A previous copy was interrupted, or the destination belongs to another sample | Choose a new `DATA_DIR`; completed folders with a manifest can be verified and reused by the guide's creation cell |
| `Unknown vendor clusters: this annotation map belongs to one frozen sample.` | The map lacks clusters that exist in the sample | Add records for every cluster (`validate_annotations(…, data_dir=…)` lists them) |
| `Gene-expression row sums differ from vendor transcript counts.` | `cells.parquet` and the matrix come from different runs or were modified | Use files from the same, unmodified XOA output |
| `Duplicate gene symbols require an explicit mapping` | Custom panel with repeated symbols | Not supported without changing the engine |
| `No tumour-candidate epithelium …` / `No region meets the declared balance/QC criteria …` | No cluster has `compartment = tumor_candidate`, or no area reaches 15 % tumour / 8 % fibroblast / 15 % immune | Review the annotation; regions are tumour-centred ([06](06_your_own_xenium_data.md#6-known-limitations)) |
| `… is not available. Choose one of [...]` (step 8) | `REGION_TO_EXPORT` misspelled, or `MANUAL` without preview | Use an ID from the table, or click *Preview manual region* first |
| `No usable region fits the total cell limit here` | Manual centre too dense, or near unusable polygons | Move the centre or lower `CORE_TARGET` |
| `Annotations changed or AUCell has not been calculated.` | Step 4 was rerun after step 5 | Rerun steps 5–8 |
| `Preserve the existing export; choose a new destination` | Same destination twice | Rerun step 8 (it adds a new timestamp) |
| Colab crashes with *"Your session crashed after using all available RAM"* | Very large section or other notebooks in the same runtime | Close other sessions, or use a High-RAM runtime |

## Part 2

| Message | Cause | Fix |
|---|---|---|
| `No simulation_input_*.zip under …` | `BUNDLE_ZIP` folder has no bundle | Point to the Part 1 run folder or the ZIP |
| `Bundle integrity check failed: …` / `Unlisted files in the input bundle.` | The ZIP was modified, re-zipped, or a file was added (macOS can add `__MACOSX/` when re-compressing) | Use the original ZIP from Part 1; never unzip and re-zip it |
| `The bundle uses a different runtime …` | The bundle was made by a different engine version | Re-export the region with this package's Part 1, or use the matching package version |
| `Unsafe path in ZIP or manifest.` / `Duplicate entry or symbolic link in ZIP.` | Not a genuine bundle | Do not use it |
| `Atlas fingerprint mismatch` / `Runtime contract mismatch` / `… ID order mismatch` | Edited bundle content | Re-export with Part 1 |

## Part 3

| Message | Cause | Fix |
|---|---|---|
| `No API key provided. No Jev request was sent.` | No secret/env var, or notebook access off | Add `TYPESAFE_API_KEY` in Secrets and enable *Notebook access* ([07](07_jev_api_key.md)) |
| `provider_failed` + `HTTP 401` / `HTTP 403` | Invalid or revoked key | Check the key in the TypeSafe console |
| `HTTP 400` | Request rejected (e.g. model name not available to your account) | Try the other `JEV_MODEL`; check the API docs |
| `HTTP 429` / `HTTP 5xx` / `Transport or response decoding failed` | Rate limit, service or network problem. **No automatic retry**; the failed batch is recorded in the journal | The same `RUN_NAME` will refuse to resend that batch (next row). Wait, then start a new `RUN_NAME`. The new run starts from t0 without the old answers, so run a short pilot with a small `MAX_HTTP_ATTEMPTS` first |
| `This journal contains a failed or unresolved attempt for a required batch. It will not be sent again automatically.` | A previous attempt for the same batch failed or its billing status is unknown | Start a **new `RUN_NAME`** (the old journal stays for audit). Raising the cap does not help |
| `Returned model differs from the requested pinned version.` | The API answered with another model | Choose that model explicitly (new `RUN_NAME`) |
| `The allowance cannot cover the first round …` | `MAX_HTTP_ATTEMPTS` < first-round batches | Raise it above the preflight `initial_http_batches` |
| `request_cap_before_round` (stop reason) | Allowance used up | Resume with the same `RUN_NAME` and a larger `MAX_HTTP_ATTEMPTS` |
| `This run is locked by another or interrupted session.` | A previous session died while running | Make sure nothing else is running it → unlock cell → rerun |
| `This run folder belongs to another configuration …` | Settings, bundle, code or package versions changed | Use a new `RUN_NAME` (or restore the exact settings) |
| `The request limit is below the number already recorded for this run.` | You lowered `MAX_HTTP_ATTEMPTS` below what was already sent | Set it ≥ the previous `http_attempts` |
| `Horizon must be a multiple of 10 minutes, between 10 minutes and 7 days.` | e.g. `HORIZON_HOURS = 0.1` | Use multiples of 1/6 h; max 168 |
| `Population ceiling must cover the starting tissue and be at most 5000 bodies.` | `MAX_POPULATION` < starting cells | Use ≥ the bundle's total cells |
| `Restored this run from …` | Not an error: a Drive copy was used | — |
| `Backup to … failed this time` | Drive copy could not be refreshed | The local run continues. Check Drive access and free space, then rerun the final copy/download cell after the run stops. A runtime loss can lose all responses since the last successful copy |
| A new runtime makes requests that were already paid for | Those responses were absent from the latest successful Drive backup | Only journaled responses in the recovered copy can be reused. Keep hourly backups on and preserve the final run folder/ZIP; see [recovery limits](05_google_colab.md#7-when-colab-disconnects) |

## Part 4

| Symptom | Fix |
|---|---|
| The replay iframe is blank in Colab | Run the cell again, or download `replay.html` and open it in a browser |
| `replay.html` does not play from Google Drive's preview | Download it first; Drive's preview does not run scripts |
| `Expected exactly one run_report.json` | Point `RESULTS_PATH` to one results ZIP or one run folder |
| Video recording does nothing | Use Chrome or Edge (WebM `MediaRecorder`), and keep the tab visible while recording |

Still stuck? Open an issue with the full error message, the notebook cell, `run_report.json` (if any) and the output of the first cell (Python and package versions). **Never paste your key.**
