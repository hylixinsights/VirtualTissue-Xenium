# Part 4 · Results and the replay viewer

**Notebooks:** section **PART 4** of `Part2-4_Run_Jev_Simulation.ipynb` (right after a run), and
`Part4_Explore_Saved_Results.ipynb` (any time later, without the API)
**Command line:** `vtx explore RESULTS.zip --plot overview.png --open`
**Python:** `virtualtissue_xenium.results`

---

## 1. Where the results are

| Location | Content |
|---|---|
| `/content/Xenium_Jev/results/<RUN_NAME>/` | The run folder on the Colab session disk (temporary) |
| `/content/Xenium_Jev/results/<RUN_NAME>_results.zip` | Everything in the run folder, zipped |
| `MyDrive/VirtualTissue/Simulation_Results/<RUN_NAME>/` and `…_results.zip` | Copy in your Drive when `SAVE_RESULTS_TO_DRIVE` is on, refreshed every simulated hour and at the end |
| Your computer | The ZIP, if `DOWNLOAD_RESULTS` is checked in the last cell |
| `/content/Xenium_Jev/results/preflight/<RUN_NAME>/` | `preflight.json`, `first_request_preview.json`, `measured_t0.png` |

Locally (`vtx run` or Jupyter), the same structure lives under the results folder you chose (`vtx_session/results/` in the notebooks).

## 2. Files in a run folder

| File | What it is |
|---|---|
| **`replay.html`** | Self-contained interactive viewer (≈ 0.3–10 MB). Opens offline in Chrome, Firefox, Edge or Safari, with no key and no network |
| `run_report.json` | `completed`, `stop_reason`, `error`, `provider_error`, `checkpoint_time_min`, `usage` (`http_attempts`, `reused_batches`, `input_tokens`, `output_tokens`), episode `summary` |
| `episode.json.gz` | The complete trajectory: `frames` (every 20 simulated min: ids, x, y, size, state of each body, fields), `bodies` (id, label, origin `roi`/`halo`/`occluded`/`born`, parent, birth time, fate), `decisions`, `events`, `metrics`, `summary`, fingerprints |
| `decisions.csv` | One row per cell question: `time_min`, `cell_id`, `provider`, `wake_reasons`, `observation` (full local record; recognition and exclusion diagnostics are omitted from the Jev payload), `options` (the menu, with gates), `answer` (choice, probabilities, confidence), `status` |
| `events.csv` | `time_min`, `kind`, `cell_id`, … for `division`, `cycle_entered`, `cytotoxic_pulse`, `apoptosis_committed`, `intact_death`, `clearance`, `corpse_degraded`, `program_started/completed`, `target_reached`, `division_blocked_population_ceiling`, `division_blocked_local_pressure`, `division_blocked_compression`, `growth_blocked_compression`, `exhaustion_transition`, `hypoxic_necrosis`, … |
| `metrics.csv` | Every 10 simulated min: `tumour_cells`, `tumour_area_um2`, `kills`, `born`, `decisions`, `cd8_tumour_contacts`, `cd8_within_20um_of_tumour`, `oxygen_tumour_mean`, `hypoxic_tumour_fraction`, `pressure_mean/max`, `pressure_by_label`, … (dict columns are JSON) |
| `jev_journal.jsonl` | Hash-chained log of every request intent and result. Used to resume. The key is never present |
| `effective_contract.json` | The contract actually used (scenario, horizon, population ceiling, model) |
| `run_config.json` | Configuration fingerprint: bundle SHA-256, settings, engine source hashes, Python and package versions |
| `input_*.json` | Provenance copied from the bundle (`SOURCE`, `aucell_provenance`, `BUNDLE_MANIFEST`, `region_selection_rules`, `score_calibration`) |
| `OUTPUT_MANIFEST.json` | SHA-256 of every output file |
| `overview.png` | Written by the notebook's summary cell |

The input bundle is **not** copied into the results. Keep it to reproduce the run.

---

## 3. The replay viewer

![Replay viewer](images/replay_R03.png)

| Control | What it does |
|---|---|
| **Play / pause** (or the space bar), **restart**, the **time slider**, **speed** (0.5–4 simulated h per second) | Navigate in time. The view starts from the **measured** segmentation and cross-fades into the simulated cells within the first 20 min |
| **Color** | `Cell type` (default) · `Current action` · `AUCell score at t0` (then pick a signature, e.g. mitotic cell cycle) |
| **Signals** | Buttons for Auto, Off and all ten fields: IFN-γ-like, CXCL9/10-like, oxygen, TGF-β-like, glucose, lactate, VEGF-like, matrix, NE-like and CGRP-like. Auto stays on oxygen throughout the run, including zero. Expand the field summary for meaning and scale. Missing data, disabled production, zero and uniform values are distinct. |
| **Highlights** | Up to three arrow callouts and clickable cards for recorded recent events and field maxima. Removed cells retain event coordinates; legacy positions are marked approximate. A peak is not a statistical or biological anomaly. |
| **Context** | Show or hide the grey context cells |
| **Interactions** | Arrows (crawling toward a target), red flashes (cytotoxic hits), cycle rings (dashed when paused; inspect the recorded cause), bursts (division), grey shrinking corpses, gold flashes (phagocytosis) |
| Mouse wheel / drag / double-click | Zoom, pan, fit to view |
| Hover a cell | Tooltip with type and current state |
| **Click a cell** | **Inspector:** its t0 AUCell percentiles as bars, its parent, and its latest six decisions with the top probabilities. Full menus, exclusions and decision history are in the saved episode and Python tables |
| Legend (cell types) | Click a type to hide or show it |
| Right panel | Live counters (tumour candidates, divisions, deaths, phagocytosed, corpses, AI decisions), a tumour-count chart and a recent-events feed |
| **● Record video** | Records the canvas as a `.webm` file while it plays |

In Colab the replay is shown in an `iframe` inside the notebook. For the smoothest experience, open `replay.html` from the downloaded ZIP or from Drive (*Download* first; Drive's preview does not run it).

---

## 4. Explore saved results later (`Part4_Explore_Saved_Results.ipynb`)

This notebook needs no API key and does not install the simulation stack. It only uses pandas and matplotlib, which Colab already has.

| Cell | ✏️ | Shows |
|---|---|---|
| Setup | `CODE_SOURCE` | — |
| Result to open | `RESULTS_SOURCE`, **`RESULTS_PATH`** (ZIP or folder in Drive) | Headline table |
| Interactive replay | — | `replay.html` inline |
| Overview | — | Four panels (living acting cells per type, cumulative decisions, cumulative events, choices by type); first/last population; choice × type table; events per 6 h |
| One cell | **`CELL_ID`** (0 = the first cell that decided) | Every decision of that cell, with the full probability distribution and its wake reasons |
| Compare runs | **`COMPARE_PATHS`** (comma-separated) | Validated matched settings and outcomes at the latest common saved time; separate achieved/requested horizons, death mechanisms and ceiling warnings |

![Overview figure](images/results_overview_local_preview.jpg)

*`plot_overview` for a 24 h **local preview** of R03 under `declared_recognition` (programmed rules, not Jev). 78 tumour cells entered the cycle at t0, but only 11 could divide: the default population ceiling of 1,000 bodies leaves 11 free slots.*

## 5. Python API

```python
from virtualtissue_xenium import results as R

res = R.load_results("R03_jev_run01_results.zip")     # or the run folder
R.summary(res)                 # dict: mode, scenario, stop_reason, simulated_hours, cell_decisions, http_attempts, tokens, …
R.population_by_type(res)      # DataFrame: living acting cells per type at every frame (index = simulated hours)
R.decision_table(res)          # one row per decision: hours, cell_id, cell_type, choice, target_id, probability, confidence
R.choice_counts(res)           # cell type × chosen action
R.cell_history(res, 11869)     # every decision of one cell with all probabilities
R.event_counts(res, bin_hours=1)
fig = R.plot_overview(res); fig.savefig("overview.png", dpi=150)
R.compare_runs(["runA_results.zip", "runB_results.zip", "runC_results.zip"])
```

`res` also gives direct access to `res["episode"]` (the full trajectory dict), `res["metrics"]`, `res["events"]`, `res["decisions"]` (pandas DataFrames), `res["report"]` and `res["replay"]` (path).

**Example: how did tumour cells with a high mitotic score decide?**

```python
import pandas as pd
d = R.decision_table(res)
tumour = d[d.cell_type == "Tumor"].copy()
# t0 percentile of each asked cell, taken from what was actually sent (missing evidence stays missing)
mitotic = {x["cell_id"]: (x["observation"].get("t0_aucell_percentile") or {}).get("MITOTIC_CELL_CYCLE")
           for x in res["episode"]["decisions"]}
tumour["mitotic_pct"] = tumour.cell_id.map(mitotic)
tumour.groupby(pd.cut(tumour.mitotic_pct, [0, .25, .5, .75, 1]), observed=False).choice.value_counts(normalize=True).unstack()
```

## 6. Interpreting results

* Probabilities are Jev's **preferences among the options offered at that moment**. They are not rates, and they change when the menu changes. Always read them together with `options` and `excluded`.
* A count such as "tumour cells at 24 h" depends on declared assumptions: the accelerated 24 h cycle, damage per hit, the population ceiling and the scenario. Compare runs that differ in **one** setting (e.g. `declared_recognition` vs `recognition_off`, same seed) rather than reading absolute numbers.
* Local-preview results are labelled `local…` in `provider` and in the folder name. Never mix them with Jev runs without saying so.

## Reading highlights and fields

Arrow labels and clickable cards show up to three recent recorded events (within 40 simulated minutes) and field maxima at the current saved frame. A maximum is not a secretion location or a biological anomaly. Field colours are approximate (8-bit display encoding); 0.1.2 also saves each frame's unquantized minimum and maximum for the legend. The scale is fixed across each run and may differ between runs. Use common scales for quantitative comparisons. Zero-field frames are explicitly labelled; inactive fields do not display the internal normalization floor as a signal scale. TNF-alpha is not implemented.

The scenario banner explains when CD8 killing is disabled. The cumulative death tile counts apoptosis commitments plus hypoxic necrosis, without recounting completion of apoptosis. Playback controls display existing results; they do not launch an intervention or call Jev.

## 7. Matched comparisons in 0.1.2

```bash
vtx compare run_declared_results.zip run_off_results.zip --shared-replays comparison_replays
```

```python
paths = ["run_declared_results.zip", "run_off_results.zip"]
table = R.compare_runs(paths)
replays = R.export_comparison_replays(paths, "comparison_replays")
```

The supported pair is `declared_recognition` versus `recognition_off`; identical-scenario sham repeats are also accepted. Bundle SHA, atlas, seed, provider/model, software, environment, decision limit and contract must match except the recognition switch. Different seeds require a separate replicate analysis. Incomplete provenance or incompatible settings raise a clear error.

Outcomes are evaluated at the latest timestamp saved in both runs. A 9.5 h endpoint and an 18.67 h endpoint are never compared as equivalent. The common saved frame may be earlier than the shorter endpoint because frames are discrete. `simulated_hours` and `comparison_time_hours` identify the comparison time; `reached_hours`, `requested_hours`, completion status and HTTP usage still describe the original run. No source recording is truncated.

Replay copies use the same scales and comparison endpoint, but their controls are not synchronized. Fields are recovered from the full episode. Opening an old standalone HTML cannot recover fields that were never exported into it.

`summary()` reports apoptosis, hypoxic necrosis, unique deaths, division-block events by cause, growth-compression blocks and the time the cumulative population ceiling was reached. Counts of block events are not counts of missed divisions. Legacy crowding events below the pressure threshold can be attributed to the ceiling; ambiguous legacy events remain explicitly unspecified. Legacy Jev comparisons warn that old prompts exposed recognition flags.
