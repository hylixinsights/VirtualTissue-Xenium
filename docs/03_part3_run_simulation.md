# Part 3 · Run the simulation

**Notebook:** `Part2-4_Run_Jev_Simulation.ipynb`, section **PART 3**
**Command line:** `vtx run BUNDLE.zip RESULTS_DIR --run-name NAME --mode jev …`
**Engine functions:** `simulation_workflow.effective_contract`, `preflight`, `run`; `jev_provider.JevProvider`

---

## 1. Experiment settings

| Setting | Default | What it controls | Notes |
|---|---|---|---|
| `RUN_NAME` | `R03_jev_run01` | Name of the results folder `<RUN_NAME>/` and `<RUN_NAME>_results.zip` | Letters, digits, `_`, `-` (≤ 80). **Keep it to resume; change it for a new experiment** |
| `RUN_MODE` | `Jev` | `Jev` = paid Choice API calls · `Local preview (no API)` = programmed rules | Local preview results get the suffix `_local_preview` |
| `SCENARIO` | `baseline_unknown` | Declared assumptions about recognition and effector competence (§2) | |
| `JEV_MODEL` | `jev-1.14.0` | Pinned model version; the run stops if the API returns another model | `jev-1.13.0` was used by the LUAD3 prototype |
| `CONTEXT_MODE` | `none` | `none`: no additional context; `custom`: your text; `bundle`: tissue from `SOURCE.json` | CLI: `--context-mode`. Recorded in the run fingerprint |
| `TISSUE` | blank | Your tissue/context description when `CONTEXT_MODE = "custom"` | CLI: `--context` (alias `--tissue`); supplying text alone selects custom. Use the same context in paired runs |
| `HORIZON_HOURS` | 24 | Simulated time | Multiple of 10 min; max 168 h (7 days) |
| `SEED` | 21 | Seeds tie-breaking, motility, ordering | Same seed + same answers → same trajectory |
| `MAX_POPULATION` | 1000 | Ceiling on allocated bodies (dead/removed bodies still count) | Must be ≥ the starting cells and ≤ 5000. R03 starts with 989, so only 11 births are possible by default |
| `MAX_HTTP_ATTEMPTS` | 600 | **Hard cap on POST requests to Jev** across all resumptions of this run, failures included | Your cost control; not a dollar cap |
| `MAX_CELL_DECISIONS` | 10000 | Cap on individual cell questions | If a round exceeds what is left, a seeded subset is asked and the run stops after it |

The cell checks the values and builds the *effective contract*: the bundle's contract with your scenario, horizon, population ceiling and model applied. The bundle itself is never modified.

Context is optional. The default does not read tissue metadata into the prompt.
For custom context, enter a single line of up to 300 characters; for bundle
context, `SOURCE.json` must contain a tissue description. Preflight shows the exact
resolved prompt. These choices affect what Jev sees, not the measured AUCell scores.

## 2. Scenarios

**For the immune-killing walkthrough, explicitly select `declared_recognition`.** The default remains `baseline_unknown` so that antigen recognition is never silently assumed. A run with that default can grow or move but cannot produce CD8 cytotoxic attacks simply by waiting longer. Necrosis from sufficiently sustained hypoxia is a separate mechanism. Compare the same bundle under `declared_recognition` and `recognition_off`; inspect active CD8 counts and growth headroom in the preflight. See the [experiment design guide](11_experiment_design.md) for the distinction between a measured starting state, a control and an intervention.

| Scenario | Recognition | Effector / clearance competence | Use it for |
|---|---|---|---|
| `baseline_unknown` | unknown (`null`) | unknown | What follows from measured state alone. Cytotoxic killing and efferocytosis are **not licensed**, because unknown never opens a gate |
| `declared_recognition` | `true` | `true` | "What if these CD8 T cells recognise the tumour and these macrophages can clear corpses?" Licenses `DEGRANULATE`, `EFFEROCYTOSE`, `HELP_CD8` when their other gates pass. It does not force killing |
| `recognition_off` | `false` | `true` | Negative control for `declared_recognition` |

All three keep motility, fields and the tumour cycle on; angiogenesis, recruitment, nerves and checkpoint blockade stay off. Recognition is **never inferred from RNA**. It is always a declared assumption, recorded in `effective_contract.json`.

## 3. Preflight: zero API calls

The *Preflight* cell:

* runs one **static** mechanical step and checks that nothing moves (`static_preflight_passed`);
* builds the **first decision round** exactly as the run would: how many cells will be asked (`initial_cell_questions`) and in how many HTTP requests (`initial_http_batches`);
* saves `first_request_preview.json`, the literal JSON body of the first request, so you can read what Jev will receive;
* reports AUCell coverage per acting type, CD8–tumour contacts/proximity, eligible questions and those deferred by the decision cap;
* warns if `MAX_HTTP_ATTEMPTS` cannot cover even the first round. In that case the run would stop before sending anything.

For the tested R03 bundle with its declared lung tissue, the `vtx-choice-3` first round is **280 questions in 25 requests**. Longer tissue descriptions or different inputs can change batching; use your preflight count. Later rounds depend on what happens: contacts form, programs finish, cycle checkpoints arrive. The first-round count is therefore not a forecast of the total. An earlier 24 h validation run on R03 with mocked answers that followed the local rules needed 122 requests. Real Jev choices can need more or fewer. Start with a small cap, look at the results, then resume with a larger cap (§6).

## 4. The Jev key and START_RUN

1. Put the key in Colab **Secrets** as `TYPESAFE_API_KEY` (or `JEV_API_KEY`) and enable **Notebook access**. Locally, use an environment variable with the same name. If no key is found, a hidden prompt appears. See [07_jev_api_key.md](07_jev_api_key.md).
2. Check **`START_RUN`** and run the cell. Checking it authorises up to `MAX_HTTP_ATTEMPTS` requests. With `START_RUN` unchecked, the cell only prints a reminder.

---

## 5. What happens during a run

```text
every 10 simulated minutes (one tick):
  kernel: fields, mechanics, programs, clocks, damage, division, clearance
  scheduler: which acting cells have a wake reason and are past their minimum interval?
  for those cells: build the menu of permitted actions (gates evaluated; failures listed in "excluded")
  cells with ≥ 2 options → one round of questions
      → batches of ≤ 20 questions / ≤ 40 KB
      → each batch: journal "intent" → POST → validate → journal "result"
      → the whole round validated → choices applied in seeded order
  every 6 ticks (1 simulated hour): checkpoint episode.json.gz (+ Drive copy in Colab)
```

**Wake reasons:** a new relevant heterotypic contact (initial contacts count as new), a field crossing a threshold bin, a finished program, a cell-cycle checkpoint, entry into the region. Quiet cells are not asked and cost nothing.

**Minimum interval between questions:** CD8 60 min (30 min while touching tumour), CD4 and macrophages 180 min, endothelium 360 min, tumour, CAF and fibroblasts 720 min; at most 24 questions per cell per simulated day. B and plasma cells have no menu and are never asked.

### What Jev receives

One request = `{"model", "state", "questions"}`. Each question is a `choice` with:

* `instructions.local_observation` — `time_min`, `primary_label`, `identity`, `t0_aucell_percentile` and `t0_aucell_raw` (the 20 signatures, missing ones absent), `t0_aucell_provenance`, `sim_state` (viability, damage, effector pulses, exhaustion, glycolysis, PD-L1, polarisation), `fields_sim` (oxygen, glucose, lactate, IFN-γ-like, CXCL9-like, VEGF-like, TGF-β-like, matrix), `local_pressure_sim`, `contacts` (id, label, viability, PD-L1, provenance). The global `recognition_sim` and `excluded` diagnostics remain omitted from `vtx-choice-3` provider payloads and retained in local decision records;
* `criteria` — the permitted options and their one-line descriptions. Targeted actions are keyed `ACTION:<target_id>`.

It never contains raw counts, gene lists, coordinates, the contract or a desired outcome. See the full real example in the [README](../README.md#43-how-one-cell-asks-jev-part-3).

### How answers are validated

The returned model must equal `JEV_MODEL`. Question IDs and option keys must match exactly. Probabilities and confidence must be finite and in [0, 1], and the probabilities must sum to 1 within a rounding tolerance. The choice must be a maximum-probability option or trail it by at most one reported rounding step (`0.01 + 1e-9`, before normalization). A compatibility copy ties an accepted near-maximum choice and is normalized for the unchanged kernel; exact reported probabilities and repairs remain in the answer and journal. See the [rounding and prompt notes](12_jev_feedback_013.md). If any answer in a round fails, **no decision of that round is applied** and the run stops with `provider_failed`. There are **no automatic retries or redirects**, and no fallback to local rules.

### Stop reasons (`run_report.json → stop_reason`)

| Stop reason | Meaning | What to do |
|---|---|---|
| `horizon_reached` | Finished `HORIZON_HOURS` | — |
| `request_cap_before_round` | The next complete round would exceed `MAX_HTTP_ATTEMPTS`; it was not sent | Resume with a larger cap (same `RUN_NAME`) |
| `decision_cap_reached` | `MAX_CELL_DECISIONS` used up (the last round may be a seeded subset) | The decision cap is part of the configuration: start a new `RUN_NAME` with a larger value |
| `provider_failed` | HTTP error, timeout, malformed answer or model drift. The message is in `provider_error` | Read the message. Fix the key, model or connectivity. An unresolved batch is **not** re-sent automatically (§6) |
| `error` set, no stop reason | Interrupted (e.g. you stopped the cell) | Rerun the same cell to resume |

In all cases the last complete checkpoint is kept, and the replay and tables are written up to that point.

---

## 6. Resuming with saved answers

Every request is appended to `jev_journal.jsonl` **before** it is sent, and its validated response is appended before it is applied. Entries are hash-chained and bound to a fingerprint of the run configuration (bundle, settings, engine source, Python and package versions, machine type).

To resume, **run the same cells with the same input ZIP, `RUN_NAME` and settings**, in the same software environment. Set `BUNDLE_ZIP` to the exact original ZIP rather than a folder that might now contain a newer export:

* the engine replays the trajectory **from t0**, reusing every matching valid response in the available journal (`reused_batches` in the report), and sends only the missing batches;
* you **may raise `MAX_HTTP_ATTEMPTS`**: the request cap is not part of the configuration, so the journal stays valid. `HORIZON_HOURS`, `SCENARIO`, `SEED`, `JEV_MODEL`, `MAX_CELL_DECISIONS` and `MAX_POPULATION` **are** part of it. Changing any of them requires a new `RUN_NAME`, and the earlier answers cannot be reused;
* if the configuration differs, the engine refuses: *"This run folder belongs to another configuration…"*. Choose a new `RUN_NAME`.

| Situation | What happens |
|---|---|
| Same Colab runtime, files still there | Rerun the run cell → resumes |
| Colab runtime was recycled, `SAVE_RESULTS_TO_DRIVE` on | The latest successful Drive copy is restored. Responses in that copy are reused; paid responses since that backup may be lost and requested again. Before the first successful copy there may be nothing to restore |
| Runtime recycled, no Drive copy | Starts from scratch and **may repeat paid calls**. Keep `SAVE_RESULTS_TO_DRIVE` on, or download results before leaving |
| A request failed (HTTP error, timeout, invalid answer) or its answer never arrived (disconnect mid-request) | The journal holds an *intent* without a valid result. The engine refuses to resend it automatically, because billing is unknown. Inspect the journal and start a new `RUN_NAME` if you accept the risk |
| `.running.lock` left by a killed session | Make sure no other session writes to this run, then use the *unlock* cell (`CONFIRM_PREVIOUS_SESSION_STOPPED`) |
| Different machine (Mac ↔ Colab) | Not resumable (the environment is part of the configuration). Start a new `RUN_NAME` |

**Backup limits:** Drive is refreshed after each completed simulated hour and at the end, not after every request. With successful hourly backups, losing the runtime can lose the responses recorded since the previous hour boundary. Failed backups can extend that gap. The recovered request count also excludes lost entries, so `MAX_HTTP_ATTEMPTS` is not an account-wide spending limit after data loss. A failed/unresolved request is blocked only when its entry survives in the available journal. Saving the notebook does not save the journal; see [Colab recovery](05_google_colab.md#7-when-colab-disconnects).

---

## 7. Local preview (free)

`RUN_MODE = "Local preview (no API)"` (CLI: `--mode local`) replaces Jev with `ScoreAwareLocalPolicy`, a small set of programmed rules that use the same menus and the same t0 AUCell evidence. For example, a tumour cell enters the cycle only if its `MITOTIC_CELL_CYCLE` percentile is ≥ 0.75 and its raw score is > 0. Use it to check the pipeline, the replay and the settings before spending credits. Its results are labelled as local everywhere (`provider = local…`, folder suffix `_local_preview`). It is a test policy, not a biological model.

---

## 8. Command-line equivalent

```bash
vtx check BUNDLE.zip --scenario declared_recognition --hours 24           # preflight only
vtx run BUNDLE.zip outputs/runs --run-name R03_local --mode local --hours 24
TYPESAFE_API_KEY=... vtx run BUNDLE.zip outputs/runs --run-name R03_decl01 --mode jev \
     --scenario declared_recognition --hours 24 --seed 21 --max-http-attempts 300 --max-population 1500
# resume later with a larger cap (same name, same settings):
TYPESAFE_API_KEY=... vtx run BUNDLE.zip outputs/runs --run-name R03_decl01 --mode jev \
     --scenario declared_recognition --hours 24 --seed 21 --max-http-attempts 600 --max-population 1500
```

`--mode jev` asks for confirmation (type `yes`) after printing the first-round estimate. Add `--yes` for unattended runs.

## 9. Python equivalent

```python
import virtualtissue_xenium as vtx
from virtualtissue_xenium.colab import get_jev_key
wf = vtx.workflow
bundle = wf.extract_bundle("simulation_input_R03_….zip", "work/extracted")
settings = dict(model="jev-1.14.0", scenario="baseline_unknown", seed=21, horizon_hours=24,
                max_decisions=10000, max_population=1000)
facts, first_request = wf.preflight(bundle, settings)
result = wf.run(bundle, settings, "work/results/R03_jev01", api_key=get_jev_key(),
                max_http_attempts=600, mode="jev")          # or mode="local_preview", api_key=""
result["report"], result["zip_path"]
```
