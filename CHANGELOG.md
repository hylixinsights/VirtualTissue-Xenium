# Changelog

## Unreleased (v3 phase 0, branch `v3/phase0-hygiene`)

* Kernel: `KILL` and `CLEAR_CORPSE` act on the nearest neighbor in spatial mode and on a seeded draw in well-mixed mode; the previous cell-ID order was a hidden bias. Events and config record the `target_rule`. Existing run fingerprints change.
* Provider: `--jev-prior hidden` omits the guidebook prior from the prompt under its own prompt version and journal identity (anchoring control). Decisions record the Jensen–Shannon divergence between Jev and the gated prior; reports aggregate it per lineage.
* `vtx2 cost` estimates HTTP calls and input tokens before any paid request (calibrated on the recorded live test) and reports how many steps a `--max-http-attempts` budget covers. README examples now use a budget that completes.
* `vtx2 control` writes null-model bundles: expression shuffled within type or globally, label permutation, coordinate shuffle, modifier removal, well-mixed equivalent. `vtx2 compare` matches a control arm to its source run and reports the control.
* Docs: `checkpoint.json` is described as the last committed state for inspection; `--resume` reconstructs from t0 using the journal cache and does not read it. The brain example declares `organism: mouse` and ships an illustrative, unreviewed mouse gene-set file; gene-set files may carry `_`-prefixed note keys.
* Kinetics (lane A): `context.kinetics.mode = hazard_per_hour` interprets action weights as per-hour hazards from the packaged `kinetics_v3.json` (draft anchors with confidence flags; verify before publication). Each step uses a first-event competing-risk discretisation (`P(any) = 1 - exp(-Σk·dt)`, shares `k/Σk`, REST = `exp(-Σk·dt)`); a gated action that leaves the menu returns its share to REST. Steps whose burden exceeds `max_step_probability` (default 0.25) are refused at preparation for always-available actions and stop a run with `kinetics_step_too_coarse` otherwise; preflight reports implied doubling times, the gated worst case and the largest admissible step. `legacy_per_decision` (default for existing contexts) reproduces 2.0.0a1 exactly. Examples now use hazard mode with 0.5 h steps. Time-step invariance is tested (`tests/test_v2_kinetics.py`): rates are under-counted by about burden/2, which the guard bounds.
* Not yet done (see the v3 plan): depth-conditioned scoring and the β effect size, timed states and motility, data-driven ontology, platform adapters, real-data runs.

## 2.0.0a1 — 2026-10-30

* Introduce the separate `vtx2` research workflow for annotated AnnData, JSON and Xenium input, explicit tissue/disease/perturbation/sampling context and panel-aware expression guidebooks.
* Add lineage-restricted actions, per-cell seeded categorical sampling, concurrent bounded Jev requests, synchronous commit barriers, abstract local signals, alive-capacity accounting and auditable state transitions. Preserve the original physical engine and `vtx` interface.
* Persist immutable API intents/results, configuration-bound cache, last-committed-state snapshots (not consumed by resume), run locks, per-cell probabilities, provenance, offline replay and matched shared-endpoint comparisons.
* Package verified skin melanoma and WT/Tg brain source manifests, bounded download/ROI/annotation adapters, executable examples, a notebook, code review and an ordered parallel/sequential development plan.
* Keep the JSON core dependency-free; scientific input support is the `xenium` extra and the unchanged pinned legacy stack is the `legacy` extra.
* Record real-data compatibility checks and a live Jev protocol test. Preferences, timescales and state transitions remain uncalibrated assumptions; biological predictive validation and learned transition models remain future work.

## 0.1.3 — 2026-10-28

* Jev prompt `vtx-choice-3` makes tissue context a user choice: none (default), custom text (`--context`, alias `--tissue`), or explicitly selected bundle metadata (`--context-mode bundle`). The notebook exposes `CONTEXT_MODE` and `TISSUE`. Concise t0 semantics and recognition-arm isolation are preserved, and the resolved context is recorded. The effect on tumour proliferation is awaiting live Jev validation; see [validation and usage](docs/12_jev_feedback_013.md).
* Accept a reported choice within one 0.01 rounding step before normalization. Preserve exact reported probabilities and repair metadata while giving the unchanged bundled kernel a normalized tie. Rejected serializable answers are journaled for audit; no automatic retry is added.
* Preflight reports AUCell coverage by acting type, CD8 access to tumour and t0 eligibility separately from the decision cap. Results report cycle probabilities and decision margins using original reported values.
* Add prompt snapshots and a bounded, dry-by-default ablation tool with fixed groups, durable request intents and refusal to resend an existing experiment. Existing input bundles remain compatible; use new run names with this version.
* Advance the simulation workflow marker so Colab requests a restart when 0.1.2 or development-candidate simulation modules are already loaded during an upgrade.

## 0.1.2 — 2026-10-26

* Results include hypoxic necrosis, unique deaths and division blocks by cause. `vtx compare` validates matching provenance/settings and compares the latest common saved frame. Original runs remain intact; optional replay copies share field scales and the comparison endpoint.
* Replay exports all ten fields, keeps Auto on oxygen, renders zero oxygen as anoxia, distinguishes disabled/missing/zero/uniform fields, and labels field maxima accurately. Field details collapse, callouts clamp to map edges, necrosis has a saved position, and population-ceiling warnings no longer claim contact inhibition.
* A reporting adapter records event coordinates and distinguishes population ceiling, local pressure and compression. The bundled biological runtime, its population-limit rule and its trajectories are unchanged. Existing input bundles retain their integrity checks and remain compatible.
* Jev prompt `vtx-choice-2` omits global recognition flags and gate-exclusion diagnostics; the latter remain in the episode. CD8 questions are batched separately from other cell types. Frozen t0 AUCell semantics are explicit. Jev results can change and require new run names; live Jev behavior was not validated for this release.
* Region feasibility is checked before whole-section AUCell. Annotation validation warns about tumour-free maps; preflight warns about few active CD8 cells and limited growth headroom. Documentation explains finite pulses, effective damage, synchronized accelerated cycles and the cumulative computational ceiling.
* Notebooks detect the current workflow version correctly on repeated setup. Resume 0.1.1 runs with the 0.1.1 checkout; do not overwrite saved fingerprints. Persistent AUCell caching and scheduled interventions remain future work.
* Validation details and limits are recorded in `tests/reports/release_0.1.2.json`. Offline JavaScript tests cover replay logic at 1400×900 and 1920×1080; browser visual inspection and live Colab/Jev are not claimed.

## 0.1.1 — 2026-10-26

* Quickstart explicitly selects declared recognition for immune-killing experiments and explains the matched recognition-off control, CD8 composition and growth ceiling. CLI and notebook preflight explain the selected scenario.
* Replay adds automatic signal selection, visible field buttons, quantitative SIM legends, zero-field messages, arrow callouts and clickable highlights of recorded events/local signal peaks. Cumulative deaths include hypoxic necrosis without double-counting completed apoptosis. Legacy flat-ring LUAD replays are readable.
* Non-tumour annotations stop preparation before AUCell and geometry with a clear CLI message, retaining the annotation audit.
* Experiment design documents the proposed paired perturbation workbench and the missing models/adapters for drugs, resistance, neutrophils, TNF, healthy tissue and other platforms. These capabilities are not advertised as implemented.
* Biological runtime and input-bundle compatibility are unchanged. Replay-source changes alter run fingerprints: resume earlier runs with their original checkout; new runs use a new name.

## 0.1.0 — 2026-10-25

First packaged release.

* The validated engine is packaged **byte-for-byte unchanged** in `src/virtualtissue_xenium/engine/`, with its 32 files hashed in `ORIGINAL_SOURCES.json`:
  * preparation `lung5k-preparation-1.1.0` (from *Xenium_Lung5K_Prepare_AUCell* v1.1);
  * simulation `xenium-jev-colab-1.1.0` (from *Xenium_Lung5K_Run_Jev_NoDrive* v1.1).
* Three notebooks, generated by `tools/build_notebooks.py`:
  * Part 1 gets the code from GitHub or Drive, can download the public example, and exports the region set in an explicit `REGION_TO_EXPORT` parameter;
  * Parts 2–4 take a bundle from Drive, upload or URL, keep a resumable Drive copy of each run, and show summary figures;
  * a new explorer notebook opens saved results.
* New helpers: `download` (HTTP byte-range download of the public example with SHA-256 verification), `own_data` (data folder and manifest for your own Xenium run, cluster-map template with marker-detection summary, validation, activation), `colab` (key lookup, Drive backup/restore), `results` (Part 4 analysis), and the `vtx` command line.
* Documentation: README and ten guides in `docs/`.
* First-user walkthroughs: explicit public/own-data routes; complete Colab and local file preparation with consistent paths, Drive mounting, a filled cluster-map example, validation and simulation handoff; local notebook path instructions; API-key access enabled in the simulation notebook.
* Recovery documentation states that only responses preserved in the available journal are reused; paid responses after the last successful Drive backup can be lost with the runtime.
* Validation: the package reproduces the earlier R03 bundle byte for byte; the original engine tests pass; all notebooks were executed end to end with a mocked Colab and a mocked Jev, including a restore from Drive that sent no new requests.
