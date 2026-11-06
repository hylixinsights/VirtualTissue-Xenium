# Changelog

## Unreleased

- Restore the original cell-body renderer as well as the dashboard: bounded weighted silhouettes, world-space sizes, nuclei, tissue context opacity and event glows. Original cells blend supplied t0 contours into the legacy simulated silhouettes; Provided outlines retains exact saved contours. Centroids remains available for inspection. These presentation modes do not alter positions, geometry evidence or simulation mechanics.
- Restore the original dark replay layout and its tissue indicators, chart, event feed and cell inspector while retaining the v3 bounded frame loader. Simulation mechanisms and rates remain unchanged.
- Add `vtx2 replay RUN_DIRECTORY [--output PATH]` to refresh a saved trajectory's presentation without rerunning the simulation or calling Jev. Scientific outputs remain unchanged.

## 3.0.0a4 — 2026-11-02 (alpha prerelease)

- Publish the accumulated v3 implementation for external software testing. The a1–a3 entries below describe internal development milestones included in this prerelease.
- Document the fresh-clone test path, required interpreters and draft, uncalibrated kinetics at the top of the README. Prepare new inputs/runs; changed fingerprints do not permit resuming older runs with this package.

- Provide a single acceptance-test entry point, with an English result page, ROI catalogue and replay. Test failures and unexpected skips prevent creation of a GitHub handoff.
- Translate the ROI interface, scientific scope and Mac launchers into English; rename the launcher files accordingly. Scientific assumptions and numerical rates remain unchanged.
- Bind the tested source, package checks and GitHub patch to an exact source snapshot and upstream base. This prepares a reviewable update without publishing or changing the real Git index.

## Internal milestone — 3.0.0a3 (2026-11-02)

- Incorporate the upstream November 2 kinetic audit into the schema-2.0 table without changing its 46 numerical hazards. Reconcile the active schema-3.0 table with the implemented timed states and action semantics; numerical rates and beta values remain assumptions.
- Distinguish eligible event hazards from cell-cycle duration, tissue turnover, net growth, secretion flux, motility speed and successful event counts. Record cross-lineage and cross-organism extrapolations rather than presenting related citations as calibrations.
- Add portable Mac launchers for installation, the offline ROI catalogue, explicit ROI preparation, expression auditing, local simulation and separately requested Jev runs. Add an object export path using existing labels and polygons before fitting scores.
- Prepared inputs bind the kinetic metadata hash. Rebuild inputs with this version; retain prior packages with earlier runs. Historical a1/a2 validation files remain unchanged.

## Internal milestone — 3.0.0a2 (2026-11-01)

- Accept existing user-supplied labels, coordinates and segmentation polygons with declared provenance. H&E and independent annotation review are optional; the revised three-claim protocol is in `docs/v3/SCIENTIFIC_SCOPE.md` and supersedes those prerequisites in the historical plan.
- Add an offline ROI catalogue showing composition, QC, expression support and geometry; the user selects the region and exports a source-bound selection. Catalogue coverage from an existing preparation is descriptive; scores must be rebuilt for the selected analysis input.
- Preserve supplied polygon vertices and their source hashes, expose invalid-contour/centroid fallbacks, and support explicit boundary-distance interactions. Runtime shapes undergo rigid translation; daughter outlines are declared synthetic inherited copies, not future segmentations or deformable mechanics. Centroid-based spatial readouts remain distinct from the boundary-distance action graph.
- Add opt-in `tiered_panel_auc`: sample/type/depth references, a consistent within-sample type fallback and a bounded weak panel-rank heuristic for rare types. The panel tier uses no invented reference cells and does not preserve the type's mean hazard; unsupported evidence remains neutral. Existing strict scoring and schema-2.0 behavior remain available.
- Add `expression-audit` for fixed-state local probability contrasts against beta-zero and expression permutations, with per-type coverage, action/hazard attribution and explicit possible Jev input routes. Changed probabilities do not guarantee a changed sampled action, trajectory effect or biological validity. Beta zero retains direct expression evidence in Jev payloads.
- Expose evidence status and reference basis in the v3 Jev prompt so empirical cell percentiles and weak panel evidence are distinguishable. No live Jev response effect is established by payload construction or local audits; changed prompt/preparation identities require new compatible runs and journals.
- Freeze practical evaluation procedures for expression effects, context-sensitive Jev preferences and held-out structural compatibility, including ablations, parameter sensitivity and prespecified relevance criteria. These are experiments to perform, not positive findings. Synthetic tests cover the new functionality; real-data claim evaluations and live Jev experiments for this version have not been performed.

## Internal milestone — 3.0.0a1 (2026-10-31)

- Add schema 3.0 measurement, platform, ontology and scoring contracts while retaining schema 2.0.
- Version 22 lineage definitions, hazards/beta, tissue/platform profiles, sourced human/mouse signatures and explicit MGI ortholog pairs. Numerical rates remain assumed; review is incomplete.
- Isolate expression references by sample; default undersized strata to null; make optional no-depth fallback a consistent partition and report its coverage.
- Add timed states, persistent movement, local crowding and daughter repulsion, corpse lifetime and abstract suppressive signals; stream simulation outputs and load offline replay frames from bounded sidecar chunks.
- Add an explicit fixed-bin occupancy kernel and purity gates with distinct readout units.
- Add ensemble, paired readouts, bootstrap summaries, valid Sobol/LHS tooling, structural statistics, ROI candidate selection, reviewer agreement, journal audit and regenerable figures.
- Add native Atera CSR Zarr, Visium HD binned and Stereo GEM adapters and annotation templates with explicit coordinate/feature-universe contracts and format-fixture tests.
- Fix hidden-prior execution to force alpha=1 and record effective settings. Compress v3 evidence payloads without converting missing scores to zero.
- Fix double-counted KILL/APOPTOSIS deaths in general v2/v3 kernels. Kernel fingerprints change; existing runs require the old checkout for resumption. Physical-engine files are unchanged.
- Enforce H5AD cell bounds before backed loading can materialize layers; reject incompatible measurement semantics and changed runtime ontology hashes.
- Compare recorded ontology changes only with an explicit forensic override, sum reported token usage, and preserve lattice units. Add paired standardized effects with explicit undefined-variance handling.
- Remove the unused checkpoint file from v3 output; resume still reconstructs from t0/journals. V2 retains its diagnostic snapshot.
- Add CI coverage/clock artifact and paper figure workflow. Real reviewed annotations, paid Jev experiments and biological/cross-platform validation remain pending; see docs/v3/HANDOFF.md and LIMITATIONS.md.


## Historical baseline — v3 phase 0 (`aa3fe9e`, branch `v3/phase0-hygiene`)

* Kernel: `KILL` and `CLEAR_CORPSE` act on the nearest neighbor in spatial mode and on a seeded draw in well-mixed mode; the previous cell-ID order was a hidden bias. Events and config record the `target_rule`. Existing run fingerprints change.
* Provider: `--jev-prior hidden` omits the guidebook prior from the prompt under its own prompt version and journal identity (anchoring control). Decisions record the Jensen–Shannon divergence between Jev and the gated prior; reports aggregate it per lineage.
* `vtx2 cost` estimates HTTP calls and input tokens before any paid request (calibrated on the recorded live test) and reports how many steps a `--max-http-attempts` budget covers. README examples now use a budget that completes.
* `vtx2 control` writes null-model bundles: expression shuffled within type or globally, label permutation, coordinate shuffle, modifier removal, well-mixed equivalent. `vtx2 compare` matches a control arm to its source run and reports the control.
* Docs: `checkpoint.json` is described as the last committed state for inspection; `--resume` reconstructs from t0 using the journal cache and does not read it. The brain example declares `organism: mouse` and ships an illustrative, unreviewed mouse gene-set file; gene-set files may carry `_`-prefixed note keys.
* Kinetics (lane A): `context.kinetics.mode = hazard_per_hour` interprets action weights as per-hour hazards from the packaged `kinetics_v3.json` (draft anchors with confidence flags; verify before publication). Each step uses a first-event competing-risk discretisation (`P(any) = 1 - exp(-Σk·dt)`, shares `k/Σk`, REST = `exp(-Σk·dt)`); a gated action that leaves the menu returns its share to REST. Steps whose burden exceeds `max_step_probability` (default 0.25) are refused at preparation for always-available actions and stop a run with `kinetics_step_too_coarse` otherwise; preflight reports implied doubling times, the gated worst case and the largest admissible step. `legacy_per_decision` (default for existing contexts) reproduces 2.0.0a1 exactly. Examples now use hazard mode with 0.5 h steps. Time-step invariance is tested (`tests/test_v2_kinetics.py`): rates are under-counted by about burden/2, which the guard bounds.
* Scoring (lane B): `context.scoring.method = depth_conditioned_auc` (`v2/scoring.py`) computes a deterministic recovery AUC over the panel, a mid-rank percentile within (cell type × depth bin) with documented fallbacks, `z = probit(p)`, and multiplies hazards by `exp(β·z)/mean_ref(exp(β·z))`, preserving the reference-group mean. `β` is global or per action (|β| ≤ 3). `min_transcripts` makes low-count cells REST-only and excludes them from references. Inspect reports signature coverage; `vtx2 control --kind beta_zero` is the label-only control. `rank_legacy` (default for old contexts) is unchanged. Examples use the new method.
* At this historical baseline, timed states and motility, data-driven ontology, signature curation, payload compression, platform adapters and real-data runs were still pending. The current status is in the implementation entry above and `docs/v3/HANDOFF.md`.

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
