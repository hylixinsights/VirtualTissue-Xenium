# VirtualTissue-Xenium 5.0.0

Released 2026-11-10. V5 is the current primary software release on `main`.

V5 keeps the v4 simulation and its replay unchanged. It improves the evidence that reaches the per-cell rulebook, what the Jev model is told about that evidence, and the tools that measure the quality of Jev's decisions. Use `vtx5`; `vtx4` and `vtx2` remain aliases and `vtx` keeps the legacy physical engine.

**Scientific status is unchanged.** Kinetic rates and beta values are assumptions. Passing software checks does not validate a biological prediction. Every evidence rule is a modelling assumption that needs biological review.

## What did not change

- Action menus, permission gates, recognition, the interaction graph.
- Movement, division, death, clearance, signalling and state clocks.
- The kinetic and beta tables, the multiplier law and the competing-risk conversion.
- The replay, the dashboard and the cell renderer (`v2/web/`, `replay.py`).
- The default behaviour of every v4 command.

Existing contexts without an evidence policy or the new optional sample-context fields give byte-identical `frames.jsonl`, `config.json`, `guidebook.json`, prepared and control fingerprints, expression-audit fingerprints and Jev request bodies, and a decision log that expands to the v4 log byte for byte. This is checked against hashes recorded on the v4.0.0 tag:

```bash
python tools/v4_equivalence.py compare tests/fixtures/v4_0_0_reference.json
```

The same comparison runs in the test suite (`tests/test_v5_v4_equivalence.py`).

## What is new

### Sample context collected in English

Use `vtx5 context-init sample_context.json --template examples/v3/context_user_object.json` before preparing your own object. The questionnaire records sample details and can add explicit lineage-specific restrictions with reasons. Free text is descriptive; it never silently changes local rates or creates actions. Mature neurons already lack division permission, whereas a brain sample can contain proliferative glia or malignant cells. `SECRETE` is an abstract inflammatory/suppressive signal; mucus production is unsupported. See [SAMPLE_CONTEXT.md](SAMPLE_CONTEXT.md).

### 1. Evidence support is reported, and it can be fixed

V4 connects six actions to six packaged marker sets. On a targeted panel those sets can have too few assayed genes; every expression multiplier is then exactly 1 and the rulebook depends on labels only. This was silent. `vtx5 inspect`, `vtx5 prepare` and every schema-3.0 run report now include `action_evidence`: for each action, the cells offered it, the cells with usable evidence, the cells whose hazard changed, and a warning when declared evidence is unusable for all of them.

An optional `evidence_policy` in a schema-3.0 context declares which evidence may modulate which action for which lineage. Evidence can be a gene set scored by the simulator or a per-cell score computed outside it and imported with a provenance manifest. See [EVIDENCE_POLICY.md](EVIDENCE_POLICY.md).

```bash
vtx5 evidence-policy my_policy.json            # the v4 table as an editable starting point
vtx5 prepare object.json context.json prepared.json \
  --lineage-map lineage_map.csv --evidence-policy my_policy.json
vtx5 expression-audit prepared.json audit.json
vtx5 evidence-decisions RUN_DIRECTORY --output association.json
```

Evidence redistributes events among the cells of a type; it does not change the type's mean rate. Expect it to change which cells act, not aggregate counts. `evidence-decisions` measures that directly from a saved run.

### 2. A prompt that states how well supported each value is

With an evidence policy, requests use prompt `vtx-context-guidebook-5.0-evidence-support`. Without one, requests use the v4 evidence format. Explicit optional sample details are forwarded in either format; existing contexts without them retain the v4 request bytes. See [JEV_DECISIONS.md](JEV_DECISIONS.md).

### 3. Measuring Jev's decisions

- `vtx5 jev-probe plan|run|analyze` asks a frozen set of t0 questions under controlled manipulations, with independent repeats, and reports divergence between arms next to divergence between identical repeats.
- Every Jev run reports `jev_event_mass`: the event probability declared by the kernel, stated by Jev and applied.
- `--jev-pooling event_share` is an opt-in pooling rule that keeps the kernel's event probability. The default is the v4 rule.

### 4. A smaller decision log

Schema-3.0 runs store frozen t0 evidence once, in `guidebook.json`, and each decision carries a reference with the SHA-256 of the evidence it stands for.

```bash
vtx5 decisions RUN_DIRECTORY --output decisions_inline.jsonl   # the v4 records, byte for byte
vtx5 run prepared.json out --decision-log inline               # write the v4 layout directly
```

In Python, `virtualtissue_xenium.v2.simulation.read_decisions(run_dir)` yields the v4 records. The reader accepts only `guidebook.json` in the run folder and fails if it is missing or if its evidence differs from the hashes recorded in the log. This catches mix-ups and accidental edits; it is not protection against someone rewriting both files. Keep `guidebook.json` with `decisions.jsonl` when moving a run. Scripts that parse `decisions.jsonl` directly and use `observation.t0_evidence` should switch to the reader or the export. Schema-2.0 logs are unchanged.

## Upgrading

```bash
git pull --ff-only
pip install -e '.[xenium,test,spatial,analysis]'
vtx5 --version
```

V4 prepared bundles and saved runs remain valid. A plain `pytest` now runs without the `legacy` extra; legacy-engine test modules are left out with a header line naming the missing packages.

## Where things live

| Path | Contents |
|---|---|
| `src/virtualtissue_xenium/engine/` | Legacy physical engine (`vtx`), vendored unchanged |
| `src/virtualtissue_xenium/v2/` | Current workflow: CLI, inputs, guidebook, scoring, kernel, provider, ensembles, replay |
| `src/virtualtissue_xenium/v3/` | Schema-3.0 data and helpers: ontology, kinetics, signatures, geometry, ROI catalogue, adapters |
| `src/virtualtissue_xenium/v5/` | Evidence policy and measured scores (`evidence.py`), provider probe (`probe.py`) |
| `examples/v5/luad_tls/` | Worked example for a LUAD TLS region; the data are not packaged |

The folder names record when a part was introduced, not which command uses it: `vtx5`, `vtx4` and `vtx2` all run `v2/cli.py`.

## Known limits

- A measured score has no panel null model. A declared type with fewer cells than `min_reference_cells` in the selected region stays neutral for that score.
- Evidence is frozen at t0. Simulated state changes do not update it.
- With the default pooling, Jev's REST preference still sets how often events happen, so the event rate depends on the time step. `jev_event_mass` quantifies this in each run; compare Jev arms only at one time step, or use `event_share`.
- `frames.json` still duplicates `frames.jsonl`.
- No live Jev request was made while preparing this release. The probe and the new prompt were exercised with an in-process provider only.
