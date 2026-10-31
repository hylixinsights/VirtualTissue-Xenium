# VirtualTissue v3 — handoff and how to continue

Status as of 2026-10-31. Read this file first in any new session; it is the single entry point.

## 1. Where everything is

| What | Where |
| --- | --- |
| Repository | https://github.com/hylixinsights/VirtualTissue-Xenium |
| Working branch | `v3/phase0-hygiene` (all v3 work so far; nothing merged to `main`) |
| Pull request (draft) | https://github.com/hylixinsights/VirtualTissue-Xenium/pull/4 |
| Full v3 plan (14 sections, task IDs, dependency diagram) | https://claude.ai/artifact/JWxubyYoAwN7eajoBdC37g (Claude Docs; private to the owner) |
| This handoff | `docs/v3/HANDOFF.md` |
| Changelog of v3 work | `CHANGELOG.md`, section "Unreleased (v3 phase 0 …)" |
| Scoring + kinetics contract (user-facing) | `docs/v2/CONTEXT_AND_GUIDEBOOK.md` §4a |
| Original critical review (v0.1.3) | `docs/v2/CODE_REVIEW_V0_1_3.md` |
| Validation plan (gates) | `docs/v2/VALIDATION_PLAN.md` |

### Source map (v3-relevant files)

| File | Role | Touched by |
| --- | --- | --- |
| `src/virtualtissue_xenium/v2/context.py` | Validates `context.json`: `kinetics`, `scoring`, modifiers | A, B |
| `src/virtualtissue_xenium/v2/guidebook.py` | Builds per-cell actions; legacy weights or per-hour hazards; QC gate | A, B |
| `src/virtualtissue_xenium/v2/scoring.py` | Depth-conditioned within-type recovery-AUC, z, normalised β multiplier | B |
| `src/virtualtissue_xenium/v2/kinetics_v3.json` | Per-hour hazards with anchors and confidence (**DRAFT**) | A |
| `src/virtualtissue_xenium/v2/simulation.py` | Kernel: graph, gated menus, competing-risk step, commit, run, preflight | A, C-1, E |
| `src/virtualtissue_xenium/v2/provider.py` | Jev payload, hidden-prior mode, cost estimate, journal | E |
| `src/virtualtissue_xenium/v2/controls.py` | Null-model bundles (`vtx2 control`) | F |
| `src/virtualtissue_xenium/v2/cli.py` | `vtx2 prepare / inspect / run / cost / control / compare / dataset` | all |
| `examples/v2/context_*.json` | Example contexts (now hazard mode, 0.5 h, depth scoring) | A, B |
| `examples/v2/gene_sets_mouse_brain_demo.json` | Illustrative, unreviewed mouse sets | I-4 |
| `tests/test_v2_kinetics.py` | Lane A, incl. Δt-invariance test (~45 s) | A |
| `tests/test_v2_scoring.py` | Lane B, incl. depth-confound test | B |
| `tests/test_v2_controls.py` | Controls + compare | F |
| `src/virtualtissue_xenium/engine/` | Preserved, hashed v0.1.3 physical engine — **do not edit** | — |

## 2. Status by task ID (IDs from the plan)

| ID | Task | Status | Commit / note |
| --- | --- | --- | --- |
| A-1 | Hazard → per-step probability (competing risks) | Done | `7da734d` |
| A-2 | Literature-anchored hazard table | **Draft — needs your review** | `kinetics_v3.json`: verify each `anchors.*.citation`, add DOIs; `assumed` rows stay in sensitivity |
| A-3 | Δt invariance test | Done (pytest) | `tests/test_v2_kinetics.py`; CI figure not done |
| A-4 | Implied doubling time / admissible Δt in inspect | Done | `vtx2 inspect` → `kinetics` |
| B-1 | β effect size, `exp(β z)` normalised to reference mean | Done | `ff60476` |
| B-2 | Depth-conditioned within-type scoring | Done | `scoring.py` |
| B-3 | Signature sets as data with provenance; add IFN, TGFβ, hypoxia, MHC-II, Treg, NK, TRM | **Not done** | Needs curated sources (MSigDB IDs) — human review |
| B-4 | Species handling | Partial | Mouse example + refusal of human defaults done; reviewed ortholog table not done |
| B-5 | QC gate `min_transcripts` | Done | REST-only, excluded from references |
| C-1 | Nearest / seeded kill and clearance target | Done | `d3b4a53` |
| C-2 | Timed activation/exhaustion, implicit deactivation | Not done | Next code lane |
| C-3 | Motility (`MOVE`) | Not done | |
| C-4 | Daughter placement + local density gate | Not done | |
| C-5 | Corpse slots, corpse lifetime | Not done | |
| C-6 | Stream frames to `frames.jsonl` | Not done | |
| D-1..D-4 | Ontology JSON, tissue profiles, `lineage_map.csv` | Not done | Lineages still in `context.py`/`guidebook.py` |
| E-1 | Hidden-prior Jev mode + JS divergence | Done | `3df7b87` |
| E-2 | Determinism replicates | Not done | Needs API key |
| E-3 | `vtx2 cost` | Done; payload compression not done | |
| E-4 | α sweep | Not done (part of F-10) | |
| E-5 | Batch-order leakage test | Not done | Needs API key |
| F-1, F-2, F-5, F-6, F-7, F-8, F-9 | `vtx2 control` kinds | Done | `7a1691f`, `ff60476` (beta_zero) |
| F-3 | Hidden prior as control | Done (via E-1) | |
| F-4 | Scrambled context for Jev | Not done | |
| F-10 | Sensitivity sweep (LHS, Sobol) | Not done | Needs `vtx2 ensemble` |
| F-11 | Δt invariance | Done (A-3) | |
| — | `vtx2 ensemble` (multi-seed runner + envelopes) | **Not done — highest-value next tool** | |
| G-1..G-5 | Real data: annotation, runs, structural validation | Not done | Needs laptop: downloads, annotation, key |
| H-1..H-5 | Platform adapters (Atera, Visium HD, Stereo-seq) | Not done | Start with H-1 `measurement` record |
| I-1, I-2, I-4 | Cost-aware README, checkpoint wording, mouse example | Done | `61bdaef` |
| I-3 | Rename `engine/xenium_runtime/v2` | **Decision needed** | Would change import paths of hashed legacy files |
| I-5..I-9 | Minor fixes, LIMITATIONS.md, paper/ folder, CI, data/ | Not done | |

## 3. What to do next, in order

Sequential spine (each step needs the previous one):

1. **Review and merge PR #4** after running the full suite on your laptop (section 5).
2. **Review `kinetics_v3.json`** (A-2): open every `anchors` paper, confirm the number, add a DOI, adjust `confidence`.
3. **`vtx2 ensemble`**: N seeds × arms → per-readout medians and 95% bootstrap envelopes; readouts from a frozen `analysis/readouts.yaml`. Every figure depends on it.
4. **Invariance gate** before any Jev spend: A-3 passes (done) and F-8 (`no_modifier_treatment`) is identical to control in local mode.
5. **Real-data local runs** (G-3), then Jev arms (shown + hidden prior), then figures (plan §13).

Parallel lanes (independent; can be separate sessions or people once step 1 is merged):

- **C** kernel (C-2 → C-3/C-4 → C-5/C-6) — `simulation.py` only.
- **D** ontology JSON — `context.py`, `guidebook.py`, new `ontology/` folder. Pair with B-3.
- **E-3 compression, F-4** — `provider.py`, `controls.py`.
- **H** adapters — new `src/virtualtissue_xenium/v2/adapters/`; H-1 first.
- **G-1 annotation** (human work) — start now; it is the longest lead time.
- **I-5..I-9** hygiene.

Collision rule: two parallel sessions must not edit the same file. C owns `simulation.py`; D owns lineage constants; E owns `provider.py`.

## 4. Decisions waiting on you

1. I-3: rename the legacy engine directory, or keep it and disambiguate in the paper?
2. A-2: accept, change or replace each hazard; which rows stay `assumed`.
3. Default β: keep 1.0 globally, or per action (e.g. only PROLIFERATE and ACTIVATE)?
4. Defaults for real data: `min_reference_cells` 30, `depth_bins` 10, `min_transcripts` (suggest 20 for Xenium 5K).
5. Venue/framing decided after the controls (plan §13 P-4).

## 5. Commands

```bash
git clone https://github.com/hylixinsights/VirtualTissue-Xenium.git
cd VirtualTissue-Xenium && git switch v3/phase0-hygiene
python3.12 -m venv .venv && source .venv/bin/activate
pip install -e '.[xenium,test]'
pytest tests/test_v2_*.py -q            # ~1 min; kinetics test ~45 s

# Quickstart (hazard mode, 0.5 h steps; 48 steps = 1 day)
vtx2 prepare examples/v2/synthetic_cells.json examples/v2/context_melanoma_control.json outputs/v2/control.json
vtx2 inspect outputs/v2/control.json      # doubling times, admissible Δt, scoring coverage
vtx2 run outputs/v2/control.json outputs/v2/control --steps 48 --seed 21

# Controls
vtx2 control outputs/v2/control.json --kind list
vtx2 control outputs/v2/control.json outputs/v2/ctl_shuf.json --kind shuffle_expression_within_type --seed 1
vtx2 run outputs/v2/ctl_shuf.json outputs/v2/ctl_shuf --steps 48 --seed 21
vtx2 compare outputs/v2/control outputs/v2/ctl_shuf

# Jev (after inspection)
vtx2 cost outputs/v2/control.json --steps 48 --max-http-attempts 100
vtx2 run outputs/v2/control.json outputs/v2/control_jev_hidden --mode jev --jev-prior hidden \
  --steps 48 --seed 21 --max-http-attempts 200
```

Notes:

- The 24-cell synthetic example has 6 cells per type, below `min_reference_cells`, so expression has **no** effect there by design (`inspect` shows `insufficient_reference`). Use real data or `tests/test_v2_scoring.py::population()` to see the effect.
- Contexts without `kinetics`/`scoring` keys reproduce 2.0.0a1 exactly (`legacy_per_decision`, `rank_legacy`).
- Run fingerprints changed because of C-1 (kill target rule). Old run directories cannot be resumed with this branch.

## 6. Prompt to start a new conversation

Paste this into a new Claude conversation (with this repository attached or linked):

> I am continuing the VirtualTissue v3 plan. Repository: https://github.com/hylixinsights/VirtualTissue-Xenium, branch `v3/phase0-hygiene`, draft PR #4. Read `docs/v3/HANDOFF.md` first; it lists what is done (by task ID), what is next, the parallel lanes and the open decisions. The full plan is at https://claude.ai/artifact/JWxubyYoAwN7eajoBdC37g. Do not change the CLI shape, bundle format, replay layout or the hashed legacy engine. Work on: **[lane / task ID]**. Run `pytest tests/test_v2_*.py -q` before committing, update `CHANGELOG.md` and the status table in `docs/v3/HANDOFF.md`, and push to the same branch.
