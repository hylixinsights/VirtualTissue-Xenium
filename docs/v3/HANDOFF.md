# VirtualTissue v3 implementation handoff

Updated 2026-11-02 for the **3.0.0a4 alpha prerelease**. Start with the repository README and `START_TEST.command`; its result page is the single entry point for software acceptance. The prepared local `V3_teste` copy has an equivalent enclosing launcher. This remains a research implementation, not a biologically validated model.

This update extends PR #4 from the historical base
`4416f9765bb09d6b678758a9cfd7fa90070ae5e5`. Its audited schema-2.0 kinetic table
has been incorporated here with all 46 numerical hazards unchanged. The local
schema-3.0 table lives separately at `v3/data/ontology/kinetics_v3.json` and has
91 hazards across 22 lineages. Its scientific review accounts for the timed
states already implemented. Earlier phase-0 handoffs describe a smaller feature
inventory; this document records the v3 implementation and its confirmed scope.
See [RATE_REVIEW.md](RATE_REVIEW.md) and the [Mac launchers](../../scripts/mac/README_MAC.md).

**Current author decision:** use existing labels and segmentation polygons from
the supplied object, with no mandatory H&E or independent annotation review.
Present ROI composition and let the user choose. The three hypotheses and their
limits are frozen in [SCIENTIFIC_SCOPE.md](SCIENTIFIC_SCOPE.md); this supersedes
the original annotation/ROI prerequisites. No claim of a proven seven-day tissue
prediction is permitted. See [ROI_WORKFLOW.md](ROI_WORKFLOW.md) and
[EXPRESSION_POLICY.md](EXPRESSION_POLICY.md).

The optional user-object profile adds `tiered_panel_auc` for rare cell types,
explicit geometric QC and boundary distances, `roi-catalog`, `select-roi` and
`expression-audit`. The original strict scoring method remains available.

## Implemented

| Plan IDs | Implementation and evidence |
|---|---|
| Phase 0, H-1/H-2 | Context/bundle 3.0 alongside 2.0, measurement/unit/organism contract, platform profiles, versioned ontology and fixed readouts |
| A-1/A-4 | Competing-risk rate conversion, step guard, gated t0 hazards and expected events/hour by lineage |
| A-3/F-11 | Original regression plus `tools/validate_clock_v3.py`: five seeds, matched horizons, analytic discrete branching expectation, MC uncertainty and figure; restricted local growth scenario |
| B-1/B-2/B-5 | Strict sample-isolated references; null default for undersized depth strata; explicit whole-type fallback uses a consistent partition; QC, beta-zero and separate fallback/informative coverage |
| B-3/B-4 | Versioned human/mouse sets with sources, exact symbols, explicit reciprocal MGI ortholog records; no implicit uppercasing in v3 |
| C-1..C-5 | Seeded target ties; timed activation and exhaustion recovery gate; persistent movement bouts and optional recent replay tracks; CD8 activation gate for division; daughter placement/soft repulsion/crowding; corpse lifetime; separate abstract suppressive channel |
| C-6 | Incremental frames.jsonl and decision logs, streamed compatibility frames.json, offline replay sidecar in 50-frame chunks with a two-chunk browser cache, and a declared maximum horizon |
| D-1..D-4 | 22 lineage definitions, rate/β tables, tissue profiles, reviewed CSV lineage adapter with hashes; comparison rejects ontology differences unless explicitly allowed and reports available recorded row/hash changes |
| E-1/E-3, F-4 | Prior-hidden forces alpha=1; compact v3 payload; donor-text scrambled context preserves local numeric mechanisms; cost estimator/usage and ensemble ledger |
| E-2/E-5 infrastructure | `jev-audit` measures matched repeated-payload/order variance from durable journals; sampling parameters recorded as unknown when unavailable; mock tests |
| F-1..F-9 | Previous controls retained; label permutations preserve sample/purity eligibility; explicit donor-context control; nulls retained |
| Ensemble, F-10 | `ensemble` runner, paired differences/dz with explicit degenerate-bootstrap handling, median bootstrap CI distinct from seed envelope; LHS exploratory and valid A/B/AB_i Sobol sensitivity runner; no local alpha sweep |
| G-1 infrastructure | Deterministic ROI candidates from reviewed input, exact-ID reviewer concordance/κ and disagreement exports |
| G-4 infrastructure | Cross-L with declared rectangular window/translation correction, cluster sizes and held-out descriptive envelopes |
| H-3 | Native Atera CSR Zarr, Visium HD binned H5/parquet and Stereo GEM readers and blank annotation templates, plus existing Xenium/H5AD; fixture tests |
| H-4 | Fixed-bin occupancy kernel; purity gate; separate occupancy readouts; rejects pooling with cell-agent runs |
| I-5/I-7/I-8/I-9 | Early H5AD count bound; symmetric v3 graph; paper figure renderer/Makefile; CI tests/coverage/clock artifact; packaged data version and dataset-registry hash |
| I-2 | Remove the unused checkpoint output in schema 3.0; keep the historical schema 2.0 diagnostic snapshot and deterministic reconstruction from t0/journals |

## Decisions implemented explicitly

- The original v3 default is sample-specific scoring with `reference_fallback: null`, minimum 30 references and platform-specific QC (Xenium: 20 transcripts). The explicit user-object profile uses `tiered_panel_auc`, including a bounded weak panel heuristic when cell references are insufficient. V2 scoring semantics remain unchanged.
- `max_step_probability` retains its compatibility name but bounds **sum(hazards) × dt**, not a probability directly. V3 rate metadata states the actual competing-risk formula.
- All default numerical biological rates remain **assumed**. All eight bibliographic identities now have verified DOI/source links; unsupported numerical descriptions were corrected or removed. Quantitative rate verification remains pending. See `RATE_REVIEW.md`.
- The old physical engine directory was retained to preserve hashes/imports. This resolves the implementation path conservatively without claiming the naming issue is scientifically significant.
- Hidden-prior Jev uses raw returned preferences, alpha=1. This prevents prior contamination but does **not** inherit local rate-clock calibration.
- An unknown bin is occupied unless occupancy=0 is explicitly supplied. Mixed/unreviewed bins are not silently promoted to cells or emptied.

## Remaining work, not claimed complete

1. **A-2 and scientific B/D review:** extract and review quantitative rate evidence, appropriateness of transfer to the target tissue, signature membership, ortholog choices, β and platform defaults with domain reviewers. Verified bibliography is not tissue calibration.
2. **Real inputs under the revised scope:** import existing object labels/polygons, document their origin and uncertainty, and record the user's ROI choice. Existing local CSVs are unknown-only engineering inputs. H&E and independent annotation review are optional, not prerequisites; supplied labels do not establish malignancy or true cell boundaries.
3. **G-3/G-4/G-5, E-2/E-4/E-5 and H-5 experiments:** real local and Jev arms, independent journal replicates, paid-call budget, held-out structural validation, alpha sweeps and real cross-platform comparisons. No live Jev calls were performed here.
4. **H adapter scope:** direct GEF import and mask-to-pseudocell aggregation remain unimplemented; use explicit vendor exports/segmented H5AD. Native new-platform readers have synthetic format tests, not public-dataset validation.
5. **Publication and historical documentation:** inherited v2 documentation remains historical; current limitations are centralized in `LIMITATIONS.md`. The figure renderer covers ensemble readouts and the clock check; final manuscript panels, spatial/sensitivity figures and study claims require the corresponding experiments.

## Replay files

V3 `replay.html` loads adjacent `replay_frames/chunk_*.js` files as needed, including when opened directly from disk. Keep the HTML and the entire sidecar folder together when copying a replay. Its browser cache holds at most two 50-frame chunks; metadata and cell identifiers remain in memory, so this is a bounded frame cache rather than constant total memory. `frames.jsonl`, `frames.json` and decisions remain complete audit files and grow with run size. V2 replay stays self-contained.

The optional `Recent tracks` toggle is off by default and draws only MOVE events
available in the two cached chunks, up to 100 recent steps. It is disabled for
well-mixed and lattice runs. Playback/loading/seek/cache/selection controls were
tested using the real JavaScript player with a DOM test double; visual browser
inspection was unavailable in this environment and is not claimed.

## Compatibility and resuming

The a3 scientific metadata changes alter kinetic hashes even where every
numerical hazard remains unchanged. Prepare new bundles for a3; keep a2 bundles
with their original a2 package and run files. Do not rewrite historical hashes
or silently reuse prior preparations. The a3 catalogue and demonstration are
regenerated separately from the preserved a1/a2 artifacts.

The a2 Jev prompt is `vtx-context-guidebook-3.1-evidence-basis` (with a hidden-prior
variant). It explicitly identifies weak panel evidence. Journals from another
prompt version cannot be reused as if they were generated with this prompt.

The plan's clock proposal was refined: the current gate uses exact discrete
branching moments, five independent simulation seeds and a continuous-time
bias bound. Equal seed labels across step sizes do not create coupled random
streams per cell-hour. The test covers homogeneous local growth; integrated
state/interaction scenarios and Jev require additional timestep experiments.
Likewise, LHS is used only as an exploratory screen; Sobol indices use an A/B
pick-freeze design, rather than incorrectly treating 200 LHS points as a Sobol
estimator. These are explicit departures from the draft statistical recipe.

Direct measured-contact-versus-capped-graph audit remains study-specific work.
The present proximity readout uses centroid distances independently of the graph
cap, even when the action graph uses polygon-boundary distance; it is not that
graph validation. Ki67/histology assessment is optional in the revised scope,
and cannot be claimed as performed when absent.

The engine under `src/virtualtissue_xenium/engine` is unchanged. Schema 2.0 remains supported. A real kernel bug was fixed in both general schemas: KILL and APOPTOSIS in the same round can no longer count one death twice. Kernel fingerprints changed to `vtx-discrete-context-2.0.1` and `vtx-cell-agent-3.0.1`; old runs are not resumable under these versions. Restore a pinned checkout to reproduce/resume old runs.

Prepared v3 bundles are bound to the packaged ontology hashes. If versioned data change, restore the original data or rebuild; an old fingerprint cannot authorize new permissions/radii. Journals remain durable and are never automatically retried after an unknown/failed attempt.

Schema 3.0 no longer writes `checkpoint.json`. `--resume` reconstructs from t0
using preserved decisions/journals; it does not restart from a saved state.
The schema 2.0 checkpoint remains a diagnostic artifact for compatibility.

## Start commands

```bash
pip install -e '.[xenium,test,spatial,analysis]'
python -m pytest tests/test_v2_*.py tests/test_v3_*.py -q
python examples/v3/make_demo.py --output outputs/v3/demo
vtx2 inspect outputs/v3/demo/control.json
vtx2 ensemble outputs/v3/demo/arms.json outputs/v3/demo/ensemble \
  --horizon-hours 24 --seeds 0 1 2 3 4
make figures ENSEMBLE=outputs/v3/demo/ensemble FIGURES=outputs/v3/demo/figures
```

The demo is synthetic and explicitly uses two depth bins to exercise reference scoring in a small population. It does not establish a biological expression effect. See `analysis/README.md`, `docs/v3/ONTOLOGY.md`, `ANNOTATION.md`, `PLATFORMS.md` and `LIMITATIONS.md` before real-data use.
