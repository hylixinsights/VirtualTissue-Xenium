# LUAD TLS region with measured GO scores

A worked v5 example on real data: one geodesic TLS candidate region of a lung adenocarcinoma Xenium section (302-gene panel), using existing labels, SPLIT-corrected expression, supplied cell boundaries and 20 saved GO Biological Process UCell scores.

**The data are not part of this repository.** The commands below assume the data folder sits beside the checkout as `../data/Xenium`, holding `LUAD3_TLS_GO20_UCell.h5ad`, `cell_boundaries.csv`, `h5ad_validation.json` and `TLS_number_of_cells.csv`. Nothing here calls Jev unless you run the last section.

## Files

| File | Purpose |
|---|---|
| `build_inputs.py` | Builds one region's raw object from the H5AD and boundaries. Verifies the H5AD checksum and the region size; changes no label; keeps uncalculated scores as null. |
| `lineage_map.csv` | Explicit mapping of `first_type` labels to model lineages. |
| `lineage_map_strict.csv` | The same with `Stromal` kept REST-only, for sensitivity. |
| `context.json` | Declared context: 15-minute step, boundary distance, `tiered_panel_auc`, recognition on. |
| `evidence_policy_go20.json` | Seven rules connecting saved GO scores to already-permitted actions. |

## Run it

From the repository root, in an environment with the `xenium` extra:

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -e '.[xenium,test,spatial,analysis]'

E=examples/v5/luad_tls
python $E/build_inputs.py --data ../data/Xenium --candidate 5 --output ../runs/B5/object.json

vtx5 prepare ../runs/B5/object.json $E/context.json ../runs/B5/prepared.json \
  --lineage-map $E/lineage_map.csv --evidence-policy $E/evidence_policy_go20.json

vtx5 expression-audit ../runs/B5/prepared.json ../runs/B5/audit.json
vtx5 run ../runs/B5/prepared.json ../runs/B5/local_seed1 --steps 24 --seed 1
open ../runs/B5/local_seed1/replay.html      # macOS; otherwise open the file in a browser
```

For a different sample, first run the English questionnaire and use its output as the context argument to `prepare`:

```bash
vtx5 context-init my_sample_context.json --template examples/v5/luad_tls/context.json
```

The questions cover tissue/organ, pathology, anatomical region, developmental stage, relevant cell states and known secretory populations. Review the template's existing treatment declarations, sampling design, recognition and numerical settings before reuse. Optional restrictions are explicit `LINEAGE:ACTION` pairs with a required reason; no restriction is inferred from an organ name. `SECRETE` means abstract inflammatory/suppressive signalling, not mucus. See [sample context and biological restrictions](../../../docs/v5/SAMPLE_CONTEXT.md).

Larger regions need `--max-cells` on `prepare` (B4 has 3,725 cells, B22 has 3,008). A 24-step local run of B5 takes about two minutes on a laptop; time grows with the number of cells.

To see what v4 does on the same region, build the object without the saved scores and prepare it without a policy. `prepare` then lists every expression-driven action as label-only:

```bash
python $E/build_inputs.py --data ../data/Xenium --candidate 5 --no-scores --output ../runs/B5/object_v4.json
vtx5 prepare ../runs/B5/object_v4.json $E/context.json ../runs/B5/prepared_v4.json --lineage-map $E/lineage_map.csv
```

### Controls and a paired ensemble

```bash
vtx5 control ../runs/B5/prepared.json ../runs/B5/beta_zero.json --kind beta_zero
vtx5 control ../runs/B5/prepared.json ../runs/B5/shuffle.json --kind shuffle_expression_within_type --seed 11
echo '{"arms": {"baseline": "prepared.json", "beta_zero": "beta_zero.json", "within_type_shuffle": "shuffle.json"}}' > ../runs/B5/arms.json
vtx5 ensemble ../runs/B5/arms.json ../runs/B5/ensemble --horizon-hours 6 --seeds 0 1 2 3 4
vtx5 evidence-decisions ../runs/B5/ensemble/runs/baseline/seed_0 ../runs/B5/ensemble/runs/beta_zero/seed_0 \
  --output ../runs/B5/association.json
```

Do not expect the aggregate readouts of this ensemble to separate the arms. The multiplier preserves each cell type's mean rate, so evidence changes which cells act, not how many. With five seeds over six simulated hours on B5, no paired difference in `paired_differences.csv` stood out from seed-to-seed variation. `evidence-decisions` shows the effect directly: in the baseline run with seed 0, B-cell decisions that selected ACTIVATE had a mean B-cell-activation z 0.89 above the decisions offered it (standard error 0.14); under `beta_zero` the difference was -0.18 (0.15).

## What to expect for B5

`prepare` prints `preflight.action_evidence`. With the policy:

| Action | Cells offered | With usable evidence | Evidence |
|---|---:|---:|---|
| PROLIFERATE | 604 | 581 | mitotic cell cycle |
| ACTIVATE | 370 | 335 | T cell activation, B cell activation, inflammatory response, antigen processing |
| SECRETE | 616 | 36 | type II interferon production (T cells only) |
| KILL | 3 | 0 | T cell mediated cytotoxicity; B5 has 7 CD8 cells, too few for a reference |
| APOPTOSIS, EXHAUST, MOVE, CLEAR_CORPSE, RECOVER | | 0 | no rule declared |

Without the policy every row has zero cells with usable evidence. In B22 (64 CD8 cells) the KILL rule is usable for all 53 CD8 cells offered it.

`expression-audit` on B5 reports 581 of 736 cells whose action probabilities differ from the label-only model.

## Decisions you should review

**Labels.** `first_type` is used as declared. A second column, `cell_type_fine`, agrees with it for only about 57% of cells at coarse level. Which column is used changes the permitted actions of many cells. `vtx5 annotation-agreement` can quantify this for a region.

**Lineage map.** `Mast_NK` and `Undefined` are REST-only because mast cells and NK cells have different permissions. `Stromal` is mapped to `fibroblast`, which grants only PROLIFERATE, APOPTOSIS and SECRETE; `lineage_map_strict.csv` keeps it REST-only.

**QC.** `min_transcripts: 20` is applied to SPLIT-corrected counts. In B5 it makes 118 of 736 cells REST-only (553 of 3,008 in B22). Lower it in `context.json` if that is not intended.

**Recognition.** `recognition: true` assumes CD8 cells recognise tumour cells. Set it to false to remove KILL.

**Rules.** Each rule in `evidence_policy_go20.json` carries its rationale and is marked `requires_biological_review`. Points to check:

- *Mitotic cell cycle → PROLIFERATE.* The 18 scored genes include cycling markers (MKI67, UBE2C, CCND1, WEE1) and broad regulators (EGFR, CTNNB1, TP53, RHOA, TGFB1). Within a cell type the score may track signalling as much as cycling.
- *No APOPTOSIS rule.* The 29 scored apoptotic-signaling genes mix pro- and anti-apoptotic regulators.
- *No EXHAUST rule.* There is no GO20 set for exhaustion and only 4 of the 6 packaged marker genes are assayed. Setting `"min_overlap": 4` in the context and adding a `gene_set` rule for `exhaustion_associated` would enable it.
- *MOVE and CLEAR_CORPSE.* Leukocyte migration and phagocytosis scores exist. A rule for them has no effect unless you also declare a beta for those actions in `scoring.beta`, because the packaged beta table does not list them.

**Region references.** Scores are standardised within the selected region. The saved UCell scores cover all 44,134 candidate cells, but a region with few cells of a type still cannot standardise that type.

## Probing Jev (costs money; needs your key)

```bash
vtx5 cost ../runs/B5/prepared.json --steps 24
vtx5 jev-probe plan ../runs/B5/probe --prepared ../runs/B5/prepared.json --cells-per-lineage 6 --repeats 3
```

For B5 the plan covers 45 cells in 5 batches: 120 requests, about 1.9 million input tokens by the package's estimate. Review `probe/payloads.json`, then:

```bash
export JEV_API_KEY=...
vtx5 jev-probe run ../runs/B5/probe --max-http-attempts 120
vtx5 jev-probe analyze ../runs/B5/probe
```

See `docs/v5/JEV_DECISIONS.md` for how to read the result.
