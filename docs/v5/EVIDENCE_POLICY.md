# Evidence policy and measured scores

An evidence policy states which experimental evidence may modulate which action. It is optional. Without it, a schema-3.0 context builds the v4 rulebook.

## Why it exists

The v4 rulebook uses one global table: PROLIFERATE follows `cell_cycle`, ACTIVATE follows `activation`, EXHAUST follows `exhaustion_associated`, APOPTOSIS follows `apoptosis_associated`, SECRETE follows `inflammation` and KILL follows `cytotoxic`. Two problems follow on real data.

1. **Coverage.** A signature needs `min_overlap` assayed genes (default 5). On a 302-gene lung panel the six sets have 2, 2, 4, 0, 3 and 2 assayed genes. Every multiplier is exactly 1, `beta_zero` and expression shuffles equal the baseline, and Jev receives null for every action-relevant signature.
2. **Lineage.** The same activation marker set is applied to T cells, B cells and myeloid cells.

Check your own bundle before anything else:

```bash
vtx5 inspect prepared.json      # read action_evidence.by_action, label_only_actions, warnings
```

## The policy

```json
{
  "policy_version": "v5-evidence-policy-1",
  "name": "my_policy_v1",
  "review_status": "requires_biological_review",
  "provider_evidence": "all",
  "rules": [
    {
      "action": "ACTIVATE",
      "lineages": ["T_CD4", "T_CD8"],
      "evidence": "GOBP_T_CELL_ACTIVATION_UCell",
      "kind": "measured_score",
      "direction": "supports",
      "rationale": "Why this evidence bears on this action for these lineages.",
      "review_status": "requires_biological_review"
    }
  ]
}
```

| Field | Meaning |
|---|---|
| `action` | A non-REST action. |
| `lineages` | A list of lineage names, or `"*"` for every lineage permitted that action. |
| `evidence` | A gene-set name (`kind: gene_set`) or a measured-score name (`kind: measured_score`). |
| `direction` | `supports` (default) or `opposes`; `opposes` uses `-beta`. |
| `rationale` | Required text. It is stored with every adjusted hazard. |
| `review_status` | Defaults to `requires_biological_review`. |
| `provider_evidence` | `"all"` or a list of evidence names shown to Jev. |

Put the policy in the context (`"evidence_policy": {...}`) or pass `--evidence-policy FILE` to `prepare` or `select-roi`.

Rules the software enforces:

- **A rule cannot add an action.** Naming a lineage that is not permitted the action is an error.
- **One evidence source per lineage and action.** Two rules that address the same pair are rejected. Scores computed from the same counts are correlated, and multiplying two of them would count the evidence twice.
- **Declared rules replace the v4 table entirely.** An action without a rule is label-only, and the report says so.
- **Every evidence name must exist** as a gene set of the bundle or as an entry of the measured-evidence manifest.

What a rule changes is only the multiplier applied to an assumed hazard:

```text
hazard = base_hazard × exp(±beta × z) / mean_reference(exp(±beta × z))
```

`beta` comes from `scoring.beta` and the packaged table, as in v4. The packaged table lists PROLIFERATE, ACTIVATE, EXHAUST, APOPTOSIS, SECRETE and KILL; an action not listed there (MOVE, CLEAR_CORPSE, RECOVER) has beta 0, so a rule for it has no effect unless you also declare a beta for it in the context. Gates, crowding and the competing-risk conversion apply afterwards, unchanged.

`vtx5 evidence-policy FILE` writes the v4 table as a policy. Unedited, and with the packaged gene sets, it gives the same hazards and the same trajectory as v4; only the guidebook's provenance records differ. One difference in strictness: every evidence name a rule uses must exist, so a custom `--gene-sets` file that lacks one of the six names is rejected, where v4 skipped the missing set silently.

## Measured scores

A measured score is a per-cell number computed outside the simulator, for example a saved UCell score.

**In a JSON object**

```json
{
  "cells": [{"cell_id": "c1", "...": "...", "measured_scores": {"MY_SCORE": 0.31, "OTHER": null}}],
  "measured_evidence": {
    "scores": {
      "MY_SCORE": {
        "method": "pyUCell 0.7.3",
        "description": "T cell activation (GO:0042110); rank score on a 0-1 scale.",
        "genes_scored": ["CD3E", "CD4", "..."],
        "n_genes_set": 627,
        "caveat": "Optional, specific to this score."
      }
    },
    "caveat": "Optional, common to all scores."
  }
}
```

**From AnnData**

```bash
vtx5 prepare roi.h5ad context.json prepared.json --cell-type-key first_type \
  --lineage-map lineage_map.csv --layer purified_raw \
  --score-key GOBP_T_CELL_ACTIVATION_UCell --score-key GOBP_MITOTIC_CELL_CYCLE_UCell \
  --measured-evidence manifest.json --evidence-policy policy.json
```

- `method` and `description` are required for every score. `genes_scored` and `n_genes_set` are optional and are shown to Jev when present.
- Give the manifest once: inside the JSON object or with `--measured-evidence`, not both.
- `null` (or NaN in AnnData) means not calculated. It is never read as zero. Such a cell is neutral for that score and is not used as a reference for other cells.
- Measured scores require an evidence policy. They are carried in the prepared bundle and its fingerprint, and are not copied into simulation frames.

### Standardisation

A measured score is standardised exactly like the internal recovery AUC:

1. Cells failing the transcript-count QC are neither scored nor used as references.
2. Each cell is compared with scored cells of the same sample, the same declared type and the same depth stratum. With `reference_fallback: "cell_type"` (always the case under `tiered_panel_auc`), a type in which any depth stratum has too few scored cells is compared within the whole type instead. With the default `reference_fallback: "null"` of `depth_conditioned_auc` there is no such fallback: a stratum with too few scored cells stays `insufficient_reference`.
3. The mid-rank percentile becomes `z = probit(percentile)`, clipped to ±3.

Because the comparison is within type, a score that differs between cell types mostly through composition contributes only its within-type variation. A reference group with fewer than `min_reference_cells` scored cells in the selected region is `insufficient_reference` and stays neutral; an imported score has no panel null model, so there is no weaker fallback for it.

## Controls and audits

- `shuffle_expression_within_type` and `shuffle_expression_global` move a cell's counts and its measured scores together.
- `beta_zero` neutralises gene-set and measured evidence alike.
- `expression-audit` reports, per action, the declared evidence and whether a probability changed through its hazard or through redistribution.
- `compare` refuses arms prepared with different policies.
- `evidence-decisions RUN [RUN ...] --output FILE` reads saved runs and asks whether the decisions that selected an action carried more of its declared evidence than the decisions that were offered it.

### What to expect from an ensemble

The multiplier preserves the mean hazard of each reference group. Evidence therefore changes **which** cells have events, not how many events a cell type has. Aggregate readouts such as cumulative divisions or the activated fraction are not expected to separate a baseline arm from `beta_zero`, and a paired ensemble that shows no such difference is not a failure of the evidence layer. Use `evidence-decisions`, and readouts that depend on where events happen, to see the effect:

```bash
vtx5 evidence-decisions ensemble/runs/baseline/seed_0 ensemble/runs/beta_zero/seed_0 --output association.json
```

A baseline run shows a positive difference for actions with a rule; `beta_zero` shows about zero. An expression-shuffle arm still shows a positive difference, because its guidebook holds the shuffled evidence that drove its hazards: the shuffle breaks the link between evidence and position, not between evidence and decision. This readout confirms that sampled decisions follow the declared evidence, in local or Jev mode. It does not show that a rule is biologically right.

## What Jev sees

See [JEV_DECISIONS.md](JEV_DECISIONS.md). The rules, their rationales and review status are not sent, so Jev's use of evidence is independent of the rulebook's mapping.

## Limits

- A rule is an assumption. A broad gene set scored on a small panel may be dominated by a few genes that do not represent the named process. Read `genes_scored` before trusting a rule.
- Gene Ontology sets can contain regulators of both signs. A score without a clear direction should not drive an action.
- Standardisation is relative to the selected region. Changing the region changes the references.
- Nothing here calibrates a rate or validates a prediction.
