# Context and guidebook specification

The v2 guidebook separates measured baseline expression from investigator assumptions and simulated state. It gives each annotated cell a compatible action menu and an uncalibrated preference distribution. Jev can reconsider those actions using the declared context and current neighborhood. Neither distribution is a measured biological transition probability.

## 1. Annotated input

Each cell needs a unique string `cell_id`, a reviewed `cell_type`, an explicit broad `lineage`, and an `expression` dictionary of gene symbols and finite nonnegative counts. Sparse dictionaries can omit zero counts; retain the complete panel separately as `assayed_genes`. Coordinates `x` and `y`, when supplied, must be in micrometres. Give separate donors/sections distinct `sample_id` values: runtime graph edges never cross samples.

```json
{
  "assayed_genes": ["ACTB", "MKI67", "TOP2A", "PCNA"],
  "cells": [
    {
      "cell_id": "donor1_cell001",
      "cell_type": "reviewed malignant melanocyte",
      "lineage": "malignant",
      "sample_id": "donor1_section1",
      "x": 120.5,
      "y": 81.0,
      "expression": {"ACTB": 15, "MKI67": 3},
      "initial_state": {"activated": false, "exhausted": false}
    }
  ]
}
```

No label automatically establishes malignancy, antigen recognition, activation, or exhaustion. The same `cell_type` must map to one lineage; split a conflicting annotation before preparation. Unmapped populations should use `lineage: "unknown"`, whose only action is `REST`. This makes unresolved coverage visible instead of assigning an unsupported behavior.

For spatial data, neighbors are the nearest centroids within `neighbor_radius_um`, capped by `max_neighbors`. This graph does not establish membrane contact. For single-cell or single-nucleus data without spatial measurements, choose `well_mixed` explicitly: every step samples an assumed neighborhood within each sample. It does not reconstruct the tissue's original anatomy. Where measured coordinates are present, they remain fixed apart from a small daughter-cell offset; v2 has no mechanics, space exclusion, diffusion solver, or pharmacokinetics.

## 2. Required experiment declaration

Start from [the brain control](../../examples/v2/context_brain_control.json) or [the melanoma treatment example](../../examples/v2/context_melanoma_anti_pd1.json), and replace the descriptions with the actual experiment.

| Field | Meaning |
|---|---|
| `schema_version` | Exactly `"2.0"`. |
| `tissue` | Anatomical origin, including compartment or region where known. |
| `organism` | Explicit species. Default signatures support `human` or `Homo sapiens`. |
| `disease` | Actual diagnosis/state; use an explicit healthy/unknown declaration as appropriate. |
| `assay` | `xenium`, `scrna`, or `snrna`. |
| `sampling.design` | `whole_tissue`, `enriched`, or `subsampled`. This is a sampling declaration, not a claim of complete biological capture. |
| `sampling.description` | Region selection, enrichment, sorting, annotation and selection details. |
| `sampling.missing_populations` | A list of excluded/absent/uncertain populations; `[]` if none are declared. |
| `perturbations` | A list of declared interventions; `[]` is a control. Each entry requires `type`, `name`, `targets` (list), and `mechanism_assumptions` (list). |
| `recognition` | Explicit boolean licensing a CD8 killing menu. Expression and an anti-PD1 label cannot establish this assumption. |
| `time_step_hours` | Positive nominal duration of one discrete decision step. It is not a fitted biological clock. |

Each perturbation additionally accepts `status: "observed"` or `"hypothetical"`. The validator inserts `"hypothetical"` when omitted, including older context examples, and exports this canonical declaration. Use `"observed"` only when the intervention was actually applied to the source sample; this label still does not establish a causal response. Optional `dose` and `duration` are descriptive text, including units and timing where known, and optional `evidence` is a list of source URLs or text references. These fields do not create pharmacokinetics, exposures or fitted response rates. Keep the measured source-sample treatment history separate from a proposed virtual intervention when both are present.

Optional graph settings are `interaction_mode: "spatial"` or `"well_mixed"`, `neighbor_radius_um` (positive, default 30) and `max_neighbors` (integer 1–100, default 20). If omitted, Xenium defaults to spatial and sc/snRNA to well mixed; declare the mode explicitly in shared experiments. The full input panel can be supplied as `assayed_genes`; the AnnData adapter returns its complete `var_names` for this purpose. Further JSON metadata, such as disease stage, duration, dose and source references, can be retained in the context, but it does not silently become a kinetic mechanism.

Context validation fails before a Jev call for missing declarations, nonfinite values, ambiguous booleans, unknown actions and invalid interaction settings. It does not certify the truth of the biological declarations.

## 3. Compatible action menus

All supported populations have `REST`. The table lists additional candidates before local state gates. These are illustrative modeling choices requiring review for each tissue, not a comprehensive cell ontology.

| Explicit lineage | Additional compatible actions |
|---|---|
| `T_CD8` | `PROLIFERATE`, `APOPTOSIS`, `ACTIVATE`, `EXHAUST`, `RECOVER`, `SECRETE`; `KILL` only with declared recognition |
| `T_CD4` | `PROLIFERATE`, `APOPTOSIS`, `ACTIVATE`, `EXHAUST`, `RECOVER`, `SECRETE` |
| `B` | `PROLIFERATE`, `APOPTOSIS`, `ACTIVATE`, `SECRETE` |
| `myeloid`, `microglia` | `PROLIFERATE`, `APOPTOSIS`, `ACTIVATE`, `SECRETE`, `CLEAR_CORPSE` |
| `fibroblast`, `endothelial`, `epithelial`, `malignant` | `PROLIFERATE`, `APOPTOSIS`, `SECRETE` |
| `astrocyte` | `PROLIFERATE`, `APOPTOSIS`, `ACTIVATE`, `SECRETE` |
| `neuron` | `APOPTOSIS`, `ACTIVATE` |
| `oligodendrocyte` | `APOPTOSIS` |
| `unknown` | None |

Here, `neuron` and `oligodendrocyte` mean mature populations; proliferating precursors need a reviewed extension rather than relabeling them to force an action. The general `myeloid` and `microglia` proliferation weights are explicit low-weight assumptions; turn them off if incompatible with the population being modeled. `ACTIVATE` in neurons is only an abstract state flag and does not represent electrophysiological firing.

The runtime gates `KILL` on CD8 lineage, declared recognition, activated/nonexhausted simulated state and a live malignant graph neighbor. `CLEAR_CORPSE` needs a dead neighbor; `EXHAUST` needs an activated nonexhausted state; `RECOVER` needs exhaustion; `ACTIVATE` is unavailable when already activated. Proliferation attempts produce at most one daughter per step and enforce a living population cap. There is no physical space check. A licensed killing action commits one target death; no antigen/receptor kinetics or measured dose-response is implied.

## 4. Expression evidence and numerical weights

The default illustrative human signature lists are in `guidebook.py` and are exported in every guidebook. They cover cell cycle, cytotoxic transcripts, activation, exhaustion-associated transcripts, apoptosis-associated transcripts, inflammation and antigen presentation. They are small descriptive panels, not a fitted pathway model.

For each cell:

1. Normalize symbols to uppercase; do not infer orthology.
2. Rank positively detected expression values within that cell, using average ranks for ties and dividing by the number of positive features.
3. For each signature, intersect its genes with the declared assayed panel. Assign undetected assayed genes a contribution of zero, and average contributions over the assayed overlap.
4. If fewer than three signature genes are assayed, record `score: null` and `status: "insufficient_panel_overlap"`.

Each evidence entry retains assayed and detected genes, panel coverage and detection fraction. A score of zero means no positive detection among the assayed overlap under this scoring definition; it does not demonstrate biological inactivity. A null score means insufficient panel support. If the panel was not supplied, the union of input expression keys is used and explicitly marked potentially incomplete. Always retain the full feature list for real analyses.

This method is **not AUCell**. It is sensitive to dropout, panel selection, RNA depth and cell composition. Positive-rank scores are invariant to multiplying all counts in a cell by a constant, but they are not calibrated across assays, tissues or panels. Activation genes such as `FOS` and `JUN`, checkpoint transcripts and apoptosis-associated genes are nonspecific evidence; they do not establish functional states or imminent death.

The guidebook starts with explicit assumed action weights, with `REST` as a positive fallback. A supported expression score multiplies its associated action weight by `1 + score`; missing evidence leaves the assumed weight unchanged. User overrides and matched context modifiers are also assumptions. The exported `actions` keep `prior_weight`, adjusted `weight`, normalized `probability`, description and provenance. The guidebook also records exact signature lists and hashes of input, context and signatures. Per-cell evidence and action preferences are summarized by annotated cell type.

At each runtime step, local gates remove unavailable actions. An abstract signal from the previous neighborhood snapshot multiplies the available activation weight by `1 + mean_neighbor_signal`; this illustrative feedback is not a measured ligand–receptor effect. Remaining weights are renormalized. Local mode samples this distribution. Jev mode combines the guidebook prior `p` with Jev's returned preference `q` using normalized `p^(1-alpha) × q^alpha`, then samples with a reproducible cell/step seed. `alpha=0` uses guidebook preferences; `alpha=1` uses Jev preferences on the same permitted menu. Neither source is an independently calibrated likelihood, so the combination is a declared pooling rule, not Bayesian biological inference.

[TypeSafe documents](https://docs.typesafe.ai/introduction) that each question is evaluated independently against shared state. V2 follows this interface contract; the shared state contains the experiment declaration, while cell-specific evidence belongs to its question. Keep the batch setting matched between experiment arms and include `--max-questions 1` in empirical batch-size ablations where useful. Concurrent HTTP dispatch can still evaluate multiple single-cell requests in parallel.

All weights apply **per decision step**, not per hour. Changing `time_step_hours` changes the simulated meaning of the same transitions; it does not automatically rescale them as hazards. Compare matched arms with the same interval, number of steps, sampling design, seed, model, graph and priors. Calibrating time-dependent division/death/exhaustion rates is future experimental work.

## 5. Declaring disease or treatment without hardcoding success

A chronic-disease label alone does not force exhaustion. An anti-PD1 label alone does not activate every T cell or license killing. These labels are available to Jev as context, so a **prompt-only** perturbation can alter Jev's preferences, while local mode retains identical guidebook weights.

An explicit modifier additionally changes a guidebook weight. For example:

```json
{
  "context_modifiers": [
    {
      "action": "ACTIVATE",
      "lineages": ["T_CD8"],
      "multiplier": 2,
      "reason": "Illustrative sensitivity assumption, not measured anti-PD1 efficacy.",
      "when": {"perturbation": "anti-PD1"}
    }
  ]
}
```

`when` supports exact matches on `perturbation` name, `tissue`, and `disease`. No synonym matching or diagnosis inference occurs. In a matched untreated context without the named perturbation, the rule remains recorded as inactive. Keep both a prompt-only arm and a modifier sensitivity arm when separating Jev context effects from explicit numerical assumptions.

Optional `action_priors` directly overrides compatible weights by lineage, for example `{"malignant": {"PROLIFERATE": 0.04, "APOPTOSIS": 0.002}}`. A zero weight removes that action. Overrides cannot add a lineage-incompatible action or bypass recognition licensing. `REST` must remain positive. After changing context or priors, rebuild the prepared bundle so its fingerprint describes the reviewed experiment.

## 6. Native AnnData and initial states

```python
import anndata as ad
from virtualtissue_xenium.v2.inputs import from_anndata, prepare, write_json

adata = ad.read_h5ad("annotated_region.h5ad")
# Required reviewed obs columns: cell_type and lineage.
# counts layer holds unscaled nonnegative counts; obsm['spatial'] is in micrometres.
cells, genes = from_anndata(
    adata,
    cell_type_key="cell_type",
    lineage_key="lineage",
    sample_key="sample_id",
    layer="counts",
    spatial_key="spatial",
    coordinate_scale=1.0,
    max_cells=1000,
)

# Optional investigator-reviewed state assumptions, separate from measured RNA.
reviewed_effector_ids = set()
reviewed_exhausted_ids = set()
for cell in cells:
    cell["initial_state"] = {
        "activated": cell["cell_id"] in reviewed_effector_ids,
        "exhausted": cell["cell_id"] in reviewed_exhausted_ids,
    }

# context is a complete validated declaration loaded from the reviewed JSON.
bundle = prepare(cells, context, assayed_genes=genes)
write_json("prepared.json", bundle)
```

`initial_state` accepts only boolean `activated` and `exhausted`; both default to false. The AnnData adapter does not automatically convert marker scores or obs labels into these flags. Do not use Python truthiness to convert strings such as `"false"`. The adapter reads one expression row at a time and does not densify the full matrix. Select a manageable region first; the default maximum is 1,000 cells.

For single nuclei, use `assay: "snrna"` and `interaction_mode: "well_mixed"` when spatial positions are unavailable. Explicitly describe nuclei selection and populations depleted by preparation. Nuclear expression remains baseline evidence for a hypothetical whole-cell state; cell shape, spatial neighbors and cytoplasmic functional activity are not measured.

For mouse or another organism, supply a reviewed `gene_sets` dictionary to `prepare(..., gene_sets=...)`. The default human lists intentionally fail for other species. Case normalization does not provide ortholog mapping. A custom signature name is scored and reported, but only the known action-signature associations in `guidebook.py` alter weights automatically; use the established names after biological review or add a reviewed mapping in code.

## 7. Review before interpreting a run

Inspect the per-type menu, panel coverage and null scores; sampling declaration and missing populations; recognized target assumptions; initial simulated states; all default/user weights and applied modifiers; graph and population limits. The guidebook is an executable hypothesis specification. Validate its outcomes with held-out donors, observed cell counts, functional assays and paired perturbation data before calling its preferences calibrated or its output predictive.
