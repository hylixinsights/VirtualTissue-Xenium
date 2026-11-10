# Sample context and biological restrictions

Collect the sample description before preparing an object:

```bash
vtx5 context-init my_context.json --template examples/v3/context_user_object.json
vtx5 prepare object.json my_context.json prepared.json --lineage-map lineage_map.csv
```

For the LUAD example, use `examples/v5/luad_tls/context.json` as the template and retain its declared evidence-policy argument when preparing measured scores. The questionnaire makes no API request and refuses to overwrite an existing output. Existing `prepare` commands with a context JSON remain noninteractive.

The English questions ask:

- Which organism is this sample from?
- Which tissue or organ was sampled?
- What disease or pathology is present (use `none` for healthy tissue)?
- How was the sample or region selected, and what material does it represent?
- Which anatomical region or compartment was sampled?
- What is the age or developmental stage of the sample?
- Which cell subtypes, injury, regeneration or tumour states affect interpretation?
- Which secretory cell populations or products are known, for example mucus-producing cells?
- Which additional lineage-specific actions should be disabled, and why?

Press Enter to retain a template value. Optional fields without a value can remain absent; use `unknown` when a fact is unavailable. Review retained values when changing the sample. Treatment declarations (`perturbations`), sampling design and missing populations, recognition, evidence policy and numerical settings come from the template. Edit those declarations explicitly when they differ. The questionnaire does not infer them from a disease name.

## How context reaches the rulebook

Required `organism`, `tissue`, `disease` and `sampling` declarations already existed. The optional `sample_context` object adds four validated text fields:

```json
{
  "anatomical_region": "adult cerebral cortex",
  "developmental_stage": "adult",
  "cell_state_notes": "Mature neurons and resident glia; no tumour declared.",
  "secretion_notes": "No mucus-producing population declared."
}
```

Place this object under `sample_context` in the full context JSON. These are investigator declarations, not RNA measurements or demonstrated functional activity. They are retained in the prepared bundle and the saved run's `guidebook.context`, bound by the bundle fingerprint in the run configuration and the full-context hash in Jev requests, and sent to Jev with that distinction stated. No retraining is needed. Existing contexts without the optional field retain the v4 request bytes and local behavior.

The local model uses the lineage ontology to define capabilities. Expression evidence can modulate the hazard of an already-permitted action. Descriptive context alone does not change that hazard. An explicit `context_modifiers` rule can further restrict an action; Jev receives only the remaining menu and cannot restore a removed action. Context-exchange controls replace or remove the original sample details when using a donor context, avoiding contradictory sample descriptions.

## Brain samples and division

The ontology already excludes `PROLIFERATE` for `neuron` and `oligodendrocyte`. It permits division for some other populations, including microglia, astrocytes and malignant cells. A brain sample therefore must be interpreted by population and condition. Adult human neocortical neuron birth dating and microglial turnover studies support distinguishing these populations; they do not calibrate the simulator's rates. See [Bhardwaj et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1567918/) and [Réu et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28746864/).

Do not represent a neural stem/progenitor population by enabling mature-neuron division. Progenitor proliferation and differentiation require a separately reviewed representation; the current model has no general neurogenesis mechanism.

For a reviewed sample assumption that further disables a permitted action, enter a pair such as `astrocyte:PROLIFERATE` in the questionnaire and provide the sample-specific reason. This is a modelling restriction, not a default biological claim about astrocytes. The equivalent declaration in the full context is:

```json
"context_modifiers": [
  {
    "action": "PROLIFERATE",
    "lineages": ["astrocyte"],
    "multiplier": 0,
    "reason": "Investigator-declared restriction for this experiment; record its evidence and scope here."
  }
]
```

The rule is sample-wide for that lineage, applies after expression modulation, and appears in each cell's weight-adjustment audit. It cannot add a capability or remove the `REST` fallback. Current restrictions do not select an individual cell subtype within a shared lineage; use an appropriate reviewed mapping or extend the ontology through a separate reviewed change when that distinction is necessary.

## Mucus and secretion

`SECRETE` emits an abstract simulated inflammatory or suppressive signal. It does not simulate mucin synthesis, mucus release, mucus accumulation or clearance. Mucus-related notes never turn it into a mucus action. Disabling generic `SECRETE` to express “no mucus” would also remove unrelated signalling.

Actual mucus production depends on anatomical region and competent secretory populations, not just an organ label. A primary human airway study illustrates the location and cell-population specificity of MUC5AC/MUC5B expression: [Okuda et al., 2019](https://doi.org/10.1164/rccm.201804-0734OC). A future mucus mechanism would need its own action semantics, subtype permissions, measured support and validation. The current implementation records that context and explicitly reports mucus as unsupported.

The seven LUAD GO20 evidence rules remain exploratory assumptions requiring biological review. Context collection does not approve their directions, the wildcard proliferation proxy, recognition, annotation mappings or assumed kinetic rates.
