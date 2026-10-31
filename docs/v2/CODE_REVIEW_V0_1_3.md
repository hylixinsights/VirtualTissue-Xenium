# Review of VirtualTissue-Xenium v0.1.3

Reviewed release: [`v0.1.3`](https://github.com/hylixinsights/VirtualTissue-Xenium/tree/v0.1.3), commit `16d0a307d118a42143a7034c81b5dc366ca2d62c`. Review date: 2026-10-30. References below point to that immutable release. This is a source review, not a claim of independent biological validation.

## Assessment

The architecture makes sense as an exploratory agent-based mechanism simulator: measured cells initialize a tissue, a contract limits their actions, Jev ranks those actions, and a deterministic kernel enforces consequences. The code already takes provenance, replay, API accounting, and the distinction between measurements and assumptions seriously. Those are useful foundations.

The current implementation is a lung tumour-interface prototype. It cannot yet accept a general annotated single-cell object and produce a biologically validated simulator for arbitrary tissues. The largest gap is the biological contract and input model, rather than the choice API itself. Version 2 should preserve the separation of responsibilities while making tissue context, sampling, action permissions, evidence, initial state, and uncertainty explicit.

| Priority | Finding | Why it matters | Version 2 response |
|---|---|---|---|
| P1 | Preparation and actions are restricted to a tumour interface | Healthy brain, skin, nasal mucosa and many nuclei datasets fail or acquire inappropriate biology | Separate generic annotated-object intake from tissue-specific action profiles |
| P1 | Jev scores are executed through a maximum-score choice | They do not constitute sampled or calibrated biological transition probabilities | Record score semantics and the explicit selection rule |
| P1 | Tissue context is a short free-text declaration | Disease, sampling, dose, timing and target specificity have no structured representation | Validate a versioned context manifest |
| P1 | Initial simulated molecular state overrides measured disease state | Chronically exhausted cells initially become simulated effectors | Require explicit mapping of measured evidence into declared initial state |
| P1 | Expression evidence stays frozen | Simulated phenotypes are not predictions of future RNA expression | Maintain separate immutable expression evidence and mutable simulated state |
| P2 | Population ceiling counts every body ever allocated | Turnover eventually stalls even if living space becomes available | Distinguish living capacity, allocated-body safety limits, and biological crowding |
| P2 | Action repertoires and physics are mechanism-specific assumptions | A name such as anti-PD1 cannot establish recognition or response | Use reviewed mechanisms and assumptions, with controls and ablations |
| P2 | AUCell reference construction limits transferability | Scores are panel- and reference-dependent; small populations can have no relative scores | Report coverage, reference sizes and missingness; never equate percentiles with transition rates |
| P2 | AUCell provenance hard-codes a 5,001-gene universe | Own-panel audit metadata can describe the wrong panel | Derive the recorded universe from the actual input |
| P2 | Existing comparison tooling supports a narrow pair of scenarios | General treatment, mutation and multi-seed experiments need another comparison contract | Validate matched inputs and changed assumptions across arbitrary declared arms |

P1 means a prerequisite for the requested general workflow; P2 means a material limitation or defect to address or expose. These are not claims of clinical safety classifications.

## Findings with code evidence

### 1. The generic-input path does not exist

[`own_data.py:40`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/own_data.py#L40) requires QC-passing tumour-candidate cells. [`xenium_prepare.py:350`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_prepare.py#L350) requires minimum tumour, fibroblast and immune fractions for a proposed region. The CLI invokes this scope check before scoring and exporting. This is an intentional scope restriction, but it blocks the user's new requirements.

The accepted runtime labels are fixed in [`own_data.py:71`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/own_data.py#L71). [`menus.py:32`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/menus.py#L32) gives B and plasma cells no decision menu; other handlers cover selected immune, tumour, fibroblast and endothelial behavior. Renaming neurons or glia to one of these labels would misrepresent biology. A general profile registry must be explicit about which cell labels are modeled and which remain context.

### 2. Choice probabilities are not biological transition probabilities

[`engine.py:17`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/engine.py#L17) requires the selected option to have a maximum probability. [`autonomous.py:291`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/autonomous.py#L291) executes that choice. The distribution is audited but is not sampled. Repeated identical cell states can therefore make identical decisions despite nonzero scores for alternatives. Changing the physics seed does not turn those scores into transition-frequency estimates.

[`jev_provider.py:52`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/jev_provider.py#L52) permits the reported choice to trail the rounded maximum by 0.01 and creates a tied kernel copy, while retaining the original distribution. This compatibility repair is transparent and tested, but future analysis must use the original model scores when studying Jev. Version 2 should distinguish reported Jev scores, guidebook prior weights, combined operational weights, and the executed action. A stochastic selection mode is useful for exploratory ensembles, provided it is labeled as a modeling assumption. Calibration requires independent observed outcomes.

### 3. Free text cannot serve as the complete experiment contract

[`jev_payload.py:20`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/jev_payload.py#L20) accepts a single-line declaration of at most 300 characters. [`simulation_workflow.py:128`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/simulation_workflow.py#L128) chooses none/custom/bundle context. The public run settings support only three recognition scenarios at [`simulation_workflow.py:106`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/simulation_workflow.py#L106).

There is no typed experiment representation for tissue origin, disease, specimen preparation, nuclei versus cells, enrichment, whole-section representativeness, perturbation target, dose, schedule, mutation, or uncertainty. Those factors should be fields in a manifest, including an explicit unknown value when unavailable. Free text can supplement that manifest; it should not silently set numerical rates or unlock an unimplemented action.

### 4. All cells start with the same simulated baseline for several key variables

[`autonomous.py:79`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/autonomous.py#L79) initializes exhaustion to zero, exhaustion state to effector, and PDL1 to zero. The expression signatures enter the decision observation, but they do not initialize these particular state variables. A chronic-disease input may consequently start with simulated effector CD8 cells regardless of the input's exhausted-cell annotations.

Version 2 needs explicit state initialization rules with provenance. An exhausted annotation, an RNA signature, measured protein, and an assumed numerical exhaustion value are different evidence types. If no validated mapping exists, expose a configurable initial-state prior and test sensitivity instead of presenting a single inferred value as measured.

### 5. Future expression is not simulated

The prompt semantics in [`jev_payload.py:9`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/jev_payload.py#L9) correctly state that t0 signatures stay frozen. [`autonomous.py:518`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/autonomous.py#L518) inherits the root observation into daughters. The mutable state contains imposed variables such as damage, exhaustion, polarization, glycolysis and PDL1.

This is an honest and useful scope boundary. However, a simulated state transition must not be described as a predicted RNA trajectory or a learned post-treatment single-cell profile. An eventual expression-transition model needs an independently trained and validated module, explicit uncertainty and held-out perturbation data.

### 6. Resource capacity and biological crowding are entangled

[`autonomous.py:543`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/autonomous.py#L543) tests `len(t.ids)` against `max_population`. Removed cells retain allocated IDs. Recruitment uses the same limit at line 577. Thus cell death does not release an allocation slot. [`recorded_simulation.py:63`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/recorded_simulation.py#L63) documents the rule, and reporting separates several blocking causes.

For turnover experiments, report current living/present population, historical allocated bodies and the relevant limit separately. Never interpret a ceiling-driven growth plateau as treatment efficacy. Use a living capacity only if that is the declared policy; retain an independent historical allocation limit to protect memory.

### 7. Biology is a declared mechanism, not a consequence established by a prompt

Cytotoxic killing requires declared recognition, effector competence, material and contact ([`menus.py:50`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/menus.py#L50)). The checkpoint effect enters an assumed damage equation ([`autonomous.py:408`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/autonomous.py#L408)). Exhaustion is a scalar accumulation-and-threshold model ([`autonomous.py:476`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/autonomous.py#L476)). Fields are dimensionless quasi-steady solutions on a reflecting rectangle ([`fields.py:1`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/fields.py#L1)).

These choices make a mechanism experiment runnable. They do not establish antigen specificity, receptor occupancy, actual cytokine concentration, validated exhaustion differentiation, or a clinically meaningful response duration. Disease should change reviewed action permissions and operational priors only through declared rules. A cancer diagnosis must not automatically make every epithelial cell malignant or every immune cell tumour-reactive. Anti-PD1 must not force efficacy; chronic disease must not force exhaustion.

### 8. Relative signature scores need coverage and reference audits

[`xenium_prepare.py:198`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_prepare.py#L198) computes panel-universe AUCell with a minimum overlap and detected-gene cutoff. Within-type depth-conditioned percentiles require at least 30 cells per group ([line 255](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_prepare.py#L255)); constant signatures and insufficient groups stay missing. This conservative missingness treatment is a strength.

Small types can therefore have no depth-conditioned scores, and the panel/reference differs between Xenium, scRNA-seq and snRNA-seq. A percentile of 0.8 means relative rank in that reference, not 80% probability of a behavior. Enrichment changes the observed cell mixture; it does not recover a whole-tissue mixture. A nuclei preparation also cannot be treated as full-cell expression without explicit assay assumptions.

There is one concrete metadata defect: [`xenium_prepare.py:270`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_prepare.py#L270) always records a universe of 5,001 assayed features, even though the own-data path can load a different panel. Derive this number and description from `sample["genes"]`.

### 9. Pairwise comparison is intentionally narrow

[`results.py:255`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/results.py#L255) compares identical scenarios or the declared-recognition/recognition-off pair. It requires matched provenance and a common saved time. Different seeds need separate analysis. These checks protect the existing demonstration, but they do not implement general treatment comparisons or uncertainty estimates.

Version 2 needs a comparison manifest identifying the one intended intervention difference, shared object/guidebook/model versions, matched seed sets, censoring and budgets. Counterfactual arms should share initial state and exogenous randomness as far as possible. Run completion must be checked before comparing endpoints.

## Strengths to preserve

- **Restricted decisions.** The provider receives permitted local actions, while the kernel retains numerical rates and physical consequences.
- **Synchronous snapshots.** [`autonomous.py:299`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/xenium_runtime/v2/autonomous.py#L299) validates the full round before state writes and uses seeded conflict reservations. Remote batching need not mean sequential biological time.
- **Fail-closed paid inference.** Model/version/option validation and durable request intents prevent silent substitution and ambiguous free retries ([`jev_provider.py`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/jev_provider.py)). Keys are omitted from persisted headers; response fields are restricted and redacted.
- **Input integrity.** [`simulation_workflow.py:43`](https://github.com/hylixinsights/VirtualTissue-Xenium/blob/v0.1.3/src/virtualtissue_xenium/engine/simulation_workflow.py#L43) validates archive paths, extraction limits, symlinks, manifests and hashes, and does not execute code bundled with data.
- **Reproducible resumption.** Run configuration, source hashes, environment, cached responses, and complete checkpoints are recorded. A seed alone does not reproduce an uncached remote model response, so retaining the response journal matters.
- **Honest readouts.** Computational ceilings, missing signature evidence, mocked transport, local previews and simulated units are labeled.

## Biological evidence informing the redesign

Anti-PD1 response depends on cell state and context. Experimental work distinguishes exhausted T-cell subsets with different proliferative and checkpoint-response properties; it does not justify converting every exhausted cell into an effector. See [Miller et al., Nature Immunology, 2019](https://www.nature.com/articles/s41590-019-0312-6).

Human melanoma single-cell data also show response-associated T-cell states, useful as independent evidence rather than a mandatory output. The corresponding [GSE120575 study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575) profiled CD45-selected cells, illustrating why enrichment must be part of the experiment contract. See [Sade-Feldman et al., Cell, 2018](https://pubmed.ncbi.nlm.nih.gov/30388456/).

Xenium measurements require panel, segmentation and transcript-contamination checks. These affect the certainty of an expression-derived guidebook; they should not silently become impossible actions. See [Resolving sensitivity, specificity and signal contamination in Xenium spatial transcriptomics, Nature Methods, 2026](https://www.nature.com/articles/s41592-026-03089-8).

## Release interpretation

A working v2 software release can support general input, explicit context, inspectable guidebooks and reproducible hypothetical interactions. It should remain labeled an exploratory simulator until prospective or held-out biological validation supports stronger claims. A live API smoke test establishes interface compatibility only.
