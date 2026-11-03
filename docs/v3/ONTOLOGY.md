# Versioned lineages, rates and signatures

The schema 3.0 guidebook reads packaged data from
`src/virtualtissue_xenium/v3/data/`. The schema 2.0 path retains its historical
constants, normalization and rate-table fingerprint. Its original metadata
formula is a historical discrepancy; it must not be copied into v3 Methods.
The v3 table records the implemented first-event competing-risk formula:
`p_a = (1 - exp(-K*dt)) * k_a/K`, `REST = exp(-K*dt)`, after gates and multipliers.

## Permissions and kinetics

`ontology/lineages_v3.json` contains 22 lineages with explicit parents,
allowed actions, motility, radius, target/killing/clearance permissions and
secretion kind. CD8 proliferation explicitly requires the activated state;
this flag is data, not an inferred property of all lymphocytes. `unknown` remains REST-only. NK has killing permission subject to
declared recognition; T_reg secretion is marked suppressive. The loader rejects
inconsistent permissions, absent rate rows, unsupported actions, parent cycles,
nonfinite rates and missing assumption/source metadata.

`ontology/kinetics_v3.json` contains per-hour hazards and β declarations.
**All default numerical rates are classified as assumed for the target tissue.**
Public studies can motivate a range but do not calibrate an individual human
Xenium ROI. All eight bibliographic identities now have verified DOI/source
links; quantitative rate extraction and biological review remain pending.
Two examples of corrected interpretations are:

- Halle et al. report **2–16 infected targets per CTL per day** in mouse viral
  infection. The earlier draft's “2–16 hours per target” was incorrect.
  [Primary article, DOI 10.1016/j.immuni.2016.01.010](https://pmc.ncbi.nlm.nih.gov/articles/PMC4846978/).
- Breart et al. estimate approximately one killed target per six hours of
  CTL–tumour interaction in a mouse adoptive-transfer experiment.
  [Primary article, DOI 10.1172/JCI34388](https://www.jci.org/articles/view/34388).

The rate file uses `bibliography_verified: true` for reference identity and
`numerical_rate_verification: pending` plus `verify: true` for the unresolved
quantitative review. See [RATE_REVIEW.md](RATE_REVIEW.md) for all eight sources,
access limitations and removed unsupported numerical descriptions. No numerical
target-tissue rate is labelled `measured`. The previous proposed activated-CD8 proliferation rate is not applied
implicitly; its selection remains a modelling decision. `MOVE` has an explicit
assumed event hazard; displacement speed/persistence are separate dynamics
parameters.

`ontology/lineages_v2_legacy.json` exports the original constants exactly for
traceability. `datasets/manifests.json` snapshots the existing dataset registry,
and its hash is included with `DATA_VERSION` in every v3 guidebook. The v2
download API and its original registry remain unchanged.
`profiles/skin_melanoma.json` and `profiles/brain_cortex.json`
declare expected lineages and support explicit rate overrides. They currently
contain no fitted overrides. A profile/organism mismatch is rejected and
unexpected lineages appear in guidebook provenance. Disease/treatment text never
selects a profile automatically.

Platform files declare initial QC, scoring-universe and geometric defaults.
Their presence does not claim that all platform file formats or binned kernels
have been validated. These defaults remain assumptions requiring sensitivity
analysis.

An explicit context `neighbor_radius_um` overrides the platform default. For
cell-agent geometry, a lineage's declared `mean_radius_um` takes precedence over
`dynamics.default_radius_um`; the latter defaults to the platform's
`mean_radius_um` and can be overridden explicitly in the context. Every current
lineage has its own radius, so the platform value is a fallback, not a global
rescaling of those lineages. Neither radius defines a binned measurement:
`measurement.unit_size_um` records that size separately, and the fixed lattice
does not use daughter placement or cell-radius repulsion.

## Signature membership and species

`signatures/human_v3.json` has 14 sets. Three are versioned Hallmark gene lists:
IFN-γ response, TGF-β signalling and hypoxia. Membership was read from the
[Broad Institute ssGSEA2.0 repository](https://github.com/broadinstitute/ssGSEA2.0/blob/master/db/msigdb/h.all.v7.0.symbols.gmt)
and pinned by the downloaded file's SHA-256 and retrieval date. The source
filename is `h.all.v7.0.symbols.gmt`; the exact content hash, rather than a moving
branch name, identifies the snapshot. Cite
[Liberzon et al., DOI 10.1016/j.cels.2015.12.004](https://pubmed.ncbi.nlm.nih.gov/26771021/)
and retain the MSigDB attribution and
[CC BY 4.0 terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).

The other sets are explicitly investigator-selected marker heuristics, including
the seven original v2 sets and MHC-II, Treg, NK and tissue-resident-memory sets.
Their metadata distinguishes these selections from complete published database
sets. CellMarker and tissue-residency references provide context; they do not
validate the particular small marker combination or its use as an action score.
Independent biological review is still required.

Mouse membership is stored explicitly in `signatures/mouse_v3.json`. The file
`human_mouse_orthologs_v3.json` contains **432 reciprocal one-to-one symbol pairs**
extracted from the complete
[MGI human/mouse homology download](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt),
with source hash, MGI and Entrez identifiers. Candidate mappings are aggregated
across homology classes before testing uniqueness; 55 required human symbols
are excluded because a reciprocal one-to-one mapping is absent. Mouse MHC-II
uses a separate explicit marker selection with exact symbols checked in MGI;
it is not an automatic HLA conversion. The computational mapping check does not
replace independent review of biological equivalence.

The runtime performs **exact, case-preserved matching** of expression, assayed
panel and species-specific signatures. No runtime uppercasing or orthology
inference occurs for schema 3.0. Missing genes reduce declared panel coverage;
insufficient overlap remains null and produces a coverage warning, rather than
being interpreted as low biological activity. Mouse signatures may lose enough
genes through conservative mapping to be unavailable under the overlap gate.

The guidebook records data hashes, species signature provenance, platform and
tissue profiles, coverage warnings, reference fallback counts and lineage rules.
These identify the assumptions used in each prepared bundle; they are not a
claim that the assumptions have been experimentally validated.
