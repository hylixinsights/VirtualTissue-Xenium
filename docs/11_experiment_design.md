# Experiments: from a measured starting tissue to a useful comparison

The user should be able to ask **what changes in this sample when I change one mechanism?** The main result is a matched control/intervention comparison, with a spatial explanation of where and when the trajectories diverge. More animation or more deaths is not the success criterion.

## Available now

Use the same bundle, model, seed, horizon and population ceiling under `declared_recognition` and `recognition_off`. The first assumes tumour recognition and effector/clearance competence; the second changes only recognition. This is a mechanism experiment, not a drug experiment. `baseline_unknown` leaves both recognition and competence unknown and answers a different question. Neither enabling recognition nor increasing the horizon guarantees killing.

Inspect active CD8 counts, contacts, excluded actions and growth headroom before running. Use `vtx compare` to validate matching and restrict outcomes to the latest common saved frame. Follow tumour survival, cytotoxic hits, apoptosis, hypoxic necrosis, division, CD8 exhaustion and clearance; keep death mechanisms separate in exported analyses. A reduced tumour count may result from less growth, more death, or both.

The replay now provides visible field buttons, automatic signal selection, scales, empty-field explanations, arrow callouts and clickable event cards. They describe saved events and local field maxima, not statistical anomalies or inferred causality. Replay plays a saved trajectory; a new intervention requires a new run.

## Proposed perturbation workbench — not implemented

1. **Starting tissue.** Show composition, measured t0 pathway scores, annotation confidence, missing signatures and spatial contacts. Separate measured evidence from assumptions about perfusion, recognition, rates and initial functional state. AUCell percentiles rank cells within type and detection-depth group; they are not absolute pathway activities or directly comparable concentration scales across samples. Daughters currently inherit t0 scores. A future dynamic pathway state needs an explicit update model.
2. **Experiment question.** Select a mechanism and retain a visible untreated/control arm. “Baseline” here means the control trajectory, not necessarily the scenario named `baseline_unknown`.
3. **Intervention.** Each control must specify what it changes, its baseline value, intensity, start time, duration and supported cell types. Distinguish an assumed relative effect from an experimentally calibrated dose. Save the intervention schedule and the effective parameters with the run.
4. **Run experiment.** Preflight checks cell types, prerequisites, model compatibility and request allowance before execution. Use matching initial states and seeds in both arms; repeat across seeds for robustness. Distinguish biological stochasticity from provider nondeterminism. Never animate an effect before the simulation produces it.
5. **Compare.** Synchronised control/intervention tissue views, common field scales and time cursor, outcome curves and differences from control. Highlight the first recorded divergence and allow inspection of its local evidence and decision history. A causal claim requires the controlled comparison and sensitivity checks; an arrow alone is insufficient.

| Desired control | Existing basis | Work still needed before offering it |
|---|---|---|
| Enable CD8 tumour recognition | Current declared-recognition scenarios | Available now; present the assumption and matched control clearly |
| Reduce hypoxia | Oxygen diffusion, endothelial perfusion and consumption | Define whether the intervention changes perfusion, oxygen supply or consumption; add bounded controls and sensitivity checks |
| Increase angiogenesis | VEGF, sprouting and delayed perfusion mechanisms in the kernel | Expose a validated scenario; vessel formation and functional perfusion must remain separate outcomes |
| Reduce PD-L1-mediated inhibition | A declared checkpoint-blockade mechanism in the kernel | Expose and test the mechanism as an assumed inhibitory effect; a named drug requires pharmacology and calibration |
| Administer drug X | No generic drug-response module | Target availability, exposure and dose/time response, target engagement, affected cell states, damage/death rules and calibration data |
| Observe acquired resistance | No acquired-resistance module | Persistent sensitive/resistant states, lineage and selection/adaptation rules, target changes and longitudinal validation. Regrowth alone does not establish resistance |
| Increase neutrophil migration | No neutrophil agent/menu | Identity mapping, recruitment source, chemotaxis, behaviours and fate rules |
| Highlight TNF-α release | No TNF field/source in this model | Implement and validate secretion, transport and response; do not label another signal as TNF |

The first perturbation milestone should expose a small set of already implemented mechanisms with matched comparisons. Named drugs and resistance follow only when their state transitions and parameter evidence exist. The interface should show unsupported mechanisms as planned capabilities, not functioning sliders.

## Healthy tissue and other technologies are separate requirements

Changing all tumour annotations to `epithelial_other` in the public cancer dataset is a negative test of the tumour-interface preparation contract. It does not create a healthy sample or validate non-cancer biology. The CLI and notebook stop before AUCell when no QC-passing tumour candidates remain. Copying the annotation JSON alone does not change its labels; retain the edited map to reproduce the test. A manual centre is not a bypass.

Non-cancer support needs a tissue/task profile with appropriate required compartments, region selection, pathway coverage and action menus. The current tumour-composition requirement cannot simply be removed and presented as a validated healthy-tissue model.

Other platforms need input adapters. Today `analysis.tar.gz` supplies cluster assignments and differential-expression tables, and the loader also requires the Xenium HDF5/parquet schemas. Removing that filename check would not make the pipeline portable. Define a platform-neutral input contract: sparse cell-by-gene counts and gene identifiers; stable cell IDs; centroids and polygons with declared physical units; reviewed labels or cluster assignments; measured panel and QC; provenance and missing-data indicators. Differential expression can be supplied or computed with a declared method rather than requiring a vendor archive. Preserve the existing Xenium adapter and validate each new adapter's joins, units, geometry and scoring coverage. Spot-based assays require a separate representation; spots should not silently become single cells.

## Release boundaries

Version 0.1.2 reports the cumulative population ceiling without changing its policy, distinguishes deaths and blocks, provides all ten fields and shared-scale replay copies, and removes unnecessary scenario information from Jev prompts. These changes do not validate Jev decisions biologically. Named drugs and intervention schedules remain outside this release. Persistent AUCell cache reuse is also deferred; geometry/region feasibility is checked before scoring, and additional regions can still be exported from the same in-memory notebook preparation.
