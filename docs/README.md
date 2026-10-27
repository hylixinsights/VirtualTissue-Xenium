# VirtualTissue-Xenium documentation

Start with the route that matches your goal:

* **Public example in Colab:** follow the [README quick start](../README.md#2-quick-start-in-google-colab); use [05](05_google_colab.md) for Drive and runtime details.
* **Your own Xenium sample:** read [06](06_your_own_xenium_data.md) first to check suitability, prepare the files and review your cluster map. It gives Colab, local Jupyter and command-line routes.
* **Your own computer:** install with [08](08_local_installation.md), then choose the public example or your own data. Replace notebook Drive paths with local paths before running input cells.
* **Saved results only:** go directly to [04](04_part4_results.md). No API key is needed.

Guides 01–04 explain the four stages in detail. The other guides are references for setup, data, keys and troubleshooting.

| # | Guide | Read it when you want to… |
|---|---|---|
| 01 | [Part 1 · Prepare a region from Xenium data](01_part1_prepare_region.md) | go from a 10x Xenium download to a simulation-ready bundle ZIP and choose the region (R01, R03, …) |
| 02 | [Part 2 · The input bundle](02_part2_input_bundle.md) | understand exactly what the simulator reads, check a bundle, or share one |
| 03 | [Part 3 · Run the simulation](03_part3_run_simulation.md) | choose a scenario, duration, seed and API limits; understand rounds, the journal, stop reasons and resuming |
| 04 | [Part 4 · Results and the replay viewer](04_part4_results.md) | use the replay, read the output files, plot, follow one cell, compare runs |
| 05 | [Google Colab, step by step](05_google_colab.md) | know every path to edit, where the key goes, which runtime to pick, how Drive is used |
| 06 | [Your own Xenium data](06_your_own_xenium_data.md) | simulate your own sample: required files, the cluster annotation map, cell-type menus |
| 07 | [The Jev API key](07_jev_api_key.md) | get, store and protect the key; understand what is sent and what it costs |
| 08 | [Local installation and the `vtx` command line](08_local_installation.md) | run everything on your own computer or a server; run the tests |
| 09 | [Architecture and maintenance](09_architecture.md) | understand why the engine is vendored unchanged; update code or notebooks safely |
| 10 | [Troubleshooting](10_troubleshooting.md) | fix an error message |

## The two files you pass between parts

```text
Part 1 ──► simulation_input_<REGION>_<HHMMSS_ffffff>.zip ──► Parts 2–3 ──► <RUN_NAME>_results.zip ──► Part 4
           (≈ 20 MB, ≤ 1,000 cells, verified by checksums)                 (replay.html + tables + journal)
```

## Glossary

| Term | Meaning |
|---|---|
| **Bundle** | The ZIP exported by Part 1 for one region. It is the only input of the simulator. |
| **Acting cell / actor** (`role = roi`) | A cell with a supported runtime label that passed QC. It asks Jev for decisions. |
| **Context** (`occluded`, `halo`) | Cells that take up space and can be pushed, but never decide. `occluded` cells sit inside the region but lack a supported identity or failed QC. `halo` cells are the 30 µm ring around the region. |
| **Runtime label** | The cell type used by the simulator: `Tumor`, `CAF`, `Fibroblast`, `Endothelial`, `Macrophage`, `T_CD8_cytotoxic`, `T_CD4`, `B_cell`, `Plasma_cell`. |
| **AUCell percentile** | Rank of a cell's AUCell score within its provisional type and RNA-depth decile, using the whole section as the reference (0 = lowest, 1 = highest). |
| **Round** | All cells that need a decision in the same 10-minute tick. A round is sent as one or more HTTP batches and is validated as a whole before it is applied. |
| **Batch / HTTP attempt** | One POST to the Jev API with up to 20 cell questions (≤ 40 KB). `MAX_HTTP_ATTEMPTS` counts these. |
| **Journal** | `jev_journal.jsonl`: an append-only, hash-chained record of every request intent and response. Resuming reuses responses preserved in this file; responses lost after the last backup may be requested again. |
| **Local preview** | The same simulation with programmed local rules instead of Jev. It is free and labelled as such everywhere. |
| **Scenario** | `baseline_unknown`, `declared_recognition` or `recognition_off`: declared assumptions about antigen recognition and effector competence. |

[11 · Experiment design](11_experiment_design.md): useful control/intervention comparisons, proposed perturbation controls, measured versus assumed baseline, and the requirements for non-cancer tissue and other platforms.

[12 · Jev feedback and 0.1.3](12_jev_feedback_013.md): declared tissue, rounding compatibility, t0 diagnostics and prompt validation.
