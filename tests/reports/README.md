# Validation reports

`release_0.1.2.json` records the 2026-10-26 checks for version 0.1.2: Python 3.12/3.13 suites, unchanged local biological trajectories, real-bundle prompt comparison, matched-run accounting, mocked notebooks and offline replay logic. It explicitly separates these checks from live service and browser validation.

The files below are the earlier 2026-10-25 records for package 0.1.0:

| File | Produced by | Shows |
|---|---|---|
| `reproduction_R03.json` | `vtx prepare`, the Part 1 notebook, the own-data route and a 24 h `vtx run --mode local`, on the public example | The package reproduces the R03 bundle of the original notebook byte for byte (39/39 files), and a 24 h simulation gives the same events as the original 24 h validation |
| `validate_r03_py313.json` | `tests/validate_r03.py` (original test) | Mocked-Jev integration on R03: budget stop before a round, deterministic resume, atomic failed round, interrupted checkpoint, unsafe and modified ZIPs rejected |
| `notebook_validation.json` | `tests/run_notebooks_mocked.py` | All three notebooks executed with a mocked Colab and Jev: Part 1 export of R03, mocked Jev run, Drive restore with no new requests, upload + local preview, explorer, no key in any output |

The final-state fingerprint of a trajectory depends on the CPU architecture (floating point). On Linux x86-64 the package and the original source give the same fingerprint (`ee2d642b…` for the 1 h R03 test). On macOS arm64 the original gave `3f7f63a6…`, with identical decisions and events.
