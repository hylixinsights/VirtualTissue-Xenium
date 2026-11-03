# Mac launcher reference

Use **START_TEST.command** next to the top-level README for the packaged test.
This internal reference describes individual launchers for advanced use. It is
not a second setup path you need to follow.

Keep the repository together. Launchers use this checkout's `src/`, even if
another version is installed, and support paths containing spaces. Python
**3.12 or 3.13** is required. Apple Silicon has been tested; Intel compatibility
is expected but has not been verified here. No Homebrew installation is automatic.

The following commands run from the repository root. Every launcher accepts
`--help`. If a copy loses executable permissions, prefix a command with `zsh`
or run `chmod +x scripts/mac/*.command`.

## Environment

```sh
scripts/mac/00_install.command
scripts/mac/03_check.command
```

Basic installation creates `.venv-mac` without network access. Launchers load
`src/` directly. Existing valid environments are reused; no directory is deleted.
To choose Python, prefix the installer with `VTX_PYTHON="/path/to/python3.12"`.
When exported, this variable takes priority over `.venv-mac`; use
`unset VTX_PYTHON` afterward if you want the launchers to use `.venv-mac`.
Manual environment activation is unnecessary.

Optional extras can download packages:

```sh
scripts/mac/00_install.command --with-h5ad
scripts/mac/00_install.command --with-tests --with-plots
scripts/mac/03_check.command --tests
```

`--with-tests` installs the `test`, `xenium`, `spatial`, and `analysis` extras.
The check runs v2/v3, Mac, and exporter tests from the correct directory. Legacy
engine tests require their separate pinned environment. Node is needed for
JavaScript replay/catalog tests, and zsh for launcher tests. Missing tools or
Python dependencies can cause skipped tests, so use the top-level release test
for the final verdict rather than relying on pytest's exit code alone.

`--offline` passes `--no-index` to pip; extras must already be available locally.
`--dry-run` only prints the configuration. The environment check verifies the
imported source, 35 engine files, 22 lineages, and 14 signatures per species.
These are engineering checks, not biological validation.

## Synthetic demonstration

```sh
scripts/mac/01_demo.command "outputs/mac/my demo"
open "outputs/mac/my demo/catalog.html"
```

The example has 96 generated cells with artificial labels, counts, and outlines.
One declared invalid outline and one low-count cell exercise QC. These are not
experimental data. No ROI is selected automatically. Choose a region in the
catalog, inspect its composition, and save the selection JSON from the browser.

## Existing objects and labels

The catalog accepts raw JSON (`cells`, panel, and metadata) or a prepared v3
bundle. Supplied polygons are preserved with QC; invalid outlines use centroid
fallback. No replacement contour or cell identity is inferred. Histology and
independent annotation review are not required to run this workflow.

```sh
scripts/mac/vtx.command roi-catalog \
  "/path/to/object.json" "/path/to/new_catalog.html" \
  --cells-per-roi 100 --grid-um 100 --max-candidates 30
```

If labels require mapping to model lineages, add
`--lineage-map "/path/to/map.csv"` with CSV columns `cell_type,lineage`. Use the
same map for catalog creation and selection.

For H5AD, install `--with-h5ad` and use the raw-object exporter:

```sh
scripts/mac/05_export_object.command --help
scripts/mac/05_export_object.command \
  "/path/to/object.h5ad" "/path/to/new_object.json" \
  --cell-type-key cell_type --lineage-key lineage \
  --polygon-key outlines --coordinate-scale 1 --max-cells 1000
```

Export preserves existing labels and does not fit scores or select an ROI.
Coordinate conversion and polygon provenance must match your object. If counts
are in a `counts` layer, add `--layer counts`; otherwise the exporter reads `X`.
Positive finite values alone do not establish that a matrix is unnormalized.
For an H5AD already restricted to the intended region, direct preparation is
available through `scripts/mac/vtx.command prepare --help`.

## Chosen ROI and local simulation

```sh
scripts/mac/02_run_roi.command \
  "outputs/mac/my demo/object.json" \
  "$HOME/Downloads/my_selection.json" \
  "outputs/mac/my demo/context.json" \
  "outputs/mac/my simulation" --steps 12 --seed 21
open "outputs/mac/my simulation/run_local/replay.html"
```

Add the same `--lineage-map` used for the catalog when applicable. The launcher
verifies selection against the source, recomputes ROI scores, saves
`prepared.json`, `inspect.json`, and `expression_audit.json`, and runs **only the
local model**. Each step's duration comes from the context's `time_step_hours`.
Use a context describing your material; the demo context declares synthetic data.

Copy the **entire `run_local` directory** when sharing replay, keeping its HTML
and JavaScript chunks together. Movement rigidly translates polygons; daughter
shapes are synthetic copies. These mechanisms are not validated cell mechanics
or proven temporal predictions. A new output directory is required. Failed runs
retain diagnostic files; no directory is erased or resumed automatically.

```sh
scripts/mac/03_check.command "outputs/mac/my simulation/prepared.json"
```

## Explicit Jev use

This estimate requires no key and **does not call the API**:

```sh
scripts/mac/04_jev.command estimate \
  "outputs/mac/my simulation/prepared.json" "outputs/mac/estimate" \
  --steps 2 --max-http-attempts 20
```

Read `cost.json`. The estimate starts from t0; population changes, failures, and
retries can change usage. `--max-http-attempts` limits HTTP attempts, **not money**.
The script does not assume service prices.

Only the following explicit `run` command makes paid calls. It requires
`TYPESAFE_API_KEY` or `JEV_API_KEY` in the environment; never put keys in arguments
or project files:

```sh
scripts/mac/04_jev.command run \
  "outputs/mac/my simulation/prepared.json" "outputs/mac/short_jev_run" \
  --steps 2 --max-http-attempts 20
```

`run` requires an explicit attempt limit. Exhausting the limit returns an error
and preserves records; no automatic increase or resume occurs. The default is
`--jev-prior shown`; `hidden` is a different experimental condition. Simulations
are not presented as proven seven-day predictions.
