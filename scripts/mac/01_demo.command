#!/bin/zsh
VTX_SCRIPT_DIR="${0:A:h}"
source "$VTX_SCRIPT_DIR/_common.zsh"
if [[ $# -gt 0 && ( "$1" == "--help" || "$1" == "-h" ) ]]; then
  cat <<'HELP'
Usage: 01_demo.command [NEW_DIRECTORY]
Generates 96 explicitly SYNTHETIC cells, polygons, context, and an HTML catalog.
Includes one declared invalid outline and one low-count cell to exercise QC.
Does not select an ROI, simulate, call Jev, or open a browser automatically.
Without an argument, uses outputs/mac/demo_DATE_TIME in the repository.
HELP
  exit 0
fi
[[ $# -le 1 ]] || vtx_fail "Usage: 01_demo.command [NEW_DIRECTORY]"
out="${1:-$VTX_REPO/outputs/mac/demo_$(date +%Y%m%d_%H%M%S)}"
[[ "$out" != --* ]] || vtx_fail "Unknown option: $out"
vtx_runtime
vtx_new_directory "$out"
"$VTX_RUNTIME_PYTHON" "$VTX_REPO/examples/v3/make_object_demo.py" \
  --output "$out/object.json" --context-output "$out/context.json" \
  --include-invalid-polygon > "$out/generation.json"
vtx_cli roi-catalog "$out/object.json" "$out/catalog.html" \
  --cells-per-roi 32 --grid-um 40 --max-candidates 12 > "$out/catalog_report.json"
print -r -- "SYNTHETIC demonstration ready: ${out:a}"
print -r -- "Open the catalog and choose an ROI: ${out:a}/catalog.html"
print -r -- "The browser saves a selection JSON. Pass that file to 02_run_roi.command."
