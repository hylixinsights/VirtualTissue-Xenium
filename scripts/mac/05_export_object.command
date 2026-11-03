#!/bin/zsh
VTX_SCRIPT_DIR="${0:A:h}"
source "$VTX_SCRIPT_DIR/_common.zsh"
vtx_runtime
vtx_require_file "$VTX_REPO/tools/export_object_v3.py"
if [[ $# -eq 0 ]]; then
  print -r -- "Export H5AD to a raw JSON object without computing scores or selecting an ROI."
  print -r -- "Requires the installer --with-h5ad extra. Exporter help:"
  set -- --help
fi
exec "$VTX_RUNTIME_PYTHON" "$VTX_REPO/tools/export_object_v3.py" "$@"
