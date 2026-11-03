#!/bin/zsh
VTX_SCRIPT_DIR="${0:A:h}"
source "$VTX_SCRIPT_DIR/_common.zsh"
if [[ $# -eq 0 || "$1" == "--help" || "$1" == "-h" ]]; then
  cat <<'HELP'
VirtualTissue — CLI for this checkout (Python 3.12/3.13).
Usage: scripts/mac/vtx.command COMMAND [arguments]

Useful commands: roi-catalog, select-roi, inspect, expression-audit, run, cost.
Example: scripts/mac/vtx.command inspect "/path/to/prepared.json"
Command help: scripts/mac/vtx.command run --help

The default simulation mode is local. Jev requires an explicit --mode jev and
an API key in the environment. Use launchers 00 through 05 for the guided workflow.
VTX_PYTHON selects an interpreter; otherwise .venv-mac takes priority.
HELP
  exit 0
fi
vtx_runtime
exec "$VTX_RUNTIME_PYTHON" -m virtualtissue_xenium.v2.cli "$@"
