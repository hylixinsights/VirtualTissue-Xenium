#!/bin/zsh
VTX_SCRIPT_DIR="${0:A:h}"
source "$VTX_SCRIPT_DIR/_common.zsh"
if [[ $# -eq 0 || "$1" == "--help" || "$1" == "-h" ]]; then
  cat <<'HELP'
Usage: 02_run_roi.command OBJECT.json SELECTION.json CONTEXT.json NEW_DIRECTORY
       [--lineage-map MAP.csv] [--steps 12] [--seed 21] [--max-population 2000]

Uses the selected ROI and saves prepared.json, inspect.json, expression_audit.json,
and run_local/replay.html. Accepts raw or prepared JSON; context is explicit.
Use the same --lineage-map used to create the catalog.
Always runs locally: no API key is read and no Jev calls are made. Existing output is refused.
HELP
  exit 0
fi
[[ $# -ge 4 ]] || vtx_fail "Object, selection, context, and a new output directory are required (see --help)."
object="$1"; selection="$2"; context="$3"; out="$4"; shift 4
steps=12; seed=21; max_population=2000
typeset -a map_args
map_args=()
while (( $# )); do
  [[ $# -ge 2 ]] || vtx_fail "Missing value for $1."
  case "$1" in
    --lineage-map) vtx_require_file "$2"; map_args=(--lineage-map "$2") ;;
    --steps) steps="$2" ;;
    --seed) seed="$2" ;;
    --max-population) max_population="$2" ;;
    *) vtx_fail "Unknown option: $1" ;;
  esac
  shift 2
done
vtx_positive_integer steps "$steps"; vtx_seed "$seed"; vtx_positive_integer max-population "$max_population"
vtx_require_file "$object"; vtx_require_file "$selection"; vtx_require_file "$context"
vtx_runtime
vtx_new_directory "$out"
vtx_cli select-roi "$object" "$selection" "$out/prepared.json" --context "$context" "${map_args[@]}" > "$out/prepare_report.json"
vtx_cli inspect "$out/prepared.json" > "$out/inspect.json"
vtx_cli expression-audit "$out/prepared.json" "$out/expression_audit.json" --seed "$seed" > /dev/null
vtx_cli run "$out/prepared.json" "$out/run_local" --mode local \
  --steps "$steps" --seed "$seed" --max-population "$max_population" > "$out/run_report.json"
print -r -- "Local simulation complete: ${out:a}/run_local/replay.html"
print -r -- "Keep replay.html and its chunk directory together when copying the output."
