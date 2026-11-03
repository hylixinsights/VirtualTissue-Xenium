#!/bin/zsh
VTX_SCRIPT_DIR="${0:A:h}"
source "$VTX_SCRIPT_DIR/_common.zsh"
if [[ $# -eq 0 || "$1" == "--help" || "$1" == "-h" ]]; then
  cat <<'HELP'
Usage: 04_jev.command estimate PREPARED.json NEW_DIRECTORY [options]
     04_jev.command run PREPARED.json NEW_DIRECTORY --max-http-attempts N [options]
Options: --steps 2 --seed 21 --model jev-1.14.0 --jev-prior shown|hidden
        --max-questions 20 --max-population 2000 --workers 1

estimate only writes cost.json; no key or network is needed. run records the
estimate and makes paid calls. It requires an explicit HTTP attempt limit and
TYPESAFE_API_KEY or JEV_API_KEY in the environment, never in arguments or project files.
The limit counts HTTP attempts, including failures; it is NOT a monetary cap.
The estimate uses t0; population changes, failures, and retries can change usage.
The estimate default is 20 attempts. Execution and resume are never automatic.
HELP
  exit 0
fi
[[ $# -ge 3 ]] || vtx_fail "Provide estimate/run, a bundle, and a new directory (see --help)."
action="$1"; bundle="$2"; out="$3"; shift 3
[[ "$action" == estimate || "$action" == run ]] || vtx_fail "Action must be estimate or run."
steps=2; seed=21; model=jev-1.14.0; prior=shown; attempts=20; explicit_budget=0
max_questions=20; max_population=2000; workers=1
while (( $# )); do
  [[ $# -ge 2 ]] || vtx_fail "Missing value for $1."
  case "$1" in
    --steps) steps="$2" ;;
    --seed) seed="$2" ;;
    --model) model="$2" ;;
    --jev-prior) prior="$2" ;;
    --max-http-attempts) attempts="$2"; explicit_budget=1 ;;
    --max-questions) max_questions="$2" ;;
    --max-population) max_population="$2" ;;
    --workers) workers="$2" ;;
    *) vtx_fail "Unknown option: $1" ;;
  esac
  shift 2
done
vtx_require_file "$bundle"
vtx_positive_integer steps "$steps"; vtx_seed "$seed"
vtx_positive_integer max-http-attempts "$attempts"; vtx_positive_integer max-questions "$max_questions"
vtx_positive_integer max-population "$max_population"; vtx_positive_integer workers "$workers"
[[ "$prior" == shown || "$prior" == hidden ]] || vtx_fail "--jev-prior must be shown or hidden."
if [[ "$action" == run ]]; then
  (( explicit_budget )) || vtx_fail "Paid execution requires an explicit --max-http-attempts."
  [[ -n "${TYPESAFE_API_KEY:-}${JEV_API_KEY:-}" ]] || vtx_fail "Set TYPESAFE_API_KEY or JEV_API_KEY in the environment before run."
fi
vtx_runtime
vtx_new_directory "$out"
typeset -a common_args
common_args=(--steps "$steps" --seed "$seed" --model "$model" --jev-prior "$prior" --max-questions "$max_questions" --max-http-attempts "$attempts")
vtx_cli cost "$bundle" "${common_args[@]}" > "$out/cost.json"
print -r -- "Estimate generated without API calls: ${out:a}/cost.json"
if [[ "$action" == run ]]; then
  vtx_cli run "$bundle" "$out/run_jev" --mode jev "${common_args[@]}" \
    --workers "$workers" --max-population "$max_population" > "$out/run_report.json"
  print -r -- "Jev execution complete: ${out:a}/run_jev/replay.html"
fi
