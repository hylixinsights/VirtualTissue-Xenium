#!/bin/zsh
VTX_SCRIPT_DIR="${0:A:h}"
source "$VTX_SCRIPT_DIR/_common.zsh"
bundle=""; run_tests=0
while (( $# )); do
  case "$1" in
    --help|-h)
      cat <<'HELP'
Usage: 03_check.command [PREPARED.json] [--tests]
Checks Python, imported code location, engine hashes, ontology, and signatures.
With a bundle, also runs inspect. --tests runs the current v2/v3/Mac and export
suite (install --with-tests). The legacy engine suite is separate.
Makes no Jev calls. These are engineering checks, not biological validation.
HELP
      exit 0 ;;
    --tests) run_tests=1 ;;
    --*) vtx_fail "Unknown option: $1" ;;
    *) [[ -z "$bundle" ]] || vtx_fail "Provide at most one bundle."; bundle="$1" ;;
  esac
  shift
done
[[ -z "$bundle" ]] || vtx_require_file "$bundle"
vtx_runtime
"$VTX_RUNTIME_PYTHON" "$VTX_SCRIPT_DIR/check_environment.py" "$VTX_REPO"
if [[ -n "$bundle" ]]; then vtx_cli inspect "$bundle"; fi
if (( run_tests )); then
  typeset -a current_tests
  current_tests=("$VTX_REPO"/tests/test_v2*.py(N) "$VTX_REPO"/tests/test_v3*.py(N)
                 "$VTX_REPO"/tests/test_mac_workflow.py(N))
  (( ${#current_tests} )) || vtx_fail "Current tests were not found; use the complete repository."
  (cd "$VTX_REPO" && "$VTX_RUNTIME_PYTHON" -m pytest "${current_tests[@]}" -q)
fi
