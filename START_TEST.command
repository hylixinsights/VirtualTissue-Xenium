#!/bin/zsh
set -eu
PROJECT_DIR="${0:A:h}"
RESULT_ROOT="${PROJECT_DIR:h}/${PROJECT_DIR:t}-test-results"
TEST_NO_OPEN=0
for ((argument_index=1; argument_index <= $#; argument_index++)); do
  if [[ "${argv[$argument_index]}" == --no-open ]]; then
    TEST_NO_OPEN=1
  fi
  if [[ "${argv[$argument_index]}" == --output-root && $argument_index -lt $# ]]; then
    RESULT_ROOT="${argv[$((argument_index+1))]}"
  fi
done
mkdir -p "$RESULT_ROOT"
print -r -- '{"status":"starting","message":"A new run has started. Previous results do not apply."}' > "$RESULT_ROOT/RESULT.json"
print -r -- '<!doctype html><html lang="en"><meta charset="utf-8"><title>VirtualTissue — starting</title><body style="font:20px system-ui;padding:40px"><h1>TEST STARTING</h1><p>Preparing the environment. A previous successful run does not apply to this run.</p></body></html>' > "$RESULT_ROOT/RESULT.html"
TRAPEXIT() {
  local test_exit_code=$?
  if [[ "$test_exit_code" -ne 0 ]] && /usr/bin/grep -q '"status":"starting"' "$RESULT_ROOT/RESULT.json"; then
    print -r -- '{"status":"failed","message":"Test setup failed. See the Terminal output; do not update GitHub from this run."}' > "$RESULT_ROOT/RESULT.json"
    print -r -- '<!doctype html><html lang="en"><meta charset="utf-8"><title>VirtualTissue — FAIL</title><body style="font:20px system-ui;padding:40px"><h1>FAIL — SETUP</h1><p>See the Terminal output. No new GitHub handoff was approved.</p></body></html>' > "$RESULT_ROOT/RESULT.html"
    if [[ "$TEST_NO_OPEN" -eq 0 ]]; then
      /usr/bin/open "$RESULT_ROOT/RESULT.html" || true
    fi
  fi
}
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
unset TYPESAFE_API_KEY JEV_API_KEY PYTEST_ADDOPTS PYTEST_PLUGINS
TEST_PYTHON="${VTX_TEST_PYTHON:-}"
if [[ -z "$TEST_PYTHON" ]]; then
  for candidate in "$PROJECT_DIR/../.environment/bin/python" "$PROJECT_DIR/.venv-test/bin/python" python3.12 python3.13; do
    if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sys; raise SystemExit(sys.version_info[:2] not in ((3,12),(3,13)))' 2>/dev/null; then
      TEST_PYTHON="$(command -v "$candidate")"
      break
    fi
  done
fi
if [[ -z "$TEST_PYTHON" ]]; then
  print -u2 'Python 3.12 or 3.13 is required. No test was run.'
  exit 1
fi
export MPLCONFIGDIR="$PROJECT_DIR/outputs/.test-cache/matplotlib"
export NUMBA_CACHE_DIR="$PROJECT_DIR/outputs/.test-cache/numba"
if ! "$TEST_PYTHON" -c 'import importlib.util; names="pytest numpy pandas scipy anndata h5py pyarrow zarr matplotlib shapely decoupler marsilea adjustText hatchling".split(); raise SystemExit(any(importlib.util.find_spec(n) is None for n in names))'; then
  print 'Preparing the isolated test environment. First-time dependency installation may use the network.'
  "$TEST_PYTHON" -m venv "$PROJECT_DIR/.venv-test"
  TEST_PYTHON="$PROJECT_DIR/.venv-test/bin/python"
  "$TEST_PYTHON" -m pip install -e "${PROJECT_DIR}[xenium,test,spatial,analysis,acceptance]"
fi
export PYTHONPATH="$PROJECT_DIR/src"
exec "$TEST_PYTHON" "$PROJECT_DIR/tools/run_acceptance.py" "$@"
