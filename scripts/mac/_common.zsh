# Shared by the macOS launchers. Source only after setting VTX_SCRIPT_DIR.
set -eu
set -o pipefail
VTX_REPO="${VTX_SCRIPT_DIR:h:h}"

vtx_fail() { print -u2 -r -- "VirtualTissue: $*"; exit 1; }
vtx_require_file() { [[ -f "$1" && -r "$1" ]] || vtx_fail "File not found or not readable: $1"; }
vtx_positive_integer() { [[ "$2" == <-> && "$2" -gt 0 ]] || vtx_fail "$1 must be a positive integer."; }
vtx_seed() { [[ "$1" == <-> ]] || vtx_fail "The seed must be a nonnegative integer."; }
vtx_new_directory() {
  [[ ! -e "$1" && ! -L "$1" ]] || vtx_fail "Output already exists; choose a new directory: $1"
  mkdir -p -- "${1:h}"
  mkdir -- "$1"
}

vtx_find_python() {
  local candidate
  if [[ -n "${VTX_PYTHON:-}" ]]; then
    candidate="$VTX_PYTHON"
    if [[ "$candidate" != */* ]]; then
      candidate="$(command -v -- "$candidate")" || vtx_fail "VTX_PYTHON was not found."
    fi
    [[ -x "$candidate" ]] || vtx_fail "VTX_PYTHON must name an executable without arguments: $candidate"
  elif [[ -x "$VTX_REPO/.venv-mac/bin/python" ]]; then
    candidate="$VTX_REPO/.venv-mac/bin/python"
  else
    candidate=""
    local name found
    for name in python3.13 python3.12 python3; do
      found="$(command -v -- "$name" 2>/dev/null)" || continue
      if "$found" -c 'import sys; raise SystemExit(sys.version_info[:2] not in ((3,12),(3,13)))' 2>/dev/null; then
        candidate="$found"
        break
      fi
    done
    [[ -n "$candidate" ]] || vtx_fail "Python 3.12 or 3.13 was not found. Install it and set VTX_PYTHON=/path/to/python3.12."
  fi
  "$candidate" -c 'import sys; raise SystemExit(sys.version_info[:2] not in ((3,12),(3,13)))' 2>/dev/null ||
    vtx_fail "Python 3.12 or 3.13 is required: $candidate"
  VTX_RUNTIME_PYTHON="${candidate:a}"
}

vtx_runtime() {
  [[ -f "$VTX_REPO/src/virtualtissue_xenium/v2/cli.py" ]] || vtx_fail "Keep scripts/mac inside the complete repository."
  vtx_find_python
  # Always use this checkout, regardless of an older pip installation/PYTHONPATH.
  export PYTHONPATH="$VTX_REPO/src"
  export PYTHONNOUSERSITE=1
  export PYTHONDONTWRITEBYTECODE=1
}

vtx_cli() { "$VTX_RUNTIME_PYTHON" -m virtualtissue_xenium.v2.cli "$@"; }
