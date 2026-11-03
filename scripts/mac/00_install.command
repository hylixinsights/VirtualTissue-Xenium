#!/bin/zsh
VTX_SCRIPT_DIR="${0:A:h}"
source "$VTX_SCRIPT_DIR/_common.zsh"
typeset -a extras pip_options
extras=(); pip_options=()
dry_run=0
while (( $# )); do
  case "$1" in
    --help|-h)
      cat <<'HELP'
Usage: 00_install.command [--with-h5ad] [--with-tests] [--with-plots] [--offline] [--dry-run]

Creates .venv-mac in this repository using Python 3.12/3.13. The basic setup
uses only the standard library, downloads no packages, and loads src/ directly.
--with-h5ad installs the xenium extra; --with-tests installs pytest and the
current test dependencies (xenium, spatial, and analysis).
--with-plots installs matplotlib. These extras use pip and may access the network.
--offline disables package indexes; extras must already be available locally.
--dry-run prints the configuration without creating or installing anything.
Set VTX_PYTHON=/path/to/python3.12 to select the interpreter.
This does not install Homebrew or change the system Python.
HELP
      exit 0 ;;
    --with-h5ad) extras+=(xenium) ;;
    --with-tests) extras+=(test xenium spatial analysis) ;;
    --with-plots) extras+=(analysis) ;;
    --offline) pip_options+=(--no-index) ;;
    --dry-run) dry_run=1 ;;
    *) vtx_fail "Unknown option: $1 (use --help)." ;;
  esac
  shift
done
vtx_runtime
venv_dir="$VTX_REPO/.venv-mac"
print -r -- "Python: $VTX_RUNTIME_PYTHON"
print -r -- "Environment: $venv_dir"
print -r -- "Extras: ${extras[*]:-none}"
if (( dry_run )); then exit 0; fi
if [[ -e "$venv_dir" || -L "$venv_dir" ]]; then
  [[ -f "$venv_dir/pyvenv.cfg" && -x "$venv_dir/bin/python" ]] || vtx_fail ".venv-mac already exists and is not a valid environment; nothing was removed."
else
  "$VTX_RUNTIME_PYTHON" -m venv "$venv_dir"
fi
"$venv_dir/bin/python" -c 'import sys; raise SystemExit(sys.version_info[:2] not in ((3,12),(3,13)))' || vtx_fail "The existing environment uses an unsupported Python version."
if (( ${#extras} )); then
  requirement_lines="$("$VTX_RUNTIME_PYTHON" -c 'import sys,tomllib; p=tomllib.load(open(sys.argv[1],"rb"))["project"]["optional-dependencies"]; print("\n".join(dict.fromkeys(r for extra in sys.argv[2:] for r in p[extra])))' "$VTX_REPO/pyproject.toml" "${extras[@]}")"
  typeset -a requirements
  requirements=("${(@f)requirement_lines}")
  PIP_DISABLE_PIP_VERSION_CHECK=1 "$venv_dir/bin/python" -m pip install "${pip_options[@]}" "${requirements[@]}"
fi
VTX_PYTHON="$venv_dir/bin/python" "$VTX_SCRIPT_DIR/03_check.command"
print -r -- "Ready. Launchers use .venv-mac when VTX_PYTHON is not set."
