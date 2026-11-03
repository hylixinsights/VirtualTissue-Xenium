"""Real offline zsh workflows, with a relocated checkout and paths with spaces."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which("zsh") is None, reason="Mac launchers require zsh")


@pytest.fixture
def mac_checkout(tmp_path):
    root = tmp_path / "checkout with spaces"
    root.mkdir()
    shutil.copytree(ROOT / "scripts", root / "scripts")
    for name in ("src", "examples", "tools", "tests"):
        (root / name).symlink_to(ROOT / name, target_is_directory=True)
    shutil.copyfile(ROOT / "pyproject.toml", root / "pyproject.toml")
    env = dict(os.environ, VTX_PYTHON=sys.executable, PYTHONPATH="/deliberately/wrong/site",
               PYTHONDONTWRITEBYTECODE="1", PIP_NO_INDEX="1", PIP_DISABLE_PIP_VERSION_CHECK="1")
    env.pop("TYPESAFE_API_KEY", None)
    env.pop("JEV_API_KEY", None)
    return root, env


def invoke(checkout, name, *args, ok=True):
    root, env = checkout
    result = subprocess.run(["zsh", str(root / "scripts/mac" / name), *map(str, args)],
                            cwd=root.parent, env=env, capture_output=True, text=True, timeout=60)
    if ok:
        assert result.returncode == 0, result.stdout + result.stderr
    else:
        assert result.returncode != 0, result.stdout + result.stderr
    return result


def snapshot(path):
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob("*") if p.is_file()}


def test_mac_help_and_interpreter_errors(mac_checkout):
    for name in ("vtx.command", "00_install.command", "01_demo.command", "02_run_roi.command",
                 "03_check.command", "04_jev.command"):
        assert "Usage:" in invoke(mac_checkout, name, "--help").stdout
    assert "export_object_v3.py" in invoke(mac_checkout, "05_export_object.command", "--help").stdout
    proof = json.loads(invoke(mac_checkout, "03_check.command").stdout)
    assert proof["status"] == "ok" and proof["lineages"] == 22
    assert proof["engine_files_verified"] == 35
    assert proof["source"] == str((ROOT / "src/virtualtissue_xenium/__init__.py").resolve())
    root, env = mac_checkout
    env["VTX_PYTHON"] = "/missing/python executable"
    failure = invoke((root, env), "03_check.command", ok=False)
    assert "VTX_PYTHON" in failure.stderr


def test_mac_installer_creates_isolated_environment_without_network(mac_checkout):
    root, env = mac_checkout
    invoke(mac_checkout, "00_install.command", "--offline", "--dry-run")
    assert not (root / ".venv-mac").exists()
    invoke(mac_checkout, "00_install.command", "--offline")
    assert (root / ".venv-mac/pyvenv.cfg").is_file()
    env.pop("VTX_PYTHON")
    proof = json.loads(invoke((root, env), "03_check.command").stdout)
    assert proof["executable"] == str(root / ".venv-mac/bin/python")
    # Reusing a valid environment is permitted; it never removes user files.
    marker = root / ".venv-mac/keep.txt"
    marker.write_text("do not replace")
    invoke((root, env), "00_install.command", "--offline")
    assert marker.read_text() == "do not replace"


def test_mac_test_runner_uses_repo_cwd_and_excludes_legacy_suite(mac_checkout):
    root, env = mac_checkout
    interpreter = root.parent / "interpreter shim"
    interpreter.write_text('''#!/bin/zsh
if [[ "$1" == -m && "$2" == pytest ]]; then
  print -r -- "$PWD" > "$VTX_TEST_LOG"
  printf '%s\\n' "$@" >> "$VTX_TEST_LOG"
  exit 0
fi
exec "$VTX_REAL_PYTHON" "$@"
''')
    interpreter.chmod(0o755)
    log = root.parent / "pytest arguments.txt"
    env.update(VTX_PYTHON=str(interpreter), VTX_REAL_PYTHON=sys.executable, VTX_TEST_LOG=str(log))
    invoke((root, env), "03_check.command", "--tests")
    lines = log.read_text().splitlines()
    assert lines[0] == str(root)
    assert any(p.endswith("test_v3_geometry.py") for p in lines)
    assert any(p.endswith("test_mac_workflow.py") for p in lines)
    assert not any(p.endswith("test_engine_integrity.py") for p in lines)
    options = invoke((root, env), "00_install.command", "--with-tests", "--offline", "--dry-run")
    assert "test xenium spatial analysis" in options.stdout


def test_mac_demo_selection_local_audit_and_cost_are_offline(mac_checkout):
    root, env = mac_checkout
    demo = root.parent / "synthetic demo with spaces"
    invoke(mac_checkout, "01_demo.command", demo)
    catalog = json.loads((demo / "catalog.catalog.json").read_text())
    assert catalog["selected_candidate"] is None
    assert len(catalog["candidates"]) == 12
    before = snapshot(demo)
    invoke(mac_checkout, "01_demo.command", demo, ok=False)
    assert snapshot(demo) == before
    exported = root.parent / "reexported raw object.json"
    invoke(mac_checkout, "05_export_object.command", demo / "object.json", exported)
    raw = json.loads(exported.read_text())
    assert len(raw["cells"]) == 96 and "guidebook" not in raw
    assert raw["annotation_provenance"]["object_export"]["scoring_status"] == "not_computed"
    exported_before = exported.read_bytes()
    invoke(mac_checkout, "05_export_object.command", demo / "object.json", exported, ok=False)
    assert exported.read_bytes() == exported_before
    selection = demo / "test selection.json"
    # Explicit test fixture choice only; the launchers never choose a candidate.
    selection.write_text(json.dumps(catalog["candidates"][0]["selection"]))
    output = root.parent / "local run with spaces"
    invoke(mac_checkout, "02_run_roi.command", demo / "object.json", selection,
           demo / "context.json", output, "--steps", "2", "--seed", "41")
    assert (output / "expression_audit.json").is_file()
    config = json.loads((output / "run_local/config.json").read_text())
    assert config["mode"] == "local"
    assert not list(output.rglob("jev_journal"))
    assert len(json.loads((output / "run_local/frames.json").read_text())) == 3
    assert (output / "run_local/replay.html").is_file()
    assert list((output / "run_local").rglob("*.js"))
    before = snapshot(output)
    invoke(mac_checkout, "02_run_roi.command", demo / "object.json", selection,
           demo / "context.json", output, "--steps", "2", ok=False)
    assert snapshot(output) == before
    estimate = root.parent / "offline Jev estimate"
    invoke(mac_checkout, "04_jev.command", "estimate", output / "prepared.json", estimate)
    assert (estimate / "cost.json").is_file()
    assert not (estimate / "run_jev").exists()
    paid = root.parent / "must not call Jev"
    missing_budget = invoke(mac_checkout, "04_jev.command", "run", output / "prepared.json", paid, ok=False)
    assert "max-http-attempts" in missing_budget.stderr and not paid.exists()
    missing_key = invoke(mac_checkout, "04_jev.command", "run", output / "prepared.json", paid,
                         "--max-http-attempts", "1", ok=False)
    assert "API_KEY" in missing_key.stderr and not paid.exists()


def test_mac_rejects_bad_inputs_before_creating_outputs(mac_checkout):
    root, _ = mac_checkout
    output = root.parent / "should not exist"
    invoke(mac_checkout, "02_run_roi.command", "missing.json", "missing.json", "missing.json", output, ok=False)
    assert not output.exists()
    invoke(mac_checkout, "02_run_roi.command", "a", "b", "c", output, "--steps", "0", ok=False)
    assert not output.exists()
    invoke(mac_checkout, "04_jev.command", "run", "a", output, "--workers", ok=False)
    assert not output.exists()
