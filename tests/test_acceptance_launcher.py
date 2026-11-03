"""Acceptance must not inherit a stale PASS or silently deselect required tests."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from run_acceptance import check_junit, clean_environment


def test_missing_interpreter_revokes_an_existing_pass(tmp_path):
    if not shutil.which("zsh"):
        pytest.skip("zsh is required for launcher validation")
    result = tmp_path / "old result"
    result.mkdir()
    (result / "RESULT.json").write_text('{"status":"passed"}')
    (result / "RESULT.html").write_text("PASS")
    proc = subprocess.run(["zsh", str(ROOT / "START_TEST.command"), "--output-root", str(result), "--no-open"],
                          env=dict(os.environ, VTX_TEST_PYTHON="/missing/test-python"), capture_output=True, text=True)
    assert proc.returncode != 0
    assert json.loads((result / "RESULT.json").read_text())["status"] == "failed"
    assert "FAIL" in (result / "RESULT.html").read_text()


def test_test_selection_environment_cannot_hide_a_failing_test(tmp_path):
    tests = tmp_path / "test_selected.py"
    tests.write_text('def test_should_run():\n    assert False, "must not be deselected"\n')
    env = clean_environment(dict(os.environ, PYTEST_ADDOPTS="-k absent_test", PYTEST_PLUGINS="nonexistent_plugin", JEV_API_KEY="sentinel"), tmp_path)
    proc = subprocess.run([sys.executable, "-m", "pytest", str(tests), "-q"], cwd=tmp_path,
                          env=env, capture_output=True, text=True)
    assert proc.returncode == 1
    assert "1 failed" in proc.stdout and "must not be deselected" in proc.stdout
    assert "JEV_API_KEY" not in env


def test_required_player_skip_rejects_otherwise_green_junit(tmp_path):
    p = tmp_path / "report.xml"
    p.write_text('<testsuites><testsuite><testcase classname="tests.test_v3_replay" name="player"><skipped message="Node unavailable"/></testcase><testcase classname="ok" name="ok"/></testsuite></testsuites>')
    with pytest.raises(RuntimeError, match="unexpected skips"):
        check_junit(p)


@pytest.mark.parametrize("content", [
    '<testsuites errors="1"><testsuite><testcase name="ok"/></testsuite></testsuites>',
    '<testsuites><testsuite failures="1"><testcase name="ok"/></testsuite></testsuites>',
    '<testsuites><testsuite><error message="collection failed"/><testcase name="ok"/></testsuite></testsuites>',
])
def test_suite_errors_cannot_look_like_a_pass(tmp_path, content):
    report = tmp_path / "suite_error.xml"
    report.write_text(content)
    with pytest.raises(RuntimeError, match="Test gate failed"):
        check_junit(report, no_skips=True)
