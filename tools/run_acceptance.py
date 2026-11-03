#!/usr/bin/env python3
"""One local acceptance workflow; exact-source GitHub handoff only after success."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import fcntl
import hashlib
import html
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time
import traceback
import uuid
import xml.etree.ElementTree as ET
import zipfile

from release_integrity import snapshot, compare_snapshots, write_github_payload

REPO = Path(__file__).resolve().parents[1]
ALLOWED_SKIPS = {
    ("tests.test_helpers.BundleRunTests", "test_cli_check_run_local_and_explore"),
    ("tests.test_release_012.RealBundleRegression", "test_reporting_preserves_one_hour_biological_state_and_legacy_bundle"),
    ("tests.test_release_013", "test_real_bundle_preflight_does_not_confuse_scheduler_with_decision_cap"),
}


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def digest_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_junit(path, *, no_skips=False, minimum=1):
    report = ET.parse(path)
    cases = list(report.iter("testcase"))
    # Collection errors may be attached to a suite, without a testcase.
    reported_errors = sum(int(node.get(key, "0")) for node in report.iter()
                          if node.tag in ("testsuites", "testsuite") for key in ("errors", "failures"))
    error_elements = list(report.iter("error")) + list(report.iter("failure"))
    failures = [c.attrib for c in cases if c.find("failure") is not None or c.find("error") is not None]
    skipped = [c for c in cases if c.find("skipped") is not None]
    unexpected = [c.attrib for c in skipped if no_skips or (c.get("classname"), c.get("name")) not in ALLOWED_SKIPS
                  or "VTX_TEST_BUNDLE" not in c.find("skipped").get("message", "")]
    passed = len(cases) - len(skipped) - len(failures)
    if failures or reported_errors or error_elements or unexpected or passed < minimum:
        raise RuntimeError(f"Test gate failed: {len(failures)} failed cases, {reported_errors} suite/root error declarations, {len(error_elements)} error elements, {len(unexpected)} unexpected skips, {passed} passes (minimum {minimum}).")
    return dict(passed=passed, skipped=len(skipped), failures=len(failures),
                skipped_tests=[dict(c.attrib, reason=c.find("skipped").get("message", "")) for c in skipped])


def page(root, result):
    state = result["status"]
    title = {"running": "TEST RUNNING", "passed": "PASS", "failed": "FAIL", "interrupted": "INCOMPLETE"}[state]
    color = {"running": "#245caf", "passed": "#146b46", "failed": "#b33131", "interrupted": "#93610b"}[state]
    links = ""
    for label, target in result.get("links", {}).items():
        links += f'<a class="button" href="{html.escape(target, quote=True)}">{html.escape(label)}</a>'
    stages = "".join(f'<li>{html.escape(item["name"])} — {html.escape(item["status"])}</li>' for item in result.get("stages", []))
    counts = result.get("tests")
    checks = f'<p><strong>{counts["passed"]} tests passed</strong>; {counts["skipped"]} explicitly listed legacy tests skipped.</p>' if counts else ""
    message = {"passed": "This tested source snapshot is ready for GitHub review. Nothing has been uploaded.",
               "running": "Please wait. A previous successful run does not apply to this run.",
               "failed": "Do not update GitHub from this run. The error and logs below explain what needs attention.",
               "interrupted": "The test did not finish. Run START_TEST.command again before updating GitHub."}[state]
    error = '<pre>' + html.escape(result.get("error", "")) + '</pre>' if result.get("error") else ""
    source = html.escape(result.get("source_tree_sha256", "not yet recorded"))
    document = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>VirtualTissue — {title}</title><style>body{{font:17px/1.6 system-ui,sans-serif;background:#f3f5f8;color:#172638;margin:0}}main{{max-width:850px;margin:40px auto;padding:32px;background:white;border-radius:16px}}h1{{color:{color};font-size:46px;margin:8px 0}}.button{{display:inline-block;margin:6px 8px 6px 0;padding:11px 17px;border-radius:8px;background:#235caf;color:white;text-decoration:none}}small,.muted{{color:#536174}}code,pre{{overflow-wrap:anywhere;white-space:pre-wrap;background:#f4f6f8;padding:10px;border-radius:6px}}ul{{padding-left:24px}}@media(max-width:600px){{main{{margin:12px;padding:20px}}}}</style>
<main><small>VirtualTissue {html.escape(result.get("version", "v3"))} · Local acceptance test</small><h1>{title}</h1><p>{message}</p>{checks}{links}{error}
<h2>What was checked</h2><ul>{stages}</ul><p class="muted">The catalogue uses synthetic data. It starts with no region selected. The automated workflow uses a separate fixed test selection; it is not a choice made by you.</p>
<p class="muted">A software PASS does not establish biological validation, real Jev behaviour or a proven prediction of future tissue. No paid API calls were made.</p>
<details><summary>Exact source and logs</summary><p>Run: {html.escape(result["run_id"])}<br>Source SHA-256: <code>{source}</code></p><a href="{html.escape(result["run_relative"], quote=True)}/validation.json">Validation record</a> · <a href="{html.escape(result["run_relative"], quote=True)}/logs/">Logs</a><p>The source package and patch in this run's GitHub handoff belong only to this snapshot. Re-run the test after changing the code.</p></details></main></html>'''
    temp = root / "RESULT.html.tmp"
    temp.write_text(document)
    temp.replace(root / "RESULT.html")


def clean_environment(source, output):
    """Do not inherit test selection or provider credentials from a user's shell."""
    env = dict(source, PYTHONPATH=str(REPO / "src"), PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1",
               MPLCONFIGDIR=str(output / "cache/matplotlib"), NUMBA_CACHE_DIR=str(output / "cache/numba"),
               PIP_DISABLE_PIP_VERSION_CHECK="1", VTX_PYTHON=sys.executable,
               PYTEST_DISABLE_PLUGIN_AUTOLOAD="1")
    for key in ("TYPESAFE_API_KEY", "JEV_API_KEY", "VTX_TEST_BUNDLE", "PYTEST_ADDOPTS", "PYTEST_PLUGINS", "PYTEST_CURRENT_TEST"):
        env.pop(key, None)
    return env


def run_workflow(root):
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    lock_handle = (root / ".test.lock").open("a+")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print("A test is already running in this folder.", file=sys.stderr)
        return 2
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:6]
    out = root / "results" / run_id
    (out / "logs").mkdir(parents=True)
    result = dict(status="running", run_id=run_id, run_relative=out.relative_to(root).as_posix(), stages=[],
                  started_at_utc=datetime.now(timezone.utc).isoformat(), paid_api_calls=0,
                  github_published=False, biological_validation="not established", commands=[])
    env = clean_environment(os.environ, out)
    def save():
        write_json(out / "validation.json", result)
        write_json(root / "RESULT.json", result)
        page(root, result)
    def stage(name):
        if result["stages"] and result["stages"][-1]["status"] == "running":
            result["stages"][-1]["status"] = "passed"
        result["stages"].append(dict(name=name, status="running"))
        print(f'[{len(result["stages"])}/8] {name}', flush=True)
        save()
    def command(name, args, *, cwd=REPO, custom_env=None):
        target = out / "logs" / (name + ".log")
        entry = dict(name=name, argv=list(map(str, args)), cwd=str(cwd), log=target.relative_to(out).as_posix())
        result["commands"].append(entry)
        started = time.monotonic()
        with target.open("w") as log:
            completed = subprocess.run(list(map(str, args)), cwd=cwd, env=custom_env or env, stdout=log, stderr=subprocess.STDOUT)
        entry.update(returncode=completed.returncode, seconds=round(time.monotonic()-started, 3))
        if completed.returncode:
            raise RuntimeError(f'{name} failed (exit {completed.returncode}). See {target}.\n' + target.read_text(errors="replace")[-3500:])
        return target
    def cli(name, *args):
        return command(name, [sys.executable, "-m", "virtualtissue_xenium.v2.cli", *args])
    save()
    try:
        stage("Environment and source integrity")
        if sys.version_info[:2] not in ((3, 12), (3, 13)):
            raise RuntimeError("Python 3.12 or 3.13 is required.")
        for executable in ("node", "zsh", "git"):
            if not shutil.which(executable):
                raise RuntimeError(f"{executable} is required; its checks must not be skipped.")
        before = snapshot(REPO)
        base_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
        result.update(source_tree_sha256=before["tree_sha256"], git_base=base_sha,
                      version=json.loads(command("environment", [sys.executable, REPO / "scripts/mac/check_environment.py", REPO]).read_text())["package_version"])
        result["environment"] = dict(python=sys.version, executable=sys.executable,
            node=subprocess.check_output(["node", "--version"], text=True).strip(),
            distributions=sorted({(d.metadata["Name"], d.version) for d in importlib.metadata.distributions()}))
        write_json(out / "source_snapshot.json", before)

        stage("Complete test suite and mandatory browser-player checks")
        command("pytest", [sys.executable, "-m", "pytest", "tests", "-q", "-ra", "--junitxml=" + str(out / "tests.xml")])
        result["tests"] = check_junit(out / "tests.xml", minimum=365)
        command("legacy_player", ["node", REPO / "tests/test_player_logic.cjs"])

        stage("Local time-step consistency benchmark")
        command("clock", [sys.executable, REPO / "tools/validate_clock_v3.py", "--output", out / "clock"])
        result["clock"] = json.loads((out / "clock/validation.json").read_text())

        stage("Object export, ROI catalogue and fixed test selection")
        demo = out / "demo"
        command("generate", [sys.executable, REPO / "examples/v3/make_object_demo.py", "--output", demo / "original.json", "--context-output", demo / "context.json", "--include-invalid-polygon"])
        command("export", [sys.executable, REPO / "tools/export_object_v3.py", demo / "original.json", demo / "object.json"])
        cli("catalog", "roi-catalog", demo / "object.json", demo / "catalog.html", "--cells-per-roi", "32", "--grid-um", "40", "--max-candidates", "12")
        catalog = json.loads((demo / "catalog.catalog.json").read_text())
        if catalog["selected_candidate"] is not None:
            raise RuntimeError("The interactive catalogue must start without a selection.")
        write_json(demo / "automated_fixture.selection.json", catalog["candidates"][0]["selection"])
        cli("select", "select-roi", demo / "object.json", demo / "automated_fixture.selection.json", demo / "prepared.json", "--context", demo / "context.json")
        cli("expression", "expression-audit", demo / "prepared.json", demo / "expression_audit.json")
        audit = json.loads((demo / "expression_audit.json").read_text())
        if audit["comparisons"]["beta_zero"]["summary"]["n_probability_changed"] < 1:
            raise RuntimeError("The assigned synthetic expression pattern did not affect local probabilities.")
        result["demo"] = dict(catalogue_candidates=len(catalog["candidates"]), interactive_selection=None,
            automated_selection="First deterministic synthetic fixture candidate; not a user's region choice.",
            beta_zero_probability_changes=audit["comparisons"]["beta_zero"]["summary"]["n_probability_changed"])

        stage("Matched local simulations and replay artifacts")
        cli("beta_zero", "control", demo / "prepared.json", demo / "beta_zero.json", "--kind", "beta_zero")
        cli("shuffle", "control", demo / "prepared.json", demo / "shuffle.json", "--kind", "shuffle_expression_within_type", "--seed", "11")
        write_json(demo / "arms.json", dict(arms=dict(baseline="prepared.json", beta_zero="beta_zero.json", within_type_shuffle="shuffle.json")))
        cli("ensemble", "ensemble", demo / "arms.json", demo / "ensemble", "--steps", "8", "--seeds", "0", "1", "2", "--bootstrap-replicates", "200")
        report = json.loads((demo / "ensemble/report.json").read_text())
        if not report["completed"] or report["n_completed_runs"] != 9 or report["cost_ledger"]["total"]["http_attempts"] != 0:
            raise RuntimeError("The local ensemble did not complete all nine runs without API calls.")
        frames_checked = 0
        for run in sorted((demo / "ensemble/runs").glob("*/*")):
            frames = [json.loads(line) for line in (run / "frames.jsonl").read_text().splitlines()]
            if frames != json.loads((run / "frames.json").read_text()):
                raise RuntimeError("Frame encodings disagree.")
            rehydrated = []
            for chunk in sorted((run / "replay_frames").glob("chunk_*.js")):
                _, block, shapes = json.loads("[" + chunk.read_text().split("(", 1)[1].rsplit(");", 1)[0] + "]")
                for frame in block:
                    for cell in frame["cells"].values():
                        if "t0_geometry_ref" in cell:
                            cell["t0_geometry"] = shapes[cell.pop("t0_geometry_ref")]
                rehydrated.extend(block)
            if rehydrated != frames:
                raise RuntimeError("Replay geometry does not match the audited frames.")
            previous = sum(c["status"] == "alive" for c in frames[0]["cells"].values())
            for frame in frames:
                events = Counter(e["kind"] for e in frame["events"])
                living = sum(c["status"] == "alive" for c in frame["cells"].values())
                if living != previous + events["division"] - events["kill"] - events["apoptosis"]:
                    raise RuntimeError("Population conservation failed.")
                previous = living
            frames_checked += len(frames)
        result["demo"].update(completed_runs=9, frames_checked=frames_checked)
        cli("cost_only", "cost", demo / "prepared.json", "--steps", "8")

        stage("Build and isolated installed-package test")
        command("build", [sys.executable, "-m", "hatchling", "build", "-t", "wheel", "-t", "sdist", "-d", out / "packages"])
        wheel = next((out / "packages").glob("*.whl"))
        source = next((out / "packages").glob("*.tar.gz"))
        target = out / "installed_package"
        command("install", [sys.executable, "-m", "pip", "install", "--no-index", "--no-deps", "--target", target, wheel], cwd=out)
        installed_env = dict(env, PYTHONPATH=str(target))
        proof = command("installed_origin", [sys.executable, "-c", "import pathlib,sys,virtualtissue_xenium as v; p=pathlib.Path(v.__file__).resolve(); assert p.is_relative_to(pathlib.Path(sys.argv[1]).resolve()); print(v.__version__,p)", target], cwd=out, custom_env=installed_env)
        command("installed_run", [sys.executable, "-m", "virtualtissue_xenium.v2.cli", "run", demo / "prepared.json", out / "installed_smoke", "--steps", "2", "--seed", "21"], cwd=out, custom_env=installed_env)
        with zipfile.ZipFile(wheel) as archive:
            for name in archive.namelist():
                if name.startswith("virtualtissue_xenium/") and not name.endswith("/"):
                    if archive.read(name) != (REPO / "src" / name).read_bytes():
                        raise RuntimeError("The wheel does not match the tested source: " + name)
        with tarfile.open(source) as archive:
            names = {"/".join(Path(m.name).parts[1:]) for m in archive.getmembers() if m.isfile()}
            for required in ("START_TEST.command", "tools/run_acceptance.py", "tools/release_integrity.py", "scripts/mac/02_run_roi.command", "docs/v3/SCIENTIFIC_SCOPE.md"):
                if required not in names:
                    raise RuntimeError("Missing source-package file: " + required)
        result["packages"] = {p.name: dict(sha256=digest_file(p), bytes=p.stat().st_size) for p in (wheel, source)}

        stage("Verify unchanged source and prepare the GitHub handoff")
        compare_snapshots(before, snapshot(REPO))
        payload = write_github_payload(REPO, out / "github", before, base_sha)
        compare_snapshots(before, snapshot(REPO))
        result["github_payload"] = payload
        body = f'''Update VirtualTissue to {result["version"]}: user-selected ROIs, supplied polygons, expression-sensitive dynamics and audited scientific assumptions, with English interfaces and one acceptance-test entry point.\n\nValidation: {result["tests"]["passed"]} tests passed; {result["tests"]["skipped"]} optional legacy real-bundle tests skipped. The local clock benchmark, nine synthetic runs, replay consistency, wheel build and isolated installation passed. No live Jev calls or biological validation were performed.\n\nBase: {base_sha}\nTested source SHA-256: {before["tree_sha256"]}\n\nApply the accompanying patch only to this base and re-run CI before merging. This folder does not publish or push automatically.\n'''
        (out / "github/PR_DESCRIPTION.md").write_text(body)
        result["handoff_files"] = {p.relative_to(out).as_posix(): digest_file(p) for p in (out / "github").iterdir() if p.is_file()}

        stage("Publish this local result page")
        prefix = out.relative_to(root).as_posix()
        result["links"] = {"Open ROI catalogue": prefix + "/demo/catalog.html",
            "Open simulation replay": prefix + "/demo/ensemble/runs/baseline/seed_0/replay.html",
            "View GitHub handoff": prefix + "/github/PR_DESCRIPTION.md", "View test log": prefix + "/logs/pytest.log"}
        result["stages"][-1]["status"] = "passed"
        result["status"] = "passed"
        print("PASS — the tested snapshot is ready for GitHub review. Nothing was uploaded.", flush=True)
        returncode = 0
    except BaseException as exc:
        result["status"] = "interrupted" if isinstance(exc, (KeyboardInterrupt, SystemExit)) else "failed"
        result["error"] = str(exc) or type(exc).__name__
        if result["stages"]:
            result["stages"][-1]["status"] = result["status"]
        (out / "logs/failure.log").write_text(traceback.format_exc())
        print(result["status"].upper() + " — " + result["error"], file=sys.stderr, flush=True)
        returncode = 1
    result["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    save()
    fcntl.flock(lock_handle, fcntl.LOCK_UN)
    lock_handle.close()
    return returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=REPO.parent / (REPO.name + "-test-results"))
    parser.add_argument("--no-open", action="store_true", help="Used by the desktop launcher for unattended verification.")
    parser.add_argument("--check-junit", type=Path)
    parser.add_argument("--allow-no-skips", action="store_true", help="Require zero skips when checking a standalone JUnit file.")
    args = parser.parse_args()
    if args.check_junit:
        print(json.dumps(check_junit(args.check_junit, no_skips=args.allow_no_skips)))
        return 0
    result = run_workflow(args.output_root)
    if not args.no_open and sys.platform == "darwin" and (args.output_root / "RESULT.html").exists():
        subprocess.run(["/usr/bin/open", str((args.output_root / "RESULT.html").resolve())], check=False)
    return result


if __name__ == "__main__":
    raise SystemExit(main())
