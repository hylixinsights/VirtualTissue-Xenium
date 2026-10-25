"""Execute every code cell of the three notebooks with a mocked Colab and a mocked Jev API.

Adapted from the original Colab_Xenium_Jev/tests/run_notebook_local.py. Only the Colab UI
(Drive mount, Secrets, file upload/download, widget manager), "/content" paths, form values and
the network transport are replaced. Code comes from this repository through the notebooks'
own CODE_SOURCE="google_drive" path, and the pinned packages must already be installed
(PIP_NO_INDEX=1 proves setup works offline).

    python tests/run_notebooks_mocked.py --output /tmp/nb_check \
        --data-dir /path/to/Xenium_Public_Lung5K            # Part 1 (needs the 10 raw files)
        --bundle /path/to/simulation_input_R03_*.zip         # Parts 2-4 (skip --data-dir to test only these)
"""
import argparse
import ast
import contextlib
import io
import json
import os
import shutil
import sys
import types
import zipfile
from pathlib import Path
from unittest.mock import patch

import nbformat

REPO = Path(__file__).resolve().parents[1]
FAKE_KEY = "FAKE_NOTEBOOK_KEY_DO_NOT_SAVE"


def mock_score_response(body):
    from validate_r03 import mock_score_response as respond
    return respond(body)


class Forms(ast.NodeTransformer):
    def __init__(self, values):
        self.values = values

    def visit_Assign(self, node):
        if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id in self.values:
            node.value = ast.Constant(self.values[node.targets[0].id])
        return node


def execute(notebook, values, content, log, calls, downloads, shown, drive_mounts, upload=None):
    nb = nbformat.read(REPO / "notebooks" / notebook, as_version=4)
    nbformat.validate(nb)
    google = types.ModuleType("google")
    colab = types.ModuleType("google.colab")
    colab.drive = types.SimpleNamespace(mount=lambda path, **kw: drive_mounts.append(path))
    colab.files = types.SimpleNamespace(upload=lambda: {upload.name: upload.read_bytes()} if upload else {},
                                        download=lambda p: downloads.append(p))
    colab.userdata = types.SimpleNamespace(get=lambda name: FAKE_KEY if name == "TYPESAFE_API_KEY" else None)
    colab.output = types.SimpleNamespace(enable_custom_widget_manager=lambda: None)
    google.colab = colab

    class Response:
        def __init__(self, body): self.raw = json.dumps(body).encode()
        def __enter__(self): return self
        def __exit__(self, *a): pass
        def read(self, n): return self.raw[:n]

    class Opener:
        def open(self, req, timeout):
            assert req.full_url == "https://api.typesafe.ai/v1/systemone"
            assert req.get_header("Authorization") == "Bearer " + FAKE_KEY
            body = json.loads(req.data)
            calls.append(body)
            return Response(mock_score_response(body))

    ns = {"display": lambda *a, **kw: shown.extend(a)}
    # Only the fake google.colab modules are added and removed: patch.dict(sys.modules) would also
    # drop numpy/h5py imported by the notebook, and compiled modules cannot be imported twice in
    # one process. (In Colab each notebook is its own process.)
    fake = {"google": google, "google.colab": colab, "google.colab.drive": colab.drive}
    saved = {k: sys.modules.get(k) for k in fake}
    sys.modules.update(fake)
    try:
        with patch("urllib.request.build_opener", return_value=Opener()), contextlib.redirect_stdout(log):
            import IPython.display
            with patch.object(IPython.display, "display", side_effect=lambda *a, **kw: shown.extend(a)):
                for i, cell in enumerate(nb.cells):
                    if cell.cell_type != "code":
                        continue
                    text = cell.source.replace("/content", str(content))
                    tree = ast.fix_missing_locations(Forms(values).visit(ast.parse(text)))
                    exec(compile(tree, f"{notebook}_cell_{i}", "exec"), ns)
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    return ns


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--data-dir", help="Public example data folder (raw/ + manifests) for Part 1.")
    parser.add_argument("--bundle", help="A simulation_input_*.zip for Parts 2-4 if Part 1 is not run.")
    parser.add_argument("--region", default="R03")
    args = parser.parse_args()
    out = Path(args.output).resolve()
    shutil.rmtree(out, ignore_errors=True)
    content = out / "content"
    drive_root = content / "drive" / "MyDrive"
    drive_root.mkdir(parents=True)
    os.environ["PIP_NO_INDEX"] = "1"
    os.environ["MPLBACKEND"] = "Agg"
    sys.path.insert(0, str(REPO / "tests"))
    report = dict(python=sys.version.split()[0])
    common = dict(CODE_SOURCE="google_drive", DRIVE_CODE_DIR=str(REPO))
    log = io.StringIO()
    calls, downloads, shown, mounts = [], [], [], []

    # ---------------------------------------------------------------- Part 1
    bundle = Path(args.bundle) if args.bundle else None
    if args.data_dir:
        values = dict(common, DATA_DIR=str(Path(args.data_dir).resolve()),
                      OUTPUT_BASE=str(drive_root / "VirtualTissue" / "Xenium_Lung5K_Prepared"),
                      REGION_TO_EXPORT=args.region, DOWNLOAD_ZIP=True)
        ns = execute("Part1_Prepare_Xenium_Region.ipynb", values, content, log, calls, downloads, shown, mounts)
        assert ns["SETUP_READY"]
        bundle = Path(ns["BUNDLE_ZIP_FOR_PART2"])
        assert bundle.is_file() and bundle.name.startswith(f"simulation_input_{args.region}_")
        assert downloads and Path(downloads[-1]) == bundle
        report["part1"] = dict(regions=list(ns["region_table"].region_id), exported=str(bundle),
                               atlas_fingerprint=ns["bundle"] and json.loads(
                                   __import__("gzip").decompress((ns["bundle"] / "atlas.json.gz").read_bytes()))["fingerprint"],
                               report={k: v for k, v in ns["report"].items() if k not in ("fidelity", "scope")})
    assert bundle, "Give --data-dir (runs Part 1) or --bundle."

    # ---------------------------------------------------------------- Parts 2-4 (Jev, mocked), Drive backup
    bundle_input = bundle.parent if args.data_dir else bundle      # a folder exercises "newest bundle inside"
    if bundle.is_dir():
        bundle = sorted(bundle.rglob("simulation_input_*.zip"), key=lambda p: p.stat().st_mtime)[-1]
    sim = dict(common, INPUT_SOURCE="google_drive", BUNDLE_ZIP=str(bundle_input),
               SAVE_RESULTS_TO_DRIVE=True, DRIVE_RESULTS_DIR=str(drive_root / "VirtualTissue" / "Simulation_Results"),
               HORIZON_HOURS=1, SCENARIO="declared_recognition", MAX_HTTP_ATTEMPTS=1000, START_RUN=True,
               RUN_NAME="notebook_validation", DOWNLOAD_RESULTS=True)
    ns = execute("Part2-4_Run_Jev_Simulation.ipynb", sim, content, log, calls, downloads, shown, mounts)
    result = ns["RESULT"]
    assert result["report"]["completed"], result["report"]
    first_calls = len(calls)
    assert first_calls == result["report"]["usage"]["http_attempts"] > 0
    backup = drive_root / "VirtualTissue" / "Simulation_Results"
    assert (backup / "notebook_validation" / "jev_journal.jsonl").is_file()
    assert (backup / "notebook_validation_results.zip").is_file()
    assert any("iframe" in getattr(x, "data", "") for x in shown if isinstance(getattr(x, "data", None), str))
    report["parts2_4_jev_mock"] = dict(http_attempts=first_calls, decisions=result["episode"]["summary"]["decisions"],
                                       simulated_minutes=result["report"]["checkpoint_time_min"],
                                       stop_reason=result["report"]["stop_reason"])

    # New runtime: session disk lost, Drive copy restored, no new paid calls.
    shutil.rmtree(content / "Xenium_Jev")
    ns = execute("Part2-4_Run_Jev_Simulation.ipynb", sim, content, log, calls, downloads, shown, mounts)
    assert len(calls) == first_calls, "Restored run must not send new requests."
    usage = ns["RESULT"]["report"]["usage"]
    assert usage["reused_batches"] == first_calls and usage["http_attempts"] == first_calls
    assert ns["RESULT"]["episode"]["final_state_fingerprint"] == result["episode"]["final_state_fingerprint"]
    report["restore_from_drive"] = dict(new_http_calls=len(calls) - first_calls, reused_batches=usage["reused_batches"])

    # Upload path + free local preview.
    local = dict(sim, INPUT_SOURCE="upload", RUN_MODE="Local preview (no API)", RUN_NAME="local_check",
                 SAVE_RESULTS_TO_DRIVE=False)
    ns = execute("Part2-4_Run_Jev_Simulation.ipynb", local, content, log, calls, downloads, shown, mounts, upload=bundle)
    assert ns["RESULT"]["report"]["completed"] and ns["RESULT"]["report"]["usage"]["http_attempts"] == 0
    assert len(calls) == first_calls
    report["local_preview_upload"] = dict(decisions=ns["RESULT"]["episode"]["summary"]["decisions"])

    # ---------------------------------------------------------------- Part 4 explorer
    explore = dict(common, RESULTS_SOURCE="google_drive",
                   RESULTS_PATH=str(backup / "notebook_validation_results.zip"), CELL_ID=0,
                   COMPARE_PATHS=f"{backup / 'notebook_validation_results.zip'},{content / 'Xenium_Jev/results/local_check_local_preview'}")
    ns = execute("Part4_Explore_Saved_Results.ipynb", explore, content, log, calls, downloads, shown, mounts)
    assert not ns["DECISIONS"].empty and ns["cell"] > 0
    report["part4"] = dict(summary=ns["vtx_results"].summary(ns["RES"]))

    # ---------------------------------------------------------------- secrets never persisted
    leaked = []
    for p in out.rglob("*"):
        if p.is_file() and p.suffix in (".json", ".jsonl", ".html", ".csv", ".gz", ".zip", ".txt"):
            if p.suffix == ".zip":
                with zipfile.ZipFile(p) as z:
                    if any(FAKE_KEY.encode() in z.read(n) for n in z.namelist()):
                        leaked.append(str(p))
            elif FAKE_KEY.encode() in p.read_bytes():
                leaked.append(str(p))
    assert not leaked, leaked
    import resource
    report["peak_memory_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    report.update(passed=True, actual_network_calls=0, mocked_http_attempts=len(calls),
                  drive_mounts=len(mounts), secret_absent_from_outputs=True)
    (out / "execution.log").write_text(log.getvalue())
    (out / "notebook_validation.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
