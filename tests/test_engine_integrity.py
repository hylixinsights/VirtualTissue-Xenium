"""The packaged engine must match the recorded hashes, including documented revisions."""
import hashlib
import json
import re
import subprocess
import sys
import tomllib
import unittest
from pathlib import Path

import virtualtissue_xenium as vtx

ROOT = Path(__file__).resolve().parents[1]


class EngineIntegrity(unittest.TestCase):
    def test_every_engine_file_matches_the_recorded_version(self):
        record = json.loads((vtx.ENGINE_DIR / "ORIGINAL_SOURCES.json").read_text())["files"]
        present = {p.relative_to(vtx.ENGINE_DIR).as_posix() for p in vtx.ENGINE_DIR.rglob("*")
                   if p.is_file() and "__pycache__" not in p.parts and p.name != "ORIGINAL_SOURCES.json"}
        self.assertEqual(present, set(record), "Engine files were added or removed.")
        for name, expected in record.items():
            digest = hashlib.sha256((vtx.ENGINE_DIR / name).read_bytes()).hexdigest()
            self.assertEqual(digest, expected["sha256"], f"{name} differs from its recorded version.")

    def test_engine_versions(self):
        import simulation_workflow
        import xenium_prepare
        self.assertEqual(xenium_prepare.VERSION, vtx.PREPARATION_VERSION)
        self.assertEqual(simulation_workflow.VERSION, vtx.SIMULATION_VERSION)
        vtx.check_engine_origin()

    def test_bundle_runtime_and_simulator_share_one_folder(self):
        # export_region copies ROOT/prepared_runtime.py and ROOT/xenium_runtime into each bundle;
        # extract_bundle compares those files with SOURCE (the simulator's folder). One folder => always equal.
        import simulation_workflow
        import xenium_prepare
        self.assertEqual(xenium_prepare.ROOT, simulation_workflow.SOURCE)
        self.assertEqual(xenium_prepare.RESOURCES, vtx.ENGINE_DIR / "resources")

    def test_pyproject_pins_equal_engine_requirements(self):
        req = [l.strip() for l in (vtx.ENGINE_DIR / "requirements-colab.txt").read_text().splitlines()
               if l.strip() and not l.startswith("#")]
        project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
        pins = project["optional-dependencies"]["legacy"]
        self.assertEqual(pins, req)

    def test_package_import_needs_only_the_standard_library(self):
        code = ("import sys, virtualtissue_xenium, virtualtissue_xenium.colab, virtualtissue_xenium.download; "
                "heavy={'numpy','pandas','scipy','h5py','anndata'} & set(sys.modules); "
                "assert not heavy, heavy; print('ok')")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stderr)


if __name__ == "__main__":
    unittest.main()
