"""Interpreter check and restart detection of both setup bootstraps.

Adapted from the original Part 1 test (Xenium_Public_Lung5K/processing/tests/test_colab_setup.py);
the test of the old embedded-code notebook was replaced by checks of the new notebooks.
"""
import ast
import contextlib
import io
import json
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import colab_setup
import colab_setup_prepare
import virtualtissue_xenium as vtx

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ["Part1_Prepare_Xenium_Region.ipynb", "Part2-4_Run_Jev_Simulation.ipynb", "Part4_Explore_Saved_Results.ipynb"]


class SetupTests(unittest.TestCase):
    def test_both_current_interpreters_are_supported(self):
        for module in (colab_setup, colab_setup_prepare):
            self.assertEqual(module.check_python((3, 12)), (3, 12))
            self.assertEqual(module.check_python((3, 13)), (3, 13))

    def test_unsupported_interpreter_reports_actual_version(self):
        for module in (colab_setup, colab_setup_prepare):
            with self.assertRaisesRegex(RuntimeError, "running Python 3.10"):
                module.check_python((3, 10))

    def test_changed_loaded_binary_needs_restart_before_import(self):
        for module in (colab_setup, colab_setup_prepare):
            fake_numpy = types.SimpleNamespace(__version__="1.26.4")
            with patch.dict(sys.modules, {"numpy": fake_numpy}), patch.object(module.metadata, "version", return_value="2.3.5"), \
                    patch.object(module.subprocess, "run") as probe, contextlib.redirect_stdout(io.StringIO()) as message:
                result = module.setup_environment(vtx.ENGINE_DIR, install_packages=False)
            self.assertFalse(result["ready"])
            self.assertTrue(result["restart_required"])
            self.assertIn("Restart session", message.getvalue())
            probe.assert_not_called()

    def test_requirements_file_is_where_setup_expects_it(self):
        self.assertTrue((vtx.ENGINE_DIR / "requirements-colab.txt").is_file())

    def test_current_simulation_does_not_trigger_spurious_restart(self):
        with patch.dict(sys.modules, {"simulation_workflow": types.SimpleNamespace(VERSION=vtx.SIMULATION_VERSION)}), \
             patch.object(colab_setup, "changed_loaded_packages", return_value=[]), \
             patch.object(colab_setup.subprocess, "run", return_value=types.SimpleNamespace(returncode=0, stdout='{"python":"3.13"}', stderr="")):
            self.assertTrue(colab_setup.setup_environment(vtx.ENGINE_DIR, install_packages=False)["ready"])

    def test_loaded_012_workflow_requires_restart_after_upgrade(self):
        with patch.dict(sys.modules, {"simulation_workflow": types.SimpleNamespace(VERSION="xenium-jev-colab-1.1.1")}), \
             patch.object(colab_setup, "changed_loaded_packages", return_value=[]), \
             patch.object(colab_setup.subprocess, "run") as probe:
            result = colab_setup.setup_environment(vtx.ENGINE_DIR, install_packages=False)
        self.assertFalse(result["ready"])
        self.assertIn("previous simulation code", result["changed_packages"])
        probe.assert_not_called()

    def test_notebooks_are_valid_english_and_have_edit_markers(self):
        for name in NOTEBOOKS:
            nb = json.loads((ROOT / "notebooks" / name).read_text())
            text = "\n".join("".join(c["source"]) for c in nb["cells"])
            for phrase in ["Preparar", "Escolha", "Células", "região", "simulação", "Resultados", "Próxima"]:
                self.assertNotIn(phrase, text)
            self.assertIn("✏️", text)
            self.assertIn("Running Python:", text)
            for c in nb["cells"]:
                if c["cell_type"] == "code":
                    ast.parse("".join(c["source"]))
            # Secrets are read from Colab Secrets / environment, never written in a cell.
            self.assertNotRegex(text, r"(?i)api_key\s*=\s*['\"][A-Za-z0-9]{12,}")

    def test_part1_region_choice_is_an_explicit_parameter(self):
        nb = json.loads((ROOT / "notebooks" / NOTEBOOKS[0]).read_text())
        text = "\n".join("".join(c["source"]) for c in nb["cells"])
        self.assertIn('REGION_TO_EXPORT = "R01"', text)
        self.assertIn("DATA_DIR = ", text)
        self.assertIn("OUTPUT_BASE = ", text)


if __name__ == "__main__":
    unittest.main()
