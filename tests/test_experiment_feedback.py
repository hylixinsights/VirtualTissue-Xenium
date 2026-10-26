"""Regression checks for unsupported tissue and experiment interpretation."""
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd

from virtualtissue_xenium import cli, own_data
from virtualtissue_xenium.experiments import experiment_notes


class PreparationScopeTests(unittest.TestCase):
    def sample(self, compartments, qc):
        return dict(obs=pd.DataFrame(dict(compartment=compartments, qc_pass=qc)),
                    output_dir=Path("/tmp/annotation-review"))

    def test_healthy_annotations_and_failed_qc_cannot_open_tumour_workflow(self):
        for sample in [self.sample(["epithelial_other", "immune"], [True, True]),
                       self.sample(["tumor_candidate", "immune"], [False, True])]:
            with self.subTest(sample=sample), self.assertRaisesRegex(own_data.PreparationScopeError, "non-tumour tissue"):
                own_data.require_tumour_interface(sample)
        own_data.require_tumour_interface(self.sample(["tumor_candidate", "immune"], [True, True]))

    def test_cli_stops_before_scoring_and_geometry_without_traceback(self):
        prep = Mock()
        prep.load_sample.return_value = self.sample(["epithelial_other"], [True])
        stderr = io.StringIO()
        with tempfile.TemporaryDirectory() as out, \
             patch.object(cli.vtx, "prepare", prep), \
             patch.object(own_data, "require_annotation_map"), \
             contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as raised:
                cli.main(["prepare", "input", out, "--list-only"])
        self.assertEqual(raised.exception.code, 2)
        prep.score_aucell.assert_not_called()
        prep.read_geometry.assert_not_called()
        self.assertIn("Preparation stopped", stderr.getvalue())
        self.assertIn("Input audit and annotation review", stderr.getvalue())
        self.assertNotIn("Traceback", stderr.getvalue())


class ExperimentNotesTests(unittest.TestCase):
    def test_recognition_is_explicit_and_no_cd8_is_explained(self):
        notes = experiment_notes(dict(scenario="declared_recognition"),
                                 dict(acting_composition={}, growth_headroom=0))
        self.assertIn("not guaranteed", notes[0])
        self.assertIn("same bundle", notes[0])
        self.assertIn("No active CD8", notes[1])
        self.assertIn("0 new bodies", notes[2])

    def test_baseline_does_not_promise_killing_with_longer_time(self):
        notes = experiment_notes(dict(scenario="baseline_unknown"))
        self.assertIn("unknown", notes[0])
        self.assertIn("increasing the duration does not enable it", notes[0])


if __name__ == "__main__":
    unittest.main()
