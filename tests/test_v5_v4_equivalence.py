"""Without a v5 evidence policy, v5 reproduces the scientific outputs of v4.0.0.

``tests/fixtures/v4_0_0_reference.json`` was recorded by ``tools/v4_equivalence.py
record`` on the v4.0.0 tag (commit 20098224). It pins, for a 270-cell synthetic
bundle and a 96-cell polygon object: prepared and control fingerprints,
expression-audit fingerprints, ``frames.jsonl``, ``config.json``, ``guidebook.json``,
the decision log in its inline form, and the hashes of every Jev request built
with the prior shown and hidden. A failure here means default behaviour changed.
"""
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_default_outputs_match_the_recorded_v4_reference(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("v4_equivalence", ROOT / "tools/v4_equivalence.py")
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setattr(sys, "path", list(sys.path))
    spec.loader.exec_module(module)
    reference = json.loads((ROOT / "tests/fixtures/v4_0_0_reference.json").read_text(encoding="utf8"))
    args = module.argparse.Namespace(h5ad=None)
    observed = module.scenarios(args, tmp_path)
    assert sorted(observed) == sorted(reference) == ["demo270", "object96"]
    assert observed == reference
