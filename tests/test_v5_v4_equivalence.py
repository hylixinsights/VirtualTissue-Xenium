"""Without a v5 evidence policy, v5 reproduces the scientific outputs of v4.0.0.

Compare candidate outputs to immutable v4 source under the same interpreter.
The historical Mac hash snapshot is retained, but is not a portable oracle:
platform math libraries can change serialized floats and content-seeded mocks.
Full byte/hash equality is required; no numerical tolerance or fixture rewrite.
"""
import json
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_default_outputs_match_pinned_v4_in_the_same_environment():
    result = subprocess.run([sys.executable, str(ROOT / "tools/v4_equivalence.py"), "compare", "--same-platform-v4"],
                            cwd=ROOT, env=dict(os.environ, PYTHONPATH=str(ROOT / "src")),
                            capture_output=True, text=True, timeout=360)
    assert result.returncode == 0, result.stdout + result.stderr
    comparison = json.loads(result.stdout)
    assert comparison == dict(identical=True, compared=["demo270", "object96"], not_compared=[], differences=[])
