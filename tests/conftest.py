"""Make the unchanged engine modules and the test helpers importable."""
import os
import sys
from pathlib import Path

import virtualtissue_xenium  # noqa: F401  (puts src/virtualtissue_xenium/engine on sys.path)

TESTS = Path(__file__).resolve().parent
if str(TESTS) not in sys.path:
    sys.path.insert(0, str(TESTS))
os.environ.setdefault("MPLBACKEND", "Agg")
