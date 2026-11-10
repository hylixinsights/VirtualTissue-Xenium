"""Make the unchanged engine modules and the test helpers importable."""
import os
import sys
from pathlib import Path

import virtualtissue_xenium  # noqa: F401  (puts src/virtualtissue_xenium/engine on sys.path)

TESTS = Path(__file__).resolve().parent
if str(TESTS) not in sys.path:
    sys.path.insert(0, str(TESTS))
os.environ.setdefault("MPLBACKEND", "Agg")

# The legacy physical engine (`vtx`) needs the `legacy` extra. The modern install
# (`.[xenium,test,spatial,analysis]`) lacks the packages below, which these modules
# need at import or run time. Leave a module out of collection when its packages
# are missing, and say so in the header, so a plain `pytest` runs everything else
# instead of aborting. Install `.[legacy,test]` to run them all.
LEGACY_ENGINE_TESTS = {
    "test_preparation.py": ("shapely", "decoupler"),
    "test_release_012.py": ("shapely",),
    "test_release_013.py": ("shapely",),
    "test_engine_integrity.py": ("shapely",),
    "test_experiment_feedback.py": ("shapely",),
    "test_helpers.py": ("shapely",),
}


def _importable(name):
    try:
        __import__(name)
    except ImportError:
        return False
    return True


_MISSING = {name: [m for m in modules if not _importable(m)] for name, modules in LEGACY_ENGINE_TESTS.items()}
collect_ignore = sorted(name for name, missing in _MISSING.items() if missing)


def pytest_report_header(config):
    if collect_ignore:
        missing = sorted({m for name in collect_ignore for m in _MISSING[name]})
        return (f"legacy engine tests not collected (missing {', '.join(missing)}): {', '.join(collect_ignore)}; "
                "install '.[legacy,test]' to run them")
