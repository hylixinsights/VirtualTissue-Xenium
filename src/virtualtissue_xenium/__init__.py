"""VirtualTissue-Xenium: AI-agent-driven virtual tissue simulations from Xenium data.

The simulation code lives in the ``engine/`` folder next to this file. The
biological kernel is unchanged; package-side reporting and payload revisions are recorded in
``engine/ORIGINAL_SOURCES.json``. Those modules use
top-level imports (``import xenium_prepare``, ``from xenium_runtime.io import ...``)
because prepared input bundles embed the same files and the simulator compares
them by SHA-256. Importing this package therefore puts ``engine/`` on ``sys.path``.

This module imports only the Python standard library, so it is safe to import
*before* the pinned scientific dependencies are installed (as the Colab
notebooks do).

Quick access (each attribute imports the engine module on first use)::

    import virtualtissue_xenium as vtx
    prep = vtx.prepare     # Part 1: xenium_prepare  (load, annotate, AUCell, regions, export)
    wf = vtx.workflow      # Parts 2-4: simulation_workflow (extract, preflight, run, replay)
    rt = vtx.runtime       # prepared_runtime (load_inputs, initialize, validate_bundle)
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

__version__ = "0.1.3"

#: Folder that holds the engine modules and their recorded resources.
ENGINE_DIR = Path(__file__).resolve().parent / "engine"

#: Folder with metadata of the public example dataset (no data files).
DATASETS_DIR = Path(__file__).resolve().parent / "datasets"

#: Versions written by the engine itself (kept here for documentation).
PREPARATION_VERSION = "lung5k-preparation-1.1.0"
SIMULATION_VERSION = "xenium-jev-colab-1.1.3"

_ENGINE_MODULES = {
    "prepare": "xenium_prepare",
    "workflow": "simulation_workflow",
    "runtime": "prepared_runtime",
    "jev_provider": "jev_provider",
    "player": "player",
}


def activate_engine() -> Path:
    """Put the engine folder first on ``sys.path`` (idempotent) and return it."""
    path = str(ENGINE_DIR)
    if path not in sys.path:
        sys.path.insert(0, path)
    return ENGINE_DIR


def engine_module_origins() -> dict:
    """Return where each already-imported engine module was loaded from.

    Useful to confirm that no other copy of ``xenium_runtime`` or
    ``prepared_runtime`` (for example from an extracted bundle folder) shadows
    the packaged engine.
    """
    names = ["xenium_prepare", "simulation_workflow", "prepared_runtime", "jev_provider",
             "player", "xenium_runtime"]
    return {n: getattr(sys.modules[n], "__file__", None) for n in names if n in sys.modules}


def check_engine_origin() -> None:
    """Raise if an imported engine module does not come from ``ENGINE_DIR``."""
    for name, origin in engine_module_origins().items():
        if origin and ENGINE_DIR not in Path(origin).resolve().parents:
            raise RuntimeError(
                f"Module {name!r} was imported from {origin}, not from {ENGINE_DIR}. "
                "Restart the Python session and import virtualtissue_xenium first.")


def __getattr__(name: str):
    if name in _ENGINE_MODULES:
        module = importlib.import_module(_ENGINE_MODULES[name])
        check_engine_origin()
        return module
    raise AttributeError(f"module 'virtualtissue_xenium' has no attribute {name!r}")


activate_engine()

__all__ = ["ENGINE_DIR", "DATASETS_DIR", "PREPARATION_VERSION", "SIMULATION_VERSION",
           "activate_engine", "engine_module_origins", "check_engine_origin", "__version__",
           *_ENGINE_MODULES]
