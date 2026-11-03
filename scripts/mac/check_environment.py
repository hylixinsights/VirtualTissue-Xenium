"""Offline environment and packaged-resource checks for the Mac launchers."""
import hashlib
import json
from pathlib import Path
import sys

import virtualtissue_xenium as vtx
from virtualtissue_xenium.v3.ontology import load_ontology, load_signatures


def main():
    repo = Path(sys.argv[1]).resolve()
    origin = Path(vtx.__file__).resolve()
    if origin != (repo / "src/virtualtissue_xenium/__init__.py").resolve():
        raise RuntimeError(f"Code imported from a different installation: {origin}")
    manifest = json.loads((vtx.ENGINE_DIR / "ORIGINAL_SOURCES.json").read_text())["files"]
    present = {p.relative_to(vtx.ENGINE_DIR).as_posix() for p in vtx.ENGINE_DIR.rglob("*")
               if p.is_file() and "__pycache__" not in p.parts and p.name != "ORIGINAL_SOURCES.json"}
    if present != set(manifest):
        raise RuntimeError("Engine files differ from the manifest.")
    for name, expected in manifest.items():
        if hashlib.sha256((vtx.ENGINE_DIR / name).read_bytes()).hexdigest() != expected["sha256"]:
            raise RuntimeError(f"Engine hash mismatch: {name}")
    ontology = load_ontology()
    counts = {species: len(load_signatures(species)[0]) for species in ("human", "mouse")}
    if len(ontology["lineages"]) != 22 or any(n != 14 for n in counts.values()):
        raise RuntimeError("This version requires 22 lineages and 14 signatures per species.")
    print(json.dumps(dict(status="ok", python=sys.version.split()[0], executable=sys.executable,
                         package_version=vtx.__version__, source=str(origin), engine_files_verified=len(manifest),
                         lineages=len(ontology["lineages"]), signatures=counts,
                         interpretation="Engineering verification; does not validate biological predictions."), indent=2))


if __name__ == "__main__":
    main()
