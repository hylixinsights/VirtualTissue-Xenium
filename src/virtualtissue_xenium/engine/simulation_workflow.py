"""Load a prepared Xenium bundle, configure an experiment, and save a replay."""
import copy
import csv
import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import stat
import tempfile
import zipfile
from collections import Counter
from pathlib import Path, PurePosixPath

from prepared_runtime import load_inputs
from xenium_runtime.io import canonical, fingerprint
from recorded_simulation import RecordedSimulation as AutonomousSimulation
from xenium_runtime.v2.episode import validate_episode
from xenium_runtime.v2.policy import ScoreAwareLocalPolicy
from jev_payload import PROMPT_VERSION, build_batches, payload
from player import write_player
from jev_provider import JevProvider

VERSION = "xenium-jev-colab-1.1.1"
SOURCE = Path(__file__).resolve().parent


def sha256(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def safe_relative(name):
    p = PurePosixPath(name)
    if not name or p.is_absolute() or "\\" in name or any(x in ("..", ".", "") for x in name.rstrip("/").split("/")):
        raise ValueError("Unsafe path in ZIP or manifest.")
    return p


def extract_bundle(zip_path, cache_dir):
    """Extract into a new cache directory; never execute Python shipped in the data ZIP."""
    zip_path = Path(zip_path).expanduser()
    digest = sha256(zip_path)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = Path(tempfile.mkdtemp(prefix="bundle_", dir=cache_dir))
    try:
        with zipfile.ZipFile(zip_path) as z:
            infos = z.infolist()
            if len(infos) > 500 or sum(i.file_size for i in infos) > 1_000_000_000:
                raise ValueError("Bundle exceeds extraction limits.")
            names = set()
            for info in infos:
                rel = safe_relative(info.filename)
                if rel.as_posix() in names or stat.S_ISLNK(info.external_attr >> 16):
                    raise ValueError("Duplicate entry or symbolic link in ZIP.")
                names.add(rel.as_posix())
                path = destination / rel
                if info.is_dir():
                    path.mkdir(parents=True, exist_ok=True)
                else:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with z.open(info) as src, path.open("wb") as dst:
                        shutil.copyfileobj(src, dst)
        manifests = list(destination.rglob("BUNDLE_MANIFEST.json"))
        if len(manifests) != 1:
            raise ValueError("Expected exactly one bundle manifest.")
        folder = manifests[0].parent
        manifest = json.loads(manifests[0].read_text())
        expected_names = {"BUNDLE_MANIFEST.json"}
        for name, expected in manifest["files"].items():
            safe_relative(name)
            expected_names.add(name)
            p = folder / name
            if not p.is_file() or p.stat().st_size != expected["bytes"] or sha256(p) != expected["sha256"]:
                raise ValueError(f"Bundle integrity check failed: {name}")
            if name.endswith(".py") and (not (SOURCE / name).is_file() or sha256(SOURCE / name) != expected["sha256"]):
                raise ValueError("The bundle uses a different runtime. Use the matching simulation notebook; bundle code is never executed.")
        actual = {p.relative_to(folder).as_posix() for p in destination.rglob("*") if p.is_file()}
        if actual != expected_names:
            raise ValueError("Unlisted files in the input bundle.")
        atlas, contract, data = load_inputs(folder)
        return dict(folder=folder, zip_path=zip_path, zip_sha256=digest, atlas=atlas, contract=contract, data=data)
    except BaseException:
        shutil.rmtree(destination)
        raise


def effective_contract(bundle, settings):
    c = copy.deepcopy(bundle["contract"])
    horizon = settings["horizon_hours"] * 60
    if not 10 <= horizon <= c["run"]["max_horizon_min"] or horizon % 10:
        raise ValueError("Horizon must be a multiple of 10 minutes, between 10 minutes and 7 days.")
    if type(settings["max_decisions"]) is not int or settings["max_decisions"] < 1:
        raise ValueError("Decision limit must be a positive integer.")
    if type(settings["seed"]) is not int or settings["seed"] < 0:
        raise ValueError("Seed must be a nonnegative integer.")
    maximum = settings["max_population"]
    if type(maximum) is not int or not len(bundle["atlas"]["cells"]) <= maximum <= 5000:
        raise ValueError("Population ceiling must cover the starting tissue and be at most 5000 bodies.")
    if not re.fullmatch(r"jev-\d+\.\d+\.\d+", settings["model"]):
        raise ValueError("A pinned Jev model version is required.")
    scenario = settings["scenario"]
    if scenario not in ("baseline_unknown", "declared_recognition", "recognition_off"):
        raise ValueError("Unknown scenario.")
    s = copy.deepcopy(c["autonomous_scenarios"]["S0"])
    s["title"] = scenario.replace("_", " ")
    if scenario != "baseline_unknown":
        s.update(recognition=scenario == "declared_recognition", effector_competence=True, clearance_competence=True)
    c["autonomous_scenarios"][scenario] = s
    c["run"]["horizon_min"] = horizon
    c["biology"]["max_population"] = maximum
    c["provider"]["pinned_model"] = settings["model"]
    return c


def create_sim(bundle, settings, policy=None):
    return AutonomousSimulation(bundle["atlas"], effective_contract(bundle, settings),
             scenario=settings["scenario"], seed=settings["seed"], policy=policy or ScoreAwareLocalPolicy(),
             max_decisions=settings["max_decisions"], stop_on_provider_failure=True)


def preflight(bundle, settings):
    sim = create_sim(bundle, settings)
    # Same seeded selection rule used by the kernel if the cell-question cap is smaller than a round.
    rows = sim.requests()
    if len(rows) > settings["max_decisions"]:
        rows = sorted(rows, key=lambda r: fingerprint(dict(seed=sim.seed, time_min=sim.time,
                    cell_id=r["cell_id"], purpose="bounded_pilot_selection")))[:settings["max_decisions"]]
    batches = build_batches(rows, settings["model"], sim.contract["provider"])
    counts = Counter(label for label, acts in zip(sim.tissue.labels, sim.tissue.acts) if acts)
    facts = dict(total_cells=len(bundle["atlas"]["cells"]), acting_cells=int(sim.tissue.acts.sum()),
                 context_cells=int((~sim.tissue.acts).sum()), acting_composition=dict(counts),
                 initial_cell_questions=len(rows), initial_http_batches=len(batches),
                 growth_headroom=settings["max_population"] - len(bundle["atlas"]["cells"]),
                 initially_cycling_or_eligible=sum(bool(s.get("cycle")) for s in sim.bio) +
                     sum("ENTER_CYCLE" in r["options"] for r in rows),
                 prompt_version=PROMPT_VERSION,
                 scenario=sim.scenario, t0_fidelity=sim.fidelity, zip_sha256=bundle["zip_sha256"])
    # A fresh static simulation checks one mechanical tick without asking any provider.
    null = AutonomousSimulation(bundle["atlas"], bundle["contract"], scenario="static", seed=settings["seed"])
    before = null.tissue.xy.copy()
    null.step()
    import numpy as np
    if null.decisions or not np.allclose(before, null.tissue.xy, rtol=0, atol=1e-8):
        raise ValueError("Static preflight failed.")
    facts["static_preflight_passed"] = True
    return facts, payload(batches[0], settings["model"]) if batches else None


def atomic_bytes(path, raw):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as f:
        f.write(raw)
        f.flush()
        os.fsync(f.fileno())
    temporary.replace(path)


def write_json(path, obj):
    atomic_bytes(path, json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False).encode())


def versions():
    packages = ["numpy", "scipy", "pandas", "shapely", "h5py", "anndata", "pyarrow"]
    return dict(python=platform.python_version(), system=platform.system(), machine=platform.machine(),
                packages={p: importlib.metadata.version(p) for p in packages})


def run_context(bundle, settings):
    sources = {p.relative_to(SOURCE).as_posix(): sha256(p) for p in sorted(SOURCE.rglob("*"))
               if p.is_file() and p.suffix in (".py", ".html", ".txt")}
    return dict(version=VERSION, prompt_version=PROMPT_VERSION, zip_sha256=bundle["zip_sha256"], atlas_fingerprint=bundle["atlas"]["fingerprint"],
                original_contract_fingerprint=fingerprint(bundle["contract"]),
                effective_contract_fingerprint=fingerprint(effective_contract(bundle, settings)),
                settings=settings, sources=sources, environment=versions())


def checkpoint(sim, atlas, folder):
    # Refresh an existing final frame, including a provider failure at t=0.
    if sim.frames and sim.frames[-1]["time_min"] == sim.time:
        sim.frames.pop()
    sim.frame()
    episode = sim.episode()
    validate_episode(episode, atlas)
    atomic_bytes(folder / "episode.json.gz", gzip.compress(canonical(episode), mtime=0))
    return episode


def export_tables(episode, folder):
    for name in ("metrics", "events", "decisions"):
        rows = episode[name]
        fields = list(dict.fromkeys(k for row in rows for k in row)) or ["time_min"]
        with (folder / f"{name}.csv").open("w", encoding="utf8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in row.items()})


def run(bundle, settings, output_dir, *, api_key="", max_http_attempts=600, transport=None,
        mode="jev", checkpoint_ticks=6, progress=print):
    """Replay cached answers from t0; only missing, authorized batches reach the network."""
    if mode not in ("jev", "local_preview"):
        raise ValueError("Unknown run mode.")
    if mode == "jev" and transport is None and not api_key:
        raise ValueError("A Jev API key is required for a live run.")
    if type(checkpoint_ticks) is not int or checkpoint_ticks < 1:
        raise ValueError("Checkpoint interval must be a positive number of ticks.")
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    # Exclusive lock prevents two notebook sessions spending from the same journal concurrently.
    lock = folder / ".running.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise RuntimeError("This run is locked by another or interrupted session. Stop any other writer, then use the notebook's unlock control before resuming.") from None
    os.close(fd)
    try:
        context = dict(run_context(bundle, settings), mode=mode,
                       transport="mock" if transport is not None else "live" if mode == "jev" else "local")
        config = folder / "run_config.json"
        if config.exists() and json.loads(config.read_text()) != context:
            raise ValueError("This run folder belongs to another configuration. Keep the original settings to resume, or choose a new RUN_NAME.")
        write_json(config, context)
        write_json(folder / "effective_contract.json", effective_contract(bundle, settings))
        for name in ("SOURCE.json", "aucell_provenance.json", "BUNDLE_MANIFEST.json", "region_selection_rules.json", "score_calibration.json"):
            shutil.copyfile(bundle["folder"] / name, folder / ("input_" + name))
        provider = (JevProvider(api_key, settings["model"], effective_contract(bundle, settings)["provider"],
                    max_http_attempts, folder / "jev_journal.jsonl", context, transport=transport,
                    pause_s=0 if transport else .25,
                    progress=lambda n, cap, hits: progress(f"Jev requests recorded: {n}/{cap}; cached batches reused: {hits}") if n % 10 == 0 else None)
                    if mode == "jev" else ScoreAwareLocalPolicy())
        sim = create_sim(bundle, settings, provider)
        ticks, error = 0, None
        # Preserve an existing valid checkpoint if replay is interrupted before reaching it.
        if not (folder / "episode.json.gz").exists():
            checkpoint(sim, bundle["atlas"], folder)
        old_time = json.loads(gzip.decompress((folder / "episode.json.gz").read_bytes()))["summary"]["time_min"]
        try:
            while not sim.stop_reason:
                sim.step()
                ticks += 1
                if ticks % checkpoint_ticks == 0 and sim.time >= old_time:
                    checkpoint(sim, bundle["atlas"], folder)
                    progress(f"Saved {sim.time / 60:.2f} simulated hours; {len(sim.decisions)} cell decisions.")
            if sim.time >= old_time:
                checkpoint(sim, bundle["atlas"], folder)
        except (KeyboardInterrupt, Exception) as e:
            # A stopped mid-tick state is not a valid checkpoint. Keep the last completed one.
            error = "Interrupted; the last valid checkpoint was preserved." if isinstance(e, KeyboardInterrupt) else f"Simulation stopped ({type(e).__name__}); the last valid checkpoint was preserved."
            progress(error)
        episode = json.loads(gzip.decompress((folder / "episode.json.gz").read_bytes()))
        validate_episode(episode, bundle["atlas"])
        usage = provider.usage() if isinstance(provider, JevProvider) else dict(http_attempts=0)
        report = dict(version=VERSION, mode=mode, test_transport=transport is not None,
                      completed=error is None and sim.stop_reason == "horizon_reached",
                      error=error, provider_error=getattr(provider, "last_error", None),
                      checkpoint_time_min=episode["summary"]["time_min"], stop_reason=sim.stop_reason,
                      usage=usage, max_http_attempts=max_http_attempts if mode == "jev" else 0,
                      summary=episode["summary"], context_fingerprint=fingerprint(context),
                      checkpoint_semantics="Episode is the last complete saved state. Journal may include later paid responses; resume by replaying from t0.")
        write_json(folder / "run_report.json", report)
        export_tables(episode, folder)
        write_player(episode, bundle["atlas"], folder / "replay.html")
        write_json(folder / "OUTPUT_MANIFEST.json", {p.name: sha256(p) for p in sorted(folder.iterdir())
                   if p.is_file() and not p.name.startswith(".") and p.name not in ("OUTPUT_MANIFEST.json",) and not p.name.endswith(".tmp")})
        result_zip = folder.with_name(folder.name + "_results.zip")
        with zipfile.ZipFile(result_zip, "w", compression=zipfile.ZIP_DEFLATED) as z:
            for p in sorted(folder.iterdir()):
                if p.is_file() and not p.name.startswith(".") and not p.name.endswith(".tmp"):
                    z.write(p, folder.name + "/" + p.name)
        progress(f"Results saved to {folder}. Stop reason: {sim.stop_reason or error}.")
        return dict(report=report, folder=folder, zip_path=result_zip, episode=episode)
    finally:
        lock.unlink(missing_ok=True)
