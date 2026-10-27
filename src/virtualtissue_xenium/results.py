"""Part 4 helpers: open a saved simulation result and summarise it.

A result is either the run folder written by ``simulation_workflow.run`` or the
``<RUN_NAME>_results.zip`` downloaded from Colab. Nothing here calls the API or
re-runs the simulation; everything is read from the saved files:

``run_report.json``   stop reason, simulated time, API usage
``episode.json.gz``   complete trajectory (frames every 20 simulated minutes, bodies, decisions, events)
``metrics.csv``       tissue metrics every 10 simulated minutes
``events.csv``        divisions, deaths, cytotoxic pulses, clearance, ...
``decisions.csv``     every cell question with its permitted options and the answer
``replay.html``       the interactive viewer (open in any browser, no network needed)

pandas and matplotlib are imported inside the functions.
"""
from __future__ import annotations

import gzip
import json
import tempfile
import zipfile
import copy
import warnings
from pathlib import Path

ACTING_ORIGINS = ("roi", "born", "recruited")  # same set as xenium_runtime.v2.mechanics


def _find_run_folder(path: Path) -> Path:
    reports = list(path.rglob("run_report.json"))
    if len(reports) != 1:
        raise ValueError(f"Expected exactly one run_report.json under {path}, found {len(reports)}.")
    return reports[0].parent


def load_results(path, extract_to=None) -> dict:
    """Open a results ZIP or run folder.

    Returns a dict with ``folder``, ``report``, ``episode``, ``metrics``, ``events``,
    ``decisions`` (pandas DataFrames) and ``replay`` (path to replay.html).
    """
    import pandas as pd

    path = Path(path)
    if path.is_file() and path.suffix == ".zip":
        target = Path(extract_to) if extract_to else Path(tempfile.mkdtemp(prefix="vtx_results_"))
        with zipfile.ZipFile(path) as z:
            for name in z.namelist():
                if name.startswith("/") or ".." in Path(name).parts:
                    raise ValueError("Unsafe path in results ZIP.")
            z.extractall(target)
        folder = _find_run_folder(target)
    elif path.is_dir():
        folder = path if (path / "run_report.json").exists() else _find_run_folder(path)
    else:
        raise FileNotFoundError(f"Not a results ZIP or run folder: {path}")
    report = json.loads((folder / "run_report.json").read_text())
    episode = json.loads(gzip.decompress((folder / "episode.json.gz").read_bytes()))
    read = lambda name: pd.read_csv(folder / name) if (folder / name).stat().st_size else pd.DataFrame()
    config = json.loads((folder / "run_config.json").read_text()) if (folder / "run_config.json").exists() else {}
    return dict(folder=folder, report=report, episode=episode, config=config, metrics=read("metrics.csv"),
                events=read("events.csv"), decisions=read("decisions.csv"), replay=folder / "replay.html")


def labels_by_body(results) -> dict:
    """Cell id → simulated cell type (daughters inherit their mother's label)."""
    return {b["id"]: b["label"] for b in results["episode"]["bodies"]}


def acting_bodies(results) -> set:
    """Ids of bodies that make decisions (initial actors, daughters, recruited cells)."""
    return {b["id"] for b in results["episode"]["bodies"] if b["origin"] in ACTING_ORIGINS}


def summary(results, at_time_min=None) -> dict:
    """Headline numbers at a saved frame; end-of-run status stays explicit."""
    from recording_analysis import event_totals, population_limit
    report, episode = results["report"], results["episode"]
    time = episode["summary"]["time_min"] if at_time_min is None else at_time_min
    frames = [f for f in episode["frames"] if f["time_min"] == time]
    if not frames:
        raise ValueError(f"No saved frame at {time} minutes; discrete events are not interpolated.")
    frame = frames[-1]
    labels, actors = labels_by_body(results), acting_bodies(results)
    tumour = sum(cid in actors and labels[cid].startswith("Tumor") and s["viability"] == "alive"
                 for cid, s in zip(frame["ids"], frame["states"]))
    limit = population_limit(episode, time)
    return dict(
        run=results["folder"].name,
        mode=report.get("mode"),
        scenario=episode.get("scenario"),
        provider=episode.get("provider"),
        completed=report.get("completed"),
        stop_reason=report.get("stop_reason") or report.get("error"),
        simulated_hours=round(time / 60, 4),
        reached_hours=round(episode["summary"]["time_min"] / 60, 4),
        requested_hours=episode.get("contract", {}).get("run", {}).get("horizon_min", 0) / 60,
        cell_decisions=sum(d["time_min"] <= time for d in episode.get("decisions", [])),
        cell_decisions_whole_run=episode["summary"].get("decisions"),
        usage_scope="HTTP usage and completion status describe the whole run; outcomes use simulated_hours.",
        invalid_answers=episode["summary"].get("invalid_answers"),
        http_attempts=report.get("usage", {}).get("http_attempts"),
        reused_batches=report.get("usage", {}).get("reused_batches"),
        input_tokens=report.get("usage", {}).get("input_tokens"),
        output_tokens=report.get("usage", {}).get("output_tokens"),
        tumour_cells_final=tumour,
        population_ceiling=limit["maximum"],
        ceiling_reached_hours=None if limit["reached_min"] is None else limit["reached_min"] / 60,
        **event_totals(episode, time),
    )


def population_by_type(results, include_context: bool = False):
    """Living cells per type at each saved frame (rows: simulated hours).

    By default only acting cells are counted (context/halo bodies are excluded).
    """
    import pandas as pd

    labels = labels_by_body(results)
    acting = acting_bodies(results)
    rows = []
    for frame in results["episode"]["frames"]:
        counts = {}
        for cid, state in zip(frame["ids"], frame["states"]):
            label = labels.get(cid, "unknown")
            if (not include_context and cid not in acting) or state.get("viability") != "alive":
                continue
            counts[label] = counts.get(label, 0) + 1
        rows.append(dict(hours=frame["time_min"] / 60, **counts))
    return pd.DataFrame(rows).set_index("hours").fillna(0).astype(int)


def decision_table(results):
    """One row per answered cell question: time, cell, type, choice, probability."""
    import pandas as pd

    labels = labels_by_body(results)
    rows = []
    for d in results["episode"]["decisions"]:
        answer = d.get("answer") or {}
        probabilities = answer.get("reported_probabilities") or answer.get("probabilities") or {}
        choice = answer.get("choice")
        target = None
        if choice and ":" in choice:
            choice, target = choice.split(":", 1)
        rows.append(dict(hours=d["time_min"] / 60, cell_id=d["cell_id"], cell_type=labels.get(d["cell_id"]),
                         choice=choice, target_id=target, probability=probabilities.get(answer.get("choice")),
                         confidence=answer.get("confidence"), n_options=len(d.get("options") or {}),
                         status=d.get("status"), provider=d.get("provider"),
                         probability_source="reported" if "reported_probabilities" in answer else "validated"))
    return pd.DataFrame(rows)


def tumour_cycle_diagnostics(results, at_time_min=None):
    """P(ENTER_CYCLE), choices and the binary decision margin, from saved results."""
    from recording_analysis import tumour_cycle_diagnostics as diagnose
    return diagnose(results["episode"], at_time_min)


def choice_counts(results):
    """Cross-table: cell type × chosen action."""
    import pandas as pd

    table = decision_table(results)
    if table.empty:
        return pd.DataFrame()
    return pd.crosstab(table.cell_type, table.choice)


def cell_history(results, cell_id: int):
    """Every decision of one cell, with the full probability distribution."""
    import pandas as pd

    rows = []
    for d in results["episode"]["decisions"]:
        if d["cell_id"] != int(cell_id):
            continue
        answer = d.get("answer") or {}
        rows.append(dict(hours=d["time_min"] / 60, choice=answer.get("choice"), confidence=answer.get("confidence"),
                         probabilities=answer.get("reported_probabilities") or answer.get("probabilities"),
                         probability_repairs=answer.get("probability_repairs", []), wake_reasons=d.get("wake_reasons"),
                         status=d.get("status")))
    return pd.DataFrame(rows)


def event_counts(results, bin_hours: float = 1.0):
    """Events per kind in time bins (rows: start of each bin, in simulated hours)."""
    import pandas as pd

    events = pd.DataFrame(results["episode"]["events"])
    if events.empty:
        return pd.DataFrame()
    events["bin"] = (events.time_min / 60 // bin_hours) * bin_hours
    return pd.crosstab(events.bin, events.kind)


def plot_overview(results, figsize=(14, 9)):
    """Four panels: living cells per type, cumulative decisions, events, top choices per type."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    population = population_by_type(results)
    if not population.empty:
        population.plot(ax=axes[0, 0], linewidth=1.8)
    axes[0, 0].set(title="Living acting cells per type", xlabel="Simulated hours", ylabel="Cells")
    axes[0, 0].legend(fontsize=8, frameon=False, ncol=2)

    metrics = results["metrics"]
    if not metrics.empty and "decisions" in metrics:
        axes[0, 1].plot(metrics.time_min / 60, metrics.decisions, color="#397f9c", linewidth=2)
    axes[0, 1].set(title="Cumulative cell decisions", xlabel="Simulated hours", ylabel="Decisions")

    events = event_counts(results)
    keep = [k for k in ("division", "cycle_entered", "cytotoxic_pulse", "apoptosis_committed",
                        "hypoxic_necrosis", "clearance", "corpse_degraded", "cd4_help",
                        "division_blocked_population_ceiling", "division_blocked_local_pressure",
                        "division_blocked_compression", "division_blocked_crowding", "growth_blocked_compression") if k in events]
    if keep:
        events[keep].cumsum().plot(ax=axes[1, 0], linewidth=1.8, marker="o", markersize=3)
        axes[1, 0].legend(fontsize=8, frameon=False)
    axes[1, 0].set(title="Cumulative events", xlabel="Simulated hours", ylabel="Events")

    counts = choice_counts(results)
    if not counts.empty:
        top = counts.sum().sort_values(ascending=False).head(10).index
        counts[top].plot(kind="barh", stacked=True, ax=axes[1, 1], width=.75, colormap="tab20")
        axes[1, 1].legend(fontsize=7, frameon=False, loc="lower right")
    axes[1, 1].set(title="Choices by cell type (top 10 actions)", xlabel="Decisions", ylabel="")
    fig.suptitle(f"{results['folder'].name} · stop reason: {results['report'].get('stop_reason')}", fontsize=12)
    fig.tight_layout()
    return fig


def _comparison_identity(result):
    """All model assumptions must match except the declared recognition switch."""
    episode, config = result["episode"], result.get("config", {})
    contract = copy.deepcopy(episode.get("contract", {}))
    scenarios = contract.pop("autonomous_scenarios", {})
    scenario = copy.deepcopy(scenarios.get(episode.get("scenario"), {}))
    for key in ("title", "recognition"):
        scenario.pop(key, None)
    if not scenario or not contract or not config.get("zip_sha256"):
        raise ValueError("Comparison provenance is incomplete: a contract, active scenario and bundle SHA are required.")
    return dict(bundle=config["zip_sha256"], atlas=episode.get("atlas_fingerprint"),
                seed=episode.get("seed"), provider=episode.get("provider"),
                kernel=episode.get("engine_fingerprint"), prompt=episode.get("prompt_version", "legacy-v1"),
                prompt_context=config.get("prompt_context") or episode.get("prompt_context"),
                reporting=episode.get("reporting_version", "legacy-v1"),
                software_sources=config.get("sources"), environment=config.get("environment"),
                decision_limit=episode.get("summary", {}).get("decision_budget", {}).get("limit"),
                contract=contract, active_scenario=scenario)


def compare_runs(paths):
    """Validate matched runs and report both at their latest common saved time.

    Supports identical scenarios (sham repeats) or the declared-recognition /
    recognition-off pair. Replicates with different seeds need a separate analysis.
    The original recordings are never truncated or overwritten.
    """
    import pandas as pd
    runs = [load_results(p) for p in paths]
    if len(runs) < 2:
        raise ValueError("Choose at least two matched runs.")
    scenarios = {r["episode"]["scenario"] for r in runs}
    if len(scenarios) > 1 and scenarios != {"declared_recognition", "recognition_off"}:
        raise ValueError("Only declared_recognition and recognition_off form a supported matched pair.")
    identities = [_comparison_identity(r) for r in runs]
    differences = sorted({k for identity in identities[1:] for k in identities[0] if identity[k] != identities[0][k]})
    if differences:
        raise ValueError("Runs are not matched: " + ", ".join(differences))
    common = set.intersection(*({f["time_min"] for f in r["episode"]["frames"]} for r in runs))
    if not common:
        raise ValueError("Runs have no common saved frame.")
    time = max(common)
    messages = []
    if any(r["episode"]["summary"]["time_min"] > time for r in runs):
        messages.append(f"Comparison restricted to the common saved time: {time / 60:g} h.")
    if time == 0:
        messages.append("Only t0 is shared; no trajectory effect can be compared.")
    if identities[0]["prompt"] == "legacy-v1" and "jev" in str(identities[0]["provider"]):
        messages.append("Legacy Jev prompts expose the recognition scenario to all cells; this comparison is potentially confounded.")
    table = pd.DataFrame([summary(r, time) for r in runs]).set_index("run")
    if table.ceiling_reached_hours.notna().any():
        messages.append("Population ceiling reached: growth is constrained by a computational limit.")
    table["comparison_time_hours"] = time / 60
    table["comparison_notes"] = " ".join(messages) or "Matched settings; common saved time."
    for message in messages:
        warnings.warn(message, UserWarning, stacklevel=2)
    table.attrs.update(common_time_min=time, warnings=messages)
    return table


def export_comparison_replays(paths, output_dir):
    """Write copies on shared field scales; comparison validation is mandatory."""
    from player import write_payload
    paths = list(paths)
    table = compare_runs(paths)
    runs = [load_results(p) for p in paths]
    payloads = []
    for result in runs:
        # A saved HTML contains polygons and pathway profiles needed by the viewer.
        import base64
        import re
        html = result["replay"].read_text()
        match = re.search(r'<script id="payload"[^>]*>\s*([A-Za-z0-9+/=]+)\s*</script>', html)
        if not match:
            raise ValueError("Could not read the saved replay payload.")
        original = json.loads(gzip.decompress(base64.b64decode(match[1])))
        # Re-export fields from the complete episode without needing the input atlas.
        payloads.append((result, original))
    from player import refresh_payload
    peaks = {}
    for result, _ in payloads:
        for frame in result["episode"]["frames"]:
            for key, values in frame.get("fields", {}).items():
                peaks[key] = max(peaks.get(key, 0), max(values, default=0))
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    destinations = []
    for index, (result, original) in enumerate(payloads, 1):
        destination = target / f"{index}_{result['folder'].name}_shared_scales.html"
        updated = refresh_payload(original, result["episode"], scales=peaks)
        updated["comparison_end_min"] = table.attrs["common_time_min"]
        updated["comparison_notes"] = table["comparison_notes"].iloc[0]
        write_payload(updated, destination)
        destinations.append(destination)
    return destinations
