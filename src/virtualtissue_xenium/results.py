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
    return dict(folder=folder, report=report, episode=episode, metrics=read("metrics.csv"),
                events=read("events.csv"), decisions=read("decisions.csv"), replay=folder / "replay.html")


def labels_by_body(results) -> dict:
    """Cell id → simulated cell type (daughters inherit their mother's label)."""
    return {b["id"]: b["label"] for b in results["episode"]["bodies"]}


def acting_bodies(results) -> set:
    """Ids of bodies that make decisions (initial actors, daughters, recruited cells)."""
    return {b["id"] for b in results["episode"]["bodies"] if b["origin"] in ACTING_ORIGINS}


def summary(results) -> dict:
    """Headline numbers of one run."""
    report, episode = results["report"], results["episode"]
    final = episode["summary"].get("final", {})
    events = episode["summary"].get("events", {})
    return dict(
        run=results["folder"].name,
        mode=report.get("mode"),
        scenario=episode.get("scenario"),
        provider=episode.get("provider"),
        completed=report.get("completed"),
        stop_reason=report.get("stop_reason") or report.get("error"),
        simulated_hours=round(report.get("checkpoint_time_min", 0) / 60, 2),
        cell_decisions=episode["summary"].get("decisions"),
        invalid_answers=episode["summary"].get("invalid_answers"),
        http_attempts=report.get("usage", {}).get("http_attempts"),
        reused_batches=report.get("usage", {}).get("reused_batches"),
        input_tokens=report.get("usage", {}).get("input_tokens"),
        output_tokens=report.get("usage", {}).get("output_tokens"),
        tumour_cells_final=final.get("tumour_cells"),
        divisions=events.get("division", 0),
        apoptosis_committed=events.get("apoptosis_committed", 0),
        cytotoxic_pulses=events.get("cytotoxic_pulse", 0),
        corpses_cleared=events.get("clearance", 0),
        corpses_degraded=events.get("corpse_degraded", 0),
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
        probabilities = answer.get("probabilities") or {}
        choice = answer.get("choice")
        target = None
        if choice and ":" in choice:
            choice, target = choice.split(":", 1)
        rows.append(dict(hours=d["time_min"] / 60, cell_id=d["cell_id"], cell_type=labels.get(d["cell_id"]),
                         choice=choice, target_id=target, probability=probabilities.get(answer.get("choice")),
                         confidence=answer.get("confidence"), n_options=len(d.get("options") or {}),
                         status=d.get("status"), provider=d.get("provider")))
    return pd.DataFrame(rows)


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
                         probabilities=answer.get("probabilities"), wake_reasons=d.get("wake_reasons"),
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
                        "clearance", "corpse_degraded", "cd4_help") if k in events]
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


def compare_runs(paths):
    """Summary table for several result ZIPs/folders (one row per run)."""
    import pandas as pd

    return pd.DataFrame([summary(load_results(p)) for p in paths]).set_index("run")
