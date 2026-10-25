"""Small helpers shared by the notebooks (Google Colab or local Jupyter).

Standard library only. Nothing here changes how the engine simulates; these
functions locate the API key, keep a copy of a run in Google Drive, and display
the replay.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

#: Secret / environment-variable names accepted for the Jev (TypeSafe) API key, in order.
KEY_NAMES = ("TYPESAFE_API_KEY", "JEV_API_KEY")


def in_colab() -> bool:
    """True when running inside Google Colab."""
    if "google.colab" in sys.modules or "COLAB_RELEASE_TAG" in os.environ:
        return True
    try:
        import importlib.util
        return importlib.util.find_spec("google.colab") is not None
    except (ImportError, ValueError):
        return False


def mount_drive(mountpoint: str = "/content/drive") -> Path | None:
    """Mount Google Drive in Colab (asks for authorization once). Returns the MyDrive path."""
    if not in_colab():
        return None
    from google.colab import drive  # type: ignore
    drive.mount(mountpoint)
    return Path(mountpoint) / "MyDrive"


def get_jev_key(names=KEY_NAMES, prompt: bool = True) -> str:
    """Return the Jev API key without printing it.

    Order: Colab Secrets (``TYPESAFE_API_KEY``, then ``JEV_API_KEY``), environment
    variables with the same names, then a hidden prompt. Returns ``""`` if nothing
    is provided.
    """
    key = ""
    if in_colab():
        try:
            from google.colab import userdata  # type: ignore
        except ImportError:
            userdata = None
        for name in names:
            if key or userdata is None:
                break
            try:
                key = userdata.get(name) or ""
            except Exception:  # SecretNotFoundError / NotebookAccessError
                key = ""
    for name in names:
        if not key:
            key = os.environ.get(name, "")
    if not key and prompt:
        from getpass import getpass
        key = getpass("Jev / TypeSafe API key (input hidden): ")
    return (key or "").strip()


def _copy_atomic(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    temporary = dst.with_name(dst.name + ".partial")
    shutil.copyfile(src, temporary)
    temporary.replace(dst)


def _run_files(run_dir: Path):
    for p in sorted(run_dir.iterdir()):
        if p.is_file() and not p.name.startswith(".") and not p.name.endswith((".tmp", ".partial")):
            yield p


def backup_run(run_dir, backup_root) -> Path:
    """Copy a run folder (and its results ZIP, if any) to ``backup_root``.

    Call it only between simulation rounds (the notebook does this right after each
    saved checkpoint and at the end), when the API journal is complete on disk.
    The lock file and temporary files are never copied.
    """
    run_dir, backup_root = Path(run_dir), Path(backup_root)
    target = backup_root / run_dir.name
    for p in _run_files(run_dir):
        _copy_atomic(p, target / p.name)
    result_zip = run_dir.with_name(run_dir.name + "_results.zip")
    if result_zip.is_file():
        _copy_atomic(result_zip, backup_root / result_zip.name)
    return target


def restore_run(backup_root, run_dir) -> bool:
    """Copy a previous backup of this run back to local disk, if the local copy is absent.

    This lets a new Colab runtime reuse the Jev answers present in the backup (the
    engine replays from t0 and sends only batches that are not in that journal).
    Responses recorded after the last successful backup may be lost and billed
    again if the original runtime disappears. The engine
    still refuses to resume if the settings, the input or the software environment
    differ from the original run (``run_config.json``).
    """
    run_dir, source = Path(run_dir), Path(backup_root) / Path(run_dir).name
    if (run_dir / "run_config.json").exists() or not (source / "run_config.json").exists():
        return False
    run_dir.mkdir(parents=True, exist_ok=True)
    for p in _run_files(source):
        _copy_atomic(p, run_dir / p.name)
    return True


def progress_with_backup(run_dir, backup_root=None, printer=print):
    """Progress callback for ``simulation_workflow.run``.

    Prints every message; when the engine reports a saved checkpoint ("Saved ...")
    or the final save ("Results saved ..."), it also refreshes the backup copy.
    """
    def progress(message: str) -> None:
        printer(message)
        if backup_root and (message.startswith("Saved ") or message.startswith("Results saved")):
            try:
                backup_run(run_dir, backup_root)
            except OSError as error:  # a Drive hiccup must not stop the simulation
                printer(f"(Backup to {backup_root} failed this time: {error})")
    return progress


def replay_iframe_html(replay_path, height: int = 850) -> str:
    """HTML for showing a saved replay inline (no external assets, no network)."""
    import html
    text = Path(replay_path).read_text(encoding="utf8")
    return ('<iframe style="width:100%;height:' + str(int(height)) + 'px;border:0" srcdoc="'
            + html.escape(text, quote=True) + '"></iframe>')


def show_replay(replay_path, height: int = 850) -> None:
    """Display ``replay.html`` inside the notebook output."""
    from IPython.display import HTML, display
    display(HTML(replay_iframe_html(replay_path, height)))
