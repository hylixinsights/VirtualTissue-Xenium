#!/usr/bin/env python3
"""Record or compare hashes of scientific outputs for fixed scenarios.

Run on the v4.0.0 tag to write a reference, then on a later checkout to confirm
that frames, expanded decisions, configurations, prepared fingerprints and Jev
request payloads are unchanged when no v5 evidence policy is declared:

  python tools/v4_equivalence.py record  REFERENCE.json
  python tools/v4_equivalence.py compare REFERENCE.json

For strict compatibility across operating systems/Python builds, compare with
the immutable v4 commit under this same interpreter (requires Git history):

  python tools/v4_equivalence.py compare --same-platform-v4

Historical hash snapshots are environment-specific: last-bit libm differences
can change serialized float hashes and therefore the content-seeded mock answers.

Optional real-data scenario (never packaged):
  --h5ad ROI.h5ad --context CONTEXT.json --lineage-map MAP.csv \
  --cell-type-key first_type --layer purified_raw

No network requests are made. The Jev scenarios use a deterministic in-process
transport, so they test payload construction, pooling and journaling only.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
V4_REFERENCE_COMMIT = "20098224d06fddcc00f2683194622a60febb2b25"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples" / "v3"))

from virtualtissue_xenium.v2.controls import apply_control  # noqa: E402
from virtualtissue_xenium.v2.inputs import digest, load_cells, prepare, read_json  # noqa: E402
from virtualtissue_xenium.v2.provider import JevProvider  # noqa: E402
from virtualtissue_xenium.v2.simulation import run  # noqa: E402
from virtualtissue_xenium.v3.expression_audit import audit_expression  # noqa: E402


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def decisions_sha256(run_dir):
    """Hash decisions in the inline v4 form, whatever the storage format."""
    try:
        from virtualtissue_xenium.v2.simulation import read_decisions
    except ImportError:  # v4: decisions are stored inline
        return file_sha256(Path(run_dir) / "decisions.jsonl")
    h = hashlib.sha256()
    for decision in read_decisions(run_dir):
        h.update((json.dumps(decision, allow_nan=False) + "\n").encode())
    return h.hexdigest()


def mock_transport(body):
    """Deterministic pseudo-answers derived from the request content only."""
    answers = {}
    for qid, question in body["questions"].items():
        actions = sorted(question["criteria"])
        seed = digest([qid, question])
        weights = [1 + int(seed[i * 4:i * 4 + 4], 16) % 97 for i in range(len(actions))]
        total = sum(weights)
        probabilities = {a: w / total for a, w in zip(actions, weights)}
        choice = max(actions, key=lambda a: (probabilities[a], a))
        answers[qid] = dict(type="choice", choice=choice, probabilities=probabilities, confidence=0.5)
    return dict(model=body["model"], answers=answers, usage=dict(input_tokens=11, output_tokens=7))


def run_record(bundle, directory, *, steps, seed, mode="local", prior_mode="shown", jev_weight=.5):
    directory = Path(directory)
    provider = None
    if mode == "jev":
        identity = dict(bundle_fingerprint=bundle["fingerprint"], steps=steps, seed=seed, prior=prior_mode)
        provider = JevProvider("", directory / "jev_journal", identity, transport=mock_transport,
                               prior_mode=prior_mode, schema_version=bundle["schema_version"], max_attempts=100000)
    report = run(bundle, directory, steps=steps, seed=seed, mode=mode, provider=provider, jev_weight=jev_weight)
    record = dict(completed=report["completed"], config_hash=report["config_hash"],
                  config_sha256=file_sha256(directory / "config.json"),
                  frames_sha256=file_sha256(directory / "frames.jsonl"),
                  guidebook_sha256=file_sha256(directory / "guidebook.json"),
                  decisions_inline_sha256=decisions_sha256(directory),
                  events=report["events"], living=report["living"])
    if mode == "jev":
        intents = sorted(p.name for p in (directory / "jev_journal").glob("*.intent.json"))
        record.update(jev_requests=len(intents), jev_request_hashes_sha256=digest(intents),
                      jev_prior_divergence=report["jev_prior_divergence"]["by_lineage"])
    return record


def scenarios(args, work):
    from make_demo import create_demo
    from make_object_demo import create_object
    result = {}
    demo = create_demo(work / "demo")
    bundle = read_json(demo / "control.json")
    result["demo270"] = dict(
        fingerprint=bundle["fingerprint"],
        beta_zero_fingerprint=read_json(demo / "beta_zero.json")["fingerprint"],
        within_type_shuffle_fingerprint=read_json(demo / "shuffle_expression_within_type.json")["fingerprint"],
        global_shuffle_fingerprint=apply_control(bundle, "shuffle_expression_global", seed=5)["fingerprint"],
        label_permutation_fingerprint=apply_control(bundle, "label_permutation", seed=5)["fingerprint"],
        coordinate_shuffle_fingerprint=apply_control(bundle, "coordinate_shuffle", seed=5)["fingerprint"],
        expression_audit_fingerprint=audit_expression(bundle, seed=7)["fingerprint"],
        local=run_record(bundle, work / "demo_local", steps=12, seed=3),
        jev_shown=run_record(bundle, work / "demo_jev_shown", steps=3, seed=3, mode="jev"),
        jev_hidden=run_record(bundle, work / "demo_jev_hidden", steps=3, seed=3, mode="jev", prior_mode="hidden"))
    obj = create_object()
    polygon = prepare(obj["cells"], obj["context"], assayed_genes=obj["assayed_genes"],
                      measurement=obj.get("measurement"), annotation_provenance=obj.get("annotation_provenance"))
    result["object96"] = dict(
        fingerprint=polygon["fingerprint"],
        expression_audit_fingerprint=audit_expression(polygon, seed=7)["fingerprint"],
        local=run_record(polygon, work / "object_local", steps=8, seed=5),
        jev_shown=run_record(polygon, work / "object_jev_shown", steps=2, seed=5, mode="jev"))
    if args.h5ad:
        from virtualtissue_xenium.v2.measurement import read_lineage_map
        cells, genes = load_cells(args.h5ad, cell_type_key=args.cell_type_key, layer=args.layer,
                                  lineage_map=read_lineage_map(args.lineage_map), max_cells=args.max_cells)
        real = prepare(cells, read_json(args.context), assayed_genes=genes, lineage_map=read_lineage_map(args.lineage_map))
        result["real_h5ad"] = dict(input_sha256=file_sha256(args.h5ad), fingerprint=real["fingerprint"],
                                   local=run_record(real, work / "real_local", steps=6, seed=4),
                                   jev_hidden=run_record(real, work / "real_jev", steps=1, seed=4, mode="jev", prior_mode="hidden"))
    return result


def same_platform_reference(args, work):
    """Run the same harness against immutable v4 source with this interpreter.

    No snapshot is regenerated from the candidate code. Missing history fails
    explicitly; numeric/hash equality is never relaxed to accommodate a host.
    """
    archive = subprocess.run(["git", "archive", "--format=tar", V4_REFERENCE_COMMIT],
                             cwd=ROOT, capture_output=True)
    if archive.returncode:
        raise ValueError("The pinned v4 baseline is unavailable. Use a Git clone with full history "
                         "(git fetch --unshallow for a shallow clone), then retry.")
    source = work / "v4_source"
    source.mkdir()
    with tarfile.open(fileobj=io.BytesIO(archive.stdout)) as saved:
        saved.extractall(source, filter="data")
    harness = source / "tools" / "v4_equivalence.py"
    harness.write_bytes(Path(__file__).read_bytes())
    reference = work / "v4_same_platform.json"
    command = [sys.executable, str(harness), "record", str(reference)]
    if args.h5ad:
        for key in ("h5ad", "context", "lineage_map"):
            command.extend(["--" + key.replace("_", "-"), str(Path(getattr(args, key)).resolve())])
        command.extend(["--cell-type-key", args.cell_type_key, "--max-cells", str(args.max_cells)])
        if args.layer:
            command.extend(["--layer", args.layer])
    environment = dict(os.environ, PYTHONPATH=str(source / "src"))
    result = subprocess.run(command, cwd=source, env=environment, capture_output=True, text=True, timeout=300)
    if result.returncode:
        raise ValueError("The pinned v4 baseline harness failed: " + result.stderr[-4000:])
    return json.loads(reference.read_text(encoding="utf8"))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("action", choices=("record", "compare"))
    p.add_argument("reference", type=Path, nargs="?")
    p.add_argument("--same-platform-v4", action="store_true",
                   help="Compare against immutable v4 source under this interpreter; requires full Git history.")
    p.add_argument("--h5ad"); p.add_argument("--context"); p.add_argument("--lineage-map")
    p.add_argument("--cell-type-key", default="cell_type"); p.add_argument("--layer")
    p.add_argument("--max-cells", type=int, default=1000)
    args = p.parse_args(argv)
    if args.same_platform_v4:
        if args.action != "compare" or args.reference is not None:
            p.error("--same-platform-v4 requires compare without a reference file.")
    elif args.reference is None:
        p.error("A reference file is required unless using compare --same-platform-v4.")
    if args.h5ad and not (args.context and args.lineage_map):
        p.error("--h5ad requires --context and --lineage-map.")
    with tempfile.TemporaryDirectory() as temporary:
        observed = scenarios(args, Path(temporary))
        if args.same_platform_v4:
            try:
                reference = same_platform_reference(args, Path(temporary))
            except (ValueError, subprocess.TimeoutExpired) as exc:
                p.error(str(exc))
    if args.action == "record":
        if args.reference.exists():
            p.error("Reference exists; choose a new path.")
        args.reference.write_text(json.dumps(observed, indent=2, sort_keys=True) + "\n", encoding="utf8")
        print(json.dumps(dict(recorded=str(args.reference), scenarios=sorted(observed)), indent=2))
        return 0
    if not args.same_platform_v4:
        reference = json.loads(args.reference.read_text(encoding="utf8"))
    differences = []

    def walk(expected, actual, path):
        if isinstance(expected, dict) and isinstance(actual, dict):
            for key in sorted(set(expected) | set(actual)):
                walk(expected.get(key), actual.get(key), path + [key])
        elif expected != actual:
            differences.append(dict(path="/".join(path), reference=expected, observed=actual))

    for name in sorted(set(reference) & set(observed)):
        walk(reference[name], observed[name], [name])
    skipped = sorted(set(reference) ^ set(observed))
    print(json.dumps(dict(identical=not differences, compared=sorted(set(reference) & set(observed)),
                          not_compared=skipped, differences=differences), indent=2))
    return 1 if differences else 0


if __name__ == "__main__":
    raise SystemExit(main())
