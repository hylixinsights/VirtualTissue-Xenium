"""Compare fixed t0 tumour questions across prompts. Dry run unless --send is set.

python tools/jev_prompt_ablation.py BUNDLE.zip OUT_DIR
python tools/jev_prompt_ablation.py BUNDLE.zip OUT_DIR --send --max-requests 24

No automatic retries. An existing journal (even an interrupted intent) refuses
another send; use a separate directory for an explicitly authorized new experiment.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import re
from pathlib import Path
from statistics import mean
import sys
import urllib.request

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import virtualtissue_xenium as vtx
from jev_payload import payload, v2_payload
from jev_provider import ENDPOINT, NoRedirect, ProviderError, validate_response
from xenium_runtime.io import canonical, fingerprint
from xenium_runtime.v2 import provider_jev as legacy


def variants(tissue):
    def changed(rows, model, *, setting=None, semantics=None, excluded=False):
        body = v2_payload(rows, model)
        if setting is not None:
            body["state"]["setting"] = setting
        if semantics is not None:
            body["state"]["semantics"] = semantics
        if excluded:
            for row in rows:
                body["questions"][f"c{row['cell_id']}"]["instructions"]["local_observation"]["excluded"] = copy.deepcopy(row["observation"].get("excluded", {}))
        return body
    return {
        "legacy_v1": legacy.payload,
        "v2": v2_payload,
        "v2_legacy_setting": lambda r, m: changed(r, m, setting="Independent local cells in an uncalibrated LUAD virtual tissue prototype."),
        "v2_legacy_semantics": lambda r, m: changed(r, m, semantics=legacy.STATE_SEMANTICS),
        "v2_with_excluded": lambda r, m: changed(r, m, excluded=True),
        "v3": lambda r, m: payload(r, m, tissue=tissue),
    }


def canary_rows(bundle, per_quartile, seed):
    settings = dict(horizon_hours=48, max_decisions=100000, seed=seed,
                    max_population=max(1000, len(bundle["atlas"]["cells"])),
                    model="jev-1.14.0", scenario="baseline_unknown")
    rows = [r for r in vtx.workflow.create_sim(bundle, settings).requests()
            if r["observation"]["primary_label"].startswith("Tumor") and "ENTER_CYCLE" in r["options"]]
    scored = [r for r in rows if (r["observation"].get("t0_aucell_percentile") or {}).get("MITOTIC_CELL_CYCLE") is not None]
    picked = []
    for lo, hi in ((0, .25), (.25, .5), (.5, .75), (.75, 1.01)):
        group = sorted((r for r in scored if lo <= r["observation"]["t0_aucell_percentile"]["MITOTIC_CELL_CYCLE"] < hi),
                       key=lambda r: fingerprint(dict(cell_id=r["cell_id"], purpose="ablation_canary", seed=seed)))
        if len(group) < per_quartile:
            raise ValueError(f"Insufficient scored tumour cells in quartile starting at {lo}; reduce --per-quartile.")
        picked += group[:per_quartile]
    return sorted(picked, key=lambda r: r["cell_id"]), len(rows), len(scored)


def groups(rows, model, makers, max_bytes, max_questions):
    """One partition satisfying both size and question limits for every variant."""
    def fits(part):
        return len(part) <= max_questions and all(len(canonical(make(part, model))) <= max_bytes for make in makers)
    out, current = [], []
    for row in rows:
        if not fits([row]):
            raise ValueError("Single question exceeds a provider request limit.")
        if current and not fits(current + [row]):
            out.append(current)
            current = []
        current.append(row)
    return out + ([current] if current else [])


def summarize(results, rows, reference=None):
    mitotic = {f"c{r['cell_id']}": r["observation"]["t0_aucell_percentile"]["MITOTIC_CELL_CYCLE"] for r in rows}
    summary = {}
    for name, per_question in results.items():
        table, pairs = {}, []
        for qid, answers in per_question.items():
            values = [(a.get("reported_probabilities") or a["probabilities"])["ENTER_CYCLE"] for a in answers]
            quartile = f"Q{min(4, int(mitotic[qid] * 4) + 1)}"
            row = table.setdefault(quartile, dict(values=[], chosen=0, repeat_spread=[]))
            row["values"] += values
            row["chosen"] += sum(a["choice"] == "ENTER_CYCLE" for a in answers)
            if len(values) > 1:
                row["repeat_spread"].append(max(values) - min(values))
            if reference and qid[1:] in reference:
                pairs.append(abs(mean(values) - reference[qid[1:]]["p_enter"]))
        summary[name] = dict(quartiles={k: dict(n=len(v["values"]), mean_p_enter=mean(v["values"]), chosen=v["chosen"],
                                  max_repeat_spread=max(v["repeat_spread"]) if v["repeat_spread"] else None)
                                       for k, v in sorted(table.items())},
                             reference_matches=len(pairs),
                             vs_reference_mean_abs_diff=mean(pairs) if pairs else None)
    return summary


def send_plan(bodies, repeats, out, key, max_requests, max_response_bytes, transport=None):
    planned = repeats * sum(len(part) for part in bodies.values())
    if repeats < 1 or not planned or planned > max_requests:
        raise ValueError("Planned requests must be positive and within --max-requests. Nothing sent.")
    def remote(body):
        req = urllib.request.Request(ENDPOINT, data=canonical(body), method="POST",
                                     headers={"Authorization": "Bearer " + key, "Content-Type": "application/json"})
        with urllib.request.build_opener(NoRedirect()).open(req, timeout=120) as response:
            raw = response.read(max_response_bytes + 1)
        if len(raw) > max_response_bytes:
            raise ProviderError("Response exceeded size limit.")
        return json.loads(raw)
    results, sent, failed = {}, 0, False
    # Exclusive creation prevents both concurrent writers and accidental repeated spending.
    with (out / "ablation_journal.jsonl").open("x", encoding="utf8") as journal:
        def record(entry):
            journal.write(canonical(entry).decode() + "\n")
            journal.flush()
            os.fsync(journal.fileno())
        for rep in range(repeats):
            for name, parts in bodies.items():
                for gi, body in enumerate(parts):
                    entry = dict(variant=name, repeat=rep, group=gi, request_hash=fingerprint(body))
                    record(dict(entry, phase="intent", request=body, usage="unknown"))
                    safe = None
                    sent += 1
                    try:
                        response = (transport or remote)(body)
                        whitelisted = {k: response.get(k) for k in ("model", "answers", "usage")} if isinstance(response, dict) else {}
                        raw = canonical(whitelisted).decode()
                        safe = json.loads(raw.replace(key, "[redacted]") if key else raw)
                        answers, repairs = validate_response(safe, body)
                        record(dict(entry, phase="result", status="valid", response=safe, repairs=repairs))
                        for qid, answer in answers.items():
                            results.setdefault(name, {}).setdefault(qid, []).append(answer)
                    except Exception:
                        # Arbitrary exception text may contain credentials or HTTP headers.
                        record(dict(entry, phase="result", status="failed", response=safe,
                                    error="Transport or validation failed; no retry made."))
                        failed = True
                    print(f"{sent}/{planned} {name} repeat {rep + 1} group {gi + 1}: {'failed' if failed else 'valid'}")
                    if failed:
                        return results, sent, failed
    return results, sent, failed


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bundle_zip"); ap.add_argument("out_dir")
    ap.add_argument("--model", default="jev-1.14.0")
    ap.add_argument("--per-quartile", type=int, default=5)
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--variants", nargs="+", choices=list(variants(None)))
    ap.add_argument("--context-mode", choices=["none", "custom", "bundle"],
                    help="Default: none; --context/--tissue alone selects custom. Bundle context is opt-in.")
    ap.add_argument("--context", "--tissue", dest="tissue", help="User-supplied context for v3.")
    ap.add_argument("--reference", help="JSON containing t0_tumour from a saved legacy run.")
    ap.add_argument("--seed", type=int, default=21)
    ap.add_argument("--send", action="store_true")
    ap.add_argument("--max-requests", type=int, default=24)
    a = ap.parse_args(argv)
    if not re.fullmatch(r"jev-\d+\.\d+\.\d+", a.model):
        ap.error("Use a pinned model, such as jev-1.14.0.")
    if min(a.per_quartile, a.repeats, a.max_requests) < 1 or a.seed < 0:
        ap.error("Counts must be positive; seed must be nonnegative.")
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if (out / "ablation_journal.jsonl").exists():
        ap.error("An ablation journal already exists. Preserve it and use a new directory for a new experiment.")
    bundle = vtx.workflow.extract_bundle(a.bundle_zip, out / "cache")
    context = vtx.workflow.prompt_context(bundle, dict(model=a.model, tissue=a.tissue, context_mode=a.context_mode))
    makers = variants(context["tissue"])
    if a.variants:
        makers = {k: makers[k] for k in dict.fromkeys(a.variants)}
    rows, n_tumour, n_scored = canary_rows(bundle, a.per_quartile, a.seed)
    limits = bundle["contract"]["provider"]
    parts = groups(rows, a.model, list(makers.values()), limits["max_payload_bytes"], limits["max_questions"])
    bodies = {name: [make(group, a.model) for group in parts] for name, make in makers.items()}
    planned = len(makers) * len(parts) * a.repeats
    reference = json.loads(Path(a.reference).read_text())["t0_tumour"] if a.reference else None
    vtx.workflow.write_json(out / "payloads.json", bodies)
    vtx.workflow.write_json(out / "plan.json", dict(bundle_sha256=bundle["zip_sha256"], prompt_context=context,
        model=a.model, seed=a.seed, repeats=a.repeats, tumour_cells=n_tumour, scored_cells=n_scored,
        canary_ids=[r["cell_id"] for r in rows], planned_requests=planned, request_cap=a.max_requests,
        payloads_sha256=fingerprint(bodies), reference_sha256=vtx.workflow.sha256(a.reference) if a.reference else None))
    print(f"Canary: {len(rows)} cells, {len(parts)} fixed groups; {planned} requests (cap {a.max_requests}).")
    if not a.send:
        print(f"Dry run: {out / 'payloads.json'}; no request sent.")
        return
    key = os.environ.get("TYPESAFE_API_KEY") or os.environ.get("JEV_API_KEY")
    if not key:
        ap.error("Set TYPESAFE_API_KEY or JEV_API_KEY. Nothing sent.")
    results, sent, failed = send_plan(bodies, a.repeats, out, key, a.max_requests, limits["max_response_bytes"])
    vtx.workflow.write_json(out / "ablation_summary.json", dict(sent=sent, failed=failed, variants=summarize(results, rows, reference)))
    if failed:
        raise SystemExit("A request failed. Remaining requests were not sent; inspect the journal.")


if __name__ == "__main__":
    main()
