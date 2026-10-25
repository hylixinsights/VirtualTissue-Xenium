"""Build the self-contained interactive replay (one HTML file, no network needed)."""
import base64
import gzip
import json
from pathlib import Path

import numpy as np

TEMPLATE = Path(__file__).resolve().parent / "web" / "player_template.html"

ACTIONS = ["", "CYCLING", "WAIT", "REST", "HOLD_POSITION", "HOLD_CONTACT", "DISENGAGE", "SCAN_RANDOM",
           "MIGRATE_UP_CXCL9", "APPROACH_TUMOR", "DEGRANULATE", "APPROACH_CORPSE", "EFFEROCYTOSE", "PATROL",
           "HELP_CD8", "UPREGULATE_PDL1", "GLYCOLYTIC_SHIFT", "OXIDATIVE_SHIFT", "SECRETE_VEGFA", "SECRETE_CXCL9",
           "SECRETE_VEGFA_TGFB", "SECRETE_TGFB", "DEPOSIT_MATRIX", "SPROUT", "PROLIFERATE", "ENTER_CYCLE", "REMAIN_QUIESCENT", "STAY_QUIESCENT"]
VIABILITY = ["alive", "apoptosis_committed", "dead_intact"]
EVENTS = ["division", "apoptosis_committed", "intact_death", "clearance", "corpse_degraded", "cytotoxic_pulse",
          "cd4_help", "target_reached", "cycle_entered", "hypoxic_necrosis", "exhaustion_transition",
          "helped_cd8_cycle_entered"]
FIELDS = {"ifng": "IFN-γ-like", "cxcl9": "CXCL9/10-like", "oxygen": "Oxygen", "tgfb": "TGF-β-like"}


def _polygon(geometry, ox, oy):
    polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
    return [[[[round(x - ox, 2), round(y - oy, 2)] for x, y in ring] for ring in polygon] for polygon in polygons]



def _u8(values, scale):
    a = np.clip(np.asarray(values, dtype=float) / scale, 0, 1)
    return base64.b64encode((a * 255).round().astype(np.uint8).tobytes()).decode("ascii")


def build_payload(episode, atlas):
    by_id = {c["cell_id"]: c for c in atlas["cells"]}
    frames0 = episode["frames"][0]
    ox, oy = float(min(frames0["x"])), float(min(frames0["y"]))
    bodies = {b["id"]: b for b in episode["bodies"]}

    def root_cell(cid):
        seen = 0
        while cid not in by_id and cid in bodies and bodies[cid]["parent"] is not None and seen < 64:
            cid, seen = bodies[cid]["parent"], seen + 1
        return by_id.get(cid)

    signatures = sorted(atlas["signatures"])
    cells = {}
    for cid, b in bodies.items():
        src = root_cell(cid)
        sig = (src or {}).get("signatures") or {}
        u = [round(float(sig[s]["percentile_within_type_depth"]), 3) if s in sig else None for s in signatures]
        ident = (src or {}).get("identity") or {}
        cells[cid] = dict(l=b["label"], o=b["origin"], p=b["parent"], b=b["born_at_min"], f=b["fate"],
                          u=u if sig else None,
                          alt=[c["label"] for c in ident.get("candidates", [])][1:2],
                          poly=_polygon(by_id[cid]["geometry"], ox, oy) if cid in by_id else None)

    frames = []
    for f in episode["frames"]:
        a, tg, v, cp, cb = [], [], [], [], []
        for s in f["states"]:
            prog = s.get("program") or None
            cyc = s.get("cycle") or None
            act = prog["action"] if prog else ("CYCLING" if cyc else "")
            a.append(ACTIONS.index(act) if act in ACTIONS else 0)
            tg.append(int(prog.get("target_id") or 0) if prog else 0)
            v.append(VIABILITY.index(s["viability"]) if s["viability"] in VIABILITY else 2)
            if cyc and cyc.get("duration"):
                cp.append(round(min(1., max(0., (f["time_min"] - cyc["start"]) / cyc["duration"])), 2))
                cb.append(1 if (cyc.get("blocked") or cyc.get("growth_blocked")) else 0)
            else:
                cp.append(-1)
                cb.append(0)
        frames.append(dict(t=f["time_min"], ids=f["ids"], x=[round(x - ox, 2) for x in f["x"]],
                           y=[round(y - oy, 2) for y in f["y"]], r=[round(r, 2) for r in f["rho"]],
                           a=a, tg=tg, v=v, cp=cp, cb=cb))

    grid = episode["field_grid"]
    field_frames, scales = [], {}
    for name in FIELDS:
        peak = max((max(f["fields"][name]) for f in episode["frames"] if f.get("fields")), default=0.)
        scales[name] = 1. if name == "oxygen" else max(peak, 1e-6)
    for i, f in enumerate(episode["frames"]):
        if f.get("fields"):
            field_frames.append(dict(t=f["time_min"], **{k: _u8(f["fields"][k], scales[k]) for k in FIELDS}))

    events = []
    for e in episode["events"]:
        if e["kind"] in EVENTS:
            events.append([e["time_min"], EVENTS.index(e["kind"]), e.get("cell_id") or 0,
                           e.get("target_id") or e.get("daughter_id") or 0])

    decisions = {}
    for d in episode["decisions"]:
        ans = d.get("answer") or {}
        probs = ans.get("probabilities") or {}
        top = sorted(probs.items(), key=lambda kv: -kv[1])[:4]
        decisions.setdefault(d["cell_id"], []).append(
            [d["time_min"], ans.get("choice") or d["status"], round(float(ans.get("confidence") or 0), 2),
             [[k, round(float(p), 3)] for k, p in top if p >= .005], d.get("wake_reasons", [])[:3], d["status"]])

    summary = episode["summary"]
    return dict(format="virtualtissue-player-1", patch=episode["patch_id"], provider=episode["provider"],
                seed=episode["seed"], scenario=episode["scenario"], hours=summary["time_min"] / 60,
                stop=summary["stop_reason"], paid_calls=episode.get("paid_calls", 0),
                model=next((a.get("model_returned") for a in episode.get("provider_audit", []) if a.get("model_returned")), None),
                origin=[ox, oy], signatures=[s.replace("GOBP_", "") for s in signatures],
                actions=ACTIONS, viability=VIABILITY, events_kinds=EVENTS,
                cells=cells, frames=frames, events=events, decisions=decisions,
                field_grid=dict(origin=[grid["origin_um"][0] - ox, grid["origin_um"][1] - oy], spacing=grid["spacing_um"],
                                shape=grid["shape"], names=FIELDS, scales=scales),
                fields=field_frames, scope=episode["scope"])


def write_player(episode, atlas, destination):
    payload = build_payload(episode, atlas)
    raw = gzip.compress(json.dumps(payload, separators=(",", ":")).encode(), mtime=0)
    data = base64.b64encode(raw).decode("ascii")
    html = TEMPLATE.read_text(encoding="utf8").replace("__PAYLOAD__", data)
    Path(destination).write_text(html, encoding="utf8")
    return Path(destination)
