"""Versioned provider payload; full gate diagnostics stay in the episode.

The original bundled payload helper is retained for input-bundle compatibility.
The paid provider and its preflight use this module instead.
"""
from xenium_runtime.io import canonical

PROMPT_VERSION = "vtx-choice-3"
V2_STATE_SEMANTICS = (
    "AUCell relative expression scores were measured at t0, conditioned on provisional type and detection depth. "
    "They stay frozen throughout the run; daughter scores are inherited, not remeasured. "
    "Unavailable evidence is missing. Scores are not protein activity or validated malignancy. "
    "SIM fields and rates are assumptions.")
STATE_SEMANTICS = (
    "AUCell relative expression scores at t0, conditioned on provisional type and detection depth; "
    "unavailable evidence is missing. Scores are not protein activity or validated malignancy. "
    "SIM fields and rates are assumptions.")


def declared_tissue(value):
    """An explicit declaration, never an inference from a filename or cell label."""
    if value is None:
        return None
    if not isinstance(value, str) or len(value) > 300 or any(ord(c) < 32 for c in value):
        raise ValueError("Tissue must be a single line of at most 300 characters.")
    return value.strip() or None


def payload(rows, model, *, tissue=None):
    tissue = declared_tissue(tissue)
    setting = (f"Independent local cells in an uncalibrated virtual tissue prototype of {tissue}."
               if tissue else "Independent local cells in an uncalibrated human tumour-interface virtual tissue prototype; tissue not declared.")
    return dict(model=model, state=dict(
        setting=setting,
        semantics=STATE_SEMANTICS, prompt_version=PROMPT_VERSION), questions={
        f"c{r['cell_id']}": dict(type="choice", instructions=dict(
            task="Choose one permitted local behavior for this cell.",
            local_observation={k: v for k, v in r["observation"].items()
                               if k not in ("recognition_sim", "excluded")}),
            criteria={k: v["description"] for k, v in r["options"].items()}) for r in rows})


def v2_payload(rows, model):
    """Frozen 0.1.2 comparator for the offline/paid ablation tool only."""
    body = payload(rows, model)
    body["state"] = dict(setting="Independent local cells in an uncalibrated Xenium virtual tissue mechanism simulation.",
                         semantics=V2_STATE_SEMANTICS, prompt_version="vtx-choice-2")
    return body


def build_batches(requests, model, p, *, tissue=None):
    # CD8 menus differ across recognition arms. Keep them out of other cells'
    # request context and use stable ordering within each group.
    batches = []
    for cd8 in (False, True):
        rows = sorted((r for r in requests if
                       (r["observation"].get("primary_label") == "T_CD8_cytotoxic") == cd8),
                      key=lambda r: r["cell_id"])
        batch = []
        for row in rows:
            candidate = batch + [row]
            if len(candidate) > p["max_questions"] or len(canonical(payload(candidate, model, tissue=tissue))) > p["max_payload_bytes"]:
                if not batch:
                    raise ValueError("Single question exceeds the request size limit")
                batches.append(batch)
                batch = [row]
                if len(canonical(payload(batch, model, tissue=tissue))) > p["max_payload_bytes"]:
                    raise ValueError("Single question exceeds the request size limit")
            else:
                batch = candidate
        if batch:
            batches.append(batch)
    return batches
