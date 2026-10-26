"""Versioned provider payload; full gate diagnostics stay in the episode.

The original bundled payload helper is retained for input-bundle compatibility.
The paid provider and its preflight use this module instead.
"""
from xenium_runtime.io import canonical

PROMPT_VERSION = "vtx-choice-2"
STATE_SEMANTICS = (
    "AUCell relative expression scores were measured at t0, conditioned on provisional type and detection depth. "
    "They stay frozen throughout the run; daughter scores are inherited, not remeasured. "
    "Unavailable evidence is missing. Scores are not protein activity or validated malignancy. "
    "SIM fields and rates are assumptions.")


def payload(rows, model):
    return dict(model=model, state=dict(
        setting="Independent local cells in an uncalibrated Xenium virtual tissue mechanism simulation.",
        semantics=STATE_SEMANTICS, prompt_version=PROMPT_VERSION), questions={
        f"c{r['cell_id']}": dict(type="choice", instructions=dict(
            task="Choose one permitted local behavior for this cell.",
            local_observation={k: v for k, v in r["observation"].items()
                               if k not in ("recognition_sim", "excluded")}),
            criteria={k: v["description"] for k, v in r["options"].items()}) for r in rows})


def build_batches(requests, model, p):
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
            if len(candidate) > p["max_questions"] or len(canonical(payload(candidate, model))) > p["max_payload_bytes"]:
                if not batch:
                    raise ValueError("Single question exceeds the request size limit")
                batches.append(batch)
                batch = [row]
                if len(canonical(payload(batch, model))) > p["max_payload_bytes"]:
                    raise ValueError("Single question exceeds the request size limit")
            else:
                batch = candidate
        if batch:
            batches.append(batch)
    return batches
