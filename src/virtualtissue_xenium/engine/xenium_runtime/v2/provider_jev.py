"""Payload size estimation only. No network provider is shipped."""
from ..io import canonical
STATE_SEMANTICS = "AUCell relative expression scores at t0, conditioned on provisional type and detection depth; unavailable evidence is missing. Scores are not protein activity or validated malignancy. SIM fields and rates are assumptions."

def payload(rows,model):
    return dict(model=model,state=dict(setting="Independent local cells in an uncalibrated LUAD virtual tissue prototype.",
                semantics=STATE_SEMANTICS),
                questions={f"c{r['cell_id']}":dict(type="choice",instructions={"task":"Choose one permitted local behavior for this cell.","local_observation":r["observation"]},
                    criteria={k:v["description"] for k,v in r["options"].items()}) for r in rows})

def build_batches(requests,model,p):
    batches, batch = [],[]
    for r in requests:
        candidate = batch+[r]
        if len(candidate) > p["max_questions"] or len(canonical(payload(candidate,model))) > p["max_payload_bytes"]:
            if not batch: raise ValueError("Single question exceeds the request size limit")
            batches.append(batch)
            batch = [r]
            if len(canonical(payload(batch,model))) > p["max_payload_bytes"]:
                raise ValueError("Single question exceeds the request size limit")
        else: batch = candidate
    if batch: batches.append(batch)
    return batches
