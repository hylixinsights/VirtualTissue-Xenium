"""Actual R03 integration: mocked Jev, replay, caps, failures, and input integrity."""
import argparse
import copy
import json
import tempfile
import time
import zipfile
from pathlib import Path

import simulation_workflow as w
from jev_provider import ProviderError
from test_provider import response
from xenium_runtime.v2.policy import ScoreAwareLocalPolicy
from xenium_runtime.v2.episode import validate_episode


def mock_score_response(body):
    rows = [dict(cell_id=int(q[1:]), observation=v["instructions"]["local_observation"],
                 options={k: dict(description=d) for k,d in v["criteria"].items()}) for q,v in body["questions"].items()]
    answers = ScoreAwareLocalPolicy().choose(rows)
    return dict(model=body["model"], answers={f"c{k}":dict(v,type="choice") for k,v in answers.items()},
                usage=dict(input_tokens=100,output_tokens=20))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--zip", default="Xenium_Lung5K_Prepared/simulation_input_R03_180834_343907.zip")
    args = parser.parse_args()
    out = Path(args.output); out.mkdir(parents=True,exist_ok=True)
    started = time.monotonic()
    bundle = w.extract_bundle(args.zip, out / "extracted")
    settings = dict(model="jev-1.14.0", scenario="declared_recognition", seed=21,
                    horizon_hours=1, max_decisions=10000, max_population=1000)
    facts, preview = w.preflight(bundle,settings)
    calls = []
    def transport(body):
        calls.append(body)
        return mock_score_response(body)
    initial = facts["initial_http_batches"]
    # First round does not fit: no calls and a valid t0 replay.
    capped = w.run(bundle,settings,out / "cap_before_round",max_http_attempts=initial-1,transport=transport)
    assert not calls and capped["report"]["stop_reason"] == "request_cap_before_round"
    assert capped["report"]["checkpoint_time_min"] == 0
    # A real trajectory with local-score decisions formatted as Jev Choice responses.
    result = w.run(bundle,settings,out / "mock_jev",max_http_attempts=1000,transport=transport)
    assert result["report"]["completed"], result["report"]
    count = len(calls)
    assert count > 0 and result["report"]["usage"]["http_attempts"] == count
    assert result["episode"]["summary"]["invalid_answers"] == 0
    # Identical configuration resumes from t0, no new transport calls.
    resumed = w.run(bundle,settings,out / "mock_jev",max_http_attempts=1200,transport=transport)
    assert len(calls) == count
    assert result["episode"]["final_state_fingerprint"] == resumed["episode"]["final_state_fingerprint"]
    assert resumed["report"]["usage"]["reused_batches"] == count
    # Even a failure after a successful first batch applies no part of the round.
    failures = []
    def partial_failure(body):
        failures.append(body)
        if len(failures) == 2:
            raise ProviderError("HTTP 503. No automatic retry was made.")
        return mock_score_response(body)
    failed = w.run(bundle,settings,out / "partial_provider_failure",max_http_attempts=1000,transport=partial_failure)
    assert failed["report"]["stop_reason"] == "provider_failed"
    assert failed["report"]["usage"]["http_attempts"] == 2
    assert not failed["episode"]["decisions"] and failed["report"]["checkpoint_time_min"] == 0
    # A mid-tick interrupt preserves the last complete checkpoint and its valid final frame.
    import unittest.mock
    original = w.AutonomousSimulation.substep
    hit = [0]
    def interrupted(self,dt):
        hit[0] += 1
        if hit[0] == 2:
            raise KeyboardInterrupt()
        return original(self,dt)
    with unittest.mock.patch.object(w.AutonomousSimulation,"substep",interrupted):
        interrupted_run = w.run(bundle,settings,out / "interrupted",max_http_attempts=1000,transport=transport)
    assert interrupted_run["report"]["error"] and interrupted_run["report"]["checkpoint_time_min"] == 0
    validate_episode(interrupted_run["episode"],bundle["atlas"])
    resumed_interrupted = w.run(bundle,settings,out / "interrupted",max_http_attempts=1000,transport=transport)
    assert resumed_interrupted["episode"]["final_state_fingerprint"] == result["episode"]["final_state_fingerprint"]
    # Scenario recognition is never inferred from AUCell, and input contract is unchanged.
    assert bundle["contract"]["autonomous_scenarios"]["S0"]["recognition"] is None
    for scenario, expected in [("baseline_unknown",None),("declared_recognition",True),("recognition_off",False)]:
        c = w.effective_contract(bundle,dict(settings,scenario=scenario))
        assert c["autonomous_scenarios"][scenario]["recognition"] is expected
    # Reject traversal and symlink entries before imports; checksum corruption also fails.
    for name in ("../escape.py", "/absolute.py", "bad\\escape.py"):
        bad = out / "unsafe.zip"
        with zipfile.ZipFile(bad,"w") as z: z.writestr(name,"pass")
        try: w.extract_bundle(bad,out / "rejected")
        except ValueError: pass
        else: raise AssertionError("Unsafe ZIP was accepted")
        bad.unlink()
    bad = out / "modified.zip"
    with zipfile.ZipFile(args.zip) as src, zipfile.ZipFile(bad,"w",zipfile.ZIP_DEFLATED) as dst:
        for info in src.infolist():
            content=src.read(info)
            if info.filename.endswith("runtime_contract.json"): content += b" "
            dst.writestr(info,content)
    try: w.extract_bundle(bad,out / "rejected")
    except ValueError: pass
    else: raise AssertionError("Modified bundle was accepted")
    bad.unlink()
    report = dict(passed=True, versions=w.versions(), zip_sha256=bundle["zip_sha256"],
                  initial_cells=facts["total_cells"],initial_questions=facts["initial_cell_questions"],
                  initial_batches=initial, mock_http_batches=count, actual_network_calls=0,
                  simulated_minutes=result["report"]["checkpoint_time_min"],
                  decisions=result["episode"]["summary"]["decisions"],events=result["episode"]["summary"]["events"],
                  deterministic_resume=True,cap_before_round=True,partial_round_failure_atomic=True,
                  interrupted_checkpoint_preserved=True,unsafe_and_modified_input_rejected=True,
                  elapsed_seconds=round(time.monotonic()-started,2),
                  final_state_fingerprint=result["episode"]["final_state_fingerprint"])
    w.write_json(out / "validation_report.json",report)
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
