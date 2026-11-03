"""Measure repeat/order sensitivity from preserved journals; never sends requests."""
from pathlib import Path
from collections import defaultdict
from statistics import mean, variance
from .inputs import read_json, digest, write_json
from .provider import validate_response


def audit_journals(paths, output=None):
    directories = [Path(p).resolve() for p in paths]
    if len(directories) < 2 or len(set(directories)) != len(directories):
        raise ValueError("Provide at least two independent journal directories, not repeated copies of one path.")
    observations = defaultdict(list)
    incomplete = []
    usage = {"input_tokens": 0, "output_tokens": 0, "unknown_attempts": 0}
    configs = set()
    for replicate, path in enumerate(directories):
        identity = read_json(path/"identity.json")["identity"]
        seen = set()
        for intent in sorted(path.glob("*.intent.json")):
            record = read_json(intent)
            body = record["request"]
            request_hash = digest(body)
            if record.get("identity") != identity or record.get("request_hash") != request_hash or intent.name != request_hash+".intent.json":
                raise ValueError("Journal intent identity/hash mismatch.")
            result_path = intent.with_name(request_hash+".result.json")
            if not result_path.exists():
                incomplete.append({"replicate": replicate, "request": request_hash, "status": "unresolved"})
                usage["unknown_attempts"] += 1
                continue
            result = read_json(result_path)
            checksum = result.pop("digest", None)
            if checksum != digest(result) or result.get("identity") != identity or result.get("request_hash") != request_hash:
                raise ValueError("Journal result identity/hash mismatch.")
            if result.get("status") != "valid":
                incomplete.append({"replicate": replicate, "request": request_hash, "status": "failed"})
                usage["unknown_attempts"] += 1
                continue
            answers = validate_response(result["response"], body)
            counts = result["response"].get("usage", {})
            if all(type(counts.get(k)) is int and counts[k] >= 0 for k in ("input_tokens", "output_tokens")):
                for k in ("input_tokens", "output_tokens"):
                    usage[k] += counts[k]
            else:
                usage["unknown_attempts"] += 1
            configs.add((body["model"], body["state"]["prompt_version"]))
            # Keep biological state and co-batched questions fixed; only ordering may differ.
            members = sorted(digest(q) for q in body["questions"].values())
            for order, (qid, question) in enumerate(body["questions"].items()):
                key = digest([body["model"], body["state"], members, question])
                if key in seen:
                    raise ValueError("A question appears repeatedly within one replicate; split independent repeats into journals.")
                seen.add(key)
                observations[key].append({"replicate": replicate, "order": order,
                                          "probabilities": answers[qid]["probabilities"]})
    if len(configs) != 1:
        raise ValueError("Replicates must share one model and prompt version.")
    questions = []
    for key, values in sorted(observations.items()):
        probabilities = {a: [v["probabilities"][a] for v in values] for a in values[0]["probabilities"]}
        questions.append({"question_sha256": key, "n_replicates": len(values),
                          "all_replicates_present": len(values) == len(directories),
                          "question_order_varied": len({v["order"] for v in values}) > 1,
                          "mean_probabilities": {a: mean(p) for a,p in probabilities.items()},
                          "probability_sample_variance": {a: variance(p) if len(p)>1 else None for a,p in probabilities.items()}})
    result = {"schema_version": "3.0", "model_and_prompt": list(next(iter(configs))),
              "n_requested_replicates": len(directories), "n_matched_questions": sum(q["all_replicates_present"] for q in questions),
              "journals": list(map(str, directories)), "questions": questions, "incomplete_requests": incomplete,
              "usage": usage, "sampling_parameters": "Not exposed by the adapter; inspect provider contract before claiming determinism.",
              "interpretation": "Empirical repeated-payload variance; zero observed variance does not establish determinism. Order effects require question_order_varied=true."}
    if output is not None:
        if Path(output).exists():
            raise ValueError("Audit output exists; choose a new path.")
        write_json(output, result)
    return result
