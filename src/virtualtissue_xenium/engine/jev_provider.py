"""Bounded Jev Choice calls with a durable, configuration-bound response journal.

No retries or redirects. A pending intent is an unknown outcome, never a free retry.
The simulation replays from t0 using successful cached answers after a disconnect.
"""
import json
import math
import os
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

from xenium_runtime.io import canonical, fingerprint
from xenium_runtime.providers import BudgetExhausted
from xenium_runtime.v2.provider_jev import build_batches, payload

ENDPOINT = "https://api.typesafe.ai/v1/systemone"


class ProviderError(RuntimeError):
    pass


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def validate_response(response, body):
    if not isinstance(response, dict) or response.get("model") != body["model"]:
        raise ProviderError("Returned model differs from the requested pinned version.")
    answers = response.get("answers")
    if not isinstance(answers, dict) or set(answers) != set(body["questions"]):
        raise ProviderError("Response question IDs differ from the request.")
    result, repairs = {}, {}
    for qid, question in body["questions"].items():
        a = answers[qid]
        if not isinstance(a, dict) or a.get("type") != "choice":
            raise ProviderError("Expected a Choice answer.")
        p = a.get("probabilities")
        if not isinstance(p, dict) or set(p) != set(question["criteria"]):
            raise ProviderError("Response options differ from the permitted menu.")
        values = [*p.values(), a.get("confidence")]
        if any(type(v) not in (int, float) or not math.isfinite(v) or not 0 <= v <= 1 for v in values):
            raise ProviderError("Invalid probability or confidence.")
        total = sum(p.values())
        if total <= 0 or abs(total - 1) > min(.02, .005 * len(p)) + 1e-9:
            raise ProviderError("Probabilities do not sum to one within rounding tolerance.")
        choice = a.get("choice")
        if choice not in p or p[choice] < max(p.values()) - 1e-6:
            raise ProviderError("Choice is absent or is not a maximum-probability option.")
        repairs[qid] = [] if abs(total - 1) < 1e-12 else ["normalized_rounded_probability_sum"]
        result[qid] = dict(choice=choice, probabilities={k: v / total for k, v in p.items()},
                           confidence=a["confidence"])
    return result, repairs


class JevProvider:
    name = "jev_choice_journal_v1"

    def __init__(self, key, model, parameters, limit, journal, context, *,
                 timeout_s=120, pause_s=.25, transport=None, progress=None):
        if not re.fullmatch(r"jev-\d+\.\d+\.\d+", model):
            raise ValueError("Use a pinned Jev version, for example jev-1.14.0.")
        if type(limit) is not int or limit < 1:
            raise ValueError("HTTP attempt limit must be a positive integer.")
        self.key, self.model, self.p, self.limit = key, model, parameters, limit
        if transport is not None:
            self.name = "mock_jev_choice_software_test"
        self.timeout_s, self.pause_s = timeout_s, pause_s
        self.journal = Path(journal)
        self.context = fingerprint(context)
        self.transport, self.progress = transport or self._remote, progress
        self.requests, self.cache_hits, self.audit = 0, 0, []
        self.cache, self.failed, self.last_error = {}, False, None
        self.chain, self.sequence = "0" * 64, 0
        self.journal.parent.mkdir(parents=True, exist_ok=True)
        if self.journal.exists():
            self._read()

    def _read(self):
        for line in self.journal.read_text(encoding="utf8").splitlines():
            try:
                entry = json.loads(line)
                digest = entry.pop("digest")
                if (entry["context"] != self.context or entry["previous"] != self.chain
                        or entry["sequence"] != self.sequence or fingerprint(entry) != digest):
                    raise ValueError()
                h = entry["request_hash"]
                if entry["phase"] == "intent":
                    if h in self.cache or fingerprint(entry["request"]) != h:
                        raise ValueError()
                    self.requests += 1
                    self.cache[h] = entry
                elif entry["phase"] == "result" and self.cache[h]["phase"] == "intent":
                    self.cache[h] = entry
                else:
                    raise ValueError()
                self.chain, self.sequence = digest, self.sequence + 1
            except (ValueError, KeyError, TypeError):
                raise ProviderError("Journal is incomplete, changed, or belongs to a different configuration. Preserve it and inspect it before starting a separate run.") from None
        if self.requests > self.limit:
            raise ProviderError("The request limit is below the number already recorded for this run.")

    def _append(self, entry):
        entry = dict(entry, context=self.context, sequence=self.sequence, previous=self.chain)
        digest = fingerprint(entry)
        with self.journal.open("a", encoding="utf8", newline="\n") as f:
            f.write(canonical(dict(entry, digest=digest)).decode() + "\n")
            f.flush()
            os.fsync(f.fileno())
        self.chain, self.sequence = digest, self.sequence + 1
        return entry

    def _remote(self, body):
        req = urllib.request.Request(ENDPOINT, data=canonical(body), method="POST",
                headers={"Authorization": "Bearer " + self.key, "Content-Type": "application/json"})
        try:
            with urllib.request.build_opener(NoRedirect()).open(req, timeout=self.timeout_s) as response:
                raw = response.read(self.p["max_response_bytes"] + 1)
            if len(raw) > self.p["max_response_bytes"]:
                raise ProviderError("API response exceeded the size limit.")
            return json.loads(raw)
        except urllib.error.HTTPError as e:
            # Never serialize headers, credentials, or an arbitrary server error body.
            raise ProviderError(f"HTTP {e.code}. No automatic retry was made.") from None
        except ProviderError:
            raise
        except Exception:
            raise ProviderError("Transport or response decoding failed; usage is unknown. No automatic retry was made.") from None

    def choose(self, requests):
        batches = build_batches(requests, self.model, self.p)
        bodies = [payload(batch, self.model) for batch in batches]
        hashes = [fingerprint(body) for body in bodies]
        for h in hashes:
            if h in self.cache and self.cache[h].get("status") != "valid":
                self.failed = True
                self.last_error = "This journal contains a failed or unresolved attempt for a required batch. It will not be sent again automatically."
                raise ProviderError(self.last_error)
        fresh = {h for h in hashes if h not in self.cache}
        if self.requests + len(fresh) > self.limit:
            raise BudgetExhausted("The remaining HTTP allowance cannot cover this complete decision round.")
        out = {}
        for batch, body, h in zip(batches, bodies, hashes):
            cached = h in self.cache
            if not cached:
                intent = self._append(dict(phase="intent", request_hash=h, request=body, usage="unknown"))
                self.cache[h] = intent
                self.requests += 1
                try:
                    response = self.transport(body)
                    # Whitelist persisted response fields; never persist arbitrary response extras.
                    safe = {k: response.get(k) for k in ("model", "answers", "usage")} if isinstance(response, dict) else {}
                    if self.key:
                        safe = json.loads(canonical(safe).decode().replace(self.key, "[redacted]"))
                    answers, repairs = validate_response(safe, body)
                    result = dict(phase="result", request_hash=h, status="valid", response=safe,
                                  answers=answers, repairs=repairs, usage=safe.get("usage", "unknown"))
                except Exception as e:
                    self.failed = True
                    # Only our own messages are safe to record.
                    self.last_error = str(e) if isinstance(e, ProviderError) else "Transport or response validation failed."
                    if self.key:
                        self.last_error = self.last_error.replace(self.key, "[redacted]")
                    self.cache[h] = self._append(dict(phase="result", request_hash=h, status="failed",
                                                     error=self.last_error, usage="unknown"))
                    self.audit.append(dict(request_hash=h, cached=False, status="failed", error=self.last_error))
                    raise ProviderError(self.last_error) from None
                self.cache[h] = self._append(result)
                if self.pause_s:
                    time.sleep(self.pause_s)
            else:
                self.cache_hits += 1
            entry = self.cache[h]
            self.audit.append(dict(request_hash=h, cached=cached, status="valid", cells=[r["cell_id"] for r in batch],
                                   model_returned=entry["response"]["model"], usage=entry["usage"], repairs=entry["repairs"]))
            for row in batch:
                out[row["cell_id"]] = entry["answers"][f"c{row['cell_id']}"]
            if self.progress:
                self.progress(self.requests, self.limit, self.cache_hits)
        return out

    def usage(self):
        totals = {"input_tokens": 0, "output_tokens": 0, "attempts_with_unknown_usage": 0}
        for entry in self.cache.values():
            usage = entry.get("usage")
            if not isinstance(usage, dict) or any(type(usage.get(k)) is not int or usage[k] < 0 for k in ("input_tokens", "output_tokens")):
                totals["attempts_with_unknown_usage"] += 1
            else:
                for k in ("input_tokens", "output_tokens"):
                    totals[k] += usage[k]
        return dict(totals, http_attempts=self.requests, reused_batches=self.cache_hits)
