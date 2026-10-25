import copy
import tempfile
import unittest
import urllib.error
from pathlib import Path
from unittest.mock import patch

from jev_provider import JevProvider, ProviderError, validate_response
from xenium_runtime.providers import BudgetExhausted
from xenium_runtime.v2.provider_jev import payload, build_batches


PARAMS = dict(max_questions=2, max_payload_bytes=40000, max_response_bytes=2000000)
MODEL = "jev-1.14.0"


def rows(n=3):
    return [dict(cell_id=i, observation=dict(time_min=0, primary_label="Tumor", t0_aucell_raw={}),
                 options={"WAIT": {"description": "Wait"}, "ENTER_CYCLE": {"description": "Enter cycle"}}) for i in range(n)]


def response(body):
    return dict(model=body["model"], usage=dict(input_tokens=100, output_tokens=20), answers={
        q: dict(type="choice", choice=next(iter(v["criteria"])),
                probabilities={k: float(i == 0) for i, k in enumerate(v["criteria"])}, confidence=.8)
        for q, v in body["questions"].items()})


class ProviderTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / "journal.jsonl"
        self.calls = []

    def tearDown(self):
        self.temp.cleanup()

    def send(self, body):
        self.calls.append(body)
        return response(body)

    def provider(self, **kwargs):
        opts = dict(key="FAKE_SECRET_NEVER_PERSIST", model=MODEL, parameters=PARAMS,
                    limit=2, journal=self.path, context={"seed": 21}, transport=self.send, pause_s=0)
        opts.update(kwargs)
        return JevProvider(**opts)

    def test_budget_is_atomic_before_round(self):
        p = self.provider(limit=1)
        with self.assertRaises(BudgetExhausted):
            p.choose(rows())
        self.assertEqual(len(self.calls), 0)
        self.assertEqual(p.requests, 0)

    def test_resume_reuses_calls_and_counts_usage_once(self):
        p = self.provider()
        expected = p.choose(rows())
        again = self.provider()
        self.assertEqual(again.choose(rows()), expected)
        self.assertEqual(len(self.calls), 2)
        self.assertEqual(again.requests, 2)
        self.assertEqual(again.usage()["input_tokens"], 200)
        self.assertNotIn("FAKE_SECRET", self.path.read_text())

    def test_changed_context_and_corruption_rejected(self):
        self.provider().choose(rows(1))
        with self.assertRaises(ProviderError):
            self.provider(context={"seed": 22})
        self.path.write_text(self.path.read_text() + '{"partial":')
        with self.assertRaises(ProviderError):
            self.provider()

    def test_interrupt_after_intent_is_never_resent(self):
        def interrupted(body):
            raise KeyboardInterrupt()
        p = self.provider(transport=interrupted)
        with self.assertRaises(KeyboardInterrupt):
            p.choose(rows(1))
        again = self.provider()
        with self.assertRaises(ProviderError):
            again.choose(rows(1))
        self.assertEqual(again.requests, 1)
        self.assertEqual(len(self.calls), 0)

    def test_http_429_is_one_attempt_without_retry(self):
        class Opener:
            calls = 0
            def open(self, *a, **kw):
                self.calls += 1
                raise urllib.error.HTTPError("https://api.typesafe.ai", 429, "limited", {}, None)
        opener = Opener()
        with patch("jev_provider.urllib.request.build_opener", return_value=opener):
            p = self.provider(transport=None)
            with self.assertRaises(ProviderError):
                p.choose(rows(1))
        self.assertEqual(opener.calls, 1)
        self.assertEqual(p.requests, 1)
        with self.assertRaises(ProviderError):
            self.provider().choose(rows(1))
        self.assertEqual(len(self.calls), 0)

    def test_response_validation(self):
        body = payload(rows(1), MODEL)
        mutations = [lambda r: r.update(model="jev-latest"),
                     lambda r: r["answers"].update(extra={}),
                     lambda r: r["answers"]["c0"].update(choice="ENTER_CYCLE"),
                     lambda r: r["answers"]["c0"].update(confidence=True),
                     lambda r: r["answers"]["c0"]["probabilities"].update(WAIT=float("nan")),
                     lambda r: r["answers"]["c0"]["probabilities"].update(WAIT=.2)]
        for mutation in mutations:
            raw = response(body)
            mutation(raw)
            with self.assertRaises(ProviderError):
                validate_response(raw, body)
        raw = response(body)
        raw["answers"]["c0"]["probabilities"] = dict(WAIT=.505, ENTER_CYCLE=.49)
        normalized, repairs = validate_response(raw, body)
        self.assertAlmostEqual(sum(normalized["c0"]["probabilities"].values()), 1)
        self.assertTrue(repairs["c0"])

    def test_arbitrary_errors_and_response_extras_do_not_leak_secrets(self):
        def bad(body):
            raise RuntimeError("Bearer FAKE_SECRET_NEVER_PERSIST")
        with self.assertRaises(ProviderError):
            self.provider(transport=bad).choose(rows(1))
        self.assertNotIn("FAKE_SECRET", self.path.read_text())

    def test_payload_size_and_question_bounds(self):
        from xenium_runtime.io import canonical
        batches = build_batches(rows(7), MODEL, PARAMS)
        self.assertEqual([len(b) for b in batches], [2, 2, 2, 1])
        for b in batches:
            self.assertLessEqual(len(canonical(payload(b, MODEL))), PARAMS["max_payload_bytes"])
        with self.assertRaises(ValueError):
            build_batches(rows(1), MODEL, dict(PARAMS, max_payload_bytes=1))


if __name__ == "__main__":
    unittest.main()
