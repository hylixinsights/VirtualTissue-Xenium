"""Declared local policy for milestone M1: a mechanics exercise, not biology.

The kernel builds each menu from its own gates; the policy only ranks the
options it is given, exactly like the future Jev adapter. Answers are one-hot
bookkeeping probabilities, never rates (H21).
"""

PREFERENCE = ("scan", "enter_cycle", "remain_quiescent", "hold")


class M1ExercisePolicy:
    name = "m1_exercise_policy"

    def choose(self, requests):
        answers = {}
        for request in requests:
            options = list(request["options"])
            choice = next((p for p in PREFERENCE if p in options), options[0])
            answers[request["cell_id"]] = dict(choice=choice, probabilities={k: float(k == choice) for k in options}, confidence=1.)
        return answers
