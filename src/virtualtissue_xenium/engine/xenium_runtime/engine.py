"""Answer validator retained from Colab_Prototype; no legacy atlas loader."""
import math

def validate_answer(answer, options):
    if not isinstance(answer, dict) or answer.get("choice") not in options:
        raise ValueError("Unknown or missing choice")
    probs = answer.get("probabilities")
    if not isinstance(probs, dict) or set(probs) != set(options):
        raise ValueError("Probability keys must exactly match submitted options")
    if any(type(v) not in (float, int) or not math.isfinite(v) or not 0 <= v <= 1 for v in probs.values()):
        raise ValueError("Invalid probability value")
    if abs(sum(probs.values()) - 1) > 1e-6:
        raise ValueError("Probabilities must sum to one")
    confidence = answer.get("confidence")
    if type(confidence) not in (int, float) or not math.isfinite(confidence) or not 0 <= confidence <= 1:
        raise ValueError("Invalid confidence")
    if probs[answer["choice"]] + 1e-6 < max(probs.values()):
        raise ValueError("Choice is not a maximum-probability option")
