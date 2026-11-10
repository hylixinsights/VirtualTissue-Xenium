"""Sample context is elicited explicitly; only declared restrictions affect menus."""
from copy import deepcopy
import json

import pytest

from virtualtissue_xenium.v2.cli import main
from virtualtissue_xenium.v2.inputs import prepare
from virtualtissue_xenium.v5.context_questions import collect_context
from v5_helpers import PANEL, cells, context, manifest, policy, rule


def answers(values, prompts):
    iterator = iter(values)

    def answer(prompt):
        prompts.append(prompt)
        return next(iterator)

    return answer


def test_questionnaire_collects_english_sample_context_and_keeps_settings():
    template = context(sample_context=dict(anatomical_region="left lower lobe"))
    original = deepcopy(template)
    prompts, messages = [], []
    result = collect_context(template, input_fn=answers(
        ["", "brain", "glioblastoma", "Tumour core and neighbouring tissue.", "", "adult", "Tumour cells and reactive glia.",
         "No mucus-producing cell population was identified.", ""], prompts), output_fn=messages.append)
    assert template == original
    assert result["organism"] == "human"
    assert result["tissue"] == "brain"
    assert result["disease"] == "glioblastoma"
    assert result["sample_context"] == dict(anatomical_region="left lower lobe", developmental_stage="adult",
        cell_state_notes="Tumour cells and reactive glia.", secretion_notes="No mucus-producing cell population was identified.")
    assert result["time_step_hours"] == template["time_step_hours"]
    assert result["recognition"] == template["recognition"]
    assert result["context_modifiers"] == []
    assert any("Which tissue or organ" in prompt for prompt in prompts)
    assert any("disease or pathology" in prompt for prompt in prompts)
    assert any("developmental stage" in prompt for prompt in prompts)
    assert any("secretory cell populations" in prompt for prompt in prompts)
    assert any("Mucus secretion is not a modelled action" in message for message in messages)
    assert any("brain glia, progenitors and tumour cells can differ" in message for message in messages)


def test_blank_answers_preserve_optional_metadata_and_existing_modifiers():
    template = context(sample_context=dict(developmental_stage="adult"), context_modifiers=[dict(
        action="SECRETE", lineages=["malignant"], multiplier=.5, reason="Declared prior assumption.")])
    result = collect_context(template, input_fn=answers([""] * 9, []), output_fn=lambda _: None)
    assert result["sample_context"] == template["sample_context"]
    assert result["context_modifiers"][0]["multiplier"] == .5
    assert len(result["context_modifiers"]) == 1


def test_explicit_restriction_removes_action_despite_positive_cycling_evidence():
    template = context(evidence_policy=policy(rule("PROLIFERATE", "S_CYCLE")))
    updated = collect_context(template, input_fn=answers(
        [""] * 8 + ["malignant:PROLIFERATE", "Sample-specific division restriction for a sensitivity analysis."], []),
        output_fn=lambda _: None)
    baseline = prepare(cells(), template, assayed_genes=PANEL, measured_evidence=manifest())
    restricted = prepare(cells(), updated, assayed_genes=PANEL, measured_evidence=manifest())
    cell_id = "tumour11"
    assert baseline["guidebook"]["cells"][cell_id]["evidence"]["S_CYCLE"]["z"] > 0
    assert "PROLIFERATE" in baseline["guidebook"]["cells"][cell_id]["actions"]
    assert "PROLIFERATE" not in restricted["guidebook"]["cells"][cell_id]["actions"]
    adjustments = restricted["guidebook"]["cells"][cell_id]["weight_adjustments"]
    assert any(row["source"] == "declared_context_modifier" and row["applied"] and row["multiplier"] == 0
               for row in adjustments)


@pytest.mark.parametrize("restriction, message", [
    ("malignant:REST", "REST fallback"),
    ("neuron:PROLIFERATE", "not a permitted action"),
    ("malignant:SECRETE_MUCUS", "not a permitted action"),
    ("goblet:SECRETE", "Unsupported restriction lineage"),
    ("malignant:PROLIFERATE, malignant:PROLIFERATE", "Duplicate restriction"),
    ("malignant", "LINEAGE:ACTION"),
])
def test_invalid_restrictions_fail_before_asking_reasons(restriction, message):
    with pytest.raises(ValueError, match=message):
        collect_context(context(), input_fn=answers([""] * 8 + [restriction], []), output_fn=lambda _: None)


def test_restriction_requires_nonempty_reason():
    with pytest.raises(ValueError, match="answer is required"):
        collect_context(context(), input_fn=answers([""] * 8 + ["malignant:PROLIFERATE", ""], []),
                        output_fn=lambda _: None)


def test_context_init_cli_writes_new_context_and_refuses_overwrite(tmp_path, monkeypatch, capsys):
    template = tmp_path / "template.json"
    template.write_text(json.dumps(context()))
    output = tmp_path / "sample.json"
    monkeypatch.setattr("builtins.input", answers([""] * 4 + ["alveolar region", "adult", "", "", ""], []))
    assert main(["context-init", str(output), "--template", str(template)]) == 0
    assert json.loads(output.read_text())["sample_context"]["anatomical_region"] == "alveolar region"
    assert '"requests_sent": 0' in capsys.readouterr().out
    original = output.read_bytes()
    monkeypatch.setattr("builtins.input", lambda _: pytest.fail("Existing output must fail before prompting."))
    assert main(["context-init", str(output), "--template", str(template)]) == 1
    assert output.read_bytes() == original
    assert "Context output exists" in capsys.readouterr().err


def test_context_init_reports_closed_stdin_without_writing(tmp_path, monkeypatch, capsys):
    template = tmp_path / "template.json"
    template.write_text(json.dumps(context()))
    output = tmp_path / "sample.json"

    def closed_input(_):
        raise EOFError

    monkeypatch.setattr("builtins.input", closed_input)
    assert main(["context-init", str(output), "--template", str(template)]) == 1
    assert not output.exists()
    assert "needs interactive input" in capsys.readouterr().err
