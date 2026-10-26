"""Explain the scientific question and limits of the existing experiment settings."""

SCENARIO_NOTES = {
    "baseline_unknown": (
        "Scenario: baseline_unknown. Immune killing and corpse engulfment are not enabled: "
        "recognition and competence are unknown. To test immune killing, explicitly choose "
        "declared_recognition; increasing the duration does not enable it."),
    "declared_recognition": (
        "Scenario: declared_recognition. Assumes tumour recognition and effector/clearance "
        "competence. Contact-dependent killing is permitted, not guaranteed. Compare with "
        "recognition_off using the same bundle, seed, duration, model and population ceiling."),
    "recognition_off": (
        "Scenario: recognition_off. Negative control: competence is declared, but tumour "
        "recognition is disabled, so CD8 cytotoxic killing is not permitted."),
}


def experiment_notes(settings, facts=None):
    notes = [SCENARIO_NOTES[settings["scenario"]]]
    if facts is not None:
        cd8 = facts["acting_composition"].get("T_CD8_cytotoxic", 0)
        if cd8 == 0:
            notes.append("No active CD8 T cells in this region; the current scenarios cannot test CD8 killing here.")
        else:
            notes.append(f"Active CD8 T cells: {cd8}. Inspect their contacts and access to tumour when interpreting killing.")
        notes.append(f"Growth headroom: {facts['growth_headroom']} new bodies before the population ceiling. "
                     "This limit includes context and removed bodies.")
    return notes
