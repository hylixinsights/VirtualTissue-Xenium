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
            if cd8 < 10:
                notes.append("Low CD8 count: outcomes may depend on a few cells and their contacts. "
                             "The warning threshold of 10 is a usability heuristic, not a validated biological cutoff.")
        notes.append(f"Growth headroom: {facts['growth_headroom']} new bodies before the population ceiling. "
                     "This limit includes context and removed bodies.")
        eligible = facts.get("initially_cycling_or_eligible", 0)
        if settings.get("horizon_hours", 0) >= 24 and eligible > facts["growth_headroom"]:
            notes.append(f"Population ceiling may constrain this horizon: {eligible} cells are initially cycling or eligible "
                         "to enter cycle. This is a capacity warning, not a prediction of divisions.")
        if "initial_http_batches" in facts:
            notes.append(f"First round needs {facts['initial_http_batches']} requests. Total demand depends on events; "
                         "equal request budgets do not ensure equal simulated time. Compare only the common saved interval.")
    return notes
