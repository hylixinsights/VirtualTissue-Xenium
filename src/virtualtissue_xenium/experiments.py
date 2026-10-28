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
        activity = facts.get("t0_activity")
        if activity and cd8:
            c = activity["cd8"]
            notes.append(f"CD8 access at t0: {c['touching_tumour']} of {c['acting']} touch a live tumour cell; "
                         f"{c['tumour_within_sensing_radius']} have one within {c['sensing_radius_um']:g} µm. "
                         f"{c['without_question_at_t0']} have no eligible question at t0; "
                         f"{c['deferred_by_decision_cap']} eligible questions are deferred by the decision cap. Cells decide only when an event "
                         "wakes them (new tumour contact, field-threshold crossing, end of a program); a CD8 cell "
                         "that is never woken cannot choose an autonomous movement program.")
        if activity:
            missing = {k: v["without_aucell"] for k, v in activity["by_type"].items() if v["without_aucell"]}
            if missing:
                total = sum(v["acting"] for v in activity["by_type"].values())
                notes.append(f"No t0 AUCell for {sum(missing.values())} of {total} acting cells "
                             f"({', '.join(f'{k} {v}/{activity['by_type'][k]['acting']}' for k, v in sorted(missing.items()))}). "
                             "Check detected-gene depth (--aucell-top-ranks), QC and gene-set overlap in the preparation audit. Missing scores gate no "
                             "action or kernel rate; Jev sees them as missing, and the local preview rules never "
                             "start a tumour cycle without a MITOTIC_CELL_CYCLE score.")
        if facts.get("prompt_context"):
            context = facts["prompt_context"]
            notes.append(f"Jev prompt {context['version']}: context mode = {context['mode']}; "
                         f"tissue = {context['tissue'] or 'not declared'} "
                         f"(source: {context['tissue_source']}). Tissue context can introduce literature priors; "
                         "the effect on ENTER_CYCLE has not yet been validated with live Jev.")
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
