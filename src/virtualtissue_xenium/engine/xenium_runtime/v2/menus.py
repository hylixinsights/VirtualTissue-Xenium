"""Licensed options and complete three-valued exclusion reasons.

This operational SIM extension does not turn RNA signatures into protein or
recognition. Identity uses the shared conservative/declared-primary evaluator.
Contact always comes from the measured-to-simulated polygon graph.
"""
from ..conditions import aggregate, compare, identity_gate
from .mechanics import TUMOUR

CD8 = "T_CD8_cytotoxic"
MACRO = "Macrophage"


def menu(sim, g):
    t, p, s = sim.tissue, sim.p, sim.bio[g]
    label, cell = t.labels[g], sim.observation_cell(g)
    options, excluded = {}, {}
    neighbours = sim.neighbours.get(g, [])
    vals = sim.local_fields[g]
    def offer(action, description, checks=(), target=None):
        key = action if target is None else f"{action}:{int(t.ids[target])}"
        gates = [("implemented_handler", action in sim.contract["handler_actions"]),
                 ("declared_functional_profile", action in sim.contract["assumed_competence_profiles"].get(sim.scenario.get("assumed_competence_profile"),{}).get("actions",[])),
                 ("identity_competence", identity_gate(cell, {label}, sim.scenario["identity_policy"])),
                 ("actor_alive", s["viability"] == "alive"), *checks]
        reports = [{"condition": name, "result": compare(value, "eq", True)} for name,value in gates]
        if aggregate([v["result"] for v in reports]) == "true":
            options[key] = dict(action=action, target_index=target, target_id=None if target is None else int(t.ids[target]), description=description, gates=reports)
        else:
            excluded[key] = reports

    if label in ("B_cell", "Plasma_cell", "occluded"):
        return options, excluded
    offer("WAIT", "Start no new program.")
    tumors = [k for k in neighbours if t.labels[k] in TUMOUR and sim.bio[k]["viability"] == "alive"]
    if label == CD8:
        offer("REST", "Rest and recover away from recognized antigen.")
        offer("HOLD_POSITION", "Remain in place for the next program interval.")
        move = [("motility_enabled", sim.scenario["motility"])]
        offer("SCAN_RANDOM", "Scan by persistent random motility.", move)
        offer("MIGRATE_UP_CXCL9", "Crawl up the local simulated chemokine gradient.", move + [("chemokine_gradient", sim.gradient_norm[g] >= p["chemotaxis_min_gradient"])])
        if not tumors:
            seen = sim.nearest(g, lambda k: t.acts[k] and t.labels[k] in TUMOUR and sim.bio[k]["viability"] == "alive", p.get("sensing_radius_um", 60))
            if seen is not None:
                offer("APPROACH_TUMOR", "Crawl toward this nearby tumor cell to make contact.", move, seen)
        if tumors:
            offer("HOLD_CONTACT", "Hold this interface without a pulse.", target=tumors[0])
            offer("DISENGAGE", "Crawl away from this tumor interface.", move, tumors[0])
        for k in tumors:
            offer("DEGRANULATE", "Begin a timed cytotoxic contact program; the kernel owns its outcome.",
                  [("target_is_acting", bool(t.acts[k])), ("recognition_declared", sim.scenario["recognition"]),
                   ("effector_competence_declared", sim.scenario["effector_competence"]),
                   ("effector_material", s["pulses"] > 0), ("cooldown_complete", sim.time >= s["cooldown_until"]),
                   ("not_terminally_exhausted", s["exhaustion"] < p["exhaustion_terminal"]),
                   ("glucose_available", vals["glucose"] >= p["glucose_gate"])], k)
        offer("PROLIFERATE", "Start a helped CD8 expansion clock, subject to space at completion.",
              [("not_cycling", s["cycle"] is None), ("helped_divisions_left", s["divisions_left"] > 0), ("recent_stimulation", s["last_stimulation"] is not None),
               ("not_terminal", s["exhaustion"] < p["exhaustion_terminal"]), ("space", sim.lp[g] <= sim.div["cycle_entry_max_local_pressure"])])
    elif label in TUMOUR:
        offer("REMAIN_QUIESCENT", "Remain out of cycle until the next checkpoint.", [("not_cycling", s["cycle"] is None)])
        offer("ENTER_CYCLE", "Start the accelerated SIM cycle; the kernel checks space again at division.",
              [("cycles_enabled", sim.scenario["tumour_cycle"]), ("not_cycling", s["cycle"] is None),
               ("oxygen", vals["oxygen"] >= p["cycle_oxygen_gate"]), ("glucose", vals["glucose"] >= p["glucose_gate"]),
               ("undamaged", s["damage_au"] == 0), ("space", sim.lp[g] <= sim.div["cycle_entry_max_local_pressure"])])
        offer("UPREGULATE_PDL1", "Increase simulated inhibitory ligand display.",
              [("ifng_signal", vals["ifng"] >= p["signal_threshold"]), ("below_maximum", s["pdl1"] < 1)])
        offer("SECRETE_VEGFA", "Release VEGF-like source during a finite program.",
              [("hypoxia_or_declared_nerve", vals["oxygen"] < p["hypoxia_threshold"] or (sim.nerve_effect("ne_tumour_vegf") and vals["ne"] >= p["signal_threshold"]))])
        offer("GLYCOLYTIC_SHIFT", "Increase this cell's glycolytic fraction by the contract increment.",
              [("metabolism_menu", sim.scenario["warburg"]), ("below_maximum", s["glycolysis"] < 1)])
        offer("OXIDATIVE_SHIFT", "Decrease this cell's glycolytic fraction by the contract increment.",
              [("metabolism_menu", sim.scenario["warburg"]), ("above_minimum", s["glycolysis"] > 0)])
    elif label == MACRO:
        if not any(sim.bio[k]["viability"] in ("apoptosis_committed", "dead_intact") for k in neighbours):
            corpse = sim.nearest(g, lambda k: sim.bio[k]["viability"] in ("apoptosis_committed", "dead_intact"), p.get("corpse_signal_radius_um", 80))
            if corpse is not None:
                offer("APPROACH_CORPSE", "Crawl toward this nearby apoptotic cell (find-me signal).",
                      [("motility_enabled", sim.scenario["motility"]), ("processing_capacity", s["processed"] < p["clearance_capacity"])], corpse)
        for k in neighbours:
            if sim.bio[k]["viability"] in ("apoptosis_committed", "dead_intact"):
                offer("EFFEROCYTOSE", "Begin finite-capacity removal of this contacted corpse.",
                      [("target_is_acting", bool(t.acts[k])), ("clearance_competence", sim.scenario["clearance_competence"]),
                       ("processing_capacity", s["processed"] < p["clearance_capacity"])], k)
        dwell = [("polarization_dwell_complete", sim.time >= s["polarized_until"])]
        offer("SECRETE_CXCL9", "Enter an inflammatory chemokine secretion program.", dwell + [("ifng_signal", vals["ifng"] >= p["signal_threshold"])])
        offer("SECRETE_VEGFA_TGFB", "Enter a reversible suppressive secretion program.", dwell + [("local_stress", vals["lactate"] >= p["signal_threshold"] or vals["oxygen"] < p["hypoxia_threshold"])])
        offer("PATROL", "Crawl through the local tissue.", [("motility_enabled", sim.scenario["motility"])])
    elif label == "Endothelial":
        offer("STAY_QUIESCENT", "Remain without a sprout.")
        offer("SPROUT", "Extend a non-perfused endothelial tip toward VEGF; connection is checked later.",
              [("angiogenesis_enabled", sim.scenario["angiogenesis"]), ("vegf_signal", vals["vegf"] >= p["signal_threshold"]),
               ("no_existing_tip", not s["tip"]), ("space", sim.lp[g] <= sim.div["cycle_entry_max_local_pressure"])])
    elif label in ("CAF", "Fibroblast"):
        gate = [("tumor_interface", bool(tumors))]
        offer("DEPOSIT_MATRIX", "Deposit matrix that increases local motility drag.", gate)
        offer("SECRETE_TGFB", "Secrete TGF-beta-like signal for a finite interval.", gate)
        if label == "CAF":
            offer("SECRETE_VEGFA", "Secrete VEGF-like signal at the tumor interface.", gate)
    elif label == "T_CD4":
        for k in neighbours:
            if t.labels[k] == CD8:
                offer("HELP_CD8", "Deliver declared CD4 help at a stimulated CD8 contact.",
                      [("help_profile_declared", sim.scenario["effector_competence"]), ("target_is_acting", bool(t.acts[k])),
                       ("target_alive", sim.bio[k]["viability"] == "alive"),
                       ("target_stimulated", sim.bio[k]["last_stimulation"] is not None)], k)
        offer("MIGRATE_UP_CXCL9", "Crawl along a local simulated chemokine gradient.",
              [("motility_enabled", sim.scenario["motility"]), ("chemokine_gradient", sim.gradient_norm[g] >= p["chemotaxis_min_gradient"])])
    if "STAY_QUIESCENT" in options or "REMAIN_QUIESCENT" in options:
        options.pop("WAIT", None)  # two synonyms for inactivity are not two behaviors
    return options, excluded
