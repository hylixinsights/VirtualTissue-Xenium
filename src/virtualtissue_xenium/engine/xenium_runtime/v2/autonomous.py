"""Autonomous SIM programs on the newly prepared Xenium initial condition.

Tick order: local fields/contact snapshot -> wake/menu/answers -> programs
start simultaneously with reservations -> mechanics -> obligatory clocks ->
new fields -> accounting/readouts -> recorded frame. Probabilities never
enter physics. No outcome target is exposed to a provider.
"""
import copy
import math
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree
from shapely import contains, distance, points, union_all

from ..engine import validate_answer
from ..io import fingerprint
from ..providers import BudgetExhausted
from .fields import FieldGrid
from .episode import engine_fingerprint
from .menus import CD8, MACRO, menu
from .mechanics import TUMOUR
from .policy import LocalPolicyV2
from .runtime import SimulationV2
from .scheduler import WakeScheduler


class AutonomousSimulation(SimulationV2):
    def __init__(self, atlas, contract, scenario="S1", seed=21, policy=None,
                 max_decisions=None, stop_on_provider_failure=False):
        if scenario not in contract["autonomous_scenarios"]:
            raise ValueError("Unknown autonomous scenario")
        if max_decisions is not None and (type(max_decisions) is not int or max_decisions < 1):
            raise ValueError("max_decisions must be a positive integer or None")
        self.max_decisions = max_decisions
        self.stop_on_provider_failure = stop_on_provider_failure
        self.decision_rounds = []
        self.pending_cap_stop = None
        super().__init__(atlas, contract, "static", seed)
        self.engine_fingerprint = engine_fingerprint()
        self.scenario_id = scenario
        self.scenario = self.contract["autonomous_scenarios"][scenario]
        self.p = self.contract["biology"]
        self.policy = policy or LocalPolicyV2()
        self.scheduler = WakeScheduler(self.contract["scheduler"])
        self.bio = [self.initial_state(label) for label in self.tissue.labels]
        self.by_id = {c["cell_id"]: c for c in atlas["cells"]}
        self.failures = []
        self.batches_estimated = []
        self.segments = []
        self.nerve = self.nerve_geometry()
        t = self.tissue
        # A nerve perturbation must not change the PDE domain of its control.
        # The external line deposits its flux onto the nearest boundary nodes.
        extent = np.vstack([t.xy-t.rho[:,None], t.xy+t.rho[:,None]])
        bounds = [*extent.min(0)-4, *extent.max(0)+4]
        self.fields = FieldGrid(bounds, self.contract["fields"]["spacing_um"], self.contract["fields"]["specs"])
        vessels = [g for g,l in enumerate(t.labels) if l == "Endothelial"]
        count = int(round(len(vessels)*self.scenario["nonperfused_fraction"]))
        off = set(self.rng.choice(vessels, count, replace=False).tolist()) if count else set()
        for g in vessels:
            self.bio[g]["perfused"] = g not in off
        for g in np.flatnonzero(t.acts):
            if t.labels[g] in TUMOUR and self.scenario["tumour_cycle"]:
                self.scheduler.wake(int(g), "initial_cycle_checkpoint")
            if t.labels[g] == "Tumor_proliferating" and self.scenario["tumour_cycle"]:
                # Prototype rule: a cell annotated as proliferating at t0 is already
                # in cycle. The kernel finishes that cycle (two daughters); every
                # later cycle entry is a provider decision. Remaining time is a seeded draw.
                remaining = float(self.rng.uniform(self.p.get("t0_cycle_min_remaining_min",120), self.div["tumour_cycle_min"]))
                self.bio[g]["cycle"] = dict(start=0., due=remaining, duration=remaining, blocked=False, kind="tumor", origin="annotated_proliferating_t0")
        self.update_fields()
        self.refresh()
        self.initial_tumor_xy = t.xy[[g for g in np.flatnonzero(t.acts) if t.labels[g] in TUMOUR]].copy()
        self.initial_nest = union_all([t.polygons[g] for g in np.flatnonzero(t.acts) if t.labels[g] in TUMOUR])
        self.frames = []
        self.frame()

    def initial_state(self, label):
        return dict(viability="alive", damage_au=0., pulses=self.p["initial_pulses"] if label == CD8 else 0,
                    exhaustion=0., exhaustion_state="effector", glycolysis=self.p["initial_glycolysis"], pdl1=0.,
                    program=None, cycle=None, cooldown_until=0., death_due=None, hypoxia_min=0.,
                    last_stimulation=None, divisions_left=0, processed=0, next_refill=self.p["refill_interval_min"],
                    polarized_until=0., polarization="unassigned", perfused=False, tip=False,
                    perfusion_due=None, next_recruit=0., recruited_here=0, initial_phase=None, dead_since=None)

    def nerve_effect(self, name):
        return self.scenario["nerve"] != "none" and self.p["nerve_effects"].get(name, False)

    def nerve_geometry(self):
        if self.scenario["nerve"] == "none":
            return np.zeros((0,2))
        t = self.tissue
        tumor = [g for g in np.flatnonzero(t.acts) if t.labels[g] in TUMOUR]
        u = (t.xy[tumor].mean(0)-t.centre) if tumor else np.array([1.,0.])
        u /= max(float(np.linalg.norm(u)), 1e-9)
        tangent = np.array([-u[1],u[0]])
        support = float(((t.xy-t.centre)@u+t.rho).max()) + self.p["nerve_radius_um"] + 4
        centre = t.centre + u*support
        return np.array([centre-tangent*40, centre+tangent*40])

    def observation_cell(self, g):
        t = self.tissue
        cid = int(t.ids[g])
        if cid in self.by_id:
            return self.by_id[cid]
        parent = t.parents[g]
        if parent in self.by_id:
            return {**self.by_id[parent], "cell_id":cid, "role":"roi"}
        # Recruit competence is explicitly assigned by the influx scenario.
        return dict(cell_id=cid, role="roi", identity=dict(agreement="SIM_assigned", candidates=[dict(label=t.labels[g], lineage="SIM", source="scenario_assigned")]), signatures={})

    def aucell_state(self, cell):
        """Compact t0 state: AUCell rank percentile within provisional type and detection-depth decile."""
        sig = cell.get("signatures") or {}
        if not sig:
            return None
        return {name.replace("GOBP_",""): round(float(v["percentile_within_type_depth"]),2) for name,v in sorted(sig.items())}

    def nearest(self, g, accept, radius):
        """Nearest present cell within ``radius`` (um) satisfying ``accept(k)``."""
        t = self.tissue
        d = np.linalg.norm(t.xy-t.xy[g], axis=1)
        order = np.argsort(d)
        for k in order.tolist():
            if d[k] > radius: break
            if k != g and t.present[k] and accept(k): return k
        return None

    def update_fields(self):
        grid, t, p = self.fields, self.tissue, self.contract["fields"]
        if not self.scenario["fields"]:
            for name in grid.values:
                grid.values[name][:] = 1. if name in ("oxygen", "glucose") else 0.
            return
        idx = np.flatnonzero(t.present)
        xy = t.xy[idx]
        alive = np.array([self.bio[g]["viability"] == "alive" for g in idx])
        area = t.A0[idx]*alive
        gly = np.array([self.bio[g]["glycolysis"] for g in idx])
        perf = np.array([self.bio[g]["perfused"] for g in idx])*area
        exchange = grid.deposit(xy, perf*p["perfusion_exchange"])
        grid.solve("oxygen", exchange, exchange+grid.deposit(xy, area*p["oxygen_consumption"]*(1-gly)))
        grid.solve("glucose", exchange, exchange+grid.deposit(xy, area*p["glucose_consumption"]*(1+gly)))
        grid.solve("lactate", grid.deposit(xy, area*gly*p["lactate_production"]))
        source = {k: np.zeros(len(idx)) for k in ("vegf","cxcl9","ifng","tgfb","matrix")}
        for j,g in enumerate(idx):
            if not alive[j] or not t.acts[g]:
                continue
            s = self.bio[g]
            action = (s["program"] or {}).get("action")
            if action == "DEGRANULATE": source["ifng"][j] = 1
            if action == "SECRETE_VEGFA": source["vegf"][j] = 1
            if action == "SECRETE_VEGFA_TGFB": source["vegf"][j] = source["tgfb"][j] = 1
            if action == "SECRETE_CXCL9": source["cxcl9"][j] = 1
            if action == "SECRETE_TGFB": source["tgfb"][j] = 1
            if action == "DEPOSIT_MATRIX": source["matrix"][j] = 1
        # IFNg exposure licenses a declared downstream chemokine response.
        exposure = grid.sample("ifng", xy)
        responsive = np.array([t.labels[g] in TUMOUR or t.labels[g] == MACRO for g in idx])
        source["cxcl9"] += alive*responsive*(exposure >= self.p["signal_threshold"])
        for name, amounts in source.items():
            rate = p["matrix_source"] if name == "matrix" else p["secretory_source"]
            deposited = grid.deposit(xy, amounts*area*rate)
            if name == "matrix":
                # Deposited matrix persists and decays; it does not disappear
                # when a secretory program ends. Implicit Euler is conservative.
                dt = self.mech["tick_min"]
                grid.solve(name, deposited+grid.values[name]/dt, np.full(grid.shape,1/dt))
            else:
                grid.solve(name, deposited)
        nerve_source = np.zeros(grid.shape)
        if len(self.nerve):
            length = float(np.linalg.norm(self.nerve[1]-self.nerve[0]))
            count = int(math.ceil(length/grid.spacing))+1
            nerve_source = grid.deposit(np.linspace(*self.nerve, count), np.full(count,length/count*p["nerve_source_per_um"]))
        for name in ("ne","cgrp"):
            grid.solve(name, nerve_source if self.scenario["nerve"] == name else np.zeros(grid.shape))

    def refresh(self):
        t = self.tissue
        self.idx, self.geo, self.lp = self.local_pressure()
        self.graph = t.contact_pairs(self.idx, self.geo)
        self.neighbours = defaultdict(list)
        for a,b in self.graph:
            self.neighbours[a].append(b)
            self.neighbours[b].append(a)
        for g in self.neighbours:
            self.neighbours[g].sort(key=lambda k:int(t.ids[k]))
        arrays = {k:self.fields.sample(k,t.xy) for k in self.fields.values}
        self.local_fields = [{k:float(v[g]) for k,v in arrays.items()} for g in range(len(t.ids))]
        self.gradients = self.fields.gradient("cxcl9",t.xy)
        self.gradient_norm = np.linalg.norm(self.gradients,axis=1)

    def requests(self):
        t = self.tissue
        requests = []
        for g in np.flatnonzero(t.acts & t.present).tolist():
            s, label, vals = self.bio[g], t.labels[g], self.local_fields[g]
            if label not in self.scheduler.p["interval_min"] or s["viability"] != "alive":
                continue
            relevant = tuple((int(t.ids[k]),self.bio[k]["viability"]) for k in self.neighbours[g]
                             if (t.labels[k] != label and (t.labels[k] in TUMOUR or t.labels[k] == CD8)) or self.bio[k]["viability"] != "alive")
            bins = tuple(vals[k] >= self.p["signal_threshold"] for k in ("cxcl9","ifng","vegf","tgfb","ne","cgrp")) + (vals["oxygen"] < self.p["hypoxia_threshold"], vals["glucose"] < self.p["glucose_gate"])
            self.scheduler.observe(g,relevant,bins)
            if not self.scheduler.eligible(g,self.time,busy=s["program"] is not None):
                continue
            options, excluded = menu(self,g)
            if len(options) < 2:
                continue
            cell = self.observation_cell(g)
            contact_rows = [[int(t.ids[k]),t.labels[k],self.bio[k]["viability"],self.bio[k]["pdl1"],self.graph[tuple(sorted((g,k)))]] for k in self.neighbours[g][:self.scheduler.p["contact_limit"]]]
            obs = dict(time_min=self.time, primary_label=label, identity=cell["identity"],
                       t0_aucell_percentile=self.aucell_state(cell),
                       t0_aucell_raw={k.replace("GOBP_", ""): v["raw"] for k, v in (cell.get("signatures") or {}).items()},
                       t0_aucell_provenance=cell.get("observation_provenance","measured at t0 (Xenium AUCell)"),
                       sim_state={k:copy.deepcopy(s[k]) for k in ("viability","damage_au","pulses","exhaustion","glycolysis","pdl1","last_stimulation","divisions_left","polarization")},
                       fields_sim={k:round(v,4) for k,v in vals.items()}, local_pressure_sim=round(float(self.lp[g]),4),
                       contacts=dict(columns=["id","primary_label","viability_sim","pdl1_sim","provenance"],rows=contact_rows),
                       recognition_sim=self.scenario["recognition"],
                       adrenergic_cycle_bias_sim=bool(self.nerve_effect("ne_cycle_bias") and vals["ne"] >= self.p["signal_threshold"]),
                       excluded={k:[c["condition"]+":"+c["result"] for c in checks if c["result"] != "true"] for k,checks in excluded.items()})
            requests.append(dict(cell_id=int(t.ids[g]),index=g,options=options,observation=obs,wake_reasons=sorted(self.scheduler.pending[g])))
        return requests

    def decide(self, lp=None):
        if not self.scenario.get("autonomy",True):
            return
        requests = self.requests()
        if not requests:
            return
        # A pilot counts individual cell questions, not time steps or batches.
        # Bound the request list before inference, including the very first tick.
        available = len(requests)
        limits = []
        if self.max_decisions is not None:
            limits.append((self.max_decisions-len(self.decisions), "decision_cap_reached"))
        # Optional one-request-per-cell providers expose their remaining budget.
        if self.stop_on_provider_failure and hasattr(self.policy, "remaining_requests"):
            limits.append((self.policy.remaining_requests, "api_call_cap_reached"))
        if limits:
            remaining, cap_reason = min(limits, key=lambda item: item[0])
            if remaining <= 0:
                self.stop_reason = cap_reason
                return
            if len(requests) > remaining:
                # Stable pseudo-random ranks are independent of cell-ID order and
                # do not consume the RNG used for the biological trajectory.
                requests = sorted(requests, key=lambda row: fingerprint(dict(
                    seed=self.seed, time_min=self.time, cell_id=row["cell_id"],
                    purpose="bounded_pilot_selection")))[:remaining]
            self.decision_rounds.append(dict(time_min=self.time, available=available,
                                            selected=len(requests), omitted=available-len(requests)))
            if len(requests) == remaining:
                self.pending_cap_stop = cap_reason
        from .provider_jev import build_batches
        batches = build_batches(requests,self.contract["provider"]["pinned_model"],self.contract["provider"])
        self.batches_estimated.append(dict(time_min=self.time,questions=len(requests),batches=len(batches)))
        try:
            answers = self.policy.choose(copy.deepcopy(requests))
        except BudgetExhausted:
            self.stop_reason = "request_cap_before_round"
            return
        except Exception:
            if self.stop_on_provider_failure:
                self.stop_reason = "provider_failed"
                self.failures.append(dict(time_min=self.time, reason="provider_failed; pilot stopped without substitution"))
                return
            answers = {}
            self.failures.append(dict(time_min=self.time,reason="provider_failed; independent clocks continue"))
        if self.stop_on_provider_failure and not isinstance(answers, dict):
            self.stop_reason = "provider_failed"
            self.failures.append(dict(time_min=self.time, reason="provider returned a non-mapping response; pilot stopped"))
            return
        validated = []
        for r in requests:
            g, cid = r["index"], r["cell_id"]
            answer = answers.get(cid)
            status = "validated"
            try:
                validate_answer(answer,r["options"])
            except (ValueError,TypeError):
                status = "invalid_answer"
                self.failures.append(dict(time_min=self.time,cell_id=cid,reason=status))
            in_contact = any(self.tissue.labels[k] in TUMOUR for k in self.neighbours[g])
            reasons = self.scheduler.consume(g,self.time,self.tissue.labels[g],in_contact)
            record = dict(time_min=self.time,cell_id=cid,provider=self.policy.name,wake_reasons=reasons,
                          observation=r["observation"],options={k:{a:b for a,b in v.items() if a != "target_index"} for k,v in r["options"].items()},
                          answer=answer,status=status)
            self.decisions.append(record)
            if status == "validated": validated.append((g,r["options"][answer["choice"]],record))
        if self.stop_on_provider_failure and (len(validated) != len(requests) or getattr(self.policy,"failed",False)):
            if len(validated) == len(requests):
                self.failures.append(dict(time_min=self.time, reason="provider reported failure; no decisions applied"))
            for _, _, record in validated:
                record["status"] = "not_applied_provider_failure"
            self.stop_reason = "provider_failed"
            return
        # Common-snapshot validation precedes all writes; tie breaking is seeded.
        reserved = {s["program"]["target"] for s in self.bio if s["program"] and s["program"].get("target") is not None and s["program"]["action"] in ("DEGRANULATE","EFFEROCYTOSE","HELP_CD8")}
        for j in self.rng.permutation(len(validated)).tolist():
            g,option,record = validated[j]
            k = option["target_index"]
            if option["action"] in ("DEGRANULATE","EFFEROCYTOSE","HELP_CD8"):
                if k in reserved or g in reserved:
                    record["status"] = "reservation_conflict"
                    continue
                reserved.update((g,k))
            self.start(g,option)

    def start(self,g,option):
        s,p,action,k = self.bio[g],self.p,option["action"],option["target_index"]
        if action in ("WAIT","REMAIN_QUIESCENT","STAY_QUIESCENT"):
            # Tumor cell-cycle checkpoints are clocks; WAIT elsewhere does not poll.
            if self.tissue.labels[g] in TUMOUR:
                s["checkpoint_due"] = self.time+self.div["quiescent_recheck_min"]
            return
        if action == "ENTER_CYCLE":
            phase = s.pop("initial_phase",None) or 0.
            s["cycle"] = dict(start=self.time-phase,due=self.time+self.div["tumour_cycle_min"]-phase,duration=self.div["tumour_cycle_min"],blocked=False,kind="tumor")
            self.event("cycle_entered",g)
            return
        if action == "PROLIFERATE":
            s["divisions_left"] -= 1
            s["cycle"] = dict(start=self.time,due=self.time+p["cd8_cycle_min"],duration=p["cd8_cycle_min"],blocked=False,kind="cd8")
            self.event("helped_cd8_cycle_entered",g)
            return
        duration = p["program_min"]
        if action == "DEGRANULATE": duration = p["contact_program_min"]
        if action == "EFFEROCYTOSE": duration = p["clearance_min"]
        if action.startswith("SECRETE") or action == "DEPOSIT_MATRIX": duration = p["secretory_program_min"]
        if action == "SPROUT": duration = p["sprout_duration_min"]
        if action in ("APPROACH_TUMOR","APPROACH_CORPSE"): duration = p.get("approach_max_min",60)
        s["program"] = dict(action=action,target=k,target_id=None if k is None else int(self.tissue.ids[k]),start=self.time,due=self.time+duration)
        if action == "SPROUT":
            s["tip"] = True
            s["sprout_origin"] = self.tissue.xy[g].tolist()
            s["sprout_was_perfused"] = s["perfused"]
            s["sprout_initial_neighbors"] = [k for k in self.neighbours[g] if self.tissue.labels[k] == "Endothelial"]
            s["perfused"] = False
        self.heading.setdefault(g,float(self.rng.uniform(0,2*math.pi)))
        self.event("program_started",g,action=action,target_id=None if k is None else int(self.tissue.ids[k]),due_min=self.time+duration)

    def substep(self,dt):
        if not self.scenario["motility"] and not self.scenario["tumour_cycle"] and not self.scenario.get("autonomy",True):
            return dict(max_step=0.,capped=0)
        t = self.tissue
        idx,geo = t.geometry()
        movers = [g for g in idx.tolist() if self.bio[g]["viability"] == "alive" and self.bio[g]["program"] and self.bio[g]["program"]["action"] in ("SCAN_RANDOM","MIGRATE_UP_CXCL9","PATROL","DISENGAGE","SPROUT","APPROACH_TUMOR","APPROACH_CORPSE")]
        self.turn(movers,dt)
        for g in movers:
            prog = self.bio[g]["program"]
            if prog["action"] == "MIGRATE_UP_CXCL9" and self.gradient_norm[g] > self.p["chemotaxis_min_gradient"]:
                self.heading[g] = math.atan2(self.gradients[g,1],self.gradients[g,0])
            elif prog["action"] in ("APPROACH_TUMOR","APPROACH_CORPSE") and t.present[prog["target"]]:
                vec = t.xy[prog["target"]]-t.xy[g]
                self.heading[g] = math.atan2(vec[1],vec[0])
            elif prog["action"] == "DISENGAGE":
                vec = t.xy[g]-t.xy[prog["target"]]
                self.heading[g] = math.atan2(vec[1],vec[0])
        active = super().propulsion(movers,idx,geo)
        for j,g in enumerate(idx.tolist()):
            if g not in movers: continue
            s, vals = self.bio[g], self.local_fields[g]
            if s["program"]["action"] == "SPROUT":
                grad = self.fields.gradient("vegf",t.xy[[g]])[0]
                norm = float(np.linalg.norm(grad))
                active[j] = grad/max(norm,1e-12)*self.p["sprout_speed_um_per_min"]
            else:
                factor = (1-.8*s["exhaustion"])/(1+self.p["matrix_drag_factor"]*vals["matrix"]+vals["lactate"])
                if self.nerve_effect("ne_motility"): factor /= 1+self.p["ne_motility_factor"]*vals["ne"]
                active[j] *= factor
        return t.relax(idx,geo,dt,active)

    def commit_death(self,g,reason):
        s = self.bio[g]
        if s["viability"] != "alive": return
        s.update(viability="apoptosis_committed",death_due=self.time+self.p["death_commitment_min"],program=None,cycle=None)
        self.tissue.fates[g] = s["viability"]
        self.event("apoptosis_committed",g,cause=reason,damage_au=s["damage_au"])
        t = self.tissue
        near = np.flatnonzero(np.linalg.norm(t.xy-t.xy[g],axis=1) <= self.p.get("corpse_signal_radius_um",80))
        for k in near.tolist():
            if t.present[k] and t.acts[k] and t.labels[k] == MACRO and self.bio[k]["viability"] == "alive":
                self.scheduler.wake(k,"apoptotic_find_me_signal")

    def complete_programs(self):
        t,p = self.tissue,self.p
        for g in range(len(self.bio)):
            s = self.bio[g]
            prog = s["program"]
            if not t.present[g] or s["viability"] != "alive" or not prog: continue
            action,k = prog["action"],prog["target"]
            contact = k is not None and t.present[k] and tuple(sorted((g,k))) in self.graph
            approach = action in ("APPROACH_TUMOR","APPROACH_CORPSE")
            if approach:
                lost = k is None or not t.present[k] or (action == "APPROACH_TUMOR" and self.bio[k]["viability"] != "alive")
                if not (contact or lost or self.time >= prog["due"]): continue
            elif self.time < prog["due"]: continue
            completed = True
            if approach:
                completed = bool(contact)
                if contact: self.event("target_reached",g,action=action,target_id=int(t.ids[k]))
            if action == "DEGRANULATE":
                valid = contact and self.bio[k]["viability"] == "alive" and s["pulses"] > 0 and s["exhaustion"] < p["exhaustion_terminal"] and self.local_fields[g]["glucose"] >= p["glucose_gate"]
                if valid:
                    vals = self.local_fields[g]
                    checkpoint = self.bio[k]["pdl1"]*(p["checkpoint_blockade_factor"] if self.scenario["blockade"] else 1)
                    damage = p["pulse_damage_au"]*(1-s["exhaustion"])/(1+p["lactate_damage_factor"]*vals["lactate"]+checkpoint)
                    s["pulses"] -= 1
                    s["cooldown_until"] = self.time+p["cooldown_min"]
                    s["last_stimulation"] = self.time
                    self.bio[k]["damage_au"] += damage
                    self.event("cytotoxic_pulse",g,target_id=int(t.ids[k]),damage_au=damage,contact_provenance=self.graph[tuple(sorted((g,k)))])
                    if self.bio[k]["damage_au"] >= p["death_threshold_au"]: self.commit_death(k,"cytotoxic_damage")
                else: completed = False
            elif action == "EFFEROCYTOSE":
                if contact and self.bio[k]["viability"] in ("apoptosis_committed","dead_intact") and s["processed"] < p["clearance_capacity"]:
                    s["processed"] += 1
                    self.bio[k].update(viability="cleared",program=None,cycle=None)
                    t.remove(k,"cleared")
                    self.event("clearance",g,target_id=int(t.ids[k]))
                else: completed = False
            elif action == "HELP_CD8":
                if contact and self.bio[k]["viability"] == "alive":
                    self.bio[k]["divisions_left"] = max(self.bio[k]["divisions_left"],p["help_divisions"])
                    self.bio[k]["exhaustion"] = max(0.,self.bio[k]["exhaustion"]-p["cd4_help_recovery"])
                    self.scheduler.wake(k,"cd4_help_completed")
                    self.event("cd4_help",g,target_id=int(t.ids[k]))
                else: completed = False
            elif action == "UPREGULATE_PDL1": s["pdl1"] = min(1.,s["pdl1"]+p["pdl1_step"])
            elif action == "GLYCOLYTIC_SHIFT": s["glycolysis"] = min(1.,s["glycolysis"]+p["glycolysis_step"])
            elif action == "OXIDATIVE_SHIFT": s["glycolysis"] = max(0.,s["glycolysis"]-p["glycolysis_step"])
            elif action in ("SECRETE_CXCL9","SECRETE_VEGFA_TGFB") and t.labels[g] == MACRO:
                s["polarization"] = "inflammatory" if action == "SECRETE_CXCL9" else "suppressive"
                s["polarized_until"] = self.time+p["macrophage_dwell_min"]
            elif action == "SPROUT":
                connected = [k for k in self.neighbours[g] if t.labels[k] == "Endothelial" and self.bio[k]["perfused"] and k not in s["sprout_initial_neighbors"]]
                segment = dict(cell_id=int(t.ids[g]),start=s["sprout_origin"],end=t.xy[g].tolist(),time_min=self.time,
                               length_um=float(np.linalg.norm(t.xy[g]-s["sprout_origin"])),connected_to=[int(t.ids[k]) for k in connected],perfused=False)
                self.segments.append(segment)
                if connected:
                    s["perfusion_due"] = self.time+p["perfusion_delay_min"]
                    s["segment_index"] = len(self.segments)-1
                    self.event("anastomosis",g,partners=segment["connected_to"],perfusion_due=s["perfusion_due"])
                else:
                    self.event("sprout_unconnected",g)
                    s["tip"] = False
            s["program"] = None
            self.event("program_completed" if completed else "program_cancelled",g,action=action,target_id=None if k is None else int(t.ids[k]))
            self.scheduler.wake(g,"program_completed")

    def advance_biology(self,tick):
        t,p = self.tissue,self.p
        for g in np.flatnonzero(t.acts & t.present).tolist():
            s,vals = self.bio[g],self.local_fields[g]
            if s["viability"] == "apoptosis_committed" and self.time >= s["death_due"]:
                s["viability"] = t.fates[g] = "dead_intact"
                s["dead_since"] = self.time
                self.event("intact_death",g)
            if s["viability"] == "dead_intact" and self.time-s["dead_since"] >= self.p.get("corpse_degradation_min",720):
                # Fallback when no phagocyte arrives: secondary necrosis and fragmentation.
                s.update(viability="degraded",program=None,cycle=None)
                t.remove(g,"degraded")
                self.event("corpse_degraded",g)
                continue
            if s["viability"] != "alive": continue
            if t.labels[g] in TUMOUR:
                s["hypoxia_min"] = s["hypoxia_min"]+tick if vals["oxygen"] < p["necrosis_threshold"] else 0.
                if s["hypoxia_min"] >= p["necrosis_delay_min"]:
                    s["removed_area_um2"] = float(self.geo.area[np.flatnonzero(self.idx == g)[0]])
                    s.update(viability="necrotic",program=None,cycle=None)
                    t.remove(g,"necrotic")
                    self.event("hypoxic_necrosis",g)
                    continue
            if t.labels[g] == CD8:
                engaged = self.scenario["recognition"] is True and any(t.labels[k] in TUMOUR and self.bio[k]["viability"] == "alive" for k in self.neighbours[g])
                pdl1 = max([self.bio[k]["pdl1"] for k in self.neighbours[g] if t.labels[k] in TUMOUR]+[0.])
                r = p["exhaustion_per_hour"]
                rate = r["antigen"]*engaged+r["pdl1"]*pdl1*engaged*(p["checkpoint_blockade_factor"] if self.scenario["blockade"] else 1)
                rate += r["lactate"]*vals["lactate"]+r["tgfb"]*vals["tgfb"]
                if self.nerve_effect("ne_exhaustion"): rate += r["ne"]*vals["ne"]
                if self.nerve_effect("cgrp_exhaustion"): rate += r["cgrp"]*vals["cgrp"]
                rest = (s["program"] or {}).get("action") == "REST" and not engaged
                s["exhaustion"] = float(np.clip(s["exhaustion"]+tick/60*(rate-p["rest_recovery_per_hour"]*rest),0,1))
                old = s["exhaustion_state"]
                s["exhaustion_state"] = "terminal_exhausted" if s["exhaustion"] >= p["exhaustion_terminal"] else "progenitor_exhausted" if s["exhaustion"] >= p["exhaustion_tpex"] else "effector"
                if old != s["exhaustion_state"]: self.event("exhaustion_transition",g,state=s["exhaustion_state"])
                if self.time >= s["next_refill"]:
                    if s["exhaustion_state"] == "effector" and not engaged and s["pulses"] < p["initial_pulses"]:
                        s["pulses"] += 1
                        self.event("effector_material_refilled",g)
                    s["next_refill"] = self.time+p["refill_interval_min"]
            if s.get("checkpoint_due") is not None and self.time >= s["checkpoint_due"]:
                self.scheduler.wake(g,"cycle_checkpoint")
                s.pop("checkpoint_due")
            if s["perfusion_due"] is not None and self.time >= s["perfusion_due"]:
                s["perfused"] = True
                s["perfusion_due"] = None
                self.segments[s["segment_index"]]["perfused"] = True
                self.event("perfusion_established",g)
        self.complete_programs()
        self.advance_autonomous_cycles()
        self.recruit()

    def add_state(self,g,parent=None):
        t = self.tissue
        s = self.initial_state(t.labels[g])
        if parent is not None:
            prior = self.bio[parent]
            for key in ("glycolysis","pdl1","exhaustion","exhaustion_state"):
                s[key] = prior[key]
            if t.labels[g] == CD8:
                # Division partitions finite effector material; it never doubles it.
                s["pulses"] = prior["pulses"]//2
                prior["pulses"] -= s["pulses"]
            root_cell = self.observation_cell(parent)
            self.by_id[int(t.ids[g])] = {**root_cell,"cell_id":int(t.ids[g]),"role":"roi","observation_provenance":"inherited from parent; not remeasured"}
        self.bio.append(s)
        self.scheduler.wake(g,"new_acting_cell")

    def advance_autonomous_cycles(self):
        t = self.tissue
        accepted=None
        if t.compression:
            proposals={g:t.A0_base[g]*(1+min(1.,(self.time-s["cycle"]["start"])/s["cycle"]["duration"]))
                       for g,s in enumerate(self.bio) if s["cycle"] and t.present[g] and s["viability"]=="alive"}
            accepted=t.grow_many(proposals)
        for g in range(len(self.bio)):
            s = self.bio[g]
            c = s["cycle"]
            if not c or not t.present[g] or s["viability"] != "alive": continue
            fraction = min(1.,(self.time-c["start"])/c["duration"])
            growth_allowed = g in accepted if accepted is not None else t.grow(g,t.A0_base[g]*(1+fraction))
            if not growth_allowed:
                c["start"]+=self.mech["tick_min"];c["due"]+=self.mech["tick_min"]
                if not c.get("growth_blocked"):
                    self.event("growth_blocked_compression",g)
                c["growth_blocked"]=True
                continue
            c["growth_blocked"]=False
            if self.time < c["due"]: continue
            if self.lp[g] > self.div["division_blocked_local_pressure"] or len(t.ids) >= self.p["max_population"]:
                if not c["blocked"]: self.event("division_blocked_crowding",g,local_pressure_sim=float(self.lp[g]))
                c["blocked"] = True
                continue
            angle = float(self.rng.uniform(0,2*math.pi))
            # Use the local long axis when the cell remains visible.
            local = np.flatnonzero(self.idx == g)
            polygon = self.geo.polygon(int(local[0])) if len(local) else None
            axis = np.array([math.cos(angle),math.sin(angle)])
            if polygon is not None and not polygon.is_empty and polygon.area > 1e-9:
                pts = np.asarray(polygon.minimum_rotated_rectangle.exterior.coords)[:3]
                edges = np.diff(pts,axis=0)
                axis = edges[np.argmax(np.linalg.norm(edges,axis=1))]
            if not t.can_divide(g,axis):
                if not c["blocked"]:self.event("division_blocked_compression",g)
                c["blocked"]=True
                continue
            daughter = t.divide(g,axis,self.time)
            self.add_state(daughter,g)
            s["cycle"] = None
            s["initial_phase"] = None
            self.scheduler.wake(g,"cycle_completed")
            self.event("division",g,daughter_id=int(t.ids[daughter]),cell_type=t.labels[g])

    def recruit(self):
        if not self.scenario["recruitment"]: return
        t,p = self.tissue,self.p
        templates = [g for g,l in enumerate(t.labels) if l == CD8 and t.origins[g] == "roi"]
        if not templates: return
        for g in np.flatnonzero(t.acts & t.present).tolist():
            s = self.bio[g]
            if t.labels[g] != "Endothelial" or not s["perfused"] or s["viability"] != "alive": continue
            if self.time < s["next_recruit"] or s["recruited_here"] >= p["recruitment_max_per_vessel"]: continue
            if g >= len(self.local_fields) or self.local_fields[g]["cxcl9"] < p["signal_threshold"] or self.lp[g] > p["recruitment_pressure_gate"]: continue
            if len(t.ids) >= p["max_population"]: break
            model = templates[0]
            angle = float(self.rng.uniform(0,2*math.pi))
            xy = t.xy[g]+4*np.array([math.cos(angle),math.sin(angle)])
            area=float(np.median(t.A0_base[templates]))
            if t.compression:
                radius=.9*(math.sqrt(area/math.pi)+math.sqrt(t.A0_base[g]/math.pi))
                candidates=[t.xy[g]+radius*np.array([math.cos(a),math.sin(a)]) for a in angle+np.arange(8)*math.pi/4]
                xy=next((position for position in candidates if np.linalg.norm(position-t.centre)<=t.fixed_radius and t.can_insert(CD8,position,area,t.kappa[model])),None)
                s["next_recruit"]=self.time+p["recruitment_interval_min"]
                if xy is None:
                    self.event("recruitment_blocked_compression",g)
                    continue
            h = t.add_body(CD8,"recruited",xy,area,t.kappa[model],None,self.time)
            self.add_state(h)
            s["next_recruit"] = self.time+p["recruitment_interval_min"]
            s["recruited_here"] += 1
            self.event("cd8_recruited",g,recruit_id=int(t.ids[h]),evidence="SIM influx, not measured")

    def accounting(self):
        t = self.tissue
        counts = Counter()
        for g in np.flatnonzero(t.acts):
            s = self.bio[g]
            key = ("alive_in_roi" if t.in_roi[g] else "alive_outside_roi") if s["viability"] == "alive" else s["viability"]
            counts[key] += 1
        return {**t.accounting(),"acting":dict(counts)}

    def step(self):
        if self.stop_reason: return
        if self.time >= self.run_p["horizon_min"]:
            self.stop_reason = "horizon_reached"
            return
        self.refresh()
        self.decide()
        if self.stop_reason: return
        tick,dt = self.mech["tick_min"],self.mech["substep_min"]
        stats = [self.substep(dt) for _ in range(int(round(tick/dt)))]
        self.time += tick
        self.refresh()
        if self.scenario.get("autonomy",True): self.advance_biology(tick)
        for g in self.tissue.update_roi():
            self.event("roi_entered" if self.tissue.in_roi[g] else "roi_exited",int(g))
            if self.tissue.in_roi[g]: self.scheduler.wake(int(g),"roi_entry")
        self.tissue.check_invariants()
        self.check_biology()
        self.update_fields()
        self.refresh()
        if self.tissue.compression:
            floors=self.tissue.floor_fractions()[self.idx]*self.tissue.A0_base[self.idx]
            if not self.tissue.valid_areas(self.geo,floors).all():raise AssertionError("Post-clock compression constraint violated")
        self.metrics.append(self.readout(self.idx,self.geo,self.lp,stats,self.accounting()))
        if self.time >= self.run_p["horizon_min"]: self.stop_reason = "horizon_reached"
        # Finish the current physics tick, then stop before any further query.
        # The final replay therefore includes the consequence of the last choices.
        if self.pending_cap_stop: self.stop_reason = self.pending_cap_stop
        if self.time % self.run_p["frame_interval_min"] == 0 or self.stop_reason: self.frame()

    def check_biology(self):
        t = self.tissue
        if len(self.bio) != len(t.ids): raise AssertionError("A physical cell has no biological state")
        for g,s in enumerate(self.bio):
            if not 0 <= s["pulses"] <= self.p["initial_pulses"] or int(s["pulses"]) != s["pulses"]:
                raise AssertionError("Effector material invariant violated")
            if not 0 <= s["processed"] <= self.p["clearance_capacity"]: raise AssertionError("Clearance capacity violated")
            for key in ("exhaustion","glycolysis","pdl1"):
                if not math.isfinite(s[key]) or not 0 <= s[key] <= 1: raise AssertionError("Invalid bounded SIM state")
            if not math.isfinite(s["damage_au"]) or s["damage_au"] < 0: raise AssertionError("Invalid damage")
            if t.present[g] != (s["viability"] not in ("cleared","necrotic","degraded")): raise AssertionError("Physical presence and fate disagree")
            if not t.acts[g] and (s["program"] or s["cycle"]): raise AssertionError("Context acquired an autonomous program")
        acc = self.accounting()
        if sum(acc["acting"].values()) != acc["expected_acting_total"]: raise AssertionError("Fate accounting does not balance")

    def readout(self,idx,geometry,lp,stats,accounting):
        base = super().readout(idx,geometry,lp,stats,accounting)
        t = self.tissue
        tumors = [g for g in idx if t.acts[g] and t.labels[g] in TUMOUR and self.bio[g]["viability"] == "alive"]
        cd8 = [g for g in idx if t.acts[g] and t.labels[g] == CD8 and self.bio[g]["viability"] == "alive"]
        areas = dict(zip(idx,geometry.area.tolist()))
        hypoxic = [g for g in tumors if self.local_fields[g]["oxygen"] < self.p["hypoxia_threshold"]]
        tumor_area = sum(areas[g] for g in tumors)
        base.update(tumour_cells=len(tumors),tumour_in_roi=sum(bool(t.in_roi[g]) for g in tumors),
                    cd8_tumour_contacts=dict(Counter(provenance for (a,b),provenance in self.graph.items() if (a in cd8 and b in tumors) or (b in cd8 and a in tumors))),
                    tumour_area_um2=float(tumor_area), hypoxic_tumour_area_um2=float(sum(areas[g] for g in hypoxic)),
                    hypoxic_tumour_fraction=len(hypoxic)/max(1,len(tumors)),
                    oxygen_tumour_mean=float(np.mean([self.local_fields[g]["oxygen"] for g in tumors])) if tumors else None,
                    glycolysis_tumour_mean=float(np.mean([self.bio[g]["glycolysis"] for g in tumors])) if tumors else None,
                    lactate_tumour_mean=float(np.mean([self.local_fields[g]["lactate"] for g in tumors])) if tumors else None,
                    cd8_exhaustion=dict(Counter(self.bio[g]["exhaustion_state"] for g in cd8)),
                    recruited=t.recruited, kills=sum(e["kind"] == "apoptosis_committed" for e in self.events),
                    vessel_length_um=sum(s["length_um"] for s in self.segments),
                    perfused_vessel_length_um=sum(s["length_um"] for s in self.segments if s["perfused"]),
                    necrotic_area_removed_um2=sum(s.get("removed_area_um2",0) for s in self.bio if s["viability"] == "necrotic"),
                    field_balance_max=max([abs(a["balance_error"]) for a in self.fields.audit.values()]+[0.]),
                    field_residual_max=max([abs(a["residual_max"]) for a in self.fields.audit.values()]+[0.]),
                    decisions=len(self.decisions),invalid_answers=len(self.failures))
        if t.compression:
            floors=t.floor_fractions()[idx]*t.A0_base[idx]
            base.update(area_min_over_floor=float((geometry.area/floors).min()),
                        compression_backtracks=t.constraint_rejections,compression_trials=t.constraint_trials,
                        compression_local_freezes=t.local_freezes,compression_global_rollbacks=t.global_rollbacks,
                        growth_steps_blocked=t.growth_blocked,insertion_trials_blocked=t.insertion_blocked)
        if cd8 and len(self.initial_tumor_xy):
            distances = cKDTree(self.initial_tumor_xy).query(t.xy[cd8])[0]
            base["cd8_near_initial_nest_20um"] = int((distances <= 20).sum())
            pts = points(t.xy[cd8])
            inside = contains(self.initial_nest,pts)
            depths = distance(pts[inside],self.initial_nest.boundary)
            base["cd8_inside_initial_tumour_footprint"] = int(inside.sum())
            base["cd8_initial_footprint_depth_q90_um"] = float(np.quantile(depths,.9)) if len(depths) else 0.
        return base

    def frame(self):
        if not hasattr(self,"bio"): return
        t = self.tissue
        idx = np.flatnonzero(t.present)
        self.frames.append(dict(time_min=self.time,ids=t.ids[idx].tolist(),x=t.xy[idx,0].tolist(),y=t.xy[idx,1].tolist(),rho=t.rho[idx].tolist(),
                                states=[{k:copy.deepcopy(self.bio[g][k]) for k in ("viability","damage_au","pulses","exhaustion","exhaustion_state","glycolysis","pdl1","perfused","program","cycle")} for g in idx],
                                fields=self.fields.frame(),segments=copy.deepcopy(self.segments),
                                metrics=copy.deepcopy(self.metrics[-1]) if self.metrics else dict(time_min=0,tumour_cells=sum(self.tissue.labels[g] in TUMOUR for g in np.flatnonzero(t.acts)),acting=self.accounting()["acting"])))

    def summary(self):
        counts = Counter(d["observation"]["primary_label"] for d in self.decisions)
        return dict(time_min=self.time,stop_reason=self.stop_reason,accounting=self.accounting(),events=dict(Counter(e["kind"] for e in self.events)),
                    decisions=len(self.decisions),decisions_by_type=dict(counts),invalid_answers=len(self.failures),
                    estimated_requests=sum(b["batches"] for b in self.batches_estimated),
                    paid_calls=getattr(self.policy,"requests",0), final=self.metrics[-1] if self.metrics else None,
                    local_inference_calls=getattr(self.policy,"inference_calls",0),
                    decision_budget=dict(limit=self.max_decisions, reached=self.max_decisions is not None and len(self.decisions)>=self.max_decisions,
                                         partial_round=any(r["omitted"] for r in self.decision_rounds),
                                         selection="seeded_hash_rank_without_replacement", rounds=copy.deepcopy(self.decision_rounds)),
                    decision_share_above_half=[k for k,n in counts.items() if n > len(self.decisions)/2],
                    uncalibrated=True)

    def episode(self):
        if not self.frames or self.frames[-1]["time_min"] != self.time: self.frame()
        t = self.tissue
        bodies = [dict(id=int(t.ids[g]),label=t.labels[g],origin=t.origins[g],parent=t.parents[g],fixed=bool(t.fixed[g]),fate=self.bio[g]["viability"],born_at_min=float(t.born_at[g])) for g in range(len(t.ids))]
        final = dict(ids=t.ids.tolist(),xy=t.xy.tolist(),present=t.present.tolist(),state=copy.deepcopy(self.bio),rng=copy.deepcopy(self.rng.bit_generator.state))
        data = dict(format="luad-episode-v2",patch_id=self.atlas["patch_id"],atlas_fingerprint=self.atlas["fingerprint"],contract=self.contract,engine_fingerprint=self.engine_fingerprint,
                    seed=self.seed,scenario=self.scenario_id,provider=self.policy.name,paid_calls=getattr(self.policy,"requests",0),
                    t0=dict(display="measured polygons",fidelity=self.fidelity),field_grid=self.fields.metadata(),nerve=self.nerve.tolist(),
                    bodies=bodies,frames=self.frames,events=self.events,decisions=self.decisions,failures=self.failures,
                    metrics=self.metrics,summary=self.summary(),batch_estimate=self.batches_estimated,provider_audit=list(self.policy.audit),
                    final_state=final,final_state_fingerprint=fingerprint(final),
                    scope="Uncalibrated SIM mechanism experiment on one measured Xenium section. Fields, recognition, rates, perfusion and virtual nerve are declared assumptions; not patient predictions.")
        return {**data,"fingerprint":fingerprint(data)}
