"""M1 runtime: moving geometry, neighbour-swap migration, growth and division.

Order within a 10-minute tick:
1. decisions for cells whose checkpoint is due (menus gated by the kernel);
2. mechanics substeps (power diagram, pressure relaxation, CD8 swaps);
3. kernel clocks (growth, division or blocked division);
4. ROI tracking, accounting, invariants, readouts, hourly frame.
Only the kernel changes geometry, clocks and fates. The policy only chooses
among options the kernel already validated.
"""
import copy
import math
from collections import Counter

import numpy as np

from ..io import fingerprint
from .fidelity import fidelity_report, unfitted_report
from .local_policy import M1ExercisePolicy
from .mechanics import TUMOUR, Tissue

CD8 = "T_CD8_cytotoxic"


class SimulationV2:
    def __init__(self, atlas, contract, scenario="m1_exercise", seed=21, policy=None):
        if scenario not in contract["scenarios"]:
            raise ValueError("Unsupported scenario")
        self.atlas = atlas
        self.contract = copy.deepcopy(contract)
        self.scenario_id = scenario
        self.scenario = contract["scenarios"][scenario]
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.policy = policy or M1ExercisePolicy()
        self.mech = contract["mechanics"]
        self.mig = contract["migration"]
        self.div = contract["division"]
        self.run_p = contract["run"]
        self.tissue = Tissue.from_atlas(atlas, contract)
        self.time = 0
        self.stop_reason = None
        self.events, self.decisions, self.metrics, self.frames = [], [], [], []
        self.fidelity = dict(fitted=fidelity_report(atlas, self.tissue), unfitted=unfitted_report(atlas, self.tissue),
                             weight_fit=self.tissue.fit)
        self.programs = {}
        self.heading = {}
        t = self.tissue
        for g in np.flatnonzero(t.acts):
            label = t.labels[g]
            if label == CD8 and self.scenario["cd8_motility"]:
                self.programs[g] = dict(kind="motility", mode=None, next_check=0)
                self.heading[g] = float(self.rng.uniform(0, 2 * math.pi))
            elif label in TUMOUR and self.scenario["tumour_cycle"]:
                if label == "Tumor_proliferating":
                    # The annotation marks proliferation; the phase is a seeded SIM draw.
                    elapsed = float(self.rng.uniform(0, self.div["tumour_cycle_min"]))
                    self.programs[g] = dict(kind="cycle", phase="cycling", start=-elapsed,
                                            due=self.div["tumour_cycle_min"] - elapsed, blocked=False, next_check=None)
                else:
                    self.programs[g] = dict(kind="cycle", phase="quiescent", next_check=0)
        self.swap_count = 0
        self.frame()

    # --------------------------------------------------------------- helpers
    def event(self, kind, g=None, **fields):
        cid = None if g is None else int(self.tissue.ids[g])
        self.events.append(dict(index=len(self.events), time_min=self.time, kind=kind, cell_id=cid, **fields))

    def local_pressure(self):
        idx, geometry = self.tissue.geometry()
        lp = np.zeros(len(self.tissue.ids))
        lp[idx] = self.tissue.local_pressure(idx, geometry)
        return idx, geometry, lp

    # ------------------------------------------------------------- decisions
    def decide(self, lp):
        t = self.tissue
        requests = []
        for g, prog in sorted(self.programs.items()):
            if not t.present[g] or prog.get("next_check") is None or self.time < prog["next_check"]:
                continue
            if prog["kind"] == "motility":
                options = {"scan": "persistent random motility at the declared speed", "hold": "stay in place"}
            elif prog["kind"] == "cycle" and prog["phase"] == "quiescent":
                gate = lp[g] <= self.div["cycle_entry_max_local_pressure"]
                prog["next_check"] = self.time + self.div["quiescent_recheck_min"]
                if not gate:
                    continue  # a single legal option is never queried
                options = {"enter_cycle": "start the accelerated SIM cycle clock", "remain_quiescent": "stay in G0"}
            else:
                continue
            requests.append(dict(cell_id=int(t.ids[g]), index=int(g), kind=prog["kind"], options=options,
                                 observation=dict(local_pressure_sim=round(float(lp[g]), 3), label=t.labels[g])))
        if not requests:
            return
        answers = self.policy.choose(copy.deepcopy(requests))
        for request in requests:
            g, answer = request["index"], answers[request["cell_id"]]
            if answer["choice"] not in request["options"]:
                raise ValueError("Policy returned an option that was not offered")
            prog = self.programs[g]
            if prog["kind"] == "motility":
                prog.update(mode=answer["choice"], next_check=None)
            elif answer["choice"] == "enter_cycle":
                prog.update(phase="cycling", start=self.time, due=self.time + self.div["tumour_cycle_min"], blocked=False, next_check=None)
                self.event("cycle_entered", g)
            self.decisions.append(dict(time_min=self.time, provider=self.policy.name, cell_id=request["cell_id"],
                                       kind=request["kind"], options=sorted(request["options"]), answer=answer))

    # ------------------------------------------------------------- mechanics
    def turn(self, movers, dt):
        for g in movers:
            if self.rng.random() < dt / self.mig["persistence_min"]:
                self.heading[g] = float(self.rng.uniform(0, 2 * math.pi))

    def substep(self, dt):
        t = self.tissue
        idx, geometry = t.geometry()
        movers = [g for g, p in sorted(self.programs.items()) if p["kind"] == "motility" and p["mode"] == "scan" and t.present[g]]
        if self.mig["mode"] == "self_propelled":
            self.turn(movers, dt)
            return t.relax(idx, geometry, dt, self.propulsion(movers, idx, geometry))
        stats = t.relax(idx, geometry, dt)
        if not movers:
            return stats
        return self.swap_substep(movers, idx, geometry, dt, stats)

    def propulsion(self, movers, idx, geometry):
        """Self-propulsion v0 * heading * traction.

        traction = sum_j L_ij * permeability_j / max(sum_j L_ij, f * 2 * pi * r_eq):
        crawling needs neighbours, so a cell exposed to a void slows down and an
        isolated cell cannot propel itself (f = substrate_perimeter_fraction).
        """
        t = self.tissue
        active = np.zeros((len(idx), 2))
        if not movers:
            return active
        local = np.full(len(t.ids), -1)
        local[idx] = np.arange(len(idx))
        i, j, L = geometry.pair_i, geometry.pair_j, geometry.pair_length
        perm = np.array([self.mig["permeability"][t.labels[g]] if not t.fixed[g] else 0. for g in idx.tolist()])
        n = len(idx)
        num = np.bincount(i, L * perm[j], n) + np.bincount(j, L * perm[i], n)
        den = np.bincount(i, L, n) + np.bincount(j, L, n)
        r_eq = np.sqrt(t.A0[idx] / math.pi)
        factor = num / np.maximum(den, self.mig["substrate_perimeter_fraction"] * 2 * math.pi * r_eq)
        v0 = self.mig["cd8_declared_speed_um_per_min"]
        for g in movers:
            k = local[g]
            active[k] = v0 * factor[k] * np.array([math.cos(self.heading[g]), math.sin(self.heading[g])])
        return active

    def swap_substep(self, movers, idx, geometry, dt, stats):
        t = self.tissue
        local = np.full(len(t.ids), -1)
        local[idx] = np.arange(len(idx))
        i, j, _ = geometry.contacts(t.g["contact_max_gap_um"])
        gi, gj = idx[i], idx[j]
        is_mover = np.zeros(len(t.ids), bool)
        is_mover[movers] = True
        sel = is_mover[gi] | is_mover[gj]
        neighbours = {g: [] for g in movers}
        for a, b in zip(gi[sel].tolist(), gj[sel].tolist()):
            if a in neighbours:
                neighbours[a].append(b)
            if b in neighbours:
                neighbours[b].append(a)
        reserved = set()
        self.turn(movers, dt)
        for g in self.rng.permutation(np.array(movers)).tolist():
            if g in reserved:
                continue
            candidates = [k for k in neighbours[g] if not t.fixed[k] and k not in reserved]
            if not candidates:
                continue
            vec = t.xy[candidates] - t.xy[g]
            dist = np.hypot(*vec.T)
            if self.rng.random() >= min(1., self.mig["cd8_declared_speed_um_per_min"] * dt / float(dist.mean())):
                continue
            heading = np.array([math.cos(self.heading[g]), math.sin(self.heading[g])])
            weight = np.exp(self.mig["direction_bias_beta"] * (vec @ heading) / dist)
            k = candidates[int(self.rng.choice(len(candidates), p=weight / weight.sum()))]
            if self.rng.random() >= self.mig["permeability"][t.labels[k]]:
                continue
            t.swap(g, k)
            reserved.update((g, k))
            self.swap_count += 1
        return stats

    # ----------------------------------------------------------------- clocks
    def advance_cycles(self, lp, geometry, idx):
        t = self.tissue
        local = {int(g): k for k, g in enumerate(idx)}
        T = self.div["tumour_cycle_min"]
        for g in sorted(self.programs):
            prog = self.programs[g]
            if prog["kind"] != "cycle" or prog["phase"] != "cycling" or not t.present[g]:
                continue
            progress = min(1., (self.time - prog["start"]) / T)
            t.A0[g] = t.A0_base[g] * (1 + progress)
            if self.time < prog["due"]:
                continue
            if lp[g] > self.div["division_blocked_local_pressure"]:
                if not prog["blocked"]:
                    self.event("division_blocked_crowding", g, local_pressure_sim=round(float(lp[g]), 3))
                prog["blocked"] = True
                prog["blocked_ticks"] = prog.get("blocked_ticks", 0) + 1
                continue
            polygon = geometry.polygon(local[g]) if g in local else None
            axis = None
            if polygon is not None and not polygon.is_empty and polygon.area > 0:
                box = np.asarray(polygon.minimum_rotated_rectangle.exterior.coords)[:3]
                e1, e2 = box[1] - box[0], box[2] - box[1]
                axis = e1 if np.hypot(*e1) >= np.hypot(*e2) else e2
            if axis is None or np.hypot(*axis) < 1e-9:
                angle = float(self.rng.uniform(0, 2 * math.pi))
                axis = np.array([math.cos(angle), math.sin(angle)])
            daughter = t.divide(g, axis, self.time)
            self.event("division", g, daughter_id=int(t.ids[daughter]), blocked_ticks=prog.get("blocked_ticks", 0))
            for h in (g, daughter):
                self.programs[h] = dict(kind="cycle", phase="quiescent", next_check=self.time + self.mech["tick_min"])

    # ------------------------------------------------------------------- loop
    def step(self):
        if self.stop_reason:
            return
        if self.time >= self.run_p["horizon_min"]:
            self.stop_reason = "horizon_reached"
            return
        _, _, lp = self.local_pressure()
        self.decide(lp)
        tick, dt = self.mech["tick_min"], self.mech["substep_min"]
        stats = [self.substep(dt) for _ in range(int(round(tick / dt)))]
        self.time += tick
        idx, geometry, lp = self.local_pressure()
        self.advance_cycles(lp, geometry, idx)
        for g in self.tissue.update_roi():
            self.event("roi_entered" if self.tissue.in_roi[g] else "roi_exited", g)
        accounting = self.tissue.check_invariants()
        idx, geometry, lp = self.local_pressure()
        self.metrics.append(self.readout(idx, geometry, lp, stats, accounting))
        if self.time % self.run_p["frame_interval_min"] == 0:
            self.frame()
        if self.time >= self.run_p["horizon_min"]:
            self.stop_reason = "horizon_reached"

    def run(self, horizon_min=None):
        if horizon_min is not None:
            if horizon_min > self.run_p["max_horizon_min"]:
                raise ValueError("Horizon exceeds the declared maximum")
            self.run_p = {**self.run_p, "horizon_min": horizon_min}
            self.contract["run"] = self.run_p
        while not self.stop_reason:
            self.step()
        return self

    # --------------------------------------------------------------- readouts
    def readout(self, idx, geometry, lp, stats, accounting):
        t = self.tissue
        P = t.pressure(idx, geometry)
        acting_present = [g for g in idx.tolist() if t.acts[g]]
        tumour = [g for g in acting_present if t.labels[g] in TUMOUR]
        cd8 = [g for g in acting_present if t.labels[g] == CD8]
        graph = t.contact_pairs(idx, geometry)
        cd8_tumour = Counter(v for (a, b), v in graph.items()
                             if {t.labels[a], t.labels[b]} & {CD8} and {t.labels[a], t.labels[b]} & TUMOUR)
        all_tumour = [g for g in idx.tolist() if t.labels[g] in TUMOUR]
        near = 0
        if cd8 and all_tumour:
            d = np.hypot(*(t.xy[cd8][:, None, :] - t.xy[all_tumour][None, :, :]).transpose(2, 0, 1))
            near = int((d.min(1) <= 20).sum())
        return dict(time_min=self.time,
                    acting=accounting["acting"], born=accounting["born"],
                    tumour_cells=len(tumour), tumour_in_roi=int(sum(t.in_roi[g] for g in tumour)),
                    tumour_area_um2=round(float(sum(geometry.area[k] for k, g in enumerate(idx.tolist()) if g in set(tumour))), 1),
                    cd8_tumour_contacts=dict(cd8_tumour), cd8_within_20um_of_tumour=near,
                    swaps_total=self.swap_count,
                    pressure_mean=round(float(P.mean()), 5), pressure_max=round(float(P.max()), 4),
                    area_total_over_target=round(float(geometry.area.sum() / t.A0[idx].sum()), 5),
                    area_min_over_target=round(float((geometry.area / t.A0[idx]).min()), 4),
                    covered_area_um2=round(float(geometry.area.sum()), 1),
                    compressed_below_30pct=dict(Counter(t.labels[g] for k, g in enumerate(idx.tolist()) if geometry.area[k] < .3 * t.A0[g])),
                    at_wall=dict(Counter(t.labels[g] for g in idx.tolist() if not t.fixed[g] and np.hypot(*(t.xy[g] - t.centre)) >= t.fixed_radius - 1e-6)),
                    hidden=int(geometry.hidden.sum()), rescues_total=t.rescues, confinements_total=t.confinements,
                    max_substep_um=round(max(s["max_step"] for s in stats), 4), capped_steps=sum(s["capped"] for s in stats),
                    pressure_by_label={label: round(float(P[[k for k, g in enumerate(idx.tolist()) if t.labels[g] == label]].mean()), 4)
                                       for label in sorted(set(t.labels[g] for g in idx.tolist()))})

    def frame(self):
        t = self.tissue
        idx = np.flatnonzero(t.present)
        self.frames.append(dict(time_min=self.time, ids=t.ids[idx].tolist(),
                                x=np.round(t.xy[idx, 0], 3).tolist(), y=np.round(t.xy[idx, 1], 3).tolist(),
                                rho=np.round(t.rho[idx], 4).tolist()))

    def displacement_by_label(self):
        """Median net displacement (um) of original, non-self-propelled bodies: an architecture-mixing diagnostic."""
        t = self.tissue
        movers = {g for g, p in self.programs.items() if p["kind"] == "motility" and p.get("mode") == "scan"}
        groups = {}
        for g in np.flatnonzero(t.present):
            if g in movers or t.origins[g] == "born" or t.fixed[g]:
                continue
            groups.setdefault(t.labels[g], []).append(float(np.hypot(*(t.xy[g] - t.xy0[g]))))
        return {k: dict(n=len(v), median=round(float(np.median(v)), 2), q90=round(float(np.quantile(v, .9)), 2)) for k, v in sorted(groups.items())}

    def summary(self):
        t = self.tissue
        movers = sorted(g for g, p in self.programs.items() if p["kind"] == "motility" and p.get("mode") == "scan" and t.present[g])
        cd8 = movers
        net = np.hypot(*(t.xy[cd8] - t.xy0[cd8]).T) if cd8 else np.zeros(0)
        hours = max(self.time / 60, 1e-9)
        return dict(time_min=self.time, stop_reason=self.stop_reason, accounting=t.accounting(),
                    events=dict(Counter(e["kind"] for e in self.events)), decisions=len(self.decisions),
                    swaps=self.swap_count,
                    cd8_path_um_per_min_mean=round(float(t.path_um[cd8].mean() / self.time), 4) if cd8 and self.time else 0.,
                    cd8_net_displacement_um=dict(median=round(float(np.median(net)), 2), q90=round(float(np.quantile(net, .9)), 2)) if cd8 else None,
                    cd8_declared_speed_um_per_min=self.mig["cd8_declared_speed_um_per_min"], migration_mode=self.mig["mode"],
                    passive_displacement_um=self.displacement_by_label(),
                    separations=t.separations,
                    final=self.metrics[-1] if self.metrics else None, hours=round(hours, 3))

    def episode(self):
        t = self.tissue
        bodies = [dict(id=int(t.ids[g]), label=t.labels[g], origin=t.origins[g], parent=t.parents[g],
                       fixed=bool(t.fixed[g]), fate=t.fates[g],
                       measured_area_um2=None if math.isnan(t.measured_area[g]) else round(float(t.measured_area[g]), 3),
                       born_at_min=float(t.born_at[g])) for g in range(len(t.ids))]
        data = dict(format="luad-episode-v2-m1", patch_id=self.atlas["patch_id"], atlas_fingerprint=self.atlas["fingerprint"],
                    contract=self.contract, seed=self.seed, scenario=self.scenario_id, provider=self.policy.name,
                    paid_calls=0, t0=dict(display="measured polygons", fidelity=self.fidelity),
                    roi=dict(centre_um=t.centre.round(4).tolist(), radius_um=round(t.roi_radius, 4), fixed_radius_um=round(t.fixed_radius, 4)),
                    bodies=bodies, frames=self.frames, events=self.events, decisions=self.decisions, metrics=self.metrics,
                    summary=self.summary(),
                    scope="M1 mechanics exercise with a declared local policy; uncalibrated SIM assumptions; not a prediction for this patient.")
        return {**data, "fingerprint": fingerprint(data)}
