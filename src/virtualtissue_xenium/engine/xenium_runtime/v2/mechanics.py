"""Mutable v2 tissue: sites, radii, target areas, fates and overdamped mechanics.

The immutable atlas is never edited. Bodies are ROI cells (they may act),
halo cells and unannotated occluded objects (mechanically active context that
never acts), and cells born or recruited during a run (they may act). The
outermost halo ring is a fixed boundary. Removed bodies stay in the arrays
with ``present = False`` and a recorded fate, so identities are never reused.
"""
import math
from collections import Counter

import numpy as np
from shapely.geometry import shape

from .geometry import power_geometry

TUMOUR = frozenset({"Tumor", "Tumor_proliferating"})
ACTING_ORIGINS = frozenset({"roi", "born", "recruited"})


class InvariantError(AssertionError):
    pass


def primary_label(cell):
    return cell["identity"]["candidates"][0]["label"] if cell.get("identity") else "occluded"


class Tissue:
    def __init__(self, ids, labels, origins, xy, measured_area, polygons, centre, roi_radius, fixed_radius, contract):
        self.contract = contract
        self.g = contract["geometry"]
        self.m = contract["mechanics"]
        n = len(ids)
        self.ids = np.asarray(ids, dtype=np.int64)
        self.labels = list(labels)
        self.origins = list(origins)
        self.parents = [None] * n
        self.xy = np.asarray(xy, dtype=float).copy()
        self.xy0 = self.xy.copy()
        self.measured_area = np.asarray(measured_area, dtype=float)
        self.polygons = polygons
        self.centre = np.asarray(centre, dtype=float)
        self.roi_radius = float(roi_radius)
        self.fixed_radius = float(fixed_radius)
        dist = np.hypot(*(self.xy - self.centre).T)
        self.fixed = np.array([o not in ACTING_ORIGINS for o in self.origins]) & (dist >= self.fixed_radius)
        self.present = np.ones(n, bool)
        self.acts = np.array([o in ACTING_ORIGINS for o in self.origins])
        self.fates = ["alive"] * n
        self.born_at = np.zeros(n)
        self.path_um = np.zeros(n)
        self.swaps = np.zeros(n, dtype=np.int64)
        self.n_initial_roi = int(sum(o == "roi" for o in self.origins))
        self.n_context = n - self.n_initial_roi
        self.born = 0
        self.recruited = 0
        self.next_id = int(self.ids.max()) + 1_000_000
        self.kappa = np.full(n, 1 / math.pi)
        self.fit = self._fit_weights()
        geometry = self.geometry()[1]
        # Target area := fitted t=0 area, so the measured configuration is a
        # declared mechanical equilibrium (zero pressure, zero force). kappa is
        # rescaled so every radius, and therefore the t=0 diagram, is unchanged.
        area0 = np.maximum(geometry.area, 1e-6)
        self.kappa = self.kappa * self.measured_area / area0
        self.A0 = area0
        self.A0_base = self.A0.copy()
        self.kappa0 = self.kappa.copy()
        self.in_roi = dist <= self.roi_radius
        self.rescues = 0
        self.anchor = self.xy.copy()
        anchoring = self.m["anchoring"]
        self.anchor_k = np.array([anchoring["stiffness_per_min"][anchoring["class_of_label"][label]] for label in self.labels])
        self.separations = 0
        self.confinements = 0
        self.compression = self.m.get("compression_constraint")
        self.constraint_rejections = 0
        self.constraint_trials = 0
        self.local_freezes = 0
        self.global_rollbacks = 0
        self.growth_blocked = 0
        self.insertion_blocked = 0
        # Pairs already closer than the hard core in the measured t=0 geometry
        # are admissible at their t=0 distance, so the null control stays static.
        idx, geometry = self.geometry()
        a, b = idx[geometry.pair_i], idx[geometry.pair_j]
        self._key_base = 1 << 32
        keys = np.minimum(a, b) * self._key_base + np.maximum(a, b)
        order = np.argsort(keys)
        self.core_keys0 = keys[order]
        self.core_d0 = geometry.pair_distance[order]

    @classmethod
    def from_atlas(cls, atlas, contract):
        cells = atlas["cells"]
        polygons = [shape(c["geometry"]) for c in cells]
        xy = np.array([[p.centroid.x, p.centroid.y] for p in polygons])
        centre = np.array(atlas["geometry"]["center_um"], dtype=float)
        roles = [c["role"] for c in cells]
        roi = np.array([r == "roi" for r in roles])
        dist = np.hypot(*(xy - centre).T)
        roi_radius = float(dist[roi].max()) + contract["mechanics"]["roi_radius_margin_um"]
        fixed_radius = roi_radius + float(atlas["geometry"]["halo_um"]) - contract["mechanics"]["fixed_ring_um"]
        tissue = cls([c["cell_id"] for c in cells], [primary_label(c) for c in cells], roles, xy,
                     [c["area_um2"] for c in cells], polygons, centre, roi_radius, fixed_radius, contract)
        measured = {tuple(sorted((c["cell_id"], e["cell_id"]))) for c in cells if c["role"] == "roi"
                    for e in c["contacts"] if e["gap_um"] <= 1}
        tissue.set_contact_memory(measured)
        return tissue

    # ----------------------------------------------------------- contact graph
    def set_contact_memory(self, measured_pairs):
        """Remember the t=0 graph so contacts hand over smoothly from measurement to simulation.

        At t=0 contact is the measured polygon graph (H04) and nothing else. A
        t=0 pair keeps its measured status (contact or no contact) while the
        distance between its two sites stays within ``contact_memory_um`` of
        its t=0 value; once either cell has moved, the simulated geometry
        decides. Only pairs with at least one acting body are tracked.
        """
        pos = {int(cid): g for g, cid in enumerate(self.ids.tolist())}
        idx, geometry = self.geometry()
        i, j, _ = geometry.contacts(self.g["contact_max_gap_um"])
        simulated = {tuple(sorted((int(idx[a]), int(idx[b])))) for a, b in zip(i.tolist(), j.tolist())}
        measured = {tuple(sorted((pos[a], pos[b]))) for a, b in measured_pairs if a in pos and b in pos}
        pairs = sorted((p for p in measured | simulated if self.acts[p[0]] or self.acts[p[1]]))
        self.memory = np.array(pairs, dtype=np.int64).reshape(-1, 2)
        self.memory_measured = np.array([p in measured for p in pairs], dtype=bool)
        self.memory_d0 = np.hypot(*(self.xy[self.memory[:, 0]] - self.xy[self.memory[:, 1]]).T) if len(pairs) else np.zeros(0)

    def contact_pairs(self, idx, geometry):
        """Current contact graph as {(g_a, g_b): provenance} over global indices (a < b)."""
        i, j, _ = geometry.contacts(self.g["contact_max_gap_um"])
        graph = {tuple(sorted((int(idx[a]), int(idx[b])))): "simulated_contact" for a, b in zip(i.tolist(), j.tolist())}
        if not hasattr(self, "memory") or not len(self.memory):
            return graph
        a, b = self.memory[:, 0], self.memory[:, 1]
        alive = self.present[a] & self.present[b]
        d = np.hypot(*(self.xy[a] - self.xy[b]).T)
        held = alive & (np.abs(d - self.memory_d0) <= self.g["contact_memory_um"])
        for (pa, pb), measured in zip(self.memory[held].tolist(), self.memory_measured[held].tolist()):
            if measured:
                graph[(pa, pb)] = "measured_contact_persisting"
            else:
                graph.pop((pa, pb), None)
        return graph

    # ------------------------------------------------------------------ geometry
    @property
    def rho(self):
        return np.sqrt(self.kappa * self._area_for_radius())

    def _area_for_radius(self):
        return getattr(self, "A0", self.measured_area)

    def geometry(self):
        idx = np.flatnonzero(self.present)
        return idx, power_geometry(self.xy[idx], self.rho[idx], self.g["dummy_sites"])

    def _fit_weights(self):
        lo, hi = (b / math.pi for b in self.g["kappa_bounds_times_pi"])
        target = self.measured_area
        for _ in range(self.g["fit_iterations"]):
            _, geometry = self.geometry()
            ratio = np.clip(target / np.maximum(geometry.area, 1e-6), .5, 2.)
            self.kappa = np.clip(self.kappa * ratio ** self.g["fit_exponent"], lo, hi)
        _, geometry = self.geometry()
        ratio = geometry.area / target
        fallback_iterations=0
        if self.g.get("recover_initial_hidden",False) and geometry.hidden.any():
            # Damped area matching at the unchanged measured sites. This is
            # initialization fitting, never a dynamic weight-inflation rescue.
            hi=32/math.pi
            for _ in range(800):
                update=np.clip(target/np.maximum(geometry.area,1e-6),.5,2.)
                self.kappa=np.clip(self.kappa*update**.15,lo,hi)
                _,geometry=self.geometry();fallback_iterations+=1
                ratio=geometry.area/target
                if not geometry.hidden.any() and np.max(abs(np.log(np.maximum(ratio,1e-12))))<.01:break
            if geometry.hidden.any():raise InvariantError("Initial geometry still hides measured cells after damped fitting")
        return dict(iterations=self.g["fit_iterations"],fallback_iterations=fallback_iterations, hidden=int(geometry.hidden.sum()),
                    area_ratio_quantiles=[round(float(q), 4) for q in np.quantile(ratio, [.05, .25, .5, .75, .95])],
                    kappa_at_bounds=int(((self.kappa <= lo * 1.0001) | (self.kappa >= hi * .9999)).sum()))

    def pressure(self, idx, geometry):
        A0 = self.A0[idx]
        return np.maximum(0., (A0 - geometry.area) / A0)

    def local_pressure(self, idx, geometry):
        """Mean one-sided pressure of each body and its simulated contacts."""
        P = self.pressure(idx, geometry)
        i, j, _ = geometry.contacts(self.g["contact_max_gap_um"])
        total = P + np.bincount(i, P[j], len(idx)) + np.bincount(j, P[i], len(idx))
        count = 1 + np.bincount(i, minlength=len(idx)) + np.bincount(j, minlength=len(idx))
        return total / count

    # ----------------------------------------------------------------- mechanics
    def relax(self, idx, geometry, dt, active=None):
        """One overdamped substep.

        v = mobility * F - k_anchor * (x - anchor) + active, with the pressure
        force F_i = -sum_j (P_i+P_j)/2 * L_ij * n_ij. ``active`` is an optional
        self-propulsion velocity (um/min) per present body (rows follow
        ``idx``); self-propelled bodies carry their anchor with them. Anchors
        creep toward the current site with the declared creep time, so a
        transient push is elastic and sustained pressure (growth) is plastic.
        A hard-core projection then keeps power-adjacent sites apart so no
        body is hidden by a heavier neighbour.
        """
        before = self.xy.copy() if self.compression else None
        old_path = self.path_um.copy() if self.compression else None
        old_anchor = self.anchor.copy() if self.compression else None
        P = self.pressure(idx, geometry)
        i, j, L, nrm = geometry.pair_i, geometry.pair_j, geometry.pair_length, geometry.pair_normal
        f = ((P[i] + P[j]) / 2 * L)[:, None] * nrm
        n = len(idx)
        F = np.c_[np.bincount(j, f[:, 0], n) - np.bincount(i, f[:, 0], n),
                  np.bincount(j, f[:, 1], n) - np.bincount(i, f[:, 1], n)]
        step = (self.m["mobility_per_min"] * F - self.anchor_k[idx, None] * (self.xy[idx] - self.anchor[idx])) * dt
        norm = np.hypot(*step.T)
        cap = self.m["max_step_um_per_substep"]
        capped = norm > cap
        step[capped] *= (cap / norm[capped])[:, None]
        if active is not None:
            step = step + active * dt
            norm = np.hypot(*step.T)
        step[self.fixed[idx]] = 0.
        self.xy[idx] += step
        self.path_um[idx] += np.hypot(*step.T)
        if active is not None:
            moving = np.hypot(*active.T) > 0
            self.anchor[idx[moving]] = self.xy[idx[moving]]
        creep = min(1., dt / self.m["anchoring"]["creep_min"])
        self.anchor[idx] += (self.xy[idx] - self.anchor[idx]) * creep
        self.separate(idx, geometry)
        self.confine()
        if self.compression:
            self.guard_motion(before,idx,geometry)
            actual = np.linalg.norm(self.xy-before,axis=1)
            self.path_um = old_path+actual
            self.anchor = old_anchor+(self.xy-old_anchor)*creep
            if active is not None:self.anchor[idx[moving]]=self.xy[idx[moving]]
            return dict(max_step=float(actual.max()),capped=int(capped.sum()))
        rescue = self.m["rescue"]
        weak = (geometry.hidden | (geometry.area < rescue["area_fraction_below"] * self.A0[idx])) & ~self.fixed[idx]
        if weak.any():
            g = idx[weak]
            self.kappa[g] = np.minimum(self.kappa[g] * rescue["kappa_factor_per_substep"],
                                       self.kappa0[g] * rescue["kappa_max_times_initial"])
            self.rescues += int(weak.sum())
        return dict(max_step=float(norm.max()) if n else 0., capped=int(capped.sum()))

    def floor_fractions(self,labels=None):
        p=self.compression
        return np.array([p["area_floor_fraction_by_label"].get(label,p["fallback_fraction"]) for label in (self.labels if labels is None else labels)])

    def valid_areas(self,geometry,floors):
        return geometry.area >= floors*(1-1e-8)

    def guard_motion(self,before,idx,previous):
        """Backtrack only the neighborhoods that violate the calibrated area floor.

        The area bound is a 2D structural proxy, not a fitted elastic modulus.
        A jam rejects displacement; it never deletes a cell or inflates weights.
        """
        floors=self.floor_fractions()[idx]*self.A0_base[idx]
        if not self.valid_areas(previous,floors).all():raise InvariantError("Compression floor already violated before motion")
        displacement=self.xy-before
        for attempt in range(self.compression["max_backtracks"]):
            self.constraint_trials+=1
            current=self.geometry()[1]
            bad=~self.valid_areas(current,floors)
            if not bad.any():return
            affected=bad.copy()
            for graph in (previous,current):
                linked=(bad[graph.pair_i]|bad[graph.pair_j]) & (graph.pair_gap<=self.g["contact_max_gap_um"])
                affected[graph.pair_i[linked]]=True;affected[graph.pair_j[linked]]=True
            displacement[idx[affected]]*=.5
            self.xy=before+displacement
            self.constraint_rejections+=1
        # A local jam must not freeze an independent part of the tissue.
        frozen=np.zeros(len(idx),bool)
        for _ in range(self.compression.get("max_freeze_passes",4)):
            current=self.geometry()[1];self.constraint_trials+=1
            bad=~self.valid_areas(current,floors)
            if not bad.any():return
            affected=bad.copy()
            for graph in (previous,current):
                linked=(bad[graph.pair_i]|bad[graph.pair_j]) & (graph.pair_gap<=self.g["contact_max_gap_um"])
                affected[graph.pair_i[linked]]=True;affected[graph.pair_j[linked]]=True
            frozen|=affected
            displacement[idx[frozen]]=0
            self.xy=before+displacement
            self.local_freezes+=int(affected.sum())
        if self.valid_areas(self.geometry()[1],floors).all():return
        # Fail closed only if even the connected-neighborhood rollback fails.
        self.xy=before
        self.global_rollbacks+=1
        self.constraint_rejections+=1

    def grow(self,g,target):
        if not self.compression:
            self.A0[g]=target;return True
        old=self.A0[g]
        self.A0[g]=target
        idx,geo=self.geometry()
        if self.valid_areas(geo,self.floor_fractions()[idx]*self.A0_base[idx]).all():return True
        self.A0[g]=old;self.growth_blocked+=1
        return False

    def grow_many(self,proposals):
        """Resolve simultaneous growth, withholding pressure-limited neighborhoods."""
        if not proposals:return set()
        old=self.A0.copy();idx=np.flatnonzero(self.present)
        changing=set(proposals)
        for g,area in proposals.items():self.A0[g]=area
        floors=self.floor_fractions()[idx]*self.A0_base[idx]
        for _ in range(self.compression["max_backtracks"]):
            geo=self.geometry()[1];bad=~self.valid_areas(geo,floors)
            if not bad.any():
                self.growth_blocked+=len(proposals)-len(changing)
                return changing
            affected=bad.copy();links=(bad[geo.pair_i]|bad[geo.pair_j]) & (geo.pair_gap<=self.g["contact_max_gap_um"])
            affected[geo.pair_i[links]]=True;affected[geo.pair_j[links]]=True
            blocked=changing.intersection(idx[affected].tolist())
            if not blocked:break
            for g in blocked:self.A0[g]=old[g]
            changing.difference_update(blocked)
        self.A0=old;self.growth_blocked+=len(proposals)
        return set()

    def can_insert(self,label,xy,area,kappa):
        if not self.compression:return True
        idx=np.flatnonzero(self.present)
        floor=self.floor_fractions()[idx]*self.A0_base[idx]
        new_floor=self.floor_fractions([label])[0]*area
        geo=power_geometry(np.vstack([self.xy[idx],xy]),np.r_[self.rho[idx],math.sqrt(area*kappa)],self.g["dummy_sites"])
        valid=self.valid_areas(geo,np.r_[floor,new_floor]).all()
        if not valid:self.insertion_blocked+=1
        return bool(valid)

    def can_divide(self,g,axis):
        if not self.compression:return True
        if self.A0[g]<2*self.A0_base[g]*(1-1e-8):return False
        idx=np.flatnonzero(self.present);j=int(np.flatnonzero(idx==g)[0])
        half=self.A0[g]/2
        axis=np.asarray(axis)/np.linalg.norm(axis)
        offset=self.contract["division"]["daughter_offset_fraction_of_radius"]*math.sqrt(self.kappa[g]*half)
        xy=self.xy[idx].copy();xy[j]-=offset*axis
        rho=self.rho[idx].copy();rho[j]=math.sqrt(self.kappa[g]*half)
        floors=self.floor_fractions()[idx]*self.A0_base[idx]
        floors[j]=self.floor_fractions([self.labels[g]])[0]*half
        geo=power_geometry(np.vstack([xy,self.xy[g]+offset*axis]),np.r_[rho,rho[j]],self.g["dummy_sites"])
        valid=self.valid_areas(geo,np.r_[floors,floors[j]]).all()
        if not valid:self.insertion_blocked+=1
        return bool(valid)

    def confine(self):
        """Mobile sites stay inside the wall at the inner edge of the fixed ring."""
        mobile = self.present & ~self.fixed
        vec = self.xy[mobile] - self.centre
        d = np.hypot(*vec.T)
        out = d > self.fixed_radius
        if out.any():
            rows = np.flatnonzero(mobile)[out]
            self.xy[rows] = self.centre + vec[out] * (self.fixed_radius / d[out])[:, None]
            self.confinements += int(out.sum())

    def separate(self, idx, geometry):
        """One Jacobi projection pushing apart power-adjacent sites closer than the hard core."""
        core = self.m["hard_core"]
        i, j = idx[geometry.pair_i], idx[geometry.pair_j]
        rho = self.rho
        w = rho * rho
        dmin = np.maximum(core["fraction_of_radii"] * (rho[i] + rho[j]), np.sqrt(np.abs(w[i] - w[j])) + core["weight_margin_um"])
        if len(self.core_keys0):
            keys = np.minimum(i, j) * self._key_base + np.maximum(i, j)
            pos = np.minimum(np.searchsorted(self.core_keys0, keys), len(self.core_keys0) - 1)
            known = self.core_keys0[pos] == keys
            dmin = np.where(known, np.minimum(dmin, self.core_d0[pos] * (1 - 1e-9)), dmin)
        vec = self.xy[j] - self.xy[i]
        d = np.hypot(*vec.T)
        close = d < dmin
        if not close.any():
            return
        i, j, vec, d, dmin = i[close], j[close], vec[close], d[close], dmin[close]
        unit = vec / np.maximum(d, 1e-9)[:, None]
        share_i = np.where(self.fixed[i], 0., np.where(self.fixed[j], 1., .5))
        share_j = np.where(self.fixed[j], 0., 1. - share_i)
        push = (dmin - d)[:, None] * unit
        n = len(self.ids)
        delta = np.zeros((n, 2))
        for axis in (0, 1):
            delta[:, axis] = np.bincount(j, push[:, axis] * share_j, n) - np.bincount(i, push[:, axis] * share_i, n)
        self.xy += delta
        self.separations += int(close.sum())

    def swap(self, a, b):
        distance = float(np.hypot(*(self.xy[a] - self.xy[b])))
        self.xy[[a, b]] = self.xy[[b, a]]
        self.anchor[[a, b]] = self.xy[[a, b]]
        self.path_um[a] += distance
        self.path_um[b] += distance
        self.swaps[a] += 1
        self.swaps[b] += 1
        return distance

    # ------------------------------------------------------ birth and removal
    def add_body(self, label, origin, xy, A0, kappa, parent, time_min):
        cid = self.next_id
        self.next_id += 1
        self.ids = np.r_[self.ids, cid]
        self.labels.append(label)
        self.origins.append(origin)
        self.parents.append(parent)
        self.xy = np.vstack([self.xy, xy])
        self.xy0 = np.vstack([self.xy0, xy])
        self.measured_area = np.r_[self.measured_area, np.nan]
        self.polygons.append(None)
        self.fixed = np.r_[self.fixed, False]
        self.present = np.r_[self.present, True]
        self.acts = np.r_[self.acts, origin in ACTING_ORIGINS]
        self.fates.append("alive")
        self.born_at = np.r_[self.born_at, time_min]
        self.path_um = np.r_[self.path_um, 0.]
        self.swaps = np.r_[self.swaps, 0]
        self.kappa = np.r_[self.kappa, kappa]
        self.kappa0 = np.r_[self.kappa0, kappa]
        self.A0 = np.r_[self.A0, A0]
        self.A0_base = np.r_[self.A0_base, A0]
        self.in_roi = np.r_[self.in_roi, np.hypot(*(np.asarray(xy) - self.centre)) <= self.roi_radius]
        self.anchor = np.vstack([self.anchor, xy])
        anchoring = self.m["anchoring"]
        self.anchor_k = np.r_[self.anchor_k, anchoring["stiffness_per_min"][anchoring["class_of_label"][label]]]
        if origin == "born":
            self.born += 1
        elif origin == "recruited":
            self.recruited += 1
        return len(self.ids) - 1

    def divide(self, g, axis, time_min):
        """Split body g along ``axis``; g keeps its identity, one daughter is born."""
        axis = np.asarray(axis, dtype=float)
        axis = axis / np.hypot(*axis)
        half = self.A0[g] / 2
        offset = self.contract["division"]["daughter_offset_fraction_of_radius"] * math.sqrt(self.kappa[g] * half)
        centre = self.xy[g].copy()
        self.xy[g] = centre - offset * axis
        self.anchor[g] = self.xy[g]
        self.A0[g] = self.A0_base[g] = half
        return self.add_body(self.labels[g], "born", centre + offset * axis, half, self.kappa[g], int(self.ids[g]), time_min)

    def remove(self, g, fate):
        if not self.present[g]:
            raise InvariantError("Removing an absent body")
        if self.fixed[g]:
            raise InvariantError("Fixed boundary bodies cannot be removed")
        self.present[g] = False
        self.fates[g] = fate

    def overlap_audit(self, quad_segs=32):
        """Exhaustive polygon check (slow): summed cell areas versus their union."""
        from shapely import union_all
        idx, geometry = self.geometry()
        polygons = [geometry.polygon(k, quad_segs) for k in range(len(idx))]
        total = float(sum(p.area for p in polygons))
        union = float(union_all(polygons).area)
        return dict(cells=len(idx), summed_area_um2=round(total, 3), union_area_um2=round(union, 3),
                    overlap_um2=round(total - union, 6), relative_overlap=(total - union) / total if total else 0.)

    # ------------------------------------------------------------ accounting
    def update_roi(self):
        """Return bodies of acting lineage that crossed the ROI boundary."""
        now = np.hypot(*(self.xy - self.centre).T) <= self.roi_radius
        crossed = np.flatnonzero((now != self.in_roi) & self.present & self.acts)
        self.in_roi = now
        return crossed

    def accounting(self):
        acting = Counter()
        for g in np.flatnonzero(self.acts):
            if self.present[g]:
                acting["alive_in_roi" if self.in_roi[g] else "alive_outside_roi"] += 1
            else:
                acting[self.fates[g]] += 1
        context = Counter("present" if self.present[g] else self.fates[g] for g in np.flatnonzero(~self.acts))
        return dict(acting=dict(acting), acting_total=sum(acting.values()),
                    expected_acting_total=self.n_initial_roi + self.born + self.recruited,
                    context=dict(context), context_total=sum(context.values()), expected_context_total=self.n_context,
                    born=self.born, recruited=self.recruited)

    def check_invariants(self):
        acc = self.accounting()
        if acc["acting_total"] != acc["expected_acting_total"] or acc["context_total"] != acc["expected_context_total"]:
            raise InvariantError(f"Individual accounting broken: {acc}")
        if len(set(self.ids.tolist())) != len(self.ids):
            raise InvariantError("Duplicate body identity")
        if not (np.isfinite(self.xy).all() and np.isfinite(self.kappa).all() and (self.A0 > 0).all()):
            raise InvariantError("Nonfinite geometry state")
        if not np.array_equal(self.xy[self.fixed], self.xy0[self.fixed]):
            raise InvariantError("A fixed boundary body moved")
        if (self.acts & self.fixed).any():
            raise InvariantError("An acting body is fixed")
        return acc
