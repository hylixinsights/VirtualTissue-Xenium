"""Reporting adapter around the unchanged, bundle-compatible biological kernel.

Only event metadata and recording metadata change. No clocks, choices, random
draws, state transitions or population-limit rules are overridden here.
"""
from xenium_runtime.io import fingerprint
from xenium_runtime.v2.autonomous import AutonomousSimulation
from jev_payload import PROMPT_VERSION

REPORTING_VERSION = "vtx-recording-2"


class RecordedSimulation(AutonomousSimulation):
    def block_reasons(self, g):
        cycle = self.bio[g].get("cycle") or {}
        if cycle.get("growth_blocked"):
            return ["compression"]
        if not cycle.get("blocked"):
            return []
        reasons = []
        if len(self.tissue.ids) >= self.p["max_population"]:
            reasons.append("population_ceiling")
        if self.lp[g] > self.div["division_blocked_local_pressure"]:
            reasons.append("local_pressure")
        return reasons or ["compression"]

    def event(self, kind, g=None, **fields):
        if g is not None:
            fields["position_um"] = self.tissue.xy[g].tolist()
            target = fields.get("target_id") or fields.get("daughter_id")
            if target is not None:
                matches = (self.tissue.ids == target).nonzero()[0]
                if len(matches):
                    fields["target_position_um"] = self.tissue.xy[int(matches[0])].tolist()
        if kind == "division_blocked_crowding":
            reasons = []
            if len(self.tissue.ids) >= self.p["max_population"]:
                reasons.append("population_ceiling")
            if self.lp[g] > self.div["division_blocked_local_pressure"]:
                reasons.append("local_pressure")
            fields.update(block_reasons=reasons, allocated_bodies=len(self.tissue.ids),
                          max_population=self.p["max_population"])
            kind = "division_blocked_" + reasons[0]
        elif kind in ("growth_blocked_compression", "division_blocked_compression"):
            fields["block_reasons"] = ["compression"]
        super().event(kind, g, **fields)

    def frame(self):
        super().frame()
        if not hasattr(self, "bio") or not self.frames:
            return
        indices = self.tissue.present.nonzero()[0]
        # Enrich the copied recording, never the live biological state.
        for g, state in zip(indices, self.frames[-1]["states"]):
            if state.get("cycle"):
                state["cycle"]["block_reasons"] = self.block_reasons(int(g))

    def episode(self):
        episode = super().episode()
        episode.update(reporting_version=REPORTING_VERSION, prompt_version=PROMPT_VERSION,
                       prompt_context=getattr(self, "prompt_context", None),
                       sample_name=self.atlas.get("provenance", {}).get("sample_id") or self.atlas["patch_id"],
                       population_limit_semantics="All allocated bodies, including context and removed bodies; computational limit.")
        episode.pop("fingerprint")
        episode["fingerprint"] = fingerprint(episode)
        return episode
