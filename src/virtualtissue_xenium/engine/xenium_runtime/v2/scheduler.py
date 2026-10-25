"""Persistent wake reasons, independent program clocks and per-cell day caps."""
from collections import defaultdict


class WakeScheduler:
    def __init__(self, parameters):
        self.p = parameters
        self.pending = defaultdict(set)
        self.previous = {}
        self.due = {}
        self.count = defaultdict(int)

    def wake(self, g, reason):
        self.pending[g].add(reason)

    def observe(self, g, contacts, bins):
        old_contacts, old_bins = self.previous.get(g, (frozenset(), bins))
        if set(contacts) - set(old_contacts):
            self.wake(g, "new_relevant_contact")
        if bins != old_bins:
            self.wake(g, "field_threshold_crossing")
        self.previous[g] = (frozenset(contacts), bins)

    def eligible(self, g, time, busy=False):
        return bool(self.pending[g]) and not busy and time >= self.due.get(g, 0) and self.count[g, time//1440] < self.p["max_per_cell_day"]

    def consume(self, g, time, label, in_contact=False):
        reasons = sorted(self.pending.pop(g, set()))
        interval = self.p["cd8_contact_interval_min"] if label == "T_CD8_cytotoxic" and in_contact else self.p["interval_min"][label]
        self.due[g] = time + interval
        self.count[g, time//1440] += 1
        return reasons
