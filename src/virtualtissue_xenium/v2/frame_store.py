"""Append committed snapshots without retaining the full trajectory in memory."""
from __future__ import annotations

import json
from pathlib import Path


class FrameStore:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.jsonl_path = self.directory / "frames.jsonl"
        self._jsonl = self.jsonl_path.open("w", encoding="utf8")
        self._array = (self.directory / "frames.json").open("w", encoding="utf8")
        self._array.write("[\n")
        self.count = 0

    def append(self, frame):
        value = json.dumps(frame, allow_nan=False)
        self._jsonl.write(value + "\n")
        self._jsonl.flush()
        if self.count:
            self._array.write(",\n")
        self._array.write(value)
        self._array.flush()
        self.count += 1

    def close(self):
        if not self._array.closed:
            self._array.write("\n]\n")
            self._array.close()
            self._jsonl.close()


def iter_frames(path):
    with Path(path).open(encoding="utf8") as f:
        for line in f:
            if line.strip():
                yield json.loads(line)
