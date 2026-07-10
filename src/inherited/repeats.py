from __future__ import annotations

from pathlib import Path


class RepeatIntervalFilter:
    """Stream repeat intervals and test whether a position falls inside one.

    Interval files are whitespace-separated rows ``chrom start end`` with
    half-open intervals ``[start, end)``.
    """

    def __init__(self, path: Path) -> None:
        self._handle = path.open(encoding="utf-8")
        self._start = 0
        self._end = 0
        self._exhausted = False
        self._load_next()

    def close(self) -> None:
        self._handle.close()

    def _load_next(self) -> None:
        for line in self._handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            self._start = int(parts[1])
            self._end = int(parts[2])
            return
        self._exhausted = True
        self._start = self._end = 2**63

    def advance_past(self, pos: int) -> None:
        """Drop intervals that end at or before ``pos``."""
        while not self._exhausted and pos >= self._end:
            self._load_next()

    def in_repeat(self, pos: int) -> bool:
        self.advance_past(pos)
        if self._exhausted:
            return False
        return self._start <= pos < self._end
