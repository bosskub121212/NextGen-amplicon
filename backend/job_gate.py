"""job_gate.py — the job queue: who may run now, in what order, and when.

- slots: how many jobs run at once (None = unlimited). The console's "Jobs at once".
- pause: seconds to wait after a job ends (finished, failed or cancelled) before the
  next one starts, so the machine has given its memory back.
- An exclusive job (shotgun + HUMAnN, ~25 GB RAM) starts only when nothing runs, and
  while it runs or waits its turn no later job starts.
- Waiting jobs keep their place, first come first served, unless move() changes the
  order. A job enters with a seq (its saved place); a job requeued after a backend
  restart therefore goes back where it was. A waiting job notices a cancel within
  `poll` seconds and gives up its turn.
"""
from __future__ import annotations

import math
import threading
import time
from typing import Callable

RUNNING = "Waiting — HUMAnN job running"
NEXT = "Waiting — HUMAnN job starts when the running job ends"


def queued_label(pos: int) -> str:
    return f"Queued — #{pos} in line"


def pause_label(pos: int, left: float) -> str:
    when = "under a minute" if left < 60 else f"{math.ceil(left / 60)} min"
    return f"Queued — #{pos} · starts in {when} (pause after the last job)"


class JobGate:
    def __init__(self, poll: float = 2.0, slots: int | None = None, pause: float = 0.0,
                 clock: Callable[[], float] = time.time):
        self._cv = threading.Condition()
        self._active: dict[str, bool] = {}               # job id → exclusive
        self._queue: list[list] = []                     # [job id, exclusive, seq], in run order
        self._poll = poll
        self._slots = slots
        self._pause = float(pause or 0)
        self._clock = clock
        self._pause_until = 0.0

    # ── settings ──────────────────────────────────────────────────────────────
    def configure(self, slots: int | None = None, pause: float | None = None) -> None:
        with self._cv:
            if slots is not None:
                self._slots = max(1, int(slots))
            if pause is not None:
                self._pause = max(0.0, float(pause))
            self._cv.notify_all()

    def pause_now(self) -> None:
        """Start a pause now (after a backend restart: an interrupted job may still be
        giving its memory back)."""
        with self._cv:
            if self._pause > 0:
                self._pause_until = self._clock() + self._pause
            self._cv.notify_all()

    # ── queue ─────────────────────────────────────────────────────────────────
    def _pos(self, job_id: str) -> int:
        for i, q in enumerate(self._queue):
            if q[0] == job_id:
                return i
        return -1

    def _blocker(self, job_id: str, exclusive: bool) -> str | None:
        """None when the job may start now, else the label to show."""
        i = self._pos(job_id)
        pos = i + 1
        ahead = self._queue[:i]
        if any(self._active.values()):
            return RUNNING
        if exclusive:
            if self._active or ahead:
                return NEXT
        elif any(q[1] for q in ahead):
            return NEXT
        if self._slots is not None:
            free = self._slots - len(self._active)
            if len(ahead) >= free:
                return queued_label(pos)
        left = self._pause_until - self._clock()
        if left > 0:
            return pause_label(pos, left)
        return None

    def enter(self, job_id: str, exclusive: bool, is_cancelled: Callable[[], bool],
              on_wait: Callable[[str], None], seq: float | None = None) -> bool:
        """Block until the job may run. False when it was cancelled while waiting."""
        with self._cv:
            seq = time.time() if seq is None else float(seq)
            k = len(self._queue)
            while k > 0 and self._queue[k - 1][2] > seq:
                k -= 1
            self._queue.insert(k, [job_id, exclusive, seq])
            try:
                while True:
                    if is_cancelled():
                        return False
                    why = self._blocker(job_id, exclusive)
                    if why is None:
                        self._active[job_id] = exclusive
                        return True
                    try:                     # every poll: the caller may add detail that changes
                        on_wait(why)         # (e.g. the running job paused at a checkpoint)
                    except Exception:
                        pass
                    self._cv.wait(self._poll)
            finally:
                self._queue = [q for q in self._queue if q[0] != job_id]
                self._cv.notify_all()

    def leave(self, job_id: str) -> None:
        with self._cv:
            if self._active.pop(job_id, None) is not None and self._pause > 0:
                self._pause_until = self._clock() + self._pause
            self._cv.notify_all()

    def move(self, job_id: str, to: int) -> bool:
        """Put a waiting job at place `to` (1 = next to run). False if it is not waiting."""
        with self._cv:
            i = self._pos(job_id)
            if i < 0:
                return False
            q = self._queue.pop(i)
            to = max(1, min(int(to), len(self._queue) + 1))
            self._queue.insert(to - 1, q)
            for n, e in enumerate(self._queue, 1):       # the new order is the saved order
                e[2] = float(n)
            self._cv.notify_all()
            return True

    # ── reading ───────────────────────────────────────────────────────────────
    def order(self) -> list[str]:
        with self._cv:
            return [q[0] for q in self._queue]

    def seqs(self) -> dict:
        with self._cv:
            return {q[0]: q[2] for q in self._queue}

    def state(self) -> dict:
        with self._cv:
            return {"active": dict(self._active), "waiting": [q[0] for q in self._queue]}

    def snapshot(self) -> dict:
        with self._cv:
            left = max(0.0, self._pause_until - self._clock())
            return {"active": list(self._active), "waiting": [q[0] for q in self._queue],
                    "exclusive": [q[0] for q in self._queue if q[1]], "slots": self._slots,
                    "pause": self._pause, "pause_left": round(left, 1)}
