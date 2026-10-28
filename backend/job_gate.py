"""job_gate.py — lets one heavy job (shotgun + HUMAnN, ~25 GB RAM) run alone.

Normal jobs share the worker slots as before. An exclusive job starts only when no
other job is running, and while it runs (or waits its turn) no later job starts.
Waiting jobs keep their place: first come, first served. A waiting job notices a
cancel within `poll` seconds and gives up its turn.
"""
from __future__ import annotations

import threading
from typing import Callable

RUNNING = "Waiting — HUMAnN job running"
NEXT = "Waiting — HUMAnN job starts when the running job ends"


class JobGate:
    def __init__(self, poll: float = 2.0):
        self._cv = threading.Condition()
        self._active: dict[str, bool] = {}        # job id → exclusive
        self._queue: list[tuple[str, bool]] = []  # waiting, in arrival order
        self._poll = poll

    def _blocker(self, job_id: str, exclusive: bool) -> str | None:
        """None when the job may start now, else the label to show."""
        if any(self._active.values()):
            return RUNNING
        ahead = []
        for jid, ex in self._queue:
            if jid == job_id:
                break
            ahead.append(ex)
        if exclusive:
            if self._active or ahead:
                return NEXT
            return None
        return NEXT if any(ahead) else None

    def enter(self, job_id: str, exclusive: bool, is_cancelled: Callable[[], bool],
              on_wait: Callable[[str], None]) -> bool:
        """Block until the job may run. False when it was cancelled while waiting."""
        with self._cv:
            self._queue.append((job_id, exclusive))
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
            self._active.pop(job_id, None)
            self._cv.notify_all()

    def state(self) -> dict:
        with self._cv:
            return {"active": dict(self._active), "waiting": [q[0] for q in self._queue]}
