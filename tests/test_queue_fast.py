#!/usr/bin/env python3
"""Fast checks for the job queue: one job at a time, a pause between jobs, queue order
that can be changed, waiting jobs kept across a backend restart.

Run alone (python3 tests/test_queue_fast.py) or from tests/test_fast.py.
"""
from __future__ import annotations

import sys
import threading
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))

FAILS = []
TESTS = []


def check(name, cond, detail=""):
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILS.append(name)


def _until(cond, secs=5.0):
    end = time.time() + secs
    while time.time() < end:
        if cond():
            return True
        time.sleep(0.02)
    return False


class Clock:
    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t


def _job(G, jid, ev, lab, excl=False, hold=None, cancel=None, seq=None):
    def body():
        ok = G.enter(jid, excl, (lambda: bool(cancel and cancel.is_set())),
                     lambda l: lab.setdefault(jid, []).append(l), seq=seq)
        ev.append(("in" if ok else "cancelled", jid))
        if ok:
            if hold:
                hold.wait(5)
            ev.append(("out", jid))
            G.leave(jid)
    t = threading.Thread(target=body, daemon=True)
    t.start()
    return t


def t_one_at_a_time(tmp):
    import job_gate as JG
    G, ev, lab = JG.JobGate(poll=0.03, slots=1, pause=0), [], {}
    ha, hb = threading.Event(), threading.Event()
    _job(G, "a", ev, lab, hold=ha)
    _until(lambda: ("in", "a") in ev)
    _job(G, "b", ev, lab, hold=hb)
    _job(G, "c", ev, lab)
    _until(lambda: "b" in lab and "c" in lab)
    check("queue: one job at a time", ("in", "b") not in ev and ("in", "c") not in ev, str(ev))
    check("queue: position shown", lab["b"][-1] == "Queued — #1 in line" and lab["c"][-1] == "Queued — #2 in line",
          str(lab))
    ha.set()
    check("queue: next starts when the running job ends", _until(lambda: ("in", "b") in ev), str(ev))
    check("queue: third still waits", ("in", "c") not in ev and _until(lambda: lab["c"][-1] == "Queued — #1 in line"),
          str(lab.get("c")))
    hb.set()
    check("queue: third runs last", _until(lambda: ("in", "c") in ev), str(ev))


def t_pause(tmp):
    import job_gate as JG
    clk = Clock()
    G, ev, lab = JG.JobGate(poll=0.03, slots=1, pause=180, clock=clk), [], {}
    ha, hb = threading.Event(), threading.Event()
    _job(G, "a", ev, lab, hold=ha)
    check("pause: none before the first job", _until(lambda: ("in", "a") in ev), str(ev))
    _job(G, "b", ev, lab, hold=hb)
    _until(lambda: "b" in lab)
    ha.set()
    _until(lambda: ("out", "a") in ev)
    check("pause: next job waits after the first ends", _until(lambda: "starts in 3 min" in lab["b"][-1])
          and ("in", "b") not in ev, str(lab["b"][-3:]))
    clk.t += 150
    check("pause: under a minute left shown", _until(lambda: "starts in under a minute" in lab["b"][-1])
          and ("in", "b") not in ev, str(lab["b"][-2:]))
    clk.t += 31
    check("pause: next job starts when the pause is over", _until(lambda: ("in", "b") in ev), str(ev))
    snap = G.snapshot()
    check("pause: snapshot", snap["slots"] == 1 and snap["pause"] == 180 and snap["active"] == ["b"], str(snap))
    hb.set()


def t_move(tmp):
    import job_gate as JG
    G, ev, lab = JG.JobGate(poll=0.03, slots=1, pause=0), [], {}
    ha = threading.Event()
    _job(G, "a", ev, lab, hold=ha)
    _until(lambda: ("in", "a") in ev)
    for j in ("b", "c", "d"):
        _job(G, j, ev, lab)
        _until(lambda: j in lab)
    check("move: arrival order", G.order() == ["b", "c", "d"], str(G.order()))
    check("move: urgent job to #2", G.move("d", 2) and G.order() == ["b", "d", "c"], str(G.order()))
    check("move: to the front", G.move("c", 1) and G.order() == ["c", "b", "d"], str(G.order()))
    check("move: position clamped", G.move("c", 99) and G.order() == ["b", "d", "c"], str(G.order()))
    check("move: running or unknown job refused", not G.move("a", 1) and not G.move("zz", 1))
    check("move: labels follow the new order", _until(lambda: lab["c"][-1] == "Queued — #3 in line"), str(lab["c"][-2:]))
    order = []
    ha.set()
    for _ in range(3):
        _until(lambda: len([e for e in ev if e[0] == "in"]) > len(order) + 1)
        order = [e[1] for e in ev if e[0] == "in"][1:]
    check("move: jobs run in the changed order", _until(lambda: [e[1] for e in ev if e[0] == "in"] == ["a", "b", "d", "c"]),
          str(ev))


def t_seq_restore(tmp):
    import job_gate as JG
    G, ev, lab = JG.JobGate(poll=0.03, slots=1, pause=0), [], {}
    ha = threading.Event()
    _job(G, "a", ev, lab, hold=ha)
    _until(lambda: ("in", "a") in ev)
    _job(G, "x", ev, lab, seq=3.0)
    _until(lambda: "x" in lab)
    _job(G, "y", ev, lab, seq=1.0)
    _until(lambda: "y" in lab)
    _job(G, "z", ev, lab, seq=2.0)
    _until(lambda: "z" in lab)
    check("restore: waiting jobs placed by their saved order", G.order() == ["y", "z", "x"], str(G.order()))
    G.move("x", 1)
    check("restore: seq() gives the order to save", [j for j, _ in sorted(G.seqs().items(), key=lambda kv: kv[1])]
          == ["x", "y", "z"], str(G.seqs()))
    ha.set()


def t_configure(tmp):
    import job_gate as JG
    G, ev, lab = JG.JobGate(poll=0.03, slots=1, pause=0), [], {}
    ha, hb, hc = threading.Event(), threading.Event(), threading.Event()
    _job(G, "a", ev, lab, hold=ha)
    _until(lambda: ("in", "a") in ev)
    _job(G, "b", ev, lab, hold=hb)
    _job(G, "c", ev, lab, hold=hc)
    _until(lambda: "b" in lab and "c" in lab)
    G.configure(slots=2)
    check("configure: one more slot lets exactly one more job in",
          _until(lambda: ("in", "b") in ev) and ("in", "c") not in ev, str(ev))
    G.configure(slots=1)
    ha.set()
    time.sleep(0.2)
    check("configure: fewer slots never stop running jobs, and no new job until under the limit",
          ("in", "c") not in ev, str(ev))
    hb.set()
    check("configure: next job after the running one ends", _until(lambda: ("in", "c") in ev), str(ev))
    hc.set()
    G2 = JG.JobGate(poll=0.03)
    check("defaults: unlimited slots, no pause (old behaviour)", G2.snapshot()["slots"] is None
          and G2.snapshot()["pause"] == 0, str(G2.snapshot()))


def t_pause_now(tmp):
    import job_gate as JG
    clk = Clock()
    G, ev, lab = JG.JobGate(poll=0.03, slots=1, pause=180, clock=clk), [], {}
    G.pause_now()
    _job(G, "a", ev, lab)
    check("restart: requeued jobs wait one pause", _until(lambda: "a" in lab and "starts in 3 min" in lab["a"][-1])
          and ("in", "a") not in ev, str(lab.get("a")))
    clk.t += 181
    check("restart: then start", _until(lambda: ("in", "a") in ev), str(ev))


def t_main_wiring(tmp):
    m = (ROOT / "backend" / "main.py").read_text()
    for k in ("GATE = JobGate(", "QUEUE_CFG", "_load_queue_cfg()", "_save_queue_cfg(", '@app.get("/queue")',
              '@app.post("/queue/move")', "pause_minutes", "def _requeue_waiting(", "GATE.pause_now()",
              "queue_seq", "Back in the queue after a server restart", "seq=seq", "starts after"):
        check(f"main.py: {k}", k in m)
    i = m.find("def load_jobs")
    check("main.py: a queued job is not marked failed on restart",
          '("running", "queued")' not in m[i:i + 1200], "")
    j = m.find("def set_config")
    check("main.py: changing slots does not rebuild the executor", "ThreadPoolExecutor(" not in m[j:j + 900], "")
    ui = (ROOT / "frontend-beta" / "index.html").read_text()
    for k in ('id="queueBtn"', 'name:"queue"', "function renderQueue(", "data-qmove", "To front", "Pause between jobs",
              '"/queue/move"', 'id="queueCount"'):
        check(f"console: {k}", k in ui)


TESTS += [t_one_at_a_time, t_pause, t_move, t_seq_restore, t_configure, t_pause_now, t_main_wiring]


def main() -> int:
    print("job queue")
    for t in TESTS:
        try:
            t(None)
        except Exception:
            check(t.__name__, False, traceback.format_exc()[-900:])
    return len(FAILS)


if __name__ == "__main__":
    n = main()
    print(f"\n{'ALL OK' if not n else f'{n} FAILED'}")
    sys.exit(1 if n else 0)
