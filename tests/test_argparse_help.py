#!/usr/bin/env python3
"""Python 3.14 argparse rejects a help string with a bare '%' at add_argument time
("badly formed help string") — on nextapp that crashed blast_taxonomy.py before it
did anything, and every 12S/CytB run silently fell back to DADA2 names. Scan every
backend script for such strings."""
import ast, re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
FAILS = []
def main():
    bad = []
    for p in (ROOT / "backend").rglob("*.py"):
        try:
            t = ast.parse(p.read_text(encoding="utf-8", errors="replace"))
        except SyntaxError:
            continue
        for n in ast.walk(t):
            if isinstance(n, ast.Call) and getattr(n.func, "attr", "") in ("add_argument", "ArgumentParser", "add_parser"):
                for k in n.keywords:
                    if k.arg in ("help", "description", "epilog") and isinstance(k.value, ast.Constant) \
                            and isinstance(k.value.value, str) \
                            and "%" in re.sub(r"%%|%\([a-z_]+\)[sd]", "", k.value.value):
                        bad.append(f"{p.relative_to(ROOT)}:{n.lineno}")
    if bad:
        print("  FAIL  argparse help with a bare % (Python 3.14 crash): " + ", ".join(bad)); FAILS.append("argparse %")
    else:
        print("  ok    no argparse help string with a bare % (Python 3.14)")
    return 1 if FAILS else 0
if __name__ == "__main__":
    sys.exit(main())
