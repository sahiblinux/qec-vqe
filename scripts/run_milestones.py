#!/usr/bin/env python3
"""End-to-end milestone driver.

Regenerates every benchmark result (M1 noiseless VQE -> M4 QEC/KPI),
writes publication-quality figures and a machine-readable summary JSON
under ``out/``, then prints the KPI report.

Usage:
    python scripts/run_milestones.py            # full run
    python scripts/run_milestones.py --m1       # noiseless baselines only
    python scripts/run_milestones.py --mitigation
    python scripts/run_milestones.py --qec
    python scripts/run_milestones.py --skip mitigation --skip qec
"""

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec_vqe.benchmark import kpi_summary, run_all  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--m1", action="store_true", help="only noiseless baselines")
    ap.add_argument("--mitigation", action="store_true",
                    help="only noise + mitigation benchmarks")
    ap.add_argument("--qec", action="store_true", help="only QEC benchmarks")
    ap.add_argument("--skip", action="append", default=[],
                    help="skip a stage: m1|mitigation|qec")
    ap.add_argument("--out", type=str, default=str(ROOT / "out"))
    args = ap.parse_args(argv)

    skip = list(args.skip)
    if args.m1:
        skip += ["mitigation", "qec"]
    if args.mitigation:
        skip += ["m1", "qec"]
    if args.qec:
        skip += ["m1", "mitigation"]
    skip = sorted(set(skip))

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print(f"[driver] running milestones into {out} "
          f"(skip={skip or 'none'}) ...", flush=True)
    results = run_all(out_dir=out, skip=skip)
    print(kpi_summary(results))
    print(f"[driver] done in {time.time() - t0:.1f}s; "
          f"figures + out/results.json written.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
