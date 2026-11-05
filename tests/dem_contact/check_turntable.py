#!/usr/bin/env python3
"""Assert the rolling-sphere-on-a-turntable orbit from test_cu_02.

test_cu_02 builds two spinning bodies: a turntable plane with w=(0,0,omega)
and a sphere launched with w=(omega*x0/R, -v/R, 0), chosen so the contact
starts without sliding. The driver prints, per output step, the measured
centre (x, y) and the closed-form orbit

    xt = x0 + vx/wr*sin(wr*t)
    yt = vx/wr - (vx/wr)*cos(wr*t),   wr = (2/7)*omega

which is the classic rolling-sphere result and exists only because the
sphere's spin couples into the contact through friction. This script turns
that comparison into a pass/fail gate so the rotational integration path has
an analytic check rather than a printed table nobody reads.

A tolerance of 1e-3 m separates the two regimes by more than an order of
magnitude at the default window: with rotation enabled the deviation is
~4e-5 m, and with rotation disabled it is ~3.5e-2 m.
"""
import argparse
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trajectory', type=Path)
    parser.add_argument('--tolerance', type=float, default=1e-3)
    parser.add_argument('--min-rows', type=int, default=1)
    args = parser.parse_args()

    rows = []
    for line in args.trajectory.read_text().splitlines()[1:]:
        fields = line.split()
        if not fields:
            continue
        t, x, y, xt, yt = (float(value) for value in fields[:5])
        rows.append((t, x, y, xt, yt))

    if len(rows) < args.min_rows:
        raise SystemExit('TURNTABLE_ROLLING_FAIL too few rows: %d < %d'
                         % (len(rows), args.min_rows))
    if not all(math.isfinite(value) for row in rows for value in row):
        raise SystemExit('TURNTABLE_ROLLING_FAIL non-finite trajectory')

    deviations = [(t, math.hypot(x-xt, y-yt)) for t, x, y, xt, yt in rows]
    worst_time, worst = max(deviations, key=lambda item: item[1])
    if worst > args.tolerance:
        raise SystemExit('TURNTABLE_ROLLING_FAIL max deviation %.6e m at t=%.4f '
                         'exceeds %.1e m; the analytic orbit requires the sphere '
                         'spin to couple through friction, so this usually means '
                         'the rotational path is not running'
                         % (worst, worst_time, args.tolerance))
    print('TURNTABLE_ROLLING_PASS rows=%d max_deviation=%.6e tolerance=%.1e'
          % (len(rows), worst, args.tolerance))


if __name__ == '__main__':
    main()
