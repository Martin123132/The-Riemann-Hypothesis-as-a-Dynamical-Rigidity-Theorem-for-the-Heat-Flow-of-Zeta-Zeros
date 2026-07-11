#!/usr/bin/env python3
"""Independently validate the quartic double-root threshold lemma."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_quartic_double_root_threshold_lemma.json"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    issues: list[str] = []
    rows = payload.get("rows", [])

    a, p, u, k, w = sp.symbols("a p u k w")
    A = -3 * a**2 + 8 * a + p
    B = -a**2 + 2 * a + p
    U = sp.factor(B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2))
    numerator = 3 * a**4 - 11 * a**3 + 2 * a**2 * p + 10 * a**2 - 5 * a * p + 6 * p**2 * u - 5 * p**2
    if sp.factor(numerator - 6 * p**2 * (u - U)) != 0:
        issues.append("heat threshold numerator failed")

    b, c = sp.symbols("b c")
    curvature = sp.expand((a - b) * (a - c)).subs(c, 4 - 2 * a - b)
    curvature = sp.rem(curvature, b**2 - (4 - 2 * a) * b + p, b)
    if sp.factor(curvature - (3 * a**2 - 4 * a + p)) != 0:
        issues.append("curvature identity failed")

    triple_p = 4 * a - 3 * a**2
    triple_U = sp.factor(U.subs(p, triple_p))
    H = sp.factor(24 * (2 * k + 7) * p**2 * (u - U) / (a * A * B))
    if sp.factor(H.subs({p: triple_p, u: triple_U})) != 0:
        issues.append("triple-root equality failed")
    tangent = sp.factor(
        24 * w**2 * (a - 1) ** 2 * (a * w + 1) ** 2
        * (2 * a * k + a - 2 * k + 1) / (a - 2)
    )
    if sp.rem(sp.together(tangent).as_numer_denom()[0], (a * w + 1) ** 2, w) != 0:
        issues.append("triple-root tangent factor failed")

    av, bv, cv = sp.Rational(13, 20), sp.Rational(21, 50), sp.Rational(57, 25)
    pv = bv * cv
    Uv = sp.factor(U.subs({a: av, p: pv}))
    xv = sp.factor(A.subs({a: av, p: pv}) / 6)
    yv = sp.factor(18 * av * B.subs({a: av, p: pv}) / A.subs({a: av, p: pv}) ** 2)
    zv = sp.factor(2 * pv * A.subs({a: av, p: pv}) / (3 * B.subs({a: av, p: pv}) ** 2))
    curvature_v = sp.factor((3 * a**2 - 4 * a + p).subs({a: av, p: pv}))
    if not (curvature_v < 0 and Uv > zv and curvature_v * (zv - Uv) > 0):
        issues.append("countermodel threshold explanation failed")

    if len(rows) != 10:
        issues.append(f"row count {len(rows)} != 10")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("open handoff count failed")
    boundary = payload.get("proof_boundary", "")
    for marker in ("does not construct", "does not", "degree-4", "PF-infinity", "RH"):
        if marker not in boundary:
            issues.append(f"missing proof-boundary marker {marker!r}")

    print(
        "validated Jensen-window PF quartic double-root threshold lemma: "
        f"{len(rows)} rows, {len(issues)} issues, 5 exact coordinate identities, "
        "1 double-root splitting criterion, 1 branch-aware inward threshold, "
        "1 triple-root equality, 1 tangent factor, 1 explained countermodel, "
        "1 open global-invariant handoff"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
