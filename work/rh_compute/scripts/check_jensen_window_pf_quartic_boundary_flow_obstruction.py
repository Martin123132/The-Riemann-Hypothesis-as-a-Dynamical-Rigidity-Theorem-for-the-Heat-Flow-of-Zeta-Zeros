#!/usr/bin/env python3
"""Independently validate the exact quartic boundary-flow obstruction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_quartic_boundary_flow_obstruction.json"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    issues: list[str] = []
    rows = payload.get("rows", [])

    w = sp.symbols("w")
    x, y, z, u, k = sp.symbols("x y z u k")
    quartic = 1 + 4 * w + 6 * x * w**2 + 4 * x**2 * y * w**3 + x**3 * y**2 * z * w**4
    discriminant = sp.factor(sp.discriminant(quartic, w))
    Q = sp.factor(discriminant / (256 * x**6 * y**2))
    F_yz = y**2 * z**2 - 6 * y * z + 4 * y + 4 * z - 3
    if sp.simplify(sp.discriminant(sp.diff(quartic, w), w) + 6912 * x**6 * y**2 * F_yz) != 0:
        issues.append("derivative cubic factorization failed")

    xp = 2 * ((2 * k - 1) * (1 - x) - (2 * k + 3) * x * (1 - y))
    yp = 2 * y * ((2 * k + 1) * (1 - y) - (2 * k + 5) * y * (1 - z))
    zp = 2 * y * z * ((2 * k + 3) * (1 - z) - (2 * k + 7) * z * (1 - u))
    Qp = sp.factor(sp.diff(Q, x) * xp + sp.diff(Q, y) * yp + sp.diff(Q, z) * zp)
    if sp.degree(Qp, u) != 1:
        issues.append("quartic heat derivative is not linear in u")

    a, b, c = sp.Rational(13, 20), sp.Rational(21, 50), sp.Rational(57, 25)
    e2 = a**2 + 2 * a * b + 2 * a * c + b * c
    e3 = a**2 * b + a**2 * c + 2 * a * b * c
    e4 = a**2 * b * c
    xv = sp.factor(e2 / 6)
    yv = sp.factor((e3 / 4) / xv**2)
    zv = sp.factor(e4 / (xv**3 * yv**2))
    uv = zv
    factorization = sp.expand((1 + a * w) ** 2 * (1 + b * w) * (1 + c * w))
    if sp.expand(quartic.subs({x: xv, y: yv, z: zv}) - factorization) != 0:
        issues.append("rational quartic factorization failed")

    F = lambda left, right: sp.factor(
        left**2 * right**2 - 6 * left * right + 4 * left + 4 * right - 3
    )
    cubic_values = [F(xv, yv), F(yv, zv), F(zv, uv)]
    if not all(value < 0 for value in cubic_values):
        issues.append("a cubic margin is not strict")
    ratio_margins = [
        xv - sp.Rational(1, 3),
        yv - sp.Rational(3, 5),
        zv - sp.Rational(5, 7),
        uv - sp.Rational(7, 9),
    ]
    if not all(value > 0 for value in ratio_margins):
        issues.append("a pointwise ratio margin is not positive")
    if not (xv < yv < zv == uv < 1):
        issues.append("contraction ordering failed")

    q_value = sp.factor(Q.subs({x: xv, y: yv, z: zv}))
    qp_value = sp.factor(Qp.subs({x: xv, y: yv, z: zv, u: uv, k: 1}))
    expected_qp = -sp.Rational(
        13108711376416987159336748097,
        20606742971316325673502124987495,
    )
    if q_value != 0 or qp_value != expected_qp or not qp_value < 0:
        issues.append("outward quartic derivative failed")

    if len(rows) != 10:
        issues.append(f"row count {len(rows)} != 10")
    if sum(row.get("role") == "forbidden_promotion" for row in rows) != 1:
        issues.append("promotion guard count failed")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("open handoff count failed")
    boundary = payload.get("proof_boundary", "")
    for marker in ("abstract and local", "does not show failure", "PF-infinity", "RH"):
        if marker not in boundary:
            issues.append(f"missing proof-boundary marker {marker!r}")

    print(
        "validated Jensen-window PF quartic boundary-flow obstruction: "
        f"{len(rows)} rows, {len(issues)} issues, 4 exact quartic identities, "
        "1 hyperbolic boundary point, 4 positive ratio margins, 3 strict cubic margins, "
        "1 negative quartic derivative, 1 blocked promotion, "
        "1 open coupled-invariant handoff"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
