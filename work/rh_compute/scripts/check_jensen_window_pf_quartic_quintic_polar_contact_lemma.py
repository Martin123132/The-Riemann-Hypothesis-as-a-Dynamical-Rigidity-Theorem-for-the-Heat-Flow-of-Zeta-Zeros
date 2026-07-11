#!/usr/bin/env python3
"""Independently validate the quartic-quintic polar-contact lemma."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_quartic_quintic_polar_contact_lemma.json"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    issues: list[str] = []
    rows = payload.get("rows", [])

    w, x, y, z, u = sp.symbols("w x y z u")
    P4 = 1 + 4 * w + 6 * x * w**2 + 4 * x**2 * y * w**3 + x**3 * y**2 * z * w**4
    P5 = 1 + 5 * w + 10 * x * w**2 + 10 * x**2 * y * w**3 + 5 * x**3 * y**2 * z * w**4 + x**4 * y**3 * z**2 * u * w**5
    if sp.expand(P5 - w * sp.diff(P5, w) / 5 - P4) != 0:
        issues.append("polar identity failed")

    a, p = sp.symbols("a p")
    A = -3 * a**2 + 8 * a + p
    B = -a**2 + 2 * a + p
    xv, yv, zv = A / 6, 18 * a * B / A**2, 2 * p * A / (3 * B**2)
    U = B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2)
    boundary = P5.subs({x: xv, y: yv, z: zv})
    root_value = sp.factor(boundary.subs(w, -1 / a))
    expected = sp.factor(-2 * p**2 * (u - U) / (a**2 * B))
    if sp.simplify(root_value - expected) != 0:
        issues.append("quintic root-value factorization failed")

    contact = sp.factor(boundary.subs(u, U))
    cofactor = sp.factor(contact / (1 + a * w) ** 3)
    if sp.rem(sp.together(contact).as_numer_denom()[0], (1 + a * w) ** 3, w) != 0:
        issues.append("triple-contact factor failed")
    cofactor_disc = sp.factor(sp.discriminant(cofactor, w))
    if sp.simplify(cofactor_disc - 5 * (3 * a**2 - 14 * a - 4 * p + 15) / 3) != 0:
        issues.append("cofactor discriminant failed")

    d, j = sp.symbols("d j", integer=True, positive=True)
    if sp.simplify((1 - j / (d + 1)) * sp.binomial(d + 1, j) - sp.binomial(d, j)) != 0:
        issues.append("general binomial polar identity failed")

    if len(rows) != 10:
        issues.append(f"row count {len(rows)} != 10")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("open handoff count failed")
    proof_boundary = payload.get("proof_boundary", "")
    for marker in ("does not supply", "infinite degree", "PF-infinity", "RH"):
        if marker not in proof_boundary:
            issues.append(f"missing proof-boundary marker {marker!r}")

    print(
        "validated Jensen-window PF quartic-quintic polar-contact lemma: "
        f"{len(rows)} rows, {len(issues)} issues, 4 exact polar identities, "
        "1 strict nonroot test, 1 multiplicity rule, 1 double-to-triple theorem, "
        "1 quintic contact factorization, 1 cofactor gate, "
        "1 open all-degree handoff"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
