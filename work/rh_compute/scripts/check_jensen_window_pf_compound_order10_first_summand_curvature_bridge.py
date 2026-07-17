#!/usr/bin/env python3
"""Independently validate the order-ten first/full curvature bridge."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order10_first_summand_curvature_bridge as bridge  # noqa: E402


@dataclass(frozen=True)
class BridgeIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> BridgeIssue:
    return BridgeIssue(section, code, str(detail))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stencil(power: int, start: int) -> Fraction:
    return Fraction(start, start - 1) ** power + 3


def independent_envelope() -> tuple[list[Fraction], Fraction]:
    a = 2 * (Fraction(321, 320) ** 12 + 3)
    ell = 4 * a
    u = 2 * a / 322 + ell * stencil(11, 322)
    v = 2 * ell / 322 + 8 * u
    w = 3 * a / 323**2 + v * stencil(10, 323)
    e = 2 * v / 323 + Fraction(5, 7) * w + ell / 323**2
    z = 4 * a / 324**3 + e * stencil(9, 324)
    y = 2 * e / 324 + v / 324**2 + Fraction(2, 3) * z
    o = 5 * a / 325**4 + y * stencil(8, 325)
    n = 2 * y / 325 + e / 325**2 + Fraction(2, 3) * o
    p = 6 * a / 326**5 + n * stencil(7, 326)
    c = 2 * n / 1249 + y / 1249**2 + Fraction(2, 3) * p
    d = 7 * a / 1250**6 + c * stencil(6, 1250)
    f = 2 * c / 1250 + n / 1250**2 + Fraction(3, 4) * d
    d7 = 8 * a / 1251**7 + f * stencil(5, 1251)
    g = 2 * f / 1251 + c / 1251**2 + 5 * d7
    transfer = g * stencil(4, 1252) / 1252**2
    return [a, ell, u, v, w, e, z, y, o, n, p, c, d, f, d7, g], transfer


def validate_source_hashes(artifact: dict) -> list[BridgeIssue]:
    issues: list[BridgeIssue] = []
    records = artifact.get("source_contract", {}).get("sources", [])
    if len(records) != 7:
        issues.append(issue("sources", "bad-count", len(records)))
        return issues
    for record in records:
        path = bridge.REPO_ROOT / str(record.get("path", ""))
        if not path.exists():
            issues.append(issue("sources", "missing", path))
            continue
        actual = sha256(path)
        if actual != record.get("sha256"):
            issues.append(issue("sources", "sha256-drift", path))
        source = json.loads(path.read_text(encoding="utf-8"))
        if source.get("kind") != record.get("kind"):
            issues.append(issue("sources", "kind-drift", path))
    return issues


def validate_floor(expression: sp.Expr, variable: sp.Symbol, start: int) -> bool:
    numerator, denominator = sp.fraction(sp.factor(expression))
    if denominator.subs(variable, start) <= 0:
        return False
    m = sp.symbols("m", integer=True, nonnegative=True)
    coefficients = sp.Poly(
        sp.expand(numerator.subs(variable, start + m)), m
    ).all_coeffs()
    return bool(coefficients) and all(value > 0 for value in coefficients)


def validate_artifact(path: Path) -> tuple[dict, list[BridgeIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[BridgeIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order10_first_summand_curvature_bridge"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))

    expected_summary = {
        "rows": 8,
        "ready_rows": 7,
        "open_rows": 1,
        "seventh_gap_floor_theorems": 1,
        "full_kernel_transfer_theorems": 1,
        "power_envelope_rows": 16,
        "transfer_constant": 10,
        "conditional_full_constant": 4211,
        "target_constant": 5500,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    rows = artifact.get("rows", [])
    if len(rows) != 8 or len({row.get("id") for row in rows}) != 8:
        issues.append(issue("rows", "bad-rows", rows))
    if sum(row.get("readiness") == "open" for row in rows) != 1:
        issues.append(issue("rows", "bad-open-count", rows))

    exact = artifact.get("exact", {})
    required_exact = {
        "continuous_target": "z_1''(t)<=4200/t^2 for every real t>=1251",
        "seventh_gap_floor": "min(W_j,W_j^(1))>1/(5*j), j>=1251",
        "full_transfer": "|Z_k-Z_k^(1)|<10/k^2 for every integer k>=1252",
        "conditional_full_ceiling": (
            "z_1''(t)<=4200/t^2 on t>=1251 => "
            "Z_k<4201/k^2+10/k^2=4211/k^2<5500/k^2, k>=1252"
        ),
        "preserved_negative_prefix": "Q_(10,n)(-100)<0 for n=0,1,2,3",
    }
    for key, expected in required_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))

    j = sp.symbols("j", integer=True, positive=True)
    first_expression = (
        8 / (2 * j + 1) - sp.Rational(4201, 1) / j**2 - 1 / (5 * j)
    )
    full_expression = (
        sp.Rational(1004, 125) / (2 * j + 1)
        - sp.Rational(4751, 1) / j**2
        - 1 / (5 * j)
    )
    if not validate_floor(first_expression, j, 1251):
        issues.append(issue("floors", "first-floor-failed", first_expression))
    if not validate_floor(full_expression, j, 1251):
        issues.append(issue("floors", "full-floor-failed", full_expression))

    constants, transfer = independent_envelope()
    stored_rows = exact.get("power_envelope", {}).get("rows", [])
    expected_starts = [
        321,
        321,
        322,
        322,
        323,
        323,
        324,
        324,
        325,
        325,
        326,
        1249,
        1250,
        1250,
        1251,
        1251,
    ]
    expected_powers = [12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6, 5, 5, 4]
    if len(stored_rows) != len(constants):
        issues.append(issue("envelope", "bad-row-count", len(stored_rows)))
    else:
        for index, (stored, expected) in enumerate(zip(stored_rows, constants)):
            if stored.get("start") != expected_starts[index]:
                issues.append(issue("envelope", f"start-drift-{index}", stored.get("start")))
            if stored.get("power") != expected_powers[index]:
                issues.append(issue("envelope", f"power-drift-{index}", stored.get("power")))
            try:
                actual = Fraction(stored.get("constant", ""))
            except Exception as exc:
                issues.append(issue("envelope", f"bad-constant-{index}", exc))
                continue
            if actual != expected:
                issues.append(issue("envelope", f"constant-drift-{index}", actual))
    stored_transfer = exact.get("power_envelope", {}).get("transfer_scaled_exact")
    try:
        actual_transfer = Fraction(stored_transfer)
    except Exception as exc:
        issues.append(issue("envelope", "bad-transfer", exc))
    else:
        if actual_transfer != transfer:
            issues.append(issue("envelope", "transfer-drift", actual_transfer))
    if transfer >= 10:
        issues.append(issue("envelope", "transfer-not-below-10", transfer))

    issues.extend(validate_source_hashes(artifact))
    try:
        canonical = bridge.build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(issue("rebuild", "exception", exc))
    else:
        if canonical != artifact:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[BridgeIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact full-kernel transfer with one open continuous",
        "min(W_j,W_j^(1))>1/(5*j)",
        "|Z_k-Z_k^(1)|<10/k^2",
        "z_1''(t)<=4200/t^2",
        "four negative endpoint rows remain negative",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-ten entry is proved",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=bridge.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=bridge.DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER10-BRIDGE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-ten first/full curvature bridge: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    scaled = artifact["exact"]["power_envelope"]["transfer_scaled_decimal"]
    print(
        "validated order-ten first/full curvature bridge: "
        f"{summary['power_envelope_rows']} exact envelope rows, "
        f"scaled transfer {scaled}<10, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
