#!/usr/bin/env python3
"""Validate the q=1 growing contiguous-terminal-prefix theorem."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_growing_prefix_finite_height_gate as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [f"ctgf_{index:02d}_{suffix}" for index, suffix in enumerate((
    "domain", "prefix", "ratio", "phase", "ratio_x", "source", "log_bounds",
    "coordinates", "leading", "individual", "aggregate", "perturbation",
    "budget", "entry", "theorem", "handoff",
), start=1)]
EXPECTED_COUNTS = {
    "rows": 16,
    "exact_ratio_logs": 1,
    "differentiated_ratio_logs": 1,
    "carrier_coordinate_errors": 4,
    "aggregate_coordinate_errors": 4,
    "growing_prefix_exponent_denominator": 18,
    "worst_height_prefix_length": 3,
    "finite_height_current_theorems": 1,
    "bulk_aggregate_closures": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated growing-prefix finite-height gate: 16 rows, 1 exact ratio, "
    "1 differentiated ratio, 4 carrier errors, 4 aggregate errors, "
    "exponent 1/18, M=3 at L=50, current <-1/400, 0 bulk closures"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing result: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid result: {exc}")
        return {}


def independent_algebra(issues: list[str]) -> None:
    h, theta, m = sp.symbols("h theta m", positive=True, real=True)
    k = theta + m
    delta = sp.log(1 - h * k) - sp.log(1 - h * theta)
    phase = (
        2 * sp.pi * delta / h**2
        + 2 * sp.pi * m * (1 / h - theta)
        + sp.pi * (m**2 + 4 * m * theta)
    )
    expected = -2 * sp.pi * (k**3 - theta**3) / 3
    if sp.simplify(sp.limit(phase / h, h, 0, dir="+") - expected) != 0:
        issues.append("independent cubic phase failed")

    t, log_a, alpha = sp.symbols("t log_a alpha", real=True)
    ell_k = sp.log(1 - h * k)
    ell_theta = sp.log(1 - h * theta)
    s = sp.Rational(1, 2) - sp.I * (
        2 * sp.pi / h**2 - sp.pi * t / 8
    )
    raw = (
        t
        * ((log_a + ell_k) ** 2 - (log_a + ell_theta) ** 2)
        / 4
        - (s + t * alpha / 2) * delta
    )
    stated = (
        -delta / 2
        - t * (alpha - log_a) * delta / 2
        + t * (ell_k + ell_theta) * delta / 4
        + sp.I * (phase - sp.pi * t * delta / 8)
    )
    log_r_m = sp.I * sp.pi * m - 4 * sp.I * sp.pi * m * theta
    branch = sp.simplify((raw - stated - log_r_m) / (2 * sp.pi * sp.I))
    branch_expected = -m * (1 / h - theta) - m * (m + 1) / 2
    if sp.simplify(branch - branch_expected) != 0:
        issues.append("independent ratio branch failed")

    ap, app, w0, w1, base = sp.symbols("ap app w0 w1 base")
    d0 = base + ap * (t / 4 + t**2 * w0**2 / 8)
    d1 = base + ap * (t / 4 + t**2 * w1**2 / 8)
    if sp.factor(d1 - d0 - ap * t**2 * (w1 - w0) * (w1 + w0) / 8) != 0:
        issues.append("independent correction difference failed")
    dx0 = app * (t / 4 + t**2 * w0**2 / 8) + t**2 * w0 * ap**2 / 4
    dx1 = app * (t / 4 + t**2 * w1**2 / 8) + t**2 * w1 * ap**2 / 4
    dx_expected = (
        app * t**2 * (w1 - w0) * (w1 + w0) / 8
        + t**2 * ap**2 * (w1 - w0) / 4
    )
    if sp.factor(dx1 - dx0 - dx_expected) != 0:
        issues.append("independent correction x-difference failed")

    c, d, e, n, ec, ed, ee, en = sp.symbols(
        "C D E N ec ed ee en", real=True
    )
    direct = (c + ec) * (n + en) - (d + ed) * (e + ee) - (c * n - d * e)
    target = c * en + n * ec + ec * en - d * ee - e * ed - ed * ee
    if sp.expand(direct - target) != 0:
        issues.append("independent current perturbation failed")


def independent_integer_budget(issues: list[str]) -> None:
    hden = gate.H_DENOMINATOR
    checks = {
        "M3": 4**18 < hden,
        "hK": 10**180 < hden**17,
        "hK3": 10**54 < hden**5,
        "current": 2_400_000**18 < hden**11,
    }
    for name, passed in checks.items():
        if not passed:
            issues.append(f"independent integer budget failed: {name}")

    aggregate = {"C": 124, "D": 300, "E": 268, "N": 940}
    linear = 4 * aggregate["N"] + aggregate["C"] + 2 * aggregate["E"] + 2 * aggregate["D"]
    quadratic = aggregate["C"] * aggregate["N"] + aggregate["D"] * aggregate["E"]
    if linear != 5020 or quadratic != 196_960:
        issues.append("independent aggregate constants failed")


def independent_majorants(issues: list[str]) -> None:
    delta = Fraction(1, 10**10)
    t_max = Fraction(1, 5000)
    phase = Fraction(44, 21) / (1 - delta)
    d_difference = 3 * t_max**2 * 14 / 8
    d_x_difference = t_max**2 * (Fraction(7, 8) + Fraction(9, 4))
    ratio = (
        1
        + 3 * t_max
        + 2 * t_max
        + phase
        + Fraction(22, 7) * t_max / 4
        + 2 * d_difference
    )
    ratio_x = Fraction(1, 2) + Fraction(1, 3000) + Fraction(1, 45)
    checks = {
        "phase": phase < 3,
        "d difference": d_difference < Fraction(1, 10**6),
        "d x difference": d_x_difference < Fraction(1, 10**6),
        "extended d": Fraction(1, 20) + delta * d_difference < Fraction(1, 10),
        "extended d x": Fraction(1, 400) + delta * d_x_difference < Fraction(1, 100),
        "ratio": ratio < 5,
        "ratio x": ratio_x < 1,
        "slope coordinate": Fraction(1, 5000) + Fraction(1, 1500) + 13 < 50,
        "slope x coordinate": (
            Fraction(3, 10000)
            + Fraction(1, 72000)
            + Fraction(1, 2)
            + Fraction(1, 3000)
            + Fraction(1, 2)
            + 9
            < 40
        ),
    }
    for name, passed in checks.items():
        if not passed:
            issues.append(f"independent majorant failed: {name}")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-07-31":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    if [row.get("id") for row in payload.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    for token in ("growing", "q=1", "no bulk", "or RH"):
        if token not in payload.get("status", ""):
            issues.append(f"status missing token: {token}")
    for token in (
        "remaining nonterminal bulk", "q>1", "Abel gap", "contact exclusion",
        "Lambda<=0", "RH",
    ):
        if token not in payload.get("proof_boundary", ""):
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        if payload.get("symbolic_certificate") != gate.symbolic_certificate():
            issues.append("symbolic certificate drifted")
        if payload.get("analytic_majorant_certificate") != gate.analytic_majorant_certificate():
            issues.append("analytic majorant certificate drifted")
        if payload.get("exact_integer_budget") != gate.exact_integer_budget():
            issues.append("integer budget drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"certificate recomputation failed: {exc}")

    exact = payload.get("exact", {})
    if exact.get("theorem") != (
        "For every admissible M>=1, the retained first-order q=1 "
        "contiguous terminal-prefix current satisfies "
        "a^2*J_(tail,M)/S_a^2<-1/400."
    ):
        issues.append("finite-height theorem drifted")

    independent_algebra(issues)
    independent_majorants(issues)
    independent_integer_budget(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Growing-Prefix Finite-Height Gate", "Date: 2026-07-31",
        "floor(a^(1/18))-1", "|W_k|<12*h*k^3", "|Delta P|<6000*h*K^7<1/400",
        "current satisfies a^2*J_(tail,M)/S_a^2<-1/400", "M=3",
        "remaining nonterminal bulk",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        validate(payload, issues)
    validate_note(issues)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
