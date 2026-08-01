#!/usr/bin/env python3
"""Validate the five-current Mangoldt normal-form gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_current_mangoldt_normal_form_gate as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"ttmg_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "correction",
            "current",
            "base",
            "transform",
            "mangoldt",
            "heat",
            "hyperbola",
            "join",
            "guard",
            "target",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 12,
    "exact_correction_transforms": 1,
    "corrected_complex_currents": 5,
    "correction_free_complex_moments": 5,
    "mangoldt_moment_orders": 4,
    "balanced_hyperbola_schemas": 1,
    "negative_prime_edge_minor_families": 4,
    "signed_type_ii_targets": 1,
    "signed_type_ii_bounds": 0,
    "phi_b_bounds": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated five-current Mangoldt normal form: 12 rows, "
    "5 correction-free moments, 4 Mangoldt orders, "
    "4 indefinite minor families, 0 signed Type-I/II bounds"
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


def independent_transform(issues: list[str]) -> None:
    l = sp.symbols("l", real=True)
    r, s, rx, sx, q = sp.symbols("r s rx sx q")
    correction = 1 + r * l + s * l**2
    correction_x = rx * l + sx * l**2
    delta = correction_x / correction
    if sp.simplify(delta * correction * q - correction_x * q) != 0:
        issues.append("independent correction-current cancellation failed")

    g = sp.symbols("g_0:5")
    matrix = sp.Matrix(
        [
            [1, r, s, 0, 0],
            [0, 1, r, s, 0],
            [0, 0, 1, r, s],
            [0, rx, sx, 0, 0],
            [0, 0, rx, sx, 0],
        ]
    )
    out = matrix * sp.Matrix(g)
    expected = sp.Matrix(
        [
            g[0] + r * g[1] + s * g[2],
            g[1] + r * g[2] + s * g[3],
            g[2] + r * g[3] + s * g[4],
            rx * g[1] + sx * g[2],
            rx * g[2] + sx * g[3],
        ]
    )
    if any(sp.expand(value) != 0 for value in out - expected):
        issues.append("independent five-current transform failed")
    determinant = sp.factor(matrix.det())
    expected_det = sp.factor(s * (-r * rx * sx + s * rx**2 + sx**2))
    if sp.expand(determinant - expected_det) != 0:
        issues.append("independent transform determinant failed")


def independent_arithmetic(issues: list[str]) -> None:
    for cutoff in (19, 43, 71):
        for n in range(2, cutoff + 1):
            factors = sp.factorint(n)
            symbols = {
                prime: sp.Symbol(f"P_{prime}") for prime in factors
            }
            log_n = sum(
                exponent * symbols[prime]
                for prime, exponent in factors.items()
            )

            def mangoldt(integer: int) -> sp.Expr:
                parts = sp.factorint(integer)
                if len(parts) != 1:
                    return sp.Integer(0)
                prime = next(iter(parts))
                return symbols.get(prime, sp.Symbol(f"P_{prime}"))

            divisors = sp.divisors(n)
            one_sided = sum(mangoldt(divisor) for divisor in divisors)
            symmetric = sum(
                (mangoldt(divisor) + mangoldt(n // divisor)) / 2
                for divisor in divisors
            )
            if sp.expand(one_sided - log_n) != 0:
                issues.append(f"independent divisor identity failed at {n}")
            if sp.expand(symmetric - log_n) != 0:
                issues.append(f"independent symmetric identity failed at {n}")

        root = int(cutoff**0.5)
        all_pairs = {
            (d, m)
            for d in range(1, cutoff + 1)
            for m in range(1, cutoff // d + 1)
        }
        square = {(d, m) for d, m in all_pairs if d <= root and m <= root}
        wings = all_pairs - square
        if any(d > root and m > root for d, m in wings):
            issues.append(f"independent hyperbola wing failed at {cutoff}")
        if square & wings or square | wings != all_pairs:
            issues.append(f"independent hyperbola partition failed at {cutoff}")

    lp = sp.symbols("lp", positive=True, real=True)
    for order in range(1, 5):
        matrix = sp.Matrix([[0, lp**order / 2], [lp**order / 2, 0]])
        if sp.factor(matrix.det()) != -lp ** (2 * order) / 4:
            issues.append(f"independent prime-edge minor failed: {order}")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-07-31":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    if sum(row.get("readiness") == "open" for row in rows) != 2:
        issues.append("open-row count drifted")

    status = payload.get("status", "")
    for token in (
        "exact terminal-tail five-current",
        "Type-I/II Phi_B bound open",
        "no bulk aggregate sign closure",
        "or RH",
    ):
        if token not in status:
            issues.append(f"status missing token: {token}")
    boundary = payload.get("proof_boundary", "")
    for token in (
        "no Type-I/II cancellation estimate",
        "signed bound on Phi_B",
        "Xi-level current theorem",
        "Lambda<=0",
        "RH",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        sources = gate.load_sources()
        if payload.get("source_audit") != gate.source_audit(sources):
            issues.append("source audit recomputation drifted")
        if payload.get("symbolic_certificate") != gate.symbolic_certificate():
            issues.append("symbolic certificate drifted")
        if payload.get("arithmetic_certificate") != gate.arithmetic_certificate():
            issues.append("arithmetic certificate drifted")
        if payload.get("exact") != gate.exact_payload():
            issues.append("exact payload drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"certificate recomputation failed: {exc}")

    if payload.get("counts", {}).get("signed_type_ii_bounds") != 0:
        issues.append("unexpected signed Type-I/II bound")

    independent_transform(issues)
    independent_arithmetic(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Terminal-Tail Five-Current Mangoldt Normal Form",
        "Date: 2026-07-31",
        "0 signed Type-I/II",
        "G_r^(B)",
        "H_2=G_2+rho_1G_3+rho_2G_4",
        "For 1<=r<=4",
        "-log(p)^8/4<0",
        "Phi_B<=1/400",
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
