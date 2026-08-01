#!/usr/bin/env python3
"""Independently validate the complex endpoint normalization corrigendum."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_complex_endpoint_source_normalization_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "c0_source_extraction": (
        REPO_ROOT
        / "work/rh_compute/scripts/"
        "jensen_window_pf_newman_polymath15_"
        "critical_RS_C1_endpoint_peeling_contract.py"
    ),
    "absolute_phase_anchor_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_absolute_phase_anchor_reduction.json"
    ),
    "direct_projection_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_direct_projection_regime_reduction.json"
    ),
    "carrier_kernel_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_carrier_kernel_abel_prefix_reduction.json"
    ),
    "first_pivot_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_endpoint_first_pivot_odd_small_ball_guard.json"
    ),
    "odd_fibre_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_endpoint_odd_fibre_"
        "correlation_feasibility_guard.json"
    ),
    "contact_transport_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_contact_signed_transport_reduction.json"
    ),
    "interior_current_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_interior_projective_current_gate.json"
    ),
}
EXPECTED_IDS = [
    "cesn_01_published_c0",
    "cesn_02_midpoint_parity",
    "cesn_03_midpoint_nonreal",
    "cesn_04_endpoint_occurrence",
    "cesn_05_endpoint_factorization",
    "cesn_06_value_projection",
    "cesn_07_jet_projection",
    "cesn_08_hermitian_product",
    "cesn_09_centered_atom",
    "cesn_10_exceptional_fibre",
    "cesn_11_odd_fibre_rotation",
    "cesn_12_endpoint_carrier_minor",
    "cesn_13_carrier_carrier_minor",
    "cesn_14_vaughan_rank",
    "cesn_15_contact_rank",
    "cesn_16_supersession",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate_source_formula(issues: list[str]) -> None:
    p = sp.symbols("p", real=True)
    a = sp.symbols("a", positive=True)
    pi = sp.pi
    c0 = (
        sp.exp(sp.I * pi * (p**2 / 2 + sp.Rational(3, 8)))
        - sp.I * sp.sqrt(2) * sp.cos(pi * p / 2)
    ) / (2 * sp.cos(pi * p))

    if sp.simplify(c0.subs(p, -p) - c0) != 0:
        issues.append("C_0 evenness failed")
    third_at_zero = sp.simplify(sp.diff(c0, p, 3).subs(p, 0))
    if third_at_zero != 0:
        issues.append("C_0'''(0) failed")

    midpoint = sp.expand_complex(c0.subs(p, 0))
    expected_real = sp.sqrt(2 - sp.sqrt(2)) / 4
    expected_imag = (
        sp.sqrt(2 + sp.sqrt(2)) / 4 - sp.sqrt(2) / 2
    )
    if sp.simplify(sp.re(midpoint) - expected_real) != 0:
        issues.append("C_0(0) real radical failed")
    if sp.simplify(sp.im(midpoint) - expected_imag) != 0:
        issues.append("C_0(0) imaginary radical failed")
    if expected_imag.is_negative is not True:
        issues.append("C_0(0) exact imaginary sign failed")

    h_midpoint = sp.simplify(
        c0.subs(p, 0)
        + sp.diff(c0, p, 3).subs(p, 0) / (12 * pi**2 * a)
    )
    if sp.simplify(h_midpoint - c0.subs(p, 0)) != 0:
        issues.append("H_a(0)=C_0(0) failed")

    source_text = SOURCE_PATHS["c0_source_extraction"].read_text(
        encoding="utf-8"
    )
    for required in (
        "sp.exp(sp.I * pi * (p**2 / 2 + sp.Rational(3, 8)))",
        "sp.I * sp.sqrt(2) * sp.cos(pi * p / 2)",
        "2 * sp.cos(pi * p)",
        "sp.diff(f, p, 3) / (12 * pi**2)",
    ):
        if required not in source_text:
            issues.append(f"local C_0 source extraction drifted: {required}")


def validate_endpoint_algebra(issues: list[str]) -> None:
    kappa, t0, hr, hi, jr, ji, c, b, u = sp.symbols(
        "kappa T_0 H_R H_I J_R J_I c b u_N", real=True
    )
    h = hr + sp.I * hi
    j = jr + sp.I * ji
    e = sp.expand_complex(kappa * (t0 + sp.I) * h)
    g = sp.expand_complex(kappa * (t0 + sp.I) * j)
    projections = (
        (sp.re(e), kappa * (t0 * hr - hi), "Re e"),
        (sp.im(e), kappa * (t0 * hi + hr), "Im e"),
        (sp.re(g), kappa * (t0 * jr - ji), "Re g"),
        (sp.im(g), kappa * (t0 * ji + jr), "Im g"),
    )
    for actual, expected, label in projections:
        if sp.simplify(actual - expected) != 0:
            issues.append(f"{label} projection failed")

    product = sp.expand_complex(g * sp.conjugate(e))
    expected_product = sp.expand_complex(
        kappa**2 * (t0**2 + 1) * j * sp.conjugate(h)
    )
    if sp.simplify(product - expected_product) != 0:
        issues.append("g*conj(e) correction failed")

    c0_atom = kappa * (t0 * hr - hi)
    d0_atom = kappa * (
        t0 * jr - ji - c * u * (t0 * hr - hi)
    )
    if sp.simplify(d0_atom - (sp.re(g) - c * u * sp.re(e))) != 0:
        issues.append("centered endpoint atom failed")

    nonzero_h_projection_zero = {
        t0: 2,
        hr: 1,
        hi: 2,
    }
    if sp.simplify((t0 * hr - hi).subs(nonzero_h_projection_zero)) != 0:
        issues.append("exceptional projection witness failed")
    if sp.simplify((hr**2 + hi**2).subs(nonzero_h_projection_zero)) == 0:
        issues.append("exceptional projection witness has H=0")


def validate_odd_rotation(issues: list[str]) -> None:
    a, b, x, y = sp.symbols("A B X Y", real=True)
    d = sp.sqrt(a**2 + b**2)
    p = (a * x + b * y) / d
    q = (-b * x + a * y) / d
    difference = sp.factor(
        (p + d) ** 2 + q**2 - ((x + a) ** 2 + (y + b) ** 2)
    )
    if sp.simplify(difference) != 0:
        issues.append("complex endpoint odd-fibre rotation failed")

    cross = sp.re(
        (x + sp.I * y) * sp.conjugate(a + sp.I * b)
    )
    if sp.simplify(cross - (a * x + b * y)) != 0:
        issues.append("odd-fibre cross projection failed")


def validate_contact_minors(issues: list[str]) -> None:
    kappa, t0, hr, hi, jr, ji, c, b, un, u_terminal = sp.symbols(
        "kappa T_0 H_R H_I J_R J_I c b u_n u_N",
        real=True,
    )
    x, y = sp.symbols("X_n Y_n", real=True)
    endpoint_value = kappa * (t0 * hr - hi)
    endpoint_slope = kappa * (
        t0 * jr
        - ji
        - c * u_terminal * (t0 * hr - hi)
    )
    carrier_slope = c * (un - u_terminal) * x - b * un * y
    actual = sp.expand(endpoint_value * carrier_slope - endpoint_slope * x)
    expected = sp.expand(
        kappa
        * (
            (t0 * hr - hi) * un * (c * x - b * y)
            - (t0 * jr - ji) * x
        )
    )
    if sp.simplify(actual - expected) != 0:
        issues.append("endpoint-carrier contact minor failed")

    positive = expected.subs(
        {
            kappa: 1,
            t0: 1,
            hr: 1,
            hi: 0,
            jr: 0,
            ji: 0,
            un: 1,
            c: 1,
            b: -1,
            x: 1,
            y: 0,
        }
    )
    negative = expected.subs(
        {
            kappa: 1,
            t0: 1,
            hr: 1,
            hi: 0,
            jr: 0,
            ji: 0,
            un: 1,
            c: 1,
            b: -1,
            x: -1,
            y: 0,
        }
    )
    if positive != 1 or negative != -1:
        issues.append("endpoint-carrier sign witnesses failed")

    xn, yn, xm, ym, um = sp.symbols(
        "X_n2 Y_n2 X_m Y_m u_m", real=True
    )
    dn = c * (un - u_terminal) * xn - b * un * yn
    dm = c * (um - u_terminal) * xm - b * um * ym
    carrier_minor = sp.expand(xn * dm - dn * xm)
    carrier_expected = sp.expand(
        c * (um - un) * xn * xm
        - b * (um * xn * ym - un * yn * xm)
    )
    if sp.simplify(carrier_minor - carrier_expected) != 0:
        issues.append("carrier-carrier contact minor failed")
    carrier_negative = carrier_expected.subs(
        {
            c: 1,
            b: -1,
            un: 2,
            um: 1,
            xn: 1,
            xm: 1,
            yn: 0,
            ym: 0,
        }
    )
    carrier_positive = carrier_expected.subs(
        {
            c: 1,
            b: -1,
            un: 2,
            um: 1,
            xn: 1,
            xm: -1,
            yn: 0,
            ym: 0,
        }
    )
    if carrier_negative != -1 or carrier_positive != 1:
        issues.append("carrier-carrier sign witnesses failed")


def validate_rank_guards(issues: list[str]) -> None:
    sigma = sp.Matrix([1, -1, -1])
    kernel = sigma * sigma.T
    expected = sp.Matrix(
        [
            [1, -1, -1],
            [-1, 1, 1],
            [-1, 1, 1],
        ]
    )
    if kernel != expected:
        issues.append("Vaughan kernel entries failed")
    if kernel.rank() != 1:
        issues.append("Vaughan kernel rank failed")
    if kernel.eigenvals() != {sp.Integer(3): 1, sp.Integer(0): 2}:
        issues.append("Vaughan kernel eigenvalues failed")

    c1, c2, c3, d1, d2, d3 = sp.symbols(
        "c1 c2 c3 d1 d2 d3", real=True
    )
    cvec = sp.Matrix([c1, c2, c3])
    dvec = sp.Matrix([d1, d2, d3])
    joined = cvec * cvec.T + dvec * dvec.T
    if joined.rank(iszerofunc=lambda value: value == 0) > 2:
        issues.append("joined contact/slope rank failed")

    x1, x2, x3 = sp.symbols("x1 x2 x3", real=True)
    xvec = sp.Matrix([x1, x2, x3])
    quadratic = sp.expand((xvec.T * joined * xvec)[0])
    expected_quadratic = sp.expand(
        (cvec.dot(xvec)) ** 2 + (dvec.dot(xvec)) ** 2
    )
    if sp.simplify(quadratic - expected_quadratic) != 0:
        issues.append("joined contact/slope pullback failed")


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for required in (
        "# Complex Endpoint Source-Normalization Corrigendum",
        "Status: exact correction gate and nonpromotion audit.",
        "C_0 and H_a are generally complex",
        "not the Vaughan low coefficient C_0(n)",
        "T_0H_R-H_I=0",
        "g*conj(e)=kappa^2(T_0^2+1)J_a*conj(H_a)",
        "rank 1 and eigenvalues 3,0,0",
        "7 historical artifacts quarantined",
        "not a proof",
    ):
        if required not in text:
            issues.append(f"note missing text: {required}")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "complex endpoint result", issues)
    for key, path in SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source {key}: {path}")
    if issues:
        return issues

    if stored.get("kind") != STEM:
        issues.append("result kind drifted")
    if stored.get("date") != "2026-07-31":
        issues.append("result date drifted")
    if [row.get("id") for row in stored.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    expected_summary = {
        "rows": 16,
        "source_nonreal_witnesses": 1,
        "corrected_cartesian_projections": 4,
        "corrected_hermitian_products": 1,
        "corrected_centered_endpoint_atoms": 1,
        "corrected_odd_fibre_rotations": 1,
        "exact_contact_minors": 2,
        "rank_guards": 2,
        "historical_artifacts_quarantined": 7,
        "signed_lower_bounds": 0,
        "abel_gaps": 0,
        "winding_bounds": 0,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    stored_hashes = (
        stored.get("source_audit", {}).get("source_sha256", {})
    )
    current_hashes = {
        key: file_hash(path) for key, path in SOURCE_PATHS.items()
    }
    if stored_hashes != current_hashes:
        issues.append("source hashes drifted")

    exact = stored.get("exact", {})
    if "generally complex" not in (
        exact.get("source_normalization", {}).get("source_semantics", "")
    ):
        issues.append("complex source semantics missing")
    if "not the Vaughan" not in (
        exact.get("source_normalization", {}).get("notation_guard", "")
    ):
        issues.append("endpoint/Vaughan notation guard missing")
    if exact.get("complex_endpoint", {}).get("hermitian_product") != (
        "g*conj(e)=kappa^2(T_0^2+1)J_a*conj(H_a)"
    ):
        issues.append("stored endpoint Hermitian product drifted")
    if exact.get("centered_endpoint_atom", {}).get("exceptional_fibre") != (
        "The division-free endpoint fibre is "
        "T_0H_R-H_I=0; H_a=0 is only a subfibre."
    ):
        issues.append("stored exceptional fibre drifted")
    if len(
        exact.get("supersession", {}).get("historical_artifacts", [])
    ) != 7:
        issues.append("historical quarantine count drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "does not prove an Xi phase law",
        "Abel-scalar gap",
        "Lambda<=0",
        "RH",
    ):
        if required not in boundary:
            issues.append(f"proof boundary missing: {required}")

    validate_source_formula(issues)
    validate_endpoint_algebra(issues)
    validate_odd_rotation(issues)
    validate_contact_minors(issues)
    validate_rank_guards(issues)
    validate_note(issues)
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated complex endpoint source-normalization gate: "
        "16 rows, 1 nonreal source witness, "
        "4 corrected Cartesian projections, 2 contact minors, "
        "2 rank guards, 7 historical artifacts quarantined, "
        "0 signed lower bounds, 0 Abel gaps, 0 winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
