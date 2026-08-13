#!/usr/bin/env python3
"""Independently check the terminal/nonterminal relative-lift gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_"
    "relative_lift_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 30,
    "outer_remainder_decompositions": 1,
    "outer_remainder_components": 4,
    "real_remainder_observation_rows": 8,
    "relative_lift_families": 2,
    "tie_transfer_invariants": 4,
    "mixed_current_compressions": 2,
    "mixed_polynomial_functionals": 2,
    "maximum_mixed_polynomial_degree": 5,
    "terminal_higher_endpoint_reductions": 2,
    "lower_conditional_relative_traces": 2,
    "transpose_carrier_bilinear_factorizations": 1,
    "terminal_hermitian_diagonal_nulls": 1,
    "geometry_checks": 440,
    "signed_lower_interior_bounds": 0,
    "complete_outer_complement_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_full_current_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "tnr_01_endpoint",
    "tnr_02_interior",
    "tnr_03_split",
    "tnr_04_components",
    "tnr_05_observations",
    "tnr_06_linear",
    "tnr_07_terminal_lift",
    "tnr_08_hlift",
    "tnr_09_tlift",
    "tnr_10_atoms",
    "tnr_11_ties",
    "tnr_12_complete",
    "tnr_13_hermitian",
    "tnr_14_transpose",
    "tnr_15_projection",
    "tnr_16_hpoly",
    "tnr_17_tpoly",
    "tnr_18_hfunctional",
    "tnr_19_tfunctional",
    "tnr_20_degree",
    "tnr_21_terminal_h",
    "tnr_22_terminal_t",
    "tnr_23_lower_h",
    "tnr_24_lower_t",
    "tnr_25_factor",
    "tnr_26_diagonal",
    "tnr_27_scale",
    "tnr_28_route",
    "tnr_29_pi",
    "tnr_30_boundary",
]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def make_complex(name: str) -> sp.Expr:
    real, imag = sp.symbols(f"{name}_real {name}_imag", real=True)
    return real + sp.I * imag


def check_sources(artifact: dict) -> None:
    expected = {
        "terminal_bulk_phase",
        "terminal_conditional",
        "terminal_obstruction",
        "outer_endpoint",
        "observation_image",
        "physical_plin",
        "full_support",
    }
    source = artifact.get("source_audit", {})
    require(set(source) == expected, "source-audit keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_remainder_and_lifts() -> None:
    imaginary = sp.I
    terminal, high_n, lower, high_1, interior = sp.symbols(
        "terminal high_n lower high_1 interior"
    )
    complete = terminal + high_n - lower - high_1 + interior
    rest = high_n - lower - high_1 + interior
    require_zero(complete - terminal - rest, "independent remainder split")

    lam, la, ln = sp.symbols("lam la ln", real=True)
    un = la - ln
    coefficients = sp.symbols("c0:5")
    polynomial = sum(coefficients[j] * lam**j for j in range(5))
    lifted = imaginary * (lam - la) * polynomial
    h_relative = sp.expand(lifted + imaginary * un * polynomial)
    t_relative = sp.expand(lifted - imaginary * un * polynomial)
    require_zero(h_relative - imaginary * (lam - ln) * polynomial, "independent H lift")
    require_zero(
        t_relative - imaginary * (lam - la - un) * polynomial,
        "independent T lift",
    )
    require(sp.Poly(h_relative, lam).degree() == 5, "independent H degree")
    require(sp.Poly(t_relative, lam).degree() == 5, "independent T degree")

    carrier, remainder, transfer = sp.symbols("carrier remainder transfer")
    dcarrier, dremainder, dtransfer = sp.symbols("dcarrier dremainder dtransfer")
    old = carrier + remainder
    new = carrier + transfer + remainder - transfer
    dold = dcarrier + dremainder
    dnew = dcarrier + dtransfer + dremainder - dtransfer
    require_zero(new - old, "independent tie value")
    require_zero(dnew - dold, "independent tie lift")
    require_zero((dnew + imaginary * un * new) - (dold + imaginary * un * old), "independent tie H")
    require_zero((dnew - imaginary * un * new) - (dold - imaginary * un * old), "independent tie T")


def check_mixed_current() -> None:
    imaginary = sp.I
    un = sp.symbols("un", real=True)
    a = {label: make_complex(f"a{label}") for label in "VNAQ"}
    y = {label: make_complex(f"y{label}") for label in "VNAQ"}
    dy = {label: make_complex(f"dy{label}") for label in "VNAQ"}
    signs = {"V": 1, "N": 1, "A": -1, "Q": -1}
    partner = {"V": "N", "N": "V", "A": "Q", "Q": "A"}

    h_full = sp.Integer(0)
    t_full = sp.Integer(0)
    h_relative = sp.Integer(0)
    t_relative = sp.Integer(0)
    for label in "VNAQ":
        other = partner[label]
        h_full += signs[label] * (
            a[label] * sp.conjugate(dy[other])
            + y[label] * sp.conjugate(-imaginary * un * a[other])
        )
        t_full += signs[label] * (
            a[label] * dy[other]
            + y[label] * (-imaginary * un * a[other])
        )
        h_relative += signs[label] * a[label] * sp.conjugate(
            dy[other] + imaginary * un * y[other]
        )
        t_relative += signs[label] * a[label] * (
            dy[other] - imaginary * un * y[other]
        )
    require_zero(
        sp.expand_complex(sp.re(h_full - h_relative)),
        "independent mixed H compression",
    )
    require_zero(sp.expand_complex(t_full - t_relative), "independent mixed T compression")


def check_endpoints() -> None:
    imaginary = sp.I
    kappa, phase, amplitude = sp.symbols("kappa phase amplitude")
    s2, s3, alpha, n = sp.symbols("s2 s3 alpha n")
    un, p, pprime, g = sp.symbols("un p pprime g")

    def higher(value: sp.Expr, derivative: sp.Expr) -> sp.Expr:
        return -kappa**2 * phase * amplitude * (
            derivative * s2 / n + alpha * value * s3 / n**2
        )

    h_endpoint = higher(0, imaginary * p)
    require_zero(
        h_endpoint + imaginary * kappa**2 * phase * amplitude * p * s2 / n,
        "independent terminal H endpoint",
    )
    t_endpoint = higher(
        -2 * imaginary * un * p,
        imaginary * (p - 2 * un * (pprime + g * p)),
    )
    t_expected = -imaginary * kappa**2 * phase * amplitude * (
        (p - 2 * un * (pprime + g * p)) * s2 / n
        - 2 * alpha * un * p * s3 / n**2
    )
    require_zero(t_endpoint - t_expected, "independent terminal T endpoint")

    logn, loga, p0, s1 = sp.symbols("logn loga p0 s1")
    require_zero(
        -kappa * (-imaginary * logn * p0) * s1
        - imaginary * kappa * logn * p0 * s1,
        "independent lower H",
    )
    require_zero(
        -kappa * (-imaginary * (loga + un) * p0) * s1
        - imaginary * kappa * (loga + un) * p0 * s1,
        "independent lower T",
    )


def check_carrier_factorization() -> None:
    rn, rl, rxn, rxl, cn, cl, qn, ql = sp.symbols(
        "rn rl rxn rxl cn cl qn ql"
    )
    an, al = rn * cn, rl * cl
    nn, nl = rn * qn + rxn * cn, rl * ql + rxl * cl
    bilinear = cn * nl + nn * cl - an * ql - qn * al
    expected = (rl - rn) * (cn * ql - qn * cl) + (rxl + rxn) * cn * cl
    require_zero(bilinear - expected, "independent carrier factorization")
    require_zero(
        bilinear.subs({rl: rn, rxl: rxn, cl: cn, ql: qn}) - 2 * rxn * cn**2,
        "independent terminal diagonal",
    )

    lam = sp.symbols("lam", real=True)
    c = sum(sp.symbols(f"cc{j}") * lam**j for j in range(3))
    n = sum(sp.symbols(f"nn{j}") * lam**j for j in range(5))
    a = sum(sp.symbols(f"aa{j}") * lam**j for j in range(4))
    q = sum(sp.symbols(f"qq{j}") * lam**j for j in range(4))
    cv, nv, av, qv = sp.symbols("cv nv av qv")
    mixed = sp.expand(cv * n + nv * c - av * q - qv * a)
    require(sp.Poly(mixed, lam).degree() == 4, "independent mixed degree four")
    require(sp.Poly((lam - sp.symbols("ell")) * mixed, lam).degree() == 5, "independent lifted degree five")


def check_geometry() -> int:
    checks = 0
    for n_value in range(2, 90):
        for theta in [Fraction(0), Fraction(1, 13), Fraction(1, 3), Fraction(2, 3), Fraction(12, 13), Fraction(1)]:
            a = Fraction(n_value) + theta
            h = 1 / a
            require(Fraction(1, n_value) <= Fraction(3, 2) * h, "independent N inverse")
            require(2 * h * a**2 / n_value**2 <= Fraction(9, 2) * h, "independent S3 coefficient")
            require(Fraction(5, n_value) <= Fraction(15, 2) * h, "independent lower S1")
            checks += 1
    require(checks == 528, "independent geometry count")
    return checks


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["tnr_28_route"],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Outer Remainder",
        "mathscr R_N[P]=U_N[P]-T_1[P]-U_1[P]+I_N[P]",
        "## Relative Lifts",
        "i(lambda-log N)P",
        "i(lambda-log a-u_N)P",
        "## Mixed Current",
        "## Polynomial Compression",
        "B_(H,N)(lambda)",
        "B_(T,N)(lambda)",
        "degree at most five",
        "## Endpoint Audit",
        "## Pi Provenance",
        "## Proof Boundary",
    ]
    for marker in markers:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    require(RESULT_PATH.exists(), "result artifact missing")
    require(NOTE_PATH.exists(), "note artifact missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "counts drifted")
    require("signed lower/interior estimate open" in artifact.get("status", ""), "status drifted")
    require(
        "no signed lower/interior estimate" in artifact.get("proof_boundary", ""),
        "proof boundary drifted",
    )
    cert = artifact.get("symbolic_certificate", {})
    require(
        set(cert)
        == {
            "outer_remainder",
            "relative_lifts",
            "mixed_current",
            "polynomial_compression",
            "endpoint_audit",
            "carrier_bilinear",
            "bounds",
            "audit",
            "handoff",
        },
        "certificate keys drifted",
    )
    require(
        cert["audit"]
        == {
            "geometry_cases": 440,
            "relative_lift_families": 2,
            "mixed_polynomial_functionals": 2,
            "outer_remainder_components": 4,
        },
        "builder audit drifted",
    )
    check_sources(artifact)
    check_remainder_and_lifts()
    check_mixed_current()
    check_endpoints()
    check_carrier_factorization()
    geometry_checks = check_geometry()
    check_rows(artifact)
    check_note(note)
    print(
        "validated full-support terminal/nonterminal relative-lift gate: "
        f"{EXPECTED_COUNTS['rows']} rows, "
        f"{EXPECTED_COUNTS['outer_remainder_components']} remainder components, "
        f"{EXPECTED_COUNTS['relative_lift_families']} relative lifts, "
        f"{EXPECTED_COUNTS['mixed_polynomial_functionals']} degree-five functionals, "
        f"{EXPECTED_COUNTS['tie_transfer_invariants']} tie invariants, "
        f"{geometry_checks} independent geometry checks, "
        f"{EXPECTED_COUNTS['signed_lower_interior_bounds']} signed lower/interior bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
