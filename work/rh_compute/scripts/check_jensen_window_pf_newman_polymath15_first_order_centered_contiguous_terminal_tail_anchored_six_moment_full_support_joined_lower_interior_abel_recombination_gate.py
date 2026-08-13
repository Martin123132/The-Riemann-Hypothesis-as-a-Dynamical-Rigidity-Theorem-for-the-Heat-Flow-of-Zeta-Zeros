#!/usr/bin/env python3
"""Independently check the joined lower/interior recombination gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_joined_lower_interior_"
    "abel_recombination_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 28,
    "physical_kernel_factorizations": 2,
    "ideal_kernel_identities": 6,
    "lower_endpoint_identities": 2,
    "modewise_recombinations": 1,
    "global_recombinations": 4,
    "tie_invariants": 2,
    "relative_recombinations": 4,
    "abel_checks": 50,
    "signed_finite_cell_defect_bounds": 0,
    "signed_completed_current_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "jli_01_hfull",
    "jli_02_hphysical",
    "jli_03_tfull",
    "jli_04_tphysical",
    "jli_05_ideal",
    "jli_06_ideal_lifts",
    "jli_07_atoms",
    "jli_08_corrections",
    "jli_09_lower",
    "jli_10_lengths",
    "jli_11_mode",
    "jli_12_outer",
    "jli_13_objects",
    "jli_14_exterior",
    "jli_15_tie",
    "jli_16_hrelative",
    "jli_17_trelative",
    "jli_18_carrier_guard",
    "jli_19_atoms",
    "jli_20_abel",
    "jli_21_habel",
    "jli_22_tabel",
    "jli_23_abel_guard",
    "jli_24_decision",
    "jli_25_route",
    "jli_26_reserve",
    "jli_27_pi",
    "jli_28_boundary",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.expand(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(
        set(source) == {"relative_lift", "finite_cell", "outer_endpoint"},
        "source-audit keys drifted",
    )
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_physical_factorizations() -> int:
    # Use a different expansion order from the builder.
    r, rn, rx, rxn = sp.symbols("r rn rx rxn")
    c, cn, d, dn, delta = sp.symbols("c cn d dn delta")
    brn, brxn, bcn, bdn, bdelta = sp.symbols("brn brxn bcn bdn bdelta")
    q = r * c + delta * c + d
    qn = rn * cn + delta * cn + dn
    bqn = brn * bcn + bdelta * bcn + bdn
    n = r * q + rx * c
    nn = rn * qn + rxn * cn
    bnn = brn * bqn + brxn * bcn

    direct_h = bcn * n + bnn * c - brn * bcn * q - bqn * r * c
    determinant_h = bcn * q - bqn * c
    factored_h = (r - brn) * determinant_h + (rx + brxn) * bcn * c
    require_zero(direct_h - factored_h, "independent Hermitian first factorization")
    expanded_h = (
        bcn
        * c
        * ((r - brn) ** 2 + (r - brn) * (delta - bdelta) + rx + brxn)
        + (r - brn) * (bcn * d - bdn * c)
    )
    require_zero(direct_h - expanded_h, "independent Hermitian physical expansion")

    direct_t = cn * n + nn * c - rn * cn * q - qn * r * c
    determinant_t = cn * q - qn * c
    factored_t = (r - rn) * determinant_t + (rx + rxn) * cn * c
    require_zero(direct_t - factored_t, "independent transpose first factorization")
    expanded_t = cn * c * ((r - rn) ** 2 + rx + rxn) + (r - rn) * (
        cn * d - dn * c
    )
    require_zero(direct_t - expanded_t, "independent transpose physical expansion")

    substitution_checks = 0
    symbols = [r, rn, rx, rxn, c, cn, d, dn, delta, brn, brxn, bcn, bdn, bdelta]
    for case in range(1, 129):
        values = {
            symbol: sp.Rational(((case + 3 * index) % 17) - 8, index + 2)
            + sp.I * sp.Rational(((2 * case + 5 * index) % 19) - 9, index + 3)
            for index, symbol in enumerate(symbols)
        }
        require_zero((direct_h - expanded_h).subs(values), f"Hermitian substitution {case}")
        require_zero((direct_t - expanded_t).subs(values), f"transpose substitution {case}")
        substitution_checks += 2
    return substitution_checks


def check_ideal_kernels() -> None:
    x, un, unx, uq = sp.symbols("x un unx uq", real=True)
    r = sp.I * x / 2
    rn = -sp.I * un / 2
    rx = -sp.I * unx / 2
    b_h = (r - sp.conjugate(rn)) ** 2 + rx + sp.conjugate(rx)
    b_t = (r - rn) ** 2 + 2 * rx
    require_zero(b_h + (x - un) ** 2 / 4, "independent ideal H")
    require_zero(b_t + (x + un) ** 2 / 4 + sp.I * unx, "independent ideal T")

    p_h = sp.I * (x + un) * b_h
    p_t = sp.I * (x - un) * b_t
    require_zero(
        p_h + sp.I * (x + un) * (x - un) ** 2 / 4,
        "independent ideal lifted H",
    )
    require_zero(
        p_t + sp.I * (x - un) * (x + un) ** 2 / 4 - unx * (x - un),
        "independent ideal lifted T",
    )
    require_zero(
        p_h.subs(x, -uq) - sp.I * (uq - un) * (uq + un) ** 2 / 4,
        "independent ideal H atom",
    )
    require_zero(
        p_t.subs(x, -uq)
        - sp.I * (uq + un) * (uq - un) ** 2 / 4
        + unx * (uq + un),
        "independent ideal T atom",
    )


def check_lower_endpoint() -> None:
    lam, length = sp.symbols("lam length")
    p0, p1, g0, kappa, alpha, s1, s2, s3 = sp.symbols(
        "p0 p1 g0 kappa alpha s1 s2 s3"
    )
    p = p0 + p1 * lam
    lifted = sp.I * (lam - length) * p
    value = lifted.subs(lam, 0)
    derivative_combo = sp.diff(lifted, lam).subs(lam, 0) + g0 * value
    lower = -kappa * value * s1 + kappa**2 * (
        derivative_combo * s2 + alpha * value * s3
    )
    expected = sp.I * kappa * length * p0 * s1 + sp.I * kappa**2 * (
        (p0 - length * (p1 + g0 * p0)) * s2 - alpha * length * p0 * s3
    )
    require_zero(lower - expected, "independent complete lower endpoint")


def check_recombination() -> None:
    tn, un, t1, u1, interior = sp.symbols("tn un t1 u1 interior")
    original = tn + un - t1 - u1 + interior
    lower = -t1 - u1 + interior
    require_zero(original - tn - un - lower, "independent mode reversal")

    f, o, m, t, tau, exterior = sp.symbols("f o m t tau exterior")
    defect = f - m
    remainder = o - t
    require_zero((remainder + defect + t).subs(o, m - f), "independent R=-D-T")
    require_zero((f + remainder - m + t).subs(o, m - f), "independent F+R")
    require_zero((defect - tau + exterior).subs(f, tau + m - exterior), "independent D=tau-E")
    require_zero(
        (remainder - exterior + tau + t)
        .subs(o, m - f)
        .subs(f, tau + m - exterior),
        "independent exterior remainder",
    )

    transfer = sp.symbols("transfer")
    require_zero((f + transfer) + (o - transfer) - f - o, "independent tie total")
    require_zero(
        ((f + transfer - m) + (o - transfer - t)) - (defect + remainder),
        "independent tie defect/remainder",
    )

    u = sp.symbols("u")
    r_h = -defect
    r_t = -defect + 2 * sp.I * u * t
    require_zero((f + r_h - m).subs(f, m + defect), "independent H mode join")
    require_zero(
        (f + r_t - m - 2 * sp.I * u * t).subs(f, m + defect),
        "independent T mode join",
    )
    require_zero(m + r_h - (m - defect), "independent H carrier join")
    require_zero(
        m + r_t - (m - defect + 2 * sp.I * u * t),
        "independent T carrier join",
    )


def check_abel() -> int:
    checks = 0
    for size in range(2, 34):
        a = [sp.Rational(((7 * j + size) % 23) - 11, 2 * j + 1) for j in range(1, size + 1)]
        k = [sp.Rational(((11 * j + size) % 29) - 14, 3 * j + 2) for j in range(1, size + 1)]
        prefix = []
        running = sp.Integer(0)
        for value in a:
            running += value
            prefix.append(running)
        lhs = sum(a[j] * k[j] for j in range(size))
        rhs = prefix[-1] * k[-1] + sum(
            prefix[j] * (k[j] - k[j + 1]) for j in range(size - 1)
        )
        require_zero(lhs - rhs, f"independent Abel full {size}")
        checks += 1

        kh = k[:-1] + [sp.Integer(0)]
        lhs_h = sum(a[j] * kh[j] for j in range(size))
        rhs_h = sum(
            prefix[j] * (kh[j] - kh[j + 1]) for j in range(size - 1)
        )
        require_zero(lhs_h - rhs_h, f"independent Abel terminal null {size}")
        checks += 1
    return checks


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["jli_25_route", "jli_26_reserve"],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Physical Kernels",
        "B_(H,N)",
        "B_(T,N)",
        "## Lower Endpoint",
        "## Exact Recombination",
        "D_N[P]=F_N[P]-M_N[P]",
        "R_N[P]=-D_N[P]-T_N[P]",
        "## Relative Channels",
        "Dropping D_N would apply Fourier inversion twice",
        "## Abel Form",
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
    require("signed finite-cell-defect estimate open" in artifact.get("status", ""), "status drifted")
    require("no signed finite-cell-defect" in artifact.get("proof_boundary", ""), "proof boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {
            "physical_kernels",
            "lower_endpoint",
            "reverse_recombination",
            "relative_channels",
            "abel",
            "handoff",
        },
        "certificate keys drifted",
    )
    check_sources(artifact)
    substitution_checks = check_physical_factorizations()
    check_ideal_kernels()
    check_lower_endpoint()
    check_recombination()
    abel_checks = check_abel()
    require(abel_checks == 64, "independent Abel count")
    check_rows(artifact)
    check_note(note)
    print(
        "validated joined lower/interior Abel recombination gate: "
        f"{EXPECTED_COUNTS['rows']} rows, "
        f"{EXPECTED_COUNTS['physical_kernel_factorizations']} physical factorizations, "
        f"{EXPECTED_COUNTS['global_recombinations']} global recombinations, "
        f"{EXPECTED_COUNTS['relative_recombinations']} relative joins, "
        f"{substitution_checks} independent physical checks, "
        f"{abel_checks} independent Abel checks, "
        f"{EXPECTED_COUNTS['signed_finite_cell_defect_bounds']} signed defect bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
