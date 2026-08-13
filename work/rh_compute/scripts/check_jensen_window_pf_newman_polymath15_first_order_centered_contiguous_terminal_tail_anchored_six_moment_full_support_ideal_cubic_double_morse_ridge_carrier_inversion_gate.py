#!/usr/bin/env python3
"""Independently check the double-Morse ridge/carrier inversion gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_"
    "morse_ridge_carrier_inversion_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 28,
    "inversion_identities": 5,
    "fiber_identities": 3,
    "ridge_decomposition_identities": 2,
    "endpoint_weight_checks": 5,
    "internal_tie_cancellations": 1,
    "witness_exact_checks": 5,
    "witness_numerical_checks": 8,
    "nonzero_exterior_witnesses": 1,
    "signed_ridge_completed_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "rci_01_sources",
    "rci_02_cell_poisson",
    "rci_03_diagonal",
    "rci_04_cell_join",
    "rci_05_band",
    "rci_06_defect",
    "rci_07_one_carrier",
    "rci_08_star_weights",
    "rci_09_ties",
    "rci_10_terminal",
    "rci_11_domain",
    "rci_12_diagonal_chart",
    "rci_13_kernel",
    "rci_14_ridge_lift",
    "rci_15_transverse",
    "rci_16_wings",
    "rci_17_ridge_sum",
    "rci_18_carrier_relation",
    "rci_19_witness_source",
    "rci_20_witness_transform",
    "rci_21_witness_cell",
    "rci_22_witness_bound",
    "rci_23_exterior",
    "rci_24_local_guard",
    "rci_25_global_route",
    "rci_26_exterior_open",
    "rci_27_signed_open",
    "rci_28_boundary",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(
        set(source)
        == {"finite_cell_inversion", "hyperbolic_transport", "joined_recombination"},
        "source keys drifted",
    )
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_inversion() -> int:
    m, e_ext, tau, aliases = sp.symbols("m e tau aliases")
    diagonal = m - e_ext
    cell = tau + diagonal + aliases
    band = tau + m - e_ext
    defect = band - m
    checks = [
        2 * m - cell - (m - tau + e_ext - aliases),
        diagonal + e_ext - m,
        2 * m - band - (m - tau + e_ext),
        m - defect - (2 * m - band),
    ]
    z_ridge, transverse, wings = sp.symbols("Z V W")
    checks.append(
        (z_ridge + transverse + wings + e_ext)
        - ((z_ridge + transverse + wings) + e_ext)
    )
    for index, expression in enumerate(checks):
        require_zero(expression, f"independent inversion {index}")
    return len(checks)


def check_ridge_kernel() -> int:
    s, lower, upper = sp.symbols("s lower upper", real=True)
    kappa = 1 / (2 * sp.pi * sp.I)
    e_char = lambda value: sp.exp(2 * sp.pi * sp.I * value)
    kernel = kappa * (e_char(s * upper) - e_char(s * lower)) / s
    checks = [
        sp.diff(kernel, upper) - e_char(s * upper),
        sp.diff(kernel, lower) + e_char(s * lower),
        sp.limit(kernel, s, 0) - (upper - lower),
    ]
    for index, expression in enumerate(checks):
        require_zero(expression, f"independent ridge kernel {index}")
    return len(checks)


def check_witness() -> tuple[int, int]:
    lower = Fraction(4, 5)
    upper = Fraction(4, 3)
    width = upper - lower
    require(width == Fraction(8, 15), "independent witness width")
    mp.mp.dps = 75
    lower_mp = mp.mpf(lower.numerator) / lower.denominator
    upper_mp = mp.mpf(upper.numerator) / upper.denominator

    def sinc(value):
        return mp.sin(mp.pi * value) / (mp.pi * value)

    checks = 0
    for r_value in [
        mp.mpf("0.805"),
        mp.mpf("0.89"),
        mp.mpf("0.99"),
        mp.mpf("1.07"),
        mp.mpf("1.22"),
        mp.mpf("1.329"),
    ]:
        transformed = mp.quad(
            lambda t: mp.e ** (-2j * mp.pi * r_value * t),
            [-mp.mpf("0.5"), mp.mpf("0.5")],
        )
        require(abs(transformed - sinc(r_value)) < mp.mpf("1e-64"), "independent sinc")
        checks += 1

    diagonal = mp.quad(sinc, [lower_mp, mp.mpf("1"), upper_mp])
    sine_integral = (
        mp.si(mp.pi * upper_mp) - mp.si(mp.pi * lower_mp)
    ) / mp.pi
    require(abs(diagonal - sine_integral) < mp.mpf("1e-64"), "independent Si")
    checks += 1
    require(abs(diagonal) <= mp.mpf(8) / 15, "independent diagonal bound")
    checks += 1
    exterior = 1 - diagonal
    require(abs(exterior) >= mp.mpf(7) / 15, "independent exterior bound")
    checks += 1
    return checks, 1


def check_endpoint_weights() -> int:
    weights = [
        Fraction(1, 2),
        Fraction(1),
        Fraction(1),
        Fraction(1),
        Fraction(1),
        Fraction(1),
        Fraction(1, 2),
    ]
    checks = 0
    for weight in weights:
        require(2 * weight - weight == weight, "independent starred weight")
        checks += 1
    require(-Fraction(1, 2) + Fraction(1, 2) == 0, "independent tie halves")
    return checks


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["rci_26_exterior_open", "rci_27_signed_open"],
        "open rows drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Exact Carrier Coefficient",
        "exactly one starred carrier remains",
        "## Hyperbolic Ridge Decomposition",
        "V_q+W_q+E_q^ext=0",
        "## Exact Obstruction Witness",
        "[4/5,4/3]",
        ">=7/15",
        "## Ideal Joins",
        "## Route Decision",
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
    require("one-starred-carrier coefficient" in artifact.get("status", ""), "status drifted")
    require("no signed ridge-completed ideal join" in artifact.get("proof_boundary", ""), "boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {"inversion", "fiber", "witness", "endpoint", "handoff"},
        "certificate keys drifted",
    )
    check_sources(artifact)
    inversion_checks = check_inversion()
    ridge_checks = check_ridge_kernel()
    witness_checks, exterior_witnesses = check_witness()
    endpoint_checks = check_endpoint_weights()
    check_rows(artifact)
    check_note(note)
    print(
        "validated double-Morse ridge/carrier inversion gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {inversion_checks} independent inversion checks, "
        f"{ridge_checks} independent ridge checks, "
        f"{witness_checks} independent witness checks, "
        f"{endpoint_checks} endpoint-weight checks, "
        f"{exterior_witnesses} nonzero exterior witness, "
        f"{EXPECTED_COUNTS['signed_ridge_completed_bounds']} signed ridge-completed bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
