#!/usr/bin/env python3
"""Independently validate the cofinal Newman boundary-degree target."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_cofinal_boundary_degree_transfer_target as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.md"
)
WINDING_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_first_jet_winding_gate.json"
)
EXHAUSTION_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json"
)
ADAPTIVE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate.json"
)
TRANSFER_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_forward_sqrt_to_corrected_rs_C1_transfer_gate.json"
)
CORRECTED_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.json"
)
Q207_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate.json"
)

EXPECTED_IDS = [
    "ncbdt_01_normalizer_orientation",
    "ncbdt_02_scaled_contact_charge",
    "ncbdt_03_boundary_rouche",
    "ncbdt_04_zero_degree_exclusion",
    "ncbdt_05_cofinal_contract",
    "ncbdt_06_q207_base",
    "ncbdt_07_boundary_only_improvement",
    "ncbdt_08_ordinary_proxy_contract",
    "ncbdt_09_corrected_proxy_contract",
    "ncbdt_10_proxy_continuity_guard",
    "ncbdt_11_endpoint_safe_model",
    "ncbdt_12_open_xi_boundary_theorem",
]


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.exists():
        issues.append(f"missing {label}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "stored boundary-degree result", issues)
    if not NOTE.exists():
        issues.append("missing rendered note")
    if not stored:
        return issues

    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored boundary-degree payload differs from reconstruction")
    if stored.get("status") != (
        "exact cofinal boundary-degree transfer target and "
        "endpoint-uniformity guard"
    ):
        issues.append("status drifted")
    rows = stored.get("rows", [])
    if len(rows) != 12:
        issues.append("expected 12 target rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_roles = {
        "exact_lemma": 2,
        "exact_composition": 5,
        "exact_reduction": 1,
        "route_decision": 1,
        "nonpromotion_gate": 1,
        "exact_countermodel": 1,
        "open_handoff": 1,
    }
    for role, expected in expected_roles.items():
        observed = sum(row.get("role") == role for row in rows)
        if observed != expected:
            issues.append(
                f"expected {expected} rows with role {role}, got {observed}"
            )

    # Re-derive positive-normalizer orientation with an arbitrary jet scale.
    amplitude, slope, ell = sp.symbols(
        "amplitude slope ell", positive=True
    )
    shear = amplitude * sp.Matrix([[1, 0], [slope / ell, 1]])
    if sp.factor(shear.det()) != amplitude**2:
        issues.append("independent normalizer determinant failed")
    for values in (
        {amplitude: 2, slope: 3, ell: 5},
        {amplitude: sp.Rational(7, 3), slope: 11, ell: 13},
    ):
        if not bool(shear.det().subs(values) > 0):
            issues.append("normalizer orientation is not positive")

    # Re-derive the scaled double-contact Jacobian, allowing variable ell.
    f_x, f_xx, f_xt = sp.symbols("f_x f_xx f_xt", real=True)
    ell_t, ell_x = sp.symbols("ell_t ell_x", real=True)
    jacobian = sp.Matrix(
        [
            [-f_xx, f_x],
            [
                f_xt / ell - f_x * ell_t / ell**2,
                f_xx / ell - f_x * ell_x / ell**2,
            ],
        ]
    )
    contact_det = sp.factor(jacobian.det().subs(f_x, 0))
    if contact_det != -f_xx**2 / ell:
        issues.append("independent scaled contact determinant failed")

    # Check the vector boundary homotopy margin on exact rational fixtures.
    proxy_vectors = [
        sp.Matrix([3, 4]),
        sp.Matrix([sp.Rational(5, 2), -sp.Rational(7, 3)]),
        sp.Matrix([-11, 2]),
    ]
    error_vectors = [
        sp.Matrix([sp.Rational(1, 4), sp.Rational(1, 3)]),
        sp.Matrix([sp.Rational(1, 5), -sp.Rational(1, 7)]),
        sp.Matrix([1, sp.Rational(1, 2)]),
    ]
    s = sp.symbols("s", real=True)
    for proxy, error in zip(proxy_vectors, error_vectors, strict=True):
        proxy_norm_squared = sp.expand(proxy.dot(proxy))
        error_norm_squared = sp.expand(error.dot(error))
        if not bool(error_norm_squared < proxy_norm_squared):
            issues.append("Rouche fixture lacks strict domination")
            continue
        homotopy_norm_squared = sp.expand(
            (proxy + s * error).dot(proxy + s * error)
        )
        minimum = sp.minimum(homotopy_norm_squared, s, sp.Interval(0, 1))
        if not bool(minimum > 0):
            issues.append("Rouche homotopy fixture crossed zero")

    # Reconstruct the endpoint-safe model and its vanishing nonuniform margin.
    t, x = sp.symbols("t x", real=True)
    epsilon, eta = sp.symbols("epsilon eta", positive=True)
    model = x**2 + epsilon - 2 * t
    if sp.simplify(sp.diff(model, t) + sp.diff(model, x, 2)) != 0:
        issues.append("independent endpoint heat equation failed")
    lower = epsilon / 2 + eta
    jet_norm = sp.factor(
        model.subs({t: lower, x: 0}) ** 2
        + sp.diff(model, x).subs({t: lower, x: 0}) ** 2
    )
    if jet_norm != 4 * eta**2:
        issues.append("independent endpoint first-jet norm failed")
    if sp.limit(jet_norm, eta, 0, dir="+") != 0:
        issues.append("endpoint margin does not tend to zero")
    if sp.simplify(
        sp.factor(sp.discriminant(model, x)) - (8 * t - 4 * epsilon)
    ) != 0:
        issues.append("endpoint root discriminant failed")

    # Audit every imported Xi input at its source boundary.
    winding = load_json(WINDING_RESULT, "first-jet winding source", issues)
    exhaustion = load_json(
        EXHAUSTION_RESULT, "diagonal exhaustion source", issues
    )
    adaptive = load_json(ADAPTIVE_RESULT, "adaptive tail source", issues)
    transfer = load_json(TRANSFER_RESULT, "ordinary transfer source", issues)
    corrected = load_json(
        CORRECTED_RESULT, "corrected C1 source", issues
    )
    q207 = load_json(Q207_RESULT, "Q207 source", issues)

    winding_exact = winding.get("exact", {})
    if "+floor(m/2)" not in winding_exact.get("multiplicity_index", ""):
        issues.append("sign-definite multiplicity index source missing")
    if "wind" not in winding_exact.get("global_winding", ""):
        issues.append("global winding source missing")
    exhaustion_exact = exhaustion.get("exact", {})
    if "R_j->infinity" not in exhaustion_exact.get(
        "abstract_exhaustion", ""
    ):
        issues.append("cofinal diagonal exhaustion source missing")
    adaptive_theorem = adaptive.get("theorem", {})
    if "exp(-3h-pi*x/8)/50" not in adaptive_theorem.get(
        "value_error", ""
    ):
        issues.append("adaptive value-error source missing")
    if "exp(-3h-pi*x/8)/50" not in adaptive_theorem.get(
        "first_derivative_error", ""
    ):
        issues.append("adaptive derivative-error source missing")
    transfer_theorem = transfer.get("theorem", {})
    if "x^(23/2)" not in transfer_theorem.get(
        "ordinary_direct_sufficient_target", ""
    ):
        issues.append("ordinary normalized target source missing")
    corrected_exact = corrected.get("exact", {})
    if "32000000*exp(-3L/2)" not in corrected_exact.get(
        "global_c1", ""
    ):
        issues.append("corrected C1 budget source missing")
    q207_summary = q207.get("summary", {})
    if not q207_summary.get("theorem_ready"):
        issues.append("Q207 source is not theorem-ready")
    if q207_summary.get("unresolved_panels") != 0:
        issues.append("Q207 source has unresolved panels")

    exact = stored.get("exact", {})
    cofinal = exact.get("cofinal_contract", {})
    if "j>=208" not in cofinal.get("open_boundary_theorem", ""):
        issues.append("open cofinal boundary index missing")
    if "implies Lambda<=0 and RH" not in cofinal.get(
        "conditional_endgame", ""
    ):
        issues.append("conditional endgame missing")
    continuity = exact.get("continuity_contract", {})
    if "no defined" not in continuity.get("warning", ""):
        issues.append("cutoff continuity warning missing")

    proof_boundary = stored.get("proof_boundary", "")
    for phrase in [
        "does not construct that cofinal proxy family",
        "does not",
        "prove Lambda<=0",
        "prove RH",
    ]:
        if phrase not in proof_boundary:
            issues.append(f"proof boundary missing phrase: {phrase}")

    if NOTE.exists():
        note = NOTE.read_text(encoding="utf-8")
        required_note = [
            "not a proof of `Lambda<=0` or RH",
            "Boundary Rouché Transfer",
            "strictly positive local index",
            "only on the selected high-frequency boundary arcs",
            "has no defined",
            "No endpoint-simplicity burden is inserted.",
            "No finite pilot is to be promoted",
        ]
        for phrase in required_note:
            if phrase not in note:
                issues.append(f"note missing required phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman cofinal boundary-degree transfer target: "
        "12 rows, 0 issues, 2 orientation identities, "
        "1 vector boundary-Rouche theorem, 1 sign-definite degree composition, "
        "1 Q207 base, 2 explicit C1 proxy budgets, "
        "1 endpoint-uniformity countermodel, 1 open cofinal Xi boundary theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
