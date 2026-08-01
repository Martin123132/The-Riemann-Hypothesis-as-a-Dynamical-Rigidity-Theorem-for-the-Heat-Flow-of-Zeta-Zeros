#!/usr/bin/env python3
"""Validate the physical q=1 saddle-phase and phase-variance reduction."""

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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_q1_saddle_phase_variance_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"spv_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "sigma",
            "frequency",
            "lock",
            "profile",
            "monotone",
            "anchor",
            "distance",
            "transfer",
            "old_guard",
            "two",
            "variance",
            "discriminant",
            "guard",
            "full",
            "target",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 16,
    "physical_q1_parameter_laws": 4,
    "common_phase_laws": 1,
    "amplitude_monotonicity_theorems": 1,
    "saddle_phase_transfer_bounds": 5,
    "two_carrier_order_signs": 1,
    "phase_variance_identities": 1,
    "ordinary_fibre_discriminants": 1,
    "ordered_profile_guards": 1,
    "signed_physical_targets": 1,
    "signed_physical_bounds": 0,
    "phi_b_bounds": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated physical q=1 saddle-phase variance reduction: 16 rows, "
    "4 parameter laws, 5 saddle transfers, 1 two-carrier sign, "
    "1 phase discriminant, 1 ordered-profile guard, "
    "0 signed physical bounds"
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


def require_zero(expression: sp.Expr, label: str, issues: list[str]) -> None:
    if sp.simplify(sp.expand(expression)) != 0:
        issues.append(label)


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")
    try:
        sources = {
            key: json.loads(path.read_text(encoding="utf-8"))
            for key, path in gate.SOURCE_PATHS.items()
        }
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"independent source load failed: {exc}")
        return
    frame = sources["critical_frame"].get("exact", {}).get(
        "critical_frame", {}
    )
    if "Re(alpha)=L/2" not in frame.get("alpha_real", ""):
        issues.append("independent alpha-real provenance failed")
    if "Im(alpha)=-(1/2)atan(x)" not in frame.get("alpha_imag", ""):
        issues.append("independent alpha-imag provenance failed")
    reciprocal = sources["reciprocal_saddle"].get("exact", {})
    if "T_0-omega" not in reciprocal.get("phase_saddle_comparison", ""):
        issues.append("independent saddle-frequency provenance failed")
    anchor = sources["phase_anchor"].get("exact", {})
    if "eta=f_1/|f_1|" not in anchor.get("unit_anchor", ""):
        issues.append("independent phase-anchor provenance failed")
    projection = sources["direct_projection"].get("exact", {})
    if "B_0=beta/(|M_t(s)|*|f_1|)>0" not in projection.get(
        "endpoint_projection", ""
    ):
        issues.append("independent B_0 provenance failed")
    endpoint = sources["endpoint_centering"].get("exact", {})
    if "kappa=(-1)^N B_0 real" not in endpoint.get("endpoint_factor", ""):
        issues.append("independent endpoint parity provenance failed")
    turan = sources["turan"]
    if turan.get("counts", {}).get("signed_joint_bounds") != 0:
        issues.append("independent upstream proof-boundary audit failed")


def independent_q1_parameter_audit(issues: list[str]) -> None:
    L, x = sp.symbols("L x", positive=True, real=True)
    t = 1 / (2 * L**2)
    alpha_r = L / 2 + sp.log(1 + x**-2) / 4 - 1 / (1 + x**2)
    alpha_i = -sp.atan(x) / 2 + 3 * x / (1 + x**2)
    sigma = sp.Rational(1, 2) + t * alpha_r / 2
    expected_sigma = (
        sp.Rational(1, 2)
        + 1 / (8 * L)
        + sp.log(1 + x**-2) / (16 * L**2)
        - 1 / (4 * L**2 * (1 + x**2))
    )
    require_zero(sigma - expected_sigma, "independent sigma identity failed", issues)

    omega = x / 2 - t * alpha_i / 2
    expected_omega = (
        x / 2
        + sp.atan(x) / (8 * L**2)
        - 3 * x / (4 * L**2 * (1 + x**2))
    )
    require_zero(
        omega - expected_omega,
        "independent frequency identity failed",
        issues,
    )
    t0 = x / 2 + sp.pi / (16 * L**2)
    epsilon = (
        sp.atan(1 / x) / (8 * L**2)
        + 3 * x / (4 * L**2 * (1 + x**2))
    )
    reciprocal_omega = expected_omega.xreplace(
        {sp.atan(x): sp.pi / 2 - sp.atan(1 / x)}
    )
    require_zero(
        reciprocal_omega - (t0 - epsilon),
        "independent saddle lock failed",
        issues,
    )

    z = sp.pi / (8 * L**2 * x)
    delta = (
        1 / (1 + x**2)
        - sp.log(1 + x**-2) / 4
        + sp.log(1 + z) / 2
    )
    log_a = L / 2 + sp.log(1 + z) / 2
    require_zero(
        delta - (log_a - alpha_r),
        "independent amplitude defect failed",
        issues,
    )
    b_a = sp.Rational(1, 2) - delta / (4 * L**2)
    u = sp.symbols("u", nonnegative=True, real=True)
    derivative = sp.diff(b_a * u + u**2 / (8 * L**2), u)
    require_zero(
        derivative - b_a - u / (4 * L**2),
        "independent amplitude derivative failed",
        issues,
    )

    # Direct high-precision checks at the weakest stated endpoint L=50.
    l_value = sp.Integer(50)
    x_value = 4 * sp.pi * sp.exp(l_value)
    epsilon_value = sp.N(epsilon.subs({L: l_value, x: x_value}), 80)
    epsilon_lower = sp.N(
        (3 / (8 * L**2 * x)).subs({L: l_value, x: x_value}), 80
    )
    epsilon_upper = sp.N(
        (7 / (8 * L**2 * x)).subs({L: l_value, x: x_value}), 80
    )
    if not epsilon_lower < epsilon_value < epsilon_upper:
        issues.append("independent epsilon endpoint bounds failed")


def independent_transfer_audit(issues: list[str]) -> None:
    distances = (
        sp.Rational(1, 7),
        sp.Rational(2, 5),
        sp.Rational(3, 4),
        sp.Rational(5, 3),
    )
    amplitudes = (
        sp.Rational(7, 5),
        sp.Rational(11, 7),
        sp.Rational(13, 9),
        sp.Rational(17, 11),
    )
    omega = sp.Rational(17, 13)
    t0 = sp.Rational(19, 14)
    epsilon = abs(t0 - omega)
    for order in range(5):
        h_omega = sum(
            amplitude * distance**order * sp.exp(-sp.I * omega * distance)
            for amplitude, distance in zip(amplitudes, distances)
        )
        h_t0 = sum(
            amplitude * distance**order * sp.exp(-sp.I * t0 * distance)
            for amplitude, distance in zip(amplitudes, distances)
        )
        bound = epsilon * sum(
            amplitude * distance ** (order + 1)
            for amplitude, distance in zip(amplitudes, distances)
        )
        difference = abs(complex(sp.N(h_omega - h_t0, 60)))
        bound_value = float(sp.N(bound, 60))
        if difference > bound_value + 1e-40:
            issues.append(f"independent order-{order} saddle transfer failed")


def independent_two_carrier_audit(issues: list[str]) -> None:
    phase_vectors = (
        (sp.Integer(1), sp.Integer(0)),
        (-sp.Integer(1), sp.Integer(0)),
        (sp.Integer(0), sp.Integer(1)),
        (sp.Integer(0), -sp.Integer(1)),
        (sp.Rational(3, 5), sp.Rational(4, 5)),
        (-sp.Rational(5, 13), sp.Rational(12, 13)),
    )
    parameter_sets = (
        (sp.Integer(5), sp.Integer(3), sp.Integer(4), sp.Integer(2)),
        (sp.Rational(7, 3), sp.Rational(11, 6), sp.Rational(5, 2), sp.Rational(4, 3)),
    )
    for a1, a2, u1, u2 in parameter_sets:
        upper = (a2 - a1) * (a1 * u1**2 - a2 * u2**2)
        if upper >= 0:
            issues.append("independent two-carrier parameter ordering failed")
            continue
        for c1, s1 in phase_vectors:
            for c2, s2 in phase_vectors:
                f = a1 * c1 + a2 * c2
                fp = a1 * u1 * s1 + a2 * u2 * s2
                fpp = -(a1 * u1**2 * c1 + a2 * u2**2 * c2)
                four_t = sp.simplify(f * fpp - fp**2)
                if four_t > upper:
                    issues.append("independent two-carrier upper bound failed")
                    return


def independent_phase_variance_audit(issues: list[str]) -> None:
    i = sp.I
    distances = (
        sp.Rational(1, 5),
        sp.Rational(3, 7),
        sp.Rational(7, 9),
        sp.Rational(11, 8),
    )
    weights = (
        1 + 2 * i / 3,
        -sp.Rational(3, 5) + i / 7,
        sp.Rational(5, 8) - 4 * i / 11,
        -sp.Rational(7, 13) - 2 * i / 9,
    )
    u_x = sp.Rational(5, 17)
    moments = tuple(
        sum(distance**order * weight for distance, weight in zip(distances, weights))
        for order in range(3)
    )
    w, m1, m2 = moments
    if sp.simplify(w) == 0:
        issues.append("independent phase-variance ordinary fibre vanished")
        return
    f = sp.re(w)
    g = sp.im(w)
    direct_four_current = sp.expand(
        -f * sp.re(m2) - sp.im(m1) ** 2 + 2 * u_x * f * g
    )
    radius2 = sp.expand(w * sp.conjugate(w))
    a_num = sp.re(m1 * sp.conjugate(w))
    c_num = (m2 * w - m1**2) * sp.conjugate(w) ** 2
    d_hat = sp.expand_complex(
        f**2 * sp.re(c_num)
        + a_num**2 * radius2
        - (2 * u_x * radius2**2 + sp.im(c_num)) * f * g
    )
    require_zero(
        d_hat + direct_four_current * radius2**2,
        "independent branch-free discriminant failed",
        issues,
    )

    mu1 = sp.simplify(m1 / w)
    mu2 = sp.simplify(m2 / w)
    variance = sp.simplify(mu2 - mu1**2)
    d_phase_times_f2 = sp.expand_complex(
        f**2 * sp.re(variance)
        + sp.re(mu1) ** 2 * radius2
        - (2 * u_x + sp.im(variance)) * f * g
    )
    require_zero(
        d_phase_times_f2 * radius2**2 - d_hat,
        "independent phase and branch-free discriminants disagree",
        issues,
    )


def independent_ordered_guard(issues: list[str]) -> None:
    d = sp.log(2)
    amplitudes = (sp.Integer(4), 2 * sp.sqrt(2), sp.sqrt(2))
    distances = (4 * d, 3 * d, d)
    weights = (-amplitudes[0], amplitudes[1], amplitudes[2])
    f = sp.simplify(sum(weights))
    fp = sp.Integer(0)
    fpp = -sp.simplify(
        sum(distance**2 * weight for distance, weight in zip(distances, weights))
    )
    current = sp.simplify((f * fpp - fp**2) / 4)
    expected = (-185 + 134 * sp.sqrt(2)) * d**2 / 2
    require_zero(
        current - expected,
        "independent limiting ordered-profile guard failed",
        issues,
    )
    if sp.N(current, 60) <= 0:
        issues.append("independent limiting guard sign failed")

    exp_x = sp.Rational(3381, 10000)
    exp_lower = 1 + exp_x + exp_x**2 / 2 + exp_x**3 / 6
    if exp_lower <= sp.Rational(7, 5):
        issues.append("independent r>7/5 certificate failed")
    f_lower = (
        sp.Rational(7071, 10000)
        * (1 - sp.Rational(343, 2000000))
        + sp.Rational(3535, 10000)
        * (1 - sp.Rational(735, 2000000))
        - 1
    )
    curvature_lower = 16 - 9 * sp.Rational(5, 7) - sp.Rational(125, 343)
    if f_lower <= sp.Rational(3, 50):
        issues.append("independent finite profile value bound failed")
    if curvature_lower <= 9:
        issues.append("independent finite profile curvature bound failed")

    # Directly test all corners of the larger certified (B,c) rectangle.
    for b_value in (sp.Rational(49, 100), sp.Rational(1, 2)):
        for c_value in (sp.Integer(0), sp.Rational(1, 20000)):
            actual = tuple(
                sp.exp(b_value * distance + c_value * distance**2)
                for distance in distances
            )
            f_value = -actual[0] + actual[1] + actual[2]
            fpp_value = (
                16 * actual[0] - 9 * actual[1] - actual[2]
            ) * d**2
            if sp.N(f_value, 60) <= 0 or sp.N(fpp_value, 60) <= 0:
                issues.append("independent finite profile corner sign failed")
                return


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Physical q=1 Saddle Phase And Variance Reduction",
        "Omega=-Im(s_*)",
        "## Phase Variance",
        "## Stronger Nonpromotion Guard",
        "0 signed physical bounds",
        "This is not a proof of RH",
        f"python work/rh_compute/scripts/check_{STEM}.py",
    ):
        if marker not in text:
            issues.append(f"note missing marker: {marker}")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-08-01":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    if sum(row.get("readiness") == "open" for row in rows) != 1:
        issues.append("open-row count drifted")
    if sum(row.get("readiness") == "guard_validated" for row in rows) != 1:
        issues.append("guard-row count drifted")

    status = payload.get("status", "")
    for token in (
        "physical q=1 saddle parameter",
        "phase-variance",
        "physical joint bound open",
        "no Phi_B bound",
        "or RH",
    ):
        if token not in status:
            issues.append(f"status missing token: {token}")
    boundary = payload.get("proof_boundary", "")
    for token in (
        "does not prove D_hat>=0",
        "bound E or R_corr",
        "upper-bound Phi_B",
        "Lambda<=0",
        "RH",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing token: {token}")

    independent_source_audit(payload, issues)
    try:
        rebuilt = gate.build_payload()
        if payload != rebuilt:
            issues.append("full builder recomputation drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"builder recomputation failed: {exc}")

    counts = payload.get("counts", {})
    if counts.get("signed_physical_bounds") != 0:
        issues.append("unexpected signed physical bound")
    if counts.get("phi_b_bounds") != 0:
        issues.append("unexpected Phi_B bound")
    if counts.get("xi_level_current_theorems") != 0:
        issues.append("unexpected Xi-level theorem")

    independent_q1_parameter_audit(issues)
    independent_transfer_audit(issues)
    independent_two_carrier_audit(issues)
    independent_phase_variance_audit(issues)
    independent_ordered_guard(issues)
    validate_note(issues)


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        validate(payload, issues)
    if issues:
        print("physical q=1 saddle-phase variance validation failed:")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
