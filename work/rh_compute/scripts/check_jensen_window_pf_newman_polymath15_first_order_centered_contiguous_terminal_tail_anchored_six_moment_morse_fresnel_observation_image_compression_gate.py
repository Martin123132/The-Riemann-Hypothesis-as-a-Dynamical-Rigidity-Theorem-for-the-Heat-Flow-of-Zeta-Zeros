#!/usr/bin/env python3
"""Validate the Morse-Fresnel observation-image compression gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_observation_image_compression_gate as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"mfic_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "source",
            "gradient",
            "residual",
            "visible",
            "normalize",
            "linearity",
            "functional",
            "base",
            "lift",
            "degree",
            "expanded",
            "plin",
            "collapse",
            "quadratic",
            "kernel",
            "tail",
            "witness",
            "route",
            "target",
            "endpoint",
            "terminal",
            "hermitian",
            "scale",
            "pi",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 25,
    "visible_residual_observations": 8,
    "linear_residual_functionals": 1,
    "maximum_polynomial_degree": 5,
    "quadratic_observation_pairings": 4,
    "tail_projection_templates": 1,
    "universal_endpoint_zeros": 0,
    "enumerated_physical_modes": 0,
    "evaluated_physical_linear_polynomials": 0,
    "coefficient_aware_interior_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
SUMMARY = (
    "validated Morse-Fresnel observation-image compression gate: 25 rows, "
    "8 visible residual observations, 1 degree-five linear functional, "
    "4 quadratic pairings, 1 tail projection template, "
    "0 interior bounds, 0 signed flow bounds"
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


def require_zero(expression: sp.Expr | sp.Matrix, label: str, issues: list[str]) -> None:
    simplified = sp.simplify(expression)
    if isinstance(simplified, sp.MatrixBase):
        if simplified != sp.zeros(*simplified.shape):
            issues.append(label)
    elif simplified != 0:
        issues.append(label)


def close(
    left: mp.mpf | mp.mpc,
    right: mp.mpf | mp.mpc,
    tolerance: str = "1e-50",
) -> bool:
    return abs(left - right) <= mp.mpf(tolerance) * max(1, abs(left), abs(right))


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    saved = payload.get("source_audit", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
            continue
        relative = str(path.relative_to(REPO_ROOT)).replace("\\", "/")
        entry = saved.get(key, {})
        if entry.get("path") != relative:
            issues.append(f"source path drifted: {key}")
        if entry.get("sha256") != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    expected = {
        "endpoint_composition": {
            "endpoint_composed_value_identities": 6,
            "observation_expansions": 1,
        },
        "flow_matrix": {
            "real_flow_observations": 8,
            "flow_rank_upper_bound": 8,
        },
        "morse_fresnel_endpoint": {"exact_transition_integrals": 6},
    }
    for source, checks in expected.items():
        path = gate.SOURCE_PATHS[source]
        if not path.is_file():
            continue
        counts = json.loads(path.read_text(encoding="utf-8")).get("counts", {})
        for key, value in checks.items():
            if counts.get(key) != value:
                issues.append(f"source count drifted: {source}.{key}")


def independent_matrix_audit(issues: list[str]) -> None:
    dimension = 7
    u0 = sp.Matrix(
        dimension,
        4,
        lambda row, column: sp.Rational((row + 3) * (column + 2) + row, row + column + 6),
    )
    u1 = sp.Matrix(
        dimension,
        4,
        lambda row, column: sp.Rational((row + 1) ** 2 + 3 * column + 1, 2 * row + column + 9),
    )
    p = sp.Matrix([sp.Rational(index + 1, index + 4) for index in range(dimension)])
    residual = sp.Matrix(
        [sp.Rational((-1) ** index * (index + 4), index + 8) for index in range(dimension)]
    )
    core = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    terminal = sp.Matrix([3, -2, 5, -7])
    matrix = u0 * core * u1.T + u1 * core * u0.T
    linear = u1 * terminal
    carrier_base = u0.T * p
    carrier_dot = u1.T * p
    error_base = u0.T * residual
    error_dot = u1.T * residual
    gradient = 2 * matrix * p + linear
    factored = u0 * (2 * core * carrier_dot) + u1 * (
        2 * core * carrier_base + terminal
    )
    require_zero(gradient - factored, "independent gradient factorization failed", issues)

    direct = (gradient.T * residual)[0] + (residual.T * matrix * residual)[0]
    compressed = (
        2 * (carrier_base.T * core * error_dot)[0]
        + 2 * (carrier_dot.T * core * error_base)[0]
        + (terminal.T * error_dot)[0]
        + 2 * (error_base.T * core * error_dot)[0]
    )
    require_zero(direct - compressed, "independent residual compression failed", issues)
    require_zero(
        (residual.T * matrix * residual)[0]
        - 2 * (error_base.T * core * error_dot)[0],
        "independent quadratic pairing failed",
        issues,
    )

    null_dimension = 12
    left = sp.Matrix.vstack(sp.eye(4), sp.zeros(null_dimension - 4, 4))
    right = sp.Matrix.vstack(sp.zeros(4, 4), sp.eye(4), sp.zeros(4, 4))
    common = sp.Matrix([0] * 8 + [1, -2, 3, -4])
    if left.T * common != sp.zeros(4, 1) or right.T * common != sp.zeros(4, 1):
        issues.append("independent invisible-kernel witness failed")


def lift(coefficients: list[mp.mpc], ell: mp.mpf) -> list[mp.mpc]:
    padded = coefficients + [mp.mpc(0)] * (6 - len(coefficients))
    result = [-1j * ell * padded[0]]
    result.extend(
        1j * (padded[index - 1] - ell * padded[index])
        for index in range(1, 6)
    )
    return result


def padded(coefficients: list[mp.mpc]) -> list[mp.mpc]:
    return coefficients + [mp.mpc(0)] * (6 - len(coefficients))


def add_scaled(target: list[mp.mpc], source: list[mp.mpc], scale: mp.mpf) -> None:
    for index in range(6):
        target[index] += scale * source[index]


def independent_polynomial_projection_audit(issues: list[str]) -> None:
    mp.mp.dps = 70
    ell = mp.log(17)
    rows = [
        [mp.mpc(1, 0), mp.mpc(2, -1) / 5, mp.mpc(-1, 2) / 7],
        [mp.mpc(2, 1) / 3, mp.mpc(-1, 1) / 5, mp.mpc(3, -2) / 7, mp.mpc(1, 1) / 11, mp.mpc(-2, 1) / 13],
        [mp.mpc(1, -1) / 4, mp.mpc(2, 3) / 7, mp.mpc(-3, 1) / 8, mp.mpc(1, -2) / 9],
        [mp.mpc(-2, 1) / 5, mp.mpc(3, 2) / 8, mp.mpc(1, -3) / 10, mp.mpc(2, 1) / 13],
    ]
    residual_moments = [
        mp.mpc(index + 2, (-1) ** index * (index + 1)) / (11 + index)
        for index in range(6)
    ]
    factor = mp.mpc(7, -3) / 13
    carrier_base = [mp.mpf(2) / 7, -mp.mpf(3) / 11, mp.mpf(5) / 13, -mp.mpf(7) / 17]
    carrier_dot = [-mp.mpf(11) / 19, mp.mpf(13) / 23, -mp.mpf(17) / 29, mp.mpf(19) / 31]
    a_tx, x_t, x_tx, a_t = mp.mpf(3) / 37, -mp.mpf(5) / 41, mp.mpf(7) / 43, -mp.mpf(11) / 47

    base_coefficients = [padded(row) for row in rows]
    lifted_coefficients = [lift(row, ell) for row in rows]

    def projection(coefficients: list[mp.mpc]) -> mp.mpf:
        return mp.re(
            factor
            * mp.fsum(
                coefficients[index] * residual_moments[index]
                for index in range(6)
            )
        )

    base_error = [projection(row) for row in base_coefficients]
    dot_error = [projection(row) for row in lifted_coefficients]
    linear = (
        (carrier_base[1] + a_tx) * dot_error[0]
        + (x_t + carrier_base[0]) * dot_error[1]
        - (carrier_base[3] + x_tx) * dot_error[2]
        - (a_t + carrier_base[2]) * dot_error[3]
        + carrier_dot[1] * base_error[0]
        + carrier_dot[0] * base_error[1]
        - carrier_dot[3] * base_error[2]
        - carrier_dot[2] * base_error[3]
    )

    p_linear = [mp.mpc(0) for _ in range(6)]
    add_scaled(p_linear, base_coefficients[0], carrier_dot[1])
    add_scaled(p_linear, base_coefficients[1], carrier_dot[0])
    add_scaled(p_linear, base_coefficients[2], -carrier_dot[3])
    add_scaled(p_linear, base_coefficients[3], -carrier_dot[2])
    add_scaled(p_linear, lifted_coefficients[0], carrier_base[1] + a_tx)
    add_scaled(p_linear, lifted_coefficients[1], x_t + carrier_base[0])
    add_scaled(p_linear, lifted_coefficients[2], -(carrier_base[3] + x_tx))
    add_scaled(p_linear, lifted_coefficients[3], -(a_t + carrier_base[2]))
    collapsed = projection(p_linear)
    if not close(linear, collapsed, "1e-60"):
        issues.append("independent one-polynomial collapse failed")

    quadratic = (
        base_error[0] * dot_error[1]
        + base_error[1] * dot_error[0]
        - base_error[2] * dot_error[3]
        - base_error[3] * dot_error[2]
    )
    core = mp.matrix(
        [
            [0, mp.mpf(1) / 2, 0, 0],
            [mp.mpf(1) / 2, 0, 0, 0],
            [0, 0, 0, -mp.mpf(1) / 2],
            [0, 0, -mp.mpf(1) / 2, 0],
        ]
    )
    matrix_quadratic = 2 * (mp.matrix(base_error).T * core * mp.matrix(dot_error))[0]
    if not close(quadratic, matrix_quadratic, "1e-60"):
        issues.append("independent four-pair quadratic failed")


def amplitude(j: int, u: mp.mpf, sigma: mp.mpf, t: mp.mpf) -> mp.mpf:
    lam = mp.log(u)
    return lam**j * mp.e ** (t * lam**2 / 4 - sigma * lam)


def endpoint_c(
    alpha: mp.mpf,
    mode: int,
    mu: mp.mpf,
    values: list[mp.mpc],
    sigma: mp.mpf,
    t: mp.mpf,
) -> mp.mpc:
    u0 = alpha / mode
    v = mu / u0
    z = mp.sign(v - 1) * mp.sqrt(2 * (v - 1 - mp.log(v)))
    jacobian = v * z / (v - 1)
    y = mp.sqrt(alpha) * z
    endpoint_amplitude = mp.fsum(
        values[index] * amplitude(index, mu, sigma, t)
        for index in range(len(values))
    )
    saddle_amplitude = mp.fsum(
        values[index] * amplitude(index, u0, sigma, t)
        for index in range(len(values))
    )
    b_endpoint = mp.sqrt(alpha) * endpoint_amplitude * jacobian / mode
    b_zero = mp.sqrt(alpha) * saddle_amplitude / mode
    return (b_endpoint - b_zero) / y


def independent_morse_linearity_audit(issues: list[str]) -> None:
    mp.mp.dps = 70
    alpha = mp.mpf("37.37")
    sigma = mp.mpf("0.501")
    t = mp.mpf(1) / 5000
    coefficients = [
        mp.mpc(2, -1) / 3,
        mp.mpc(-1, 2) / 5,
        mp.mpc(3, 1) / 7,
        mp.mpc(-2, -1) / 11,
        mp.mpc(1, 3) / 13,
        mp.mpc(-3, 2) / 17,
    ]
    for mode in (2, 7, 19, 38, 71):
        for mu in (mp.mpf(1), mp.mpf(5)):
            composed = endpoint_c(alpha, mode, mu, coefficients, sigma, t)
            separated = mp.fsum(
                coefficients[index]
                * endpoint_c(
                    alpha,
                    mode,
                    mu,
                    [mp.mpc(0)] * index + [mp.mpc(1)],
                    sigma,
                    t,
                )
                for index in range(6)
            )
            if not close(composed, separated, "1e-58"):
                issues.append(f"Morse amplitude linearity failed at r={mode}, mu={mu}")


def independent_lift_and_witness_audit(payload: dict, issues: list[str]) -> None:
    lam, ell = sp.symbols("lambda ell", real=True)
    coefficients = sp.symbols("c0:5")
    polynomial = sum(coefficients[index] * lam**index for index in range(5))
    lifted = -sp.I * ell * coefficients[0]
    lifted += sum(
        sp.I * (coefficients[index - 1] - ell * coefficients[index]) * lam**index
        for index in range(1, 5)
    )
    lifted += sp.I * coefficients[4] * lam**5
    require_zero(lifted - sp.I * (lam - ell) * polynomial, "independent polynomial lift failed", issues)

    p_v = 1 + lam**2
    p_n = 1 + lam**4
    p_a = lam**3
    p_q = 1 + lam**3
    witness = sp.expand(
        p_v
        + 2 * p_n
        - 3 * p_a
        - 4 * p_q
        + sp.I * (lam - 2) * (5 * p_v + 6 * p_n - 7 * p_a - 8 * p_q)
    )
    if sp.degree(witness, lam) != 5 or sp.LC(sp.Poly(witness, lam)) != 6 * sp.I:
        issues.append("independent degree-five witness failed")
    if witness.subs(lam, 2) != -21:
        issues.append("independent no-lift-factor witness failed")
    if payload.get("polynomial_certificate", {}).get("witness_polynomial") != str(witness):
        issues.append("saved witness polynomial drifted")


def independent_tail_template_audit(payload: dict, issues: list[str]) -> None:
    expected = [6, 14, 55, 336, 2738, 27936]
    if payload.get("tail_certificate", {}).get("constants") != expected:
        issues.append("tail constants drifted")
    coefficients = [2, -3, 5, -7, 11, -13]
    component_bounds = [2 * value for value in expected]
    direct_triangle = sum(abs(coefficients[index]) * component_bounds[index] for index in range(6))
    template = 2 * sum(abs(coefficients[index]) * expected[index] for index in range(6))
    if direct_triangle != template:
        issues.append("tail projection template arithmetic failed")


def structural_audit(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("kind drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append("row ids drifted")
    if len(ids) != len(set(ids)):
        issues.append("duplicate row ids")
    readiness = [row.get("readiness") for row in rows]
    if readiness.count("proved") != 24 or readiness.count("open") != 1:
        issues.append("row readiness boundary drifted")
    for key in (
        "universal_endpoint_zeros",
        "enumerated_physical_modes",
        "evaluated_physical_linear_polynomials",
        "coefficient_aware_interior_bounds",
        "quadratic_remainder_bounds",
        "signed_flow_bounds",
        "phi_b_bounds",
    ):
        if payload.get("counts", {}).get(key) != 0:
            issues.append(f"{key} was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "does not evaluate P_lin",
        "bound its c-prime interior",
        "signed flow",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    try:
        rebuilt = gate.build_payload()
    except Exception as exc:  # pragma: no cover
        issues.append(f"builder replay failed: {exc}")
    else:
        if payload != rebuilt:
            issues.append("saved result differs from a fresh builder replay")


def note_audit(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Morse-Fresnel Observation-Image Compression Gate",
        "g_p^Tr=Re E[P_lin]",
        "closure through degree five",
        "x_Vz_N+x_Nz_V-x_Az_Q-x_Qz_A",
        "|T[P]|<2h^2",
        "0 coefficient-aware interior bounds",
        "this is not a proof",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
    ):
        if marker not in text:
            issues.append(f"note marker missing: {marker}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        independent_source_audit(payload, issues)
        independent_matrix_audit(issues)
        independent_polynomial_projection_audit(issues)
        independent_morse_linearity_audit(issues)
        independent_lift_and_witness_audit(payload, issues)
        independent_tail_template_audit(payload, issues)
        structural_audit(payload, issues)
    note_audit(issues)
    if issues:
        for issue in issues:
            print(f"FAIL: {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
