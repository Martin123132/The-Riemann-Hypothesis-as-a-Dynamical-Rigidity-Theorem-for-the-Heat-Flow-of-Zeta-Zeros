#!/usr/bin/env python3
"""Validate the six-moment flow-matrix and reciprocal-phase reduction."""

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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_flow_matrix_phase_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"fmpr_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "shift",
            "lift",
            "observations",
            "flow",
            "matrix",
            "factor",
            "rank",
            "inertia",
            "fibre",
            "normalize",
            "rotation",
            "pi",
            "poisson",
            "endpoint",
            "hermitian",
            "null",
            "transpose",
            "window",
            "handoff",
            "guard",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 22,
    "correction_free_moments": 6,
    "real_moment_coordinates": 12,
    "real_flow_observations": 8,
    "flow_rank_upper_bound": 8,
    "generic_flow_rank": 8,
    "generic_positive_inertia": 4,
    "generic_negative_inertia": 4,
    "generic_zero_inertia": 4,
    "reciprocal_phase_families": 3,
    "hermitian_reciprocal_diagonal_nulls": 1,
    "stationary_phase_remainder_bounds": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
SUMMARY = (
    "validated six-moment flow matrix/phase reduction: 22 rows, "
    "rank <=8, generic inertia (4,4,4), 3 reciprocal phase families, "
    "1 Hermitian diagonal null, 0 signed flow bounds"
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


def require_zero(
    expression: sp.Expr | sp.Matrix,
    label: str,
    issues: list[str],
) -> None:
    if isinstance(expression, sp.MatrixBase):
        if any(sp.simplify(sp.expand(entry)) != 0 for entry in expression):
            issues.append(label)
        return
    if sp.simplify(sp.expand(expression)) != 0:
        issues.append(label)


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    saved = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif saved.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")
    try:
        sources = {
            key: json.loads(path.read_text(encoding="utf-8"))
            for key, path in gate.SOURCE_PATHS.items()
        }
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"independent source load failed: {exc}")
        return
    if sources["five_moment_quadratic"].get("counts", {}).get(
        "quadratic_rank_upper_bound"
    ) != 4:
        issues.append("independent rank-four source audit failed")
    if sources["six_moment_flow"].get("counts", {}).get(
        "max_moment_order"
    ) != 5:
        issues.append("independent six-moment source audit failed")
    if "standard Poisson phase" not in sources["reciprocal_saddle"].get(
        "exact", {}
    ).get("pi_provenance", ""):
        issues.append("independent reciprocal-pi source audit failed")


def independent_shift_audit(issues: list[str]) -> None:
    log_a = sp.symbols("A", real=True)
    real = sp.symbols("x0:6", real=True)
    imag = sp.symbols("y0:6", real=True)
    coefficient_real = sp.symbols("p0:5", real=True)
    coefficient_imag = sp.symbols("q0:5", real=True)
    moments = [real[index] + sp.I * imag[index] for index in range(6)]
    coefficients = [
        coefficient_real[index] + sp.I * coefficient_imag[index]
        for index in range(5)
    ]
    direct = sum(
        coefficients[index]
        * sp.I
        * (moments[index + 1] - log_a * moments[index])
        for index in range(5)
    )
    lifted_coefficients = [-sp.I * log_a * coefficients[0]]
    lifted_coefficients.extend(
        sp.I * (coefficients[index - 1] - log_a * coefficients[index])
        for index in range(1, 5)
    )
    lifted_coefficients.append(sp.I * coefficients[4])
    lifted = sum(
        lifted_coefficients[index] * moments[index] for index in range(6)
    )
    require_zero(direct - lifted, "independent complex shift failed", issues)


def independent_matrix_audit(issues: list[str]) -> None:
    dimension = 7
    u0 = sp.Matrix(
        dimension,
        4,
        lambda row, column: sp.Rational(
            (row + 2) * (column + 3) + (row == column),
            row + column + 5,
        ),
    )
    u1 = sp.Matrix(
        dimension,
        4,
        lambda row, column: sp.Rational(
            (row + 1) ** 2 + 2 * column + 1,
            2 * row + column + 7,
        ),
    )
    state = sp.Matrix(
        [sp.Rational((-1) ** index * (index + 2), index + 3) for index in range(dimension)]
    )
    core = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    tail = sp.Matrix([2, 3, -5, -7])
    matrix = u0 * core * u1.T + u1 * core * u0.T
    linear = u1 * tail
    observation = u0.T * state
    derivative = u1.T * state
    lhs = (state.T * matrix * state)[0] + (linear.T * state)[0]
    rhs = 2 * (observation.T * core * derivative)[0] + (
        tail.T * derivative
    )[0]
    require_zero(lhs - rhs, "independent flow matrix failed", issues)
    if matrix != matrix.T:
        issues.append("independent flow symmetry failed")
    factor = sp.Matrix.hstack(u0, u1)
    zero = sp.zeros(4)
    doubled_core = sp.Matrix.vstack(
        sp.Matrix.hstack(zero, core),
        sp.Matrix.hstack(core, zero),
    )
    require_zero(
        matrix - factor * doubled_core * factor.T,
        "independent rank-eight factorization failed",
        issues,
    )


def _real_row(vector: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        [sp.expand(entry).as_real_imag()[0] for entry in vector]
        + [-sp.expand(entry).as_real_imag()[1] for entry in vector]
    )


def independent_rank_witness(issues: list[str]) -> None:
    i = sp.I
    log_n = sp.Integer(2)
    log_a = sp.Integer(3)
    rho_1, rho_2 = 1 + i, 2 - i
    rho_1_x, rho_2_x = 1 - 2 * i, -1 + i
    s_prime, s_second = 1 - i / 2, 2 + i
    kappa = 1 + i
    b, b_x, u_n, u_n_x = -sp.Rational(1, 2), 1, 1, 1
    v0 = sp.Matrix([1, rho_1, rho_2, 0, 0])
    v1 = sp.Matrix(
        [
            log_n,
            log_n * rho_1 - 1,
            log_n * rho_2 - rho_1,
            -rho_2,
            0,
        ]
    )
    v2 = sp.Matrix(
        [
            log_n**2,
            log_n**2 * rho_1 - 2 * log_n,
            log_n**2 * rho_2 - 2 * log_n * rho_1 + 1,
            rho_1 - 2 * log_n * rho_2,
            rho_2,
        ]
    )
    e0 = sp.Matrix([0, rho_1_x, rho_2_x, 0, 0])
    e1 = sp.Matrix(
        [0, log_n * rho_1_x, log_n * rho_2_x - rho_1_x, -rho_2_x, 0]
    )
    q0 = kappa * v0 + s_prime * v1 + e0
    q1 = kappa * v1 + s_prime * v2 + e1
    a0 = s_prime * v1 + i * b * u_n * v0
    n0 = (
        s_second * v1
        + s_prime * q1
        + i * (b_x * u_n + b * u_n_x) * v0
        + i * b * u_n * q0
    )
    base = sp.Matrix.hstack(
        _real_row(v0), _real_row(n0), _real_row(a0), _real_row(q0)
    )
    embedding = sp.zeros(10, 12)
    shift = sp.zeros(10, 12)
    for index in range(5):
        embedding[index, index] = 1
        embedding[5 + index, 6 + index] = 1
        shift[index, index + 1] = 1
        shift[5 + index, 6 + index + 1] = 1
    identity = sp.eye(5)
    zero = sp.zeros(5)
    complex_structure = sp.Matrix.vstack(
        sp.Matrix.hstack(zero, -identity),
        sp.Matrix.hstack(identity, zero),
    )
    derivative = complex_structure * (shift - log_a * embedding)
    factor = sp.Matrix.hstack(embedding.T * base, derivative.T * base)
    minor = sp.factor(factor[list(range(8)), :].det())
    if factor.rank() != 8 or minor != sp.Rational(51095, 16):
        issues.append("independent rank-eight witness failed")
    saved = load_json(RESULT, issues).get("rank_certificate", {})
    if saved.get("witness_minor") != "51095/16":
        issues.append("saved rank witness drifted")
    if saved.get("flow_inertia") != [4, 4, 4]:
        issues.append("saved flow inertia drifted")


def independent_common_phase_audit(issues: list[str]) -> None:
    matrix = sp.Matrix(
        [
            [2, 3, 5, 7],
            [3, 11, 13, 17],
            [5, 13, 19, 23],
            [7, 17, 23, 29],
        ]
    )
    identity = sp.eye(2)
    zero = sp.zeros(2)
    j = sp.Matrix.vstack(
        sp.Matrix.hstack(zero, -identity),
        sp.Matrix.hstack(identity, zero),
    )
    hermitian = (matrix - j * matrix * j) / 2
    transpose = (matrix + j * matrix * j) / 2
    require_zero(hermitian * j - j * hermitian, "independent Hermitian projection failed", issues)
    require_zero(transpose * j + j * transpose, "independent transpose projection failed", issues)
    require_zero(
        j * matrix - matrix * j + 2 * transpose * j,
        "independent common-rotation cancellation failed",
        issues,
    )


def independent_reciprocal_phase_audit(issues: list[str]) -> None:
    a, n, m, k, r = sp.symbols("a n m k r", positive=True)
    n_star, m_star = a**2 / k, a**2 / r
    endpoint = a**2 * sp.log(n) - k * n
    hermitian = a**2 * sp.log(n) - a**2 * sp.log(m) - k * n + r * m
    transpose = a**2 * sp.log(n) + a**2 * sp.log(m) - k * n - r * m
    for label, expression, variable, saddle in (
        ("endpoint", endpoint, n, n_star),
        ("Hermitian n", hermitian, n, n_star),
        ("Hermitian m", hermitian, m, m_star),
        ("transpose n", transpose, n, n_star),
        ("transpose m", transpose, m, m_star),
    ):
        require_zero(
            sp.diff(expression, variable).subs(variable, saddle),
            f"independent {label} saddle failed",
            issues,
        )
    phase = sp.expand_log(
        hermitian.subs({n: n_star, m: m_star}), force=True
    )
    multiplier = sp.expand_log(sp.log(n_star / m_star), force=True)
    require_zero(
        phase - a**2 * multiplier,
        "independent Hermitian phase lock failed",
        issues,
    )
    require_zero(
        phase.subs(r, k),
        "independent Hermitian diagonal phase failed",
        issues,
    )
    require_zero(
        multiplier.subs(r, k),
        "independent Hermitian diagonal multiplier failed",
        issues,
    )


def independent_finite_carrier_audit(issues: list[str]) -> None:
    mp.mp.dps = 70
    log_a = mp.log(17)
    logarithms = [mp.log(value) for value in (2, 3, 5, 7, 11)]
    amplitudes = [mp.mpf(value) / 13 for value in (17, 19, 23, 29, 31)]
    rows = [
        [mp.mpc(j + 2, 2 * j - 1) / (11 + j) for j in range(5)],
        [mp.mpc(3 * j - 2, j + 1) / (17 + j) for j in range(5)],
        [mp.mpc(j**2 + 1, 4 - j) / (23 + j) for j in range(5)],
        [mp.mpc(5 - j, 2 * j + 3) / (29 + j) for j in range(5)],
    ]
    x_t, a_tx, a_t, x_tx = (
        mp.mpf(2) / 7,
        -mp.mpf(3) / 11,
        mp.mpf(5) / 13,
        -mp.mpf(7) / 17,
    )
    xi = mp.mpf(19) / 7

    def moments(value: mp.mpf) -> list[mp.mpc]:
        carriers = [
            amplitudes[index]
            * mp.e ** (-mp.j * value * (log_a - logarithms[index]))
            for index in range(5)
        ]
        return [
            sum(
                logarithms[index] ** order * carriers[index]
                for index in range(5)
            )
            for order in range(6)
        ]

    def observation(row: list[mp.mpc], values: list[mp.mpc]) -> mp.mpf:
        return mp.re(sum(row[index] * values[index] for index in range(5)))

    def derivative_observation(
        row: list[mp.mpc], values: list[mp.mpc]
    ) -> mp.mpf:
        coefficients = [-mp.j * log_a * row[0]]
        coefficients.extend(
            mp.j * (row[index - 1] - log_a * row[index])
            for index in range(1, 5)
        )
        coefficients.append(mp.j * row[4])
        return mp.re(sum(coefficients[index] * values[index] for index in range(6)))

    def psi(value: mp.mpf) -> mp.mpf:
        values = moments(value)
        v, n_value, a_value, q = [observation(row, values) for row in rows]
        return (
            v * n_value
            - a_value * q
            + x_t * n_value
            + a_tx * v
            - a_t * q
            - x_tx * a_value
        )

    values = moments(xi)
    v, n_value, a_value, q = [observation(row, values) for row in rows]
    dv, dn, da, dq = [derivative_observation(row, values) for row in rows]
    exact = (
        dv * n_value
        + v * dn
        - da * q
        - a_value * dq
        + x_t * dn
        + a_tx * dv
        - a_t * dq
        - x_tx * da
    )
    numerical = mp.diff(psi, xi)
    if abs(exact - numerical) > mp.mpf("1e-60"):
        issues.append("independent finite-carrier flow audit failed")


def structural_audit(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("kind drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids drifted")
    if len({row.get("id") for row in rows}) != len(rows):
        issues.append("duplicate row ids")
    counts = payload.get("counts", {})
    for key in (
        "stationary_phase_remainder_bounds",
        "signed_flow_bounds",
        "phi_b_bounds",
    ):
        if counts.get(key) != 0:
            issues.append(f"{key} was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "no Poisson summation remainder bound",
        "signed flow estimate",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")


def note_audit(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Six-Moment Flow Matrix and Reciprocal Phase Reduction",
        "G_j'(xi)=-i*log(a)G_j(xi)+iG_(j+1)(xi)",
        "rank(M_xi)<=8",
        "inertia (4,4,4)",
        "stationary phase a^2log(r/k)",
        "zero-phase reciprocal diagonal",
        "0 signed flow bounds",
        "this is not a proof",
        f"python work/rh_compute/scripts/check_{STEM}.py",
    ):
        if marker not in text:
            issues.append(f"note marker missing: {marker}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        independent_source_audit(payload, issues)
        independent_shift_audit(issues)
        independent_matrix_audit(issues)
        independent_rank_witness(issues)
        independent_common_phase_audit(issues)
        independent_reciprocal_phase_audit(issues)
        independent_finite_carrier_audit(issues)
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
