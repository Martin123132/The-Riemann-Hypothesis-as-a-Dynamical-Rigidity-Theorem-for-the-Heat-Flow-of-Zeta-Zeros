#!/usr/bin/env python3
"""Build the six-moment flow-matrix and reciprocal-phase reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_six_moment_flow_matrix_"
    "phase_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "five_moment_quadratic": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_explicit_"
        "phi_quadratic_reduction.json"
    ),
    "six_moment_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_saddle_flow_"
        "reduction.json"
    ),
    "reciprocal_saddle": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "reciprocal_saddle_self_duality_gate.json"
    ),
}


@dataclass(frozen=True)
class FlowMatrixRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def require_zero(expression: sp.Expr | sp.Matrix, label: str) -> None:
    if isinstance(expression, sp.MatrixBase):
        if any(sp.simplify(sp.expand(entry)) != 0 for entry in expression):
            raise RuntimeError(label)
        return
    if sp.simplify(sp.expand(expression)) != 0:
        raise RuntimeError(label)


def _real(expression: sp.Expr) -> sp.Expr:
    return sp.expand(expression).as_real_imag()[0]


def _imag(expression: sp.Expr) -> sp.Expr:
    return sp.expand(expression).as_real_imag()[1]


def _real_row(vector: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        [_real(entry) for entry in vector]
        + [-_imag(entry) for entry in vector]
    )


def _complex_structure(size: int) -> sp.Matrix:
    identity = sp.eye(size)
    zero = sp.zeros(size)
    return sp.Matrix.vstack(
        sp.Matrix.hstack(zero, -identity),
        sp.Matrix.hstack(identity, zero),
    )


def _embedding_and_shift() -> tuple[sp.Matrix, sp.Matrix]:
    embedding = sp.zeros(10, 12)
    shift = sp.zeros(10, 12)
    for index in range(5):
        embedding[index, index] = 1
        embedding[5 + index, 6 + index] = 1
        shift[index, index + 1] = 1
        shift[5 + index, 6 + index + 1] = 1
    return embedding, shift


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    quadratic = payloads["five_moment_quadratic"]
    symbolic = quadratic.get("symbolic_certificate", {})
    if "rank(M)<=4" not in symbolic.get("rank_factorization", ""):
        raise RuntimeError("five-moment rank source drifted")
    if quadratic.get("counts", {}).get("real_bulk_observations") != 4:
        raise RuntimeError("five-moment observation count drifted")

    flow = payloads["six_moment_flow"]
    if flow.get("counts", {}).get("correction_free_moments") != 6:
        raise RuntimeError("six-moment source drifted")
    if "Psi'(xi)" not in flow.get("full_current_flow_certificate", {}).get(
        "derivative", ""
    ):
        raise RuntimeError("full-flow derivative source drifted")

    reciprocal = payloads["reciprocal_saddle"].get("exact", {})
    if "a_omega^2=omega/(2*pi)" not in reciprocal.get(
        "pi_provenance", ""
    ):
        raise RuntimeError("reciprocal pi provenance drifted")
    if "u_nu=a_omega^2/nu" not in reciprocal.get("model_phase", ""):
        raise RuntimeError("reciprocal saddle source drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_signed_flow_bounds": 0,
        "imported_stationary_phase_bounds": 0,
    }


def moment_shift_certificate() -> dict:
    log_a = sp.symbols("L_a", real=True)
    coefficient_r = sp.symbols("r0_R:5", real=True)
    coefficient_i = sp.symbols("r0_I:5", real=True)
    moment_r = sp.symbols("G0_R:6", real=True)
    moment_i = sp.symbols("G0_I:6", real=True)
    coefficient = sp.Matrix(
        [coefficient_r[index] + sp.I * coefficient_i[index] for index in range(5)]
    )
    moments = sp.Matrix(
        [moment_r[index] + sp.I * moment_i[index] for index in range(6)]
    )

    padded = sp.Matrix((*coefficient, 0))
    shifted = sp.Matrix((0, *coefficient))
    lifted = sp.I * (shifted - log_a * padded)
    direct = _real(
        sum(
            coefficient[index]
            * sp.I
            * (moments[index + 1] - log_a * moments[index])
            for index in range(5)
        )
    )
    lifted_value = _real(sum(lifted[index] * moments[index] for index in range(6)))
    require_zero(direct - lifted_value, "complex derivative lift failed")

    embedding, shift = _embedding_and_shift()
    complex_structure = _complex_structure(5)
    derivative = complex_structure * (shift - log_a * embedding)
    real_moments = sp.Matrix((*moment_r, *moment_i))
    expected_derivative = sp.Matrix(
        [
            log_a * moment_i[index] - moment_i[index + 1]
            for index in range(5)
        ]
        + [
            moment_r[index + 1] - log_a * moment_r[index]
            for index in range(5)
        ]
    )
    require_zero(
        derivative * real_moments - expected_derivative,
        "real moment-shift matrix failed",
    )
    require_zero(
        derivative.T * _real_row(coefficient) - _real_row(lifted),
        "real observation derivative lift failed",
    )

    return {
        "moment_flow": (
            "For 0<=j<=4, G_j'(xi)=-i*log(a)G_j(xi)+iG_(j+1)(xi). "
            "Thus the derivative of a degree-four observation uses exactly "
            "G_0,...,G_5 and no G_6."
        ),
        "coefficient_lift": (
            "For r=(r_0,...,r_4), pad(r)=(r_0,...,r_4,0), "
            "S(r)=(0,r_0,...,r_4), and "
            "delta_a(r)=i[S(r)-log(a)pad(r)]. Then "
            "d_xi Re(r dot G_[0:4])=Re(delta_a(r) dot G_[0:5])."
        ),
        "real_shift": (
            "For y=(Re G_0,...,Re G_5,Im G_0,...,Im G_5)^T, let E "
            "select G_0,...,G_4, P select G_1,...,G_5, and "
            "J_5=[[0,-I_5],[I_5,0]]. Then x=Ey and "
            "x'=D_ay with D_a=J_5(P-log(a)E)."
        ),
        "symbolic_differences": {
            "complex_lift": str(sp.simplify(direct - lifted_value)),
            "real_shift": "zero_10_vector",
            "real_row_lift": "zero_12_vector",
        },
        "input_complex_moments": 6,
        "input_real_coordinates": 12,
    }


def observation_flow_certificate() -> dict:
    v, n, a_obs, q = sp.symbols("V N A Q", real=True)
    dv, dn, da, dq = sp.symbols("V_xi N_xi A_xi Q_xi", real=True)
    x_t, a_tx, a_t, x_tx = sp.symbols(
        "X_T A_T_x A_T X_T_x", real=True
    )
    direct = sp.expand(
        dv * n
        + v * dn
        - da * q
        - a_obs * dq
        + x_t * dn
        + a_tx * dv
        - a_t * dq
        - x_tx * da
    )
    grouped = sp.expand(
        dv * (n + a_tx)
        + (x_t + v) * dn
        - da * (q + x_tx)
        - (a_t + a_obs) * dq
    )
    require_zero(direct - grouped, "four-observation flow failed")

    dimension = 5
    u0_symbols = sp.symbols(f"u0:{dimension * 4}", real=True)
    u1_symbols = sp.symbols(f"d0:{dimension * 4}", real=True)
    state = sp.Matrix(sp.symbols(f"y0:{dimension}", real=True))
    u0 = sp.Matrix(dimension, 4, u0_symbols)
    u1 = sp.Matrix(dimension, 4, u1_symbols)
    core = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    tail = sp.Matrix([a_tx, x_t, -x_tx, -a_t])
    observation = u0.T * state
    derivative_observation = u1.T * state
    flow_matrix = u0 * core * u1.T + u1 * core * u0.T
    flow_linear = u1 * tail
    matrix_value = sp.expand(
        (state.T * flow_matrix * state)[0]
        + (flow_linear.T * state)[0]
    )
    observation_value = sp.expand(
        2 * (observation.T * core * derivative_observation)[0]
        + (tail.T * derivative_observation)[0]
    )
    require_zero(
        matrix_value - observation_value,
        "flow-matrix factorization failed",
    )
    if flow_matrix != flow_matrix.T:
        raise RuntimeError("flow matrix is not symmetric")

    contact = sp.expand(grouped.subs({v: -x_t, a_obs: -a_t}))
    expected_contact = sp.expand(dv * (n + a_tx) - da * (q + x_tx))
    require_zero(contact - expected_contact, "double-contact flow fibre failed")

    return {
        "observation_rows": (
            "Retain the Section 11.159 complex coefficient rows "
            "v_0,n_0,a_0,q_0 and put dot(r)=delta_a(r). For y in R^12, "
            "V=tilde(v)^Ty, N=tilde(n)^Ty, A=tilde(a)^Ty, Q=tilde(q)^Ty "
            "and V_xi=dot(v)^Ty, N_xi=dot(n)^Ty, "
            "A_xi=dot(a)^Ty, Q_xi=dot(q)^Ty."
        ),
        "explicit_flow": (
            "Psi'=V_xi*N+V*N_xi-A_xi*Q-A*Q_xi+X_T*N_xi+"
            "A_(T,x)*V_xi-A_T*Q_xi-X_(T,x)*A_xi. Equivalently, "
            "Psi'=V_xi[N+A_(T,x)]+[X_T+V]N_xi-"
            "A_xi[Q+X_(T,x)]-[A_T+A]Q_xi."
        ),
        "real_quadratic": (
            "Let U_0=[tilde(v),tilde(n),tilde(a),tilde(q)], "
            "U_1=[dot(v),dot(n),dot(a),dot(q)], "
            "J=(1/2)diag_offdiag(+1,-1), and "
            "t_T=(A_(T,x),X_T,-X_(T,x),-A_T)^T. Then "
            "Psi'=y^T M_xi y+ell_xi^T y, where "
            "M_xi=U_0JU_1^T+U_1JU_0^T and ell_xi=U_1t_T."
        ),
        "eight_observation_factor": (
            "With W=[U_0,U_1] and K=[[0,J],[J,0]], "
            "M_xi=WKW^T. Hence rank(M_xi)<=8, ell_xi lies in range(U_1), "
            "and Psi' depends on at most eight real observations."
        ),
        "kernel": (
            "Every y satisfying U_0^Ty=U_1^Ty=0 has Psi'(xi)=0. "
            "The common flow-invisible kernel has dimension at least four."
        ),
        "double_contact_fibre": (
            "On X_T+V=A_T+A=0, Psi'=V_xi[N+A_(T,x)]-"
            "A_xi[Q+X_(T,x)]. This is division-free but has no algebraic sign."
        ),
        "symbolic_difference": str(sp.simplify(matrix_value - observation_value)),
    }


def _witness_vectors() -> tuple[sp.Matrix, sp.Expr]:
    log_n = sp.Integer(2)
    rho_1 = 1 + sp.I
    rho_2 = 2 - sp.I
    rho_1_x = 1 - 2 * sp.I
    rho_2_x = -1 + sp.I
    s_prime = 1 - sp.I / 2
    s_second = 2 + sp.I
    kappa = 1 + sp.I
    b = -sp.Rational(1, 2)
    b_x = sp.Integer(1)
    u_n = sp.Integer(1)
    u_n_x = sp.Integer(1)
    log_a = sp.Integer(3)

    v_0 = sp.Matrix([1, rho_1, rho_2, 0, 0])
    v_1 = sp.Matrix(
        [
            log_n,
            log_n * rho_1 - 1,
            log_n * rho_2 - rho_1,
            -rho_2,
            0,
        ]
    )
    v_2 = sp.Matrix(
        [
            log_n**2,
            log_n**2 * rho_1 - 2 * log_n,
            log_n**2 * rho_2 - 2 * log_n * rho_1 + 1,
            rho_1 - 2 * log_n * rho_2,
            rho_2,
        ]
    )
    e_0 = sp.Matrix([0, rho_1_x, rho_2_x, 0, 0])
    e_1 = sp.Matrix(
        [
            0,
            log_n * rho_1_x,
            log_n * rho_2_x - rho_1_x,
            -rho_2_x,
            0,
        ]
    )
    q_0 = kappa * v_0 + s_prime * v_1 + e_0
    q_1 = kappa * v_1 + s_prime * v_2 + e_1
    a_0 = s_prime * v_1 + sp.I * b * u_n * v_0
    n_0 = (
        s_second * v_1
        + s_prime * q_1
        + sp.I * (b_x * u_n + b * u_n_x) * v_0
        + sp.I * b * u_n * q_0
    )
    base = sp.Matrix.hstack(
        _real_row(v_0), _real_row(n_0), _real_row(a_0), _real_row(q_0)
    )
    return base, log_a


def rank_certificate() -> dict:
    base, log_a = _witness_vectors()
    embedding, shift = _embedding_and_shift()
    derivative = _complex_structure(5) * (shift - log_a * embedding)
    u0 = embedding.T * base
    u1 = derivative.T * base
    factor = sp.Matrix.hstack(u0, u1)
    minor = sp.factor(factor[list(range(8)), :].det())
    if base.rank() != 4:
        raise RuntimeError("rank witness base observations failed")
    if factor.rank() != 8 or minor != sp.Rational(51095, 16):
        raise RuntimeError("rank-eight flow witness failed")

    core = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    zero = sp.zeros(4)
    doubled_core = sp.Matrix.vstack(
        sp.Matrix.hstack(zero, core),
        sp.Matrix.hstack(core, zero),
    )
    eigenvalues = doubled_core.eigenvals()
    positive = sum(
        multiplicity
        for value, multiplicity in eigenvalues.items()
        if value.is_positive
    )
    negative = sum(
        multiplicity
        for value, multiplicity in eigenvalues.items()
        if value.is_negative
    )
    if (positive, negative) != (4, 4):
        raise RuntimeError("doubled core inertia failed")
    flow_matrix = factor * doubled_core * factor.T
    if flow_matrix.rank() != 8:
        raise RuntimeError("rank-eight flow matrix failed")

    return {
        "generic_rank_guard": (
            "An exact algebraic specialization log(N)=2, log(a)=3, "
            "rho_1=1+i, rho_2=2-i, rho_(1,x)=1-2i, "
            "rho_(2,x)=-1+i, s_*'=1-i/2, s_*''=2+i, "
            "chi_N=1+i, b=-1/2, b_x=1, and u_N=u_(N,x)=1 has "
            "rank(W)=8 and leading 8-by-8 minor 51095/16."
        ),
        "generic_inertia_guard": (
            "Since K=[[0,J],[J,0]] has inertia (4,4), the witness flow "
            "matrix has inertia (4,4,4). Therefore the rank-eight bound is "
            "generically sharp and the flow is not coefficient-blind "
            "semidefinite. The witness is algebraic, not a physical Xi state."
        ),
        "witness_minor": str(minor),
        "factor_rank": factor.rank(),
        "flow_rank": flow_matrix.rank(),
        "flow_inertia": [4, 4, 4],
    }


def common_phase_certificate() -> dict:
    size = 2
    log_a = sp.symbols("L_a", real=True)
    entries = sp.symbols("m00 m01 m02 m11 m12 m22 m03 m13 m23 m33", real=True)
    m00, m01, m02, m11, m12, m22, m03, m13, m23, m33 = entries
    matrix = sp.Matrix(
        [
            [m00, m01, m02, m03],
            [m01, m11, m12, m13],
            [m02, m12, m22, m23],
            [m03, m13, m23, m33],
        ]
    )
    complex_structure = _complex_structure(size)
    hermitian_part = sp.expand(
        (matrix - complex_structure * matrix * complex_structure) / 2
    )
    transpose_part = sp.expand(
        (matrix + complex_structure * matrix * complex_structure) / 2
    )
    require_zero(
        hermitian_part * complex_structure
        - complex_structure * hermitian_part,
        "Hermitian commuting projection failed",
    )
    require_zero(
        transpose_part * complex_structure
        + complex_structure * transpose_part,
        "transpose anticommuting projection failed",
    )
    common_flow = sp.expand(
        log_a * (complex_structure * matrix - matrix * complex_structure)
    )
    require_zero(
        common_flow + 2 * log_a * transpose_part * complex_structure,
        "common-phase flow projection failed",
    )

    xi, omega, ell_a, ell = sp.symbols(
        "xi Omega ell_a ell", real=True
    )
    c_xi = sp.exp(sp.I * (omega - xi) * ell_a)
    q_xi = sp.exp(sp.I * xi * ell)
    carrier = c_xi * q_xi
    require_zero(
        sp.diff(carrier, xi) + sp.I * (ell_a - ell) * carrier,
        "external common-phase derivative failed",
    )

    return {
        "normalized_moments": (
            "Write G_j(xi)=c_xi*Ghat_j(xi), where "
            "c_xi=exp[i(Omega-xi)log(a)] and Ghat_j contains the full "
            "Mangoldt or balanced divisor sum but no c_xi. Then "
            "Ghat_j'(xi)=iGhat_(j+1)(xi), while c_xi'=-i log(a)c_xi."
        ),
        "family_separation": (
            "The endpoint family carries c_xi, the Hermitian family carries "
            "c_xi*conj(c_xi)=1, and the transpose family carries c_xi^2. "
            "Hence log(a) cancels exactly from the Hermitian flow and remains "
            "only through c_xi' and (c_xi^2)' in the endpoint and transpose "
            "families. No common phase enters a divisor sum."
        ),
        "matrix_projection": (
            "For J_5=[[0,-I],[I,0]], decompose a real symmetric current "
            "matrix M as M_H=(M-J_5MJ_5)/2 and M_T=(M+J_5MJ_5)/2. "
            "Then M_HJ_5=J_5M_H, M_TJ_5=-J_5M_T, and the common-rotation "
            "quadratic flow is log(a)(J_5M-MJ_5)=-2log(a)M_TJ_5. "
            "The Hermitian part contributes zero."
        ),
        "flow_multipliers": (
            "In logarithmic coordinates ell_n=log(n), the three exact "
            "multipliers are i(ell_n-log(a)), i(ell_n-ell_m), and "
            "i(ell_n+ell_m-2log(a)) for endpoint, Hermitian, and transpose "
            "respectively."
        ),
        "symbolic_differences": {
            "hermitian_commutator": "zero_4_matrix",
            "transpose_anticommutator": "zero_4_matrix",
            "common_flow": "zero_4_matrix",
            "external_phase": "0",
        },
    }


def reciprocal_phase_certificate() -> dict:
    a, n, m, k, r = sp.symbols("a n m k r", positive=True)
    endpoint = a**2 * sp.log(n) - k * n
    hermitian = a**2 * sp.log(n) - a**2 * sp.log(m) - k * n + r * m
    transpose = a**2 * sp.log(n) + a**2 * sp.log(m) - k * n - r * m
    n_star = a**2 / k
    m_star = a**2 / r

    require_zero(
        sp.diff(endpoint, n).subs(n, n_star),
        "endpoint reciprocal saddle failed",
    )
    require_zero(
        sp.diff(hermitian, n).subs(n, n_star),
        "Hermitian n saddle failed",
    )
    require_zero(
        sp.diff(hermitian, m).subs(m, m_star),
        "Hermitian m saddle failed",
    )
    require_zero(
        sp.diff(transpose, n).subs(n, n_star),
        "transpose n saddle failed",
    )
    require_zero(
        sp.diff(transpose, m).subs(m, m_star),
        "transpose m saddle failed",
    )

    endpoint_curvature = sp.simplify(sp.diff(endpoint, n, 2).subs(n, n_star))
    hermitian_hessian = sp.diag(
        sp.simplify(sp.diff(hermitian, n, 2).subs(n, n_star)),
        sp.simplify(sp.diff(hermitian, m, 2).subs(m, m_star)),
    )
    transpose_hessian = sp.diag(
        sp.simplify(sp.diff(transpose, n, 2).subs(n, n_star)),
        sp.simplify(sp.diff(transpose, m, 2).subs(m, m_star)),
    )
    require_zero(
        endpoint_curvature + k**2 / a**2,
        "endpoint reciprocal curvature failed",
    )
    require_zero(
        hermitian_hessian - sp.diag(-k**2 / a**2, r**2 / a**2),
        "Hermitian reciprocal Hessian failed",
    )
    require_zero(
        transpose_hessian - sp.diag(-k**2 / a**2, -r**2 / a**2),
        "transpose reciprocal Hessian failed",
    )

    hermitian_at_saddle = sp.expand_log(
        hermitian.subs({n: n_star, m: m_star}), force=True
    )
    expected_hermitian = sp.expand_log(a**2 * sp.log(r / k), force=True)
    require_zero(
        hermitian_at_saddle - expected_hermitian,
        "Hermitian stationary phase failed",
    )
    hermitian_multiplier = sp.expand_log(
        sp.log(n_star / m_star), force=True
    )
    expected_multiplier = sp.expand_log(sp.log(r / k), force=True)
    require_zero(
        hermitian_multiplier - expected_multiplier,
        "Hermitian stationary multiplier failed",
    )
    require_zero(
        hermitian_at_saddle / a**2 - hermitian_multiplier,
        "Hermitian phase-multiplier lock failed",
    )
    require_zero(
        hermitian_at_saddle.subs(r, k),
        "Hermitian reciprocal diagonal phase failed",
    )
    require_zero(
        hermitian_multiplier.subs(r, k),
        "Hermitian reciprocal diagonal multiplier failed",
    )

    return {
        "pi_provenance": (
            "Use e(t)=exp(2*pi*i*t). The T_0=2*pi*a^2 factor is inherited "
            "from the completed-zeta/Riemann-Siegel saddle. Dividing its "
            "radian phase by the standard Poisson 2*pi gives a^2. These are "
            "the two appearances of the same normalization; no circle or "
            "fitted pi is introduced."
        ),
        "poisson_convention": (
            "At xi=T_0 and with Poisson character e(-kn), discard only "
            "n-independent unit phases when locating saddles. The endpoint "
            "cycle phase is a^2log(n)-kn. For a two-variable mode, use "
            "a^2log(n)-a^2log(m)-kn+rm in the Hermitian family and "
            "a^2log(n)+a^2log(m)-kn-rm in the transpose family, with k,r>0."
        ),
        "endpoint_saddle": (
            "The endpoint mode has n_*=a^2/k, phase "
            "a^2[log(a^2/k)-1], curvature -k^2/a^2, and flow multiplier "
            "log(a/k). Its retained-support window 1<=n_*<=B is exactly "
            "a^2/B<=k<=a^2."
        ),
        "hermitian_saddle": (
            "The Hermitian mode has (n_*,m_*)=(a^2/k,a^2/r), stationary "
            "phase a^2log(r/k), Hessian diag(-k^2/a^2,+r^2/a^2), and flow "
            "multiplier log(r/k). Thus the multiplier is exactly the "
            "stationary phase divided by a^2."
        ),
        "hermitian_diagonal_null": (
            "On the reciprocal diagonal k=r, both the Hermitian stationary "
            "phase and its differentiated-flow multiplier vanish exactly. "
            "The zero-phase reciprocal diagonal channel is therefore absent from "
            "Psi'_H, although near-diagonal modes still require an estimate."
        ),
        "transpose_saddle": (
            "The transpose mode has (n_*,m_*)=(a^2/k,a^2/r), phase "
            "a^2[log(a^4/(kr))-2], Hessian "
            "diag(-k^2/a^2,-r^2/a^2), and flow multiplier "
            "log(a^2/(kr))."
        ),
        "support_signs": (
            "Because B<a, every interior retained saddle has k,r>a. Hence "
            "the endpoint multiplier log(a/k) and transpose multiplier "
            "log(a^2/(kr)) are negative there. They multiply complex "
            "coefficients and unit phases, so this is not a current-sign "
            "theorem. The Hermitian multiplier is antisymmetric under k<->r."
        ),
        "external_phases": (
            "Since Omega=T_0-epsilon, c_(T_0)=exp[-i epsilon log(a)]. "
            "The endpoint Poisson family retains this unit phase, the "
            "Hermitian family has none, and the transpose family retains "
            "exp[-2i epsilon log(a)]."
        ),
        "active_mode_window": (
            "For each primal variable in [1,B], the only continuous "
            "stationary modes satisfy a^2/B<=mode<=a^2. Modes outside this "
            "window are nonstationary, but a discrete Poisson formula, hard-"
            "cutoff boundary terms, and uniform integration-by-parts constants "
            "remain to be proved."
        ),
        "symbolic": {
            "endpoint_curvature": str(endpoint_curvature),
            "hermitian_hessian": str(hermitian_hessian),
            "transpose_hessian": str(transpose_hessian),
            "hermitian_stationary_phase": str(hermitian_at_saddle),
            "hermitian_stationary_multiplier": str(hermitian_multiplier),
        },
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use the q=1, L>=50, fixed-N chart, B=N-M-1, and the exact "
            "auxiliary frequency segment T_0-epsilon<=xi<=T_0 from Section "
            "11.164. All physical tail jets and coefficient rows are held "
            "fixed along xi."
        ),
        "arithmetic_handoff": (
            "Insert the balanced Mangoldt square-and-wing formula for each "
            "Ghat_0,...,Ghat_5 into the eight observations "
            "(V,N,A,Q,V_xi,N_xi,A_xi,Q_xi). Keep c_xi outside the sums. "
            "Apply one- and two-variable Poisson/B-process transforms before "
            "taking absolute values, preserving endpoint, Hermitian, and "
            "transpose terms in one signed expression."
        ),
        "new_cancellation_clue": (
            "The canonical zero phase-difference channel is the Hermitian "
            "reciprocal diagonal, and the differentiated kernel kills it "
            "exactly. The next estimate should exploit the resulting "
            "log(r/k) factor by dual antisymmetry or summation by parts. The "
            "endpoint and transpose saddles remain and must be composed with "
            "the certified terminal recurrence rather than bounded separately."
        ),
        "exact_target": (
            "Prove Psi(T_0)-epsilon*integral_0^1 "
            "Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800, "
            "by an endpoint-complete signed reciprocal-Poisson estimate for "
            "the eight-observation normal form."
        ),
        "falsification_program": (
            "Any candidate estimate must retain the hard-cutoff boundary, "
            "adjacent saddle recurrence, c_xi, correction derivative D, "
            "M_0=0, X_T+V=0, double contact, W_0=0, p=+-1, removable "
            "p=+-1/2, prime edges, square/cube transports, n=64, and q=1."
        ),
        "proof_boundary": (
            "This proves an exact six-moment/12-real flow matrix, rank at "
            "most eight, generic inertia (4,4,4), common-phase separation, "
            "three reciprocal stationary-phase skeletons, and the Hermitian "
            "reciprocal-diagonal null. It proves no Poisson summation remainder "
            "bound, hard-cutoff boundary estimate, near-diagonal bilinear "
            "gain, signed flow estimate, saddle-proxy or Phi_B upper bound, "
            "contact exclusion, retained aggregate or Xi theorem, Q209, "
            "cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or "
            "prize-level conclusion."
        ),
    }


def build_rows(
    shift: dict,
    observations: dict,
    rank: dict,
    phase: dict,
    reciprocal: dict,
    exact: dict,
) -> list[FlowMatrixRow]:
    return [
        FlowMatrixRow("fmpr_01_domain", "fixed-chart domain", "proved", "The flow normal form uses the certified auxiliary segment.", exact["domain"], "Intermediate xi is auxiliary, not a new physical source state."),
        FlowMatrixRow("fmpr_02_shift", "six-moment shift", "proved", "Frequency differentiation is one exact shifted moment operator.", shift["moment_flow"], "Only derivatives of G_0 through G_4 are used."),
        FlowMatrixRow("fmpr_03_lift", "coefficient lift", "proved", "Every degree-four complex observation has one explicit degree-five derivative row.", shift["coefficient_lift"], shift["real_shift"], shift["symbolic_differences"]),
        FlowMatrixRow("fmpr_04_observations", "eight observations", "proved", "The full flow uses four values and their four frequency derivatives.", observations["observation_rows"], "No denominator or relative moment is introduced."),
        FlowMatrixRow("fmpr_05_flow", "explicit flow polynomial", "proved", "Psi' is one displayed division-free polynomial in eight observations.", observations["explicit_flow"], "The endpoint, Hermitian, and transpose terms have not been separately bounded."),
        FlowMatrixRow("fmpr_06_matrix", "12-real flow matrix", "proved", "The six complex moments give one symmetric real quadratic plus a linear term.", observations["real_quadratic"], "This is an identity, not a sign estimate."),
        FlowMatrixRow("fmpr_07_factor", "rank-eight factorization", "proved", "The flow matrix factors through eight real observations.", observations["eight_observation_factor"], observations["kernel"]),
        FlowMatrixRow("fmpr_08_rank", "generic rank witness", "proved", "The rank-eight upper bound is algebraically sharp.", rank["generic_rank_guard"], "The witness is not claimed to occur on the Xi source.", {"minor": rank["witness_minor"], "rank": rank["flow_rank"]}),
        FlowMatrixRow("fmpr_09_inertia", "generic inertia guard", "proved", "Coefficient-blind semidefiniteness is impossible.", rank["generic_inertia_guard"], "Physical phase cancellation remains possible.", {"inertia": rank["flow_inertia"]}),
        FlowMatrixRow("fmpr_10_fibre", "exceptional fibres", "proved", "The matrix form remains valid on the contact and zero-moment fibres.", observations["double_contact_fibre"], "No contact exclusion or flow sign is inferred."),
        FlowMatrixRow("fmpr_11_normalize", "external common phase", "proved", "The common phase can be removed from every arithmetic divisor sum.", phase["normalized_moments"], "The phase remains in the endpoint and transpose prefactors."),
        FlowMatrixRow("fmpr_12_rotation", "common-rotation cancellation", "proved", "The log(a) rotation cancels exactly from the Hermitian flow.", phase["family_separation"] + " " + phase["matrix_projection"], phase["flow_multipliers"], phase["symbolic_differences"]),
        FlowMatrixRow("fmpr_13_pi", "pi provenance", "proved", "The reciprocal scale introduces no unexplained pi.", reciprocal["pi_provenance"], "The standard Poisson character fixes the normalization."),
        FlowMatrixRow("fmpr_14_poisson", "Poisson phase convention", "proved", "All three flow families have explicit cycle phases at T_0.", reciprocal["poisson_convention"], "This locates continuous saddles only."),
        FlowMatrixRow("fmpr_15_endpoint", "endpoint saddle", "proved", "Every active endpoint mode has one reciprocal saddle.", reciprocal["endpoint_saddle"], reciprocal["external_phases"]),
        FlowMatrixRow("fmpr_16_hermitian", "Hermitian saddle", "proved", "The phase-difference family has an indefinite reciprocal saddle.", reciprocal["hermitian_saddle"], "No two-variable stationary-phase remainder is supplied."),
        FlowMatrixRow("fmpr_17_null", "Hermitian diagonal null", "proved", "The zero-phase reciprocal diagonal is killed by the flow multiplier.", reciprocal["hermitian_diagonal_null"], "Near-diagonal modes remain an open bilinear problem."),
        FlowMatrixRow("fmpr_18_transpose", "transpose saddle", "proved", "The phase-sum family has a negative-definite reciprocal saddle.", reciprocal["transpose_saddle"], reciprocal["support_signs"]),
        FlowMatrixRow("fmpr_19_window", "active dual window", "proved", "The exact retained saddle window is a^2/B<=mode<=a^2.", reciprocal["active_mode_window"], "Boundary and nonstationary constants remain open."),
        FlowMatrixRow("fmpr_20_handoff", "signed arithmetic handoff", "open", "The eight-observation form must be estimated after reciprocal composition.", exact["arithmetic_handoff"] + " " + exact["new_cancellation_clue"], exact["exact_target"]),
        FlowMatrixRow("fmpr_21_guard", "falsification guard", "open", "Every proposed cancellation estimate must survive the existing exceptional cases.", exact["falsification_program"], "No guard may be removed by genericity."),
        FlowMatrixRow("fmpr_22_boundary", "proof boundary", "proved", "The exact algebra is separated from the missing analytic estimate.", exact["proof_boundary"], "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    shift = payload["moment_shift_certificate"]
    observations = payload["observation_flow_certificate"]
    rank = payload["rank_certificate"]
    phase = payload["common_phase_certificate"]
    reciprocal = payload["reciprocal_phase_certificate"]
    exact = payload["exact"]
    return f"""# Six-Moment Flow Matrix and Reciprocal Phase Reduction

Date: 2026-08-01

Status: exact matrix/phase reduction with `0 signed flow bounds`, `0
stationary-phase remainder bounds`, and `0 Phi_B bounds`; this is not a proof
of RH.

## Domain

{exact['domain']}

## Six-Moment Shift

{shift['moment_flow']}

{shift['coefficient_lift']}

{shift['real_shift']}

## Eight-Observation Flow

{observations['observation_rows']}

{observations['explicit_flow']}

{observations['real_quadratic']}

{observations['eight_observation_factor']}

{observations['kernel']}

{observations['double_contact_fibre']}

## Rank and Inertia

{rank['generic_rank_guard']}

{rank['generic_inertia_guard']}

## Common Phase

{phase['normalized_moments']}

{phase['family_separation']}

{phase['matrix_projection']}

{phase['flow_multipliers']}

## Reciprocal Phases

{reciprocal['pi_provenance']}

{reciprocal['poisson_convention']}

{reciprocal['endpoint_saddle']}

{reciprocal['hermitian_saddle']}

{reciprocal['hermitian_diagonal_null']}

{reciprocal['transpose_saddle']}

{reciprocal['support_signs']}

{reciprocal['external_phases']}

{reciprocal['active_mode_window']}

## Handoff

{exact['arithmetic_handoff']}

{exact['new_cancellation_clue']}

{exact['exact_target']}

{exact['falsification_program']}

## Proof Boundary

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    sources = load_sources()
    shift = moment_shift_certificate()
    observations = observation_flow_certificate()
    rank = rank_certificate()
    phase = common_phase_certificate()
    reciprocal = reciprocal_phase_certificate()
    exact = exact_payload()
    rows = build_rows(shift, observations, rank, phase, reciprocal, exact)
    return {
        "kind": STEM,
        "date": "2026-08-01",
        "status": (
            "exact six-moment flow matrix, common-phase separation, and "
            "reciprocal-phase skeleton complete; signed analytic estimate open"
        ),
        "source_audit": source_audit(sources),
        "moment_shift_certificate": shift,
        "observation_flow_certificate": observations,
        "rank_certificate": rank,
        "common_phase_certificate": phase,
        "reciprocal_phase_certificate": reciprocal,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
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
        },
        "proof_boundary": exact["proof_boundary"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "built six-moment flow matrix/phase reduction: 22 rows, "
        "rank <=8, generic inertia (4,4,4), 3 reciprocal phase families, "
        "1 Hermitian diagonal null, 0 signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
