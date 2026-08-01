#!/usr/bin/env python3
"""Build the endpoint, Hermitian, and transpose two-carrier Phi_B kernel."""

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
    "contiguous_terminal_tail_anchored_five_moment_two_carrier_"
    "kernel_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "explicit_phi": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_explicit_phi_"
        "quadratic_reduction.json"
    ),
    "mangoldt": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_current_mangoldt_"
        "normal_form_gate.json"
    ),
    "growing_tail": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_growing_prefix_finite_height_gate.json"
    ),
    "q1_finite_edge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "q1_finite_height_real_edge_remainder_gate.json"
    ),
    "carrier_rotation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
}


@dataclass(frozen=True)
class KernelRow:
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


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    explicit = payloads["explicit_phi"]
    counts = explicit.get("counts", {})
    if counts.get("real_bulk_observations") != 4:
        raise RuntimeError("explicit-Phi observation source drifted")
    if counts.get("signed_phi_bounds") != 0:
        raise RuntimeError("explicit-Phi proof boundary drifted")

    mangoldt = payloads["mangoldt"].get("counts", {})
    if mangoldt.get("mangoldt_moment_orders") != 4:
        raise RuntimeError("Mangoldt order source drifted")
    if mangoldt.get("signed_type_ii_bounds") != 0:
        raise RuntimeError("Mangoldt proof boundary drifted")

    tail = payloads["growing_tail"].get("exact", {})
    for marker in (
        "|s_*'+i/2|<h^2/6000",
        "|s_*''|<h^4/40000",
        "u_x=h^2/(8*pi)",
    ):
        if marker not in tail.get("source_bounds", ""):
            raise RuntimeError(f"growing-tail source marker missing: {marker}")
    finite_edge = payloads["q1_finite_edge"]
    finite_edge_text = json.dumps(finite_edge, sort_keys=True)
    if "|W_x|<2h^2" not in finite_edge_text:
        raise RuntimeError("terminal logarithmic derivative source drifted")
    if "Q_x/Q=-i*h*theta/2" not in finite_edge_text:
        raise RuntimeError("terminal phase derivative source drifted")

    rotation = payloads["carrier_rotation"].get("exact", {})
    relative = rotation.get("relative_carrier_rotation", "")
    for marker in (
        "|d_n|<2189/x<1/2",
        "|d_(n,x)|<4223/x^2",
        "d_(n,x)/(1+d_n)-d_(1,x)/(1+d_1)",
    ):
        if marker not in relative:
            raise RuntimeError(f"carrier correction marker missing: {marker}")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_phi_bounds": 0,
        "imported_type_ii_bounds": 0,
    }


def _complex_symbols(prefix: str) -> tuple[sp.Symbol, sp.Symbol, sp.Expr]:
    real, imag = sp.symbols(f"{prefix}_R {prefix}_I", real=True)
    return real, imag, real + sp.I * imag


def symbolic_certificate() -> dict:
    log_n, ell_n, ell_m = sp.symbols("L ell_n ell_m", real=True)
    u_terminal, b, b_x, u_x = sp.symbols("u_N b b_x u_N_x", real=True)
    _, _, rho_1 = _complex_symbols("rho_1")
    _, _, rho_2 = _complex_symbols("rho_2")
    _, _, rho_1_x = _complex_symbols("rho_1_x")
    _, _, rho_2_x = _complex_symbols("rho_2_x")
    _, _, s_prime = _complex_symbols("s_prime")
    _, _, s_second = _complex_symbols("s_second")
    _, _, chi_n = _complex_symbols("chi_N")

    def distance(ell: sp.Expr) -> sp.Expr:
        return log_n - ell

    def correction(ell: sp.Expr) -> sp.Expr:
        return 1 + rho_1 * ell + rho_2 * ell**2

    def correction_x(ell: sp.Expr) -> sp.Expr:
        return rho_1_x * ell + rho_2_x * ell**2

    def radial(ell: sp.Expr) -> sp.Expr:
        return s_prime * distance(ell) + sp.I * b * u_terminal

    def radial_x(ell: sp.Expr) -> sp.Expr:
        return (
            s_second * distance(ell)
            + sp.I * (b_x * u_terminal + b * u_x)
        )

    def value_rate(ell: sp.Expr) -> sp.Expr:
        return (
            (chi_n + s_prime * distance(ell)) * correction(ell)
            + correction_x(ell)
        )

    def scalar_rate(ell: sp.Expr) -> sp.Expr:
        return (
            radial(ell) * value_rate(ell)
            + radial_x(ell) * correction(ell)
        )

    c_n = correction(ell_n)
    c_m = correction(ell_m)
    d_n = correction_x(ell_n)
    d_m = correction_x(ell_m)
    r_n = radial(ell_n)
    r_m = radial(ell_m)
    r_n_x = radial_x(ell_n)
    r_m_x = radial_x(ell_m)
    q_n = value_rate(ell_n)
    q_m = value_rate(ell_m)
    n_n = scalar_rate(ell_n)
    n_m = scalar_rate(ell_m)

    expanded_n_m = (
        s_second * distance(ell_m) * c_m
        + s_prime
        * (
            chi_n * distance(ell_m) * c_m
            + s_prime * distance(ell_m) ** 2 * c_m
            + distance(ell_m) * d_m
        )
        + sp.I * (b_x * u_terminal + b * u_x) * c_m
        + sp.I * b * u_terminal * q_m
    )
    if sp.simplify(n_m - expanded_n_m) != 0:
        raise RuntimeError("scalar-rate carrier factorization failed")

    defect = sp.expand(r_n * c_n - q_n)
    expected_defect = sp.expand(
        (sp.I * b * u_terminal - chi_n) * c_n - d_n
    )
    if sp.simplify(defect - expected_defect) != 0:
        raise RuntimeError("radial-subtracted defect failed")

    x_t, a_t, x_t_x, a_t_x = sp.symbols(
        "X_T A_T X_T_x A_T_x", real=True
    )
    endpoint = (
        x_t * n_n + a_t_x * c_n - a_t * q_n - x_t_x * r_n * c_n
    )
    endpoint_factored = (
        (x_t * r_n - a_t) * q_n
        + (x_t * r_n_x + a_t_x - x_t_x * r_n) * c_n
    )
    if sp.simplify(endpoint - endpoint_factored) != 0:
        raise RuntimeError("endpoint-linear kernel factorization failed")

    hermitian = c_n * sp.conjugate(n_m) - r_n * c_n * sp.conjugate(q_m)
    transpose = c_n * n_m - r_n * c_n * q_m
    hermitian_factored = c_n * (
        (sp.conjugate(r_m) - r_n) * sp.conjugate(q_m)
        + sp.conjugate(r_m_x) * sp.conjugate(c_m)
    )
    transpose_factored = c_n * (
        (r_m - r_n) * q_m + r_m_x * c_m
    )
    if sp.simplify(hermitian - hermitian_factored) != 0:
        raise RuntimeError("Hermitian kernel factorization failed")
    if sp.simplify(transpose - transpose_factored) != 0:
        raise RuntimeError("transpose kernel factorization failed")

    hermitian_symmetric = sp.expand(
        (hermitian + sp.conjugate(
            c_m * sp.conjugate(n_n) - r_m * c_m * sp.conjugate(q_n)
        ))
        / 2
    )
    hermitian_symmetric_expected = sp.expand(
        (
            (sp.conjugate(r_m) - r_n)
            * (c_n * sp.conjugate(q_m) - sp.conjugate(c_m) * q_n)
            + c_n
            * sp.conjugate(c_m)
            * (sp.conjugate(r_m_x) + r_n_x)
        )
        / 2
    )
    if sp.simplify(
        hermitian_symmetric - hermitian_symmetric_expected
    ) != 0:
        raise RuntimeError("Hermitian symmetrization failed")

    transpose_reverse = c_m * n_n - r_m * c_m * q_n
    transpose_symmetric = sp.expand((transpose + transpose_reverse) / 2)
    delta_lambda = distance(ell_m) - distance(ell_n)
    transpose_symmetric_expected = sp.expand(
        (
            s_prime**2 * delta_lambda**2 * c_n * c_m
            + s_prime
            * delta_lambda
            * (c_n * d_m - c_m * d_n)
            + c_n * c_m * (r_n_x + r_m_x)
        )
        / 2
    )
    if sp.simplify(
        transpose_symmetric - transpose_symmetric_expected
    ) != 0:
        raise RuntimeError("transpose symmetrization failed")

    hermitian_difference = sp.expand(
        c_n * sp.conjugate(q_m) - sp.conjugate(c_m) * q_n
    )
    hermitian_difference_expected = sp.expand(
        c_n
        * sp.conjugate(c_m)
        * (
            sp.conjugate(chi_n)
            - chi_n
            + sp.conjugate(s_prime) * distance(ell_m)
            - s_prime * distance(ell_n)
        )
        + c_n * sp.conjugate(d_m)
        - sp.conjugate(c_m) * d_n
    )
    if sp.simplify(
        hermitian_difference - hermitian_difference_expected
    ) != 0:
        raise RuntimeError("Hermitian common-rate cancellation failed")

    leading = {
        symbol: 0
        for symbol in (
            *sp.re(rho_1).free_symbols,
            *sp.im(rho_1).free_symbols,
            *sp.re(rho_2).free_symbols,
            *sp.im(rho_2).free_symbols,
            *sp.re(rho_1_x).free_symbols,
            *sp.im(rho_1_x).free_symbols,
            *sp.re(rho_2_x).free_symbols,
            *sp.im(rho_2_x).free_symbols,
        )
    }
    leading.update(
        {
            sp.re(s_prime): 0,
            sp.im(s_prime): -sp.Rational(1, 2),
            sp.re(s_second): 0,
            sp.im(s_second): 0,
            sp.re(chi_n): 0,
            sp.im(chi_n): -u_terminal / 2,
            b: -sp.Rational(1, 2),
            b_x: 0,
        }
    )
    # The expressions above use component symbols rather than re/im nodes;
    # build the substitution by symbol name to keep it exact.
    all_symbols = set().union(
        hermitian_symmetric.free_symbols,
        transpose_symmetric.free_symbols,
    )
    by_name = {str(symbol): symbol for symbol in all_symbols}
    leading = {
        by_name[name]: value
        for name, value in {
            "rho_1_R": 0,
            "rho_1_I": 0,
            "rho_2_R": 0,
            "rho_2_I": 0,
            "rho_1_x_R": 0,
            "rho_1_x_I": 0,
            "rho_2_x_R": 0,
            "rho_2_x_I": 0,
            "s_prime_R": 0,
            "s_prime_I": -sp.Rational(1, 2),
            "s_second_R": 0,
            "s_second_I": 0,
            "chi_N_R": 0,
            "chi_N_I": -u_terminal / 2,
            "b": -sp.Rational(1, 2),
            "b_x": 0,
        }.items()
        if name in by_name
    }
    lambda_n = distance(ell_n)
    lambda_m = distance(ell_m)
    leading_hermitian = sp.simplify(hermitian_symmetric.subs(leading))
    leading_transpose = sp.simplify(transpose_symmetric.subs(leading))
    expected_leading_hermitian = -(
        lambda_n + lambda_m + 2 * u_terminal
    ) ** 2 / 8
    expected_leading_transpose = -(
        lambda_n - lambda_m
    ) ** 2 / 8 - sp.I * u_x / 2
    if sp.simplify(
        leading_hermitian - expected_leading_hermitian
    ) != 0:
        raise RuntimeError("leading Hermitian kernel failed")
    if sp.simplify(
        leading_transpose - expected_leading_transpose
    ) != 0:
        raise RuntimeError("leading transpose kernel failed")

    return {
        "carrier_polynomials": (
            "For ell=log(n), lambda=L-ell=log(N/n), put "
            "C(ell)=1+rho_1*ell+rho_2*ell^2, "
            "D(ell)=rho_(1,x)*ell+rho_(2,x)*ell^2, "
            "R(ell)=s_*'*lambda+i*b*u_N, "
            "Q(ell)=(chi_N+s_*'*lambda)C(ell)+D(ell), "
            "R_x(ell)=s_*''*lambda+i(b_x*u_N+b*u_(N,x)), and "
            "N(ell)=R(ell)Q(ell)+R_x(ell)C(ell)."
        ),
        "observation_sums": (
            "For w_n=(eta/S_a)exp[t(log n)^2/4-s_*log n], "
            "mathscr V=sum C_nw_n, mathscr Q=sum Q_nw_n, "
            "mathscr A=sum R_nC_nw_n, mathscr N=sum N_nw_n; "
            "the real observations are their real parts."
        ),
        "radial_defect": (
            "R(ell)C(ell)-Q(ell)="
            "[i*b*u_N-chi_N]C(ell)-D(ell)."
        ),
        "endpoint_kernel": (
            "L_n=X_T*N_n+A_(T,x)C_n-A_TQ_n-X_(T,x)R_nC_n="
            "(X_TR_n-A_T)Q_n+"
            "[X_TR_(n,x)+A_(T,x)-X_(T,x)R_n]C_n."
        ),
        "unsymmetrized_kernels": (
            "H_(n,m)=C_n*conj(N_m)-R_nC_n*conj(Q_m)="
            "C_n[(conj(R_m)-R_n)conj(Q_m)+conj(R_(m,x))conj(C_m)]; "
            "T_(n,m)=C_nN_m-R_nC_nQ_m="
            "C_n[(R_m-R_n)Q_m+R_(m,x)C_m]."
        ),
        "full_two_carrier_identity": (
            "h^2*Phi_B=Re sum_(n<=B)L_nw_n+"
            "(1/2)Re sum_(n,m<=B)[H_(n,m)w_nconj(w_m)+"
            "T_(n,m)w_nw_m]."
        ),
        "hermitian_symmetrization": (
            "H^+_(n,m)=(H_(n,m)+conj(H_(m,n)))/2="
            "(1/2){[conj(R_m)-R_n][C_nconj(Q_m)-conj(C_m)Q_n]+"
            "C_nconj(C_m)[conj(R_(m,x))+R_(n,x)]}."
        ),
        "transpose_symmetrization": (
            "For Delta=lambda_m-lambda_n, T^+_(n,m)="
            "(T_(n,m)+T_(m,n))/2=(1/2){s_*'^2Delta^2C_nC_m+"
            "s_*'Delta(C_nD_m-C_mD_n)+"
            "C_nC_m[R_(n,x)+R_(m,x)]}."
        ),
        "common_rate_cancellation": (
            "T^+ is independent of chi_N. In H^+, "
            "C_nconj(Q_m)-conj(C_m)Q_n="
            "C_nconj(C_m)[conj(chi_N)-chi_N+"
            "conj(s_*')lambda_m-s_*'lambda_n]+"
            "C_nconj(D_m)-conj(C_m)D_n, so Re(chi_N) cancels."
        ),
        "leading_q1_kernels": (
            "At C=1,D=0,s_*'=-i/2,s_*''=0,b=-1/2,b_x=0,"
            "chi_N=i*b*u_N, one has R_n=Q_n=-i*u_n/2 with "
            "u_n=lambda_n+u_N, H^+_(n,m)=-(u_n+u_m)^2/8, and "
            "T^+_(n,m)=-(lambda_n-lambda_m)^2/8-i*u_(N,x)/2."
        ),
    }


def finite_sum_certificate() -> dict:
    i = sp.I
    log_n = sp.Integer(5)
    rho_1 = (1 + i) / 10
    rho_2 = (-1 + 2 * i) / 50
    rho_1_x = (2 - i) / 100
    rho_2_x = (1 + i) / 200
    s_prime = sp.Rational(1, 20) - i / 2
    s_second = sp.Rational(1, 300) + i / 400
    chi_n = sp.Rational(1, 7) + 2 * i / 9
    b = -sp.Rational(1, 2)
    b_x = sp.Rational(1, 400)
    u_terminal = sp.Rational(1, 5)
    u_x = sp.Rational(1, 100)
    ells = (sp.Integer(0), sp.Integer(1), sp.Integer(3))
    weights = (1 + i / 2, -sp.Rational(2, 3) + i / 5, sp.Rational(3, 7) - 2 * i / 9)

    def lam(ell: sp.Expr) -> sp.Expr:
        return log_n - ell

    def c(ell: sp.Expr) -> sp.Expr:
        return 1 + rho_1 * ell + rho_2 * ell**2

    def d(ell: sp.Expr) -> sp.Expr:
        return rho_1_x * ell + rho_2_x * ell**2

    def r(ell: sp.Expr) -> sp.Expr:
        return s_prime * lam(ell) + i * b * u_terminal

    def rx(ell: sp.Expr) -> sp.Expr:
        return s_second * lam(ell) + i * (b_x * u_terminal + b * u_x)

    def q(ell: sp.Expr) -> sp.Expr:
        return (chi_n + s_prime * lam(ell)) * c(ell) + d(ell)

    def nn(ell: sp.Expr) -> sp.Expr:
        return r(ell) * q(ell) + rx(ell) * c(ell)

    complex_v = sum(c(ell) * weight for ell, weight in zip(ells, weights))
    complex_q = sum(q(ell) * weight for ell, weight in zip(ells, weights))
    complex_a = sum(r(ell) * c(ell) * weight for ell, weight in zip(ells, weights))
    complex_n = sum(nn(ell) * weight for ell, weight in zip(ells, weights))
    bulk_direct = sp.re(complex_v) * sp.re(complex_n) - sp.re(complex_a) * sp.re(complex_q)

    double_sum = 0
    for ell_left, weight_left in zip(ells, weights):
        for ell_right, weight_right in zip(ells, weights):
            hermitian = (
                c(ell_left) * sp.conjugate(nn(ell_right))
                - r(ell_left) * c(ell_left) * sp.conjugate(q(ell_right))
            )
            transpose = (
                c(ell_left) * nn(ell_right)
                - r(ell_left) * c(ell_left) * q(ell_right)
            )
            double_sum += (
                hermitian * weight_left * sp.conjugate(weight_right)
                + transpose * weight_left * weight_right
            )
    bulk_kernel = sp.re(double_sum) / 2
    if sp.simplify(bulk_direct - bulk_kernel) != 0:
        raise RuntimeError("finite exact two-carrier audit failed")

    h = sp.symbols("h", positive=True, real=True)
    positive_witness = sp.simplify(
        -sp.Rational(1, 4)
        * (1 - sp.Rational(1, 2))
        * (h**2 - sp.Rational(1, 2) * (2 * h) ** 2)
    )
    negative_witness = -h**2 / 4
    if positive_witness != h**2 / 8 or negative_witness != -h**2 / 4:
        raise RuntimeError("leading pair-sign guard failed")

    h_cap = sp.Rational(1, 72000000000)
    if not 16892 * h_cap**2 < 1:
        raise RuntimeError("physical chi_N bound comparison failed")

    return {
        "carrier_count": len(ells),
        "exact_bulk_difference": str(sp.simplify(bulk_direct - bulk_kernel)),
        "physical_chi_identity": (
            "At n=N, hat(z_N)=C_Nw_N=-Q(p)exp(W_theta), "
            "Q_x/Q=-i*h*theta/2, and C_(N,x)/C_N=delta_N; hence "
            "chi_N=-i*h*theta/2+W_(theta,x)-delta_N."
        ),
        "physical_chi_bound": (
            "Since delta_N=d_(N,x)/(1+d_N)-d_(1,x)/(1+d_1), "
            "the bounds |d_n|<1/2 and |d_(n,x)|<4223/x^2 give "
            "|delta_N|<16892/x^2. Together with "
            "|W_(theta,x)|<2h^2, x>a^2=h^-2, and 16892h^2<1, "
            "this gives |chi_N+i*h*theta/2|<3h^2."
        ),
        "physical_radial_defect_bound": (
            "For u_N=-log(1-h*theta), |s_*'+i/2|<h^2/6000 "
            "and 0<=u_N-h*theta<h^2 imply "
            "|i*b*u_N-chi_N|<4h^2. The full coefficient defect also "
            "contains -D(ell), which must remain in the bulk estimate."
        ),
        "positive_pair_guard": "u_1=h,u_2=2h,w_1=1,w_2=-1/2 gives P_bulk=h^2/8",
        "negative_single_guard": "u_1=h,w_1=1 gives P_bulk=-h^2/4",
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use the fixed physical q=1, L>=50 terminal/bulk chart with "
            "B=N-M-1. Write ell_n=log(n), lambda_n=log(N/n), and use "
            "the correction-free normalized carrier w_n in G_r=sum "
            "ell_n^r w_n. The centered bulk rate is chi_N and is distinct "
            "from the established endpoint factor kappa_N."
        ),
        "kernel_interpretation": (
            "The Hermitian kernel multiplies w_n conjugate(w_m) and carries "
            "phase differences. The transpose kernel multiplies w_nw_m and "
            "carries phase sums. Both occur in the real projected current; "
            "controlling only one correlation family is incomplete."
        ),
        "two_variable_gain": (
            "Expanding the rank-four moment quadratic directly produces one "
            "ordered two-carrier sum, not a four-variable product of two "
            "separately absolute-valued Mangoldt sums. The degree-four "
            "polynomials may be rejoined with the exact divisor identities "
            "while preserving both oscillatory kernels."
        ),
        "physical_gain": (
            "The terminal phase anchor makes chi_N=O(h), with an O(h^2) "
            "error around -i*h*theta/2, and makes i*b*u_N-chi_N=O(h^2). "
            "Thus the real observations A and Q have a source-specific "
            "near-coincidence. This is genuine structure absent from the "
            "generic inertia witness, but no bound on the summed D(ell) "
            "defect or the two oscillatory kernels is yet proved."
        ),
        "leading_sign_guard": (
            "In the ideal q=1 leading kernel every single positive-distance "
            "carrier is clockwise, but two real opposite-phase carriers "
            "with distances h and 2h and amplitudes 1 and 1/2 give positive "
            "bulk current h^2/8. This is an algebraic kernel guard, not an "
            "attained Xi configuration or aggregate counterexample."
        ),
        "sharp_target": (
            "Prove that the endpoint-linear sum plus half the real part of "
            "the complete Hermitian and transpose two-carrier sums is at "
            "most h^2/400, preferably h^2/800, for the actual Xi carriers."
        ),
        "type_ii_handoff": (
            "Apply balanced Type-I/II or Vaughan decomposition jointly to "
            "the phase-difference and phase-sum kernels. Retain the complete "
            "square, both wings, correction derivative D(ell), terminal "
            "anchor, and adjacent-cutoff transport. A theorem for only the "
            "transpose-symmetric Mangoldt form cannot close Phi_B."
        ),
        "proof_boundary": (
            "This proves an exact endpoint-linear plus Hermitian/transpose "
            "two-carrier decomposition, common-rate cancellations, the "
            "physical chi_N anchor bound, leading kernels, and a pair-sign "
            "guard. It proves no signed Type-I/II or Vaughan estimate, upper "
            "bound on Phi_B, contact exclusion, retained aggregate sign, Xi "
            "residual transfer, Q209, cofinal descendant theorem, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
    }


def build_rows(exact: dict, symbolic: dict, finite: dict) -> list[KernelRow]:
    return [
        KernelRow("ttkc_01_domain", "fixed-chart domain", "proved", "The kernel uses the exact certified terminal/bulk chart and disambiguated rate.", exact["domain"], "No cutoff derivative is inserted."),
        KernelRow("ttkc_02_polynomials", "carrier polynomials", "proved", "Five moments collapse to four explicit degree-at-most-four carrier polynomials.", symbolic["carrier_polynomials"], "Every correction derivative is retained."),
        KernelRow("ttkc_03_observations", "observation sums", "proved", "The four real observations are real parts of four carrier sums.", symbolic["observation_sums"], "No real projection is divided out."),
        KernelRow("ttkc_04_defect", "radial defect", "proved", "The A-Q coefficient defect loses every logarithmic-distance term exactly.", symbolic["radial_defect"], finite["physical_radial_defect_bound"]),
        KernelRow("ttkc_05_endpoint", "endpoint-linear kernel", "proved", "All tail/bulk cross terms form one exact endpoint-linear carrier sum.", symbolic["endpoint_kernel"], "The recurrent endpoint remains composed."),
        KernelRow("ttkc_06_unsymmetrized", "two-carrier kernels", "proved", "The bulk quadratic is one Hermitian plus one transpose ordered pair sum.", symbolic["unsymmetrized_kernels"], symbolic["full_two_carrier_identity"]),
        KernelRow("ttkc_07_hermitian", "Hermitian symmetrization", "proved", "The phase-difference kernel has an exact Hermitian symmetrization.", symbolic["hermitian_symmetrization"], exact["kernel_interpretation"]),
        KernelRow("ttkc_08_transpose", "transpose symmetrization", "proved", "The phase-sum kernel has an exact transpose symmetrization.", symbolic["transpose_symmetrization"], exact["kernel_interpretation"]),
        KernelRow("ttkc_09_rates", "common-rate cancellation", "proved", "The transpose kernel loses chi_N and the Hermitian kernel loses Re(chi_N).", symbolic["common_rate_cancellation"], "Im(chi_N) remains physical data."),
        KernelRow("ttkc_10_anchor", "physical terminal anchor", "proved", "The physical q=1 terminal carrier constrains chi_N sharply.", finite["physical_chi_identity"] + " " + finite["physical_chi_bound"], "This bound is q=1 and uses the certified terminal chart."),
        KernelRow("ttkc_11_leading", "leading q=1 kernels", "proved", "The ideal source has closed Hermitian and transpose kernels.", symbolic["leading_q1_kernels"], "Entrywise signs do not imply aggregate sign."),
        KernelRow("ttkc_12_guard", "pair-sign guard", "proved", "Single-carrier clockwise signs do not survive a generic two-carrier join.", finite["positive_pair_guard"] + "; " + finite["negative_single_guard"], exact["leading_sign_guard"]),
        KernelRow("ttkc_13_target", "joint signed target", "open", "Both oscillatory kernel families and the endpoint must satisfy one inequality.", exact["sharp_target"], "No signed constant is proved."),
        KernelRow("ttkc_14_handoff", "Type-I/II handoff", "open", "The two-variable kernel is ready for a joint balanced decomposition.", exact["type_ii_handoff"], exact["proof_boundary"]),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    finite = payload["finite_sum_certificate"]
    return f"""# Terminal-Tail Five-Moment Two-Carrier Kernel Reduction

Date: 2026-07-31

Status: exact endpoint/Hermitian/transpose reduction with `0 signed Type-I/II
bounds`. This is not a proof artifact; `Phi_B` and RH remain open.

## Carrier Polynomials

{exact['domain']}

```text
{symbolic['carrier_polynomials']}
```

```text
{symbolic['observation_sums']}
```

The exact radial-subtracted defect is

```text
{symbolic['radial_defect']}
```

## Endpoint And Pair Kernels

```text
{symbolic['endpoint_kernel']}

{symbolic['unsymmetrized_kernels']}

{symbolic['full_two_carrier_identity']}
```

The two exact symmetrizations are

```text
{symbolic['hermitian_symmetrization']}

{symbolic['transpose_symmetrization']}
```

{exact['kernel_interpretation']}

{symbolic['common_rate_cancellation']}

{exact['two_variable_gain']}

## Physical Anchor

```text
{finite['physical_chi_identity']}

{finite['physical_chi_bound']}

{finite['physical_radial_defect_bound']}
```

{exact['physical_gain']}

## Leading Kernel And Guard

```text
{symbolic['leading_q1_kernels']}
```

```text
{finite['negative_single_guard']}
{finite['positive_pair_guard']}
```

{exact['leading_sign_guard']}

## Open Arithmetic Target

```text
{exact['sharp_target']}
```

{exact['type_ii_handoff']}

## Boundary

No new `pi` occurs in the kernel decomposition. The `pi` in
`Q_x/Q=-i*h*theta/2` comes from the inherited Riemann--Siegel terminal phase.

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    payloads = load_sources()
    source = source_audit(payloads)
    symbolic = symbolic_certificate()
    finite = finite_sum_certificate()
    exact = exact_payload()
    rows = build_rows(exact, symbolic, finite)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact terminal-tail endpoint-linear and Hermitian/transpose "
            "two-carrier Phi_B kernel; signed Type-I/II bound open; no "
            "retained aggregate sign, Xi-level theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "finite_sum_certificate": finite,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "carrier_polynomials": 4,
            "endpoint_linear_kernels": 1,
            "hermitian_kernels": 1,
            "transpose_kernels": 1,
            "exact_kernel_symmetrizations": 2,
            "common_rate_cancellations": 2,
            "physical_chi_anchor_identities": 1,
            "physical_chi_bounds": 1,
            "leading_q1_kernels": 2,
            "pair_sign_guards": 1,
            "signed_type_ii_targets": 1,
            "signed_type_ii_bounds": 0,
            "phi_b_bounds": 0,
            "xi_level_current_theorems": 0,
        },
        "rows": [asdict(row) for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "wrote five-moment two-carrier kernel reduction: "
        f"{payload['counts']['rows']} rows, 1 Hermitian kernel, "
        "1 transpose kernel, 1 physical chi_N bound, "
        "0 signed Type-I/II bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
