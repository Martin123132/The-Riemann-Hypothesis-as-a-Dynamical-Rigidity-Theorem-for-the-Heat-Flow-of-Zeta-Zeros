#!/usr/bin/env python3
"""Compile explicit arbitrary-N m=9 bounds for the modular theta tail."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import sys

try:
    import psutil

    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception:
    pass


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb
import sympy as sp


STEM = "jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate"
DEFAULT_OUT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
MODULAR_SOURCE = RESULT_ROOT / "jensen_window_pf_newman_theta_modular_blend_gate.json"
ENVELOPE_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.json"
)
ARBITRARY_N_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate.json"
)
OUTER_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.json"
)
FULL_BUDGET_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_full_derivative_budget_certificate.json"
)
DATE = "2026-07-24"
PRECISION_BITS = 256
MAX_ORDER = 9
M_ORDER = 9
T_EXACT = Fraction(1, 5)
MIN_FIRST_OMITTED = 3
K_VALUES = (3, 4, 5, 8, 11, 17, 26, 65)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def rational_arb(value: Fraction | sp.Rational | int) -> arb:
    if isinstance(value, int):
        return arb(value)
    return arb(int(value.p if isinstance(value, sp.Rational) else value.numerator)) / int(
        value.q if isinstance(value, sp.Rational) else value.denominator
    )


def arb_power_rational(base: arb, exponent: Fraction | sp.Rational) -> arb:
    exponent_ball = rational_arb(exponent)
    return (exponent_ball * base.log()).exp()


def serialize_positive(value: arb) -> dict:
    if not value.is_finite() or value.lower() <= 0:
        raise RuntimeError(f"nonpositive or nonfinite enclosure: {value}")
    logarithm = value.log() / arb(10).log()
    return {
        "enclosure": value.str(70, more=True),
        "upper": value.upper().str(70),
        "log10_enclosure": logarithm.str(55, more=True),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def polynomial_rows() -> list[dict]:
    x = sp.symbols("X", nonnegative=True)
    current = 2 * x - 3
    rows: list[dict] = []
    for order in range(MAX_ORDER + 1):
        polynomial = sp.Poly(sp.expand(current), x)
        if polynomial.degree() != order + 1:
            raise RuntimeError(f"kernel polynomial degree drift at q={order}")
        rows.append(
            {
                "q": order,
                "P_q": str(polynomial.as_expr()),
                "degree": polynomial.degree(),
                "A_q": int(sum(abs(value) for value in polynomial.all_coeffs())),
            }
        )
        current = sp.expand(
            (5 - 4 * x) * polynomial.as_expr()
            + 4 * x * sp.diff(polynomial.as_expr(), x)
        )
    return rows


def laurent_data(expression: sp.Expr, y: sp.Symbol) -> tuple[sp.Rational, int]:
    coefficient_norm = sp.Rational(0)
    maximum_degree: int | None = None
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient, exponent = term.as_coeff_exponent(y)
        if coefficient.has(y) or not exponent.is_integer:
            raise RuntimeError(f"not a Laurent monomial: {term}")
        coefficient_norm += abs(coefficient)
        degree = int(exponent)
        maximum_degree = degree if maximum_degree is None else max(maximum_degree, degree)
    if maximum_degree is None:
        raise RuntimeError("empty switch Laurent polynomial")
    return coefficient_norm, maximum_degree


def switch_derivative_rows() -> list[dict]:
    y = sp.symbols("y", positive=True)
    z = sp.Rational(3, 2) * (y - y**-1)
    z_prime = 4 * y * sp.diff(z, y)
    current = sp.expand(-z_prime)
    rows = [
        {
            "a": 0,
            "laurent_polynomial": None,
            "coefficient_norm": "1/2",
            "maximum_degree": 0,
            "normalization": "1",
            "bound": "0<s(u)<=0.5*exp(-h(u))",
        }
    ]
    for order in range(1, MAX_ORDER + 1):
        coefficient_norm, maximum_degree = laurent_data(current, y)
        if maximum_degree != 2 * order - 1:
            raise RuntimeError(f"switch degree drift at a={order}")
        rows.append(
            {
                "a": order,
                "laurent_polynomial": str(current),
                "coefficient_norm": str(coefficient_norm),
                "maximum_degree": maximum_degree,
                "normalization": "1/sqrt(pi)",
                "bound": (
                    "|s^(a)(u)|<=C_a*y^h_a*exp(-h(u))/sqrt(pi), "
                    "h_a=2a-1"
                ),
            }
        )
        current = sp.expand(4 * y * sp.diff(current, y) - 2 * z * z_prime * current)
    return rows


def forward_tail(
    p: int,
    kernel_power: int,
    beta: sp.Rational,
    first: int,
) -> arb:
    pi = arb.pi()
    degree = 2 * kernel_power - 2
    gamma_power = beta - sp.Rational(3, 4)
    if gamma_power < 0:
        gamma_power = sp.Rational(0)
    denominator_shift = gamma_power + sp.Rational(1, 80)
    theta = 1 - rational_arb(denominator_shift) / (pi * first**2)
    if theta.lower() <= 0:
        raise RuntimeError("forward denominator margin is not positive")
    rho = (
        arb(degree) / first - pi * (2 * first + 1)
    ).exp()
    if rho.upper() >= 1:
        raise RuntimeError("forward arithmetic ratio is not below one")
    return (
        factorial(p)
        / arb(4)
        * (arb(1) / 80).exp()
        * pi ** (kernel_power - 1)
        / theta
        * first**degree
        * (-pi * first**2).exp()
        / (1 - rho)
    )


def reflected_compact_tail(
    p: int,
    kernel_power: int,
    beta: sp.Rational,
    first: int,
) -> arb:
    pi = arb.pi()
    sqrt_two = arb(2).sqrt()
    length = arb(2).log() / 8
    degree = 2 * kernel_power
    rho = (
        arb(degree) / first
        - pi / sqrt_two * (2 * first + 1)
    ).exp()
    if rho.upper() >= 1:
        raise RuntimeError("compact reflected arithmetic ratio is not below one")
    beta_factor = (
        rational_arb(beta / 2) * arb(2).log()
    ).exp()
    return (
        length ** (p + 1)
        * (length**2 / 5).exp()
        * beta_factor
        * pi**kernel_power
        * first**degree
        * (-pi * first**2 / sqrt_two).exp()
        / (1 - rho)
    )


def safe_phase_constant() -> arb:
    pi = arb.pi()
    c_zero = (
        arb(3)
        / 4
        * arb_power_rational(arb(3), Fraction(2, 3))
        * arb_power_rational(pi, Fraction(2, 3))
    )
    return c_zero / 4


def stretched_exponential_sum(kernel_power: int, first: int) -> arb:
    exponent = 2 * kernel_power
    c_safe = safe_phase_constant()
    boundary = arb(first - 1)
    mode = arb_power_rational(
        arb(3 * exponent) / (4 * c_safe),
        Fraction(3, 4),
    )
    if mode.upper() <= boundary:
        supremum_point = boundary
    elif mode.lower() >= boundary:
        supremum_point = mode
    else:
        raise RuntimeError("stretched-tail mode decision is not rigorous")
    supremum = (
        supremum_point**exponent
        * (
            -c_safe
            * arb_power_rational(supremum_point, Fraction(4, 3))
        ).exp()
    )
    shape = Fraction(3 * (exponent + 1), 4)
    gamma_argument = (
        c_safe * arb_power_rational(boundary, Fraction(4, 3))
    )
    integral = (
        arb(3)
        / 4
        * arb_power_rational(c_safe, -shape)
        * gamma_argument.gamma_upper(rational_arb(shape))
    )
    return supremum + integral


def reflected_large_tail(
    p: int,
    kernel_power: int,
    beta: sp.Rational,
    first: int,
) -> arb:
    pi = arb.pi()
    quadratic = sp.Rational(27, 64)
    remaining_quadratic = sp.Rational(27, 128)
    heat_linear = sp.Rational(1, 80)
    delta = beta - sp.Rational(3, 4)
    shape = (delta + 1) / 2
    if shape <= 0:
        raise RuntimeError("large reflected gamma shape is not positive")
    young_constant = (
        rational_arb(heat_linear**2 / (2 * quadratic))
    ).exp()
    gamma_argument = rational_arb(2 * remaining_quadratic)
    y_integral = (
        factorial(p)
        / arb(8)
        * young_constant
        * arb_power_rational(rational_arb(remaining_quadratic), -shape)
        * gamma_argument.gamma_upper(rational_arb(shape))
    )
    return (
        y_integral
        * pi**kernel_power
        * stretched_exponential_sum(kernel_power, first)
    )


def monomial_budget(
    p: int,
    derivative_order: int,
    first: int,
    polynomials: list[dict],
    switches: list[dict],
) -> dict[str, arb]:
    components = {
        "forward": arb(0),
        "reflected_compact": arb(0),
        "reflected_large": arb(0),
    }
    pi = arb.pi()
    for switch_order in range(derivative_order + 1):
        kernel_order = derivative_order - switch_order
        polynomial_constant = arb(polynomials[kernel_order]["A_q"])
        kernel_power = kernel_order + 2
        multiplicity = comb(derivative_order, switch_order)

        if switch_order == 0:
            forward_constant = multiplicity * polynomial_constant
            forward_beta = sp.Rational(kernel_order) + sp.Rational(9, 4)
            reflected_constant = (
                multiplicity
                * sp.Rational(1, 2)
                * polynomials[kernel_order]["A_q"]
                * 2 ** (kernel_order + 1)
            )
            reflected_constant_ball = rational_arb(reflected_constant)
            reflected_beta = sp.Rational(0)
        else:
            switch_constant = rational_arb(
                sp.Rational(switches[switch_order]["coefficient_norm"])
            ) / pi.sqrt()
            switch_degree = switches[switch_order]["maximum_degree"]
            forward_constant = (
                multiplicity * switch_constant * polynomial_constant
            )
            forward_beta = (
                sp.Rational(switch_degree)
                + kernel_order
                + sp.Rational(9, 4)
            )
            reflected_constant_ball = (
                multiplicity
                * switch_constant
                * polynomial_constant
                * 2 ** (kernel_order + 1)
            )
            reflected_beta = max(
                sp.Rational(switch_degree) - sp.Rational(5, 4),
                sp.Rational(0),
            )

        components["forward"] += forward_constant * forward_tail(
            p,
            kernel_power,
            forward_beta,
            first,
        )
        components["reflected_compact"] += (
            reflected_constant_ball
            * reflected_compact_tail(
                p,
                kernel_power,
                reflected_beta,
                first,
            )
        )
        components["reflected_large"] += (
            reflected_constant_ball
            * reflected_large_tail(
                p,
                kernel_power,
                reflected_beta,
                first,
            )
        )
    components["total"] = sum(components.values(), arb(0))
    return components


def heat_terms(order: int) -> list[tuple[int, Fraction]]:
    terms: list[tuple[int, Fraction]] = []
    for ell in range(order // 2 + 1):
        power = order - 2 * ell
        coefficient = Fraction(
            factorial(order) * 2**power,
            factorial(ell)
            * factorial(power)
            * 5 ** (order - ell),
        )
        terms.append((power, coefficient))
    return terms


def empty_components() -> dict[str, arb]:
    return {
        "forward": arb(0),
        "reflected_compact": arb(0),
        "reflected_large": arb(0),
        "total": arb(0),
    }


def add_scaled(
    destination: dict[str, arb],
    source: dict[str, arb],
    scale: Fraction | int,
) -> None:
    scale_ball = rational_arb(scale)
    for key in destination:
        destination[key] += scale_ball * source[key]


def serialize_components(components: dict[str, arb]) -> dict:
    return {key: serialize_positive(value) for key, value in components.items()}


def compose_budget(
    first: int,
    polynomials: list[dict],
    switches: list[dict],
) -> dict:
    cache: dict[tuple[int, int], dict[str, arb]] = {}

    def integral(power: int, order: int) -> dict[str, arb]:
        key = (power, order)
        if key not in cache:
            cache[key] = monomial_budget(
                power,
                order,
                first,
                polynomials,
                switches,
            )
        return cache[key]

    value = empty_components()
    first_jet_u = empty_components()
    first_jet_lower = empty_components()

    for heat_order in range(M_ORDER + 1):
        derivative_order = M_ORDER - heat_order
        leibniz = comb(M_ORDER, heat_order)
        for power, heat_coefficient in heat_terms(heat_order):
            scale = leibniz * heat_coefficient
            add_scaled(value, integral(power, derivative_order), scale)
            add_scaled(
                first_jet_u,
                integral(power + 1, derivative_order),
                scale,
            )

    for heat_order in range(M_ORDER):
        derivative_order = M_ORDER - 1 - heat_order
        leibniz = M_ORDER * comb(M_ORDER - 1, heat_order)
        for power, heat_coefficient in heat_terms(heat_order):
            add_scaled(
                first_jet_lower,
                integral(power, derivative_order),
                leibniz * heat_coefficient,
            )

    first_jet = {
        key: first_jet_u[key] + first_jet_lower[key]
        for key in first_jet_u
    }
    return {
        "N": first - 1,
        "K": first,
        "monomial_integrals_reused": len(cache),
        "d0_m9_t1_5": serialize_components(value),
        "d1_u_m9_t1_5": serialize_components(first_jet_u),
        "d1_lower_m9_t1_5": serialize_components(first_jet_lower),
        "d1_m9_t1_5": serialize_components(first_jet),
    }


def witness_rows(polynomials: list[dict], switches: list[dict]) -> list[dict]:
    rows = [compose_budget(first, polynomials, switches) for first in K_VALUES]
    for previous, current in zip(rows, rows[1:]):
        for field in ("d0_m9_t1_5", "d1_m9_t1_5"):
            previous_value = arb(previous[field]["total"]["enclosure"])
            current_value = arb(current[field]["total"]["enclosure"])
            if current_value.upper() >= previous_value.lower():
                raise RuntimeError(f"{field} witness bounds are not strictly decreasing")
    return rows


def condition_diagnostics(
    polynomials: list[dict],
    switches: list[dict],
) -> dict:
    forward_betas: list[sp.Rational] = []
    reflected_betas: list[sp.Rational] = []
    for order in range(MAX_ORDER + 1):
        for switch_order in range(order + 1):
            kernel_order = order - switch_order
            if switch_order == 0:
                forward_betas.append(
                    sp.Rational(kernel_order) + sp.Rational(9, 4)
                )
                reflected_betas.append(sp.Rational(0))
            else:
                degree = switches[switch_order]["maximum_degree"]
                forward_betas.append(
                    sp.Rational(degree)
                    + kernel_order
                    + sp.Rational(9, 4)
                )
                reflected_betas.append(
                    max(
                        sp.Rational(degree) - sp.Rational(5, 4),
                        sp.Rational(0),
                    )
                )
    max_forward_beta = max(forward_betas)
    max_gamma_power = max_forward_beta - sp.Rational(3, 4)
    pi = arb.pi()
    theta_min = (
        1
        - rational_arb(max_gamma_power + sp.Rational(1, 80))
        / (pi * MIN_FIRST_OMITTED**2)
    )
    forward_rho_max = (
        arb(2 * (MAX_ORDER + 2) - 2) / MIN_FIRST_OMITTED
        - pi * (2 * MIN_FIRST_OMITTED + 1)
    ).exp()
    reflected_rho_max = (
        arb(2 * (MAX_ORDER + 2)) / MIN_FIRST_OMITTED
        - pi / arb(2).sqrt() * (2 * MIN_FIRST_OMITTED + 1)
    ).exp()
    if theta_min.lower() <= 0:
        raise RuntimeError("global forward denominator diagnostic failed")
    if max(forward_rho_max.upper(), reflected_rho_max.upper()) >= 1:
        raise RuntimeError("global arithmetic ratio diagnostic failed")
    return {
        "max_forward_beta_exact": str(max_forward_beta),
        "max_forward_gamma_power_exact": str(max_gamma_power),
        "max_reflected_beta_exact": str(max(reflected_betas)),
        "minimum_forward_theta_at_K3": serialize_positive(theta_min),
        "maximum_forward_rho_at_K3": serialize_positive(forward_rho_max),
        "maximum_compact_reflected_rho_at_K3": serialize_positive(
            reflected_rho_max
        ),
        "c0_exact": "3*3^(2/3)*pi^(2/3)/4",
        "c_safe_exact": "3*3^(2/3)*pi^(2/3)/16",
        "c_safe": serialize_positive(safe_phase_constant()),
        "switch_coefficient_norms": [
            {
                "a": row["a"],
                "C_a": row["coefficient_norm"],
                "h_a": row["maximum_degree"],
            }
            for row in switches
        ],
        "kernel_coefficient_norms": [
            {"q": row["q"], "A_q": row["A_q"]} for row in polynomials
        ],
    }


def source_audit(
    polynomials: list[dict],
    switches: list[dict],
) -> dict:
    modular = load(MODULAR_SOURCE)
    envelope = load(ENVELOPE_SOURCE)
    arbitrary = load(ARBITRARY_N_SOURCE)
    outer = load(OUTER_SOURCE)
    full = load(FULL_BUDGET_SOURCE)

    tail_identity = (
        modular.get("exact", {})
        .get("decaying_tail_enclosure", {})
        .get("tail_kernel")
    )
    if tail_identity != "r_N(u)=sum_(n>N)b_n(u)=Phi(u)-sum_(n<=N)b_n(u)":
        raise RuntimeError("modular tail identity drifted")
    safe_exponent = (
        envelope.get("exact", {})
        .get("derivative_tail_theorem", {})
        .get("safe_exponent")
    )
    if safe_exponent != "c_*=c_0/4":
        raise RuntimeError("reflected safe exponent drifted")
    if arbitrary.get("polynomials") != polynomials:
        raise RuntimeError("arbitrary-N kernel polynomials drifted")

    outer_switches = outer.get("switch_derivatives", [])
    for expected, observed in zip(switches, outer_switches, strict=True):
        for key in (
            "a",
            "laurent_polynomial",
            "coefficient_norm",
            "maximum_degree",
            "normalization",
        ):
            if expected[key] != observed.get(key):
                raise RuntimeError(f"outer switch row drifted at a={expected['a']}")
    parameters = full.get("parameters", {})
    if parameters.get("m_order") != M_ORDER or parameters.get("t_cap_exact") != "1/5":
        raise RuntimeError("full derivative-budget parameters drifted")

    return {
        "modular_source_sha256": file_hash(MODULAR_SOURCE),
        "envelope_source_sha256": file_hash(ENVELOPE_SOURCE),
        "arbitrary_n_source_sha256": file_hash(ARBITRARY_N_SOURCE),
        "outer_source_sha256": file_hash(OUTER_SOURCE),
        "full_budget_source_sha256": file_hash(FULL_BUDGET_SOURCE),
        "tail_identity": tail_identity,
        "safe_exponent": safe_exponent,
        "m_order": M_ORDER,
        "t_cap_exact": "1/5",
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    polynomials = polynomial_rows()
    switches = switch_derivative_rows()
    audit = source_audit(polynomials, switches)
    diagnostics = condition_diagnostics(polynomials, switches)
    witnesses = witness_rows(polynomials, switches)

    rows = [
        GateRow(
            id="ntsdecg_01_exact_kernel_recurrences",
            role="exact_identity",
            readiness="proved",
            claim=(
                "Kernel and positive-order switch derivatives are generated "
                "by exact polynomial and Laurent-polynomial recurrences."
            ),
            formula=(
                "D^r phi_n=Q*y^(5/4)*P_r(Qy)*exp(-Qy); "
                "s^(a)=L_a(y)*exp(-h)/sqrt(pi), a>=1"
            ),
            proof_boundary="Exact through derivative order nine.",
        ),
        GateRow(
            id="ntsdecg_02_direct_tail_leibniz_envelope",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Every derivative of b_n=(1-s)phi_n+s*phi_n(-u) has an "
                "explicit positive four-template Leibniz envelope."
            ),
            formula=(
                "D^k b_n=sum_(a=0)^k binom(k,a)"
                "[D^a(1-s)D^(k-a)phi_n+D^a(s)D^(k-a)phi_n(-u)]"
            ),
            proof_boundary=(
                "Absolute values discard cancellation but introduce no "
                "finite retained-count cap."
            ),
        ),
        GateRow(
            id="ntsdecg_03_forward_monomial_integral",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "The forward weighted monomial integral is reduced to a "
                "Gaussian arithmetic tail with an explicit denominator."
            ),
            formula=(
                "F_(p,R,beta,K)=p!*exp(1/80)*pi^(R-1)"
                "*K^(2R-2)*exp(-pi*K^2)"
                "/[4*theta*(1-rho_F)]"
            ),
            proof_boundary=(
                "Uses u^p<=p!*exp(u), u^2<=exp(4u)/16, and "
                "log(y)<=y-1."
            ),
            diagnostics={
                "theta": (
                    "1-(beta-3/4+1/80)/(pi*K^2), "
                    "with beta-3/4>=0"
                ),
                "rho_F": "exp((2R-2)/K-pi*(2K+1))",
            },
        ),
        GateRow(
            id="ntsdecg_04_reflected_compact_integral",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "On 1<=y<=sqrt(2), the reflected weighted integral is an "
                "explicit Gaussian arithmetic tail."
            ),
            formula=(
                "C_(p,R,beta,K)=L^(p+1)*exp(L^2/5)*2^(beta/2)"
                "*pi^R*K^(2R)*exp(-pi*K^2/sqrt(2))/(1-rho_C), "
                "L=log(2)/8"
            ),
            proof_boundary="The switch phase is dropped only on this compact interval.",
        ),
        GateRow(
            id="ntsdecg_05_reflected_large_integral",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "On y>=sqrt(2), one quarter of the reflected phase gives "
                "uniform n^(4/3) decay and the remaining phase is integrated "
                "by an upper incomplete gamma function."
            ),
            formula=(
                "exp(-Q/y-h)<=exp(-c_safe*n^(4/3))"
                "*exp(-27*y^2/64), "
                "c_safe=3*3^(2/3)*pi^(2/3)/16"
            ),
            proof_boundary=(
                "The discrete stretched-exponential tail is bounded by "
                "one unimodal supremum plus its directed gamma integral."
            ),
        ),
        GateRow(
            id="ntsdecg_06_arbitrary_n_monomial_theorem",
            role="exact_theorem",
            readiness="proved",
            claim=(
                "For every N>=2, p<=10, and k<=9 used below, the positive-"
                "half-line integral of u^p*exp(u^2/5)*|r_N^(k)| has an "
                "explicit finite directed upper bound."
            ),
            formula=(
                "I_(p,k)(N)<=sum_(a=0)^k binom(k,a)"
                "[forward_(a)+reflected_compact_(a)+reflected_large_(a)]"
            ),
            proof_boundary=(
                "Termwise summation is justified by the displayed positive "
                "summable majorants."
            ),
        ),
        GateRow(
            id="ntsdecg_07_full_m9_d0_d1_compiler",
            role="exact_composition",
            readiness="proved",
            claim=(
                "The explicit monomial bounds compose into arbitrary-N "
                "m=9 value and first-jet budgets at T=1/5."
            ),
            formula=(
                "d0<=sum_b binom(9,b)sum_ell H_(b,ell)I_(b-2ell,9-b); "
                "d1<=sum_b binom(9,b)sum_ell H_(b,ell)I_(b-2ell+1,9-b)"
                "+9*sum_(b<=8)binom(8,b)sum_ell H_(b,ell)I_(b-2ell,8-b)"
            ),
            proof_boundary=(
                "The witness table samples the all-N formula; it is not a "
                "finite calibration substituted for the theorem."
            ),
        ),
        GateRow(
            id="ntsdecg_08_three_quarter_scale_instantiation",
            role="exact_consequence",
            readiness="proved",
            claim=(
                "Substitution N=ceil(kappa*(1+x)^(3/4)) gives an explicit "
                "computable exponential-in-x absolute-error budget."
            ),
            formula=(
                "d_(N,j,9)(1/5)<=poly(N)*exp(-c_safe*(N+1)^(4/3)); "
                "N=ceil(kappa*(1+x)^(3/4))"
            ),
            proof_boundary=(
                "This controls only the omitted tail; it does not lower-bound "
                "the retained value or derivative."
            ),
        ),
        GateRow(
            id="ntsdecg_09_cofinal_separation_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "The next obligation is to compare these all-N upper budgets "
                "with a rigorous retained J/J' lower margin on every adaptive "
                "transition cell."
            ),
            formula=(
                "|J-J_N|<=16*d0/x^5; "
                "|J'-J_N'|<=64*d0/x^6+16*d1/x^5"
            ),
            proof_boundary=(
                "No retained first-jet separation, terminating cofinal cover, "
                "strict Laguerre positivity, Lambda<=0, RH, or Clay-prize "
                "proof is supplied."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "explicit all-N modular switch-tail constants and m=9 "
            "value/first-jet compiler"
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "max_derivative_order": MAX_ORDER,
            "m_order": M_ORDER,
            "t_cap_exact": "1/5",
            "minimum_retained_count": MIN_FIRST_OMITTED - 1,
            "minimum_first_omitted": MIN_FIRST_OMITTED,
            "retained_count_upper_bound": None,
            "witness_K": list(K_VALUES),
        },
        "source_audit": audit,
        "kernel_polynomials": polynomials,
        "switch_derivatives": switches,
        "exact_majorants": {
            "variables": (
                "y=exp(4u), Q=pi*n^2, h=9*sinh(4u)^2, "
                "r_N=sum_(n>N)b_n"
            ),
            "forward_kernel": (
                "|D^r phi_n(u)|<=A_r*Q^(r+2)*y^(r+9/4)*exp(-Qy)"
            ),
            "reflected_kernel": (
                "|D^r phi_n(-u)|<=A_r*2^(r+1)*Q^(r+2)"
                "*y^(-5/4)*exp(-Q/y), n>=3"
            ),
            "switch": (
                "s<=exp(-h)/2; "
                "|s^(a)|<=C_a*y^(2a-1)*exp(-h)/sqrt(pi), a>=1"
            ),
            "elementary_weight_bounds": (
                "u^p<=p!*exp(u); u^2<=exp(4u)/16; "
                "y^g<=exp(g*(y-1)), g>=0"
            ),
            "reflected_split": "1<=y<=sqrt(2) and y>=sqrt(2)",
            "large_phase": (
                "Q/y+h>=c0*n^(4/3); "
                "h>=9y^2/16 for y>=sqrt(2)"
            ),
            "unimodal_sum_lemma": (
                "sum_(n>=K)f(n)<=sup_(x>=K-1)f(x)"
                "+integral_(K-1)^infinity f(x)dx for nonnegative unimodal f"
            ),
            "heat_coefficients": (
                "H_(b,ell)=b!*2^(b-2ell)*(1/5)^(b-ell)"
                "/(ell!*(b-2ell)!)"
            ),
        },
        "condition_diagnostics": diagnostics,
        "budgets": witnesses,
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact replaces the formerly merely effective switch-tail "
            "constant by an explicit directed formula for every N>=2 and "
            "assembles full m=9 d0/d1 budgets at T=1/5. It does not prove a "
            "retained first-jet lower separation, does not certify an unbounded "
            "adaptive transition cover for x>38, and does not prove strict "
            "Laguerre positivity, Lambda<=0, RH, or a Clay-prize result."
        ),
    }


def render_note(artifact: dict) -> str:
    lines = [
        "# Newman Theta Switch-Defect Explicit Constant Gate",
        "",
        f"Date: {DATE}",
        "",
        "Status: explicit all-`N` modular-tail constants and `m=9`",
        "value/first-jet compiler. This is not a retained-separation theorem",
        "and not a proof of `Lambda<=0`, RH, or a Clay-prize result.",
        "",
        "## Direct Tail",
        "",
        "For `N>=2`, set `K=N+1`, `y=exp(4u)`, `Q=pi*n^2`, and",
        "`s=erfc(3*sinh(4u))/2`. The proof bounds the direct tail",
        "",
        "```text",
        "r_N(u)=sum_(n>=K)[(1-s(u))*phi_n(u)+s(u)*phi_n(-u)].",
        "```",
        "",
        "The derivative recurrence and switch Laurent recurrence give four",
        "positive templates for each Leibniz index: forward order zero,",
        "forward positive switch order, reflected order zero, and reflected",
        "positive switch order.",
        "",
        "## Integral Split",
        "",
        "The forward template uses",
        "",
        "```text",
        "u^p<=p!*exp(u),  u^2<=exp(4u)/16,  log(y)<=y-1,",
        "```",
        "",
        "and becomes an explicit Gaussian arithmetic tail. The reflected",
        "template is split at `y=sqrt(2)`. Its compact part is another",
        "Gaussian arithmetic tail. On the large part,",
        "",
        "```text",
        "Q/y+h >= c0*n^(4/3),",
        "h >= 9*y^2/16,",
        "c_safe=c0/4=3*3^(2/3)*pi^(2/3)/16.",
        "```",
        "",
        "One quarter of the phase supplies the discrete stretched",
        "exponential. The remaining three quarters produce a directed upper",
        "incomplete-gamma integral.",
        "",
        "## m=9 Witnesses",
        "",
        "The following rows sample the symbolic all-`N` formula.",
        "",
        "| N | K | log10 d0 upper | log10 d1 upper |",
        "|---:|---:|---:|---:|",
    ]
    for row in artifact["budgets"]:
        lines.append(
            f"| {row['N']} | {row['K']} | "
            f"`{row['d0_m9_t1_5']['total']['log10_enclosure']}` | "
            f"`{row['d1_m9_t1_5']['total']['log10_enclosure']}` |"
        )
    lines.extend(
        [
            "",
            "Every stored number is an Arb enclosure at 256-bit precision.",
            "The independent checker regenerates the kernel and switch",
            "recurrences and recompiles every witness.",
            "",
            "## Remaining Boundary",
            "",
            "The switch-tail constants are no longer an effective-constant",
            "placeholder. The main cofinal obligation is now a rigorous",
            "retained value-or-derivative lower margin across all adaptive",
            "transition cells, strong enough to dominate these upper budgets.",
            "",
        ]
    )
    return "\n".join(lines)


def write_artifact(artifact: dict, out: Path, note: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_artifact(artifact, args.out, args.note)
    print(
        "wrote Newman theta switch-defect explicit constant gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['budgets'])} arbitrary-N witnesses"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
