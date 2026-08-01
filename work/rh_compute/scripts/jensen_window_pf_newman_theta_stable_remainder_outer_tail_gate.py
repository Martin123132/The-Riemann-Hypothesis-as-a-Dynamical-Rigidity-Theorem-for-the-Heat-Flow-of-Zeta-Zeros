#!/usr/bin/env python3
"""Build explicit Arb bounds for the stable remainder on u > 11/5."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
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


STEM = "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate"
DEFAULT_OUT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
TAIL_GATE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_forward_remainder_tail_gate.json"
)
MODULAR_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
DATE = "2026-07-24"
PRECISION_BITS = 256
M_ORDER = 9
N_VALUES = tuple(range(4, 11))
U_EXACT = "11/5"
T_EXACT = "1/5"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def polynomial_rows() -> list[dict]:
    x = sp.symbols("X", nonnegative=True)
    current = 2 * x - 3
    rows: list[dict] = []
    for q in range(M_ORDER + 1):
        polynomial = sp.Poly(sp.expand(current), x)
        rows.append(
            {
                "q": q,
                "P_q": str(polynomial.as_expr()),
                "A_q": int(
                    sum(
                        abs(value)
                        for value in polynomial.all_coeffs()
                    )
                ),
            }
        )
        current = sp.expand(
            (5 - 4 * x) * polynomial.as_expr()
            + 4 * x * sp.diff(polynomial.as_expr(), x)
        )
    return rows


def laurent_data(expression: sp.Expr, y: sp.Symbol) -> tuple[sp.Expr, int]:
    coefficient_norm = sp.Rational(0)
    maximum_degree: int | None = None
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient, exponent = term.as_coeff_exponent(y)
        if coefficient.has(y) or not exponent.is_integer:
            raise RuntimeError(f"not a Laurent monomial: {term}")
        coefficient_norm += abs(coefficient)
        degree = int(exponent)
        maximum_degree = (
            degree
            if maximum_degree is None
            else max(maximum_degree, degree)
        )
    if maximum_degree is None:
        raise RuntimeError("empty Laurent polynomial")
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
            "bound": (
                "0<s(u)<=0.5*exp(-(9/16)*exp(8u))"
            ),
        }
    ]
    for order in range(1, M_ORDER + 1):
        coefficient_norm, maximum_degree = laurent_data(current, y)
        rows.append(
            {
                "a": order,
                "laurent_polynomial": str(current),
                "coefficient_norm": str(coefficient_norm),
                "maximum_degree": maximum_degree,
                "normalization": "1/sqrt(pi)",
                "bound": (
                    "|s^(a)(u)|<=C_a*exp(4*h_a*u)"
                    "*exp(-(9/16)*exp(8u))/sqrt(pi)"
                ),
            }
        )
        current = sp.expand(
            4 * y * sp.diff(current, y)
            - 2 * z * z_prime * current
        )
    return rows


def matrix_specs() -> list[dict]:
    specs: list[dict] = []
    for b in range(M_ORDER + 1):
        specs.append(
            {
                "component": "d0",
                "p": 0,
                "b": b,
                "q": M_ORDER - b,
                "coefficient": comb(M_ORDER, b),
            }
        )
        specs.append(
            {
                "component": "d1_u",
                "p": 1,
                "b": b,
                "q": M_ORDER - b,
                "coefficient": comb(M_ORDER, b),
            }
        )
    for b in range(M_ORDER):
        specs.append(
            {
                "component": "d1_lower",
                "p": 0,
                "b": b,
                "q": M_ORDER - 1 - b,
                "coefficient": M_ORDER * comb(M_ORDER - 1, b),
            }
        )
    return specs


def sympy_rational_to_arb(value: str) -> arb:
    rational = sp.Rational(value)
    return arb(int(rational.p)) / int(rational.q)


def heat_polynomial(order: int, t: arb, u: arb) -> arb:
    value = arb(0)
    for ell in range(order // 2 + 1):
        value += (
            factorial(order)
            * (2 * u) ** (order - 2 * ell)
            * t ** (order - ell)
            / (factorial(ell) * factorial(order - 2 * ell))
        )
    return value


def serialize_positive(value: arb) -> dict:
    if not value.is_finite() or value.lower() <= 0:
        raise RuntimeError(f"nonpositive or nonfinite bound: {value}")
    log10_value = value.log() / arb(10).log()
    return {
        "enclosure": value.str(65, more=True),
        "upper": value.upper().str(65),
        "log10_enclosure": log10_value.str(50, more=True),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def forward_outer_bound(
    retained: int,
    p: int,
    b: int,
    q: int,
    polynomials: list[dict],
) -> tuple[arb, dict]:
    u0 = arb(11) / 5
    t = arb(1) / 5
    y4 = (4 * u0).exp()
    first = retained + 1
    degree = 2 * q + 4
    linear = 4 * q + 9
    a_q = polynomials[q]["A_q"]
    pi = arb.pi()
    rho = (
        arb(degree) / first
        - pi * y4 * (2 * first + 1)
    ).exp()
    if rho.upper() >= 1:
        raise RuntimeError("forward arithmetic ratio is not below one")
    first_value = (
        a_q
        * pi ** (q + 2)
        * first**degree
        * (linear * u0 - pi * first**2 * y4).exp()
    )
    arithmetic_sum = first_value / (1 - rho)
    weight_at_split = (
        u0**p
        * heat_polynomial(b, t, u0)
        * (t * u0**2).exp()
    )
    decay_rate = (
        4 * pi * first**2 * y4
        - 2 * t * u0
        - linear
        - arb(p + b) / u0
    )
    if decay_rate.lower() <= 0:
        raise RuntimeError("forward outer decay rate is not positive")
    bound = weight_at_split * arithmetic_sum / decay_rate
    return bound, {
        "first_forward_index": first,
        "degree": degree,
        "linear_exponent": linear,
        "rho_enclosure": rho.str(45, more=True),
        "decay_rate_enclosure": decay_rate.str(45, more=True),
    }


def defect_outer_bound(
    retained: int,
    p: int,
    b: int,
    q: int,
    polynomials: list[dict],
    switch_rows: list[dict],
) -> tuple[arb, int]:
    u0 = arb(11) / 5
    t = arb(1) / 5
    y4 = (4 * u0).exp()
    y8 = (8 * u0).exp()
    pi = arb.pi()
    c = arb(9) / 16
    polynomial_weight = u0**p * heat_polynomial(b, t, u0)
    total = arb(0)
    term_count = 0

    for index in range(1, retained + 1):
        x_minus = pi * index**2 / y4
        if x_minus.upper() >= 1:
            raise RuntimeError(
                f"reflected X is not below one at n={index}"
            )
        for switch_order in range(q + 1):
            kernel_order = q - switch_order
            polynomial = polynomials[kernel_order]
            a_k = polynomial["A_q"]
            forward_linear = 4 * kernel_order + 9
            forward_at_split = (
                a_k
                * pi ** (kernel_order + 2)
                * index ** (2 * kernel_order + 4)
                * (
                    forward_linear * u0
                    - pi * index**2 * y4
                ).exp()
            )
            forward_decay_rate = (
                4 * pi * index**2 * y4 - forward_linear
            )
            if forward_decay_rate.lower() <= 5:
                raise RuntimeError(
                    "forward summand does not decay at rate five"
                )
            difference_coefficient = (
                a_k * pi * index**2
                + forward_at_split * (5 * u0).exp()
            )

            switch = switch_rows[switch_order]
            switch_degree = int(switch["maximum_degree"])
            switch_constant = sympy_rational_to_arb(
                switch["coefficient_norm"]
            )
            if switch_order > 0:
                switch_constant /= pi.sqrt()
            switch_linear = 4 * switch_degree - 5
            initial = (
                comb(q, switch_order)
                * switch_constant
                * difference_coefficient
                * polynomial_weight
                * (
                    t * u0**2
                    + switch_linear * u0
                    - c * y8
                ).exp()
            )
            decay_rate = (
                8 * c * y8
                - 2 * t * u0
                - switch_linear
                - arb(p + b) / u0
            )
            if decay_rate.lower() <= 0:
                raise RuntimeError(
                    "switch-defect decay rate is not positive"
                )
            contribution = initial / decay_rate
            total += contribution
            term_count += 1
    return total, term_count


def build_outer_rows(
    polynomials: list[dict], switch_rows: list[dict]
) -> tuple[list[dict], list[dict]]:
    entries: list[dict] = []
    aggregates: list[dict] = []
    for retained in N_VALUES:
        components = {
            "d0": arb(0),
            "d1_u": arb(0),
            "d1_lower": arb(0),
        }
        for spec in matrix_specs():
            forward, forward_diagnostic = forward_outer_bound(
                retained,
                spec["p"],
                spec["b"],
                spec["q"],
                polynomials,
            )
            defect, defect_term_count = defect_outer_bound(
                retained,
                spec["p"],
                spec["b"],
                spec["q"],
                polynomials,
                switch_rows,
            )
            total = forward + defect
            components[spec["component"]] += (
                spec["coefficient"] * total
            )
            entries.append(
                {
                    **spec,
                    "N": retained,
                    "forward_outer": serialize_positive(forward),
                    "retained_defect_outer": serialize_positive(defect),
                    "total_outer": serialize_positive(total),
                    "forward_diagnostic": forward_diagnostic,
                    "defect_term_count": defect_term_count,
                }
            )
        d1 = components["d1_u"] + components["d1_lower"]
        aggregates.append(
            {
                "N": retained,
                "outer_d0": serialize_positive(components["d0"]),
                "outer_d1_u": serialize_positive(
                    components["d1_u"]
                ),
                "outer_d1_lower": serialize_positive(
                    components["d1_lower"]
                ),
                "outer_d1": serialize_positive(d1),
            }
        )
    return entries, aggregates


def source_audit() -> dict:
    modular = json.loads(MODULAR_SOURCE.read_text(encoding="utf-8"))
    switch = modular["exact"]["modular_switch"]
    expected = (
        "(1-omega(u))/omega(u)<=erfc(3*sinh(4u))"
        "<=exp(-9*sinh(4u)^2), u>=0"
    )
    if switch.get("tail_bound") != expected:
        raise RuntimeError("upstream switch-tail bound drifted")
    tail = json.loads(TAIL_GATE.read_text(encoding="utf-8"))
    if tail.get("kind") != (
        "jensen_window_pf_newman_theta_forward_remainder_tail_gate"
    ):
        raise RuntimeError("forward-remainder tail gate kind drifted")
    return {
        "modular_source_sha256": file_hash(MODULAR_SOURCE),
        "tail_gate_sha256": file_hash(TAIL_GATE),
        "switch_tail_bound": expected,
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    polynomials = polynomial_rows()
    switch_rows = switch_derivative_rows()
    entries, aggregates = build_outer_rows(polynomials, switch_rows)
    rows = [
        GateRow(
            id="ntsrotg_01_stable_outer_split",
            role="exact_identity",
            readiness="proved",
            claim=(
                "The modular remainder is the full forward tail above N "
                "plus finitely many retained switch defects."
            ),
            formula=(
                "r_N=sum_(n>N)phi_n+sum_(n=1)^N "
                "s(u)*(phi_n(u)-phi_n(-u)), "
                "s=erfc(3*sinh(4u))/2"
            ),
            proof_boundary="Exact for 4<=N<=10 and all u>=0.",
        ),
        GateRow(
            id="ntsrotg_02_forward_outer_envelope",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Every weighted forward derivative tail has an explicit "
                "geometric arithmetic factor and positive outer decay rate."
            ),
            formula=(
                "F_q<=A_q*pi^(q+2)*n^(2q+4)"
                "*exp((4q+9)u-pi*n^2*exp(4u))"
            ),
            proof_boundary=(
                "Uses X>=1, log(1+1/n)<=1/n, and u>=11/5."
            ),
        ),
        GateRow(
            id="ntsrotg_03_switch_derivative_envelope",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Each positive-order switch derivative is a Laurent "
                "polynomial times exp(-9*sinh(4u)^2)/sqrt(pi)."
            ),
            formula=(
                "s^(a)=L_a(exp(4u))*exp(-9*sinh(4u)^2)/sqrt(pi), "
                "|L_a(y)|<=C_a*y^h_a"
            ),
            proof_boundary=(
                "The a=0 row separately uses erfc(z)<=exp(-z^2)."
            ),
        ),
        GateRow(
            id="ntsrotg_04_reflected_kernel_envelope",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "For retained n<=10, reflected kernel derivatives have a "
                "uniform exp(-5u) envelope on the outer interval."
            ),
            formula=(
                "|D^k(phi_n(u)-phi_n(-u))|<=G_(k,n)*exp(-5u)"
            ),
            proof_boundary=(
                "Uses pi*n^2*exp(-4u)<1 and explicit forward decay."
            ),
        ),
        GateRow(
            id="ntsrotg_05_heat_absorption",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "The heat exponential and polynomial weight are absorbed "
                "without numerical quadrature."
            ),
            formula=(
                "u^p*P_b(T,u)<=U^p*P_b(T,U)"
                "*exp((p+b)(u-U)/U); "
                "u^2-U^2<=U*(exp(ku)-exp(kU))/(k*exp(kU)/2), "
                "k in {4,8}"
            ),
            proof_boundary=(
                "Applied only at U=11/5, T=1/5, p<=1, b<=9."
            ),
        ),
        GateRow(
            id="ntsrotg_06_outer_matrix",
            role="directed_rounding_certificate",
            readiness="proved",
            claim=(
                "The 29 weighted outer integrals are explicitly bounded "
                "for every retained N from 4 through 10."
            ),
            formula=(
                "outer_entry<=forward_initial/lambda_4"
                "+sum_(n<=N,a<=q)defect_initial/lambda_8"
            ),
            proof_boundary=(
                "Arb encloses constants only; no quadrature or sampled "
                "tail cutoff is used."
            ),
        ),
        GateRow(
            id="ntsrotg_07_compact_handoff",
            role="proof_search_target",
            readiness="not_ready_to_apply",
            claim=(
                "Add these outer d0/d1 bounds to the independently checked "
                "stable compact matrix."
            ),
            formula="d_j<=compact_d_j+outer_d_j",
            proof_boundary=(
                "Retained J/J' lower separation and transition-cell "
                "coverage remain open."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": "explicit full-arithmetic outer-u derivative theorem",
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "m_order": M_ORDER,
            "n_values": list(N_VALUES),
            "u_split_exact": U_EXACT,
            "t_cap_exact": T_EXACT,
            "forward_phase": "pi*n^2*exp(4u)",
            "switch_phase_lower": "(9/16)*exp(8u)",
        },
        "source_audit": source_audit(),
        "derivative_polynomials": polynomials,
        "switch_derivatives": switch_rows,
        "outer_entries": entries,
        "outer_aggregates": aggregates,
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact bounds the complete arithmetic remainder on "
            "u>11/5 for m=9, j=0,1, T=1/5 and N=4..10. It does not "
            "supply the compact matrix, retained first-jet J/J' lower balls, "
            "transition-cell coverage, Lambda<=0, RH, or a Clay-prize "
            "conclusion."
        ),
        "versions": {
            "python_flint": getattr(flint, "__version__", "unknown"),
            "sympy": sp.__version__,
        },
    }


def render_note(artifact: dict) -> str:
    lines = [
        "# Newman Theta Stable-Remainder Outer-Tail Gate",
        "",
        f"Date: {DATE}",
        "",
        "Status: explicit rigorous full-arithmetic outer-u theorem.",
        "This is not a retained first-jet separation theorem and not a proof",
        "of `Lambda<=0`, RH, or a Clay-prize result.",
        "",
        "## Exact Split",
        "",
        "```text",
        "s(u)=erfc(3*sinh(4u))/2",
        "delta_n=s(u)*(phi_n(u)-phi_n(-u))",
        "r_N=sum_(n>N)phi_n+sum_(n=1)^N delta_n",
        "```",
        "",
        "For `u>=11/5`, the forward phase is",
        "`pi*n^2*exp(4u)` and the switch phase is at least",
        "`(9/16)*exp(8u)`. Positive heat-polynomial growth is absorbed",
        "into those phases, leaving a positive explicit decay rate.",
        "",
        "## Certified Outer Budgets",
        "",
        "| N | outer d0 upper | outer d1 upper |",
        "|---:|---:|---:|",
    ]
    for row in artifact["outer_aggregates"]:
        lines.append(
            f"| {row['N']} | `{row['outer_d0']['upper']}` | "
            f"`{row['outer_d1']['upper']}` |"
        )
    lines.extend(
        [
            "",
            "All constants are Arb enclosures at 256-bit precision.",
            "No improper numerical quadrature or sampled tail cutoff is used.",
            "",
            "## Remaining Handoff",
            "",
            "Add these bounds to the stable compact matrix. The retained",
            "`J/J'` lower separation and transition-cell theorem remain open.",
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
        "wrote Newman theta stable-remainder outer-tail gate: "
        f"{len(artifact['outer_entries'])} entries, "
        f"{len(artifact['outer_aggregates'])} retained counts"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
