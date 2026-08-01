#!/usr/bin/env python3
"""Build the Edrei Hankel boundary-flux and escaping-order gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_edrei_hankel_boundary_flux_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_edrei_hankel_boundary_flux_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def poly_mul(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] += x * y
    return out


def annihilator(betas: list[Fraction]) -> list[Fraction]:
    out = [Fraction(1)]
    for beta in betas:
        out = poly_mul(out, [-beta, Fraction(1)])
    return out


def derivative_value(coefficients: list[Fraction], x: Fraction) -> Fraction:
    return sum(
        i * coefficients[i] * x ** (i - 1)
        for i in range(1, len(coefficients))
    )


def moment_values(
    atoms: list[tuple[Fraction, Fraction]], maximum: int
) -> list[Fraction]:
    return [
        sum(rho * beta**n for beta, rho in atoms)
        for n in range(maximum + 1)
    ]


def flow_rhs(values: list[Fraction], n: int) -> Fraction:
    convolution = sum(values[k] * values[n - k] for k in range(n + 1))
    return (
        -2 * (n + 1) * (2 * n + 3) * values[n + 1]
        + 4 * (n + 1) * convolution
    )


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    sign = 1
    out = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        out *= pivot_value
        for row in range(column + 1, len(work)):
            ratio = work[row][column] / pivot_value
            for entry in range(column + 1, len(work)):
                work[row][entry] -= ratio * work[column][entry]
    return sign * out


def boundary_audit() -> dict:
    cases = (
        ([Fraction(3, 2)], [1]),
        ([Fraction(3, 2)], [2]),
        ([Fraction(3, 2), Fraction(2, 3)], [1, 3]),
        (
            [Fraction(5, 4), Fraction(3, 4), Fraction(1, 5)],
            [1, 2, 4],
        ),
    )
    quadratic_checks = 0
    determinant_checks = 0
    samples: list[dict] = []
    for betas, multiplicities in cases:
        atoms = [
            (beta, Fraction(multiplicity) * beta)
            for beta, multiplicity in zip(betas, multiplicities)
        ]
        coefficients = annihilator(betas)
        degree = len(coefficients) - 1
        for shift in range(4):
            values = moment_values(atoms, shift + 2 * degree + 1)
            derivatives = [
                flow_rhs(values, n) for n in range(shift + 2 * degree + 1)
            ]
            observed = sum(
                coefficients[i]
                * coefficients[j]
                * derivatives[shift + i + j]
                for i in range(degree + 1)
                for j in range(degree + 1)
            )
            expected = 8 * sum(
                rho
                * beta ** (shift + 2)
                * (rho - beta)
                * derivative_value(coefficients, beta) ** 2
                for beta, rho in atoms
            )
            if observed != expected:
                raise RuntimeError("finite-rank boundary flux audit failed")
            quadratic_checks += 1
        if degree == 0:
            continue
        for shift in range(1, 4):
            values = moment_values(atoms, shift + 2 * degree + 1)
            derivatives = [
                flow_rhs(values, n) for n in range(shift + 2 * degree + 1)
            ]
            size = degree + 1
            hankel = [
                [values[shift + i + j] for j in range(size)]
                for i in range(size)
            ]
            hankel_dot = [
                [derivatives[shift + i + j] for j in range(size)]
                for i in range(size)
            ]
            observed = Fraction(0)
            for row in range(size):
                replaced = [entries[:] for entries in hankel]
                replaced[row] = hankel_dot[row][:]
                observed += determinant(replaced)
            tau = Fraction(1)
            for beta, rho in atoms:
                tau *= rho * beta**shift
            for i, beta_i in enumerate(betas):
                for beta_j in betas[i + 1 :]:
                    tau *= (beta_i - beta_j) ** 2
            boundary_form = 8 * sum(
                rho
                * beta ** (shift + 2)
                * (rho - beta)
                * derivative_value(coefficients, beta) ** 2
                for beta, rho in atoms
            )
            expected = tau * boundary_form
            if observed != expected:
                raise RuntimeError("finite-rank determinant flux audit failed")
            determinant_checks += 1
            if len(samples) < 3:
                samples.append(
                    {
                        "betas": [str(beta) for beta in betas],
                        "multiplicities": multiplicities,
                        "shift": shift,
                        "determinant_derivative": str(observed),
                    }
                )
    fractional_beta = Fraction(1)
    fractional_rho = Fraction(1, 2)
    fractional_flux = (
        8
        * fractional_rho
        * fractional_beta**3
        * (fractional_rho - fractional_beta)
    )
    if fractional_flux != -2:
        raise RuntimeError("fractional-residue guard failed")
    return {
        "quadratic_form_checks": quadratic_checks,
        "determinant_checks": determinant_checks,
        "fractional_residue_flux": str(fractional_flux),
        "samples": samples,
    }


def symbolic_two_atom_audit() -> dict:
    shift = sp.symbols("s", integer=True, nonnegative=True)
    q = sp.symbols("q", positive=True)
    multiplicity = sp.symbols("m", integer=True, positive=True)
    tail_multiplicity = sp.symbols("n", integer=True, positive=True)

    def moment(index: sp.Expr) -> sp.Expr:
        return multiplicity + tail_multiplicity * q ** (index + 1)

    def convolution(index: sp.Expr) -> sp.Expr:
        return (
            (index + 1) * multiplicity**2
            + 2
            * multiplicity
            * tail_multiplicity
            * q
            * (1 - q ** (index + 1))
            / (1 - q)
            + (index + 1)
            * tail_multiplicity**2
            * q ** (index + 2)
        )

    def flow(index: sp.Expr) -> sp.Expr:
        return (
            -2 * (index + 1) * (2 * index + 3) * moment(index + 1)
            + 4 * (index + 1) * convolution(index)
        )

    hankel_det = sp.factor(
        moment(shift) * moment(shift + 2) - moment(shift + 1) ** 2
    )
    hankel_det_dot = sp.factor(
        flow(shift) * moment(shift + 2)
        + moment(shift) * flow(shift + 2)
        - 2 * moment(shift + 1) * flow(shift + 1)
    )
    simple_ratio = sp.factor(
        (hankel_det_dot / hankel_det).subs(
            {multiplicity: 1, tail_multiplicity: 1}
        )
    )
    repeated_half_ratio = sp.factor(
        (hankel_det_dot / hankel_det).subs(
            {
                multiplicity: 2,
                tail_multiplicity: 1,
                q: sp.Rational(1, 2),
            }
        )
    )
    expected_simple = (
        -2
        * (q + 1)
        * (
            q**2 * shift
            + 3 * q**2
            - 2 * q * shift
            - 14 * q
            + shift
            + 3
        )
        / (q - 1) ** 2
    )
    expected_repeated = 128 * 2**shift + 4 * shift**2 + 29 * shift + 131
    if sp.simplify(simple_ratio - expected_simple) != 0:
        raise RuntimeError("simple two-atom ratio identity failed")
    if sp.simplify(repeated_half_ratio - expected_repeated) != 0:
        raise RuntimeError("repeated two-atom ratio identity failed")
    return {
        "measure": "nu=m*delta_1+n*q*delta_q, 0<q<1",
        "hankel_determinant": str(hankel_det),
        "simple_log_derivative": str(simple_ratio),
        "repeated_q_half_log_derivative": str(repeated_half_ratio),
        "simple_root_limit": "1",
        "repeated_root_limit": "2",
    }


def exact_statements() -> dict:
    return {
        "moment_flow": "a_n'=-2(n+1)(2n+3)a_(n+1)+4(n+1)sum_(k=0)^n a_k*a_(n-k)",
        "quadratic_form": "Q_s(P)=sum_(i,j)c_i*c_j*a_(s+i+j)=integral x^s*P(x)^2 dnu(x)",
        "operator_identity": "Q_s'(P)=-2 integral x(E+1)(2E+3)[x^s P^2]dnu+4 double_integral [h_s(x)-h_s(y)]/(x-y)dnu(x)dnu(y), h_s(t)=t*d/dt[t^(s+1)P(t)^2]",
        "edrei_measure": "For F in LP+: dnu=gamma*delta_0+sum_j rho_j*delta_(beta_j), rho_j=m_j*beta_j, m_j positive integers",
        "null_flux": "If s>=1 and P(beta_j)=0 on every positive atom, then Q_s'(P)=8 sum_j rho_j*beta_j^(s+2)*(rho_j-beta_j)*P'(beta_j)^2",
        "integer_flux": "For Edrei residues rho_j=m_j*beta_j: Q_s'(P)=8 sum_j m_j*(m_j-1)*beta_j^(s+4)*P'(beta_j)^2>=0",
        "determinant_flux": "For r positive atoms, monic P=product_j(t-beta_j), and D_(r,s)=det(a_(s+i+j))_(i,j=0)^r: D_(r,s)'=tau_(r,s)*Q_s'(P), tau_(r,s)=product_j(rho_j*beta_j^s)*product_(i<j)(beta_i-beta_j)^2",
        "fractional_guard": "For nu=(1/2)delta_1, P(t)=t-1, and s=1, Q_s'(P)=-2: Stieltjes positivity without the integer Edrei residue condition does not orient the wall forward",
        "infinite_support": "An LP+ Edrei measure with infinitely many distinct positive atoms makes every finite H_(d,s)=(a_(s+i+j)) positive definite for s=0 and s=1",
        "canonical_escape": "If a continuous family exits LP+ through an infinite-support boundary, every outside sequence approaching the boundary requires canonical Stieltjes witness order d->infinity in one of the shifts s=0 or s=1",
        "orthogonal_quotient": "For the monic x^s*dnu-orthogonal polynomial pi_(r,s): h_(r,s)=Q_s(pi_(r,s))=D_(r,s)/D_(r-1,s) and partial_lambda log D_(r,s)=partial_lambda log D_(r-1,s)+Q_s'(pi_(r,s))/Q_s(pi_(r,s))",
        "orthogonal_convergence": "If P_r(t)=product_(j=1)^r(t-beta_j), then pi_(r,s)(beta_i)=O((beta_(r+1)/beta_i)^s) for i<=r and pi_(r,s)=P_r+O((beta_(r+1)/beta_r)^s) coefficientwise",
        "collision_sensor": "If beta_1>beta_2>... and r is the first index with m_r>=2, then lim_(s->infinity)(partial_lambda log D_(r,s))^(1/s)=beta_r/beta_(r+1)>1",
        "newman_consequence": "Lambda<=0 iff, for every fixed r and every lambda in (0,1/5] for which H_lambda has only real zeros, limsup_(s->infinity)max(1,partial_lambda log D_(r,s))^(1/s)<=1",
        "open_target": "Prove directly from the Phi transform that for every fixed r and every lambda in (0,1/5] for which H_lambda has only real zeros, limsup_(s->infinity)max(1,partial_lambda log D_(r,s))^(1/s)<=1",
    }


def build_payload() -> dict:
    exact = exact_statements()
    rows = [
        GateRow(
            "ehbf_01_moment_flow",
            "exact_input",
            "available_exact",
            "The Edrei moments obey the radial heat hierarchy.",
            exact["moment_flow"],
            "Inherited from Lemma 11.22C; no positivity follows yet.",
        ),
        GateRow(
            "ehbf_02_quadratic_form",
            "exact_identity",
            "available_exact",
            "Every shifted Hankel quadratic form is a polynomial square norm when an LP+ representing measure exists.",
            exact["quadratic_form"],
            "The representation is available only on the LP+ side or boundary.",
        ),
        GateRow(
            "ehbf_03_operator_flux",
            "exact_identity",
            "available_exact",
            "The nonlinear moment hierarchy has an exact polynomial divided-difference flux.",
            exact["operator_identity"],
            "The divided difference is interpreted diagonally as h_s'(x).",
        ),
        GateRow(
            "ehbf_04_integer_edrei_measure",
            "exact_structure",
            "available_exact",
            "Entireness quantizes every positive pole residue of the Stieltjes logarithmic derivative.",
            exact["edrei_measure"],
            "This is the Sokal/Edrei LP+ representation, not a property of arbitrary Stieltjes measures.",
        ),
        GateRow(
            "ehbf_05_null_form_flux",
            "exact_boundary_identity",
            "available_exact",
            "At finite support the null-form flux separates atom by atom.",
            exact["null_flux"],
            "The polynomial must vanish on the complete positive support.",
        ),
        GateRow(
            "ehbf_06_integer_residue_orientation",
            "exact_boundary_theorem",
            "available_exact",
            "The finite-rank Edrei wall is nonnegative forward and is strict exactly at repeated factors.",
            exact["integer_flux"],
            "Simple factors are tangent at first order.",
        ),
        GateRow(
            "ehbf_07_vandermonde_determinant_flux",
            "exact_boundary_theorem",
            "available_exact",
            "The first vanishing Hankel determinant has the same multiplicity flux times an explicit positive Vandermonde factor.",
            exact["determinant_flux"],
            "Stated for s>=1 so the optional exponential atom at zero is invisible.",
        ),
        GateRow(
            "ehbf_08_rank_one_recovery",
            "exact_consistency_check",
            "available_exact",
            "The r=1 specialization recovers the previous repeated-zero 2x2 formula.",
            "D_(1,s)'=8*m^2*(m-1)*beta^(2s+5)",
            "This is a consistency check, not a new closing estimate.",
        ),
        GateRow(
            "ehbf_09_fractional_residue_guard",
            "exact_countermodel",
            "guard_validated",
            "A positive Stieltjes atom with non-Edrei residue can point out of the cone even under forward heat.",
            exact["fractional_guard"],
            "Its logarithmic primitive is not entire; the example isolates why residue quantization matters.",
        ),
        GateRow(
            "ehbf_10_infinite_support_strictness",
            "exact_structure",
            "available_exact",
            "Infinite LP+ support makes every fixed finite shifted Hankel block strictly positive.",
            exact["infinite_support"],
            "Strict finite blocks do not supply a uniform all-order margin.",
        ),
        GateRow(
            "ehbf_11_canonical_order_escape",
            "exact_nonpromotion_gate",
            "available_exact",
            "Approach to an infinite-support Newman boundary cannot be detected at one fixed canonical Hankel order.",
            exact["canonical_escape"],
            "Finite order-eleven and finite-shift certificates remain bounded evidence.",
        ),
        GateRow(
            "ehbf_12_cauchy_binet_expansion",
            "exact_identity",
            "available_exact",
            "For s>=1, Cauchy-Binet expands every shifted determinant as a positive sum of weighted squared Vandermondes.",
            "For s>=1: D_(d,s)=sum_(|J|=d+1) product_(j in J)(rho_j*beta_j^s)*V(beta_J)^2",
            "At s=0 the optional gamma*delta_0 atom adds nonnegative Cauchy-Binet terms; the displayed positive-atom sum still proves strictness.",
        ),
        GateRow(
            "ehbf_13_high_shift_collision_sensor",
            "exact_boundary_asymptotic",
            "available_exact",
            "The first repeated atom forces exponential growth of one fixed-size shifted determinant logarithmic derivative.",
            exact["orthogonal_quotient"] + "; " + exact["collision_sensor"],
            "The Schur quotient and orthogonality remove the coefficient-derivative term; the ratio beta_r/beta_(r+1) is greater than one because the distinct atoms are strictly ordered.",
        ),
        GateRow(
            "ehbf_14_simple_two_atom_guard",
            "exact_countermodel",
            "guard_validated",
            "Two simple atoms have only polynomial shifted determinant logarithmic growth.",
            "For nu=delta_1+q*delta_q: D_(1,s)'/D_(1,s)=-2(q+1)[q^2*s+3q^2-2q*s-14q+s+3]/(q-1)^2",
            "Shows that the exponential signature is multiplicity-sensitive.",
        ),
        GateRow(
            "ehbf_15_repeated_two_atom_witness",
            "exact_boundary_witness",
            "guard_validated",
            "A repeated leading atom gives the predicted exponential signature exactly.",
            "For nu=2delta_1+(1/2)delta_(1/2): D_(1,s)'/D_(1,s)=128*2^s+4s^2+29s+131",
            "A finite LP+ witness, not the Xi function.",
        ),
        GateRow(
            "ehbf_16_positive_newman_consequence",
            "exact_equivalence",
            "available_exact",
            "The fixed-size subexponential high-shift determinant-flux criterion is exactly equivalent to Lambda<=0.",
            exact["newman_consequence"],
            "The forward direction uses positive-time simplicity; the reverse direction uses positive-boundary attainment and the collision sensor.",
        ),
        GateRow(
            "ehbf_17_subexponential_phi_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "A subexponential fixed-size determinant-flux theorem would exclude every positive Newman collision.",
            exact["open_target"],
            "No such Xi/Phi-specific estimate is proved here.",
        ),
        GateRow(
            "ehbf_18_proof_boundary",
            "proof_boundary",
            "guard_validated",
            "The gate separates exact flux, exact asymptotics, countermodels, and the missing Phi estimate.",
            "finite-rank flux and collision sensor proved; Xi/Phi subexponential bound open",
            "Nothing in this artifact proves PF-infinity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_edrei_hankel_boundary_flux_gate",
        "date": "2026-07-24",
        "status": (
            "exact finite-rank Hankel boundary flux, infinite-support order escape, "
            "and high-shift collision sensor with one open Xi/Phi determinant-growth gate"
        ),
        "exact": exact,
        "boundary_audit": boundary_audit(),
        "two_atom_audit": symbolic_two_atom_audit(),
        "rows": [asdict(row) for row in rows],
        "sources": [
            {
                "title": "Sokal, When does a hypergeometric function belong to LP+?, Proposition 6",
                "url": "https://arxiv.org/abs/2204.01045",
            },
            {
                "title": "D. H. J. Polymath, Effective approximation of heat flow evolution of the Riemann xi function",
                "url": "https://arxiv.org/abs/1904.12438",
            },
            {
                "title": "Local Edrei heat-flow boundary gate",
                "path": "outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md",
            },
            {
                "title": "Local positive Newman boundary and strict Laguerre correlation results",
                "path": "outputs/formal_core.md",
            },
        ],
        "proof_boundary": (
            "This artifact proves the polynomial Hankel-flux identity, the finite-rank "
            "integer-residue orientation and determinant formula, strictness and "
            "canonical-order escape at an infinite-support LP+ boundary, and the "
            "high-shift determinant signature forced by a repeated atom. It does not "
            "prove the required Xi/Phi subexponential determinant-growth estimate, "
            "all-order Stieltjes positivity at lambda zero, PF-infinity, Jensen "
            "hyperbolicity for zeta, RH, or Lambda<=0."
        ),
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    audit = payload["boundary_audit"]
    two_atom = payload["two_atom_audit"]
    return "\n".join(
        [
            "# Jensen-Window PF Edrei Hankel Boundary-Flux Gate",
            "",
            "Date: 2026-07-24",
            "",
            "Status: exact finite-rank Hankel boundary flux, infinite-support",
            "canonical-order escape, and high-shift collision sensor with one open",
            "Xi/Phi determinant-growth gate. This is not a proof of PF-infinity,",
            "Jensen hyperbolicity for zeta, RH, or `Lambda <= 0`.",
            "",
            "Artifact kind: `jensen_window_pf_edrei_hankel_boundary_flux_gate`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_edrei_hankel_boundary_flux_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_edrei_hankel_boundary_flux_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_edrei_hankel_boundary_flux_gate.py",
            "```",
            "",
            "## Polynomial Flux Identity",
            "",
            "Retain the exact moment hierarchy from Lemma 11.22C:",
            "",
            "```text",
            exact["moment_flow"],
            "```",
            "",
            "At a time where the logarithmic derivative has a Stieltjes measure",
            "`nu`, write `P(t)=sum_i c_i t^i`, put `E=t*d/dt`, and define",
            "",
            "```text",
            exact["quadratic_form"],
            "f_s(t)=t^s*P(t)^2",
            "g_s(t)=t*f_s(t)",
            "h_s(t)=t*g_s'(t)",
            "```",
            "",
            "Coefficient extraction and",
            "",
            "```text",
            "(n+1)sum_(k=0)^n x^k*y^(n-k)",
            " =[x*d/dx(x^(n+1))-y*d/dy(y^(n+1))]/(x-y)",
            "```",
            "",
            "give the exact divided-difference identity",
            "",
            "```text",
            exact["operator_identity"],
            "```",
            "",
            "with the diagonal value interpreted as `h_s'(x)`. This is an identity",
            "for every polynomial `P`; no sign has been imposed.",
            "",
            "## Finite-Rank Edrei Boundary",
            "",
            "Sokal's logarithmic-derivative criterion supplies the quantized measure",
            "",
            "```text",
            exact["edrei_measure"],
            "```",
            "",
            "Suppose `s>=1` and `P` vanishes on every positive atom. Off-diagonal",
            "divided differences vanish. At one atom `beta`, direct differentiation",
            "gives",
            "",
            "```text",
            "-2*beta*(E+1)*(2E+3)[t^s P(t)^2]_(t=beta)",
            " =-8*beta^(s+3)*P'(beta)^2",
            "h_s'(beta)=2*beta^(s+2)*P'(beta)^2.",
            "```",
            "",
            "The linear and diagonal quadratic pieces therefore combine to",
            "",
            "```text",
            exact["null_flux"],
            exact["integer_flux"],
            "```",
            "",
            "This proves forward nonnegativity at every finite-rank Edrei wall.",
            "It is strict precisely when a factor is repeated; simple factors are",
            "tangent at first order.",
            "",
            "For `r` positive atoms and the monic annihilator",
            "`P(t)=product_j(t-beta_j)`, the shifted moment matrix has rank `r`.",
            "Its adjugate is the outer product of the coefficient vector of `P`",
            "times its leading Vandermonde cofactor. Hence",
            "",
            "```text",
            exact["determinant_flux"],
            "```",
            "",
            "At `r=1` this recovers",
            "",
            "```text",
            "D_(1,s)'=8*m^2*(m-1)*beta^(2s+5).",
            "```",
            "",
            f"The builder checks {audit['quadratic_form_checks']} exact rational null-form",
            f"instances and {audit['determinant_checks']} exact rational determinant",
            "instances.",
            "",
            "## Residue-Quantization Guard",
            "",
            "The sign is not a generic consequence of Stieltjes positivity:",
            "",
            "```text",
            exact["fractional_guard"],
            "```",
            "",
            "The formal logarithmic primitive has a fractional power and is not",
            "entire. The integer Edrei residue, not positivity of `nu` alone, is the",
            "structural input behind the forward orientation.",
            "",
            "## Infinite-Support Boundary And Escaping Order",
            "",
            "For an LP+ endpoint with distinct atoms",
            "`beta_1>beta_2>...>0` and `s>=1`, Cauchy-Binet gives",
            "",
            "```text",
            "D_(d,s)=sum_(|J|=d+1) product_(j in J)(rho_j*beta_j^s)",
            "          *V(beta_J)^2>0.",
            "```",
            "",
            "For `s=0`, the same positive-atom sum is strictly positive and the",
            "optional `gamma*delta_0` atom contributes only additional nonnegative",
            "Cauchy-Binet terms. Thus every fixed block in both canonical shifts",
            "`s=0,1` is positive definite.",
            "",
            "Thus every fixed finite shifted block is positive definite when the",
            "support is infinite. If a continuous family approaches such a boundary",
            "from outside LP+, each fixed canonical block at shifts zero and one",
            "stays positive in a neighborhood of the boundary. Sokal's criterion",
            "then proves",
            "",
            "```text",
            exact["canonical_escape"],
            "```",
            "",
            "This is stronger than saying that a large finite computation has not",
            "yet found the right row: no fixed canonical order can be the boundary",
            "detector.",
            "",
            "## High-Shift Collision Sensor",
            "",
            "Assume the first repeated atom is `beta_r`, so",
            "`m_1=...=m_(r-1)=1` and `m_r>=2`. Let `pi_(r,s)` be the monic",
            "degree-`r` orthogonal polynomial for `x^s*dnu`, and put",
            "`h_(r,s)=Q_s(pi_(r,s))`. Schur complementation and differentiation",
            "at fixed leading coefficient give the exact identities",
            "",
            "```text",
            exact["orthogonal_quotient"],
            "```",
            "",
            "Indeed, `partial_lambda pi_(r,s)` has degree at most `r-1`, so its",
            "cross term with `pi_(r,s)` vanishes by orthogonality. Cauchy-Binet",
            "then gives",
            "",
            "```text",
            "D_(r,s)=C_D*(beta_1*...*beta_r*beta_(r+1))^s*(1+o(1)).",
            "D_(r-1,s)=C_0*(beta_1*...*beta_r)^s*(1+o(1)).",
            "h_(r,s)=C_H*beta_(r+1)^s*(1+o(1)),",
            "```",
            "",
            "To control the numerator, use the Lagrange polynomial `L_i` on the",
            "first `r` atoms in the orthogonality equation. Since `a_0<infinity`,",
            "the tail is uniformly summable and",
            "",
            "```text",
            exact["orthogonal_convergence"],
            "```",
            "",
            "Inserting this into the exact divided-difference flux shows that every",
            "simple leading atom and every tail term is",
            "`O(s^2*beta_(r+1)^s)`, while the repeated atom contributes",
            "`8*m_r*(m_r-1)*beta_r^(s+4)*P_r'(beta_r)^2*(1+o(1))`.",
            "The lower logarithmic derivative is only `O(s^2)`, because its",
            "leading `r` distinct atoms already give full rank. Consequently",
            "",
            "```text",
            "partial_lambda D_(r,s)",
            " =C_N*(beta_1*...*beta_(r-1)*beta_r^2)^s*(1+o(1)),",
            exact["collision_sensor"],
            "```",
            "",
            "where `C_D,C_N>0`. In the exact two-atom benchmark",
            f"`{two_atom['measure']}`, the simple leading atom gives",
            "",
            "```text",
            two_atom["simple_log_derivative"],
            "```",
            "",
            "whose `s`th-root growth is one. Replacing the leading atom by a double",
            "factor and taking `q=1/2` gives",
            "",
            "```text",
            two_atom["repeated_q_half_log_derivative"],
            "```",
            "",
            "whose `s`th-root growth is two, exactly `beta_1/beta_2`.",
            "",
            "## Newman Reduction",
            "",
            "If `Lambda<=0`, Lemma 11.9 says every zero is simple for every",
            "`0<lambda<=1/5`. For each fixed `r`, the same orthogonal-quotient",
            "argument now has no repeated atom among the annihilated leading support:",
            "both `Q_s'(pi_(r,s))` and `Q_s(pi_(r,s))` have the",
            "`beta_(r+1)^s` scale up to polynomial factors. Recursing through the",
            "lower quotients gives `partial_lambda log D_(r,s)=O(s^2)`, so the",
            "subexponential condition below holds.",
            "",
            "Conversely, the audited positive-boundary theorem says that `Lambda>0`",
            "produces a finite multiple real zero of `H_Lambda` with",
            "`0<Lambda<=1/5`. The transform has infinitely many distinct zero",
            "locations, so its Edrei measure has a next atom. The first repeated",
            "descending reciprocal-square atom then violates the condition by the",
            "collision sensor. Therefore the following criterion is exactly",
            "equivalent to `Lambda<=0`:",
            "",
            "```text",
            exact["open_target"],
            "```",
            "",
            "This equivalence does not prove its open direction. A proof must use the",
            "actual Phi family. The repeated quadratic from",
            "Lemma 11.22C violates this bound, so it cannot follow from LP+,",
            "Stieltjes positivity, radial heat, or positive scale mixing alone.",
            "",
            "## Literature Audit",
            "",
            "Sokal's Proposition 6 supplies the discrete integer-residue measure.",
            "Standard Gram/Cauchy-Binet identities supply the positive Vandermonde",
            "expansion. The searched Hankel/Toda and orthogonal-polynomial machinery",
            "organizes moment determinants but does not provide the required",
            "Phi-specific subexponential heat-flux estimate. The Polymath positive-time",
            "asymptotics supply infinite high zeros and boundary localization, not",
            "the missing determinant bound.",
            "",
            "## Proof Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built Edrei Hankel boundary-flux gate: "
        f"{len(payload['rows'])} rows, "
        f"{payload['boundary_audit']['quadratic_form_checks']} null-form audits, "
        f"{payload['boundary_audit']['determinant_checks']} determinant audits, "
        "1 open Xi/Phi determinant-growth gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
