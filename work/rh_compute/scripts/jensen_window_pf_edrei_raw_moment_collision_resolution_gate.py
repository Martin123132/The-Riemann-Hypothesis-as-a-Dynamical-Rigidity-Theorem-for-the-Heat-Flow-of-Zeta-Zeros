#!/usr/bin/env python3
"""Build the Edrei raw-moment collision-resolution gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_edrei_raw_moment_collision_resolution_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def polynomial_from_atoms(atoms: list[Fraction]) -> list[Fraction]:
    coefficients = [Fraction(1)]
    for atom in atoms:
        extended = [Fraction(0)] * (len(coefficients) + 1)
        for index, value in enumerate(coefficients):
            extended[index] += value
            extended[index + 1] += atom * value
        coefficients = extended
    return coefficients


def log_moments_from_taylor(
    coefficients: list[Fraction], maximum: int
) -> list[Fraction]:
    padded = coefficients + [Fraction(0)] * (maximum + 2)
    moments: list[Fraction] = []
    for n in range(maximum + 1):
        lower = sum(
            (-1) ** k * moments[k] * padded[n - k] for k in range(n)
        )
        signed = (n + 1) * padded[n + 1] - lower
        moments.append((-1) ** n * signed)
    return moments


def edrei_moments(
    atoms: list[tuple[Fraction, Fraction]], maximum: int
) -> list[Fraction]:
    return [
        sum(residue * beta**n for beta, residue in atoms)
        for n in range(maximum + 1)
    ]


def flow_rhs(values: list[Fraction], n: int) -> Fraction:
    return (
        -2 * (n + 1) * (2 * n + 3) * values[n + 1]
        + 4
        * (n + 1)
        * sum(values[k] * values[n - k] for k in range(n + 1))
    )


def shifted_two_by_two(values: list[Fraction], shift: int) -> Fraction:
    return (
        values[shift] * values[shift + 2] - values[shift + 1] ** 2
    )


def shifted_two_by_two_derivative(
    values: list[Fraction], derivatives: list[Fraction], shift: int
) -> Fraction:
    return (
        derivatives[shift] * values[shift + 2]
        + values[shift] * derivatives[shift + 2]
        - 2 * values[shift + 1] * derivatives[shift + 1]
    )


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    sign = 1
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (
                row
                for row in range(column, len(work))
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, len(work)):
            ratio = work[row][column] / pivot_value
            for entry in range(column + 1, len(work)):
                work[row][entry] -= ratio * work[column][entry]
    return sign * result


def shifted_determinant(
    values: list[Fraction], rank: int, shift: int
) -> Fraction:
    return determinant(
        [
            [values[shift + row + column] for column in range(rank + 1)]
            for row in range(rank + 1)
        ]
    )


def vandermonde_square(atoms: list[Fraction]) -> Fraction:
    result = Fraction(1)
    for i, left in enumerate(atoms):
        for right in atoms[i + 1 :]:
            result *= (left - right) ** 2
    return result


def cauchy_binet_sum(
    weighted_atoms: list[tuple[Fraction, Fraction]],
    rank: int,
    shift: int,
) -> Fraction:
    total = Fraction(0)
    for indices in combinations(range(len(weighted_atoms)), rank + 1):
        selected = [weighted_atoms[index] for index in indices]
        atoms = [atom for atom, _ in selected]
        weight = Fraction(1)
        for atom, residue in selected:
            weight *= residue * atom**shift
        total += weight * vandermonde_square(atoms)
    return total


def simple_atom_velocities(atoms: list[Fraction]) -> list[Fraction]:
    velocities = []
    for i, beta_i in enumerate(atoms):
        interaction = sum(
            atoms[j] / (beta_i - atoms[j])
            for j in range(len(atoms))
            if j != i
        )
        velocities.append(beta_i**2 * (-2 + 8 * interaction))
    return velocities


def pair_term(beta_i: Fraction, beta_j: Fraction, shift: int) -> Fraction:
    return (beta_i * beta_j) ** (shift + 1) * (beta_i - beta_j) ** 2


def pair_term_derivative(
    beta_i: Fraction,
    beta_j: Fraction,
    velocity_i: Fraction,
    velocity_j: Fraction,
    shift: int,
) -> Fraction:
    term = pair_term(beta_i, beta_j, shift)
    logarithmic_derivative = (
        (shift + 1)
        * (velocity_i / beta_i + velocity_j / beta_j)
        + 2 * (velocity_i - velocity_j) / (beta_i - beta_j)
    )
    return term * logarithmic_derivative


def repeated_flux_formula(
    beta: Fraction, q: Fraction, shift: int
) -> Fraction:
    x = q / beta
    polynomial = (
        2 * x**3 * shift
        + 6 * x**3
        - 4 * x**2 * shift**2
        - 18 * x**2 * shift
        - 50 * x**2
        + 8 * x * shift**2
        + 38 * x * shift
        + 18 * x
        - 4 * shift**2
        - 22 * shift
        - 30
    )
    return beta * (16 * x ** (-shift - 1) - polynomial) / (1 - x) ** 2


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def transfer_audit() -> dict:
    cases = (
        (Fraction(1), Fraction(1, 2), Fraction(1, 4)),
        (Fraction(5, 4), Fraction(1, 2), Fraction(1, 8)),
        (Fraction(3, 2), Fraction(2, 3), Fraction(1, 6)),
        (Fraction(7, 5), Fraction(3, 5), Fraction(1, 10)),
    )
    triangular_checks = 0
    samples = []
    for beta, q, epsilon in cases:
        collision_coefficients = polynomial_from_atoms([beta, beta, q])
        split_coefficients = polynomial_from_atoms(
            [beta + epsilon, beta - epsilon, q]
        )
        expected_difference = [
            Fraction(0),
            Fraction(0),
            -epsilon**2,
            -q * epsilon**2,
        ]
        observed_difference = [
            split_coefficients[index] - collision_coefficients[index]
            for index in range(4)
        ]
        if observed_difference != expected_difference:
            raise RuntimeError("quadratic coefficient-closeness audit failed")
        for atoms, coefficients in (
            ([beta, beta, q], collision_coefficients),
            (
                [beta + epsilon, beta - epsilon, q],
                split_coefficients,
            ),
        ):
            recovered = log_moments_from_taylor(coefficients, 12)
            expected = [
                sum(atom ** (n + 1) for atom in atoms) for n in range(13)
            ]
            if recovered != expected:
                raise RuntimeError("triangular Taylor-to-Edrei audit failed")
            triangular_checks += len(expected)
        samples.append(
            {
                "beta": fraction_text(beta),
                "q": fraction_text(q),
                "epsilon": fraction_text(epsilon),
                "coefficient_difference": [
                    fraction_text(value) for value in observed_difference
                ],
            }
        )
    return {
        "cases": len(cases),
        "triangular_checks": triangular_checks,
        "samples": samples,
    }


def collision_audit() -> dict:
    beta = Fraction(1)
    q = Fraction(1, 2)
    epsilons = (Fraction(1, 4), Fraction(1, 8), Fraction(1, 32))
    determinant_checks = 0
    flux_checks = 0
    repeated_checks = 0
    samples = []

    repeated_cases = (
        (Fraction(1), Fraction(1, 2)),
        (Fraction(3, 2), Fraction(1, 2)),
        (Fraction(5, 4), Fraction(2, 3)),
    )
    for repeated_beta, repeated_q in repeated_cases:
        repeated_atoms = [
            (repeated_beta, 2 * repeated_beta),
            (repeated_q, repeated_q),
        ]
        repeated_values = edrei_moments(repeated_atoms, 15)
        repeated_derivatives = [
            flow_rhs(repeated_values, n) for n in range(15)
        ]
        for shift in range(1, 13):
            determinant = shifted_two_by_two(repeated_values, shift)
            expected_determinant = (
                2
                * (repeated_beta * repeated_q) ** (shift + 1)
                * (repeated_beta - repeated_q) ** 2
            )
            if determinant != expected_determinant:
                raise RuntimeError("repeated determinant audit failed")
            derivative = shifted_two_by_two_derivative(
                repeated_values, repeated_derivatives, shift
            )
            ratio = derivative / determinant
            expected_ratio = repeated_flux_formula(
                repeated_beta, repeated_q, shift
            )
            if ratio != expected_ratio:
                raise RuntimeError(
                    "general repeated logarithmic-flux audit failed"
                )
            if (
                repeated_beta == 1
                and repeated_q == Fraction(1, 2)
                and ratio
                != 128 * 2**shift + 4 * shift**2 + 29 * shift + 131
            ):
                raise RuntimeError(
                    "benchmark repeated logarithmic-flux audit failed"
                )
            repeated_checks += 1

    for epsilon in epsilons:
        simple_atoms = [beta + epsilon, beta - epsilon, q]
        weighted_atoms = [(atom, atom) for atom in simple_atoms]
        values = edrei_moments(weighted_atoms, 15)
        derivatives = [flow_rhs(values, n) for n in range(15)]
        velocities = simple_atom_velocities(simple_atoms)
        for n in range(13):
            velocity_derivative = sum(
                (n + 1) * atom**n * velocity
                for atom, velocity in zip(simple_atoms, velocities)
            )
            if velocity_derivative != derivatives[n]:
                raise RuntimeError("simple-root velocity audit failed")
            flux_checks += 1
        for shift in range(1, 13):
            determinant = shifted_two_by_two(values, shift)
            pair_sum = sum(
                pair_term(simple_atoms[i], simple_atoms[j], shift)
                for i in range(3)
                for j in range(i + 1, 3)
            )
            if determinant != pair_sum:
                raise RuntimeError("split Cauchy-Binet audit failed")
            determinant_checks += 1
            derivative = shifted_two_by_two_derivative(
                values, derivatives, shift
            )
            pair_derivative = sum(
                pair_term_derivative(
                    simple_atoms[i],
                    simple_atoms[j],
                    velocities[i],
                    velocities[j],
                    shift,
                )
                for i in range(3)
                for j in range(i + 1, 3)
            )
            if derivative != pair_derivative:
                raise RuntimeError("split determinant-flux audit failed")
            flux_checks += 1
        samples.append(
            {
                "epsilon": fraction_text(epsilon),
                "top_pair_product": fraction_text(
                    (beta + epsilon) * (beta - epsilon)
                ),
                "next_pair_product": fraction_text((beta + epsilon) * q),
                "shift_12_logarithmic_flux": fraction_text(
                    shifted_two_by_two_derivative(
                        values, derivatives, 12
                    )
                    / shifted_two_by_two(values, 12)
                ),
            }
        )
    return {
        "beta": fraction_text(beta),
        "q": fraction_text(q),
        "split_epsilons": [fraction_text(value) for value in epsilons],
        "repeated_cases": [
            [fraction_text(left), fraction_text(right)]
            for left, right in repeated_cases
        ],
        "determinant_checks": determinant_checks,
        "flux_checks": flux_checks,
        "repeated_checks": repeated_checks,
        "samples": samples,
    }


def conditioning_audit() -> dict:
    beta = Fraction(1)
    q = Fraction(1, 2)
    values = edrei_moments([(beta, 2 * beta), (q, q)], 43)
    rows = []
    for shift in (4, 8, 12, 20, 40):
        determinant = shifted_two_by_two(values, shift)
        scale = (
            values[shift] * values[shift + 2]
            + values[shift + 1] ** 2
        )
        condition_number = scale / determinant
        rows.append(
            {
                "shift": shift,
                "condition_number": fraction_text(condition_number),
                "condition_over_2_to_s": float(
                    condition_number / 2**shift
                ),
            }
        )
    if not 31.9 < rows[-1]["condition_over_2_to_s"] < 32.1:
        raise RuntimeError("condition-number asymptotic audit failed")
    return {
        "definition": "kappa_s=(a_s*a_(s+2)+a_(s+1)^2)/D_(1,s)",
        "asymptotic": (
            "kappa_s~[4*beta^2/(beta-q)^2]"
            "*(beta/q)^(s+1)"
        ),
        "benchmark_limit": "kappa_s/2^s -> 32 for beta=1, q=1/2",
        "rows": rows,
    }


def general_rank_audit() -> dict:
    cases = (
        ([], Fraction(1), Fraction(1, 2), Fraction(1, 8)),
        (
            [Fraction(3, 2)],
            Fraction(1),
            Fraction(1, 2),
            Fraction(1, 8),
        ),
        (
            [Fraction(2), Fraction(3, 2)],
            Fraction(1),
            Fraction(1, 2),
            Fraction(1, 10),
        ),
    )
    determinant_checks = 0
    cauchy_binet_checks = 0
    ratio_checks = 0
    samples = []
    for prefix, beta, q, epsilon in cases:
        rank = len(prefix) + 1
        collision_atoms = (
            [(atom, atom) for atom in prefix]
            + [(beta, 2 * beta), (q, q)]
        )
        split_atoms = (
            [(atom, atom) for atom in prefix]
            + [
                (beta + epsilon, beta + epsilon),
                (beta - epsilon, beta - epsilon),
                (q, q),
            ]
        )
        maximum = 8 + 2 * rank
        collision_values = edrei_moments(collision_atoms, maximum)
        split_values = edrei_moments(split_atoms, maximum)
        prefix_product = math.prod(prefix, start=Fraction(1))
        collision_support = prefix + [beta, q]
        collision_vandermonde = vandermonde_square(collision_support)
        top_split_support = prefix + [beta + epsilon, beta - epsilon]
        top_split_term_factor = vandermonde_square(top_split_support)
        for shift in range(1, 9):
            collision_determinant = shifted_determinant(
                collision_values, rank, shift
            )
            expected_collision = (
                2
                * (prefix_product * beta * q) ** (shift + 1)
                * collision_vandermonde
            )
            if collision_determinant != expected_collision:
                raise RuntimeError(
                    "general-r collision determinant audit failed"
                )
            determinant_checks += 1
            split_determinant = shifted_determinant(
                split_values, rank, shift
            )
            split_cauchy_binet = cauchy_binet_sum(
                split_atoms, rank, shift
            )
            if split_determinant != split_cauchy_binet:
                raise RuntimeError(
                    "general-r split Cauchy-Binet audit failed"
                )
            cauchy_binet_checks += 1
            top_split_term = (
                (
                    prefix_product
                    * (beta**2 - epsilon**2)
                )
                ** (shift + 1)
                * top_split_term_factor
            )
            expected_ratio = (
                Fraction(1, 2)
                * (
                    (beta**2 - epsilon**2) / (beta * q)
                )
                ** (shift + 1)
                * top_split_term_factor
                / collision_vandermonde
            )
            if top_split_term / collision_determinant != expected_ratio:
                raise RuntimeError("general-r crossover-ratio audit failed")
            ratio_checks += 1
        asymptotic_prefactor = (
            Fraction(2)
            * epsilon**2
            / (beta - q) ** 2
            * math.prod(
                (
                    (atom - beta) / (atom - q)
                ) ** 2
                for atom in prefix
            )
        )
        samples.append(
            {
                "rank": rank,
                "prefix": [fraction_text(atom) for atom in prefix],
                "beta": fraction_text(beta),
                "q": fraction_text(q),
                "epsilon": fraction_text(epsilon),
                "small_epsilon_ratio_prefactor": fraction_text(
                    asymptotic_prefactor
                ),
            }
        )
    return {
        "rank_cases": len(cases),
        "determinant_checks": determinant_checks,
        "cauchy_binet_checks": cauchy_binet_checks,
        "ratio_checks": ratio_checks,
        "samples": samples,
        "general_crossover": (
            "For any fixed larger simple prefix, "
            "s_cross=2*log(beta/epsilon)/log(beta/q)+O(1)"
        ),
    }


def crossover_audit() -> dict:
    beta = 1.0
    q = 0.5
    epsilons = (2.0**-4, 2.0**-8, 2.0**-12, 2.0**-16)
    rows = []
    for epsilon in epsilons:
        base = (beta**2 - epsilon**2) / (beta * q)
        equality_shift = (
            math.log((beta - q) ** 2 / (2 * epsilon**2))
            / math.log(base)
            - 1
        )
        asymptotic_shift = (
            2 * math.log(beta / epsilon) / math.log(beta / q)
        )
        rows.append(
            {
                "epsilon": epsilon,
                "exact_benchmark_shift": equality_shift,
                "leading_asymptotic_shift": asymptotic_shift,
                "difference": equality_shift - asymptotic_shift,
            }
        )
    if abs(rows[-1]["difference"] - rows[-2]["difference"]) > 0.02:
        raise RuntimeError("crossover asymptotic audit failed")
    return {
        "ratio": (
            "T_split/D_collision="
            "[2*epsilon^2/(beta-q)^2]"
            "*[(beta^2-epsilon^2)/(beta*q)]^(s+1)"
        ),
        "asymptotic_shift": (
            "s_cross=2*log(beta/epsilon)/log(beta/q)+O(1)"
        ),
        "rows": rows,
    }


def exact_statements() -> dict:
    return {
        "raw_moments": (
            "M_n(lambda)=integral_R exp(lambda*u^2)"
            "*Phi(u)*u^(2n)du"
        ),
        "taylor_map": (
            "H_lambda(z)=sum_(n>=0)d_n z^n, "
            "d_n=M_n/((2n)!*M_0)"
        ),
        "triangular_map": (
            "(n+1)d_(n+1)=sum_(k=0)^n(-1)^k*a_k*d_(n-k); "
            "a_n=(-1)^n[(n+1)d_(n+1)"
            "-sum_(k=0)^(n-1)(-1)^k*a_k*d_(n-k)]"
        ),
        "collision_family": (
            "H_0(z)=(1+beta*z)^2(1+q*z), "
            "H_epsilon(z)=(1+(beta+epsilon)z)"
            "*(1+(beta-epsilon)z)*(1+q*z)"
        ),
        "coefficient_gap": (
            "H_epsilon-H_0=-epsilon^2*z^2*(1+q*z)"
        ),
        "power_gap": (
            "a_n(epsilon)-a_n(0)="
            "(beta+epsilon)^(n+1)+(beta-epsilon)^(n+1)"
            "-2*beta^(n+1)"
        ),
        "collision_determinant": (
            "D_(1,s)(0)=2*(beta*q)^(s+1)*(beta-q)^2"
        ),
        "condition_number": (
            "kappa_s=(a_s*a_(s+2)+a_(s+1)^2)/D_(1,s)"
            "~[4*beta^2/(beta-q)^2]*(beta/q)^(s+1)"
        ),
        "split_determinant": (
            "D_(1,s)(epsilon)="
            "4*epsilon^2*(beta^2-epsilon^2)^(s+1)"
            "+[((beta+epsilon)*q)^(s+1)"
            "*(beta+epsilon-q)^2]"
            "+[((beta-epsilon)*q)^(s+1)"
            "*(beta-epsilon-q)^2]"
        ),
        "crossover": (
            "s_cross=2*log(beta/epsilon)/log(beta/q)+O(1)"
        ),
        "general_rank_collision": (
            "For B={b_1,...,b_(r-1)} above beta>q and "
            "P_B=product B: D_(r,s)(0)="
            "2*(P_B*beta*q)^(s+1)*V(B,beta,q)^2"
        ),
        "general_rank_top_split": (
            "T_top(epsilon)="
            "[P_B*(beta^2-epsilon^2)]^(s+1)"
            "*V(B,beta+epsilon,beta-epsilon)^2"
        ),
        "repeated_flux": (
            "With x=q/beta, partial_lambda log D_(1,s)(0)="
            "beta*[16*x^(-s-1)-P_s(x)]/(1-x)^2, "
            "P_s(x)=2*x^3*s+6*x^3-4*x^2*s^2-18*x^2*s"
            "-50*x^2+8*x*s^2+38*x*s+18*x-4*s^2-22*s-30"
        ),
        "noncommuting_limits": (
            "I_r(H)=limsup_(s->infinity)"
            "max(1,partial_lambda log D_(r,s)(H))^(1/s); "
            "after any fixed larger simple prefix, "
            "I_r(H_epsilon)=1 for epsilon>0 and "
            "I_r(H_0)=beta/q>1"
        ),
        "open_target": (
            "Control the Xi/Phi determinant flux by a structural identity "
            "that preserves the rank-one cancellation, or obtain "
            "gap-adapted exponentially small remainders; fixed algebraic "
            "raw-moment asymptotics alone do not resolve the collision."
        ),
    }


def build_payload() -> dict:
    exact = exact_statements()
    rows = [
        GateRow(
            "ermc_01_phi_raw_moment_map",
            "exact_identity",
            "available_exact",
            "The Phi transform supplies normalized Taylor coefficients by exact even raw moments.",
            exact["raw_moments"] + "; " + exact["taylor_map"],
            "This is an identity, not an asymptotic estimate.",
        ),
        GateRow(
            "ermc_02_triangular_log_map",
            "exact_identity",
            "available_exact",
            "The Taylor coefficients determine every Edrei power moment through a triangular Newton recurrence.",
            exact["triangular_map"],
            "Computing a high-shift determinant requires the recurrence through the same high order.",
        ),
        GateRow(
            "ermc_03_collision_family",
            "exact_model",
            "available_exact",
            "A symmetric split gives an LP+ family converging to a repeated factor at any fixed support rank.",
            exact["collision_family"],
            "Any fixed string of larger simple factors can be adjoined; the family is a conditioning model, not the Xi function.",
        ),
        GateRow(
            "ermc_04_quadratic_coefficient_gap",
            "exact_identity",
            "available_exact",
            "The split and collision are quadratically close in coefficient and compact-open topology.",
            exact["coefficient_gap"],
            "Multiplication by a common normalized LP+ tail preserves the O(epsilon^2) compact-open gap.",
        ),
        GateRow(
            "ermc_05_power_moment_gap",
            "exact_identity",
            "available_exact",
            "The logarithmic moments amplify the small symmetric split at high order.",
            exact["power_gap"],
            "For fixed n the gap is O(epsilon^2), but the estimate is not uniform as n grows.",
        ),
        GateRow(
            "ermc_06_collision_determinant",
            "exact_identity",
            "available_exact",
            "At collision the leading rank-one terms cancel and expose the next atom.",
            exact["collision_determinant"],
            "The formula is Cauchy-Binet for residues 2*beta and q.",
        ),
        GateRow(
            "ermc_07_exponential_condition_number",
            "exact_obstruction",
            "guard_validated",
            "Entrywise reconstruction of the collision determinant is exponentially ill-conditioned.",
            exact["condition_number"],
            "Correlated structural formulas can evade this black-box condition number; independent algebraic remainders cannot.",
        ),
        GateRow(
            "ermc_08_split_cauchy_binet",
            "exact_identity",
            "available_exact",
            "At every fixed rank the split determinant contains a tiny top-subset term with a larger exponential base.",
            exact["split_determinant"] + "; " + exact["general_rank_top_split"],
            "Every Cauchy-Binet subset term is positive.",
        ),
        GateRow(
            "ermc_09_crossover_scale",
            "exact_asymptotic",
            "available_exact",
            "At every fixed rank the high shift needed to resolve a gap diverges logarithmically as the gap closes.",
            exact["crossover"],
            "The larger simple prefix changes only the O(1) term, not the beta/q exponential base.",
        ),
        GateRow(
            "ermc_10_repeated_flux",
            "exact_identity",
            "available_exact",
            "The repeated model has an explicit exponentially growing determinant logarithmic flux.",
            exact["repeated_flux"],
            "Its sth-root limit is beta/q.",
        ),
        GateRow(
            "ermc_11_simple_split_flux",
            "exact_asymptotic",
            "available_exact",
            "For each fixed nonzero split, the largest distinct pair dominates and the logarithmic flux is only affine in s up to exponentially small terms.",
            "partial_lambda log D_(1,s)(epsilon)=A_epsilon*s+B_epsilon+o(1)",
            "The constants diverge as epsilon tends to zero, so this estimate is not uniform at collision.",
        ),
        GateRow(
            "ermc_12_noncommuting_limits",
            "exact_obstruction",
            "guard_validated",
            "Every fixed-rank all-shift collision indicator is discontinuous under compact-open LP+ convergence.",
            exact["noncommuting_limits"],
            "This blocks promotion from any fixed collection of shifts or any nonuniform coefficient asymptotic.",
        ),
        GateRow(
            "ermc_13_phi_precision_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "A Phi-specific proof must preserve the determinant cancellation structurally or resolve it exponentially.",
            exact["open_target"],
            "No such uniform Xi/Phi estimate is proved here.",
        ),
        GateRow(
            "ermc_14_proof_boundary",
            "proof_boundary",
            "guard_validated",
            "The gate separates an exact conditioning obstruction from the still-open zeta estimate.",
            "exact transfer and collision barrier proved; Xi/Phi determinant-flux bound open",
            "Nothing here proves PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_edrei_raw_moment_collision_resolution_gate",
        "date": "2026-07-24",
        "status": (
            "exact Phi-moment transfer and LP+ collision-resolution "
            "obstruction with one open Xi/Phi precision handoff"
        ),
        "exact": exact,
        "transfer_audit": transfer_audit(),
        "collision_audit": collision_audit(),
        "conditioning_audit": conditioning_audit(),
        "general_rank_audit": general_rank_audit(),
        "crossover_audit": crossover_audit(),
        "rows": [asdict(row) for row in rows],
        "sources": [
            {
                "title": "Local Edrei Hankel boundary-flux gate",
                "path": "outputs/jensen_window_pf_edrei_hankel_boundary_flux_gate.md",
            },
            {
                "title": "Local formal core, Lemma 11.22C and Corollary 11.22C.1",
                "path": "outputs/formal_core.md",
            },
            {
                "title": "Local cancellation-preserving strict-Laguerre correlation target",
                "path": "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            },
        ],
        "proof_boundary": (
            "This artifact proves the exact raw-Phi-moment to Taylor to "
            "logarithmic-moment transfer, the symmetric LP+ split identities, "
            "the exponential condition number of the high-shift determinant "
            "at a repeated leading atom, the same logarithmic crossover scale "
            "at every fixed collision rank, and the noncommuting collision "
            "limits. It does not prove a "
            "Xi/Phi-specific cancellation-preserving estimate, PF-infinity, "
            "Jensen hyperbolicity for zeta, RH, or Lambda<=0."
        ),
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    transfer = payload["transfer_audit"]
    collision = payload["collision_audit"]
    conditioning = payload["conditioning_audit"]
    general_rank = payload["general_rank_audit"]
    crossover = payload["crossover_audit"]
    return "\n".join(
        [
            "# Jensen-Window PF Edrei Raw-Moment Collision-Resolution Gate",
            "",
            "Date: 2026-07-24",
            "",
            "Status: exact Phi-moment transfer and LP+ collision-resolution",
            "obstruction with one open Xi/Phi precision handoff. This is not a",
            "proof of PF-infinity, Jensen hyperbolicity for zeta, RH, or",
            "`Lambda <= 0`.",
            "",
            "Artifact kind:",
            "`jensen_window_pf_edrei_raw_moment_collision_resolution_gate`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py",
            "```",
            "",
            "## Exact Phi-Moment Transfer",
            "",
            "Put",
            "",
            "```text",
            exact["raw_moments"],
            "F_lambda(z)=sum_(n>=0)M_n(lambda)*z^n/(2n)!",
            exact["taylor_map"],
            "```",
            "",
            "For",
            "`H_lambda'/H_lambda=sum_(n>=0)(-1)^n*a_n*z^n`,",
            "coefficient comparison in `H'=H*(H'/H)` gives",
            "",
            "```text",
            exact["triangular_map"],
            "```",
            "",
            "Thus `D_(r,s)` depends on the raw moments through at least order",
            "`s+2r+1`; it is not a local fixed-order functional of a saddle",
            "expansion. The builder performs",
            f"{transfer['triangular_checks']} exact rational triangular checks.",
            "",
            "## Symmetric Collision Family",
            "",
            "Fix `0<q<beta` and `0<epsilon<beta-q`. Consider",
            "",
            "```text",
            exact["collision_family"],
            "```",
            "",
            "Both functions are normalized LP+ polynomials. Direct multiplication",
            "gives the exact compact-open gap",
            "",
            "```text",
            exact["coefficient_gap"],
            "```",
            "",
            "Multiplying both functions by one common normalized LP+ tail whose",
            "atoms lie below `q` gives the same identity times that tail, so the",
            "obstruction is not a finite-support artifact. Their Edrei moments",
            "differ by",
            "",
            "```text",
            exact["power_gap"],
            " =2*sum_(j>=1)binom(n+1,2j)",
            "    *beta^(n+1-2j)*epsilon^(2j).",
            "```",
            "",
            "The gap is quadratic for each fixed `n`, but not uniformly small at",
            "the high orders needed by the shifted determinant.",
            "",
            "## Exponential Cancellation Condition Number",
            "",
            "At `epsilon=0`, the distinct Edrei atoms have residues `2*beta` and",
            "`q`. Cauchy-Binet gives",
            "",
            "```text",
            exact["collision_determinant"],
            "```",
            "",
            "The two products subtracted in",
            "`D_(1,s)=a_s*a_(s+2)-a_(s+1)^2`",
            "are each of order `beta^(2s)`. With the natural componentwise",
            "subtraction condition number",
            "",
            "```text",
            exact["condition_number"],
            "```",
            "",
            "a worst-case entrywise relative remainder must therefore be",
            "`o((q/beta)^s)` to resolve the collision determinant uniformly.",
            "Any fixed algebraic remainder `O(s^-N)` loses this comparison.",
            "For `beta=1,q=1/2`, the exact audit confirms",
            f"`{conditioning['benchmark_limit']}`.",
            "",
            "This is a black-box reconstruction barrier. A correlated integral",
            "identity, sum-of-squares formula, or operator inequality that preserves",
            "the cancellation can bypass it; mere independent error bars cannot.",
            "",
            "## Split-Pair Crossover",
            "",
            "For `epsilon>0`, Cauchy-Binet gives",
            "",
            "```text",
            exact["split_determinant"],
            "```",
            "",
            "The first term has the largest exponential base but coefficient",
            "`4*epsilon^2`. Relative to the collision determinant,",
            "",
            "```text",
            crossover["ratio"],
            "```",
            "",
            "so its benchmark crossover is",
            "",
            "```text",
            exact["crossover"],
            "```",
            "",
            f"The builder checks {len(crossover['rows'])} decreasing gaps.",
            "Every finite shift range can therefore be made collision-like by",
            "taking a sufficiently small but nonzero split.",
            "",
            "### Arbitrary Fixed Rank",
            "",
            "Let `B={b_1,...,b_(r-1)}` be any fixed string of larger simple",
            "atoms, `b_1>...>b_(r-1)>beta>q`, and put",
            "`P_B=product_(b in B)b`. Multiplying both collision polynomials",
            "by `product_(b in B)(1+b*z)` moves the repeated factor to rank",
            "`r`. Cauchy-Binet gives",
            "",
            "```text",
            exact["general_rank_collision"],
            exact["general_rank_top_split"],
            "```",
            "",
            "The top split term divided by the collision determinant is",
            "",
            "```text",
            "(1/2)*[(beta^2-epsilon^2)/(beta*q)]^(s+1)",
            " *V(B,beta+epsilon,beta-epsilon)^2/V(B,beta,q)^2.",
            "```",
            "",
            "As `epsilon->0`, the Vandermonde ratio is",
            "",
            "```text",
            "4*epsilon^2/(beta-q)^2",
            " *product_(b in B)[(b-beta)/(b-q)]^2*(1+O(epsilon^2)).",
            "```",
            "",
            "The prefix therefore changes only the bounded prefactor. The",
            "crossover remains",
            "`2*log(beta/epsilon)/log(beta/q)+O(1)` at every fixed rank.",
            f"The builder checks {general_rank['rank_cases']} ranks,",
            f"{general_rank['determinant_checks']} exact collision determinants,",
            f"{general_rank['cauchy_binet_checks']} split Cauchy-Binet",
            f"identities, and {general_rank['ratio_checks']} exact top-term",
            "ratios.",
            "",
            "## Noncommuting Flux Limits",
            "",
            "The exact repeated-model logarithmic flux is",
            "",
            "```text",
            exact["repeated_flux"],
            "```",
            "",
            "For `beta=1,q=1/2` this reduces to",
            "",
            "```text",
            "partial_lambda log D_(1,s)",
            " =128*2^s+4*s^2+29*s+131.",
            "```",
            "",
            "For every fixed `epsilon>0`, the atoms are distinct. Differentiating",
            "their positive Cauchy-Binet pair sum under radial heat shows that the",
            "largest pair contributes an affine function of `s`, while all other",
            "pairs are exponentially smaller. Hence",
            "",
            "```text",
            "partial_lambda log D_(1,s)(epsilon)",
            " =A_epsilon*s+B_epsilon+o(1),",
            exact["noncommuting_limits"],
            "```",
            "",
            "For each fixed `s`, the coefficient and heat-flow formulas are",
            "polynomial in `epsilon^2`, so the split flux tends to the repeated",
            "flux. The limits `epsilon->0` and `s->infinity` do not commute.",
            f"The builder checks {collision['determinant_checks']} exact split",
            f"determinants, {collision['flux_checks']} exact velocity/flux",
            f"identities, and {collision['repeated_checks']} repeated benchmarks.",
            "",
            "## Proof Consequence",
            "",
            "This artifact proves:",
            "",
            "- the exact `Phi raw moments -> Taylor coefficients -> Edrei moments`",
            "  transfer;",
            "- an exponentially ill-conditioned collision determinant;",
            "- a logarithmically escaping split-resolution shift;",
            "- a compact-open LP+ family with noncommuting collision limits.",
            "",
            "It does not prove:",
            "",
            "- the Xi/Phi subexponential determinant-flux estimate;",
            "- all-order Stieltjes positivity at `lambda=0`;",
            "- PF-infinity, Jensen hyperbolicity for zeta, RH, or `Lambda<=0`.",
            "",
            "The surviving handoff is precise:",
            "",
            "```text",
            exact["open_target"],
            "```",
            "",
            "The existing cancellation-preserving Phi formulation is the strict",
            "Laguerre correlation target",
            "",
            "```text",
            "L_t(x)=H_t'(x)^2-H_t(x)H_t''(x)",
            "      =Fourier[K_(1,t)](2x).",
            "```",
            "",
            "It is recorded in",
            "`outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md`.",
            "The present gate explains why returning to that correlated kernel is",
            "structurally preferable to extending entrywise raw-moment expansions;",
            "it does not prove `L_t(x)>0`.",
            "",
            "A finite Poincare saddle expansion with an algebraic remainder is",
            "not, by itself, strong enough at a possible collision.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Edrei raw-moment collision-resolution gate: "
        f"{len(payload['rows'])} rows, "
        f"{payload['transfer_audit']['triangular_checks']} transfer audits, "
        f"{payload['collision_audit']['determinant_checks']} determinant audits, "
        f"{payload['collision_audit']['flux_checks']} flux audits, "
        f"{payload['general_rank_audit']['determinant_checks']} "
        "general-r determinant audits"
    )


if __name__ == "__main__":
    main()
