#!/usr/bin/env python3
"""Independently validate the Edrei raw-moment collision-resolution gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_edrei_raw_moment_collision_resolution_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.md"
)
DEFAULT_FORMAL = REPO_ROOT / "outputs/formal_core.md"

REQUIRED_IDS = {
    f"ermc_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "phi_raw_moment_map"),
        (2, "triangular_log_map"),
        (3, "collision_family"),
        (4, "quadratic_coefficient_gap"),
        (5, "power_moment_gap"),
        (6, "collision_determinant"),
        (7, "exponential_condition_number"),
        (8, "split_cauchy_binet"),
        (9, "crossover_scale"),
        (10, "repeated_flux"),
        (11, "simple_split_flux"),
        (12, "noncommuting_limits"),
        (13, "phi_precision_handoff"),
        (14, "proof_boundary"),
    )
}

REQUIRED_NOTE_TEXT = (
    "d_n=M_n/((2n)!*M_0)",
    "(n+1)d_(n+1)=sum_(k=0)^n(-1)^k*a_k*d_(n-k)",
    "H_epsilon-H_0=-epsilon^2*z^2*(1+q*z)",
    "D_(1,s)(0)=2*(beta*q)^(s+1)*(beta-q)^2",
    "kappa_s=(a_s*a_(s+2)+a_(s+1)^2)/D_(1,s)",
    "s_cross=2*log(beta/epsilon)/log(beta/q)+O(1)",
    "### Arbitrary Fixed Rank",
    "V(B,beta+epsilon,beta-epsilon)^2/V(B,beta,q)^2",
    "128*2^s+4*s^2+29*s+131",
    "The limits `epsilon->0` and `s->infinity` do not commute.",
    "This artifact proves:",
    "It does not prove:",
    "fixed algebraic remainder",
)

REQUIRED_FORMAL_TEXT = (
    "#### Corollary 11.22C.2: Raw-Moment Collision-Resolution Barrier",
    "d_n(lambda)=M_n(lambda)/((2n)!M_0(lambda))",
    "H_epsilon(z)-H_0(z)=-epsilon^2 z^2(1+qz).",
    "kappa_s",
    "s_cross",
    "I_1(H_epsilon)=1  for every epsilon>0",
    "T_top(epsilon)/D_(r,s)(0)",
    "commute at any fixed collision rank",
    "outputs/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.md",
)


def polynomial_from_atoms(atoms: list[Fraction]) -> list[Fraction]:
    coefficients = [Fraction(1)]
    for atom in atoms:
        out = [Fraction(0)] * (len(coefficients) + 1)
        for index, value in enumerate(coefficients):
            out[index] += value
            out[index + 1] += atom * value
        coefficients = out
    return coefficients


def recover_log_moments(
    coefficients: list[Fraction], maximum: int
) -> list[Fraction]:
    padded = coefficients + [Fraction(0)] * (maximum + 2)
    recovered: list[Fraction] = []
    for n in range(maximum + 1):
        lower = sum(
            (-1) ** k * recovered[k] * padded[n - k]
            for k in range(n)
        )
        recovered.append(
            (-1) ** n * ((n + 1) * padded[n + 1] - lower)
        )
    return recovered


def moments(
    weighted_atoms: list[tuple[Fraction, Fraction]], maximum: int
) -> list[Fraction]:
    return [
        sum(residue * atom**n for atom, residue in weighted_atoms)
        for n in range(maximum + 1)
    ]


def flow(values: list[Fraction], n: int) -> Fraction:
    return (
        -2 * (n + 1) * (2 * n + 3) * values[n + 1]
        + 4
        * (n + 1)
        * sum(values[k] * values[n - k] for k in range(n + 1))
    )


def determinant(values: list[Fraction], shift: int) -> Fraction:
    return (
        values[shift] * values[shift + 2] - values[shift + 1] ** 2
    )


def matrix_determinant(matrix: list[list[Fraction]]) -> Fraction:
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


def rank_determinant(
    values: list[Fraction], rank: int, shift: int
) -> Fraction:
    return matrix_determinant(
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


def cauchy_binet(
    weighted_atoms: list[tuple[Fraction, Fraction]],
    rank: int,
    shift: int,
) -> Fraction:
    total = Fraction(0)
    for selected_indices in combinations(
        range(len(weighted_atoms)), rank + 1
    ):
        selected = [weighted_atoms[index] for index in selected_indices]
        selected_atoms = [atom for atom, _ in selected]
        weight = Fraction(1)
        for atom, residue in selected:
            weight *= residue * atom**shift
        total += weight * vandermonde_square(selected_atoms)
    return total


def determinant_derivative(
    values: list[Fraction], derivatives: list[Fraction], shift: int
) -> Fraction:
    return (
        derivatives[shift] * values[shift + 2]
        + values[shift] * derivatives[shift + 2]
        - 2 * values[shift + 1] * derivatives[shift + 1]
    )


def velocities(atoms: list[Fraction]) -> list[Fraction]:
    return [
        atom_i**2
        * (
            -2
            + 8
            * sum(
                atoms[j] / (atom_i - atoms[j])
                for j in range(len(atoms))
                if j != i
            )
        )
        for i, atom_i in enumerate(atoms)
    ]


def pair_term(left: Fraction, right: Fraction, shift: int) -> Fraction:
    return (left * right) ** (shift + 1) * (left - right) ** 2


def pair_derivative(
    left: Fraction,
    right: Fraction,
    left_velocity: Fraction,
    right_velocity: Fraction,
    shift: int,
) -> Fraction:
    return pair_term(left, right, shift) * (
        (shift + 1)
        * (left_velocity / left + right_velocity / right)
        + 2 * (left_velocity - right_velocity) / (left - right)
    )


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


def transfer_checks() -> tuple[list[str], int]:
    issues: list[str] = []
    checks = 0
    cases = (
        (Fraction(1), Fraction(2, 5), Fraction(1, 7)),
        (Fraction(4, 3), Fraction(1, 2), Fraction(1, 9)),
        (Fraction(7, 5), Fraction(3, 5), Fraction(1, 10)),
        (Fraction(3, 2), Fraction(3, 4), Fraction(1, 8)),
    )
    for beta, q, epsilon in cases:
        collision = polynomial_from_atoms([beta, beta, q])
        split = polynomial_from_atoms(
            [beta + epsilon, beta - epsilon, q]
        )
        difference = [
            split[index] - collision[index] for index in range(4)
        ]
        if difference != [
            Fraction(0),
            Fraction(0),
            -epsilon**2,
            -q * epsilon**2,
        ]:
            issues.append(
                f"coefficient gap failed for beta={beta}, q={q}, "
                f"epsilon={epsilon}"
            )
        for atoms, coefficients in (
            ([beta, beta, q], collision),
            ([beta + epsilon, beta - epsilon, q], split),
        ):
            recovered = recover_log_moments(coefficients, 12)
            expected = [
                sum(atom ** (n + 1) for atom in atoms)
                for n in range(13)
            ]
            for n, (observed, target) in enumerate(
                zip(recovered, expected)
            ):
                if observed != target:
                    issues.append(
                        f"triangular recurrence failed at n={n}, "
                        f"atoms={atoms}: {observed} != {target}"
                    )
                checks += 1
    return issues, checks


def collision_checks() -> tuple[list[str], int, int, int]:
    issues: list[str] = []
    determinant_checks = 0
    flux_checks = 0
    repeated_checks = 0
    beta = Fraction(1)
    q = Fraction(1, 2)

    for repeated_beta, repeated_q in (
        (Fraction(1), Fraction(1, 2)),
        (Fraction(4, 3), Fraction(2, 5)),
        (Fraction(7, 5), Fraction(3, 5)),
    ):
        repeated_values = moments(
            [
                (repeated_beta, 2 * repeated_beta),
                (repeated_q, repeated_q),
            ],
            15,
        )
        repeated_flow = [flow(repeated_values, n) for n in range(15)]
        for shift in range(1, 13):
            observed_determinant = determinant(repeated_values, shift)
            expected_determinant = (
                2
                * (repeated_beta * repeated_q) ** (shift + 1)
                * (repeated_beta - repeated_q) ** 2
            )
            if observed_determinant != expected_determinant:
                issues.append(
                    f"repeated determinant failed at beta={repeated_beta}, "
                    f"q={repeated_q}, s={shift}"
                )
            observed_ratio = (
                determinant_derivative(
                    repeated_values, repeated_flow, shift
                )
                / observed_determinant
            )
            expected_ratio = repeated_flux_formula(
                repeated_beta, repeated_q, shift
            )
            if observed_ratio != expected_ratio:
                issues.append(
                    f"repeated flux failed at beta={repeated_beta}, "
                    f"q={repeated_q}, s={shift}: "
                    f"{observed_ratio} != {expected_ratio}"
                )
            if (
                repeated_beta == 1
                and repeated_q == Fraction(1, 2)
                and observed_ratio
                != 128 * 2**shift
                + 4 * shift**2
                + 29 * shift
                + 131
            ):
                issues.append(
                    f"benchmark repeated flux failed at s={shift}"
                )
            repeated_checks += 1

    for epsilon in (
        Fraction(1, 5),
        Fraction(1, 10),
        Fraction(1, 40),
    ):
        atoms = [beta + epsilon, beta - epsilon, q]
        values = moments([(atom, atom) for atom in atoms], 15)
        derivatives = [flow(values, n) for n in range(15)]
        root_velocities = velocities(atoms)
        for n in range(13):
            differentiated = sum(
                (n + 1) * atoms[index] ** n * root_velocities[index]
                for index in range(3)
            )
            if differentiated != derivatives[n]:
                issues.append(
                    f"root-velocity moment flow failed at epsilon={epsilon}, "
                    f"n={n}"
                )
            flux_checks += 1
        for shift in range(1, 13):
            pair_sum = sum(
                pair_term(atoms[i], atoms[j], shift)
                for i in range(3)
                for j in range(i + 1, 3)
            )
            if determinant(values, shift) != pair_sum:
                issues.append(
                    f"split determinant failed at epsilon={epsilon}, "
                    f"s={shift}"
                )
            determinant_checks += 1
            pair_flux = sum(
                pair_derivative(
                    atoms[i],
                    atoms[j],
                    root_velocities[i],
                    root_velocities[j],
                    shift,
                )
                for i in range(3)
                for j in range(i + 1, 3)
            )
            direct_flux = determinant_derivative(
                values, derivatives, shift
            )
            if pair_flux != direct_flux:
                issues.append(
                    f"split pair flux failed at epsilon={epsilon}, "
                    f"s={shift}"
                )
            flux_checks += 1
        if not atoms[0] * atoms[1] > atoms[0] * atoms[2]:
            issues.append("split leading-pair product is not dominant")
    return issues, determinant_checks, flux_checks, repeated_checks


def conditioning_checks() -> tuple[list[str], int]:
    issues: list[str] = []
    beta = Fraction(1)
    q = Fraction(1, 2)
    values = moments([(beta, 2 * beta), (q, q)], 43)
    shifts = (4, 8, 12, 20, 40)
    scaled = []
    for shift in shifts:
        det = determinant(values, shift)
        scale = (
            values[shift] * values[shift + 2]
            + values[shift + 1] ** 2
        )
        condition = scale / det
        scaled.append(float(condition / 2**shift))
    if any(right >= left for left, right in zip(scaled, scaled[1:])):
        issues.append("scaled condition numbers do not decrease to the limit")
    if abs(scaled[-1] - 32.0) > 0.1:
        issues.append(
            f"scaled condition-number limit is {scaled[-1]}, expected 32"
        )
    return issues, len(shifts)


def general_rank_checks() -> tuple[list[str], int, int, int]:
    issues: list[str] = []
    determinant_checks = 0
    cauchy_checks = 0
    ratio_checks = 0
    cases = (
        ([], Fraction(6, 5), Fraction(2, 5), Fraction(1, 9)),
        (
            [Fraction(7, 4)],
            Fraction(6, 5),
            Fraction(2, 5),
            Fraction(1, 9),
        ),
        (
            [Fraction(9, 4), Fraction(7, 4)],
            Fraction(6, 5),
            Fraction(2, 5),
            Fraction(1, 10),
        ),
    )
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
        collision_values = moments(collision_atoms, maximum)
        split_values = moments(split_atoms, maximum)
        prefix_product = Fraction(1)
        for atom in prefix:
            prefix_product *= atom
        collision_support = prefix + [beta, q]
        collision_vandermonde = vandermonde_square(collision_support)
        top_support = prefix + [beta + epsilon, beta - epsilon]
        top_vandermonde = vandermonde_square(top_support)
        for shift in range(1, 9):
            observed_collision = rank_determinant(
                collision_values, rank, shift
            )
            expected_collision = (
                2
                * (prefix_product * beta * q) ** (shift + 1)
                * collision_vandermonde
            )
            if observed_collision != expected_collision:
                issues.append(
                    f"general-r collision determinant failed at "
                    f"r={rank}, s={shift}"
                )
            determinant_checks += 1
            observed_split = rank_determinant(
                split_values, rank, shift
            )
            expected_split = cauchy_binet(
                split_atoms, rank, shift
            )
            if observed_split != expected_split:
                issues.append(
                    f"general-r Cauchy-Binet failed at r={rank}, "
                    f"s={shift}"
                )
            cauchy_checks += 1
            top_term = (
                (
                    prefix_product
                    * (beta**2 - epsilon**2)
                )
                ** (shift + 1)
                * top_vandermonde
            )
            expected_ratio = (
                Fraction(1, 2)
                * (
                    (beta**2 - epsilon**2) / (beta * q)
                )
                ** (shift + 1)
                * top_vandermonde
                / collision_vandermonde
            )
            if top_term / observed_collision != expected_ratio:
                issues.append(
                    f"general-r top-term ratio failed at r={rank}, "
                    f"s={shift}"
                )
            ratio_checks += 1
    return issues, determinant_checks, cauchy_checks, ratio_checks


def crossover_checks() -> tuple[list[str], int]:
    issues: list[str] = []
    beta = 1.0
    q = 0.5
    epsilons = (2.0**-4, 2.0**-8, 2.0**-12, 2.0**-16)
    differences = []
    for epsilon in epsilons:
        base = (beta**2 - epsilon**2) / (beta * q)
        exact_shift = (
            math.log((beta - q) ** 2 / (2 * epsilon**2))
            / math.log(base)
            - 1
        )
        leading_shift = (
            2 * math.log(beta / epsilon) / math.log(beta / q)
        )
        differences.append(exact_shift - leading_shift)
    if abs(differences[-1] - differences[-2]) > 0.02:
        issues.append("crossover correction does not approach a constant")
    return issues, len(epsilons)


def artifact_checks(payload: dict, note: str, formal: str) -> list[str]:
    issues: list[str] = []
    if (
        payload.get("kind")
        != "jensen_window_pf_edrei_raw_moment_collision_resolution_gate"
    ):
        issues.append("wrong artifact kind")
    rows = payload.get("rows")
    if not isinstance(rows, list):
        return issues + ["rows is not a list"]
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_IDS:
        issues.append(
            f"row IDs differ: missing={sorted(REQUIRED_IDS - ids)}, "
            f"extra={sorted(ids - REQUIRED_IDS)}"
        )
    open_rows = [
        row
        for row in rows
        if row.get("readiness") == "not_ready_to_apply"
    ]
    if len(open_rows) != 1:
        issues.append(f"expected one open handoff, found {len(open_rows)}")
    elif open_rows[0].get("id") != "ermc_13_phi_precision_handoff":
        issues.append("wrong row marked as the open handoff")
    for snippet in REQUIRED_NOTE_TEXT:
        if snippet not in note:
            issues.append(f"missing note text: {snippet!r}")
    for snippet in REQUIRED_FORMAL_TEXT:
        if snippet not in formal:
            issues.append(f"missing formal-core text: {snippet!r}")
    proof_boundary = payload.get("proof_boundary", "")
    for guard in ("does not prove", "RH", "Lambda<=0"):
        if guard not in proof_boundary:
            issues.append(f"proof boundary missing {guard!r}")
    transfer = payload.get("transfer_audit", {})
    collision = payload.get("collision_audit", {})
    if transfer.get("triangular_checks") != 104:
        issues.append("builder transfer-audit count is not 104")
    if collision.get("determinant_checks") != 36:
        issues.append("builder determinant-audit count is not 36")
    if collision.get("flux_checks") != 75:
        issues.append("builder flux-audit count is not 75")
    if collision.get("repeated_checks") != 36:
        issues.append("builder repeated-audit count is not 36")
    general_rank = payload.get("general_rank_audit", {})
    if general_rank.get("rank_cases") != 3:
        issues.append("builder general-r rank count is not 3")
    if general_rank.get("determinant_checks") != 24:
        issues.append("builder general-r determinant count is not 24")
    if general_rank.get("cauchy_binet_checks") != 24:
        issues.append("builder general-r Cauchy-Binet count is not 24")
    if general_rank.get("ratio_checks") != 24:
        issues.append("builder general-r ratio count is not 24")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--formal", type=Path, default=DEFAULT_FORMAL)
    args = parser.parse_args()

    payload = json.loads(args.artifact.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    formal = args.formal.read_text(encoding="utf-8")
    issues = artifact_checks(payload, note, formal)

    transfer_issues, transfer_count = transfer_checks()
    (
        collision_issues,
        determinant_count,
        flux_count,
        repeated_count,
    ) = collision_checks()
    condition_issues, condition_count = conditioning_checks()
    (
        general_rank_issues,
        general_rank_determinants,
        general_rank_cauchy,
        general_rank_ratios,
    ) = general_rank_checks()
    crossover_issues, crossover_count = crossover_checks()
    issues.extend(transfer_issues)
    issues.extend(collision_issues)
    issues.extend(condition_issues)
    issues.extend(general_rank_issues)
    issues.extend(crossover_issues)

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        raise SystemExit(1)

    print(
        "validated Edrei raw-moment collision-resolution gate: "
        f"{len(payload['rows'])} rows, 0 issues, "
        f"{transfer_count} triangular transfer audits, "
        f"{determinant_count} split determinant audits, "
        f"{flux_count} simple-root/flux audits, "
        f"{repeated_count} repeated-flux audits, "
        f"{condition_count} conditioning audits, "
        f"{general_rank_determinants} general-r collision determinant audits, "
        f"{general_rank_cauchy} general-r Cauchy-Binet audits, "
        f"{general_rank_ratios} general-r crossover-ratio audits, "
        f"{crossover_count} crossover audits, "
        "1 noncommuting-limit obstruction, "
        "1 open Xi/Phi precision handoff"
    )


if __name__ == "__main__":
    main()
