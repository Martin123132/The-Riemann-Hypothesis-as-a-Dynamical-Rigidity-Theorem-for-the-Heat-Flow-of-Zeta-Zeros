#!/usr/bin/env python3
"""Independently validate the Edrei Hankel boundary-flux gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_edrei_hankel_boundary_flux_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_edrei_hankel_boundary_flux_gate.md"
)

REQUIRED_IDS = {
    f"ehbf_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "moment_flow"),
        (2, "quadratic_form"),
        (3, "operator_flux"),
        (4, "integer_edrei_measure"),
        (5, "null_form_flux"),
        (6, "integer_residue_orientation"),
        (7, "vandermonde_determinant_flux"),
        (8, "rank_one_recovery"),
        (9, "fractional_residue_guard"),
        (10, "infinite_support_strictness"),
        (11, "canonical_order_escape"),
        (12, "cauchy_binet_expansion"),
        (13, "high_shift_collision_sensor"),
        (14, "simple_two_atom_guard"),
        (15, "repeated_two_atom_witness"),
        (16, "positive_newman_consequence"),
        (17, "subexponential_phi_target"),
        (18, "proof_boundary"),
    )
}

REQUIRED_NOTE_TEXT = (
    "Q_s'(P)=-2 integral x(E+1)(2E+3)",
    "Q_s'(P)=8 sum_j rho_j*beta_j^(s+2)*(rho_j-beta_j)",
    "m_j*(m_j-1)*beta_j^(s+4)",
    "optional `gamma*delta_0` atom contributes only additional nonnegative",
    "canonical Stieltjes witness order d->infinity",
    "h_(r,s)=Q_s(pi_(r,s))=D_(r,s)/D_(r-1,s)",
    "pi_(r,s)(beta_i)=O((beta_(r+1)/beta_i)^s)",
    "lim_(s->infinity)(partial_lambda log D_(r,s))^(1/s)",
    "128*2**s + 4*s**2 + 29*s + 131",
    "Xi/Phi subexponential determinant-growth estimate",
    "This artifact proves",
    "It does not prove",
)


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


def moments(
    atoms: list[tuple[Fraction, Fraction]], maximum: int
) -> list[Fraction]:
    return [
        sum(rho * beta**n for beta, rho in atoms)
        for n in range(maximum + 1)
    ]


def flow(values: list[Fraction], n: int) -> Fraction:
    return (
        -2 * (n + 1) * (2 * n + 3) * values[n + 1]
        + 4
        * (n + 1)
        * sum(values[k] * values[n - k] for k in range(n + 1))
    )


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    sign = 1
    result = Fraction(1)
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
        result *= pivot_value
        for row in range(column + 1, len(work)):
            ratio = work[row][column] / pivot_value
            for entry in range(column + 1, len(work)):
                work[row][entry] -= ratio * work[column][entry]
    return sign * result


def determinant_derivative(
    matrix: list[list[Fraction]], derivative: list[list[Fraction]]
) -> Fraction:
    total = Fraction(0)
    for row_index in range(len(matrix)):
        replaced = [row[:] for row in matrix]
        replaced[row_index] = derivative[row_index][:]
        total += determinant(replaced)
    return total


def solve_linear(
    matrix: list[list[Fraction]], rhs: list[Fraction]
) -> list[Fraction]:
    size = len(matrix)
    work = [matrix[row][:] + [rhs[row]] for row in range(size)]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            raise ValueError("singular orthogonality system")
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    work[row][entry] - factor * work[column][entry]
                    for entry in range(size + 1)
                ]
    return [work[row][-1] for row in range(size)]


def null_form_checks() -> tuple[list[str], int]:
    issues: list[str] = []
    checks = 0
    beta_sets = (
        [Fraction(7, 5)],
        [Fraction(7, 5), Fraction(4, 5)],
        [Fraction(7, 5), Fraction(4, 5), Fraction(1, 3)],
    )
    multiplicity_sets = ((1,), (2,), (1, 3), (2, 1), (1, 2, 4))
    for betas in beta_sets:
        for multiplicities in multiplicity_sets:
            if len(multiplicities) != len(betas):
                continue
            atoms = [
                (beta, Fraction(multiplicity) * beta)
                for beta, multiplicity in zip(betas, multiplicities)
            ]
            coefficients = annihilator(betas)
            degree = len(betas)
            for shift in range(1, 6):
                values = moments(atoms, shift + 2 * degree + 1)
                derivatives = [
                    flow(values, n) for n in range(shift + 2 * degree + 1)
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
                    issues.append(
                        f"null-form mismatch betas={betas} m={multiplicities} "
                        f"s={shift}: {observed} != {expected}"
                    )
                if expected < 0:
                    issues.append("integer-residue boundary flux is negative")
                checks += 1
    return issues, checks


def determinant_checks() -> tuple[list[str], int]:
    issues: list[str] = []
    checks = 0
    cases = (
        ([Fraction(7, 5)], [2]),
        ([Fraction(7, 5), Fraction(4, 5)], [1, 3]),
        (
            [Fraction(7, 5), Fraction(4, 5), Fraction(1, 3)],
            [1, 2, 4],
        ),
    )
    for betas, multiplicities in cases:
        atoms = [
            (beta, Fraction(multiplicity) * beta)
            for beta, multiplicity in zip(betas, multiplicities)
        ]
        coefficients = annihilator(betas)
        rank = len(betas)
        for shift in range(1, 5):
            values = moments(atoms, shift + 2 * rank + 1)
            derivatives = [
                flow(values, n) for n in range(shift + 2 * rank + 1)
            ]
            size = rank + 1
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
            expected_form = 8 * sum(
                rho
                * beta ** (shift + 2)
                * (rho - beta)
                * derivative_value(coefficients, beta) ** 2
                for beta, rho in atoms
            )
            expected = tau * expected_form
            if determinant(hankel) != 0:
                issues.append("finite-rank boundary determinant is nonzero")
            if observed != expected:
                issues.append(
                    f"determinant-flux mismatch rank={rank} s={shift}: "
                    f"{observed} != {expected}"
                )
            checks += 1
    return issues, checks


def cauchy_binet_checks() -> tuple[list[str], int]:
    issues: list[str] = []
    checks = 0
    atoms = [
        (Fraction(5, 4), Fraction(5, 4)),
        (Fraction(3, 4), Fraction(3, 2)),
        (Fraction(2, 5), Fraction(6, 5)),
        (Fraction(1, 6), Fraction(2, 3)),
    ]
    from itertools import combinations

    for degree in (0, 1, 2):
        size = degree + 1
        for shift in range(5):
            values = moments(atoms, shift + 2 * degree)
            observed = determinant(
                [
                    [values[shift + i + j] for j in range(size)]
                    for i in range(size)
                ]
            )
            expected = Fraction(0)
            for subset in combinations(atoms, size):
                weight = Fraction(1)
                betas = []
                for beta, rho in subset:
                    betas.append(beta)
                    weight *= rho * beta**shift
                for i, beta_i in enumerate(betas):
                    for beta_j in betas[i + 1 :]:
                        weight *= (beta_i - beta_j) ** 2
                expected += weight
            if observed != expected:
                issues.append(
                    f"Cauchy-Binet mismatch d={degree} s={shift}: "
                    f"{observed} != {expected}"
                )
            if observed <= 0:
                issues.append(f"positive determinant failed d={degree} s={shift}")
            checks += 1
    return issues, checks


def orthogonal_quotient_checks() -> tuple[list[str], int]:
    issues: list[str] = []
    checks = 0
    atom_sets = (
        [
            (Fraction(5, 4), Fraction(5, 4)),
            (Fraction(3, 4), Fraction(3, 2)),
            (Fraction(2, 5), Fraction(6, 5)),
            (Fraction(1, 6), Fraction(2, 3)),
        ],
        [
            (Fraction(3, 2), Fraction(3)),
            (Fraction(5, 6), Fraction(5, 6)),
            (Fraction(1, 2), Fraction(3, 2)),
            (Fraction(1, 5), Fraction(1, 5)),
            (Fraction(1, 10), Fraction(1, 5)),
        ],
    )
    for atoms in atom_sets:
        for rank in (1, 2, 3):
            for shift in range(1, 5):
                maximum = shift + 2 * rank + 1
                values = moments(atoms, maximum)
                derivatives = [flow(values, n) for n in range(maximum)]

                matrix = [
                    [values[shift + i + j] for j in range(rank + 1)]
                    for i in range(rank + 1)
                ]
                matrix_dot = [
                    [derivatives[shift + i + j] for j in range(rank + 1)]
                    for i in range(rank + 1)
                ]
                previous = [
                    [values[shift + i + j] for j in range(rank)]
                    for i in range(rank)
                ]
                previous_dot = [
                    [derivatives[shift + i + j] for j in range(rank)]
                    for i in range(rank)
                ]

                determinant_value = determinant(matrix)
                determinant_dot = determinant_derivative(matrix, matrix_dot)
                previous_value = determinant(previous)
                previous_derivative = determinant_derivative(
                    previous, previous_dot
                )

                lower_coefficients = solve_linear(
                    previous,
                    [-values[shift + rank + i] for i in range(rank)],
                )
                coefficients = lower_coefficients + [Fraction(1)]
                quadratic_value = sum(
                    coefficients[i]
                    * coefficients[j]
                    * values[shift + i + j]
                    for i in range(rank + 1)
                    for j in range(rank + 1)
                )
                quadratic_dot = sum(
                    coefficients[i]
                    * coefficients[j]
                    * derivatives[shift + i + j]
                    for i in range(rank + 1)
                    for j in range(rank + 1)
                )

                expected_value = determinant_value / previous_value
                expected_dot = (
                    determinant_dot * previous_value
                    - determinant_value * previous_derivative
                ) / previous_value**2
                if quadratic_value != expected_value:
                    issues.append(
                        f"orthogonal quotient mismatch rank={rank} s={shift}: "
                        f"{quadratic_value} != {expected_value}"
                    )
                if quadratic_dot != expected_dot:
                    issues.append(
                        f"orthogonal flux mismatch rank={rank} s={shift}: "
                        f"{quadratic_dot} != {expected_dot}"
                    )
                checks += 1
    return issues, checks


def two_atom_checks() -> tuple[list[str], int]:
    issues: list[str] = []
    checks = 0
    q = Fraction(1, 2)
    for leading_multiplicity in (1, 2):
        atoms = [
            (Fraction(1), Fraction(leading_multiplicity)),
            (q, q),
        ]
        for shift in range(21):
            values = moments(atoms, shift + 3)
            derivatives = [flow(values, n) for n in range(shift + 3)]
            det_value = (
                values[shift] * values[shift + 2] - values[shift + 1] ** 2
            )
            det_dot = (
                derivatives[shift] * values[shift + 2]
                + values[shift] * derivatives[shift + 2]
                - 2 * values[shift + 1] * derivatives[shift + 1]
            )
            observed = det_dot / det_value
            if leading_multiplicity == 1:
                expected = 39 - 3 * shift
            else:
                expected = 128 * 2**shift + 4 * shift**2 + 29 * shift + 131
            if observed != expected:
                issues.append(
                    f"two-atom mismatch m={leading_multiplicity} s={shift}: "
                    f"{observed} != {expected}"
                )
            checks += 1
    fractional_flux = (
        8
        * Fraction(1, 2)
        * Fraction(1) ** 3
        * (Fraction(1, 2) - Fraction(1))
    )
    if fractional_flux != -2:
        issues.append("fractional-residue guard mismatch")
    checks += 1
    return issues, checks


def validate(artifact: Path, note: Path) -> tuple[list[str], dict]:
    issues: list[str] = []
    counts: dict[str, int] = {}
    if not artifact.exists():
        return [f"missing artifact: {artifact}"], counts
    if not note.exists():
        return [f"missing note: {note}"], counts
    payload = json.loads(artifact.read_text(encoding="utf-8"))
    if payload.get("kind") != "jensen_window_pf_edrei_hankel_boundary_flux_gate":
        issues.append("unexpected artifact kind")
    rows = payload.get("rows", [])
    if len(rows) != 18:
        issues.append(f"expected 18 rows, found {len(rows)}")
    ids = {row.get("id") for row in rows}
    if ids != REQUIRED_IDS:
        issues.append(
            f"row id mismatch: missing={sorted(REQUIRED_IDS-ids)} "
            f"extra={sorted(ids-REQUIRED_IDS)}"
        )
    row_by_id = {row.get("id"): row for row in rows}
    if (
        row_by_id.get("ehbf_17_subexponential_phi_target", {}).get("readiness")
        != "not_ready_to_apply"
    ):
        issues.append("Xi/Phi determinant-growth target is not open")
    if (
        row_by_id.get("ehbf_09_fractional_residue_guard", {}).get("readiness")
        != "guard_validated"
    ):
        issues.append("fractional-residue guard is not validated")
    for key in (
        "operator_identity",
        "null_flux",
        "integer_flux",
        "determinant_flux",
        "canonical_escape",
        "orthogonal_quotient",
        "orthogonal_convergence",
        "collision_sensor",
        "open_target",
    ):
        if not payload.get("exact", {}).get(key):
            issues.append(f"missing exact statement: {key}")

    for name, checker in (
        ("null_form", null_form_checks),
        ("determinant", determinant_checks),
        ("cauchy_binet", cauchy_binet_checks),
        ("orthogonal_quotient", orthogonal_quotient_checks),
        ("two_atom", two_atom_checks),
    ):
        new_issues, count = checker()
        issues.extend(new_issues)
        counts[name] = count

    text = note.read_text(encoding="utf-8")
    for marker in REQUIRED_NOTE_TEXT:
        if marker not in text:
            issues.append(f"note missing marker: {marker}")
    for forbidden in (
        "therefore RH",
        "proves RH",
        "proves `Lambda <= 0`",
        "the Xi/Phi determinant-growth estimate is proved",
    ):
        if forbidden.lower() in text.lower():
            issues.append(f"forbidden promotion language: {forbidden}")
    return issues, counts


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, counts = validate(args.artifact, args.note)
    for issue in issues:
        print(f"EDREI-HANKEL-FLUX issue: {issue}")
    print(
        "validated Edrei Hankel boundary-flux gate: "
        f"18 rows, {len(issues)} issues, "
        f"{counts.get('null_form', 0)} null-form audits, "
        f"{counts.get('determinant', 0)} determinant audits, "
        f"{counts.get('cauchy_binet', 0)} Cauchy-Binet audits, "
        f"{counts.get('orthogonal_quotient', 0)} orthogonal-quotient audits, "
        f"{counts.get('two_atom', 0)} countermodel/collision audits, "
        "1 open Xi/Phi determinant-growth gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
