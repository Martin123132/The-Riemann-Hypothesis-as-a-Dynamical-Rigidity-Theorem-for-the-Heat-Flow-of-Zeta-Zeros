#!/usr/bin/env python3
"""Validate the planar edge-Gram and Vaughan-symmetrization handoff."""

from __future__ import annotations

import cmath
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.py"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"


def close(left: float, right: float, tolerance: float = 2e-10) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def mobius_values(limit: int) -> list[int]:
    values = [1] * (limit + 1)
    prime = [True] * (limit + 1)
    values[0] = 0
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        for multiple in range(p, limit + 1, p):
            prime[multiple] = False
            values[multiple] *= -1
        square = p * p
        for multiple in range(square, limit + 1, square):
            values[multiple] = 0
    return values


def divisors(value: int) -> list[int]:
    result: list[int] = []
    for divisor in range(1, math.isqrt(value) + 1):
        if value % divisor:
            continue
        result.append(divisor)
        if divisor * divisor != value:
            result.append(value // divisor)
    return sorted(result)


def type_i_brute(
    value: int,
    upper_u: int,
    upper_v: int,
    mu: list[int],
) -> int:
    return sum(
        mu[left] * mu[right]
        for left in range(1, upper_u + 1)
        for right in range(1, upper_v + 1)
        if value % (left * right) == 0
    )


def type_ii_brute(
    value: int,
    upper_u: int,
    upper_v: int,
    mu: list[int],
) -> int:
    return sum(
        mu[left] * mu[right]
        for left in range(upper_u + 1, value + 1)
        for right in range(upper_v + 1, value // left + 1)
        if value % (left * right) == 0
    )


def a_coefficient(
    value: int,
    upper_u: int,
    upper_v: int,
    mu: list[int],
) -> int:
    return sum(
        mu[left] * mu[value // left]
        for left in divisors(value)
        if left <= upper_u and value // left <= upper_v
    )


def b_coefficient(value: int, upper_v: int, mu: list[int]) -> int:
    return sum(mu[divisor] for divisor in divisors(value) if divisor > upper_v)


def type_i_grouped(
    value: int,
    upper_u: int,
    upper_v: int,
    mu: list[int],
) -> int:
    return sum(
        a_coefficient(divisor, upper_u, upper_v, mu)
        for divisor in divisors(value)
        if divisor <= upper_u * upper_v
    )


def type_ii_grouped(
    value: int,
    upper_u: int,
    upper_v: int,
    mu: list[int],
) -> int:
    return sum(
        b_coefficient(divisor, upper_v, mu) * mu[value // divisor]
        for divisor in divisors(value)
        if divisor > upper_v and value // divisor > upper_u
    )


def all_edges(size: int) -> list[tuple[int, int]]:
    return [
        (left, right)
        for left in range(size + 1, 4 * size)
        for right in range(left + 1, 4 * size)
    ]


def threshold_edges(
    edges: list[tuple[int, int]],
    base_cutoff: int,
    shift_cutoff: int,
) -> list[tuple[int, int]]:
    return [
        edge
        for edge in edges
        if edge[0] <= base_cutoff
        and edge[1] - edge[0] <= shift_cutoff
    ]


def directed_pair(
    edges: list[tuple[int, int]],
    left: dict[int, int],
    right: dict[int, int],
) -> int:
    return sum(left[p] * right[q] for p, q in edges)


def symmetric_pair(
    edges: list[tuple[int, int]],
    left: dict[int, int],
    right: dict[int, int],
) -> float:
    return 0.5 * (
        directed_pair(edges, left, right)
        + directed_pair(edges, right, left)
    )


def check_schema(
    payload: dict,
    note_text: str,
    formal_text: str,
    issues: list[str],
) -> None:
    if payload.get("kind") != (
        "jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff"
    ):
        issues.append("result kind mismatch")
    if payload.get("date") != "2026-07-24":
        issues.append("result date mismatch")
    if payload.get("status") != (
        "exact edge-Gram diagonal reduction and Vaughan "
        "symmetrization with one open signed quartic gate"
    ):
        issues.append("result status mismatch")

    rows = payload.get("rows", [])
    summary = payload.get("summary", {})
    if len(rows) != 38 or summary.get("row_count") != 38:
        issues.append("expected 38 ledger rows")
    ids = [item.get("id") for item in rows]
    if len(ids) != len(set(ids)):
        issues.append("duplicate row id")
    if summary.get("exact_reduction_count") != 28:
        issues.append("exact reduction count mismatch")
    if summary.get("proof_guard_count") != 5:
        issues.append("proof guard count mismatch")
    if summary.get("literature_guard_count") != 3:
        issues.append("literature guard count mismatch")
    if summary.get("finite_validation_count") != 1:
        issues.append("finite validation count mismatch")
    if summary.get("open_edge_gate_count") != 1:
        issues.append("open gate count mismatch")
    for key in (
        "edge_gram_proved",
        "diagonal_target_proved",
        "vaughan_symmetrization_proved",
        "unit_cutoff_cancellation_proved",
    ):
        if summary.get(key) is not True:
            issues.append(f"{key} must be true")
    for key in (
        "offdiagonal_edge_gate_proved",
        "curvature_energy_proved",
        "full_burnol_bound_proved",
        "rh_proved",
        "pf_infinity_proved",
        "lambda_le_zero_proved",
    ):
        if summary.get(key) is not False:
            issues.append(f"{key} must be false")

    required_note_markers = (
        "Status: exact edge-Gram and Vaughan-symmetrization reduction",
        "not a proof of RH",
        "The Entire Diagonal Closes",
        "D_(alpha,K)",
        "(27/2)pi^2 K(1+log(2R))",
        "Exact Vaughan Symmetrization",
        "Both cross terms remain",
        "U=V=1",
        "Fourier And Theorem-Fit Audit",
        "Lewko and Lewko",
        "Tao and Teravainen",
        "Matomaki, Radziwill, and Tao",
        "No estimate of (PEGV.19)",
        "Lambda<=0",
    )
    for marker in required_note_markers:
        if marker not in note_text:
            issues.append(f"note missing marker {marker!r}")

    required_formal_markers = (
        "Corollary 11.22Z.6",
        "(11.22Z.100)",
        "(11.22Z.101)",
        "(11.22Z.102)",
        "(11.22Z.103)",
        "(11.22Z.104)",
        "(11.22Z.105)",
        "(11.22Z.106)",
        "(11.22Z.107)",
        "(11.22Z.108)",
        "(11.22Z.109)",
        "jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.md",
    )
    for marker in required_formal_markers:
        if marker not in formal_text:
            issues.append(f"formal core missing marker {marker!r}")


def check_pointwise_vaughan(issues: list[str]) -> int:
    mu = mobius_values(100)
    audits = 0
    for upper_u, upper_v in ((1, 1), (2, 3), (3, 2), (4, 4)):
        for value in range(2, 81):
            type_i = type_i_brute(value, upper_u, upper_v, mu)
            type_ii = type_ii_brute(value, upper_u, upper_v, mu)
            grouped_i = type_i_grouped(value, upper_u, upper_v, mu)
            grouped_ii = type_ii_grouped(value, upper_u, upper_v, mu)
            if (
                value > max(upper_u, upper_v)
                and mu[value] != -type_i + type_ii
            ):
                issues.append(
                    f"pointwise Vaughan mismatch n={value}, "
                    f"U={upper_u}, V={upper_v}"
                )
            if type_i != grouped_i:
                issues.append(
                    f"grouped Type I mismatch n={value}, "
                    f"U={upper_u}, V={upper_v}"
                )
            if type_ii != grouped_ii:
                issues.append(
                    f"grouped Type II mismatch n={value}, "
                    f"U={upper_u}, V={upper_v}"
                )
            audits += 3
    return audits


def check_band_and_gram(issues: list[str]) -> tuple[int, int]:
    mu = mobius_values(40)
    threshold_audits = 0
    gram_audits = 0
    cross_witness = False

    for size in range(2, 7):
        vertices = range(size + 1, 4 * size)
        edges = all_edges(size)
        mu_vector = {value: mu[value] for value in vertices}
        one_vector = {value: 1 for value in vertices}
        type_i = {
            value: type_i_brute(value, 2, 2, mu)
            for value in vertices
        }
        type_ii = {
            value: type_ii_brute(value, 2, 2, mu)
            for value in vertices
        }
        unit_i = {value: 1 for value in vertices}
        unit_ii = {value: 1 + mu[value] for value in vertices}

        thresholds: list[tuple[int, int]] = []
        nu: dict[tuple[int, int], int] = {}
        for base_cutoff in range(size + 1, 4 * size - 1):
            for shift_cutoff in range(1, 3 * size - 1):
                threshold = (base_cutoff, shift_cutoff)
                thresholds.append(threshold)
                nu[threshold] = (
                    1 + (base_cutoff + 2 * shift_cutoff + size) % 7
                )
                active = threshold_edges(
                    edges,
                    base_cutoff,
                    shift_cutoff,
                )
                signed_prefix = directed_pair(
                    active,
                    mu_vector,
                    mu_vector,
                )
                if signed_prefix != symmetric_pair(
                    active,
                    mu_vector,
                    mu_vector,
                ):
                    issues.append(
                        f"band symmetrization mismatch K={size}, "
                        f"x={base_cutoff}, H={shift_cutoff}"
                    )

                four_term = (
                    directed_pair(active, type_i, type_i)
                    + directed_pair(active, type_ii, type_ii)
                    - directed_pair(active, type_i, type_ii)
                    - directed_pair(active, type_ii, type_i)
                )
                three_term = (
                    symmetric_pair(active, type_i, type_i)
                    + symmetric_pair(active, type_ii, type_ii)
                    - 2.0 * symmetric_pair(active, type_i, type_ii)
                )
                if signed_prefix != four_term:
                    issues.append(
                        f"directed four-term mismatch K={size}, "
                        f"x={base_cutoff}, H={shift_cutoff}"
                    )
                if not close(signed_prefix, three_term):
                    issues.append(
                        f"symmetric three-term mismatch K={size}, "
                        f"x={base_cutoff}, H={shift_cutoff}"
                    )

                edge_count = len(active)
                left_linear = directed_pair(
                    active,
                    one_vector,
                    mu_vector,
                )
                right_linear = directed_pair(
                    active,
                    mu_vector,
                    one_vector,
                )
                if left_linear != right_linear:
                    cross_witness = True
                unit_terms = (
                    directed_pair(active, unit_i, unit_i),
                    directed_pair(active, unit_ii, unit_ii),
                    directed_pair(active, unit_i, unit_ii),
                    directed_pair(active, unit_ii, unit_i),
                )
                expected_terms = (
                    edge_count,
                    edge_count
                    + left_linear
                    + right_linear
                    + signed_prefix,
                    edge_count + left_linear,
                    edge_count + right_linear,
                )
                if unit_terms != expected_terms:
                    issues.append(
                        f"unit-cutoff term mismatch K={size}, "
                        f"x={base_cutoff}, H={shift_cutoff}"
                    )
                if (
                    unit_terms[0]
                    + unit_terms[1]
                    - unit_terms[2]
                    - unit_terms[3]
                    != signed_prefix
                ):
                    issues.append(
                        f"unit-cutoff cancellation mismatch K={size}, "
                        f"x={base_cutoff}, H={shift_cutoff}"
                    )

                same_base_direct = 0
                by_base: dict[int, list[int]] = {}
                for left, right in active:
                    by_base.setdefault(left, []).append(right)
                for left, tips in by_base.items():
                    for right in tips:
                        for other in tips:
                            if right != other:
                                same_base_direct += (
                                    mu[left] ** 2
                                    * mu[right]
                                    * mu[other]
                                )
                same_base_formula = sum(
                    mu[left] ** 2
                    * (
                        sum(mu[right] for right in tips) ** 2
                        - sum(mu[right] ** 2 for right in tips)
                    )
                    for left, tips in by_base.items()
                )
                if same_base_direct != same_base_formula:
                    issues.append(
                        f"same-base collision mismatch K={size}, "
                        f"x={base_cutoff}, H={shift_cutoff}"
                    )
                if edge_count > math.comb(3 * size - 1, 2):
                    issues.append(f"edge-count bound failed K={size}")
                threshold_audits += 8

        tail: dict[tuple[int, int], int] = {}
        for base_cutoff in range(size + 1, 4 * size):
            for shift_cutoff in range(1, 3 * size):
                tail[(base_cutoff, shift_cutoff)] = sum(
                    weight
                    for (x_value, h_value), weight in nu.items()
                    if x_value >= base_cutoff
                    and h_value >= shift_cutoff
                )

        edge_labels = {
            edge: mu[edge[0]] * mu[edge[1]]
            for edge in edges
        }
        direct_energy = 0
        direct_diagonal = 0
        for threshold in thresholds:
            active = threshold_edges(edges, *threshold)
            prefix = sum(edge_labels[edge] for edge in active)
            direct_energy += nu[threshold] * prefix * prefix
            direct_diagonal += nu[threshold] * sum(
                edge_labels[edge] ** 2 for edge in active
            )

        gram_energy = 0
        gram_diagonal = 0
        collision = 0
        disjoint = 0
        for edge in edges:
            for other in edges:
                gram_value = tail[
                    (
                        max(edge[0], other[0]),
                        max(
                            edge[1] - edge[0],
                            other[1] - other[0],
                        ),
                    )
                ]
                term = (
                    edge_labels[edge]
                    * edge_labels[other]
                    * gram_value
                )
                gram_energy += term
                if edge == other:
                    gram_diagonal += term
                elif len(set(edge + other)) == 3:
                    collision += term
                elif len(set(edge + other)) == 4:
                    disjoint += term
                else:
                    issues.append(
                        f"unexpected edge-pair stratum K={size}: "
                        f"{edge}, {other}"
                    )
        if direct_energy != gram_energy:
            issues.append(f"edge-Gram energy mismatch K={size}")
        if direct_diagonal != gram_diagonal:
            issues.append(f"edge-Gram diagonal mismatch K={size}")
        if gram_energy != gram_diagonal + collision + disjoint:
            issues.append(f"edge-pair partition mismatch K={size}")
        if gram_diagonal > math.comb(3 * size - 1, 2) * sum(nu.values()):
            issues.append(f"diagonal edge-count bound mismatch K={size}")
        gram_audits += 4

    if not cross_witness:
        issues.append("no directed cross-term noninterchangeability witness")
    return threshold_audits, gram_audits


def check_fourier_identity(issues: list[str]) -> int:
    mu = mobius_values(50)
    audits = 0
    sample_count = 256
    for size in range(3, 9):
        values = list(range(size + 1, 4 * size))
        for shift_cutoff in range(1, 3 * size - 1):
            direct = sum(
                mu[left] * mu[left + shift]
                for shift in range(1, shift_cutoff + 1)
                for left in values
                if left + shift < 4 * size
            )
            quadrature = 0.0
            for index in range(sample_count):
                theta = index / sample_count
                fourier = sum(
                    mu[value] * cmath.exp(2j * math.pi * value * theta)
                    for value in values
                )
                kernel = sum(
                    math.cos(2.0 * math.pi * shift * theta)
                    for shift in range(1, shift_cutoff + 1)
                )
                quadrature += abs(fourier) ** 2 * kernel
            quadrature /= sample_count
            if not close(float(direct), quadrature, tolerance=2e-9):
                issues.append(
                    f"Fourier identity mismatch K={size}, "
                    f"H={shift_cutoff}: {direct} vs {quadrature}"
                )
            audits += 1
    return audits


def check_rebuild(
    payload: dict,
    note_text: str,
    issues: list[str],
) -> None:
    with tempfile.TemporaryDirectory() as directory:
        result = Path(directory) / "result.json"
        note = Path(directory) / "note.md"
        completed = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--result",
                str(result),
                "--note",
                str(note),
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        if completed.returncode:
            issues.append(
                "builder reproduction failed: "
                + (completed.stderr or completed.stdout).strip()
            )
            return
        rebuilt = json.loads(result.read_text(encoding="utf-8"))
        rebuilt_note = note.read_text(encoding="utf-8")
        if rebuilt != payload:
            issues.append("builder reproduction JSON mismatch")
        if rebuilt_note != note_text:
            issues.append("builder reproduction note mismatch")


def main() -> int:
    issues: list[str] = []
    if not RESULT.exists():
        print(f"missing result: {RESULT}")
        return 1
    if not NOTE.exists():
        print(f"missing note: {NOTE}")
        return 1
    if not FORMAL_CORE.exists():
        print(f"missing formal core: {FORMAL_CORE}")
        return 1

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    note_text = NOTE.read_text(encoding="utf-8")
    formal_text = FORMAL_CORE.read_text(encoding="utf-8")
    check_schema(payload, note_text, formal_text, issues)
    vaughan_audits = check_pointwise_vaughan(issues)
    threshold_audits, gram_audits = check_band_and_gram(issues)
    fourier_audits = check_fourier_identity(issues)
    check_rebuild(payload, note_text, issues)

    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        print(
            "validated planar edge-Gram/Vaughan handoff: "
            f"{len(payload.get('rows', []))} rows, {len(issues)} issues"
        )
        return 1

    print(
        "validated planar edge-Gram/Vaughan handoff: "
        f"{len(payload['rows'])} rows, 0 issues, "
        f"{vaughan_audits} Vaughan audits, "
        f"{threshold_audits} band audits, "
        f"{gram_audits} Gram audits, "
        f"{fourier_audits} Fourier audits, "
        "1 open signed edge gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
