#!/usr/bin/env python3
"""Validate the planar incidence and signed anti-diagonal handoff."""

from __future__ import annotations

import cmath
from fractions import Fraction
import itertools
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
    "jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.py"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"


def close(left: float, right: float, tolerance: float = 2e-9) -> bool:
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
        (left, right)
        for left, right in edges
        if left <= base_cutoff and right - left <= shift_cutoff
    ]


def edge_sum(
    edges: list[tuple[int, int]],
    values: dict[int, int],
) -> int:
    return sum(values[left] * values[right] for left, right in edges)


def edge_diagonal(
    edges: list[tuple[int, int]],
    values: dict[int, int],
) -> int:
    return sum(
        values[left] ** 2 * values[right] ** 2
        for left, right in edges
    )


def collision_disjoint(
    edges: list[tuple[int, int]],
    values: dict[int, int],
) -> tuple[int, int]:
    collision = 0
    disjoint = 0
    for edge in edges:
        for other in edges:
            if edge == other:
                continue
            product = (
                values[edge[0]]
                * values[edge[1]]
                * values[other[0]]
                * values[other[1]]
            )
            endpoint_count = len(set(edge + other))
            if endpoint_count == 3:
                collision += product
            elif endpoint_count == 4:
                disjoint += product
            else:
                raise AssertionError("unexpected edge intersection")
    return collision, disjoint


def neighbor_collision(
    vertices: list[int],
    edges: list[tuple[int, int]],
    values: dict[int, int],
) -> int:
    adjacent = {vertex: [] for vertex in vertices}
    for left, right in edges:
        adjacent[left].append(right)
        adjacent[right].append(left)
    return sum(
        values[vertex] ** 2
        * (
            sum(values[other] for other in adjacent[vertex]) ** 2
            - sum(values[other] ** 2 for other in adjacent[vertex])
        )
        for vertex in vertices
    )


def directed_pair(
    edges: list[tuple[int, int]],
    left: dict[int, int],
    right: dict[int, int],
) -> int:
    return sum(left[p] * right[q] for p, q in edges)


def sinc(value: float) -> float:
    if value == 0.0:
        return 1.0
    return math.sin(value) / value


def check_schema(
    payload: dict,
    note_text: str,
    formal_text: str,
    issues: list[str],
) -> None:
    expected_kind = (
        "jensen_window_pf_mertens_planar_"
        "incidence_antidiagonal_handoff"
    )
    if payload.get("kind") != expected_kind:
        issues.append("result kind mismatch")
    if payload.get("date") != "2026-07-24":
        issues.append("result date mismatch")
    if payload.get("status") != (
        "exact endpoint incidence, random-cut bilinear, and signed "
        "current-interior anti-diagonal reductions with one open "
        "projection gate"
    ):
        issues.append("result status mismatch")

    rows = payload.get("rows", [])
    summary = payload.get("summary", {})
    if len(rows) != 43 or summary.get("row_count") != 43:
        issues.append("expected 43 ledger rows")
    row_ids = [item.get("id") for item in rows]
    if len(row_ids) != len(set(row_ids)):
        issues.append("duplicate row id")
    for key, expected in (
        ("exact_reduction_count", 32),
        ("proof_guard_count", 7),
        ("literature_guard_count", 1),
        ("countermodel_count", 1),
        ("finite_validation_count", 1),
        ("open_projection_gate_count", 1),
    ):
        if summary.get(key) != expected:
            issues.append(f"{key} mismatch")
    for key in (
        "endpoint_expansion_proved",
        "perfect_stratum_cancellation_witness_proved",
        "cut_bilinear_equivalence_proved",
        "partial_prefix_fourier_identity_proved",
        "antidiagonal_projection_proved",
    ):
        if summary.get(key) is not True:
            issues.append(f"{key} must be true")
    for key in (
        "projection_estimate_proved",
        "remaining_blocks_proved",
        "curvature_energy_proved",
        "rh_proved",
        "pf_infinity_proved",
        "lambda_le_zero_proved",
    ):
        if summary.get(key) is not False:
            issues.append(f"{key} must be false")

    for source in payload.get("source", {}).values():
        if not (REPO_ROOT / source).exists():
            issues.append(f"missing source {source}")
    if payload.get("literature", {}).get("menon_2026") != (
        "https://arxiv.org/abs/2607.15574"
    ):
        issues.append("missing Menon 2026 literature guard")

    required_note_markers = (
        "Status: exact incidence, cut-bilinear",
        "This is not a proof of RH",
        "Complete Endpoint Incidence Expansion",
        "C_3=210, C_4=-210",
        "Exact Random-Cut Bilinearization",
        "4 E_chi[J_t(chi)^2]=S_t^2+Q_t",
        "P_K(z)=D+z C_3+z^2 C_4",
        "Partial-Prefix Fourier Identity",
        "Signed Current-Interior Anti-Diagonal",
        "sum_(s=0)^(2N-1)kappa_(K,R)(s)^2",
        "sum_(r=1)^R|Y_(K,r)|^2",
        "Recent Short-Interval Literature Guard",
        "arXiv:2607.15574",
        "anchored quadratic edge aggregate",
        "transition,",
        "No",
        "Lambda<=0",
    )
    for marker in required_note_markers:
        if marker not in note_text:
            issues.append(f"note missing marker {marker!r}")

    required_formal_markers = (
        "Corollary 11.22Z.7",
        "(11.22Z.110)",
        "(11.22Z.111)",
        "(11.22Z.112)",
        "(11.22Z.113)",
        "(11.22Z.114)",
        "(11.22Z.115)",
        "(11.22Z.116)",
        "(11.22Z.117)",
        "(11.22Z.118)",
        "(11.22Z.119)",
        "jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.md",
    )
    for marker in required_formal_markers:
        if marker not in formal_text:
            issues.append(f"formal core missing marker {marker!r}")


def check_endpoint_incidence(issues: list[str]) -> int:
    mu = mobius_values(40)
    audits = 0
    for size in range(2, 7):
        vertices = list(range(size + 1, 4 * size))
        edges = all_edges(size)
        vectors = (
            {vertex: mu[vertex] for vertex in vertices},
            {
                vertex: ((vertex * vertex + 3 * size) % 5) - 2
                for vertex in vertices
            },
        )
        thresholds = [
            (base, shift)
            for base in range(size + 1, 4 * size - 1)
            for shift in range(1, 3 * size - 1)
        ]
        weights = {
            threshold: 1 + (threshold[0] + 2 * threshold[1] + size) % 5
            for threshold in thresholds
        }

        def tail(base: int, shift: int) -> int:
            return sum(
                weight
                for (x_value, h_value), weight in weights.items()
                if x_value >= base and h_value >= shift
            )

        for values in vectors:
            integrated_collision = 0
            integrated_disjoint = 0
            for threshold, weight in weights.items():
                active = threshold_edges(edges, *threshold)
                signed = edge_sum(active, values)
                diagonal = edge_diagonal(active, values)
                collision, disjoint = collision_disjoint(active, values)
                if signed * signed - diagonal != collision + disjoint:
                    issues.append(
                        f"offdiagonal split mismatch K={size}, "
                        f"t={threshold}"
                    )
                if collision != neighbor_collision(
                    vertices,
                    active,
                    values,
                ):
                    issues.append(
                        f"neighbor collision mismatch K={size}, "
                        f"t={threshold}"
                    )
                direct_ordered = 0
                for edge in active:
                    for other in active:
                        if edge == other:
                            continue
                        direct_ordered += (
                            values[edge[0]]
                            * values[edge[1]]
                            * values[other[0]]
                            * values[other[1]]
                        )
                if direct_ordered != collision + disjoint:
                    issues.append(
                        f"ordered-pair mismatch K={size}, t={threshold}"
                    )
                integrated_collision += weight * collision
                integrated_disjoint += weight * disjoint
                audits += 3

            triple_formula = 0
            for first, second, third in itertools.combinations(
                vertices,
                3,
            ):
                triple_formula += 2 * (
                    values[first] ** 2
                    * values[second]
                    * values[third]
                    * tail(first, third - first)
                    + values[first]
                    * values[second] ** 2
                    * values[third]
                    * tail(
                        second,
                        max(second - first, third - second),
                    )
                    + values[first]
                    * values[second]
                    * values[third] ** 2
                    * tail(second, third - first)
                )
            quadruple_formula = 0
            for first, second, third, fourth in itertools.combinations(
                vertices,
                4,
            ):
                quadruple_formula += (
                    2
                    * values[first]
                    * values[second]
                    * values[third]
                    * values[fourth]
                    * (
                        tail(
                            third,
                            max(second - first, fourth - third),
                        )
                        + tail(
                            second,
                            max(third - first, fourth - second),
                        )
                        + tail(second, fourth - first)
                    )
                )
            if integrated_collision != triple_formula:
                issues.append(f"sorted-triple mismatch K={size}")
            if integrated_disjoint != quadruple_formula:
                issues.append(f"sorted-quadruple mismatch K={size}")
            audits += 2

    size = 4
    vertices = list(range(5, 16))
    values = dict(
        zip(
            vertices,
            (1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1),
        )
    )
    active = threshold_edges(all_edges(size), 11, 10)
    signed = edge_sum(active, values)
    diagonal = edge_diagonal(active, values)
    collision, disjoint = collision_disjoint(active, values)
    witness = (
        len(active),
        signed,
        diagonal,
        signed * signed - diagonal,
        collision,
        disjoint,
    )
    if witness != (49, 7, 49, 0, 210, -210):
        issues.append(f"cancellation witness mismatch: {witness}")
    audits += 1
    return audits


def check_cut_identities(issues: list[str]) -> int:
    mu = mobius_values(20)
    audits = 0
    for size in (2, 3):
        vertices = list(range(size + 1, 4 * size))
        values = {vertex: mu[vertex] for vertex in vertices}
        edges = all_edges(size)
        for base in range(size + 1, 4 * size - 1):
            for shift in range(1, 3 * size - 1):
                active = threshold_edges(edges, base, shift)
                signed = edge_sum(active, values)
                diagonal = edge_diagonal(active, values)
                cut_squares = []
                for bits in itertools.product((0, 1), repeat=len(vertices)):
                    colors = dict(zip(vertices, bits))
                    cut_sum = sum(
                        values[left] * values[right]
                        for left, right in active
                        if colors[left] != colors[right]
                    )
                    left_vector = {
                        vertex: values[vertex]
                        if colors[vertex] == 0
                        else 0
                        for vertex in vertices
                    }
                    right_vector = {
                        vertex: values[vertex]
                        if colors[vertex] == 1
                        else 0
                        for vertex in vertices
                    }
                    bilinear = directed_pair(
                        active,
                        left_vector,
                        right_vector,
                    ) + directed_pair(
                        active,
                        right_vector,
                        left_vector,
                    )
                    if cut_sum != bilinear:
                        issues.append(
                            f"cut bilinear mismatch K={size}, "
                            f"x={base}, H={shift}"
                        )
                    cut_squares.append(cut_sum * cut_sum)
                    audits += 1
                cut_average = Fraction(
                    sum(cut_squares),
                    len(cut_squares),
                )
                if 4 * cut_average != signed * signed + diagonal:
                    issues.append(
                        f"cut moment mismatch K={size}, "
                        f"x={base}, H={shift}"
                    )
                audits += 1

                if size == 2:
                    collision, disjoint = collision_disjoint(
                        active,
                        values,
                    )
                    mean = Fraction(1, 3)
                    plus_probability = (1 + mean) / 2
                    minus_probability = (1 - mean) / 2
                    biased_moment = Fraction(0)
                    for signs in itertools.product(
                        (-1, 1),
                        repeat=len(vertices),
                    ):
                        probability = Fraction(1)
                        signed_mask = dict(zip(vertices, signs))
                        for sign in signs:
                            probability *= (
                                plus_probability
                                if sign == 1
                                else minus_probability
                            )
                        masked_sum = sum(
                            values[left]
                            * values[right]
                            * signed_mask[left]
                            * signed_mask[right]
                            for left, right in active
                        )
                        biased_moment += (
                            probability * masked_sum * masked_sum
                        )
                    expected = (
                        diagonal
                        + mean**2 * collision
                        + mean**4 * disjoint
                    )
                    if biased_moment != expected:
                        issues.append(
                            f"biased chaos mismatch x={base}, H={shift}"
                        )
                    audits += 1
    return audits


def check_partial_fourier(issues: list[str]) -> int:
    mu = mobius_values(30)
    audits = 0
    sample_count = 257
    roots = [
        cmath.exp(2j * math.pi * sample / sample_count)
        for sample in range(sample_count)
    ]
    for size in range(2, 6):
        vertices = list(range(size + 1, 4 * size))
        values = {vertex: mu[vertex] for vertex in vertices}
        edges = all_edges(size)
        for base in range(size + 1, 4 * size - 1):
            for shift in range(1, 3 * size - 1):
                expected = edge_sum(
                    threshold_edges(edges, base, shift),
                    values,
                )
                integral = 0j
                for root in roots:
                    partial = sum(
                        values[value] * root**value
                        for value in range(size + 1, base + 1)
                    )
                    full = sum(
                        values[value] * root**value
                        for value in vertices
                    )
                    shift_kernel = sum(
                        root**offset
                        for offset in range(1, shift + 1)
                    )
                    integral += (
                        partial * full.conjugate() * shift_kernel
                    )
                integral /= sample_count
                if not close(integral.real, expected, 2e-8):
                    issues.append(
                        f"partial Fourier real mismatch K={size}, "
                        f"x={base}, H={shift}"
                    )
                if abs(integral.imag) > 2e-8:
                    issues.append(
                        f"partial Fourier imaginary residue K={size}, "
                        f"x={base}, H={shift}"
                    )
                audits += 1
    return audits


def kappa_spectral(size: int, rank: int, index: int) -> float:
    normalizer = 2 * size + 1
    return (
        4
        * math.pi**2
        / normalizer**3
        * sum(
            sinc((2 * mode - 1) * math.pi / (2 * normalizer))
            * sinc((2 * mode - 1) * math.pi / normalizer)
            * math.cos(
                (2 * mode - 1)
                * (index + 1.5)
                * math.pi
                / normalizer
            )
            for mode in range(1, rank + 1)
        )
    )


def kappa_integral(size: int, rank: int, index: int) -> float:
    normalizer = 2 * size + 1
    first_width = math.pi / normalizer
    second_width = 2 * math.pi / normalizer
    origin = index * math.pi / normalizer
    total = 0.0
    for mode in range(1, rank + 1):
        frequency = 2 * mode - 1
        integral = (
            math.cos(frequency * (origin + second_width))
            + math.cos(frequency * (origin + first_width))
            - math.cos(
                frequency
                * (origin + first_width + second_width)
            )
            - math.cos(frequency * origin)
        ) / frequency**2
        total += 2 * integral / normalizer
    return total


def check_antidiagonal(issues: list[str]) -> int:
    mu = mobius_values(60)
    audits = 0
    for size in range(5, 13):
        normalizer = 2 * size + 1
        indices = list(range(3, 2 * size - 4))
        values = {
            value: mu[value]
            for value in range(size + 1, 2 * size - 1)
        }
        line_aggregate: dict[int, int] = {}
        for index in indices:
            direct = 0
            for base_offset in range(1, size - 2):
                shift = index - 2 * base_offset
                if shift < 1 or base_offset + shift > size - 2:
                    continue
                active = [
                    (left, right)
                    for left in range(size + 1, size + base_offset + 1)
                    for right in range(left + 1, 2 * size - 1)
                    if right - left <= shift
                ]
                direct += edge_sum(active, values)

            multiplicity_formula = 0
            for edge_base in range(1, size - 2):
                for gap in range(1, size - 1 - edge_base):
                    lower = max(
                        edge_base,
                        index - size + 2,
                    )
                    upper = (index - gap) // 2
                    multiplicity = max(0, upper - lower + 1)
                    multiplicity_formula += (
                        values[size + edge_base]
                        * values[size + edge_base + gap]
                        * multiplicity
                    )
            if direct != multiplicity_formula:
                issues.append(
                    f"edge multiplicity mismatch K={size}, s={index}"
                )
            line_aggregate[index] = direct
            audits += 1

        for rank in range(1, 5):
            for index in indices:
                spectral = kappa_spectral(size, rank, index)
                integral = kappa_integral(size, rank, index)
                if not close(spectral, integral, 2e-11):
                    issues.append(
                        f"kappa formula mismatch K={size}, "
                        f"R={rank}, s={index}"
                    )
                audits += 1

            full_square = sum(
                kappa_spectral(size, rank, index) ** 2
                for index in range(2 * normalizer)
            )
            sigma_square = sum(
                (
                    sinc(
                        (2 * mode - 1)
                        * math.pi
                        / (2 * normalizer)
                    )
                    * sinc(
                        (2 * mode - 1)
                        * math.pi
                        / normalizer
                    )
                )
                ** 2
                for mode in range(1, rank + 1)
            )
            expected_square = (
                16
                * math.pi**4
                / normalizer**5
                * sigma_square
            )
            if not close(full_square, expected_square, 2e-10):
                issues.append(
                    f"kappa Parseval mismatch K={size}, R={rank}"
                )

            direct_pairing = sum(
                kappa_spectral(size, rank, index)
                * line_aggregate[index]
                for index in indices
            )
            projections = []
            for mode in range(1, rank + 1):
                frequency = 2 * mode - 1
                projections.append(
                    sum(
                        line_aggregate[index]
                        * math.cos(
                            frequency
                            * (index + 1.5)
                            * math.pi
                            / normalizer
                        )
                        for index in indices
                    )
                )
            spectral_pairing = (
                4
                * math.pi**2
                / normalizer**3
                * sum(
                    sinc(
                        (2 * mode - 1)
                        * math.pi
                        / (2 * normalizer)
                    )
                    * sinc(
                        (2 * mode - 1)
                        * math.pi
                        / normalizer
                    )
                    * projections[mode - 1]
                    for mode in range(1, rank + 1)
                )
            )
            if not close(direct_pairing, spectral_pairing, 2e-10):
                issues.append(
                    f"spectral projection mismatch K={size}, R={rank}"
                )

            cauchy_bound = (
                16
                * math.pi**4
                * rank
                / normalizer**6
                * sum(value * value for value in projections)
            )
            if direct_pairing * direct_pairing > (
                cauchy_bound + 2e-10
            ):
                issues.append(
                    f"projection Cauchy failed K={size}, R={rank}"
                )

            bessel_left = sum(value * value for value in projections)
            bessel_right = normalizer * sum(
                value * value for value in line_aggregate.values()
            )
            if bessel_left > bessel_right + 2e-9:
                issues.append(
                    f"odd-frequency Bessel failed K={size}, R={rank}"
                )
            audits += 4
    return audits


def check_reproduction(issues: list[str]) -> None:
    with tempfile.TemporaryDirectory() as temp_directory:
        temp_root = Path(temp_directory)
        temp_result = temp_root / "result.json"
        temp_note = temp_root / "note.md"
        process = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--output",
                str(temp_result),
                "--note",
                str(temp_note),
            ],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        if process.returncode:
            issues.append(
                "builder reproduction failed: "
                + (process.stderr.strip() or process.stdout.strip())
            )
            return
        if temp_result.read_bytes() != RESULT.read_bytes():
            issues.append("builder result reproduction mismatch")
        if temp_note.read_bytes() != NOTE.read_bytes():
            issues.append("builder note reproduction mismatch")


def main() -> int:
    issues: list[str] = []
    try:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        note_text = NOTE.read_text(encoding="utf-8")
        formal_text = FORMAL_CORE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as error:
        print(f"PEIA [load-failure] {error}")
        return 1

    check_schema(payload, note_text, formal_text, issues)
    incidence_audits = check_endpoint_incidence(issues)
    cut_audits = check_cut_identities(issues)
    fourier_audits = check_partial_fourier(issues)
    antidiagonal_audits = check_antidiagonal(issues)
    check_reproduction(issues)

    for issue in issues:
        print(f"PEIA [issue] {issue}")
    print(
        "validated planar incidence/anti-diagonal handoff: "
        f"{len(payload.get('rows', []))} rows, "
        f"{len(issues)} issues, "
        f"{incidence_audits} incidence audits, "
        f"{cut_audits} cut audits, "
        f"{fourier_audits} Fourier audits, "
        f"{antidiagonal_audits} anti-diagonal audits, "
        f"{payload.get('summary', {}).get('open_projection_gate_count')} "
        "open projection gate"
    )
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
