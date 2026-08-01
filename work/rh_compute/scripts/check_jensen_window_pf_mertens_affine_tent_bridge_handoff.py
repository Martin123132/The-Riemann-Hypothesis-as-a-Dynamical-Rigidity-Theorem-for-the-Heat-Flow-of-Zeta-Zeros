#!/usr/bin/env python3
"""Check the affine-tent and Brownian-bridge unified Vaughan handoff."""

from __future__ import annotations

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
    "jensen_window_pf_mertens_affine_tent_bridge_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_affine_tent_bridge_handoff.md"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_affine_tent_bridge_handoff.py"
)

EXPECTED_KIND = "jensen_window_pf_mertens_affine_tent_bridge_handoff"
EXPECTED_STATUS = (
    "exact affine scale-localization and unified local Vaughan "
    "reduction with one open signed gate"
)
EXPECTED_ROW_COUNT = 53
EXPECTED_EXACT_COUNT = 44
EXPECTED_GUARD_COUNT = 4
EXPECTED_OPEN_COUNT = 1

REQUIRED_NOTE_MARKERS = (
    "u_K-(1/2)u_(2K)",
    "rho_alpha:=2^(-(1+alpha)/2)",
    "W_K=w_K-(1/2)w_(2K)",
    "(1-rho_alpha)||d||_2",
    "M_alpha<infinity",
    "0<=D_(aff,K)<=3K^(-alpha)",
    "0<=G_(K,h)(n)<=(5/4)K^(-1-alpha)",
    "D_(loc,K)",
    "(19/6)K^(-alpha)",
    "R_alpha<infinity",
    "O_(loc,K)=-TI_(loc,K)+TII_(loc,K)",
    "K^(alpha-1)i*j",
    "Z_K/K+sum_(i=t)^(K-1)q_(K+i)",
    "(MATBH.34)",
    "beta_K",
    "(1/(4sqrt(2)))*K^(-s)||B_K||_2",
    "sum_(K dyadic)K^(-2-alpha)",
    "(MATBH.39)",
    "H_K=O(K^3)",
    "theta>0.55",
    "exceptional sets of base points",
    "fails at `alpha=0`",
    "-TI_(loc,K)+TII_(loc,K)",
    "one full power",
    "all remain open",
)

REQUIRED_FORMAL_MARKERS = (
    "Lemma 11.22Y: Mertens Affine-Tent/Bridge Unification",
    "v_K:=u_K-(1/2)u_(2K)",
    "rho_alpha:=2^(-(1+alpha)/2)",
    "0<=G_(K,h)(n)<=(5/4)K^(-1-alpha)",
    "O_(loc,K)=-TI_(loc,K)+TII_(loc,K)",
    "K^(alpha-1)i*j",
    "Z_K/K+sum_(i=t)^(K-1)q_(K+i)",
    "beta_K",
    "(1/(4sqrt(2)))*K^(-s)||B_K||_2",
    "sum_(K dyadic)K^(-2-alpha)",
    "(11.22Y.31)",
)


def close(left: float, right: float, tolerance: float = 1e-11) -> bool:
    return abs(left - right) <= tolerance * max(
        1.0,
        abs(left),
        abs(right),
    )


def mobius_values(limit: int) -> list[int]:
    mu = [0] * (limit + 1)
    mu[1] = 1
    primes: list[int] = []
    composite = [False] * (limit + 1)
    for n in range(2, limit + 1):
        if not composite[n]:
            primes.append(n)
            mu[n] = -1
        for prime in primes:
            value = n * prime
            if value > limit:
                break
            composite[value] = True
            if n % prime == 0:
                mu[value] = 0
                break
            mu[value] = -mu[n]
    return mu


def tent(size: int, n: int) -> float:
    if size < n <= 2 * size:
        return float(n - size)
    if 2 * size < n < 4 * size:
        return (4.0 * size - n) / 2.0
    return 0.0


def combined_kernel(
    alpha: float,
    size: int,
    shift: int,
    n: int,
) -> float:
    right = n + shift
    value = 0.0
    if size < n < right < 4 * size:
        value += (
            size ** (alpha - 1.0)
            * (n * right) ** (-1.0 - alpha)
            * tent(size, n)
            * tent(size, right)
        )
    if size < n < right < 2 * size:
        value += (
            size ** (-3.0 - alpha)
            * (n - size)
            * (2 * size - n - shift)
        )
    return value


def check_schema(payload: dict, note_text: str, issues: list[str]) -> None:
    if payload.get("kind") != EXPECTED_KIND:
        issues.append("unexpected result kind")
    if payload.get("date") != "2026-07-23":
        issues.append("unexpected result date")
    if payload.get("status") != EXPECTED_STATUS:
        issues.append("unexpected result status")

    rows = payload.get("rows")
    if not isinstance(rows, list):
        issues.append("rows must be a list")
        return
    if len(rows) != EXPECTED_ROW_COUNT:
        issues.append(
            f"expected {EXPECTED_ROW_COUNT} rows, found {len(rows)}"
        )
    row_ids = [entry.get("id") for entry in rows if isinstance(entry, dict)]
    if len(row_ids) != len(set(row_ids)):
        issues.append("row ids are not unique")
    for index, entry in enumerate(rows):
        if not isinstance(entry, dict):
            issues.append(f"row {index} is not an object")
            continue
        for field in ("id", "role", "status", "statement", "proof_boundary"):
            value = entry.get(field)
            if not isinstance(value, str) or not value.strip():
                issues.append(f"row {index} has invalid {field}")

    audit = payload.get("audit", {})
    expected_scalars = {
        "row_count": EXPECTED_ROW_COUNT,
        "exact_reduction_count": EXPECTED_EXACT_COUNT,
        "literature_guard_count": 4,
        "proof_guard_count": EXPECTED_GUARD_COUNT,
        "open_signed_gate_count": EXPECTED_OPEN_COUNT,
    }
    for key, expected in expected_scalars.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key} expected {expected}, found {audit.get(key)}"
            )
    expected_true = (
        "affine_infinite_tail_localized",
        "scale_filter_invertible",
        "affine_bridge_kernel_unified",
        "combined_diagonal_summable",
        "two_interval_vaughan_handoff_proved",
        "positive_tail_lattice_square_proved",
        "ordinary_mobius_anchored_criterion_proved",
    )
    for key in expected_true:
        if audit.get(key) is not True:
            issues.append(f"audit flag {key} must be true")
    expected_false = (
        "averaged_chowla_closes_gate",
        "signed_type_i_ii_gain_proved",
        "full_burnol_bound_proved",
        "rh_proved",
        "pf_infinity_proved",
        "lambda_le_zero_proved",
        "alpha_zero_passage_allowed",
    )
    for key in expected_false:
        if audit.get(key) is not False:
            issues.append(f"audit flag {key} must be false")

    for marker in REQUIRED_NOTE_MARKERS:
        if marker not in note_text:
            issues.append(f"note missing marker: {marker}")


def check_builder_reproduction(
    payload: dict,
    note_text: str,
    issues: list[str],
) -> None:
    with tempfile.TemporaryDirectory(prefix="matbh_check_") as temp_dir:
        result_path = Path(temp_dir) / "result.json"
        note_path = Path(temp_dir) / "note.md"
        process = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--result",
                str(result_path),
                "--note",
                str(note_path),
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        if process.returncode != 0:
            issues.append(
                "builder reproduction failed: "
                + (process.stderr.strip() or process.stdout.strip())
            )
            return
        rebuilt = json.loads(result_path.read_text(encoding="utf-8"))
        rebuilt_note = note_path.read_text(encoding="utf-8")
        if rebuilt != payload:
            issues.append("builder reproduction changed the JSON payload")
        if rebuilt_note != note_text:
            issues.append("builder reproduction changed the note")


def check_arbitrary_tail_filter(issues: list[str]) -> None:
    alpha_values = (0.09, 0.37, 0.79)
    support = 1024
    coefficients = [0.0] * (support + 1)
    for n in range(1, support + 1):
        sign = (-1.0, 0.0, 1.0)[(17 * n + n * n) % 3]
        coefficients[n] = sign * n ** (-1.0)

    for alpha in alpha_values:
        q = [0.0] + [
            coefficients[n] * n ** (-alpha)
            for n in range(1, support + 1)
        ]
        tails = [0.0] * (support + 1)
        for n in range(support - 1, -1, -1):
            tails[n] = tails[n + 1] + q[n + 1]

        dyadic_sizes = []
        size = 1
        while size <= support:
            dyadic_sizes.append(size)
            size *= 2
        d_values: list[float] = []
        e_values: list[float] = []
        p = (alpha - 1.0) / 2.0
        rho = 2.0 ** (-(1.0 + alpha) / 2.0)

        for size in dyadic_sizes:
            if 4 * size <= support:
                u_value = sum(tails[n] for n in range(size, 2 * size))
                kernel_value = sum(
                    q[n] * min(max(n - size, 0), size)
                    for n in range(1, support + 1)
                )
                if not close(u_value, kernel_value, tolerance=5e-12):
                    issues.append(
                        f"tail-count kernel failed at alpha={alpha}, "
                        f"K={size}"
                    )

                u_double = sum(
                    tails[n] for n in range(2 * size, 4 * size)
                )
                filtered = u_value - 0.5 * u_double
                tent_value = sum(
                    q[n] * tent(size, n)
                    for n in range(1, support + 1)
                )
                if not close(filtered, tent_value, tolerance=5e-12):
                    issues.append(
                        f"tent localization failed at alpha={alpha}, "
                        f"K={size}"
                    )

            u_value = sum(
                tails[n]
                for n in range(size, min(2 * size, support + 1))
            )
            u_double = sum(
                tails[n]
                for n in range(
                    2 * size,
                    min(4 * size, support + 1),
                )
            )
            d_values.append(size**p * u_value)
            e_values.append(size**p * (u_value - 0.5 * u_double))

        for index in range(len(d_values) - 1):
            if not close(
                e_values[index],
                d_values[index] - rho * d_values[index + 1],
                tolerance=5e-12,
            ):
                issues.append(
                    f"normalized scale filter failed at alpha={alpha}, "
                    f"index={index}"
                )

        d_norm = math.sqrt(sum(value * value for value in d_values))
        e_norm = math.sqrt(sum(value * value for value in e_values))
        if e_norm > (1.0 + rho) * d_norm + 2e-11:
            issues.append(f"filter upper norm bound failed at alpha={alpha}")
        if e_norm + 2e-11 < (1.0 - rho) * d_norm:
            issues.append(f"filter lower norm bound failed at alpha={alpha}")


def check_signed_expansions(issues: list[str]) -> None:
    mu = mobius_values(2048)
    for alpha in (0.13, 0.47, 0.83):
        for size in (4, 8, 16, 31, 64, 96):
            e_value = (
                size ** ((alpha - 1.0) / 2.0)
                * sum(
                    mu[n] * n ** (-1.0 - alpha) * tent(size, n)
                    for n in range(size + 1, 4 * size)
                )
            )
            affine_diagonal = (
                size ** (alpha - 1.0)
                * sum(
                    mu[n] ** 2
                    * n ** (-2.0 - 2.0 * alpha)
                    * tent(size, n) ** 2
                    for n in range(size + 1, 4 * size)
                )
            )
            affine_offdiagonal = (
                size ** (alpha - 1.0)
                * sum(
                    mu[n]
                    * mu[right]
                    * (n * right) ** (-1.0 - alpha)
                    * tent(size, n)
                    * tent(size, right)
                    for n in range(size + 1, 4 * size)
                    for right in range(n + 1, 4 * size)
                )
            )
            if not close(
                e_value * e_value,
                affine_diagonal + 2.0 * affine_offdiagonal,
                tolerance=2e-11,
            ):
                issues.append(
                    f"affine signed expansion failed at alpha={alpha}, "
                    f"K={size}"
                )
            if affine_diagonal > 3.0 * size ** (-alpha) * (1.0 + 1e-12):
                issues.append(
                    f"affine diagonal bound failed at alpha={alpha}, "
                    f"K={size}"
                )

            path = [0.0]
            for index in range(1, size):
                path.append(path[-1] + mu[size + index])
            mean = sum(path) / size
            bridge_energy = (
                size ** (-2.0 - alpha)
                * sum((value - mean) ** 2 for value in path)
            )
            bridge_diagonal = (
                size ** (-2.0 - alpha)
                * sum(
                    mu[size + index] ** 2
                    * index
                    * (size - index)
                    / size
                    for index in range(1, size)
                )
            )
            bridge_offdiagonal = (
                size ** (-2.0 - alpha)
                * sum(
                    mu[size + index]
                    * mu[size + index + shift]
                    * index
                    * (size - index - shift)
                    / size
                    for shift in range(1, size - 1)
                    for index in range(1, size - shift)
                )
            )
            if not close(
                bridge_energy,
                bridge_diagonal + 2.0 * bridge_offdiagonal,
                tolerance=2e-11,
            ):
                issues.append(
                    f"bridge signed expansion failed at alpha={alpha}, "
                    f"K={size}"
                )

            weighted_path = [0.0]
            for index in range(1, size):
                weighted_path.append(
                    weighted_path[-1]
                    + mu[size + index]
                    * (size + index) ** (-1.0 - alpha)
                )
            weighted_mean = sum(weighted_path) / size
            weighted_centered = sum(
                (value - weighted_mean) ** 2
                for value in weighted_path
            )
            future_anchor = (
                size
                * mu[2 * size]
                * (2 * size) ** (-1.0 - alpha)
                + sum(
                    (4 * size - n)
                    / 2.0
                    * mu[n]
                    * n ** (-1.0 - alpha)
                    for n in range(2 * size + 1, 4 * size)
                )
            )
            z_value = future_anchor / size
            s_value = 1.0 + alpha
            beta_value = (2.0 * size) ** s_value * z_value
            beta_explicit = mu[2 * size] + sum(
                ((2.0 * size) / n) ** s_value
                * (4 * size - n)
                / (2.0 * size)
                * mu[n]
                for n in range(2 * size + 1, 4 * size)
            )
            if not close(beta_value, beta_explicit, tolerance=3e-11):
                issues.append(
                    f"future anchor normalization failed at "
                    f"alpha={alpha}, K={size}"
                )

            weighted_suffixes = [
                z_value
                + sum(
                    mu[size + index]
                    * (size + index) ** (-s_value)
                    for index in range(t, size)
                )
                for t in range(1, size + 1)
            ]
            ordinary_suffixes = [
                beta_value
                + sum(
                    mu[size + index]
                    for index in range(t, size)
                )
                for t in range(1, size + 1)
            ]
            weighted_norm = math.sqrt(
                sum(value * value for value in weighted_suffixes)
            )
            ordinary_norm = math.sqrt(
                sum(value * value for value in ordinary_suffixes)
            )
            if (
                weighted_norm
                > math.sqrt(6.0)
                * size ** (-s_value)
                * ordinary_norm
                + 3e-11
            ):
                issues.append(
                    f"ordinary suffix upper transfer failed at "
                    f"alpha={alpha}, K={size}"
                )
            if (
                weighted_norm + 3e-11
                < size ** (-s_value)
                * ordinary_norm
                / (4.0 * math.sqrt(2.0))
            ):
                issues.append(
                    f"ordinary suffix lower transfer failed at "
                    f"alpha={alpha}, K={size}"
                )

            tail_lattice = size**alpha * sum(
                (
                    z_value
                    + sum(
                        mu[size + index]
                        * (size + index) ** (-1.0 - alpha)
                        for index in range(t, size)
                    )
                )
                ** 2
                for t in range(1, size + 1)
            )
            if not close(
                e_value * e_value + size**alpha * weighted_centered,
                tail_lattice,
                tolerance=3e-11,
            ):
                issues.append(
                    f"positive tail-lattice square failed at "
                    f"alpha={alpha}, K={size}"
                )
            for left in range(1, size):
                for right in range(1, size):
                    completed = (
                        size ** (alpha - 1.0) * left * right
                        + size**alpha
                        * (
                            min(left, right)
                            - left * right / size
                        )
                    )
                    target = size**alpha * min(left, right)
                    if not close(completed, target, tolerance=2e-12):
                        issues.append(
                            f"constant-direction completion failed at "
                            f"alpha={alpha}, K={size}"
                        )
                        break

            combined_offdiagonal = sum(
                mu[n]
                * mu[n + shift]
                * combined_kernel(alpha, size, shift, n)
                for shift in range(1, 3 * size)
                for n in range(size + 1, 4 * size - shift)
            )
            if not close(
                combined_offdiagonal,
                affine_offdiagonal + bridge_offdiagonal,
                tolerance=3e-11,
            ):
                issues.append(
                    f"combined kernel failed at alpha={alpha}, K={size}"
                )
            for shift in range(1, 3 * size):
                for n in range(size + 1, 4 * size - shift):
                    kernel = combined_kernel(alpha, size, shift, n)
                    if kernel < -1e-16:
                        issues.append(
                            f"negative combined kernel at alpha={alpha}, "
                            f"K={size}"
                        )
                    if (
                        kernel
                        > 1.25
                        * size ** (-1.0 - alpha)
                        * (1.0 + 2e-12)
                    ):
                        issues.append(
                            f"combined amplitude failed at alpha={alpha}, "
                            f"K={size}"
                    )


def check_full_suffix_abel_transfer(issues: list[str]) -> None:
    for alpha in (0.01, 0.37, 0.99):
        s_value = 1.0 + alpha
        for size in (2, 3, 8, 17, 64, 129):
            weights = [0.0] + [
                (size + index) ** (-s_value)
                for index in range(1, size)
            ] + [(2 * size) ** (-s_value)]
            values = [0.0] + [
                ((17 * index + 3 * index * index) % 19 - 9) / 7.0
                for index in range(1, size + 1)
            ]
            weighted = [0.0] * (size + 1)
            ordinary = [0.0] * (size + 1)
            for t in range(size, 0, -1):
                weighted[t] = (
                    weighted[t + 1] if t < size else 0.0
                ) + weights[t] * values[t]
                ordinary[t] = (
                    ordinary[t + 1] if t < size else 0.0
                ) + values[t]

            for t in range(1, size + 1):
                forward = weights[t] * ordinary[t] + sum(
                    (weights[index] - weights[index - 1])
                    * ordinary[index]
                    for index in range(t + 1, size + 1)
                )
                inverse = weighted[t] / weights[t] + sum(
                    (
                        1.0 / weights[index]
                        - 1.0 / weights[index - 1]
                    )
                    * weighted[index]
                    for index in range(t + 1, size + 1)
                )
                if not close(forward, weighted[t], tolerance=3e-11):
                    issues.append(
                        f"forward suffix Abel identity failed at "
                        f"alpha={alpha}, K={size}, t={t}"
                    )
                    break
                if not close(inverse, ordinary[t], tolerance=3e-11):
                    issues.append(
                        f"inverse suffix Abel identity failed at "
                        f"alpha={alpha}, K={size}, t={t}"
                    )
                    break

            weighted_norm = math.sqrt(
                sum(value * value for value in weighted[1:])
            )
            ordinary_norm = math.sqrt(
                sum(value * value for value in ordinary[1:])
            )
            if (
                weighted_norm
                > math.sqrt(6.0)
                * size ** (-s_value)
                * ordinary_norm
                + 3e-11
            ):
                issues.append(
                    f"arbitrary suffix upper norm failed at "
                    f"alpha={alpha}, K={size}"
                )
            if (
                weighted_norm + 3e-11
                < size ** (-s_value)
                * ordinary_norm
                / (4.0 * math.sqrt(2.0))
            ):
                issues.append(
                    f"arbitrary suffix lower norm failed at "
                    f"alpha={alpha}, K={size}"
                )


def vaughan_interval(
    mu: list[int],
    base: int,
    test,
) -> tuple[float, float]:
    u_cut = max(1, math.isqrt(base))
    v_cut = u_cut
    a_coefficients = {
        d: sum(
            mu[b] * mu[c]
            for b in range(1, u_cut + 1)
            for c in range(1, v_cut + 1)
            if b * c == d
        )
        for d in range(1, u_cut * v_cut + 1)
    }
    b_coefficients = {
        d: sum(
            mu[c]
            for c in range(v_cut + 1, d + 1)
            if d % c == 0
        )
        for d in range(v_cut + 1, (2 * base) // u_cut + 1)
    }
    type_i = 0.0
    for d, coefficient in a_coefficients.items():
        type_i += coefficient * sum(
            test(d * w)
            for w in range(base // d + 1, (2 * base) // d + 1)
        )
    type_ii = 0.0
    for d, coefficient in b_coefficients.items():
        type_ii += coefficient * sum(
            mu[w] * test(d * w)
            for w in range(
                max(u_cut, base // d) + 1,
                (2 * base) // d + 1,
            )
        )
    return type_i, type_ii


def check_two_interval_vaughan(issues: list[str]) -> None:
    mu = mobius_values(4096)
    for alpha in (0.19, 0.61):
        for size in (4, 8, 16, 24, 32):
            direct = 0.0
            type_i = 0.0
            type_ii = 0.0
            for shift in range(1, 3 * size):
                for interval_base in (size, 2 * size):
                    def test(n: int) -> float:
                        return (
                            mu[n + shift]
                            * combined_kernel(alpha, size, shift, n)
                        )

                    direct += sum(
                        mu[n] * test(n)
                        for n in range(
                            interval_base + 1,
                            2 * interval_base + 1,
                        )
                    )
                    interval_i, interval_ii = vaughan_interval(
                        mu,
                        interval_base,
                        test,
                    )
                    type_i += interval_i
                    type_ii += interval_ii
            if not close(
                direct,
                -type_i + type_ii,
                tolerance=3e-10,
            ):
                issues.append(
                    f"two-interval Vaughan failed at alpha={alpha}, "
                    f"K={size}: direct={direct}, "
                    f"rhs={-type_i + type_ii}"
                )


def check_formal_core(issues: list[str]) -> None:
    if not FORMAL_CORE.exists():
        issues.append("formal core is missing")
        return
    text = FORMAL_CORE.read_text(encoding="utf-8")
    for marker in REQUIRED_FORMAL_MARKERS:
        if marker not in text:
            issues.append(f"formal core missing marker: {marker}")


def main() -> int:
    issues: list[str] = []
    if not RESULT.exists():
        issues.append(f"missing result: {RESULT}")
    if not NOTE.exists():
        issues.append(f"missing note: {NOTE}")
    if not BUILDER.exists():
        issues.append(f"missing builder: {BUILDER}")
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    note_text = NOTE.read_text(encoding="utf-8")
    check_schema(payload, note_text, issues)
    check_builder_reproduction(payload, note_text, issues)
    check_arbitrary_tail_filter(issues)
    check_signed_expansions(issues)
    check_full_suffix_abel_transfer(issues)
    check_two_interval_vaughan(issues)
    check_formal_core(issues)

    audit = payload.get("audit", {})
    print(
        "validated Mertens affine-tent/bridge handoff: "
        f"{len(payload.get('rows', []))} rows, "
        f"{len(issues)} issues, "
        f"{audit.get('exact_reduction_count')} exact reductions, "
        f"{audit.get('proof_guard_count')} proof guards, "
        f"{audit.get('open_signed_gate_count')} open gate"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
