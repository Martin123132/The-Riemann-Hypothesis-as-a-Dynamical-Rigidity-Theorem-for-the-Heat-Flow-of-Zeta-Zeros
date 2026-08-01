#!/usr/bin/env python3
"""Check the Mertens shift-kernel variation and joint-cancellation handoff."""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_shift_kernel_variation_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_shift_kernel_variation_handoff.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_shift_kernel_variation_handoff.py"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"
KIND = "jensen_window_pf_mertens_shift_kernel_variation_handoff"
STATUS = (
    "exact shift decomposition, two-axis bounded variation, "
    "and joint-cancellation handoff with one open arithmetic gate"
)
REQUIRED_NOTE_MARKERS = (
    "# Jensen-Window PF Mertens Shift-Kernel Variation Handoff",
    "O_(alpha,K)=sum_(h=1)^(3K-2)Q_(K,h)",
    "Var_(n in I_(K,h))W_(K,h)(n)<6rho_K",
    "Var_h W_(K,h)(n)<pi^2/N+rho_K",
    "sum_(K<n<m<4K)G(n,m)<pi^2*K/2",
    "Why Axiswise Cancellation Is Insufficient",
    "<u_I,G u_I>+<u_II,G u_II>-2<u_I,G u_II>",
    "The weighted joint Mobius gain",
)
SPECIAL_ROLES = {
    "proof_guard",
    "literature_guard",
    "conditional_calibration",
    "finite_validation",
    "open_gate",
}


def close(left: float, right: float, tolerance: float = 1e-11) -> bool:
    return abs(left - right) <= tolerance * max(
        1.0, abs(left), abs(right)
    )


def mobius_values(limit: int) -> list[int]:
    values = [1] * (limit + 1)
    prime = [True] * (limit + 1)
    primes: list[int] = []
    values[0] = 0
    for n in range(2, limit + 1):
        if prime[n]:
            primes.append(n)
            values[n] = -1
        for p in primes:
            if n * p > limit:
                break
            prime[n * p] = False
            if n % p == 0:
                values[n * p] = 0
                break
            values[n * p] = -values[n]
    return values


def angle(size: int, mode: int) -> float:
    return (2 * mode - 1) * math.pi / (2 * size + 1)


def eigenvector_value(size: int, mode: int, index: int) -> float:
    return (
        2.0
        / math.sqrt(2 * size + 1)
        * math.sin(index * angle(size, mode))
    )


def kernel(
    size: int,
    cutoff: int,
    left: int,
    right: int,
) -> float:
    return sum(
        eigenvector_value(size, mode, left)
        * eigenvector_value(size, mode, right)
        / (2 * mode - 1) ** 2
        for mode in range(1, cutoff + 1)
    )


def future_weight(alpha: float, size: int, n: int) -> float:
    if not 2 * size <= n < 4 * size:
        raise ValueError("future weight outside support")
    return (
        (2.0 * size / n) ** (1.0 + alpha)
        * (4.0 * size - n)
        / (2.0 * size)
    )


def position_weight(
    alpha: float,
    size: int,
    n: int,
) -> tuple[int, float]:
    if size < n < 2 * size:
        return n - size, 1.0
    if 2 * size <= n < 4 * size:
        return size, future_weight(alpha, size, n)
    raise ValueError("ordinary position outside support")


def ordinary_kernel(
    alpha: float,
    size: int,
    cutoff: int,
    left: int,
    right: int,
) -> float:
    left_position, left_weight = position_weight(alpha, size, left)
    right_position, right_weight = position_weight(alpha, size, right)
    return (
        left_weight
        * right_weight
        * kernel(size, cutoff, left_position, right_position)
    )


def interval_values(lower: int, upper: int) -> list[int]:
    if lower > upper:
        return []
    return list(range(lower, upper + 1))


def block_supports(
    size: int,
    shift: int,
) -> tuple[list[int], list[int], list[int]]:
    current_current = interval_values(
        size + 1,
        2 * size - 1 - shift,
    )
    current_future = interval_values(
        max(size + 1, 2 * size - shift),
        min(2 * size - 1, 4 * size - 1 - shift),
    )
    future_future = interval_values(
        2 * size,
        4 * size - 1 - shift,
    )
    return current_current, current_future, future_future


def variation(values: list[float]) -> float:
    return sum(
        abs(right - left)
        for left, right in zip(values, values[1:])
    )


def check_schema(payload: dict, note_text: str, issues: list[str]) -> None:
    if payload.get("kind") != KIND:
        issues.append("result kind mismatch")
    if payload.get("date") != "2026-07-24":
        issues.append("result date mismatch")
    if payload.get("status") != STATUS:
        issues.append("result status mismatch")
    rows = payload.get("rows")
    if not isinstance(rows, list) or len(rows) != 45:
        issues.append("expected exactly 45 rows")
        return
    row_ids = [item.get("id") for item in rows]
    expected_ids = [f"skv_{index:02d}_" for index in range(1, 46)]
    if len(set(row_ids)) != 45:
        issues.append("row ids are not unique")
    for prefix in expected_ids:
        if not any(
            isinstance(row_id, str) and row_id.startswith(prefix)
            for row_id in row_ids
        ):
            issues.append(f"missing row prefix {prefix}")
    exact_count = sum(
        item.get("role") not in SPECIAL_ROLES for item in rows
    )
    if exact_count != 39:
        issues.append(f"expected 39 exact rows, found {exact_count}")
    role_counts = {
        role: sum(item.get("role") == role for item in rows)
        for role in SPECIAL_ROLES
    }
    expected_role_counts = {
        "proof_guard": 2,
        "literature_guard": 1,
        "conditional_calibration": 1,
        "finite_validation": 1,
        "open_gate": 1,
    }
    if role_counts != expected_role_counts:
        issues.append(
            f"special role counts mismatch: {role_counts}"
        )
    audit = payload.get("audit", {})
    expected_audit = {
        "row_count": 45,
        "exact_reduction_count": 39,
        "literature_guard_count": 1,
        "proof_guard_count": 2,
        "conditional_calibration_count": 1,
        "finite_validation_count": 1,
        "open_signed_gate_count": 1,
        "shift_decomposition_proved": True,
        "fixed_shift_variation_proved": True,
        "fixed_base_variation_proved": True,
        "mass_bounds_proved": True,
        "axiswise_routes_rejected": True,
        "joint_mobius_gain_proved": False,
        "full_burnol_bound_proved": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(f"audit mismatch for {key}")
    for marker in REQUIRED_NOTE_MARKERS:
        if marker not in note_text:
            issues.append(f"note missing marker: {marker}")
    sources = payload.get("source_anchors", [])
    for source in (
        "https://doi.org/10.2140/ant.2015.9.2167",
        "https://doi.org/10.1017/fmp.2023.28",
        "https://doi.org/10.5802/aif.2401",
    ):
        if source not in sources:
            issues.append(f"missing primary source {source}")


def check_builder_reproduction(
    payload: dict,
    note_text: str,
    issues: list[str],
) -> None:
    spec = importlib.util.spec_from_file_location("skv_builder", BUILDER)
    if spec is None or spec.loader is None:
        issues.append("could not load builder")
        return
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    if module.build_payload() != payload:
        issues.append("builder payload does not reproduce result")
    if module.render_note() != note_text:
        issues.append("builder note does not reproduce output")
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        result_path = temp / "result.json"
        note_path = temp / "note.md"
        completed = subprocess.run(
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
            check=False,
            timeout=30,
        )
        if completed.returncode != 0:
            issues.append(
                "builder subprocess failed: "
                f"{completed.stderr.strip()}"
            )
        elif (
            json.loads(result_path.read_text(encoding="utf-8")) != payload
            or note_path.read_text(encoding="utf-8") != note_text
        ):
            issues.append("builder subprocess output mismatch")


def check_shift_geometry(issues: list[str]) -> None:
    for size in (2, 3, 5, 8, 13):
        ambient = 2 * size + 1
        rho = math.pi**2 * size / ambient**2
        mu = mobius_values(4 * size)
        for alpha in (0.2, 0.55, 0.85):
            cutoff = min(
                size,
                math.ceil(size ** (1.0 - alpha / 2.0)),
            )
            direct_offdiagonal = 0.0
            shifted_offdiagonal = 0.0
            total_mass = 0.0
            weights: dict[tuple[int, int], float] = {}
            for left in range(size + 1, 4 * size):
                for right in range(left + 1, 4 * size):
                    value = ordinary_kernel(
                        alpha,
                        size,
                        cutoff,
                        left,
                        right,
                    )
                    weights[(left, right - left)] = value
                    total_mass += value
                    direct_offdiagonal += (
                        mu[left] * mu[right] * value
                    )
            if not total_mass < math.pi**2 * size / 2.0 + 1e-12:
                issues.append(
                    f"total mass bound failed at K={size}, alpha={alpha}"
                )
            for shift in range(1, 3 * size - 1):
                support = list(
                    range(size + 1, 4 * size - shift)
                )
                blocks = block_supports(size, shift)
                flattened = [value for block in blocks for value in block]
                if flattened != support:
                    issues.append(
                        f"block partition failed at K={size}, h={shift}"
                    )
                    continue
                values = [weights[(n, shift)] for n in support]
                if any(
                    value <= 0.0 or value > rho + 1e-11
                    for value in values
                ):
                    issues.append(
                        f"height bound failed at K={size}, h={shift}"
                    )
                if not variation(values) < 6.0 * rho + 1e-11:
                    issues.append(
                        f"fixed-shift variation failed at "
                        f"K={size}, h={shift}"
                    )
                current_current, current_future, future_future = blocks
                cc_values = [weights[(n, shift)] for n in current_current]
                if any(
                    right <= left
                    for left, right in zip(cc_values, cc_values[1:])
                ):
                    issues.append(
                        f"cc monotonicity failed at K={size}, h={shift}"
                    )
                if not variation(cc_values) < rho + 1e-11:
                    issues.append(
                        f"cc variation failed at K={size}, h={shift}"
                    )
                cf_values = [weights[(n, shift)] for n in current_future]
                if not variation(cf_values) < 2.0 * rho + 1e-11:
                    issues.append(
                        f"cf variation failed at K={size}, h={shift}"
                    )
                ff_values = [weights[(n, shift)] for n in future_future]
                if any(
                    right >= left
                    for left, right in zip(ff_values, ff_values[1:])
                ):
                    issues.append(
                        f"ff monotonicity failed at K={size}, h={shift}"
                    )
                if not variation(ff_values) < rho + 1e-11:
                    issues.append(
                        f"ff variation failed at K={size}, h={shift}"
                    )
                for n in current_current:
                    index = n - size
                    expected = kernel(
                        size,
                        cutoff,
                        index,
                        index + shift,
                    )
                    if not close(weights[(n, shift)], expected):
                        issues.append(
                            f"cc formula failed at K={size}, h={shift}"
                        )
                for n in current_future:
                    index = n - size
                    expected = (
                        kernel(size, cutoff, index, size)
                        * future_weight(alpha, size, n + shift)
                    )
                    if not close(weights[(n, shift)], expected):
                        issues.append(
                            f"cf formula failed at K={size}, h={shift}"
                        )
                for n in future_future:
                    expected = (
                        kernel(size, cutoff, size, size)
                        * future_weight(alpha, size, n)
                        * future_weight(alpha, size, n + shift)
                    )
                    if not close(weights[(n, shift)], expected):
                        issues.append(
                            f"ff formula failed at K={size}, h={shift}"
                        )
                prefix: list[int] = []
                running = 0
                direct_shift = 0.0
                for n, weight in zip(support, values):
                    running += mu[n] * mu[n + shift]
                    prefix.append(running)
                    direct_shift += mu[n] * mu[n + shift] * weight
                abel_shift = prefix[-1] * values[-1] + sum(
                    prefix[index] * (values[index] - values[index + 1])
                    for index in range(len(values) - 1)
                )
                if not close(direct_shift, abel_shift):
                    issues.append(
                        f"fixed-shift Abel identity failed at "
                        f"K={size}, h={shift}"
                    )
                maximum = max(abs(value) for value in prefix)
                if abs(direct_shift) > 7.0 * rho * maximum + 1e-10:
                    issues.append(
                        f"fixed-shift Abel bound failed at "
                        f"K={size}, h={shift}"
                    )
                if sum(values) >= 13.0 * math.pi**2 / 24.0 + 1e-10:
                    issues.append(
                        f"per-shift mass failed at K={size}, h={shift}"
                    )
                shifted_offdiagonal += direct_shift
            if not close(direct_offdiagonal, shifted_offdiagonal):
                issues.append(
                    f"shift decomposition failed at "
                    f"K={size}, alpha={alpha}"
                )
            for n in range(size + 1, 4 * size - 1):
                shifts = list(range(1, 4 * size - n))
                values = [weights[(n, shift)] for shift in shifts]
                bound = math.pi**2 / ambient + rho
                if not variation(values) < bound + 1e-11:
                    issues.append(
                        f"fixed-base variation failed at K={size}, n={n}"
                    )
                prefix = []
                running = 0
                direct_base = 0.0
                for shift, weight in zip(shifts, values):
                    running += mu[n] * mu[n + shift]
                    prefix.append(running)
                    direct_base += mu[n] * mu[n + shift] * weight
                abel_base = prefix[-1] * values[-1] + sum(
                    prefix[index] * (values[index] - values[index + 1])
                    for index in range(len(values) - 1)
                )
                if not close(direct_base, abel_base):
                    issues.append(
                        f"fixed-base Abel identity failed at K={size}, n={n}"
                    )
                local_maximum = max(abs(value) for value in prefix)
                if abs(direct_base) > (
                    math.pi**2 / size * local_maximum + 1e-10
                ):
                    issues.append(
                        f"fixed-base Abel bound failed at K={size}, n={n}"
                    )


def check_mass_calibration(issues: list[str]) -> None:
    for size in (1, 2, 3, 5, 8, 13, 21):
        for alpha in (0.1, 0.5, 0.9):
            weights = [
                future_weight(alpha, size, n)
                for n in range(2 * size, 4 * size)
            ]
            first = sum(weights)
            second = sum(value * value for value in weights)
            if first > (2 * size + 1) / 2.0 + 1e-12:
                issues.append(
                    f"future l1 bound failed at K={size}, alpha={alpha}"
                )
            expected_second = (
                (2 * size + 1)
                * (4 * size + 1)
                / (12.0 * size)
            )
            if second > expected_second + 1e-12:
                issues.append(
                    f"future l2 bound failed at K={size}, alpha={alpha}"
                )
        if not (
            math.pi**2 / (2 * size + 1)
            + 2 * math.pi**2 * size / (2 * size + 1) ** 2
            < math.pi**2 / size
        ):
            issues.append(f"fixed-base constant failed at K={size}")


def check_formal_core(issues: list[str]) -> None:
    text = FORMAL_CORE.read_text(encoding="utf-8")
    for marker in (
        "Corollary 11.22Z.4",
        "(11.22Z.64)",
        "(11.22Z.85)",
        "fixed-shift variation",
        "joint cancellation across shifts",
    ):
        if marker not in text:
            issues.append(f"formal core missing marker: {marker}")


def main() -> int:
    issues: list[str] = []
    if not RESULT.exists():
        issues.append(f"missing result {RESULT}")
    if not NOTE.exists():
        issues.append(f"missing note {NOTE}")
    if not BUILDER.exists():
        issues.append(f"missing builder {BUILDER}")
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    note_text = NOTE.read_text(encoding="utf-8")
    check_schema(payload, note_text, issues)
    check_builder_reproduction(payload, note_text, issues)
    check_shift_geometry(issues)
    check_mass_calibration(issues)
    check_formal_core(issues)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    print(
        "validated Mertens shift-kernel variation handoff: "
        "45 rows, 0 issues, 39 exact reductions, "
        "2 proof guards, 1 open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
