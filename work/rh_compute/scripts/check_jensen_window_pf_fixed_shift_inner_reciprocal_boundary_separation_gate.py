#!/usr/bin/env python3
"""Validate the fixed-shift inner/reciprocal-boundary separation gate."""

from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_fixed_shift_inner_"
    "reciprocal_boundary_separation_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_fixed_shift_inner_"
    "reciprocal_boundary_separation_gate.md"
)

REQUIRED_IDS = {
    f"fsirbsg_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "parameters"),
        (2, "toy_completed_function"),
        (3, "symmetry"),
        (4, "zero_quartet"),
        (5, "fixed_shift_quotient"),
        (6, "boundary_zero_cancellation"),
        (7, "quotient_simplification"),
        (8, "pole_location"),
        (9, "boundary_unimodularity"),
        (10, "rational_inner"),
        (11, "causal_difference"),
        (12, "exact_hardy_energy"),
        (13, "reciprocal_boundary_energy"),
        (14, "shifted_line_zero"),
        (15, "reciprocal_divergence"),
        (16, "fixed_shift_nonpromotion"),
        (17, "cofinal_distinction"),
        (18, "open_arithmetic_target"),
    )
}

ROLE_STATUS = {
    "exact_definition": "available_exact",
    "exact_identity": "available_exact",
    "exact_consequence": "available_exact",
    "proof_guard": "guard_validated",
    "open_target": "open_target",
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Fixed-Shift Inner/Reciprocal-Boundary Separation Gate",
    "F(z-omega)/F(z+omega)",
    "boundary-zero factor cancels",
    "|Q(it)|=1",
    "rational inner",
    "integral_R |Q(it)-1|^2 dt=16*pi*omega",
    "I_(omega,T)=infinity",
    "same boundary zero is removable",
    "universal norm multiplier",
    "exceptional cancellation shift",
    "sum_N P_(alpha/2,N)/N<infinity",
    "one open arithmetic target",
    "not a proof of RH",
)


def close(
    left: complex | float,
    right: complex | float,
    tolerance: float = 2e-10,
) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def toy_f(z: complex, omega: float, height: float) -> complex:
    return (
        ((z - omega) ** 2 + height * height)
        * ((z + omega) ** 2 + height * height)
    )


def raw_q(z: complex, omega: float, height: float) -> complex:
    return toy_f(z - omega, omega, height) / toy_f(
        z + omega, omega, height
    )


def reduced_q(z: complex, omega: float, height: float) -> complex:
    return (
        ((z - 2.0 * omega) ** 2 + height * height)
        / ((z + 2.0 * omega) ** 2 + height * height)
    )


def check_schema(payload: dict, note: str, issues: list[str]) -> None:
    expected_kind = (
        "jensen_window_pf_fixed_shift_inner_"
        "reciprocal_boundary_separation_gate"
    )
    if payload.get("kind") != expected_kind:
        issues.append("bad kind")
    rows = payload.get("rows", [])
    if len(rows) != 18:
        issues.append(f"expected 18 rows, found {len(rows)}")
    ids = {item.get("id") for item in rows}
    if ids != REQUIRED_IDS:
        issues.append(
            "row id mismatch: "
            f"missing={sorted(REQUIRED_IDS-ids)}, "
            f"extra={sorted(ids-REQUIRED_IDS)}"
        )
    for item in rows:
        role = item.get("role")
        expected_status = ROLE_STATUS.get(role)
        if expected_status is None:
            issues.append(f"{item.get('id')}: bad role {role!r}")
        elif item.get("status") != expected_status:
            issues.append(
                f"{item.get('id')}: role/status mismatch "
                f"{role!r}/{item.get('status')!r}"
            )
        for field in ("statement", "proof_boundary"):
            if not item.get(field):
                issues.append(f"{item.get('id')}: missing {field}")

    expected_audit = {
        "row_count": 18,
        "exact_reduction_count": 15,
        "proof_guard_count": 2,
        "open_arithmetic_target_count": 1,
        "rational_inner_proved": True,
        "causal_energy_finite": True,
        "reciprocal_boundary_energy_finite": False,
        "fixed_shift_promotion_valid": False,
        "cofinal_full_burnol_implication_rejected": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    audit = payload.get("audit", {})
    for key, value in expected_audit.items():
        if audit.get(key) != value:
            issues.append(
                f"audit {key}: expected {value!r}, "
                f"found {audit.get(key)!r}"
            )
    for needle in REQUIRED_NOTE:
        if needle not in note:
            issues.append(f"note missing {needle!r}")


def check_polynomial_and_quotient(issues: list[str]) -> None:
    for omega, height in ((0.21, 1.7), (0.63, 3.2)):
        samples = (
            complex(0.3, 0.4),
            complex(1.2, -2.1),
            complex(-0.7, 1.1),
        )
        for z in samples:
            value = toy_f(z, omega, height)
            if not close(toy_f(-z, omega, height), value):
                issues.append("even symmetry failed")
            if not close(
                toy_f(z.conjugate(), omega, height),
                value.conjugate(),
            ):
                issues.append("real symmetry failed")
            if not close(
                raw_q(z, omega, height),
                reduced_q(z, omega, height),
            ):
                issues.append("fixed-shift quotient cancellation failed")
            expected_difference = (
                -8.0
                * omega
                * z
                / ((z + 2.0 * omega) ** 2 + height * height)
            )
            if not close(
                reduced_q(z, omega, height) - 1.0,
                expected_difference,
            ):
                issues.append("causal quotient-difference identity failed")

        for real_sign in (-1.0, 1.0):
            for imag_sign in (-1.0, 1.0):
                root = complex(real_sign * omega, imag_sign * height)
                if abs(toy_f(root, omega, height)) > 2e-11:
                    issues.append(f"toy root failed at {root}")
                derivative_probe = toy_f(
                    root + 1e-7, omega, height
                ) / 1e-7
                if abs(derivative_probe) < 1e-5:
                    issues.append(f"toy root is not simple at {root}")

        for index in range(-120, 121):
            t = index * height / 31.0
            modulus = abs(reduced_q(complex(0.0, t), omega, height))
            if not close(modulus, 1.0, tolerance=2e-11):
                issues.append(
                    f"boundary unimodularity failed at omega={omega}, t={t}"
                )
                break

        for x in (0.05, 0.3, 1.0, 4.0):
            for y in (-2.0 * height, -0.4, 0.0, 0.7, 1.5 * height):
                modulus = abs(reduced_q(complex(x, y), omega, height))
                if modulus >= 1.0 + 3e-12:
                    issues.append(
                        f"right-half-plane inner bound failed at "
                        f"omega={omega}, z={x}+{y}i"
                    )

        for pole in (
            complex(-2.0 * omega, height),
            complex(-2.0 * omega, -height),
        ):
            denominator = (pole + 2.0 * omega) ** 2 + height * height
            if abs(denominator) > 2e-12 or pole.real >= 0.0:
                issues.append(f"pole-location check failed at {pole}")


def transformed_boundary_integrand(
    u: float,
    omega: float,
    height: float,
) -> float:
    endpoint = math.pi / 2.0
    if abs(abs(u) - endpoint) < 1e-15:
        return 64.0 * omega * omega
    t = math.tan(u)
    d_minus = (t - height) ** 2 + 4.0 * omega * omega
    d_plus = (t + height) ** 2 + 4.0 * omega * omega
    value = 64.0 * omega * omega * t * t / (d_minus * d_plus)
    return value / (math.cos(u) ** 2)


def simpson_integral(
    omega: float,
    height: float,
    panels: int = 20000,
) -> float:
    left = -math.pi / 2.0
    right = math.pi / 2.0
    step = (right - left) / panels
    total = (
        transformed_boundary_integrand(left, omega, height)
        + transformed_boundary_integrand(right, omega, height)
    )
    for index in range(1, panels):
        weight = 4.0 if index % 2 else 2.0
        total += weight * transformed_boundary_integrand(
            left + index * step,
            omega,
            height,
        )
    return total * step / 3.0


def check_hardy_energy(issues: list[str]) -> None:
    for omega, height in ((0.17, 0.9), (0.41, 2.6), (0.82, 5.0)):
        numerical = simpson_integral(omega, height)
        expected = 16.0 * math.pi * omega
        if not close(numerical, expected, tolerance=2e-8):
            issues.append(
                f"Hardy boundary integral failed at omega={omega}: "
                f"numerical={numerical}, expected={expected}"
            )

        for t in (-3.1, -height / 2.0, 0.4, 1.7 * height):
            d_minus = (t - height) ** 2 + 4.0 * omega * omega
            d_plus = (t + height) ** 2 + 4.0 * omega * omega
            left = t * t / (d_minus * d_plus)
            right = t / (4.0 * height) * (
                1.0 / d_minus - 1.0 / d_plus
            )
            if not close(left, right, tolerance=3e-12):
                issues.append("Cauchy-kernel partial fraction failed")


def reciprocal_integrand(
    t: float,
    omega: float,
    height: float,
) -> float:
    value = toy_f(complex(omega, t), omega, height)
    return 1.0 / ((1.0 + t * t) * abs(value) ** 2)


def check_reciprocal_divergence(issues: list[str]) -> None:
    for omega, height in ((0.23, 1.4), (0.71, 3.7)):
        for sign in (-1.0, 1.0):
            root = complex(omega, sign * height)
            if abs(toy_f(root, omega, height)) > 2e-11:
                issues.append("shifted-line zero check failed")

        expected_constant = 1.0 / (
            64.0
            * omega
            * omega
            * height
            * height
            * (omega * omega + height * height)
            * (1.0 + height * height)
        )
        scaled_values = []
        for epsilon in (2e-3, 1e-3, 5e-4, 2.5e-4):
            scaled_values.append(
                epsilon
                * epsilon
                * reciprocal_integrand(
                    height + epsilon, omega, height
                )
            )
        if abs(scaled_values[-1] / expected_constant - 1.0) > 4e-3:
            issues.append(
                f"reciprocal pole asymptotic failed at omega={omega}: "
                f"ratio={scaled_values[-1] / expected_constant}"
            )
        if not all(
            scaled_values[index + 1] > 0.95 * scaled_values[index]
            for index in range(len(scaled_values) - 1)
        ):
            issues.append("reciprocal scaled-pole sequence is unstable")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []
    check_schema(payload, note, issues)
    check_polynomial_and_quotient(issues)
    check_hardy_energy(issues)
    check_reciprocal_divergence(issues)
    for issue in issues:
        print(f"FIXED-SHIFT-INNER-RECIPROCAL {issue}")
    audit = payload.get("audit", {})
    print(
        "validated fixed-shift inner/reciprocal-boundary separation gate: "
        f"{len(payload.get('rows', []))} rows, {len(issues)} issues, "
        f"{audit.get('exact_reduction_count', 0)} exact reductions, "
        f"{audit.get('proof_guard_count', 0)} proof guards, "
        f"{audit.get('open_arithmetic_target_count', 0)} open target"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
