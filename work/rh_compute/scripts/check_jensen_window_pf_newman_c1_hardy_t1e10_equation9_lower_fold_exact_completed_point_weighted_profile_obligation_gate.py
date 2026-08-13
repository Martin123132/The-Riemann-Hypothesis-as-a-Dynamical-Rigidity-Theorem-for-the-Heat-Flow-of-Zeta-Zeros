#!/usr/bin/env python3
"""Independently check the completed-point/weighted-profile obligation."""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_completed_point_weighted_profile_obligation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_ball(text: str) -> tuple[Decimal, Decimal]:
    body = text.strip().removeprefix("[").removesuffix("]").strip()
    if "+/-" not in body:
        return Decimal(body), Decimal(0)
    midpoint, radius = body.split("+/-", 1)
    return Decimal(midpoint.strip()), Decimal(radius.strip())


def exact_nonmatching_counterprofile() -> None:
    """Use different rationals from the builder and sparse Fourier algebra."""

    center = 39_895
    mode_coefficients = {
        center + 198: Fraction(43, 47),
        center + 199: Fraction(-43, 47),
    }
    require(sum(mode_coefficients.values(), Fraction()) == 0, "s=0 endpoint changed")
    require(sum(mode_coefficients.values(), Fraction()) == 0, "s=1 endpoint changed")
    require(39_894 not in mode_coefficients, "event-zero coefficient changed")

    weights = {40_093: Fraction(-13, 37), 40_094: Fraction(17, 41)}
    direct = sum((weights[mode] * value for mode, value in mode_coefficients.items()), Fraction())
    expected = Fraction(43, 47) * (weights[40_093] - weights[40_094])
    require(direct == expected != 0, "independent weighted counterprofile shift failed")

    # Independently check the elementary L-infinity/L1 transfer on a finite
    # rational step model; this tests the inequality without using builder code.
    profile_values = [Fraction(2, 7), Fraction(-3, 11), Fraction(5, 13)]
    kernel_values = [Fraction(-7, 17), Fraction(11, 19), Fraction(3, 23)]
    pairing = sum((abs(g * p) for g, p in zip(profile_values, kernel_values)), Fraction(0))
    bound = max(map(abs, profile_values)) * sum(map(abs, kernel_values), Fraction(0))
    require(pairing <= bound, "independent Linfinity-L1 transfer failed")


def main() -> None:
    getcontext().prec = 90
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")

    exact_nonmatching_counterprofile()
    certificate = artifact["certificate"]
    weight_difference = parse_ball(certificate["last_two_weight_difference_absolute_ball"])
    kernel_l1 = parse_ball(certificate["two_piece_weighted_kernel_L1_bound_ball"])
    completed_value = parse_ball(certificate["uniform_completed_scalar_reduced_error_ball"])
    spacing = parse_ball(certificate["detuning_spacing_h_ball"])
    require(weight_difference[0] - weight_difference[1] > Decimal("1.917e-7"), "weight distinction failed")
    require(kernel_l1[0] + kernel_l1[1] < Decimal("3.840283e-5"), "kernel L1 failed")
    require(completed_value[0] + completed_value[1] < Decimal("9.497e-12"), "completed point error failed")
    spacing_lower = spacing[0] - spacing[1]
    spacing_upper = spacing[0] + spacing[1]
    require(spacing_lower > Decimal("0.054"), "spacing interval failed")
    independent_scale_lower = Decimal(1) / spacing_upper
    independent_scale_upper = Decimal(1) / spacing_lower
    stored_scale = [Decimal(value) for value in certificate["canonical_profile_scale_interval"]]
    require(stored_scale[0] == independent_scale_lower, "profile scale lower drift")
    require(stored_scale[1] == independent_scale_upper, "profile scale upper drift")
    canonical_profile_upper = (completed_value[0] + completed_value[1]) * independent_scale_upper
    raw_transfer_upper = (kernel_l1[0] + kernel_l1[1]) * independent_scale_upper
    require(canonical_profile_upper == Decimal(certificate["uniform_completed_canonical_profile_value_upper_bound"]), "canonical profile bound drift")
    require(raw_transfer_upper == Decimal(certificate["raw_fold_profile_transfer_coefficient_upper_bound"]), "raw transfer drift")
    require(canonical_profile_upper < Decimal("1.759e-10"), "canonical profile bound failed")
    require(raw_transfer_upper < Decimal("7.112e-4"), "raw transfer bound failed")

    decision = artifact["decision"]
    require(decision["completed_scalar_value_controls_weighted_profile_correction"] is False, "point insufficiency drift")
    require(decision["conditional_whole_kernel_transfer_exact"] is True, "conditional transfer drift")
    require(decision["physical_exact_common_profile_operator_identity_proved"] is False, "operator proof boundary drift")
    require(decision["uniform_exact_minus_beta4_profile_bound_proved"] is False, "profile proof boundary drift")
    require(decision["rh_implication"] is False, "RH boundary drift")
    print("independently checked exact completed-point versus weighted-profile obligation", flush=True)


if __name__ == "__main__":
    main()
