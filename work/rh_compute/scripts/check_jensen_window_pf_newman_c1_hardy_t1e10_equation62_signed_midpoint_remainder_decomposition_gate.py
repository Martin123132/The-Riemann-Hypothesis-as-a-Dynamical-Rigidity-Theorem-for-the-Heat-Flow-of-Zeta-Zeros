#!/usr/bin/env python3
"""Independently check the signed midpoint-remainder decomposition."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation62_signed_midpoint_remainder_decomposition_gate as gate


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing signed decomposition result")
    require(gate.NOTE.is_file(), "missing signed decomposition note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_equation62_signed_midpoint_remainder_decomposition_gate", "kind drift")
    require(artifact["status"] == "finite_diagnostic_midpoint_residual_resolved_as_opposite_sign_channels", "status drift")
    require(bool(artifact["passed"]), "signed decomposition gate failed")
    require(artifact["scope"]["classical_collar"] == list(gate.dgate.COLLAR), "collar roster drift")
    require(len(artifact["output_rows"]) == gate.EXPECTED_OUTPUTS, "output roster drift")

    ctx.dps = CHECK_PRECISION
    ctx.threads = 1
    complementary: list[arb] = []
    D_boundary: list[arb] = []
    midpoint: list[arb] = []
    triangles: list[arb] = []
    retained: list[arb] = []
    cancelled: list[arb] = []
    for output_index, saved in enumerate(artifact["output_rows"], start=1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        fresh_cells = gate.cells.output_atlas(offset, CHECK_PRECISION)
        fresh_D = gate.dgate.output_transport(target_t, CHECK_PRECISION)
        final = fresh_cells["prefix_rows"][-1]

        source_remainder = arb(final["prefix_minus_source_endpoint_target_ball"])
        direct_midpoint = arb(final["prefix_minus_midpoint_ceil_target_ball"])
        exact_collar = arb(final["midpoint_ceil_target_minus_source_endpoint_target_ball"])
        transported_D = arb(fresh_D["transformed_D_collar_real_ball"])
        require(exact_collar.overlaps(arb(fresh_D["exact_theta_classical_collar_ball"])), "fresh collar handoff drift")

        fresh_complementary = source_remainder - transported_D
        fresh_D_boundary = transported_D - exact_collar
        fresh_midpoint = fresh_complementary + fresh_D_boundary
        fresh_triangle = abs(fresh_complementary) + abs(fresh_D_boundary)
        fresh_retained = abs(fresh_midpoint) / fresh_triangle
        fresh_cancelled = 1 - fresh_retained

        require(saved["output_index"] == output_index, "saved output index drift")
        require(saved["target_t"] == target_t, "saved target height drift")
        require(fresh_complementary < 0, "fresh complementary sign drift")
        require(fresh_D_boundary > 0, "fresh D-boundary sign drift")
        require(fresh_midpoint > 0, "fresh midpoint sign drift")
        require(fresh_midpoint.overlaps(direct_midpoint), "fresh signed identity drift")
        require(abs(fresh_midpoint) < gate.REQUESTED_TOLERANCE, "fresh midpoint target failed")
        require(fresh_triangle > gate.REQUESTED_TOLERANCE, "fresh triangle unexpectedly closes")

        for key, value in (
            ("complementary_channel_source_minus_D_ball", fresh_complementary),
            ("D_minus_exact_classical_collar_ball", fresh_D_boundary),
            ("signed_channel_sum_ball", fresh_midpoint),
            ("direct_midpoint_remainder_ball", direct_midpoint),
            ("independent_absolute_triangle_ball", fresh_triangle),
            ("signed_remainder_retained_fraction_ball", fresh_retained),
            ("opposite_sign_cancellation_fraction_ball", fresh_cancelled),
        ):
            require(arb(saved[key]).overlaps(value), f"fresh saved-row drift: {key}")

        complementary.append(fresh_complementary)
        D_boundary.append(fresh_D_boundary)
        midpoint.append(fresh_midpoint)
        triangles.append(fresh_triangle)
        retained.append(fresh_retained)
        cancelled.append(fresh_cancelled)

    aggregate = artifact["aggregate"]
    require(aggregate["output_count"] == gate.EXPECTED_OUTPUTS, "aggregate output count drift")
    require(aggregate["opposite_sign_channel_count"] == gate.EXPECTED_OUTPUTS, "aggregate sign count drift")
    require(aggregate["signed_midpoint_below_0p005_count"] == gate.EXPECTED_OUTPUTS, "aggregate signed count drift")
    require(aggregate["independent_triangle_above_0p005_count"] == gate.EXPECTED_OUTPUTS, "aggregate triangle count drift")
    for key, value in (
        ("minimum_complementary_channel_ball", min(complementary, key=lambda x: x.lower())),
        ("maximum_complementary_channel_ball", max(complementary, key=lambda x: x.upper())),
        ("minimum_D_boundary_channel_ball", min(D_boundary, key=lambda x: x.lower())),
        ("maximum_D_boundary_channel_ball", max(D_boundary, key=lambda x: x.upper())),
        ("minimum_signed_midpoint_remainder_ball", min(midpoint, key=lambda x: x.lower())),
        ("maximum_signed_midpoint_remainder_ball", max(midpoint, key=lambda x: x.upper())),
        ("minimum_independent_triangle_ball", min(triangles, key=lambda x: x.lower())),
        ("maximum_independent_triangle_ball", max(triangles, key=lambda x: x.upper())),
        ("minimum_signed_remainder_retained_fraction_ball", min(retained, key=lambda x: x.lower())),
        ("maximum_signed_remainder_retained_fraction_ball", max(retained, key=lambda x: x.upper())),
        ("minimum_opposite_sign_cancellation_fraction_ball", min(cancelled, key=lambda x: x.lower())),
        ("maximum_opposite_sign_cancellation_fraction_ball", max(cancelled, key=lambda x: x.upper())),
    ):
        require(arb(aggregate[key]).overlaps(value), f"fresh aggregate drift: {key}")

    decision = artifact["decision"]
    require(bool(decision["source_remainder_minus_D_strict_negative_all_outputs"]), "negative sign result lost")
    require(bool(decision["D_minus_exact_classical_collar_strict_positive_all_outputs"]), "positive sign result lost")
    require(bool(decision["signed_sum_equals_midpoint_remainder_all_outputs"]), "signed identity result lost")
    require(bool(decision["midpoint_remainder_below_0p005_all_outputs"]), "midpoint target result lost")
    require(not bool(decision["separate_absolute_triangle_bound_closes_target"]), "triangle bound overpromotion")
    require(not bool(decision["selector_split_changes_diagnostic_midpoint_residual"]), "selector split falsely changes diagnostic residual")
    require(not bool(decision["channel_sign_uniformity_is_source_aligned_next_theorem"]), "auxiliary channel signs promoted to source-aligned theorem")
    require(not bool(decision["diagnostic_midpoint_is_a_published_cutoff_rule"]), "diagnostic midpoint promoted to source rule")
    require(not bool(decision["height_uniform_sign_pairing_proved"]), "height-uniform overpromotion")
    require(not bool(decision["finite_gate_has_rh_implication"]), "RH overpromotion")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file: {path}")
            require(file_hash(path) == record["sha256"], f"hash drift: {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "Pure algebra", "strictly negative", "strictly positive", "above `0.005` everywhere", "split is auxiliary", "paper's intrinsic error", "no implication"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated signed midpoint remainder decomposition independently: "
        f"outputs={gate.EXPECTED_OUTPUTS}, opposite-sign={len(midpoint)}, "
        f"signed<0.005={len(midpoint)}, triangle>0.005={len(triangles)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
