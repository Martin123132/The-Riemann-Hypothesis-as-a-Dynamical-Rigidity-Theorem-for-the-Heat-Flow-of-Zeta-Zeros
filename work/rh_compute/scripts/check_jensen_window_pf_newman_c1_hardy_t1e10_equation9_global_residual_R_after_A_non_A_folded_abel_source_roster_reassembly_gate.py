#!/usr/bin/env python3
"""Independently check the Abel-zero source-roster reassembly gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_folded_abel_source_roster_reassembly_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()

A = 159_577
B = 5_122_421
L = 2_481_422
FIRST = 622
ORDINARY_LAST = 39_852
A_FIRST = 39_853
TARGET_LAST = 39_894
LAST = 39_936


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_hashes(artifact: dict[str, Any]) -> None:
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")


def independent_source_replay(artifact: dict[str, Any]) -> None:
    source = artifact["source_roster_certificate"]
    require(B - A == 2 * L, "source endpoint arithmetic failed")
    require(source["source_count"] == L + 1, "source count mismatch")

    # Rebuild the endpoint coefficient vector on a different surrogate length.
    length = 11
    h = [Fraction(0) for _ in range(length + 1)]
    h[0] = h[-1] = Fraction(1, 2)
    folded = [Fraction(0) for _ in range(length + 1)]
    for index in range(length):
        folded[index] += Fraction(1, 2)
    for index in range(1, length + 1):
        folded[index] += Fraction(1, 2)
    require(all(left + right == 1 for left, right in zip(h, folded)), "endpoint vector replay failed")

    h0 = Fraction(A + B, 2)
    f0 = L * A + L * (L - 1)
    f1 = L * A + L * (L + 1)
    full = sum((A, B)) * (L + 1) // 2
    require(h0 + Fraction(f0 + f1, 2) == full, "x=0 source identity failed")
    recorded = source["x_zero"]
    require(recorded["H_0"] == h0.numerator, "recorded H_0 mismatch")
    require(recorded["F_0_0"] == f0 and recorded["F_0_1"] == f1, "recorded folded sums mismatch")
    require(recorded["complete_source_sum"] == full, "recorded source sum mismatch")


def independent_ownership_replay(artifact: dict[str, Any]) -> None:
    ownership = artifact["ownership_reassembly_certificate"]
    # Dictionary replay is deliberately different from the builder tuple replay.
    minus_i = {"P_m": -1, "A_m": -1, "A_-m": 0, "B_m": -1}
    low_row = {"P_m": 0, "A_m": 1, "A_-m": 0, "B_m": 1}
    high_row = {"P_m": 0, "A_m": 0, "A_-m": -1, "B_m": 1}
    low_net = {key: minus_i[key] + low_row[key] for key in minus_i}
    high_net = {key: minus_i[key] + high_row[key] for key in minus_i}
    require(low_net == {"P_m": -1, "A_m": 0, "A_-m": 0, "B_m": 0}, "low row failed")
    require(
        high_net == {"P_m": -1, "A_m": -1, "A_-m": -1, "B_m": 0},
        "high row failed",
    )
    require(ownership["ordinary_net_vector"] == [-1, 0, 0, 0], "recorded low vector mismatch")
    require(ownership["A_window_net_vector"] == [-1, -1, -1, 0], "recorded high vector mismatch")
    require(ORDINARY_LAST - FIRST + 1 == 39_231, "ordinary count failed")
    require(LAST - A_FIRST + 1 == 84, "A-window count failed")
    partition = ownership["post_A_positive_notch_route_partition"]
    require(sum(row["count"] for row in partition) == 39_315, "route partition failed")
    require(ownership["independent_norm_permission_created"] is False, "independent norm guard failed")


def independent_transition_replay(artifact: dict[str, Any]) -> None:
    transitions = artifact["transition_schedule_certificate"]
    require(math_gcd(A, B) == 1, "endpoint gcd drift")

    # Within-sequence spacings are fixed.  For cross-sequence spacings, test
    # the two nearest A numerators to every B event without sorting the roster.
    best: tuple[int, int, int] | None = None
    for b_mode in range(FIRST, LAST + 1):
        quotient, _ = divmod(b_mode * A, B)
        for a_mode in (quotient, quotient + 1):
            if FIRST <= a_mode <= TARGET_LAST:
                numerator = abs(b_mode * A - a_mode * B)
                candidate = (numerator, b_mode, a_mode)
                if best is None or candidate < best:
                    best = candidate
    require(best == (441, 20_223, 630), "minimum cross-event numerator drift")
    cross_gap = Fraction(2 * best[0], A * B)
    minimum = min(Fraction(2, B), Fraction(2, A), cross_gap)
    require(minimum == Fraction(882, A * B), "minimum event gap mismatch")
    recorded = transitions["minimum_adjacent_gap"]
    require(Fraction(recorded["numerator"], recorded["denominator"]) == minimum, "recorded gap mismatch")
    require(transitions["total_distinct_events"] == 39_315 + 39_273, "event count mismatch")
    require(A // 4 == TARGET_LAST, "midpoint exit floor mismatch")
    require(LAST - TARGET_LAST == 42, "midpoint-persistent roster mismatch")
    classification = transitions["saved_height_classification"]
    require(classification["interior_joint_saddles"]["mode_range"] == [622, 39_852], "interior class drift")
    require(classification["lower_A_transition"]["mode_range"] == [39_853, 39_894], "lower A class drift")
    require(classification["reflected_A_transition"]["mode_range"] == [39_895, 39_936], "reflected class drift")


def math_gcd(left: int, right: int) -> int:
    while right:
        left, right = right, left % right
    return abs(left)


def main() -> int:
    require(RESULT.is_file() and BUILDER.is_file(), "missing gate files")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact did not pass")
    verify_hashes(artifact)
    independent_source_replay(artifact)
    independent_ownership_replay(artifact)
    independent_transition_replay(artifact)

    decision = artifact["decision"]
    require(decision["half_current_plus_Abel_boundary_layer_equals_complete_source_roster"] is True, "source decision drift")
    require(decision["raw_39315_Fresnel_physical_quadrature_selected"] is False, "raw pilot reselected")
    require(decision["source_roster_minus_owned_bulk_and_A_face_route_selected"] is True, "compressed route lost")
    require(decision["physical_x_quadrature_completed"] is False, "physical quadrature overclaim")
    require(decision["non_A_bound_proved"] is False, "non-A overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked folded Abel source-roster reassembly and transition schedule", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
