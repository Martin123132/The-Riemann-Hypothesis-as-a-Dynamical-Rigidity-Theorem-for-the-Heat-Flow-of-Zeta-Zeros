#!/usr/bin/env python3
"""Independent high-precision replay of the local height event atlas."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_local_height_event_atlas_gate as gate


CHECK_PRECISION_BITS = 512


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(num: int, den: int) -> arb:
    return arb(num) / den


def q_value(height: arb, y: arb) -> arb:
    return y + height / (2 * arb.pi() * y)


def q_prime(height: arb, y: arb) -> arb:
    return 1 - height / (2 * arb.pi() * y * y)


def margin(height: arb, y: arb, q: arb) -> arb:
    detuning = q - q_value(height, y)
    return detuning * detuning - 64 * abs(q_prime(height, y))


def interval_ball(lower: arb, upper: arb) -> arb:
    midpoint = (lower + upper) / 2
    return arb(midpoint, ((upper - lower) / 2).upper())


def overlap_record(record: dict[str, str], value: arb, label: str) -> None:
    require(arb(record["ball"]).overlaps(value), f"saved {label} misses independent replay")


def main() -> None:
    ctx.prec = CHECK_PRECISION_BITS
    require(gate.RESULT.is_file(), f"missing artifact: {gate.RESULT}")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["passed"], "artifact is not passed")
    require(artifact["kind"] == gate.STEM, "artifact kind drift")

    require(artifact["sources"]["builder"]["sha256"] == file_hash(gate.BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(Path(__file__)), "checker hash drift")
    for name, path in gate.DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {path}")
        saved = artifact["dependencies"][name]
        require(saved["sha256"] == file_hash(path), f"dependency hash drift: {name}")

    c = artifact["certificate"]
    cell = c["primary_open_cell"]
    guards = c["uniform_combinatorial_guards"]
    route = c["uniform_route_guards"]

    pi = arb.pi()
    t0 = arb(gate.HEIGHT)
    a = rational(gate.LOWER_A_CELL_NUM, gate.LOWER_A_CELL_DEN)
    U = rational(gate.UPPER_A_CELL_NUM, gate.UPPER_A_CELL_DEN)
    L = rational(gate.L_NUM, gate.L_DEN)
    core_split = rational(gate.CORE_SPLIT_NUM, gate.CORE_SPLIT_DEN)
    far_split = rational(gate.FAR_SPLIT_NUM, gate.FAR_SPLIT_DEN)
    q_minus = rational(gate.Q_MINUS_NUM, gate.Q_MINUS_DEN)
    q_plus = rational(gate.Q_PLUS_NUM, gate.Q_PLUS_DEN)

    # Independent formulas use fixed-cell Q contacts instead of tau_m(alpha).
    lower = 2 * pi * U * gate.ORDINARY_END
    upper = 2 * pi * a * gate.TRANSITION_END
    lower_tau = pi * gate.ORDINARY_END * (gate.A - 2 * gate.ORDINARY_END)
    upper_tau = pi * gate.TRANSITION_END * (gate.A - 2 * gate.TRANSITION_END)
    require((lower - lower_tau).contains(0), "lower alternative event formula failed")
    require((upper - upper_tau).contains(0), "upper alternative event formula failed")
    require(lower < t0 < upper, "saved height left independently reconstructed cell")
    overlap_record(cell["lower_ball"], lower, "lower event")
    overlap_record(cell["upper_ball"], upper, "upper event")
    overlap_record(cell["saved_height_distance_above_lower_ball"], t0 - lower, "lower distance")
    overlap_record(cell["saved_height_distance_below_upper_ball"], upper - t0, "upper distance")
    overlap_record(cell["width_ball"], upper - lower, "cell width")

    def roots(alpha: int, height: arb) -> tuple[arb, arb]:
        discriminant = (arb(alpha) ** 2 - 8 * height / pi).sqrt()
        return (arb(alpha) - discriminant) / 4, (arb(alpha) + discriminant) / 4

    rlm, rlp = roots(gate.A, lower)
    rum, rup = roots(gate.A, upper)
    require((rlm - gate.ORDINARY_END).contains(0), "lower-wall lower root failed")
    require((rlp - U).contains(0), "lower-wall upper root failed")
    require((rum - a).contains(0), "upper-wall lower root failed")
    require((rup - gate.TRANSITION_END).contains(0), "upper-wall upper root failed")

    center_lo = (lower / (2 * pi)).sqrt()
    center_hi = (upper / (2 * pi)).sqrt()
    require(center_lo > gate.CLASSICAL_END and center_hi < gate.CLASSICAL_END + 1, "classical floor drift")
    overlap_record(guards["classical_center_at_lower_ball"], center_lo, "lower center")
    overlap_record(guards["classical_center_at_upper_ball"], center_hi, "upper center")

    bll, bul = roots(gate.B, lower)
    blh, buh = roots(gate.B, upper)
    require(arb(621) < bll < blh < arb(622), "B lower floor drift")
    require(arb(2_560_588) < buh < bul < arb(2_560_589), "B upper floor drift")
    overlap_record(guards["B_lower_saddle_at_lower_ball"], bll, "B lower saddle at lower wall")
    overlap_record(guards["B_lower_saddle_at_upper_ball"], blh, "B lower saddle at upper wall")
    overlap_record(guards["B_upper_saddle_at_lower_ball"], bul, "B upper saddle at lower wall")
    overlap_record(guards["B_upper_saddle_at_upper_ball"], buh, "B upper saddle at upper wall")

    existing = q_value(lower, L) - (q_plus + 229)
    next_gap = q_plus + 230 - q_value(upper, L)
    lower_gap = q_value(lower, a) - (q_minus - 1)
    core_gap = q_plus - q_value(upper, core_split)
    require(existing > 0 and next_gap > 0, "230-label census is not stable")
    require(lower_gap > 0, "lower complement gained a stationary label")
    require(core_gap > 0, "first upper saddle escaped the core split")
    overlap_record(guards["existing_230th_label_margin_at_lower_ball"], existing, "230th-label margin")
    overlap_record(guards["next_231st_label_gap_at_upper_ball"], next_gap, "231st-label gap")
    overlap_record(guards["lower_complement_nonstationary_gap_at_lower_ball"], lower_gap, "lower gap")
    overlap_record(guards["first_upper_saddle_before_core_split_gap_at_upper_ball"], core_gap, "core gap")

    require((q_value(lower, U) - q_minus).contains(0), "lower Q(U) contact failed")
    require((q_value(upper, a) - q_minus).contains(0), "upper Q(a) contact failed")
    require(lower / (2 * pi) > a * a, "Q changed monotonicity on physical interval")

    t_cell = interval_ball(lower, upper)
    finite = margin(t_cell, far_split, q_plus)
    remote = margin(t_cell, L, q_plus + 752)
    prev_grid = margin(t_cell, L + rational(2, 16), q_plus)
    prev_count = margin(t_cell, L, q_plus + 736)
    require(finite.lower() > 0 and remote.lower() > 0, "selected 752 launches failed")
    require(prev_grid.upper() < 0 and prev_count.upper() < 0, "752 minimality predecessor changed sign")
    overlap_record(route["finite_752_eight_width_margin_ball"], finite, "finite 752 margin")
    overlap_record(route["remote_752_eight_width_margin_ball"], remote, "remote 752 margin")
    overlap_record(route["previous_grid_eight_width_margin_ball"], prev_grid, "previous-grid margin")
    overlap_record(route["previous_count_eight_width_margin_ball"], prev_count, "previous-count margin")

    require(math.gcd(9, 16) == math.gcd(11, 16) == 1, "primitive order-16 arithmetic failed")
    require(240 % 16 == 0 and 752 % 16 == 0, "closed-count divisibility failed")
    require(240 // 16 == 15 and 752 // 16 == 47, "closed-count quotient failed")

    sixth_root = (t_cell.log() / 6).exp()
    lambda_star = (sixth_root + (sixth_root * sixth_root - 4).sqrt()) / 2
    y_star = (t_cell / (2 * pi)).sqrt() / lambda_star
    require(rational(1719, 2) < y_star < rational(1721, 2), "diagnostic route split drift")
    overlap_record(route["diagnostic_unequal_truncation_y_star_ball"], y_star, "diagnostic y star")

    require(guards["ordinary_roster"] == [622, 39852], "saved ordinary roster drift")
    require(guards["transition_roster"] == [39853, 39936], "saved transition roster drift")
    require(guards["transition_half_counts"] == [42, 42], "saved half counts drift")
    require(guards["upper_complement_stationary_count"] == 230, "saved upper census drift")
    require(not c["identity_scope"]["uniform_numerical_enclosures"], "packet bounds were overpromoted")
    require(not c["identity_scope"]["saved_height_sign_transported"], "saved sign was overpromoted")
    require(not artifact["decision"]["all_height_theorem"], "all-height theorem was overpromoted")
    require(not artifact["decision"]["rh_implication"], "RH implication was overpromoted")

    print(
        "independently checked local height event atlas: maximal same-roster cell, "
        "230-label census, rational closures, and 752 launch signs"
    )


if __name__ == "__main__":
    main()
