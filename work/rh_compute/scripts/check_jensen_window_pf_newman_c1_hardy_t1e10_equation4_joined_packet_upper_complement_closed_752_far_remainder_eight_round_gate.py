#!/usr/bin/env python3
"""Independently check the closed 752-label far-remainder enclosure."""

from __future__ import annotations

import hashlib
import json
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

from flint import acb, arb, ctx

import check_jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as old_check


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_far_remainder_eight_round_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()
TRANSITION_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_label_exterior_transition_layer_gate.json"

CHECK_ORDER = 9
HEIGHT = arb("10000000000")
L = arb("621.5")
C = arb("621.6875")
A = arb("39852.5")
Q_PLUS = arb("2561211.5")
Q_REMOTE = Q_PLUS + 752


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def saved_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def endpoint_sum(family: str, y: arb, p: arb, power: int) -> acb:
    qy = old_check.q_curve(y, p)
    if family == "F_C":
        delta = Q_PLUS - qy
        if y == C:
            value = old_check.finite_sum(delta, 16, 11, 752, power)
        else:
            value = old_check.finite_sum(delta, 2, 1, 752, power)
        return old_check.phase(Q_PLUS, y) * value
    if family == "I_L":
        delta = Q_REMOTE - qy
        return old_check.phase(Q_REMOTE, y) * old_check.eta(delta, power)
    raise RuntimeError(f"unknown family: {family}")


def endpoint_value(
    level: dict[int, dict[tuple[int, int], object]],
    family: str,
    p: arb,
    y: arb,
) -> acb:
    value = acb(0)
    for denominator_power, terms in level.items():
        value += old_check.coefficient_value(terms, p, y) * endpoint_sum(
            family, y, p, denominator_power + 1
        )
    return old_check.carrier(HEIGHT, y) * value


def endpoint_partial(
    table: list[dict[int, dict[tuple[int, int], object]]],
    family: str,
    p: arb,
    left: arb,
) -> acb:
    value = acb(0)
    scale = acb(0, 2 * arb.pi())
    for level in range(CHECK_ORDER):
        boundary = endpoint_value(table[level], family, p, A) - endpoint_value(
            table[level], family, p, left
        )
        value += (-1) ** level * scale ** (-(level + 1)) * boundary
    return value


def changed_ball(
    table: list[dict[int, dict[tuple[int, int], object]]],
    family: str,
    p: arb,
) -> acb:
    if family == "F_C":
        left, q0, count = C, Q_PLUS, 752
    elif family == "I_L":
        left, q0, count = L, Q_REMOTE, None
    else:
        raise RuntimeError(f"unknown family: {family}")
    partial = endpoint_partial(table, family, p, left)
    remainder = old_check.production_remainder(
        table[CHECK_ORDER],
        p,
        left,
        A,
        q0,
        "upper",
        count,
        CHECK_ORDER,
    )
    return old_check.add_error(partial, remainder)


def main() -> int:
    resource_mode = old_check.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 448
    ctx.threads = 1
    for path in (RESULT, NOTE, BUILDER, CHECKER, TRANSITION_RESULT):
        require(path.is_file(), f"missing artifact: {path}")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production gate not passed")
    require(
        artifact["decision"]["complete_upper_complement_complex_ball_certified"] is True,
        "upper complement was not certified",
    )
    require(
        artifact["decision"]["lower_complementary_tail_enclosed"] is False,
        "lower tail overpromoted",
    )
    for relpath, expected in artifact["sources"].items():
        path = ROOT / relpath
        require(path.is_file() and file_hash(path) == expected, f"source hash drift: {path}")

    c = artifact["certificate"]
    require(c["integration_rounds"] == 8, "production order drift")
    require(c["transition_label_count"] == 752, "transition count drift")
    require(
        c["endpoint_families"]["F_C_finite_752"]["endpoint_C"]["root_numerator"] == 11,
        "C root numerator drift",
    )
    require(
        c["endpoint_families"]["F_C_finite_752"]["endpoint_C"]["cycles"] == 47,
        "C cycle count drift",
    )
    require(
        c["endpoint_families"]["F_C_finite_752"]["endpoint_a"]["cycles"] == 376,
        "a cycle count drift",
    )

    p = HEIGHT / (2 * arb.pi())
    require(A * A < p, "physical interval crossed the Q minimum")
    table = old_check.coefficient_table(CHECK_ORDER)
    changed: dict[str, acb] = {}
    for family in ("F_C", "I_L"):
        changed[family] = changed_ball(table, family, p)
        saved = saved_complex(c["families"][family]["tail_integral_ball"])
        require(changed[family].overlaps(saved), f"changed-order {family} replay misses")

    transition = json.loads(TRANSITION_RESULT.read_text(encoding="utf-8"))
    g752 = saved_complex(transition["certificate"]["grouped_752_label_transition_ball"])
    changed_upper = g752 + changed["F_C"] + changed["I_L"]
    saved_upper = saved_complex(c["complete_upper_complement_ball"])
    require(changed_upper.overlaps(saved_upper), "changed-order upper join misses production")

    altered_enclosure, altered_direct = old_check.altered_sign_witness(
        [arb("7.5") + index for index in range(4)]
    )
    require(altered_enclosure.overlaps(altered_direct), "altered recurrence witness misses")
    altered_width = abs(altered_enclosure - altered_direct)

    text = " ".join(NOTE.read_text(encoding="utf-8").split())
    for fragment in (
        "T_upper=G_752+F_C+I_L",
        "F_C remainder <=",
        "I_L remainder <=",
        "complete upper complement",
        "No circle, polygon, fitted constant, or visual pattern supplies",
        "No lower complementary-tail",
    ):
        require(fragment in text, f"missing note fragment: {fragment}")

    print(
        "independently checked complete upper complement; "
        f"altered recurrence overlap width {altered_width.str(8, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
