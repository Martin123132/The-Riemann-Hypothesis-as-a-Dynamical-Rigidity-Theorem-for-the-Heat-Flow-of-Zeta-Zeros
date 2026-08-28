#!/usr/bin/env python3
"""Independently check the signed-coefficient lower-complement enclosure."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_complement_signed_coefficient_interval_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()
UPPER_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_far_remainder_eight_round_gate.json"
COMPLEMENT_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate.json"

CHECK_ORDER = 4
CHECK_SLABS = 20000
CHECK_CLUSTER_POWER = 5
HEIGHT = arb("10000000000")
L = arb("621.5")
A = arb("39852.5")
Q_LOW = arb("79787.5")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def saved_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def clustered_right(index: int) -> arb:
    width = A - L
    denominator = CHECK_SLABS**CHECK_CLUSTER_POWER
    return A - width * arb((CHECK_SLABS - index) ** CHECK_CLUSTER_POWER) / denominator


def signed_remainder(
    level: dict[int, dict[tuple[int, int], object]], p: arb
) -> arb:
    total = arb(0)
    for index in range(CHECK_SLABS):
        slab_left = clustered_right(index)
        slab_right = clustered_right(index + 1)
        midpoint = (slab_left + slab_right) / 2
        radius = (slab_right - slab_left) / 2
        y_ball = arb(midpoint, radius)
        delta = (old_check.q_curve(slab_right, p) - Q_LOW).lower()
        require(delta > 0, f"changed lower gap lost in slab {index}")
        integrand = arb(0)
        for power, terms in level.items():
            coefficient = abs(old_check.coefficient_value(terms, p, y_ball)).upper()
            integrand += coefficient * old_check.series_bound(delta, power, None)
        total += (slab_right - slab_left) * integrand
    return (total / (2 * arb.pi()) ** CHECK_ORDER).upper()


def gamma_radius(complement: dict[str, object]) -> arb:
    log_ball = arb(
        complement["certificate"]["ordinary_Gamma_defect_packet_log10_upper_ball"][
            "ball"
        ]
    )
    return (log_ball.upper() * arb(10).log()).exp().upper()


def main() -> int:
    resource_mode = old_check.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 448
    ctx.threads = 1
    for path in (RESULT, NOTE, BUILDER, CHECKER, UPPER_RESULT, COMPLEMENT_RESULT):
        require(path.is_file(), f"missing artifact: {path}")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production gate not passed")
    require(
        artifact["decision"]["complete_ordinary_packet_complex_ball_certified"] is True,
        "ordinary packet was not certified",
    )
    require(
        artifact["decision"]["complete_later_joined_packet_enclosed"] is False,
        "later join overpromoted",
    )
    for relpath, expected in artifact["sources"].items():
        path = ROOT / relpath
        require(path.is_file() and file_hash(path) == expected, f"source hash drift: {path}")

    c = artifact["certificate"]
    require(c["integration_rounds"] == 3, "production order drift")
    require(c["remainder_slabs"] == 16384, "production slab count drift")
    require(len(c["coefficient_cancellation_diagnostics"]) == 4, "coefficient census drift")
    require(c["selected_audited_order"] == 3, "audited order selection drift")
    require(len(c["audited_order_remainder_bounds"]) == 8, "order sweep drift")

    p = HEIGHT / (2 * arb.pi())
    require(A * A < p, "physical interval crossed Q minimum")
    table = old_check.coefficient_table(CHECK_ORDER)
    partial = old_check.endpoint_partial(
        table, "R_lower", p, L, A, CHECK_ORDER
    )
    remainder = signed_remainder(table[CHECK_ORDER], p)
    changed_lower = old_check.add_error(partial, remainder)
    saved_lower = saved_complex(c["lower_complementary_tail_ball"])
    require(changed_lower.overlaps(saved_lower), "changed-order lower tail misses production")

    upper = json.loads(UPPER_RESULT.read_text(encoding="utf-8"))
    upper_ball = saved_complex(upper["certificate"]["complete_upper_complement_ball"])
    complement = json.loads(COMPLEMENT_RESULT.read_text(encoding="utf-8"))
    changed_ordinary = old_check.add_error(
        -changed_lower - upper_ball, gamma_radius(complement)
    )
    saved_ordinary = saved_complex(c["complete_ordinary_packet_ball"])
    require(changed_ordinary.overlaps(saved_ordinary), "changed ordinary join misses production")

    altered_enclosure, altered_direct = old_check.altered_sign_witness(
        [arb("2.5") - index for index in range(4)]
    )
    require(altered_enclosure.overlaps(altered_direct), "altered negative-D witness misses")
    altered_width = abs(altered_enclosure - altered_direct)

    text = " ".join(NOTE.read_text(encoding="utf-8").split())
    for fragment in (
        "sup_(y in Y_j)|c_k(y)|",
        "sum_(n>=0)(delta+n)^(-k)",
        "O_join=-T_lower-T_upper",
        "complete ordinary packet",
        "No complete later joined packet",
        "no circle, polygon, fitted constant, or visual pattern supplies",
    ):
        require(fragment in text, f"missing note fragment: {fragment}")

    print(
        "independently checked lower complement and ordinary packet; "
        f"altered negative-D overlap width {altered_width.str(8, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
