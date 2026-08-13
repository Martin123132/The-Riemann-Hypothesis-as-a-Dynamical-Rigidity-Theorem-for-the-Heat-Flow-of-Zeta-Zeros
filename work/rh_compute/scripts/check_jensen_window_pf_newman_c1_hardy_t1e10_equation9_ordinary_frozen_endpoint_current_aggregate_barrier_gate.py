#!/usr/bin/env python3
"""Independent replay of the frozen endpoint-current aggregate barrier."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_frozen_endpoint_current_aggregate_barrier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()
T = 10_000_000_000
A = 159_577
B = 5_122_421
SEGMENTS = (
    ("B_collar", 622, 999),
    ("low_bulk", 1_000, 9_999),
    ("middle_bulk", 10_000, 29_999),
    ("upper_bulk", 30_000, 38_999),
    ("A_collar", 39_000, 39_694),
    ("fold_overlap", 39_695, 39_852),
    ("lower_transition", 39_853, 39_894),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresnel(z: arb, pi: arb, imaginary_unit: acb) -> acb:
    return (imaginary_unit * pi / 4).exp() / arb(2).sqrt() * (
        ((-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt() * z).erf()
    )


def ratio(mode: int, pi: arb, imaginary_unit: acb) -> acb:
    m = arb(mode)
    t = arb(T)
    x = 2 * pi * m**2 / (t + 2 * pi * m**2)
    q_a = (x / 2).sqrt() * (arb(A) - 2 * m / x)
    q_b = (x / 2).sqrt() * (arb(B) - 2 * m / x)
    i0 = fresnel(q_b, pi, imaginary_unit) - fresnel(q_a, pi, imaginary_unit)
    i1 = ((imaginary_unit * pi * q_b**2 / 2).exp() - (imaginary_unit * pi * q_a**2 / 2).exp()) / (imaginary_unit * pi)
    return i0 / (1 + imaginary_unit) + x.sqrt() * i1 / (m * arb(2).sqrt() * (1 + imaginary_unit))


def replay_segment(start: int, end: int, theta_zero: arb, theta_exact: arb, pi: arb, imaginary_unit: acb) -> tuple[acb, acb, acb, arb]:
    total = acb(0)
    endpoint = acb(0)
    phase = acb(0)
    triangle = arb(0)
    for mode in range(start, end + 1):
        m = arb(mode)
        r = ratio(mode, pi, imaginary_unit)
        c0 = (imaginary_unit * (theta_zero - arb(T) * m.log())).exp() / m.sqrt()
        ce = (imaginary_unit * (theta_exact - arb(T) * m.log())).exp() / m.sqrt()
        defect = c0 * r - ce
        total += defect
        endpoint += c0 * (r - 1)
        phase += c0 - ce
        triangle += abs(defect)
    return total, endpoint, phase, triangle


def stored_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def overlap_complex(left: acb, right: acb) -> bool:
    return left.real.overlaps(right.real) and left.imag.overlaps(right.imag)


def main() -> int:
    ctx.dps = 125
    ctx.threads = 1
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "aggregate barrier artifact is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    t = arb(T)
    theta_exact = acb(arb("0.25"), t / 2).lgamma().imag - t * pi.log() / 2
    theta_zero = t / 2 * ((t / (2 * pi)).log() - 1) - pi / 8
    stored_phase = artifact["certificate"]["phase"]
    require(theta_exact.overlaps(arb(stored_phase["theta_exact_ball"])), "exact theta replay mismatch")
    require(theta_zero.overlaps(arb(stored_phase["theta_zero_ball"])), "theta-zero replay mismatch")

    stored_rows = artifact["certificate"]["segments"]
    require([(row["name"], *row["range"]) for row in stored_rows] == list(SEGMENTS), "segment roster drift")
    grand_total = acb(0)
    grand_endpoint = acb(0)
    grand_phase = acb(0)
    grand_triangle = arb(0)
    for row, (name, start, end) in zip(stored_rows, SEGMENTS):
        total, endpoint, phase, triangle = replay_segment(start, end, theta_zero, theta_exact, pi, imaginary_unit)
        require(overlap_complex(total, stored_complex(row["complex_frozen_minus_exact_ball"])), f"total replay mismatch: {name}")
        require(overlap_complex(endpoint, stored_complex(row["complex_endpoint_current_only_ball"])), f"endpoint replay mismatch: {name}")
        require(overlap_complex(phase, stored_complex(row["complex_phase_only_ball"])), f"phase replay mismatch: {name}")
        require(triangle.overlaps(arb(row["termwise_complex_triangle_ball"])), f"triangle replay mismatch: {name}")
        grand_total += total
        grand_endpoint += endpoint
        grand_phase += phase
        grand_triangle += triangle

    aggregate = artifact["certificate"]["aggregate"]
    require(overlap_complex(grand_total, stored_complex(aggregate["complex_frozen_minus_exact_ball"])), "aggregate total mismatch")
    require(overlap_complex(grand_endpoint, stored_complex(aggregate["complex_endpoint_current_only_ball"])), "aggregate endpoint mismatch")
    require(overlap_complex(grand_phase, stored_complex(aggregate["complex_phase_only_ball"])), "aggregate phase mismatch")
    require(grand_triangle.overlaps(arb(aggregate["termwise_complex_triangle_ball"])), "aggregate triangle mismatch")
    require(2 * grand_total.real < arb("-0.096"), "aggregate obstruction replay failed")
    require(abs(2 * grand_phase.real) < arb("0.00000001"), "phase-only replay unexpectedly large")
    require(arb(2).sqrt() / pi * abs(grand_total) > arb("0.028"), "physical-scale replay failed")
    require(artifact["decision"]["endpoint_retained_nonfrozen_morse_method_rejected"] is False, "method overclaim")
    require(artifact["decision"]["complete_x_integrated_error_lower_bound_proved"] is False, "full-integral overclaim")
    require(artifact["decision"]["complete_T_upper_proved"] is False, "T_upper overclaim")

    print("checked frozen endpoint-current aggregate barrier over 39273 modes at 125 digits", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
