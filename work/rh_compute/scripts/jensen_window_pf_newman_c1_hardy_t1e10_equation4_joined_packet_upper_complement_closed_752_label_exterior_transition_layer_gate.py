#!/usr/bin/env python3
"""Certify the closed 752-label upper exterior-transition regrouping."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_endpoint_saddle_core_gate as base


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_label_exterior_transition_layer_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

BASE_BUILDER = Path(base.__file__).resolve()
BASE_CHECKER = BASE_BUILDER.with_name("check_" + BASE_BUILDER.name)
COLLAR_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_rational_240_label_collar_gate.json"
ENDPOINT_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate.json"
ENDPOINT_BUILDER = ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate.py"
ENDPOINT_CHECKER = ENDPOINT_BUILDER.with_name("check_" + ENDPOINT_BUILDER.name)
OBSTRUCTION_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate.json"
OBSTRUCTION_BUILDER = ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate.py"
OBSTRUCTION_CHECKER = OBSTRUCTION_BUILDER.with_name("check_" + OBSTRUCTION_BUILDER.name)

PRECISION_BITS = 384
PANELS = 192
HEIGHT = arb("10000000000")
LOWER = arb("621.5")
SPLIT = arb("621.6875")
PREVIOUS_SPLIT = arb("621.625")
UPPER = arb("39852.5")
Q_PLUS = arb("2561211.5")
STATIONARY_COUNT = 230
PACKET_COUNT = 752
PREVIOUS_PACKET_COUNT = 736
ROOT_ORDER = 16
CHI_TARGET_SQUARED = 64


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def q_value(y: arb, p: arb) -> arb:
    return y + p / y


def q_prime_absolute(y: arb, p: arb) -> arb:
    derivative = 1 - p / (y * y)
    require(derivative < 0, "endpoint escaped the lower saddle branch")
    return -derivative


def scaled_detuning(delta: arb, q_prime_abs: arb) -> arb:
    require(delta > 0, "detuning must be positive")
    return delta / q_prime_abs.sqrt()


def source_paths() -> tuple[Path, ...]:
    return (
        COLLAR_RESULT,
        ENDPOINT_RESULT,
        ENDPOINT_BUILDER,
        ENDPOINT_CHECKER,
        OBSTRUCTION_RESULT,
        OBSTRUCTION_BUILDER,
        OBSTRUCTION_CHECKER,
        BASE_BUILDER,
        BASE_CHECKER,
        BUILDER,
        CHECKER,
    )


def build_certificate() -> dict[str, Any]:
    pi = arb.pi()
    p = HEIGHT / (2 * pi)

    require(math.gcd(11, ROOT_ORDER) == 1, "split Fourier step is not primitive order 16")
    require(PACKET_COUNT == 47 * ROOT_ORDER, "packet factorization drift")
    require(PACKET_COUNT % ROOT_ORDER == 0, "packet does not close at C")
    require(PACKET_COUNT % 2 == 0, "packet does not close at L")
    require(PREVIOUS_PACKET_COUNT == PACKET_COUNT - ROOT_ORDER, "previous packet drift")
    require(UPPER * UPPER < p, "physical interval crosses the minimum of Q")

    first_root = base.lower_saddle(Q_PLUS, p)
    last_stationary = Q_PLUS + STATIONARY_COUNT - 1
    last_root = base.lower_saddle(last_stationary, p)
    require(last_root > LOWER, "last stationary root escaped below L")
    require(first_root < SPLIT, "first stationary root escaped above C")

    q_lower = q_value(LOWER, p)
    q_split = q_value(SPLIT, p)
    q_previous_split = q_value(PREVIOUS_SPLIT, p)
    split_derivative_abs = q_prime_absolute(SPLIT, p)
    previous_split_derivative_abs = q_prime_absolute(PREVIOUS_SPLIT, p)
    lower_derivative_abs = q_prime_absolute(LOWER, p)

    finite_detuning = Q_PLUS - q_split
    previous_finite_detuning = Q_PLUS - q_previous_split
    remote_first_q = Q_PLUS + PACKET_COUNT
    previous_remote_first_q = Q_PLUS + PREVIOUS_PACKET_COUNT
    remote_detuning = remote_first_q - q_lower
    previous_remote_detuning = previous_remote_first_q - q_lower

    finite_margin = finite_detuning**2 - CHI_TARGET_SQUARED * split_derivative_abs
    previous_finite_margin = (
        previous_finite_detuning**2
        - CHI_TARGET_SQUARED * previous_split_derivative_abs
    )
    remote_margin = remote_detuning**2 - CHI_TARGET_SQUARED * lower_derivative_abs
    previous_remote_margin = (
        previous_remote_detuning**2
        - CHI_TARGET_SQUARED * lower_derivative_abs
    )

    require(finite_detuning > 0, "finite continuation is not nonstationary at C")
    require(previous_finite_detuning > 0, "previous grid point is not yet nonstationary")
    require(remote_detuning > 0, "remote tail is not nonstationary at L")
    require(previous_remote_detuning > 0, "previous remote packet is not yet nonstationary")
    require(finite_margin > 0, "finite continuation does not exceed eight widths")
    require(previous_finite_margin < 0, "previous 1/16-grid split already exceeds eight widths")
    require(remote_margin > 0, "remote tail does not exceed eight widths")
    require(previous_remote_margin < 0, "previous multiple-of-16 packet already exceeds eight widths")

    width = SPLIT - LOWER
    integral, panels = base.integrate_panels(
        base.grouped_integrand_factory(
            t=HEIGHT,
            lower=LOWER,
            q0=Q_PLUS,
            count=PACKET_COUNT,
        ),
        width,
        PANELS,
    )
    require(integral.is_finite(), "752-label grouped transition integral is not finite")

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "panels": PANELS,
        "t": HEIGHT.str(50, more=True),
        "p_equals_t_over_2pi_ball": base.arb_record(p),
        "lower_endpoint": LOWER.str(30, more=True),
        "transition_split": SPLIT.str(30, more=True),
        "physical_upper_endpoint": UPPER.str(30, more=True),
        "physical_interval_below_sqrt_p": True,
        "transition_width": width.str(30, more=True),
        "previous_grid_split": PREVIOUS_SPLIT.str(30, more=True),
        "stationary_label_count": STATIONARY_COUNT,
        "transition_label_count": PACKET_COUNT,
        "previous_multiple_of_16_count": PREVIOUS_PACKET_COUNT,
        "root_of_unity_order_at_split": ROOT_ORDER,
        "packet_factorization_exact": "752=47*16",
        "left_endpoint_zero_exact": "sum_(j=0)^751 (-1)^j=0",
        "right_endpoint_zero_exact": (
            "sum_(j=0)^751 exp(2*pi*i*j*(621+11/16))=0; gcd(11,16)=1"
        ),
        "first_stationary_root_ball": base.arb_record(first_root),
        "last_stationary_root_ball": base.arb_record(last_root),
        "all_230_stationary_roots_in_transition_interval": True,
        "q_at_lower_ball": base.arb_record(q_lower),
        "q_at_split_ball": base.arb_record(q_split),
        "q_at_previous_split_ball": base.arb_record(q_previous_split),
        "finite_continuation_detuning_at_split_ball": base.arb_record(finite_detuning),
        "finite_continuation_abs_qprime_at_split_ball": base.arb_record(split_derivative_abs),
        "finite_continuation_scaled_detuning_ball": base.arb_record(
            scaled_detuning(finite_detuning, split_derivative_abs)
        ),
        "finite_continuation_eight_width_margin_ball": base.arb_record(finite_margin),
        "previous_grid_eight_width_margin_ball": base.arb_record(previous_finite_margin),
        "remote_tail_first_q": remote_first_q.str(30, more=True),
        "remote_tail_detuning_at_lower_ball": base.arb_record(remote_detuning),
        "remote_tail_abs_qprime_at_lower_ball": base.arb_record(lower_derivative_abs),
        "remote_tail_scaled_detuning_ball": base.arb_record(
            scaled_detuning(remote_detuning, lower_derivative_abs)
        ),
        "remote_tail_eight_width_margin_ball": base.arb_record(remote_margin),
        "previous_multiple_eight_width_margin_ball": base.arb_record(previous_remote_margin),
        "grouped_752_label_transition_ball": base.acb_record(integral),
        "panel_rows": panels,
        "exact_decomposition": (
            "T_upper=G_752+F_C+I_L, where G_752 sums q_+ through q_++751 on [L,C], "
            "F_C continues those 752 labels on [C,a], and I_L sums q>=q_++752 on [L,a]"
        ),
        "eight_width_criterion": (
            "chi=delta/sqrt(abs(Q'(y))); chi>8 is certified by "
            "delta^2-64*abs(Q'(y))>0"
        ),
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Closed 752-label exterior transition layer

Date: 2026-08-27

Status: the exact dual-endpoint regrouping, its rigorous transition-packet
complex ball, and both eight-width far-launch tests are certified; the two far
integrals remain open; this is not a proof of the complete joined packet or RH.

The failed unsplit eight-round integration-by-parts route begins too close to
the endpoint saddles.  Replace that split by

```text
C=L+3/16=621.6875=621+11/16,
N=752=47*16.
```

At `L=621.5` the Fourier step is `-1`.  At `C` it is the primitive order-16
root `exp(2*pi*i*11/16)`, because `gcd(11,16)=1`.  Thus the same roster
vanishes exactly at both endpoints:

```text
sum_(j=0)^751 (-1)^j=0,
sum_(j=0)^751 exp(2*pi*i*j*C)=0.                 (ET1)
```

All 230 stationary labels lie in `[L,C]`; their extreme lower roots are
`{c['last_stationary_root_ball']['ball']}` and
`{c['first_stationary_root_ball']['ball']}`.  The packet

```text
G_752=sum_(q=q_+)^(q_++751)
      integral_L^C y^(-s)exp(-i*pi*y^2+2*pi*i*q*y)dy             (ET2)
```

is evaluated with one grouped geometric amplitude.  At 384 bits on 192 Arb
panels,

```text
G_752 = {c['grouped_752_label_transition_ball']['real_ball']}
      + i {c['grouped_752_label_transition_ball']['imag_ball']},

|G_752| = {c['grouped_752_label_transition_ball']['absolute_ball']}.
```

For `Q(y)=y+p/y`, `p=t/(2*pi)`, define

```text
chi=delta/sqrt(abs(Q'(y))),
delta=q-Q(y).
```

The comparison `chi>8` is made without decimal threshold fitting by proving
`delta^2-64*abs(Q'(y))>0`.  For the finite continuation at `C`, the certified
scaled detuning is
`{c['finite_continuation_scaled_detuning_ball']['ball']}` and its squared
margin is `{c['finite_continuation_eight_width_margin_ball']['ball']}`.
The preceding 1/16-grid point `L+2/16` fails, with margin
`{c['previous_grid_eight_width_margin_ball']['ball']}`.

For the remote tail beginning at `q_++752` from `L`, the scaled detuning is
`{c['remote_tail_scaled_detuning_ball']['ball']}` and its squared margin is
`{c['remote_tail_eight_width_margin_ball']['ball']}`.  The preceding multiple
of 16, `736`, fails, with margin
`{c['previous_multiple_eight_width_margin_ball']['ball']}`.  Monotonicity of
the positive detuning on the lower saddle branch makes these the first
1/16-grid split and the least multiple-of-16 roster meeting the threshold.

The exact upper-complement decomposition is therefore

```text
T_upper=G_752+F_C+I_L,                              (ET3)

F_C=sum_(q=q_+)^(q_++751) integral_C^a (...),
I_L=sum_(q>=q_++752)       integral_L^a (...).
```

Both retained pieces are genuinely far nonstationary candidates for the
already-certified endpoint-Hurwitz algebra and repeated integration by parts.
They are not bounded by this gate.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
Fourier characters, or `p=t/(2*pi)`.  The split `3/16` comes from the first
rational Fourier grid point meeting the stated squared-detuning inequality;
`752` is the least multiple of its exact root order 16 meeting the remote
inequality.  No circle, polygon, fitted constant, or visual pattern supplies
`pi` here.

Proof boundary: exact dual-endpoint cancellation, the exact decomposition,
stationary-root census, rigorous `G_752` complex enclosure, and the two
eight-width launch inequalities only.  No quantitative enclosure of `F_C`,
`I_L`, the lower complementary tail, complete ordinary or joined packet,
`J_Z`, or `D_K` is proved, nor is any non-A, all-height, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    resource_mode = base.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1

    for path in source_paths():
        require(path.is_file(), f"missing dependency: {path}")
    require(load_json(COLLAR_RESULT).get("passed") is True, "240-label collar gate not passed")
    endpoint = load_json(ENDPOINT_RESULT)
    require(endpoint.get("passed") is True, "endpoint-Hurwitz gate not passed")
    require(
        endpoint.get("status")
        == "production_nonstationary_complement_endpoint_sums_through_power_16_certified",
        "endpoint-Hurwitz status drift",
    )
    obstruction = load_json(OBSTRUCTION_RESULT)
    require(obstruction.get("passed") is True, "eight-round obstruction gate not passed")
    require(
        obstruction.get("status")
        == "three_tail_balls_certified_unsplit_eight_round_absolute_route_rejected",
        "eight-round obstruction status drift",
    )

    certificate = build_certificate()
    artifact = {
        "kind": "rh_c1_hardy_equation4_upper_complement_closed_752_label_exterior_transition_layer_gate",
        "date": "2026-08-27",
        "status": "exact_closed_752_label_exterior_transition_layer_and_far_separation_certified",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "unsplit_eight_round_route_remains_rejected": True,
            "all_230_stationary_labels_captured": True,
            "both_transition_endpoints_vanish_exactly": True,
            "grouped_752_label_transition_enclosed": True,
            "both_far_pieces_exceed_eight_scaled_widths": True,
            "finite_continuation_enclosed": False,
            "remote_tail_enclosed": False,
            "complete_upper_complement_enclosed": False,
        },
        "sources": {relative(path): file_hash(path) for path in source_paths()},
        "proof_boundary": (
            "Exact dual-endpoint cancellation and regrouping, stationary-root census, a rigorous "
            "complex enclosure of G_752, and two eight-width launch inequalities only. No F_C, "
            "I_L, lower complementary tail, complete ordinary or joined packet, J_Z, D_K, non-A, "
            "all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    base.atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    base.atomic_write(NOTE, note_text(artifact))
    print("certified closed 752-label exterior transition layer")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
