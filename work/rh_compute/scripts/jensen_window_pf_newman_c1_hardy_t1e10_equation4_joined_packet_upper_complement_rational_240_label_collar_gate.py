#!/usr/bin/env python3
"""Certify the rationally closed 240-label upper-complement collar."""

from __future__ import annotations

import hashlib
import json
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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_rational_240_label_collar_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
BASE_BUILDER = Path(base.__file__).resolve()
BASE_CHECKER = BASE_BUILDER.with_name("check_" + BASE_BUILDER.name)
BASE_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_endpoint_saddle_core_gate.json"

PRECISION_BITS = 384
PANELS = 32
HEIGHT = arb("10000000000")
LOWER = arb("621.5")
SPLIT = arb("621.5625")
Q_PLUS = arb("2561211.5")
STATIONARY_COUNT = 230
COLLAR_COUNT = 240
ROOT_ORDER = 16


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_certificate() -> dict[str, Any]:
    pi = arb.pi()
    p = HEIGHT / (2 * pi)
    q_s_lower = LOWER + p / LOWER
    q_s_split = SPLIT + p / SPLIT
    last_stationary = Q_PLUS + STATIONARY_COUNT - 1
    first_remote = Q_PLUS + COLLAR_COUNT
    last_stationary_root = base.lower_saddle(last_stationary, p)
    first_stationary_root = base.lower_saddle(Q_PLUS, p)

    require(COLLAR_COUNT % ROOT_ORDER == 0, "collar does not close at B")
    require(COLLAR_COUNT % 2 == 0, "collar does not close at L")
    require(COLLAR_COUNT - STATIONARY_COUNT == 10, "collar padding drift")
    require(last_stationary_root > LOWER, "last stationary root escaped below L")
    require(first_stationary_root < SPLIT, "first stationary root escaped above B")
    require(Q_PLUS - q_s_split > 0, "collar continuation is not nonstationary at B")
    require(first_remote - q_s_lower > 10, "remote tail launch gap is not above ten")

    integral, panels = base.integrate_panels(
        base.grouped_integrand_factory(
            t=HEIGHT,
            lower=LOWER,
            q0=Q_PLUS,
            count=COLLAR_COUNT,
        ),
        SPLIT - LOWER,
        PANELS,
    )
    require(integral.is_finite(), "240-label collar integral is not finite")

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "panels": PANELS,
        "t": HEIGHT.str(50, more=True),
        "p_equals_t_over_2pi_ball": base.arb_record(p),
        "lower_endpoint": LOWER.str(30, more=True),
        "rational_split": SPLIT.str(30, more=True),
        "split_width": (SPLIT - LOWER).str(30, more=True),
        "stationary_label_count": STATIONARY_COUNT,
        "collar_label_count": COLLAR_COUNT,
        "nonstationary_padding_labels": COLLAR_COUNT - STATIONARY_COUNT,
        "root_of_unity_order_at_split": ROOT_ORDER,
        "least_closed_collar_count": COLLAR_COUNT,
        "left_endpoint_zero_exact": "sum_(j=0)^239 (-1)^j=0",
        "right_endpoint_zero_exact": "sum_(j=0)^239 exp(2*pi*i*j*621.5625)=0",
        "first_stationary_root_ball": base.arb_record(first_stationary_root),
        "last_stationary_root_ball": base.arb_record(last_stationary_root),
        "collar_continuation_gap_at_split_ball": base.arb_record(Q_PLUS - q_s_split),
        "remote_tail_first_q": first_remote.str(30, more=True),
        "remote_tail_launch_gap_at_lower_ball": base.arb_record(first_remote - q_s_lower),
        "grouped_240_label_collar_ball": base.acb_record(integral),
        "panel_rows": panels,
        "exact_decomposition": (
            "T_upper=P_240+R_B+R_L, where P_240 sums q_+ through q_++239 on [L,B], "
            "R_B continues those 240 labels on [B,a], and R_L sums q>=q_++240 on [L,a]"
        ),
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Rationally closed 240-label upper-complement collar

Date: 2026-08-27

Status: the least rationally closed endpoint collar containing all 230 upper
saddles is rigorously enclosed; both remaining upper tails are nonstationary;
not a proof of the complete joined packet or RH.

At `L=621.5`, every even consecutive half-integer roster has zero grouped
amplitude.  At `B=L+1/16=621.5625`, the Fourier step has exact order 16.
Therefore the least multiple of 16 containing the 230 stationary labels is

```text
240=15*16.
```

The extended collar has exact zeros at both endpoints:

```text
sum_(j=0)^239 (-1)^j=0,
sum_(j=0)^239 exp(2*pi*i*j*621.5625)=0.            (RC1)
```

It contains the 230 stationary labels and ten additional nonstationary labels.
The grouped collar is

```text
P_240=sum_(q=q_+)^(q_++239)
      integral_L^B y^(-s)exp(-i*pi*y^2+2*pi*i*q*y)dy.           (RC2)
```

Thirty-two Arb panels at 384 bits give

```text
P_240 = {c['grouped_240_label_collar_ball']['real_ball']}
      + i {c['grouped_240_label_collar_ball']['imag_ball']},

|P_240| = {c['grouped_240_label_collar_ball']['absolute_ball']}.
```

The exact remainder split is

```text
T_upper=P_240+R_B+R_L,                            (RC3)

R_B=sum_(q=q_+)^(q_++239) integral_B^a (...),
R_L=sum_(q>=q_++240) integral_L^a (...).
```

Every phase in both remainders is nonstationary.  The minimum `R_B` gap at
`B` is `{c['collar_continuation_gap_at_split_ball']['ball']}`.  The minimum
`R_L` launch gap at `L` is
`{c['remote_tail_launch_gap_at_lower_ball']['ball']}`, replacing the previous
near-endpoint gap `0.320323...` by a gap greater than ten.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
half-integer Fourier spacing, or `p=t/(2*pi)`.  The numbers `1/16`, `16`, and
`240` come exactly from the rational endpoint and its root-of-unity order; no
fitted geometric constant is used.

Proof boundary: exact dual-endpoint cancellation, least closed collar length,
production gap census, and rigorous grouped complex enclosure of `P_240` only.
No quantitative enclosure of `R_B`, `R_L`, the lower complementary tail,
complete ordinary or joined packet, `J_Z`, or `D_K` is proved, nor is any
non-A, all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    resource_mode = base.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    for path in (BASE_RESULT, BASE_BUILDER, BASE_CHECKER, CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    require(load_json(BASE_RESULT).get("passed") is True, "230-label core gate not passed")

    certificate = build_certificate()
    artifact = {
        "kind": "rh_c1_hardy_equation4_upper_complement_rational_240_label_collar_gate",
        "date": "2026-08-27",
        "status": "dual_endpoint_zero_240_label_collar_complex_ball_certified",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "least_closed_collar_is_240_labels": True,
            "both_collar_endpoints_vanish_exactly": True,
            "all_230_stationary_labels_captured": True,
            "both_upper_remainders_nonstationary": True,
            "grouped_240_label_collar_enclosed": True,
            "upper_remainders_enclosed": False,
        },
        "sources": {
            relative(BASE_RESULT): file_hash(BASE_RESULT),
            relative(BASE_BUILDER): file_hash(BASE_BUILDER),
            relative(BASE_CHECKER): file_hash(BASE_CHECKER),
            relative(BUILDER): file_hash(BUILDER),
            relative(CHECKER): file_hash(CHECKER),
        },
        "proof_boundary": (
            "Exact dual-endpoint cancellation, least closed collar length, production gap "
            "census, and rigorous grouped complex enclosure of P_240 only. No R_B, R_L, "
            "lower complementary tail, complete ordinary or joined packet, J_Z, D_K, "
            "non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    base.atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    base.atomic_write(NOTE, note_text(artifact))
    print("certified rationally closed 240-label upper collar")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
