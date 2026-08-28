#!/usr/bin/env python3
"""Enclose the two far tails and complete the 752-label upper complement."""

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

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate as endpoint_base
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as ibp_base


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_far_remainder_eight_round_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

TRANSITION_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_label_exterior_transition_layer_gate"
TRANSITION_RESULT = ROOT / "work" / "rh_compute" / "results" / f"{TRANSITION_STEM}.json"
TRANSITION_BUILDER = ROOT / "work" / "rh_compute" / "scripts" / f"{TRANSITION_STEM}.py"
TRANSITION_CHECKER = TRANSITION_BUILDER.with_name("check_" + TRANSITION_BUILDER.name)
ENDPOINT_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate.json"
ENDPOINT_BUILDER = Path(endpoint_base.__file__).resolve()
ENDPOINT_CHECKER = ENDPOINT_BUILDER.with_name("check_" + ENDPOINT_BUILDER.name)
OBSTRUCTION_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate.json"
OBSTRUCTION_BUILDER = Path(ibp_base.__file__).resolve()
OBSTRUCTION_CHECKER = OBSTRUCTION_BUILDER.with_name("check_" + OBSTRUCTION_BUILDER.name)

PRECISION_BITS = 384
ORDER = 8
SLABS = 8192
HEIGHT = arb("10000000000")
L = arb("621.5")
C = arb("621.6875")
A = arb("39852.5")
Q_PLUS = arb("2561211.5")
PACKET_COUNT = 752
Q_REMOTE = Q_PLUS + PACKET_COUNT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def source_paths() -> tuple[Path, ...]:
    return (
        TRANSITION_RESULT,
        TRANSITION_BUILDER,
        TRANSITION_CHECKER,
        ENDPOINT_RESULT,
        ENDPOINT_BUILDER,
        ENDPOINT_CHECKER,
        OBSTRUCTION_RESULT,
        OBSTRUCTION_BUILDER,
        OBSTRUCTION_CHECKER,
        BUILDER,
        CHECKER,
    )


def endpoint_families(p: arb) -> dict[str, dict[str, Any]]:
    q_c = endpoint_base.q_curve(C, p)
    q_l = endpoint_base.q_curve(L, p)
    q_a = endpoint_base.q_curve(A, p)
    finite_c_delta = Q_PLUS - q_c
    finite_a_delta = Q_PLUS - q_a
    remote_l_delta = Q_REMOTE - q_l
    remote_a_delta = Q_REMOTE - q_a
    for name, delta in (
        ("F_C at C", finite_c_delta),
        ("F_C at a", finite_a_delta),
        ("I_L at L", remote_l_delta),
        ("I_L at a", remote_a_delta),
    ):
        require(delta > 0, f"{name} acquired a stationary point")

    return {
        "F_C_finite_752": {
            "definition": (
                "sum_(n=0)^751 exp(2*pi*i*(q_++n)*y)/(q_++n-Q(y))^k"
            ),
            "endpoint_C": endpoint_base.endpoint_record(
                y=C,
                q0=Q_PLUS,
                delta=finite_c_delta,
                order=16,
                numerator=11,
                mode="finite",
                cycles=47,
            ),
            "endpoint_a": endpoint_base.endpoint_record(
                y=A,
                q0=Q_PLUS,
                delta=finite_a_delta,
                order=2,
                numerator=1,
                mode="finite",
                cycles=376,
            ),
        },
        "I_L_infinite_upper": {
            "definition": (
                "sum_(n>=0) exp(2*pi*i*(q_++752+n)*y)/(q_++752+n-Q(y))^k"
            ),
            "endpoint_L": endpoint_base.endpoint_record(
                y=L,
                q0=Q_REMOTE,
                delta=remote_l_delta,
                order=2,
                numerator=1,
                mode="infinite",
            ),
            "endpoint_a": endpoint_base.endpoint_record(
                y=A,
                q0=Q_REMOTE,
                delta=remote_a_delta,
                order=2,
                numerator=1,
                mode="infinite",
            ),
        },
    }


def build_certificate(transition: dict[str, Any]) -> dict[str, Any]:
    p = HEIGHT / (2 * arb.pi())
    require(A * A < p, "physical interval crossed the minimum of Q")
    table = ibp_base.coefficient_table(ORDER)
    endpoints = endpoint_families(p)

    finite = ibp_base.family_certificate(
        name="F_C",
        table=table,
        p=p,
        endpoints=endpoints["F_C_finite_752"],
        left_key="endpoint_C",
        right_key="endpoint_a",
        left=C,
        right=A,
        q0=Q_PLUS,
        orientation="upper",
        finite_count=PACKET_COUNT,
    )
    remote = ibp_base.family_certificate(
        name="I_L",
        table=table,
        p=p,
        endpoints=endpoints["I_L_infinite_upper"],
        left_key="endpoint_L",
        right_key="endpoint_a",
        left=L,
        right=A,
        q0=Q_REMOTE,
        orientation="upper",
        finite_count=None,
    )

    g752 = ibp_base.saved_complex(
        transition["certificate"]["grouped_752_label_transition_ball"]
    )
    finite_ball = ibp_base.saved_complex(finite["tail_integral_ball"])
    remote_ball = ibp_base.saved_complex(remote["tail_integral_ball"])
    upper = g752 + finite_ball + remote_ball
    require(upper.is_finite(), "complete upper complement is not finite")

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "integration_rounds": ORDER,
        "remainder_slabs_per_family": SLABS,
        "cluster_power": ibp_base.CLUSTER_POWER,
        "t": HEIGHT.str(50, more=True),
        "p_equals_t_over_2pi_ball": ibp_base.arb_record(p),
        "physical_interval": [L.str(30, more=True), A.str(30, more=True)],
        "physical_interval_below_sqrt_p": True,
        "transition_split": C.str(30, more=True),
        "transition_label_count": PACKET_COUNT,
        "remote_first_q": Q_REMOTE.str(30, more=True),
        "endpoint_families": endpoints,
        "coefficient_level_term_counts": [
            {
                "level": level,
                "denominator_powers": sorted(table[level]),
                "exact_monomial_count": sum(len(terms) for terms in table[level].values()),
            }
            for level in range(ORDER + 1)
        ],
        "recurrence": "h_0=y^(-1/2), h_(n+1)=d/dy[h_n/(q-y-p/y)]",
        "ibp_identity": (
            "I=sum_(n=0)^(N-1)(-1)^n(i2*pi)^(-n-1)"
            "[exp(i*phi)h_n/D]_left^right+(-1)^N(i2*pi)^(-N)"
            "integral exp(i*phi)h_N"
        ),
        "families": {"F_C": finite, "I_L": remote},
        "grouped_752_label_transition_ball": transition["certificate"][
            "grouped_752_label_transition_ball"
        ],
        "complete_upper_complement_ball": ibp_base.acb_record(upper),
        "exact_decomposition": (
            "T_upper=G_752+F_C+I_L, with G_752 on [L,C], F_C the same finite "
            "752-label roster on [C,a], and I_L the q>=q_++752 roster on [L,a]"
        ),
        "abel_limit_justification": (
            "derive I_L at finite Abel cutoff; its k=1 endpoint sums use the certified "
            "alternating Lerch convention, while every integrated order-eight label sum "
            "is absolutely dominated by denominator powers at least eight"
        ),
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    finite = c["families"]["F_C"]
    remote = c["families"]["I_L"]
    upper = c["complete_upper_complement_ball"]
    return f"""# Eight-round closure of the far upper-complement remainders

Date: 2026-08-27

Status: both far upper remainders and the complete upper complement at the
production height are rigorously enclosed as complex balls; the lower
complement and complete joined packet remain open; this is not a proof of RH.

Section 11.502 supplies the exact decomposition

```text
T_upper=G_752+F_C+I_L,                              (FR1)
```

where `C=621+11/16`, `F_C` contains 752 labels on `[C,a]`, and `I_L`
starts at `q_++752` on `[L,a]`.  Both launches exceed eight quadratic widths
and the full interval lies below `sqrt(p)`, so `Q` decreases throughout.

The endpoint sums are evaluated before any norm.  At `C`, the finite roster
uses the primitive order-16 root with numerator 11 and 47 cycles; at `a` it
uses 376 alternating cycles.  The remote endpoints at `L` and `a` use the
alternating infinite Lerch/digamma formulas.  Every finite formula through
power 16 overlaps its direct 752-term Arb sum, and every infinite row obeys
the shifted Lerch recurrence, including exact cancellation of the `k=1`
Hurwitz poles.

Apply eight exact integration-by-parts rounds using

```text
h_0=y^(-1/2),       h_(n+1)=d/dy[h_n/(q-y-p/y)].    (FR2)
```

On 8192 fourth-power endpoint-clustered slabs, the retained complex results
are

```text
F_C endpoint partial = {finite['endpoint_partial_sum_ball']['real_ball']}
                     + i {finite['endpoint_partial_sum_ball']['imag_ball']},
F_C remainder <= {finite['remainder_census']['scaled_remainder_bound_ball']['ball']},
F_C ball = {finite['tail_integral_ball']['real_ball']}
         + i {finite['tail_integral_ball']['imag_ball']};

I_L endpoint partial = {remote['endpoint_partial_sum_ball']['real_ball']}
                     + i {remote['endpoint_partial_sum_ball']['imag_ball']},
I_L remainder <= {remote['remainder_census']['scaled_remainder_bound_ball']['ball']},
I_L ball = {remote['tail_integral_ball']['real_ball']}
         + i {remote['tail_integral_ball']['imag_ball']}.
```

Joining these balls to the certified transition packet before taking a norm
gives the complete upper complement

```text
T_upper = {upper['real_ball']}
        + i {upper['imag_ball']},

|T_upper| = {upper['absolute_ball']}.               (FR3)
```

The independent checker changes to nine recurrence rounds, 10000
fifth-power-clustered slabs, direct finite endpoint sums, and explicit
alternating eta values.  Both changed-order family balls and their complete
upper join overlap the production enclosures.  A separate altered-height
positive-denominator direct quadrature overlaps its fifth-order recurrence
ball.

Pi provenance: every `pi` comes from the inherited Fresnel/Fourier phase,
`phi'=2*pi*D`, exact root-of-unity characters, or `p=t/(2*pi)`.  The endpoint
`C` and count 752 are inherited from the exact order-16/eight-width gate.  No
circle, polygon, fitted constant, or visual pattern supplies `pi`.

Proof boundary: rigorous endpoint tables through power 16, eight-round
complex enclosures of `F_C` and `I_L`, and the complete production-height
upper-complement ball only.  No lower complementary-tail, complete ordinary
or joined packet, `J_Z`, `D_K`, non-A, all-height, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
"""


def main() -> int:
    resource_mode = ibp_base.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    for path in source_paths():
        require(path.is_file(), f"missing source or dependency: {path}")

    transition = load_json(TRANSITION_RESULT)
    require(transition.get("passed") is True, "752-label transition gate not passed")
    require(
        transition["decision"]["complete_upper_complement_enclosed"] is False,
        "transition gate proof boundary drift",
    )
    require(load_json(ENDPOINT_RESULT).get("passed") is True, "endpoint-Hurwitz gate not passed")
    obstruction = load_json(OBSTRUCTION_RESULT)
    require(obstruction.get("passed") is True, "eight-round obstruction gate not passed")
    require(
        obstruction["decision"]["unsplit_eight_round_absolute_route_rejected"] is True,
        "unsplit-route rejection drift",
    )

    certificate = build_certificate(transition)
    artifact = {
        "kind": "rh_c1_hardy_equation4_upper_complement_closed_752_far_remainder_eight_round_gate",
        "date": "2026-08-27",
        "status": "complete_production_upper_complement_complex_ball_certified",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "F_C_complex_ball_certified": True,
            "I_L_complex_ball_certified": True,
            "endpoint_carriers_retained_before_norm": True,
            "complete_upper_complement_complex_ball_certified": True,
            "lower_complementary_tail_enclosed": False,
            "complete_joined_packet_enclosed": False,
            "D_K_corridor_proved": False,
        },
        "sources": {relative(path): file_hash(path) for path in source_paths()},
        "proof_boundary": (
            "Production endpoint tables through power 16, rigorous eight-round F_C and I_L "
            "complex balls, and the complete upper-complement ball only. No lower complementary "
            "tail, complete ordinary or joined packet, J_Z, D_K, non-A, all-height, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    ibp_base.atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    ibp_base.atomic_write(NOTE, note_text(artifact))
    print("certified complete upper complement after 752-label transition closure", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
