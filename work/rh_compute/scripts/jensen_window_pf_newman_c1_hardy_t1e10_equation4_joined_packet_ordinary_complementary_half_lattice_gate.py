#!/usr/bin/env python3
"""Reduce the ordinary Gamma-subtracted packet to complementary half-lattice tails."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

LOWER_CELL = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_cell_reversed_roster_endpoint_gate.json"
COMMON_CARRIER = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_common_carrier_subtraction_reduction_gate.json"
TRANSITION_PACKET = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_alternating_boundary_jet_packet_gate.json"
EXACT_H = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_continuation_defect_target_gate.json"
GAMMA_BULK = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json"

HEIGHT = 10_000_000_000
LOWER_BOUNDARY = arb("621.5")
UPPER_BOUNDARY = arb("39852.5")
Q_MINUS = arb("79788.5")
Q_PLUS = arb("2561211.5")
ORDINARY_FIRST = 622
ORDINARY_LAST = 39852
SOURCE_COUNT = 2_481_423
PRECISION_BITS = 384


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def arb_record(value: arb, digits: int = 80) -> dict[str, Any]:
    return {
        "ball": value.str(digits, more=True),
        "lower_float": float(value.lower()),
        "mid_float": float(value.mid()),
        "upper_float": float(value.upper()),
    }


def complex_from_record(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def lower_stationary_root(q: arb, p: arb) -> arb:
    discriminant = q**2 - 4 * p
    require(discriminant.lower() > 0, "stationary root is not real")
    root = discriminant.sqrt()
    return 2 * p / (q + root)


def build_certificate() -> dict[str, Any]:
    pi = arb.pi()
    t = arb(HEIGHT)
    p = t / (2 * pi)
    central = 2 * p.sqrt()
    q_s_lower = LOWER_BOUNDARY + p / LOWER_BOUNDARY
    q_s_upper = UPPER_BOUNDARY + p / UPPER_BOUNDARY
    upper_stationary_span = q_s_lower - Q_PLUS
    upper_first_nonstationary_gap = Q_PLUS + 230 - q_s_lower
    lower_nearest_gap = q_s_upper - (Q_MINUS - 1)

    require(Q_PLUS - Q_MINUS == SOURCE_COUNT, "finite roster endpoint drift")
    require(arb(229) < upper_stationary_span and upper_stationary_span < arb(230), "upper stationary count drift")
    require(upper_first_nonstationary_gap.lower() > 0, "upper nonstationary tail gap lost")
    require(lower_nearest_gap.lower() > 0, "lower complementary tail acquired a saddle")
    require(Q_MINUS - 1 < central and central < Q_MINUS, "fold placement drift")
    require(q_s_upper < Q_MINUS and Q_MINUS - 1 < q_s_upper, "upper boundary fold placement drift")

    first_stationary_q = Q_PLUS
    last_stationary_q = Q_PLUS + 229
    first_stationary_y = lower_stationary_root(first_stationary_q, p)
    last_stationary_y = lower_stationary_root(last_stationary_q, p)
    require(LOWER_BOUNDARY < last_stationary_y, "last stationary label fell below L")
    require(first_stationary_y < arb("622.5"), "stationary packet escaped the first ordinary cell")

    # Exact C_G=1/sqrt(1+exp(-2*pi*t)); retain the exponentially small
    # correction symbolically because its size lies far below working precision.
    gamma_x_log10 = -2 * pi * t / arb(10).log()
    gamma_delta_log10_upper = gamma_x_log10 - arb(2).log() / arb(10).log()
    mode_absolute_bound = 2 * (arb(ORDINARY_LAST).sqrt() - arb(ORDINARY_FIRST - 1).sqrt())
    gamma_packet_log10_upper = gamma_delta_log10_upper + mode_absolute_bound.log() / arb(10).log()

    transition = load_json(TRANSITION_PACKET)
    saved_cg = complex_from_record(transition["joined_packet"]["common_Gamma_coefficient_ball"])
    require(saved_cg.real.contains(1) and saved_cg.imag.contains(0), "saved common Gamma ball misses exact reduction")

    return {
        "precision_bits": ctx.prec,
        "p_equals_t_over_2pi_ball": arb_record(p),
        "central_fold_q_equals_2sqrtp_ball": arb_record(central),
        "q_s_at_L_ball": arb_record(q_s_lower),
        "q_s_at_a_ball": arb_record(q_s_upper),
        "upper_stationary_span_ball": arb_record(upper_stationary_span),
        "upper_stationary_label_count": 230,
        "upper_stationary_q_range": ["2561211.5", "2561440.5"],
        "upper_first_nonstationary_q": "2561441.5",
        "upper_first_nonstationary_gap_ball": arb_record(upper_first_nonstationary_gap),
        "lower_nearest_omitted_q": "79787.5",
        "lower_nearest_gap_ball": arb_record(lower_nearest_gap),
        "stationary_y_range": {
            "q_plus_root_ball": arb_record(first_stationary_y),
            "q_plus_plus_229_root_ball": arb_record(last_stationary_y),
            "containing_cell": ["621.5", "622.5"],
        },
        "common_Gamma_coefficient_exact": "C_G(t)=(1+exp(-2*pi*t))^(-1/2)",
        "Gamma_coefficient_defect_bound": "0<1-C_G<exp(-2*pi*t)/2",
        "Gamma_coefficient_defect_log10_upper_ball": arb_record(gamma_delta_log10_upper),
        "ordinary_mode_absolute_sum_bound_ball": arb_record(mode_absolute_bound),
        "ordinary_Gamma_defect_packet_log10_upper_ball": arb_record(gamma_packet_log10_upper),
        "ordinary_packet_identity": (
            "O_join=-T_lower-T_upper+(1-C_G)*sum_(m=622)^39852 m^(-s), where "
            "T_lower=sum_(q in Z+1/2,q<=79787.5) integral_L^a y^(-s)e^(-i*pi*y^2+2*pi*i*q*y)dy "
            "and T_upper=sum_(q in Z+1/2,q>=2561211.5) of the same integral, in the Abel sense"
        ),
        "passed": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Complementary half-lattice reduction of the ordinary packet

Date: 2026-08-27

Status: exact ordinary-packet complement and stationary-label census certified;
quantitative complementary-tail enclosure remains open.

For `622<=m<=39852`, put

```text
w_m(u)=(1+u/m)^(-s)exp(-i*pi*u^2),
beta_m=integral_(-1/2)^(1/2) w_m(u)
       sum_(q=q_-)^(q_+-1) exp[2*pi*i(q-m)u]du,       (CL1)

q_-=79788.5,       q_+=2561211.5.
```

The complete half-integer Fourier lattice has the Abel/Poisson identity

```text
sum_(q in Z+1/2) exp[2*pi*i(q-m)u]
 =sum_(k in Z)(-1)^k delta(u-k).                    (CL2)
```

Only `k=0` lies in the cell, and `w_m(0)=1`.  Consequently the complete
lattice contributes exactly one, so the finite coefficient defect is exactly
the negative of the two missing half-lattice tails.  Reassembling the finitely
many cells before taking the Abel limit gives

```text
O_join=sum_(m=622)^39852 m^(-s)[beta_m-C_G]
      =-T_lower-T_upper
       +(1-C_G)sum_(m=622)^39852 m^(-s),             (CL3)

T_lower=sum_(q in Z+1/2, q<=79787.5)
        integral_L^a y^(-s)exp(-i*pi*y^2+2*pi*i*q*y)dy,

T_upper=sum_(q in Z+1/2, q>=2561211.5)
        integral_L^a y^(-s)exp(-i*pi*y^2+2*pi*i*q*y)dy,

L=621.5,       a=39852.5.                            (CL4)
```

No source label is duplicated: the finite roster, lower complement, and upper
complement partition `Z+1/2` exactly.

The Gamma factor also simplifies exactly.  The phase of `B=I_bulk/I_G` is
`theta-theta_0`; combining its modulus with the exact `H(t)` normalization
gives

```text
C_G(t)=[1+exp(-2*pi*t)]^(-1/2),
0<1-C_G<exp(-2*pi*t)/2.                              (CL5)
```

At `t=10^10`, the entire final term in (CL3) has

```text
log10 absolute upper bound
 < {c['ordinary_Gamma_defect_packet_log10_upper_ball']['ball']}.
```

It is retained exactly; it is not silently replaced by zero.

For one complementary label the phase is

```text
Phi_q(y)=-t log(y)-pi*y^2+2*pi*q*y,
Phi_q'(y)=2*pi[q-y-p/y],       p=t/(2*pi).           (CL6)
```

Since `y+p/y` decreases on `[L,a]`, the lower complement has no stationary
label.  Its nearest gap is

```text
{c['lower_nearest_gap_ball']['ball']}.
```

The upper complement has exactly 230 stationary labels,

```text
q=2561211.5, 2561212.5, ..., 2561440.5,             (CL7)
```

because

```text
y+p/y at L minus q_+
 = {c['upper_stationary_span_ball']['ball']}.
```

Their lower saddle roots all lie in the first ordinary cell:

```text
{c['stationary_y_range']['q_plus_plus_229_root_ball']['ball']}
 <= y_q <=
{c['stationary_y_range']['q_plus_root_ball']['ball']}.
```

The next upper label is nonstationary with endpoint gap
`{c['upper_first_nonstationary_gap_ball']['ball']}`.  Thus the old 39,231-cell
ordinary wall is reduced to one 230-label endpoint-saddle packet plus two
nonstationary complementary tails and the explicit negligible Gamma defect.

Pi provenance: each `pi` comes from the inherited Fresnel phase, the exact
half-integer Fourier spacing, or `p=t/(2*pi)`.  No fitted geometric constant is
introduced.

Proof boundary: exact Abel/Poisson half-lattice completion, exact
complementary-tail identity, exact Gamma coefficient reduction and bound, and
production stationary-label census only.  No quantitative enclosure of the
230-label packet or nonstationary tails, complete ordinary or joined packet,
`J_Z`, or `D_K` is proved, nor is any non-A, all-height, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    dependencies = {
        "lower_cell": LOWER_CELL,
        "common_carrier": COMMON_CARRIER,
        "transition_packet": TRANSITION_PACKET,
        "exact_H": EXACT_H,
        "gamma_bulk": GAMMA_BULK,
    }
    for path in dependencies.values():
        require(path.is_file(), f"missing dependency: {path}")
        require(load_json(path).get("passed", True) is True, f"dependency not passed: {path}")
    require(CHECKER.is_file(), f"missing checker: {CHECKER}")

    certificate = build_certificate()
    artifact: dict[str, Any] = {
        "kind": "rh_c1_hardy_equation4_joined_packet_ordinary_complementary_half_lattice_gate",
        "date": "2026-08-27",
        "status": "ordinary_packet_exact_complementary_half_lattice_and_230_label_census_certified",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "complete_half_lattice_cell_value_is_one": True,
            "finite_defect_equals_missing_tail_pair": True,
            "Gamma_defect_retained_explicitly": True,
            "lower_complement_nonstationary": True,
            "upper_stationary_labels_exactly_230": True,
            "all_230_stationary_labels_in_first_ordinary_cell": True,
            "ordinary_packet_numerically_enclosed": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in dependencies.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Construct an endpoint-uniform Fresnel/erfc evaluator for the 230 upper complementary labels "
            "in the first cell, and repeated-integration-by-parts Hurwitz bounds for both remaining "
            "nonstationary tails. Add those complex balls to the retained lower-cell, transition, and upper-arc balls before any norm."
        ),
        "proof_boundary": (
            "Exact complementary half-lattice reduction and stationary census only. No quantitative ordinary "
            "or complete joined packet, J_Z or D_K enclosure, non-A, all-height, Lambda<=0, PF-infinity, RH, "
            "or prize-level conclusion is proved."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified ordinary complementary half-lattice reduction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
