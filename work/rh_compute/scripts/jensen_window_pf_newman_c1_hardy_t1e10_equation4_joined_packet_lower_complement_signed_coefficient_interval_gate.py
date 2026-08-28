#!/usr/bin/env python3
"""Enclose the lower complement with signed coefficient interval bounds."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as ibp_base


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_complement_signed_coefficient_interval_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

UPPER_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_far_remainder_eight_round_gate"
UPPER_RESULT = ROOT / "work" / "rh_compute" / "results" / f"{UPPER_STEM}.json"
UPPER_BUILDER = ROOT / "work" / "rh_compute" / "scripts" / f"{UPPER_STEM}.py"
UPPER_CHECKER = UPPER_BUILDER.with_name("check_" + UPPER_BUILDER.name)
ENDPOINT_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate.json"
ENDPOINT_BUILDER = ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate.py"
ENDPOINT_CHECKER = ENDPOINT_BUILDER.with_name("check_" + ENDPOINT_BUILDER.name)
OBSTRUCTION_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate.json"
OBSTRUCTION_BUILDER = Path(ibp_base.__file__).resolve()
OBSTRUCTION_CHECKER = OBSTRUCTION_BUILDER.with_name("check_" + OBSTRUCTION_BUILDER.name)
COMPLEMENT_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate.json"
COMPLEMENT_BUILDER = ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate.py"
COMPLEMENT_CHECKER = COMPLEMENT_BUILDER.with_name("check_" + COMPLEMENT_BUILDER.name)

PRECISION_BITS = 384
ORDER = 3
MAX_AUDIT_ORDER = 9
SLABS = 16384
CLUSTER_POWER = 4
HEIGHT = arb("10000000000")
L = arb("621.5")
A = arb("39852.5")
Q_LOW = arb("79787.5")


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
        UPPER_RESULT,
        UPPER_BUILDER,
        UPPER_CHECKER,
        ENDPOINT_RESULT,
        ENDPOINT_BUILDER,
        ENDPOINT_CHECKER,
        OBSTRUCTION_RESULT,
        OBSTRUCTION_BUILDER,
        OBSTRUCTION_CHECKER,
        COMPLEMENT_RESULT,
        COMPLEMENT_BUILDER,
        COMPLEMENT_CHECKER,
        BUILDER,
        CHECKER,
    )


def clustered_right(index: int) -> arb:
    width = A - L
    denominator = SLABS**CLUSTER_POWER
    return A - width * arb((SLABS - index) ** CLUSTER_POWER) / denominator


def signed_coefficient_remainder(
    level: dict[int, dict[tuple[int, int], object]], p: arb, *, order: int
) -> tuple[arb, dict[str, Any]]:
    require(A * A < p, "physical interval crossed the Q minimum")
    total = arb(0)
    minimum_gap: arb | None = None
    for index in range(SLABS):
        slab_left = clustered_right(index)
        slab_right = clustered_right(index + 1)
        midpoint = (slab_left + slab_right) / 2
        radius = (slab_right - slab_left) / 2
        y_ball = arb(midpoint, radius)
        delta = (ibp_base.q_curve(slab_right, p) - Q_LOW).lower()
        require(delta > 0, f"lower-tail gap lost in slab {index}")
        if minimum_gap is None or delta < minimum_gap:
            minimum_gap = delta
        slab_majorant = arb(0)
        for power, terms in level.items():
            signed_coefficient = ibp_base.evaluate_coefficient(terms, p, y_ball)
            coefficient_bound = abs(signed_coefficient).upper()
            slab_majorant += coefficient_bound * ibp_base.infinite_denominator_bound(
                delta, power
            )
        total += (slab_right - slab_left) * slab_majorant
    require(minimum_gap is not None and total.is_finite(), "signed remainder census failed")
    scaled = (total / (2 * arb.pi()) ** order).upper()
    return scaled, {
        "slabs": SLABS,
        "cluster_side": "right",
        "cluster_power": CLUSTER_POWER,
        "integration_rounds": order,
        "coefficient_bound": (
            "for each denominator power, interval-evaluate the complete signed exact "
            "coefficient on the slab before taking its absolute value"
        ),
        "minimum_denominator_gap_ball": ibp_base.arb_record(minimum_gap),
        "unscaled_integral_majorant_ball": ibp_base.arb_record(total.upper()),
        "scaled_remainder_bound_ball": ibp_base.arb_record(scaled),
    }


def coefficient_cancellation_diagnostics(
    level: dict[int, dict[tuple[int, int], object]], p: arb
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for power, terms in sorted(level.items()):
        signed = abs(ibp_base.evaluate_coefficient(terms, p, A))
        triangle = ibp_base.coefficient_absolute_upper(terms, p, A)
        require(signed > 0, f"zero endpoint coefficient at power {power}")
        rows.append(
            {
                "denominator_power": power,
                "signed_coefficient_abs_at_a_ball": ibp_base.arb_record(signed),
                "old_monomial_triangle_at_a_ball": ibp_base.arb_record(triangle),
                "triangle_to_signed_ratio_ball": ibp_base.arb_record(triangle / signed),
            }
        )
    return rows


def gamma_defect_radius(complement: dict[str, Any]) -> tuple[arb, arb]:
    log_ball = arb(
        complement["certificate"]["ordinary_Gamma_defect_packet_log10_upper_ball"][
            "ball"
        ]
    )
    log_upper = log_ball.upper()
    radius = (log_upper * arb(10).log()).exp().upper()
    require(radius > 0, "Gamma defect radius vanished")
    return log_ball, radius


def build_certificate(
    endpoint: dict[str, Any],
    obstruction: dict[str, Any],
    upper: dict[str, Any],
    complement: dict[str, Any],
) -> dict[str, Any]:
    p = HEIGHT / (2 * arb.pi())
    table = ibp_base.coefficient_table(MAX_AUDIT_ORDER)
    endpoints = endpoint["certificate"]["families"]["R_lower_infinite"]
    partial, endpoint_terms = ibp_base.endpoint_expansion(
        table,
        p,
        HEIGHT,
        endpoints["endpoint_L"],
        endpoints["endpoint_a"],
        ORDER,
    )
    sweep_raw: dict[int, tuple[arb, dict[str, Any]]] = {}
    for order in range(2, MAX_AUDIT_ORDER + 1):
        sweep_raw[order] = signed_coefficient_remainder(
            table[order], p, order=order
        )
    selected_order = 2
    for order in range(3, MAX_AUDIT_ORDER + 1):
        if sweep_raw[order][0] < sweep_raw[selected_order][0]:
            selected_order = order
    require(selected_order == ORDER, "audited lower-tail order selection drift")
    remainder, census = sweep_raw[ORDER]
    lower_ball = ibp_base.add_complex_error(partial, remainder)
    upper_ball = ibp_base.saved_complex(
        upper["certificate"]["complete_upper_complement_ball"]
    )
    gamma_log10, gamma_radius = gamma_defect_radius(complement)
    ordinary_without_gamma = -lower_ball - upper_ball
    ordinary = ibp_base.add_complex_error(ordinary_without_gamma, gamma_radius)

    old_remainder = arb(
        obstruction["certificate"]["families"]["R_lower"]["remainder_census"][
            "scaled_remainder_bound_ball"
        ]["ball"]
    )
    require(remainder < old_remainder, "signed coefficient bound did not improve")

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "integration_rounds": ORDER,
        "audited_order_range": [2, MAX_AUDIT_ORDER],
        "audited_order_remainder_bounds": [
            {
                "integration_rounds": order,
                "scaled_remainder_bound_ball": ibp_base.arb_record(
                    sweep_raw[order][0]
                ),
            }
            for order in range(2, MAX_AUDIT_ORDER + 1)
        ],
        "selected_audited_order": selected_order,
        "remainder_slabs": SLABS,
        "cluster_power": CLUSTER_POWER,
        "t": HEIGHT.str(50, more=True),
        "p_equals_t_over_2pi_ball": ibp_base.arb_record(p),
        "physical_interval": [L.str(30, more=True), A.str(30, more=True)],
        "recurrence": "h_0=y^(-1/2), h_(n+1)=d/dy[h_n/(q-y-p/y)]",
        "coefficient_cancellation_diagnostics": coefficient_cancellation_diagnostics(
            table[ORDER], p
        ),
        "rejected_order_eight_coefficient_cancellation_diagnostics": (
            coefficient_cancellation_diagnostics(table[8], p)
        ),
        "lower_endpoint_terms": endpoint_terms,
        "lower_endpoint_partial_sum_ball": ibp_base.acb_record(partial),
        "lower_remainder_census": census,
        "lower_complementary_tail_ball": ibp_base.acb_record(lower_ball),
        "old_monomial_triangle_remainder_bound_ball": ibp_base.arb_record(old_remainder),
        "remainder_improvement_factor_ball": ibp_base.arb_record(old_remainder / remainder),
        "complete_upper_complement_ball": upper["certificate"][
            "complete_upper_complement_ball"
        ],
        "ordinary_Gamma_defect_log10_upper_ball": ibp_base.arb_record(gamma_log10),
        "ordinary_Gamma_defect_radius_ball": ibp_base.arb_record(gamma_radius),
        "ordinary_packet_without_Gamma_defect_ball": ibp_base.acb_record(
            ordinary_without_gamma
        ),
        "complete_ordinary_packet_ball": ibp_base.acb_record(ordinary),
        "ordinary_packet_identity": (
            "O_join=-T_lower-T_upper+Gamma_defect, with the defect bounded by "
            "(1-C_G)sum_(m=622)^39852 m^(-s)"
        ),
        "abel_limit_justification": (
            "derive at finite lower-label cutoff; the k=1 endpoint current uses the "
            "certified alternating Lerch convention, and the integrated order-three "
            "label sums are absolutely dominated by powers at least three"
        ),
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    lower = c["lower_complementary_tail_ball"]
    ordinary = c["complete_ordinary_packet_ball"]
    ratios = c["rejected_order_eight_coefficient_cancellation_diagnostics"]
    sweep = c["audited_order_remainder_bounds"]
    return f"""# Signed-coefficient interval enclosure of the lower complement

Date: 2026-08-27

Status: the lower complementary tail and complete ordinary packet at the
production height are rigorously enclosed as complex balls; the later full
joined packet remains open; this is not a proof of RH.

The old lower-tail remainder bound expanded every order-eight exact recurrence
coefficient into absolute monomials.  At the limiting endpoint `a`, that
inflates the signed coefficient by factors ranging from
`{ratios[0]['triangle_to_signed_ratio_ball']['ball']}` through
`{ratios[-1]['triangle_to_signed_ratio_ball']['ball']}`.  The largest losses
occur in denominator powers 13 through 16, exactly where the recurrence
contains strong internal cancellation.

Keep the exact collected coefficient `c_k(y)` for each denominator power and,
on every slab `Y_j`, use the rigorous interval bound

```text
sup_(y in Y_j)|c_k(y)|
 <= abs(c_k(Y_j)).upper.                            (LC1)
```

Only after (LC1) is formed is it multiplied by the positive denominator sum

```text
sum_(n>=0)(delta+n)^(-k)
 <= delta^(-k)+delta^(1-k)/(k-1).                  (LC2)
```

This preserves cancellation internal to each exact coefficient without
assuming cancellation between labels, powers, or slabs.  A rigorous sweep of
the first admissible orders gives

```text
order 2: {sweep[0]['scaled_remainder_bound_ball']['ball']},
order 3: {sweep[1]['scaled_remainder_bound_ball']['ball']},
order 4: {sweep[2]['scaled_remainder_bound_ball']['ball']},
order 5: {sweep[3]['scaled_remainder_bound_ball']['ball']},
order 6: {sweep[4]['scaled_remainder_bound_ball']['ball']},
order 7: {sweep[5]['scaled_remainder_bound_ball']['ball']},
order 8: {sweep[6]['scaled_remainder_bound_ball']['ball']},
order 9: {sweep[7]['scaled_remainder_bound_ball']['ball']}.
```

Order one is excluded because its integrated infinite label sum still has a
denominator-power-one term and is not absolutely summable.  Thus order three
is the certified minimum over orders 2 through 9 under this slab rule.  Using
16384 fourth-power right-endpoint-clustered slabs gives

```text
T_lower endpoint partial =
  {c['lower_endpoint_partial_sum_ball']['real_ball']}
 +i{c['lower_endpoint_partial_sum_ball']['imag_ball']},

T_lower remainder <= {c['lower_remainder_census']['scaled_remainder_bound_ball']['ball']},

T_lower = {lower['real_ball']}
         +i{lower['imag_ball']}.                    (LC3)
```

The former monomial-triangle remainder bound was
`{c['old_monomial_triangle_remainder_bound_ball']['ball']}`; the certified
improvement factor is `{c['remainder_improvement_factor_ball']['ball']}`.

The exact complementary identity is

```text
O_join=-T_lower-T_upper
       +(1-C_G)sum_(m=622)^39852 m^(-s).            (LC4)
```

The last packet has log10 absolute upper bound
`{c['ordinary_Gamma_defect_log10_upper_ball']['ball']}`.  Joining it as a
complex error ball to the certified complete upper complement and (LC3)
before taking a norm gives

```text
O_join = {ordinary['real_ball']}
        +i{ordinary['imag_ball']},

|O_join| = {ordinary['absolute_ball']}.             (LC5)
```

The independent checker changes to four recurrence rounds, 20000
fifth-power-clustered slabs, and explicit alternating eta endpoint values.
Its lower-tail and ordinary-packet balls overlap (LC3) and (LC5).  A separate
altered-height negative-denominator direct quadrature overlaps the corresponding
fifth-order recurrence enclosure.

Pi provenance: every `pi` comes from the inherited Fresnel/Fourier phase,
`phi'=2*pi*D`, alternating Fourier characters, or `p=t/(2*pi)`.  The signed
coefficient collection follows exact differentiation; no circle, polygon,
fitted constant, or visual pattern supplies `pi`.

Proof boundary: rigorous signed-coefficient order-three lower-tail enclosure
and the resulting complete production-height ordinary-packet complex ball
only.  No complete later joined packet, `J_Z`, `D_K`, non-A, all-height,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    resource_mode = ibp_base.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    for path in source_paths():
        require(path.is_file(), f"missing source or dependency: {path}")

    endpoint = load_json(ENDPOINT_RESULT)
    obstruction = load_json(OBSTRUCTION_RESULT)
    upper = load_json(UPPER_RESULT)
    complement = load_json(COMPLEMENT_RESULT)
    for name, dependency in (
        ("endpoint", endpoint),
        ("obstruction", obstruction),
        ("upper", upper),
        ("complement", complement),
    ):
        require(dependency.get("passed") is True, f"{name} dependency not passed")
    require(
        upper["decision"]["complete_upper_complement_complex_ball_certified"] is True,
        "upper complement proof boundary drift",
    )
    require(
        obstruction["decision"]["lower_complementary_tail_complex_ball_certified"] is True,
        "old lower-tail obstruction missing",
    )

    certificate = build_certificate(endpoint, obstruction, upper, complement)
    artifact = {
        "kind": "rh_c1_hardy_equation4_lower_complement_signed_coefficient_interval_gate",
        "date": "2026-08-27",
        "status": "lower_complement_and_complete_ordinary_packet_complex_balls_certified",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "signed_coefficient_interval_bound_certified": True,
            "order_three_minimal_over_audited_orders_two_through_nine": True,
            "lower_complementary_tail_complex_ball_certified": True,
            "complete_upper_complement_dependency_retained": True,
            "complete_ordinary_packet_complex_ball_certified": True,
            "complete_later_joined_packet_enclosed": False,
            "D_K_corridor_proved": False,
        },
        "sources": {relative(path): file_hash(path) for path in source_paths()},
        "proof_boundary": (
            "Rigorous signed-coefficient order-three lower-complement enclosure and the resulting "
            "complete production-height ordinary-packet complex ball only. No complete later "
            "joined packet, J_Z, D_K, non-A, all-height, Lambda<=0, PF-infinity, RH, or "
            "prize-level conclusion is proved."
        ),
    }
    ibp_base.atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    ibp_base.atomic_write(NOTE, note_text(artifact))
    print("certified lower complement and complete ordinary packet", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
