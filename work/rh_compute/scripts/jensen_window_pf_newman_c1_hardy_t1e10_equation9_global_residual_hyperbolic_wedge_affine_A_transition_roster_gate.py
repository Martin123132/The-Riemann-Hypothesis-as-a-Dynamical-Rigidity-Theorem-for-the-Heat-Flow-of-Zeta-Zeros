#!/usr/bin/env python3
"""Certify the carrier-oriented affine tangent wedge on the full A roster."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_transition_roster_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "orientation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate.json",
    "affine_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate.json",
    "null_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_null_coordinate_unfolding_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
LOWER_START = 39_853
TARGET_END = 39_894
OUTER_START = 39_895
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))
RAY_CUTOFF = 20


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def full_lower_indented(A_quad: arb, B_linear: arb) -> acb:
    require(not A_quad.contains(0), "quadratic sign is unresolved")
    sign = 1 if A_quad.lower() > 0 else -1
    direction = acb(1, sign) / arb(2).sqrt()
    argument = direction * B_linear / (2 * abs(A_quad).sqrt())
    return acb(0, 1) * arb.pi() * (1 + argument.erf())


def outward_ray(
    A_quad: arb,
    B_linear: arb,
    r0: arb,
    cutoff: arb,
    side: int,
) -> tuple[acb, arb, arb]:
    """Integrate from r0 toward right (side=+1) or left (side=-1)."""
    require(side in (-1, 1), "invalid ray side")
    require(not A_quad.contains(0), "quadratic sign is unresolved")
    sign = 1 if A_quad.lower() > 0 else -1
    magnitude = abs(A_quad)
    sqrt_two = arb(2).sqrt()
    direction = acb(1, sign) / sqrt_two
    imaginary = acb(0, 1)

    def integrand(y: acb, analytic: bool) -> acb:
        del analytic
        R = acb(r0) + arb(side) * direction * y
        # The left real interval has reversed y orientation, so its
        # parameterized integral also carries +direction.
        return direction * (imaginary * (A_quad * R**2 + B_linear * R)).exp() / R

    derivative_at_r0 = 2 * A_quad * r0 + B_linear
    linear_growth = -arb(side * sign) * derivative_at_r0 / sqrt_two
    require(linear_growth.upper() <= 0, "ray does not point away from the quadratic saddle")
    require((cutoff - r0).lower() > 0, "ray cutoff is too short")
    value = acb.integral(
        integrand,
        arb(0),
        cutoff,
        abs_tol=arb("1e-30"),
        rel_tol=arb("1e-30"),
        eval_limit=400_000,
        depth_limit=55,
    )
    decay_rate = 2 * magnitude * cutoff - linear_growth
    require(decay_rate.lower() > 0, "ray tail is not decreasing")
    tail = (-magnitude * cutoff**2 + linear_growth * cutoff).exp() / (
        decay_rate * (cutoff - r0)
    )
    return value + acb(arb(0, tail), arb(0, tail)), tail, linear_growth


def canonical_wedge(a: arb, rho: arb, h: arb, cutoff: arb) -> tuple[acb, dict[str, Any]]:
    delta = 1 + rho
    require(not delta.contains(0), "delta sign is unresolved")
    r0 = a + delta * h
    require(r0.lower() > 0, "A-roster null coordinate is not positive")
    A1, B1 = arb("0.5"), -h
    A2, B2 = arb("0.5") - 1 / delta, a / delta

    if delta.upper() < 0:
        right1, tail1, growth1 = outward_ray(A1, B1, r0, cutoff, +1)
        right2, tail2, growth2 = outward_ray(A2, B2, r0, cutoff, +1)
        value = (
            full_lower_indented(A1, B1)
            - right1
            - full_lower_indented(A2, B2)
            + right2
        ) / (2 * arb.pi() * acb(0, 1))
        branch = "delta_negative_full_lines_minus_monotone_right_rays"
    else:
        # Here h>r0, so the lower A1 interval lies wholly on the side of its
        # saddle that points down the direct left steepest ray. This avoids
        # subtracting a right ray that first grows toward the saddle.
        left1, tail1, growth1 = outward_ray(A1, B1, r0, cutoff, -1)
        right2, tail2, growth2 = outward_ray(A2, B2, r0, cutoff, +1)
        value = (left1 + right2) / (2 * arb.pi() * acb(0, 1))
        branch = "delta_positive_lower_indented_left_ray_includes_half_jump_plus_right_ray"
    return value, {
        "delta": delta,
        "r0": r0,
        "A2": A2,
        "B2": B2,
        "tail1": tail1,
        "tail2": tail2,
        "linear_growth1": growth1,
        "linear_growth2": growth2,
        "branch": branch,
    }


def wedge_derivatives(a: arb, rho: arb, h: arb) -> tuple[acb, acb]:
    pi, imaginary = arb.pi(), acb(0, 1)
    c = rho**2 - 1
    require(not c.contains(0), "c sign is unresolved")
    sign = 1 if c.lower() > 0 else -1
    u = (abs(c) / pi).sqrt() * (h + a * rho / c)
    rotation = (imaginary * arb(sign) * pi / 4).exp()
    argument = (-imaginary * arb(sign) * pi / 4).exp() * (pi / 2).sqrt() * u
    fresnel_tail = rotation / arb(2).sqrt() * (1 - argument.erf())
    C_a = (-imaginary * a**2 / (2 * c)).exp() * fresnel_tail / (2 * (pi * abs(c)).sqrt())

    face_at_h = a + rho * h
    inner_argument = (-imaginary * pi / 4).exp() * face_at_h / arb(2).sqrt()
    inner_fresnel = (imaginary * pi / 4).exp() * (pi / 2).sqrt() * (1 + inner_argument.erf())
    C_h = -(-imaginary * h**2 / 2).exp() * inner_fresnel / (2 * pi)
    return C_a, C_h


def mode_geometry(mode_int: int) -> dict[str, Any]:
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    mode = arb(mode_int)
    denominator = 2 * pi * mode**2 + t
    alpha = 2 * mode + t / (pi * mode)
    a = (pi * mode / alpha).sqrt() * (endpoint - alpha)
    rho = -(2 * t).sqrt() * (pi * endpoint * mode + denominator) / (
        2 * denominator ** arb("1.5")
    )
    v_half = 2 * pi * mode**2 / t
    defect = v_half - 1 - v_half.log()
    root = (2 * defect).sqrt()
    h = (t / 2).sqrt() * (-root if (v_half - 1).upper() < 0 else root)
    return {"mode": mode, "a": a, "rho": rho, "h": h}


def evaluate_roster(cutoff: int, precision: int = 100) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    imaginary = acb(0, 1)
    cutoff_ball = arb(cutoff)
    normalizer = (pi / (32 * t)) ** (arb(1) / 4)
    paper_rotation = (-imaginary * pi / 8).exp()
    mu_S = arb(5) / 12 * (arb(2) / t).sqrt()

    rows: list[dict[str, Any]] = []
    scalar_sum = arb(0)
    affine_sum = arb(0)
    correction_sum = arb(0)
    absolute_affine_sum = arb(0)
    normalized_half_sum = acb(0)
    lower_affine_sum = arb(0)
    upper_affine_sum = arb(0)
    max_tail = arb(0)
    min_r0: arb | None = None
    max_r0: arb | None = None

    for mode_int in MODES:
        geometry = mode_geometry(mode_int)
        mode, a, rho, h = geometry["mode"], geometry["a"], geometry["rho"], geometry["h"]
        C, metadata = canonical_wedge(a, rho, h, cutoff_ball)
        C_a, C_h = wedge_derivatives(a, rho, h)
        K = pi * endpoint * mode
        scale = (pi * mode * endpoint).sqrt()
        lambda_P = acb(
            (3 * K**2 + 1) / (2 * (K**2 + 1) * scale),
            K / ((K**2 + 1) * scale),
        )
        C_aff = C + imaginary * mu_S * C_h - imaginary * (lambda_P + mu_S * rho) * C_a
        correction = C_aff - C
        r = t / (2 * pi * mode**2)
        W0 = acb(K, -1) * (4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt()))
        phase = t * ((t / (2 * pi)).log() - 1) / 2 - t * mode.log()
        carrier = (imaginary * phase).exp()
        raw_scalar = -carrier * W0 * C
        raw_affine = -carrier * W0 * C_aff
        raw_correction = -carrier * W0 * correction
        normalized_half = normalizer * raw_affine
        physical_scalar = 2 * normalizer * (paper_rotation * raw_scalar).real
        physical_affine = 2 * normalizer * (paper_rotation * raw_affine).real
        physical_correction = 2 * normalizer * (paper_rotation * raw_correction).real
        target_indicator = 1 if mode_int <= TARGET_END else 0

        scalar_sum += physical_scalar
        affine_sum += physical_affine
        correction_sum += physical_correction
        absolute_affine_sum += abs(physical_affine)
        normalized_half_sum += normalized_half
        if target_indicator:
            lower_affine_sum += physical_affine
        else:
            upper_affine_sum += physical_affine
        max_tail = max(max_tail, metadata["tail1"], metadata["tail2"])
        min_r0 = metadata["r0"] if min_r0 is None else min(min_r0, metadata["r0"])
        max_r0 = metadata["r0"] if max_r0 is None else max(max_r0, metadata["r0"])

        rows.append(
            {
                "mode": mode_int,
                "target_indicator": target_indicator,
                "positive_full_line_bulk_projector_coefficient": 1 - target_indicator,
                "paired_endpoint_projector_coefficient": 1,
                "branch": metadata["branch"],
                "a_face_detuning_ball": a.str(65, more=True),
                "delta_one_plus_rho_ball": (1 + rho).str(65, more=True),
                "h_half_boundary_ball": h.str(65, more=True),
                "r0_ball": metadata["r0"].str(65, more=True),
                "canonical_wedge_C": complex_record(C),
                "affine_wedge_C_aff": complex_record(C_aff),
                "normalized_affine_half_mode": complex_record(normalized_half),
                "physical_scalar_A_tangent_wedge_ball": physical_scalar.str(65, more=True),
                "physical_affine_A_tangent_wedge_ball": physical_affine.str(65, more=True),
                "physical_affine_correction_ball": physical_correction.str(65, more=True),
                "ray_tail_absolute_bounds": [metadata["tail1"].str(55, more=True), metadata["tail2"].str(55, more=True)],
            }
        )

    require(len(rows) == 84, "A transition roster count drift")
    require(sum(row["target_indicator"] for row in rows) == 42, "lower target roster count drift")
    require(sum(1 - row["target_indicator"] for row in rows) == 42, "upper roster count drift")
    require(all(row["paired_endpoint_projector_coefficient"] == 1 for row in rows), "endpoint orientation drift")
    require(all(arb(row["delta_one_plus_rho_ball"]).upper() < 0 for row in rows[:42]), "lower delta sign drift")
    require(all(arb(row["delta_one_plus_rho_ball"]).lower() > 0 for row in rows[42:]), "upper delta sign drift")
    require(max_tail.upper() < arb("1e-70"), "roster ray tail target failed")
    require(absolute_affine_sum.lower() > abs(affine_sum).upper(), "signed roster sum did not improve on termwise norm")
    cancellation_ratio = abs(affine_sum) / absolute_affine_sum
    return {
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": len(rows),
        "target_side_range": [LOWER_START, TARGET_END],
        "target_side_count": 42,
        "outer_side_range": [OUTER_START, UPPER_END],
        "outer_side_count": 42,
        "ray_cutoff": cutoff,
        "paper_normalizer_ball": normalizer.str(70, more=True),
        "minimum_r0_ball": min_r0.str(65, more=True) if min_r0 is not None else "",
        "maximum_r0_ball": max_r0.str(65, more=True) if max_r0 is not None else "",
        "maximum_ray_tail_absolute_bound": max_tail.str(55, more=True),
        "rows": rows,
        "signed_sums": {
            "target_side_physical_affine_A_tangent_wedge_ball": lower_affine_sum.str(75, more=True),
            "outer_side_physical_affine_A_tangent_wedge_ball": upper_affine_sum.str(75, more=True),
            "complete_physical_scalar_A_tangent_wedge_ball": scalar_sum.str(75, more=True),
            "complete_physical_affine_A_tangent_wedge_ball": affine_sum.str(75, more=True),
            "complete_physical_affine_correction_ball": correction_sum.str(75, more=True),
            "termwise_absolute_physical_affine_sum_ball": absolute_affine_sum.str(75, more=True),
            "signed_to_termwise_absolute_ratio_ball": cancellation_ratio.str(65, more=True),
            "complete_normalized_affine_half_sum": complex_record(normalized_half_sum),
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    s = c["signed_sums"]
    return f"""# Carrier-oriented affine A-transition roster

Date: 2026-08-13

Status: rigorous 84-mode canonical tangent-wedge sum; not an exact curved-face
or transformed-amplitude estimate

The exact A-transition roster at `t=10^10` is

```text
target side: 39853..39894, 42 modes,
outer side:  39895..39936, 42 modes.                   (AT1)
```

Every mode uses the common odd-endpoint carrier and `epsilon_A=-1` from
Section 11.431.  The paired endpoint coefficient remains `+1` across the
target edge; only the positive full-line bulk projector coefficient changes.

For `delta=1+rho<0`, the null-coordinate formula is evaluated as two
lower-indented full lines minus monotone right steepest rays.  For
`delta>0`, the first lower interval is evaluated directly on the left
steepest ray and the second upper interval directly on the right ray:

```text
C=(L_left^-(A1,B1;r0)+L_right(A2,B2;r0))/(2pi i).     (AT2)
```

This is the same Abel value as the full-minus-right formula, but its rays
point away from their quadratic saddles.  The superscript minus records that
the direct left ray passes below the pole and already contains the Abel
half-jump; no additional `1/2` is inserted.  Along each ray the modulus is
`exp(-|A|y^2+Ly)` with certified `L<=0`.  With cutoff `{c['ray_cutoff']}` the
largest omitted ray tail is

```text
{c['maximum_ray_tail_absolute_bound']}.               (AT3)
```

The null join remains numerically tight throughout:

```text
min r0={c['minimum_r0_ball']},
max r0={c['maximum_r0_ball']}.                         (AT4)
```

After restoring `2(pi/(32t))^(1/4)`, `exp(-i*pi/8)`, the common phase, and
the A endpoint sign, summing in deterministic mode order before taking any
norm gives

```text
target-side affine sum={s['target_side_physical_affine_A_tangent_wedge_ball']},
outer-side affine sum ={s['outer_side_physical_affine_A_tangent_wedge_ball']},
complete scalar sum   ={s['complete_physical_scalar_A_tangent_wedge_ball']},
complete affine sum   ={s['complete_physical_affine_A_tangent_wedge_ball']},
affine-minus-scalar   ={s['complete_physical_affine_correction_ball']}. (AT5)
```

The termwise absolute affine sum and signed cancellation ratio are

```text
sum |q_A,m^aff|={s['termwise_absolute_physical_affine_sum_ball']},
|sum q_A,m^aff|/sum |q_A,m^aff|
 ={s['signed_to_termwise_absolute_ratio_ball']}.       (AT6)
```

Equations (AT5)--(AT6) measure cancellation only inside the canonical affine
tangent-wedge A roster.  They are not a bound for the exact A block or for
`R_Dir`: the curved-face strip and transformed-amplitude remainder remain
unsummed, and the positive full-line projector jump plus the B, zero,
negative, half-current, and outer sectors must stay in the final common
assembly.

Pi provenance: every `pi` in (AT1)--(AT6) comes from equation (9), the exact
bi-Morse/null-coordinate changes, odd-endpoint parity, Gaussian contour
rotation, and the paper normalizer.  No fitted constant is introduced.

Proof boundary: rigorous canonical scalar and affine tangent-wedge values,
ray tails, and signed sums over modes 39853..39936 only.  No exact curved-face
or transformed-amplitude remainder, complete A endpoint theorem, `R_Dir`
estimate, complete `Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["orientation"]["decision"]["paired_endpoint_coefficient_continuous_across_target_edge"] is True, "orientation dependency drift")
    require(dependencies["affine_remainder"]["decision"]["affine_wedge_reduced_to_C_and_boundary_derivatives"] is True, "affine dependency drift")
    require(dependencies["null_reduction"]["decision"]["near_null_wedge_reduced_to_one_dimensional_Abel_integrals"] is True, "null dependency drift")

    certificate = evaluate_roster(RAY_CUTOFF)
    artifact = {
        "kind": STEM,
        "status": "complete_84_mode_carrier_oriented_affine_A_tangent_wedge_roster_certified_exact_remainders_open",
        "passed": True,
        "certificate": certificate,
        "decision": {
            "complete_84_mode_A_transition_roster_evaluated": True,
            "target_and_outer_sides_each_have_42_modes": True,
            "delta_negative_and_positive_branches_joined_without_c_denominator_instability": True,
            "all_steepest_rays_point_away_from_their_quadratic_saddles": True,
            "common_phase_and_paired_endpoint_orientation_retained": True,
            "signed_affine_tangent_wedge_sum_certified_before_norms": True,
            "termwise_absolute_value_shown_to_lose_cancellation": True,
            "curved_face_remainder_bound_proved": False,
            "transformed_amplitude_remainder_bound_proved": False,
            "complete_A_endpoint_block_proved": False,
            "R_Dir_bound_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
            "deterministic_mode_order": [LOWER_START, UPPER_END],
        },
        "next_obligation": "Express the exact curved A face as R_A(S)=eta_A plus its cubic null unfolding on a finite characteristic box. Bound the oriented face strip and exact-minus-affine amplitude there, then treat the complement by hyperbolic integration by parts. Preserve the 84-mode signed carrier through that assembly.",
        "proof_boundary": "Rigorous carrier-oriented canonical affine tangent-wedge roster and signed sum over modes 39853..39936 only. No exact curved-face or transformed-amplitude remainder, complete A endpoint theorem, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified carrier-oriented affine A-transition roster", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
