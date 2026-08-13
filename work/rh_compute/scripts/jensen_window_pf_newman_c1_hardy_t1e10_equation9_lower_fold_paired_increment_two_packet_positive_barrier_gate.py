#!/usr/bin/env python3
"""Certify that the selector fold increment is a large positive source term."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "common_profile": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "exact_profile": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.json",
    "selector_increment": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate.json",
    "companion_decomposition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate.json",
    "inner_outer_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_inner_outer_corrected_defect_projection_gate.json",
}

C = 159_577
LOW_LENGTH = 200
HIGH_FIRST = 79_590
HIGH_LAST = 79_789
MIDPOINT_PANELS = 65_536
DERIVATIVE_PANELS = 2_048
PRECISION = 70


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def beta4_weights(lam: arb, y: arb, beta: arb) -> tuple[arb, arb, arb, arb, arb, arb]:
    epsilon = beta**-2
    u1 = -(13 * lam + 3 * y) / 60
    v1 = (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    n_u2 = 448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    n_v2 = 40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27
    u = 1 + epsilon * u1 - epsilon**2 * n_u2 / 50400
    v = epsilon * v1 + epsilon**2 * n_v2 / 1680
    n_u2_y = 840 * lam**2 * y**2 - 420 * lam * y**3 + 30 * lam + 315 * y**4 - 810 * y
    n_v2_y = -20 * lam**2 + 2 * lam * y - 27 * y**2
    u_y = -epsilon / 20 - epsilon**2 * n_u2_y / 50400
    v_y = epsilon * (-4 * lam + 6 * y) / 60 + epsilon**2 * n_v2_y / 1680
    n_u2_yy = 1680 * lam**2 * y - 1260 * lam * y**2 + 1260 * y**3 - 810
    n_v2_yy = 2 * lam - 54 * y
    u_yy = -epsilon**2 * n_u2_yy / 50400
    v_yy = epsilon / 10 + epsilon**2 * n_v2_yy / 1680
    return u, v, u_y, v_y, u_yy, v_yy


def profile_derivatives(lam: arb, y: arb, beta: arb, y_max: arb) -> tuple[acb, acb, acb]:
    x = lam + y
    ai, ai_prime, _, _ = acb(-x).airy()
    a = ai
    a_x = -ai_prime
    u, v, u_y, v_y, u_yy, v_yy = beta4_weights(lam, y, beta)
    h = u * a + v * a_x
    h_y = (u_y - x * v) * a + (u + v_y) * a_x
    h_yy = (u_yy - v - x * u - 2 * x * v_y) * a + (2 * u_y + v_yy - x * v) * a_x
    psi = y**2 / (4 * beta) - 3 * arb.pi() * y / (2 * y_max)
    psi_y = y / (2 * beta) - 3 * arb.pi() / (2 * y_max)
    psi_yy = 1 / (2 * beta)
    phase = acb(0, psi).exp()
    g = y_max * phase * h
    g_s = y_max**2 * phase * (h_y + acb(0, psi_y) * h)
    g_ss = y_max**3 * phase * (
        h_yy + 2 * acb(0, psi_y) * h_y + (acb(0, psi_yy) - psi_y**2) * h
    )
    return g, g_s, g_ss


def derivative_atlas(lam: arb, beta: arb, y_max: arb, panels: int = DERIVATIVE_PANELS) -> dict[str, arb]:
    maxima = [arb(0), arb(0), arb(0)]
    half_width = y_max / (2 * panels)
    for index in range(panels):
        midpoint = y_max * arb(2 * index + 1) / (2 * panels)
        y = arb(midpoint, half_width)
        values = profile_derivatives(lam, y, beta, y_max)
        maxima = [max(old, abs(value).upper()) for old, value in zip(maxima, values)]
    return {"g": maxima[0], "g_s": maxima[1], "g_ss": maxima[2]}


def low_packet_value(lam: arb, beta: arb, y_max: arb, panels: int = MIDPOINT_PANELS) -> acb:
    total = acb(0)
    two_pi_i = acb(0, 2 * arb.pi())
    for index in range(panels):
        s = arb(2 * index + 1) / (2 * panels)
        y = y_max * s
        g, _, _ = profile_derivatives(lam, y, beta, y_max)
        z = (two_pi_i * s).exp()
        packet = z * (1 - z**LOW_LENGTH) / (1 - z)
        total += g * packet
    return total / panels


def low_packet_midpoint_error(derivatives: dict[str, arb], panels: int = MIDPOINT_PANELS) -> arb:
    count = arb(LOW_LENGTH)
    sum_k = arb(LOW_LENGTH * (LOW_LENGTH + 1) // 2)
    sum_k2 = arb(LOW_LENGTH * (LOW_LENGTH + 1) * (2 * LOW_LENGTH + 1) // 6)
    p0 = count
    p1 = 2 * arb.pi() * sum_k
    p2 = (2 * arb.pi()) ** 2 * sum_k2
    second = derivatives["g_ss"] * p0 + 2 * derivatives["g_s"] * p1 + derivatives["g"] * p2
    return (second / (24 * panels**2)).upper()


def high_packet_by_parts(lam: arb, beta: arb, y_max: arb, g_ss_bound: arb) -> dict[str, Any]:
    g0, gs0, _ = profile_derivatives(lam, arb(0), beta, y_max)
    g1, gs1, _ = profile_derivatives(lam, y_max, beta, y_max)
    first_sum = acb(0)
    second_sum = arb(0)
    for k in range(HIGH_FIRST, HIGH_LAST + 1):
        n = -k
        omega = 2 * arb.pi() * n
        first_sum += 1 / acb(0, -omega)
        second_sum += 1 / omega**2
    center = (g1 - g0) * first_sum + (gs1 - gs0) * second_sum
    remainder = (g_ss_bound * second_sum).upper()
    return {
        "center": center,
        "remainder": remainder,
        "first_boundary_sum": first_sum,
        "second_boundary_sum": second_sum,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Paired selector-increment two-packet positive barrier

Date: 2026-08-13

Status: rigorous positive lower bound for the selector fold increment; this
rejects the signed-increment route and is not a fixed-state residual theorem

For the exact common-profile coefficient

```text
G_ex,(39895+n)=integral_0^1 g_ex(s)e^(-2pi*i*n*s)ds,
```

the paired fold increment has exactly two 200-frequency packets:

```text
D_fold/kappa_C=sum_(n=-200)^-1 G_ex,(39895+n)
              +sum_(n=-79789)^-79590 G_ex,(39895+n),  (PB1)
kappa_C=e^(i beta^3)2sqrt(2).
```

The first packet is evaluated directly from the beta-minus-four profile by
a `{MIDPOINT_PANELS}`-panel Arb midpoint rule.  A `{DERIVATIVE_PANELS}`-panel
rectangle atlas bounds the profile and its first two `s` derivatives, giving

```text
{c['low_beta4_real_lower_ball']} < Re low beta4
 < {c['low_beta4_real_upper_ball']},
low midpoint error <{c['low_midpoint_error_ball']}.    (PB2)
```

The remote packet is integrated by parts twice.  Integer endpoint phases are
one, so its two boundary currents are explicit and the remaining integral is
bounded by `sup|g4''| sum(2pi n)^-2`:

```text
{c['high_beta4_real_lower_ball']} < Re high beta4
 < {c['high_beta4_real_upper_ball']},
high second-current remainder <{c['high_remainder_ball']}. (PB3)
```

Finally the exact Kummer-ODE theorem gives a uniform profile error.  The two
length-200 Dirichlet packets have the same `L1` majorant, so the complete
exact-minus-beta-minus-four packet error is

```text
canonical packet error <{c['exact_profile_two_packet_error_ball']}.      (PB4)
```

Combining (PB1)--(PB4), restoring the carrier, and applying the exact paired
projection proves uniformly on the ordinary top corridor

```text
20<{c['physical_exact_fold_increment_lower_ball']}<Delta Q_fold
 <{c['physical_exact_fold_increment_upper_ball']}<35.                  (PB5)
```

Thus the selector increment is one large removed source contribution, not a
small signed correction.  Its positive scale is more than
`{c['ratio_to_inner_corrected_margin_lower_ball']}` times the certified inner
corrected margin.  Telescoping these increments without simultaneous target
and ownership updates cannot close the fixed-state `Q_K-T` residual.

Pi provenance: every frequency and projection factor comes from the exact
Kummer/Fourier character, `beta^3=pi*C^2/8`, and odd-square half-domain
reflection.  No fitted constant is used.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No fixed-state fold residual, moving-target telescope, complete `Q_K-T` or
`T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    ctx.dps = PRECISION
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in DEPENDENCIES.items()
    }
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["selector_increment"]["decision"]["fold_increment_equals_signed_strip_pair_block"] is True,
        "paired increment identity drift",
    )
    require(
        dependencies["common_profile"]["decision"]["exact_transformed_strip_common_profile_identity_proved"] is True,
        "common-profile identity drift",
    )

    c = arb(C)
    tstar = arb.pi() * c * c / 8
    beta = tstar ** (arb(1) / 3)
    lambda_max = arb.pi() / (16 * beta)
    lam = arb(lambda_max / 2, lambda_max / 2)
    y_max = arb.pi() * c / (2 * beta)
    derivatives = derivative_atlas(lam, beta, y_max)
    low = low_packet_value(lam, beta, y_max)
    low_error = low_packet_midpoint_error(derivatives)
    low_real_lower = (low.real.lower() - low_error).lower()
    low_real_upper = (low.real.upper() + low_error).upper()
    high = high_packet_by_parts(lam, beta, y_max, derivatives["g_ss"])
    high_real_lower = (high["center"].real.lower() - high["remainder"]).lower()
    high_real_upper = (high["center"].real.upper() + high["remainder"]).upper()

    exact_profile_sup = arb(
        dependencies["exact_profile"]["interval_certificate"]["uniform_canonical_exact_minus_beta4_profile_bound"]
    ).upper()
    packet_l1 = (1 + (arb.pi() * LOW_LENGTH / 2).log()).upper()
    exact_packet_error = (2 * packet_l1 * exact_profile_sup).upper()
    exact_real_lower = (low_real_lower + high_real_lower - exact_packet_error).lower()
    exact_real_upper = (low_real_upper + high_real_upper + exact_packet_error).upper()
    physical_lower = (4 * arb(2).sqrt() * exact_real_lower).lower()
    physical_upper = (4 * arb(2).sqrt() * exact_real_upper).upper()
    require(physical_lower > 20, "paired selector increment did not clear positive barrier 20")
    require(physical_upper < 35, "paired selector increment did not clear upper barrier 35")

    inner_margin = arb(
        dependencies["inner_outer_projection"]["certificate"]["sufficient_paired_companion_real_upper_target"]
    ).lower()
    ratio = (physical_lower / inner_margin).lower()
    require(ratio > 290_000, "positive barrier ratio below 290000")

    certificate = {
        "height_interval": "t*-pi/16<=t<=t*",
        "beta_ball": beta.str(PRECISION, more=True),
        "Y_ball": y_max.str(PRECISION, more=True),
        "profile_modulus_bound_ball": derivatives["g"].str(PRECISION, more=True),
        "profile_first_derivative_bound_ball": derivatives["g_s"].str(PRECISION, more=True),
        "profile_second_derivative_bound_ball": derivatives["g_ss"].str(PRECISION, more=True),
        "low_beta4_midpoint_ball": {"real_ball": low.real.str(PRECISION, more=True), "imag_ball": low.imag.str(PRECISION, more=True)},
        "low_midpoint_error_ball": low_error.str(PRECISION, more=True),
        "low_beta4_real_lower_ball": low_real_lower.str(PRECISION, more=True),
        "low_beta4_real_upper_ball": low_real_upper.str(PRECISION, more=True),
        "high_beta4_boundary_center_ball": {"real_ball": high["center"].real.str(PRECISION, more=True), "imag_ball": high["center"].imag.str(PRECISION, more=True)},
        "high_remainder_ball": high["remainder"].str(PRECISION, more=True),
        "high_beta4_real_lower_ball": high_real_lower.str(PRECISION, more=True),
        "high_beta4_real_upper_ball": high_real_upper.str(PRECISION, more=True),
        "single_packet_L1_majorant_ball": packet_l1.str(PRECISION, more=True),
        "exact_profile_two_packet_error_ball": exact_packet_error.str(PRECISION, more=True),
        "canonical_exact_fold_increment_real_lower_ball": exact_real_lower.str(PRECISION, more=True),
        "canonical_exact_fold_increment_real_upper_ball": exact_real_upper.str(PRECISION, more=True),
        "physical_exact_fold_increment_lower_ball": physical_lower.str(PRECISION, more=True),
        "physical_exact_fold_increment_upper_ball": physical_upper.str(PRECISION, more=True),
        "inner_corrected_margin_ball": inner_margin.str(PRECISION, more=True),
        "ratio_to_inner_corrected_margin_lower_ball": ratio.str(PRECISION, more=True),
    }
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "rigorous_positive_selector_fold_increment_barrier_certified",
        "passed": True,
        "exact_packet_identity": {
            "profile_coefficient": "G_ex,(39895+n)=int_0^1 g_ex(s)exp(-2pi*i*n*s)ds",
            "low_packet": [-200, -1],
            "remote_packet": [-79_789, -79_590],
            "packet_length_each": 200,
            "fold_increment": "D_fold=kappa_C*(sum_low G_ex,n+sum_remote G_ex,n)",
            "physical_projection": "Delta Q_fold=4sqrt(2)Re(sum of exact packets)",
        },
        "certificate": certificate,
        "decision": {
            "exact_two_packet_fold_increment_identity_proved": True,
            "low_packet_rigorous_midpoint_enclosure_proved": True,
            "remote_packet_two_current_enclosure_proved": True,
            "uniform_exact_profile_error_applied_to_both_packets": True,
            "selector_fold_increment_uniformly_positive_above_20": True,
            "selector_fold_increment_uniformly_below_35": True,
            "selector_increment_is_small_residual_route": False,
            "signed_negative_selector_increment_proved": False,
            "fixed_state_fold_residual_bounded": False,
            "complete_T_upper_proved": False,
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
            "workers": 1,
            "process_priority": priority,
            "midpoint_panels": MIDPOINT_PANELS,
            "derivative_panels": DERIVATIVE_PANELS,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Abandon the fixed-target selector-increment sign route. Return to one fixed selector state and derive the "
            "cancellation-preserving A-endpoint paired residual against its allocated target carrier, or formulate a "
            "moving-target telescope in which the large removed source term is canceled explicitly at every step."
        ),
        "proof_boundary": (
            "Rigorous positive lower bound for one fixed-B selector fold increment on one ordinary top corridor only. "
            "No fixed-state fold residual, moving-target telescope, complete Q_K-T or T_upper, height-uniform theorem, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print(f"certified selector fold increment > {physical_lower}; fixed-target increment route rejected", flush=True)


if __name__ == "__main__":
    main()
