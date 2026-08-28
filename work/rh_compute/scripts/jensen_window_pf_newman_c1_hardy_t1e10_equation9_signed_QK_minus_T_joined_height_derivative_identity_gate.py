#!/usr/bin/env python3
"""Certify an A-free joined height-derivative identity for Q_K-T."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
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
import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_gamma_balanced_mellin_hurwitz_finite_difference_gate as mellin_gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as interval_base


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)


def result_path(stem: str) -> Path:
    return ROOT / "work" / "rh_compute" / "results" / f"{stem}.json"


DEPENDENCIES = {
    "event_atlas": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_local_height_event_atlas_gate"
    ),
    "saved_closure": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_one_sided_closure_gate"
    ),
    "signed_JZ": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_complete_signed_JZ_DK_ownership_gate"
    ),
    "physical_ownership": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_physical_transform_ownership_ledger_gate"
    ),
    "Gamma_insertion": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate"
    ),
    "A_transition": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate"
    ),
    "lower_cell": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_cell_reversed_roster_endpoint_gate"
    ),
    "quarter_arc": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_quarter_disk_vertical_arc_tail_gate"
    ),
    "ordinary_packet": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_complement_signed_coefficient_interval_gate"
    ),
    "transition_packet": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_alternating_boundary_jet_packet_gate"
    ),
    "exact_H": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_continuation_defect_target_gate"
    ),
    "QK_bridge": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_QK_exact_hardy_truncation_bridge_target_gate"
    ),
    "mellin_coordinate": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_gamma_balanced_mellin_hurwitz_finite_difference_gate"
    ),
    "lower_fold_transport": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_nonzero_height_cell_gate"
    ),
    "B_derivative_transport": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_derivative_transport_gate"
    ),
}

PRECISION_BITS = 384
HEIGHT = 10_000_000_000
TARGET_START = 622
TARGET_END = 39_894
LOWER_TRANSITION_START = 39_853
EXTRA_START = 39_895
EXTENDED_END = 39_936


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


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
            available = process.cpu_affinity()
            process.cpu_affinity([available[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        available = process.cpu_affinity()
        process.cpu_affinity([available[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def complex_record(value: acb, digits: int = 90) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def ball_record(value: arb, digits: int = 90) -> dict[str, str]:
    return {"ball": value.str(digits, more=True)}


def hardy_projection(theta: arb, value: acb) -> arb:
    return 2 * (acb(0, theta).exp() * value).real


def direct_log_H(t: arb) -> arb:
    pi = arb.pi()
    thermal = (1 + (-2 * pi * t).exp()).log()
    return (
        arb(5) * arb(2).log() / 4
        + t.log() / 4
        + 3 * pi.log() / 2
        - 3 * pi * t / 4
        - thermal
        - acb(arb(1) / 4, t / 2).lgamma().real
        - 2 * acb(arb(3) / 4, t / 2).lgamma().real
    )


def duplication_log_H(t: arb) -> arb:
    pi = arb.pi()
    log_two = arb(2).log()
    thermal = (1 + (-2 * pi * t).exp()).log()
    log_cosh = pi * t - log_two + thermal
    return (
        3 * log_two / 4
        + t.log() / 4
        + pi.log() / 2
        - 3 * pi * t / 4
        + log_cosh / 2
        - thermal
        - acb(arb(3) / 4, t / 2).lgamma().real
    )


def direct_log_H_derivative(t: arb) -> arb:
    pi = arb.pi()
    thermal = (-2 * pi * t).exp()
    psi_quarter = acb(arb(1) / 4, t / 2).digamma()
    psi_three_quarters = acb(arb(3) / 4, t / 2).digamma()
    return (
        1 / (4 * t)
        - 3 * pi / 4
        + 2 * pi * thermal / (1 + thermal)
        + psi_quarter.imag / 2
        + psi_three_quarters.imag
    )


def duplication_log_H_derivative(t: arb) -> arb:
    pi = arb.pi()
    thermal = (-2 * pi * t).exp()
    psi_three_quarters = acb(arb(3) / 4, t / 2).digamma()
    return (
        1 / (4 * t)
        - 3 * pi / 4
        + pi * (pi * t).tanh() / 2
        + 2 * pi * thermal / (1 + thermal)
        + psi_three_quarters.imag / 2
    )


def theta_and_derivative(t: arb) -> tuple[arb, arb]:
    pi = arb.pi()
    argument = acb(arb(1) / 4, t / 2)
    theta = argument.lgamma().imag - t * pi.log() / 2
    theta_prime = argument.digamma().real / 2 - pi.log() / 2
    return theta, theta_prime


def finite_dirichlet(t: arb, start: int, end: int, reverse: bool = False) -> tuple[acb, acb]:
    imaginary = acb(0, 1)
    value = acb(0)
    derivative = acb(0)
    modes = range(end, start - 1, -1) if reverse else range(start, end + 1)
    for mode in modes:
        log_mode = arb(mode).log()
        term = (-imaginary * t * log_mode).exp() / arb(mode).sqrt()
        value += term
        derivative += -imaginary * log_mode * term
    return value, derivative


def upper_tail_radius(quarter: dict[str, Any]) -> arb:
    log_ball = arb(
        quarter["arb_certificate"]["error_budgets"]["positive_real_tail_log10_upper"]["ball"]
    )
    return (log_ball.upper() * arb(10).log()).exp().upper()


def symbolic_certificate() -> dict[str, Any]:
    qk, g_target, g_extra, a_endpoint, classical = sp.symbols(
        "Q_K G_target G_extra A_endpoint T"
    )
    a_transition = a_endpoint + g_extra
    jz = qk - g_target - a_transition
    gamma_to_classical = g_target - classical
    require(
        sp.expand(jz + a_transition + gamma_to_classical - (qk - classical)) == 0,
        "physical A/projector cancellation failed",
    )

    h, h_prime, theta_prime = sp.symbols("H H_prime theta_prime", real=True)
    source, target = sp.symbols("S_W D_T")
    hardy_source, hardy_target = sp.symbols("Hardy_S Hardy_D")
    require(
        sp.expand((hardy_source - h * hardy_target) / h - (hardy_source / h - hardy_target))
        == 0,
        "direct classical-target lift failed",
    )

    theta = sp.symbols("theta", real=True)
    yr, yi, yrp, yip = sp.symbols("y_r y_i y_r_prime y_i_prime", real=True)
    log_h_prime = sp.symbols("ell_H_prime", real=True)
    c, s = sp.cos(theta), sp.sin(theta)
    direct_derivative = 2 * (
        -s * theta_prime * yr
        + c * yrp
        - c * theta_prime * yi
        - s * yip
        - log_h_prime * (c * yr - s * yi)
    ) / h
    kernel_real = yrp - log_h_prime * yr - theta_prime * yi
    kernel_imag = yip - log_h_prime * yi + theta_prime * yr
    kernel_derivative = 2 * (c * kernel_real - s * kernel_imag) / h
    require(sp.simplify(direct_derivative - kernel_derivative) == 0, "Hardy quotient derivative failed")

    source_prime, target_prime = sp.symbols("S_prime D_prime")
    imaginary = sp.I
    joined = source - h * target
    joined_prime = source_prime - h_prime * target - h * target_prime
    collapsed = sp.expand(joined_prime + (imaginary * theta_prime - h_prime / h) * joined)
    expected = sp.expand(
        source_prime
        + (imaginary * theta_prime - h_prime / h) * source
        - h * (target_prime + imaginary * theta_prime * target)
    )
    require(sp.simplify(collapsed - expected) == 0, "H-prime target cancellation failed")

    log_x = sp.symbols("log_x", real=True)
    direct_weight = 3 * sp.pi / 4 - imaginary * log_x
    mellin_ray_weight = sp.pi / 2 - imaginary * (log_x + imaginary * sp.pi / 4)
    require(sp.simplify(direct_weight - mellin_ray_weight) == 0, "Mellin derivative weight drift")

    return {
        "physical_collapse": (
            "J_Z+A_transition+(G-T)=Q_K-T, with "
            "A_transition=A_endpoint+G_extra"
        ),
        "A_free_join": "Q_K-T=Hardy_t[S_W-H*D_T]/H",
        "target_roster": "D_T=sum_(m=622)^39894 m^(-1/2-it)",
        "joined_derivative": (
            "(Q_K-T)'=Hardy_t[K_T]/H, "
            "K_T=S_W'+(i*theta'-H'/H)S_W-H(D_T'+i*theta'*D_T)"
        ),
        "finite_target_derivative": "D_T'=-i sum_(m=622)^39894 log(m)m^(-1/2-it)",
        "half_line_source_derivative": (
            "S_W'=K_0 integral_0^infinity (3pi/4-i log x)"
            "x^(-s)e^(-pi*x^2)G_W(x)dx"
        ),
        "mellin_source_derivative": (
            "S_W'=E_s integral_0^infinity [pi/2-i log(u/(2pi))]"
            "u^(-s)e^(iu^2/(4pi))Delta_q(u)du"
        ),
        "Mellin_ray_weight_matches_half_line_weight": True,
        "H_prime_target_terms_cancel_inside_joined_kernel": True,
    }


def direct_source_derivative(
    s: mp.mpc,
    first_label: int,
    count: int,
    low_cuts: tuple[mp.mpf, ...],
    high_cuts: tuple[mp.mpf, ...],
) -> mp.mpc:
    t = mp.im(s)
    k0 = mp.exp(3 * mp.pi * t / 4 + 3j * mp.pi / 8)
    integral = mellin_gate.integrate_zero_to_infinity(
        s,
        lambda x: (3 * mp.pi / 4 - 1j * mp.log(x))
        * mp.exp(-mp.pi * x * x)
        * mellin_gate.geometric_roster(x, first_label, count),
        low_cuts,
        high_cuts,
    )
    return k0 * integral


def mellin_source_derivative(
    s: mp.mpc,
    first_label: int,
    count: int,
    low_cuts: tuple[mp.mpf, ...],
    high_cuts: tuple[mp.mpf, ...],
) -> mp.mpc:
    q_minus = mp.mpf(first_label) / 2
    integral = mellin_gate.integrate_zero_to_infinity(
        s,
        lambda u: (mp.pi / 2 - 1j * mp.log(u / (2 * mp.pi)))
        * mp.exp(1j * u * u / (4 * mp.pi))
        * mellin_gate.geometric_difference(u, q_minus, count),
        low_cuts,
        high_cuts,
    )
    return mellin_gate.scaled_prefactor(s) * integral


def altered_derivative_fixture() -> dict[str, str | int]:
    mp.mp.dps = 82
    s = mp.mpc("0.5", "3.25")
    first_label = 5
    count = 5
    low_cuts = (mp.mpf(0), mp.mpf(1), mp.mpf(3), mp.mpf(7), mp.mpf(15), mp.inf)
    high_cuts = tuple(mp.mpf(item) for item in (1, 2, 4, 8, 16, 32, 64))
    direct = direct_source_derivative(s, first_label, count, low_cuts, high_cuts)
    mellin = mellin_source_derivative(s, first_label, count, low_cuts, high_cuts)
    discrepancy = abs(direct - mellin)
    require(discrepancy < mp.mpf("1e-48"), "altered derivative coordinates do not agree")
    return {
        "height": "3.25",
        "first_odd_label": first_label,
        "label_count": count,
        "direct_half_line_derivative": mp.nstr(direct, 55),
        "positive_real_mellin_derivative": mp.nstr(mellin, 55),
        "absolute_discrepancy": mp.nstr(discrepancy, 30),
    }


def production_certificate(deps: dict[str, dict[str, Any]]) -> dict[str, Any]:
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    t = arb(HEIGHT)
    theta, theta_prime = theta_and_derivative(t)
    log_h_direct = direct_log_H(t)
    log_h_duplication = duplication_log_H(t)
    h_log_prime_direct = direct_log_H_derivative(t)
    h_log_prime_duplication = duplication_log_H_derivative(t)
    require(log_h_direct.overlaps(log_h_duplication), "H formulas do not overlap")
    require(
        h_log_prime_direct.overlaps(h_log_prime_duplication),
        "H logarithmic derivative formulas do not overlap",
    )
    H = log_h_direct.exp()
    saved_H = arb(deps["exact_H"]["exact_H_certificate"]["high_precision_H_ball"])
    require(H.overlaps(saved_H), "recomputed H misses saved H")
    saved_theta = arb(deps["signed_JZ"]["certificate"]["theta_ball"]["ball"])
    require(theta.overlaps(saved_theta), "recomputed theta misses saved theta")

    D_target, D_target_prime = finite_dirichlet(t, TARGET_START, TARGET_END)
    D_lower_transition, _ = finite_dirichlet(t, LOWER_TRANSITION_START, TARGET_END)
    D_extra, _ = finite_dirichlet(t, EXTRA_START, EXTENDED_END)
    target_value = hardy_projection(theta, D_target)
    target_derivative = hardy_projection(
        theta, D_target_prime + acb(0, theta_prime) * D_target
    )
    bridge_target = arb(deps["QK_bridge"]["bridge_target_certificate"]["classical_upper_main_ball"])
    require(target_value.overlaps(bridge_target), "finite classical target misses bridge target")

    lower_cell = interval_base.saved_complex(
        deps["lower_cell"]["certificate"]["lower_cell_integral_ball"]
    )
    upper_arc = interval_base.saved_complex(deps["quarter_arc"]["arb_certificate"]["actual_arc"])
    tail_radius = upper_tail_radius(deps["quarter_arc"])
    upper_tail = acb(arb(0, tail_radius), arb(0, tail_radius))
    unowned = lower_cell + upper_arc + upper_tail
    ordinary = interval_base.saved_complex(
        deps["ordinary_packet"]["certificate"]["complete_ordinary_packet_ball"]
    )
    transition = deps["transition_packet"]
    transition_source = interval_base.saved_complex(
        transition["source_certificate"]["complete_transition_source_integral_ball"]
    )
    saved_C_G = interval_base.saved_complex(
        transition["joined_packet"]["common_Gamma_coefficient_ball"]
    )
    C_G = (1 + (-2 * arb.pi() * t).exp()) ** (-arb(1) / 2)
    require(acb(C_G).overlaps(saved_C_G), "exact C_G misses transition coefficient")

    joined_target = (
        unowned
        + ordinary
        + transition_source
        - C_G * D_lower_transition
        + (C_G - H) * D_target
    )
    direct_QKT = hardy_projection(theta, joined_target) / H

    signed = deps["signed_JZ"]["certificate"]
    J_argument = interval_base.saved_complex(signed["complete_preprojection_argument_ball"])
    natural_A = interval_base.saved_complex(
        transition["joined_packet"]["grouped_natural_A_lift_ball"]
    )
    alternate_target = J_argument + natural_A + C_G * D_extra + (C_G - H) * D_target
    require(joined_target.overlaps(alternate_target), "two A-free preprojection assemblies miss")
    alternate_QKT = hardy_projection(theta, alternate_target) / H
    require(direct_QKT.overlaps(alternate_QKT), "two A-free Q_K-T assemblies miss")

    saved_QKT = arb(deps["saved_closure"]["certificate"]["Q_K_minus_T_ball"]["ball"])
    require(direct_QKT.overlaps(saved_QKT), "A-free Q_K-T misses saved closure")
    require(direct_QKT.upper() < 0, "A-free Q_K-T assembly lost negativity")

    margin = arb(
        deps["saved_closure"]["certificate"]["target_audit"][
            "Q_K_minus_T_negativity_margin_ball"
        ]["ball"]
    )
    transport_budget = margin.lower() / 2
    radius_rows = []
    for radius_text in ("0.01", "0.001", "0.0001", "0.00001"):
        radius = arb(radius_text)
        radius_rows.append(
            {
                "radius": radius_text,
                "maximum_joined_derivative_if_half_margin_reserved": (
                    transport_budget / radius
                ).str(70, more=True),
            }
        )

    return {
        "height": HEIGHT,
        "precision_bits": PRECISION_BITS,
        "phase_and_prefactor": {
            "theta_ball": ball_record(theta),
            "theta_prime_ball": ball_record(theta_prime),
            "H_ball": ball_record(H),
            "H_logarithmic_derivative_direct_ball": ball_record(h_log_prime_direct),
            "H_logarithmic_derivative_duplication_ball": ball_record(
                h_log_prime_duplication
            ),
            "two_H_formulas_overlap": True,
            "two_H_logarithmic_derivatives_overlap": True,
        },
        "classical_target": {
            "mode_range": [TARGET_START, TARGET_END],
            "mode_count": TARGET_END - TARGET_START + 1,
            "D_T_ball": complex_record(D_target),
            "D_T_prime_ball": complex_record(D_target_prime),
            "T_ball": ball_record(target_value),
            "T_prime_ball": ball_record(target_derivative),
            "bridge_target_overlap": True,
        },
        "A_free_saved_height_assembly": {
            "packet_identity": (
                "Y_T=U_unowned+O_join+I_transition-C_G*D_(39853..39894)"
                "+(C_G-H)*D_T=S_W-H*D_T"
            ),
            "unowned_packet_ball": complex_record(unowned),
            "ordinary_packet_ball": complex_record(ordinary),
            "transition_source_ball": complex_record(transition_source),
            "C_G_ball": ball_record(C_G),
            "joined_classical_target_preprojection_ball": complex_record(joined_target),
            "alternate_preprojection_ball": complex_record(alternate_target),
            "direct_Q_K_minus_T_ball": ball_record(direct_QKT),
            "alternate_Q_K_minus_T_ball": ball_record(alternate_QKT),
            "saved_Q_K_minus_T_overlap": True,
            "A_endpoint_coefficient": 0,
            "A_transition_coefficient": 0,
            "B_trace_or_outer_input_used": False,
        },
        "first_subcell_budget": {
            "saved_negativity_margin_ball": ball_record(margin),
            "reserved_transport_budget_ball": ball_record(transport_budget),
            "reservation_rule": (
                "Reserve one half of the saved upper-endpoint negativity margin for height "
                "transport and retain one half as the final strict-sign margin."
            ),
            "candidate_radius_rows": radius_rows,
            "radius_selected": False,
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    phase = c["phase_and_prefactor"]
    target = c["classical_target"]
    joined = c["A_free_saved_height_assembly"]
    budget = c["first_subcell_budget"]
    rows = "\n".join(
        f"| `{row['radius']}` | `{row['maximum_joined_derivative_if_half_margin_reserved']}` |"
        for row in budget["candidate_radius_rows"]
    )
    return f"""# A-free joined height derivative for Q_K-T

Date: 2026-08-28

Status: exact derivative identity and saved-height A-free reassembly
certified; no nonzero-radius sign transport yet

The finite-source transform and the classical target use the same Hardy
projection.  With

```text
D_T(t)=sum_(m=622)^39894 m^(-1/2-it),
Hardy_t[X]=2 Re[exp(i theta(t))X],
```

the exact identities `Q_K=Hardy_t[S_W]/H` and `T=Hardy_t[D_T]` give

```text
Q_K-T=Hardy_t[S_W-H D_T]/H.                         (HD1)
```

This is also obtained from the previously certified route

```text
J_Z=Q_K-G-A_transition,
A_transition=A_endpoint+G_extra,
Q_K-T=J_Z+A_transition+(G-T).                       (HD2)
```

Thus the endpoint coefficient is exactly zero.  The signed B trace and outer
term belong to a valid alternate residual decomposition, but they are not
inputs to (HD1).

Write `h=H'/H`.  Differentiating (HD1) before splitting the joined object
gives

```text
(Q_K-T)'=Hardy_t[K_T]/H,                            (HD3)

K_T=S_W'+(i theta'-h)S_W
    -H[D_T'+i theta'D_T],

D_T'=-i sum_(m=622)^39894 log(m)m^(-1/2-it).        (HD4)
```

The apparent `H'D_T` terms cancel exactly inside (HD3).  The phase and
prefactor derivatives at the saved height are

```text
theta'={phase['theta_prime_ball']['ball']},
H'/H ={phase['H_logarithmic_derivative_direct_ball']['ball']}.        (HD5)
```

The second value is independently reproduced from the Gamma-duplication
formula.  Its scale is about `6.25e-32`; it was not obtained by subtracting
floating-point approximations.

For the source, differentiation is legal in two exact coordinates:

```text
S_W'=K_0 integral_0^infinity (3pi/4-i log x)
              x^(-s)e^(-pi x^2)G_W(x)dx,            (HD6)

S_W'=E_s integral_0^infinity [pi/2-i log(u/(2pi))]
              u^(-s)e^(iu^2/(4pi))Delta_q(u)du.     (HD7)
```

Near zero the dominating factor is `u^(-1/2)(1+|log u|)`, which is
integrable.  At infinity the fixed finite difference has exponential decay;
the half-line coordinate has Gaussian decay.  These majorants are uniform on
compact subcells of the Section 11.507 event cell.  On the rotated ray
`u=2pi exp(i*pi/4)x`, the weight in (HD7) is exactly the weight in (HD6).
An altered five-label fixture evaluates both derivatives independently and
finds absolute discrepancy
`{artifact['altered_derivative_fixture']['absolute_discrepancy']}`.

At `t=10^10`, direct 384-bit summation gives

```text
T ={target['T_ball']['ball']},
T'={target['T_prime_ball']['ball']}.                 (HD8)
```

The A-free packet assembly

```text
U_unowned+O_join+I_transition-C_G D_(39853..39894)
 +(C_G-H)D_T=S_W-H D_T                              (HD9)
```

produces

```text
Q_K-T={joined['direct_Q_K_minus_T_ball']['ball']},   (HD10)
```

overlapping the saved closure and remaining strictly negative.  A second
assembly adds the natural A lift back to `J_arg`, restores the extra Gamma
roster, and reaches the same ball; this checks every cancellation sign.

Half of the saved negativity margin is reserved for transport:

| radius | allowed supremum of `|(Q_K-T)'|` |
|---:|---:|
{rows}

No radius is selected by this gate.  The next certificate must bound the
single joined kernel `K_T` on one row of this ladder; separately norming
`Q_K'` and `T'` would discard the cancellation that (HD3) preserves.

Pi provenance: every `pi` in (HD1)--(HD9) comes from the inherited
Riemann--Siegel phase, exact finite Fourier/Gaussian source, Mellin quarter
turn, or Gamma normalization.  No circle, polygon, fitted constant, or
visual pattern supplies `pi`.

Proof boundary: exact A-free object identity, exact joined first derivative,
rigorous saved-height target value/derivative, and a saved-height A-free
negative reassembly only.  No bound for `K_T` on a nonzero-radius interval,
off-height sign, event-wall handoff, all-height theorem, equation-(4) global
infinite-series identity, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    for path in (*DEPENDENCIES.values(), BUILDER, CHECKER):
        require(path.is_file(), f"missing source or dependency: {path}")
    deps = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(dep.get("passed") is True for dep in deps.values()), "dependency failure")
    require(
        deps["event_atlas"]["decision"]["maximal_same_roster_cell_certified"] is True,
        "event-cell dependency drift",
    )
    require(
        deps["physical_ownership"]["identity_certificate"]["joined_source_owned_transform"]
        == "J_Z=P_t[Z]=Q_K-G-A_transition=R_KGamma-A_transition",
        "physical ownership drift",
    )
    require(
        deps["Gamma_insertion"]["decision"]["Gamma_residual_identified_as_QK_minus_G"]
        is True,
        "Gamma residual ownership drift",
    )
    require(
        deps["A_transition"]["symbolic_certificate"]["projector_completed_transition"]
        == "A_transition=A_endpoint+sum_(m=39895)^39936 G_m",
        "A-transition ownership drift",
    )
    require(
        deps["mellin_coordinate"]["decision"]["gamma_balanced_mellin_coordinate_selected"]
        is True,
        "Mellin coordinate dependency drift",
    )

    artifact = {
        "kind": STEM,
        "date": "2026-08-28",
        "status": "A_free_joined_QK_minus_T_height_derivative_identity_and_saved_height_reassembly_certified",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "altered_derivative_fixture": altered_derivative_fixture(),
        "certificate": production_certificate(deps),
        "transport_inventory": {
            "joined_source_classical_target_kernel": "selected",
            "A_endpoint_height_derivative_required_for_direct_QK_minus_T_route": False,
            "signed_B_trace_height_derivative_required_for_direct_QK_minus_T_route": False,
            "outer_term_height_derivative_required_for_direct_QK_minus_T_route": False,
            "lower_fold_height_cell_role": (
                "donor transport technology only; it is not identified with S_W or K_T and is not "
                "inserted numerically"
            ),
            "completed_B_derivative_role": (
                "retained for the alternate residual route; bypassed, not invalidated, by the direct "
                "Q_K-T identity"
            ),
            "remaining_quantitative_object": "K_T on a nonzero-radius subcell",
        },
        "decision": {
            "A_endpoint_cancels_exactly_from_QK_minus_T": True,
            "extra_Gamma_roster_rejoins_exactly_to_target_roster": True,
            "QK_minus_T_reduced_to_single_A_free_Hardy_join": True,
            "exact_joined_first_height_derivative_proved": True,
            "source_derivative_coordinates_match_on_altered_fixture": True,
            "classical_target_value_and_point_derivative_rigorously_enclosed": True,
            "saved_height_A_free_QK_minus_T_negative": True,
            "nonzero_radius_joined_derivative_bound_proved": False,
            "off_height_QK_minus_T_sign_proved": False,
            "all_height_transport_theorem_proved": False,
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
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": (
            "Differentiate the cancellation-preserving packet form in (HD9), not Q_K and T "
            "separately. Enclose U_unowned', O_join', and I_transition' together with the explicit "
            "finite Dirichlet corrections on the first radius whose joined bound fits the stored "
            "half-margin budget, then certify Q_K-T<0 throughout that subcell."
        ),
        "proof_boundary": (
            "Exact A-free Q_K-T identity, exact joined first-height derivative, production target "
            "value/derivative, and saved-height A-free negative reassembly only. No nonzero-radius "
            "joined derivative bound, off-height sign, event-wall handoff, all-height theorem, "
            "equation-(4) global infinite-series identity, Lambda<=0, PF-infinity, RH, or prize-level "
            "conclusion is proved."
        ),
    }
    interval_base.atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    interval_base.atomic_write(NOTE, render_note(artifact))
    print(
        "certified A-free joined Q_K-T height derivative identity and transport budget",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
