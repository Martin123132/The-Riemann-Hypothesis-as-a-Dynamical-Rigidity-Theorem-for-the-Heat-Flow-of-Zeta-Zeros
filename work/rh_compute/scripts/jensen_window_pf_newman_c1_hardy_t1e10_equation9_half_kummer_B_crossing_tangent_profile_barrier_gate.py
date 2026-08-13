#!/usr/bin/env python3
"""Certify the exact tangent half-plane model for the B-endpoint crossing."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_B_crossing_tangent_profile_barrier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "half_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
    "universal_morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
    "Gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "symmetric_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.json",
}

PRECISION = 110
T = 10_000_000_000
B = 5_122_421
LOW_END = 621
TARGET_START = 622
TARGET_END = 39_894
NORMALIZED_TARGET = arb("0.000019")
PHYSICAL_TARGET = arb("0.0000086")
WITNESS_MODES = (1, 256, 257, 258, 620, 621, 622, 623, 624, 1000, 39_894)
SEGMENTS = (
    ("positive_delta_low", 1, 257),
    ("negative_delta_low", 258, 599),
    ("B_crossing_collar", 600, 650),
    ("target_near", 651, 999),
    ("target_bulk", 1000, 39_894),
)


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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def symbolic_certificate() -> dict[str, str]:
    m, t, endpoint, v = sp.symbols("m t B v", positive=True, real=True)
    pi = sp.pi
    r = t / (2 * pi * m**2)
    x = 1 / (1 + r * v)
    q = sp.sqrt(x / 2) * (endpoint - 2 * m / x)
    kappa = sp.factor(2 * sp.diff(q, v).subs(v, 1) / sp.sqrt(t))
    expected_kappa = -sp.sqrt(t) * (pi * endpoint * m + 2 * pi * m**2 + t) / (
        sp.sqrt(pi) * (2 * pi * m**2 + t) ** sp.Rational(3, 2)
    )
    require(sp.simplify(kappa - expected_kappa) == 0, "B-boundary tangent slope failed")

    y, rho, k = sp.symbols("y rho k", real=True)
    delta = 1 - pi * k**2 / 2
    quadratic = -y**2 + pi * (rho + k * y) ** 2 / 2
    completed = -delta * (y - pi * k * rho / (2 * delta)) ** 2 + pi * rho**2 / (2 * delta)
    require(sp.simplify(quadratic - completed) == 0, "two-dimensional quadratic completion failed")

    center_y = pi * k * rho / (2 * delta)
    centered_q = sp.simplify(rho + k * center_y)
    require(sp.simplify(centered_q - rho / delta) == 0, "q-moment center failed")

    return {
        "outer_scale": "y=sqrt(t)s/2, so the universal Morse Gaussian is exp(-i*y^2)",
        "linear_B_boundary": "q_B(s)=q_B*+kappa_B*y in the tangent model",
        "tangent_slope": str(kappa),
        "quadratic_defect": "Delta_B=1-pi*kappa_B^2/2",
        "half_plane_phase": "-y^2+(pi/2)(rho+kappa*y)^2=-Delta(y-pi*kappa*rho/(2Delta))^2+pi*rho^2/(2Delta)",
        "inner_density": "After factoring the classical m-density, the exact grouped tangent density is 1+c_m*q with c_m=sqrt(x_m)/(m*sqrt(2)).",
        "rotation": "chi(Delta)=1 for Delta>0 and i for Delta<0",
        "scalar_tail": "U_0=chi*T_signDelta(q*/sqrt(|Delta|))/(1+i)",
        "q_density_tail": "U_1=-c_m*chi*exp(i*pi*q*^2/(2Delta))/[sqrt(|Delta|)*i*pi*(1+i)]",
        "crossing_residual": "R_B,m=1_(m<=621)-(U_0+U_1)",
        "stable_small_tail_form": "Because q_B*<0 through m=621 and q_B*>0 from m=622, R_0=-sign(q_B*)chi*T_signDelta(|q_B*|/sqrt(|Delta|))/(1+i).",
        "classical_carrier": "exp(i(theta_0-t log m))/sqrt(m)",
        "interpretation": "This evaluates the tangent half-plane model exactly; it does not bound the nonlinear Morse amplitude or the remaining symmetric Poisson completion.",
    }


def fresnel_tail_positive(u: arb, imaginary_unit: acb, pi: arb) -> acb:
    full = (1 + imaginary_unit) / 2
    argument = (-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt() * u
    primitive = full * argument.erf()
    return full - primitive


def mode_profile(mode: int, pi: arb, imaginary_unit: acb) -> dict[str, Any]:
    m = arb(mode)
    t = arb(T)
    endpoint = arb(B)
    x = 2 * pi * m**2 / (t + 2 * pi * m**2)
    alpha_saddle = 2 * m + t / (pi * m)
    q = m * (pi / (t + 2 * pi * m**2)).sqrt() * (endpoint - alpha_saddle)
    kappa = -t.sqrt() * (pi * endpoint * m + 2 * pi * m**2 + t) / (
        pi.sqrt() * (2 * pi * m**2 + t) ** (arb(3) / 2)
    )
    delta = 1 - pi * kappa**2 / 2
    require(not delta.contains(0), f"integer tangent defect contains zero at mode {mode}")
    require((mode <= LOW_END and q < 0) or (mode >= TARGET_START and q > 0), f"B crossing sign drift at mode {mode}")

    abs_delta = abs(delta)
    u = abs(q) / abs_delta.sqrt()
    positive_tail = fresnel_tail_positive(u, imaginary_unit, pi)
    if delta > 0:
        rotation = acb(1)
        signed_tail = positive_tail
        delta_sign = 1
    else:
        rotation = imaginary_unit
        signed_tail = positive_tail.conjugate()
        delta_sign = -1

    sign_q = -1 if q < 0 else 1
    scalar_profile = -sign_q * rotation * signed_tail / (1 + imaginary_unit)
    coefficient = x.sqrt() / (m * arb(2).sqrt())
    effective_phase = (imaginary_unit * pi * q**2 / (2 * delta)).exp()
    q_density_profile = (
        coefficient
        * rotation
        * effective_phase
        / (abs_delta.sqrt() * imaginary_unit * pi * (1 + imaginary_unit))
    )
    profile = scalar_profile + q_density_profile
    return {
        "mode": mode,
        "x": x,
        "q": q,
        "kappa": kappa,
        "delta": delta,
        "delta_sign": delta_sign,
        "u": u,
        "scalar_profile": scalar_profile,
        "q_density_profile": q_density_profile,
        "profile": profile,
    }


def aggregate_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    t = arb(T)
    theta_zero = t / 2 * ((t / (2 * pi)).log() - 1) - pi / 8
    theta_exact = acb(arb("0.25"), t / 2).lgamma().imag - t * pi.log() / 2
    theta_shift = theta_exact - theta_zero

    total_scalar = acb(0)
    total_q_density = acb(0)
    triangle = arb(0)
    segment_sums = {name: acb(0) for name, _, _ in SEGMENTS}
    segment_triangles = {name: arb(0) for name, _, _ in SEGMENTS}
    witnesses: list[dict[str, Any]] = []
    minimum_abs_delta = None
    minimum_abs_delta_mode = None

    for mode in range(1, TARGET_END + 1):
        row = mode_profile(mode, pi, imaginary_unit)
        carrier = (imaginary_unit * (theta_zero - t * arb(mode).log())).exp() / arb(mode).sqrt()
        scalar_term = carrier * row["scalar_profile"]
        q_term = carrier * row["q_density_profile"]
        term = scalar_term + q_term
        total_scalar += scalar_term
        total_q_density += q_term
        triangle += abs(term)
        for name, start, end in SEGMENTS:
            if start <= mode <= end:
                segment_sums[name] += term
                segment_triangles[name] += abs(term)
                break
        abs_delta = abs(row["delta"])
        if minimum_abs_delta is None or abs_delta < minimum_abs_delta:
            minimum_abs_delta = abs_delta
            minimum_abs_delta_mode = mode
        if mode in WITNESS_MODES:
            witnesses.append(
                {
                    "mode": mode,
                    "q_B_star_ball": row["q"].str(PRECISION, more=True),
                    "kappa_B_ball": row["kappa"].str(PRECISION, more=True),
                    "Delta_B_ball": row["delta"].str(PRECISION, more=True),
                    "scaled_tail_argument_ball": row["u"].str(PRECISION, more=True),
                    "scalar_profile_ball": complex_record(row["scalar_profile"]),
                    "q_density_profile_ball": complex_record(row["q_density_profile"]),
                    "grouped_profile_ball": complex_record(row["profile"]),
                }
            )

    total_theta_zero = total_scalar + total_q_density
    phase_rotation = (imaginary_unit * theta_shift).exp()
    total_exact_theta = phase_rotation * total_theta_zero
    scalar_exact_theta = phase_rotation * total_scalar
    q_exact_theta = phase_rotation * total_q_density
    two_real = 2 * total_exact_theta.real
    scalar_two_real = 2 * scalar_exact_theta.real
    q_two_real = 2 * q_exact_theta.real
    physical_two_real = arb(2).sqrt() / pi * two_real
    target_ratio = two_real / NORMALIZED_TARGET
    retained_triangle_ratio = abs(total_exact_theta) / triangle

    require(minimum_abs_delta_mode == 257, "nearest integer tangent defect mode drift")
    require(minimum_abs_delta > arb("0.0009"), "integer tangent defect too close to zero")
    require(scalar_two_real > 0 and scalar_two_real < NORMALIZED_TARGET, "scalar-only diagnostic no longer lies below target")
    require(q_two_real > arb("0.000013"), "q-density channel unexpectedly small")
    require(two_real > NORMALIZED_TARGET and two_real < arb("0.000023"), "grouped tangent barrier interval failed")
    require(physical_two_real > PHYSICAL_TARGET, "physical tangent barrier failed")
    require(triangle > arb("0.0027"), "tangent profile triangle drift")
    require(abs(theta_shift) < arb("2.1e-12"), "theta shift drift")

    return {
        "height": T,
        "upper_endpoint": B,
        "mode_range": [1, TARGET_END],
        "mode_count": TARGET_END,
        "step_boundary": [LOW_END, TARGET_START],
        "precision_decimal_digits": PRECISION,
        "minimum_integer_abs_Delta_B_mode": minimum_abs_delta_mode,
        "minimum_integer_abs_Delta_B_ball": minimum_abs_delta.str(PRECISION, more=True),
        "theta_exact_minus_theta_zero_ball": theta_shift.str(PRECISION, more=True),
        "scalar_only_exact_theta_aggregate_ball": complex_record(scalar_exact_theta),
        "scalar_only_two_real_ball": scalar_two_real.str(PRECISION, more=True),
        "q_density_exact_theta_aggregate_ball": complex_record(q_exact_theta),
        "q_density_two_real_ball": q_two_real.str(PRECISION, more=True),
        "grouped_exact_theta_aggregate_ball": complex_record(total_exact_theta),
        "grouped_two_real_ball": two_real.str(PRECISION, more=True),
        "grouped_physical_two_real_ball": physical_two_real.str(PRECISION, more=True),
        "normalized_target": str(NORMALIZED_TARGET),
        "physical_target": str(PHYSICAL_TARGET),
        "grouped_to_normalized_target_ratio_ball": target_ratio.str(PRECISION, more=True),
        "termwise_complex_triangle_ball": triangle.str(PRECISION, more=True),
        "signed_complex_to_triangle_ratio_ball": retained_triangle_ratio.str(PRECISION, more=True),
        "segments": [
            {
                "name": name,
                "range": [start, end],
                "count": end - start + 1,
                "grouped_sum_ball": complex_record(segment_sums[name] * phase_rotation),
                "termwise_triangle_ball": segment_triangles[name].str(PRECISION, more=True),
            }
            for name, start, end in SEGMENTS
        ],
        "witnesses": witnesses,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Tangent half-plane model for the B crossing

Date: 2026-08-13

Status: countermodel gate; not a proof or bound for the nonlinear half-Kummer residual

On the half-domain from Section 11.352, scale the universal Morse coordinate
by

```text
y=sqrt(t)s/2,
q_B(s)=q_B*+kappa_B y                              (tangent model),
Delta_B=1-pi*kappa_B^2/2.
```

After factoring the classical full-line density, the complete grouped alpha
current has tangent density

```text
1+c_m q,       c_m=sqrt(x_m)/(m sqrt(2)).                (BT1)
```

The second term is the exact current corresponding to the sine/pivot-
derivative channel.  It must not be dropped.

Put `q=q_B*+kappa_B y` and then `rho=q-kappa_B y`.  The joint phase completes
exactly:

```text
-y^2+(pi/2)(rho+kappa_B y)^2
 =-Delta_B[y-pi kappa_B rho/(2Delta_B)]^2
   +pi rho^2/(2Delta_B).                                (BT2)
```

Let `chi=1` for `Delta_B>0`, `chi=i` for `Delta_B<0`, and let

```text
T_eps(u)=integral_u^infinity exp(i eps pi v^2/2)dv.
```

The excluded upper half-plane has normalized tangent profile

```text
U_0=chi T_sign(Delta)(q_B*/sqrt(|Delta|))/(1+i),

U_1=-c_m chi exp[i pi q_B*^2/(2Delta)]
     /[sqrt(|Delta|) i pi(1+i)].                        (BT3)
```

Relative to the sharp classical step at `621|622`, the exact tangent-model
residual is

```text
R_B,m=1_(m<=621)-(U_0+U_1).                             (BT4)
```

No integer mode hits the tangent degeneracy.  The nearest is `m=257`, with

```text
|Delta_B|={c['minimum_integer_abs_Delta_B_ball']}>0.0009.
```

Rigorous Arb summation of (BT4) over all positive lower-branch modes
`1..39894`, using the exact Riemann--Siegel theta carrier, gives

```text
2 Re sum carrier*R_B,m
 ={c['grouped_two_real_ball']},

(sqrt(2)/pi) times this
 ={c['grouped_physical_two_real_ball']}.                (BT5)
```

Both exceed their respective targets `0.000019` and `0.0000086`; the
normalized ratio is

```text
{c['grouped_to_normalized_target_ratio_ball']}.
```

The mechanism is instructive.  If the `c_m q` channel is illegally omitted,
the same signed diagnostic would be

```text
{c['scalar_only_two_real_ball']},                       (BT6)
```

which lies below the normalized target.  The retained q-density contributes

```text
{c['q_density_two_real_ball']}.                         (BT7)
```

Thus dropping the sine/current term creates a false pass.  The complete
termwise triangle is `{c['termwise_complex_triangle_ball']}`, so the result
also depends on substantial signed cancellation.

This is a barrier only to promoting the tangent crossing by itself.  It is
not a lower bound for the exact half-Kummer error: nonlinear Morse amplitude,
the `m<=0` and `m>=39895` nonstationary completion, the endpoint half-current,
and the characteristic A/half-boundary chart remain outside (BT4).  They must
be assembled before a final sign or magnitude is claimed.

Pi provenance: `pi` in (BT1)--(BT7) comes from the equation-(9) Kummer
quadratic, the Fourier-Poisson character, Gaussian completion, and
Riemann--Siegel normalization.  No fitted constant is used.

Proof boundary: exact tangent half-plane transform and a saved-height finite
aggregate of that model only.  No nonlinear B-crossing remainder, complete
symmetric half-domain completion, A-fold splice, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["half_reflection"]["decision"]["half_domain_has_one_positive_stationary_branch"] is True, "half-reflection dependency drift")
    require(dependencies["universal_morse"]["decision"]["grouped_boundary_fresnel_current_preserved"] is True, "grouped-current dependency drift")
    require(dependencies["Gamma_bulk"]["decision"]["ordinary_full_line_bulk_closed_at_saved_height"] is True, "Gamma dependency drift")
    require(dependencies["symmetric_residual"]["decision"]["symmetric_residual_limit_proved"] is True, "residual dependency drift")

    artifact = {
        "kind": STEM,
        "status": "exact_B_crossing_tangent_half_plane_profile_aggregated_and_grouped_target_barrier_proved",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "certificate": aggregate_certificate(),
        "decision": {
            "tangent_half_plane_transform_proved": True,
            "scalar_and_q_density_channels_retained": True,
            "integer_tangent_defect_nonzero_on_1_to_39894": True,
            "scalar_only_counterfactual_below_normalized_target": True,
            "grouped_tangent_profile_below_normalized_target": False,
            "grouped_tangent_profile_below_physical_target": False,
            "dropping_q_density_channel_admissible": False,
            "tangent_profile_is_exact_half_Kummer_residual": False,
            "nonlinear_B_crossing_remainder_bounded": False,
            "symmetric_nonstationary_completion_bounded": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Compose the tangent B profile with the m<=0 and m>=39895 half-domain nonstationary Poisson completion and endpoint half-current before estimating nonlinear Morse corrections; test whether their exact x=1/2 boundary currents cancel the q-density excess in (BT7).",
        "proof_boundary": "Exact tangent half-plane transform and saved-height finite aggregate only. No nonlinear B-crossing remainder, symmetric half-domain completion, A-fold splice, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified grouped B tangent-profile barrier: "
        f"two-real={arb(artifact['certificate']['grouped_two_real_ball'])}, "
        f"scalar-only={arb(artifact['certificate']['scalar_only_two_real_ball'])}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
