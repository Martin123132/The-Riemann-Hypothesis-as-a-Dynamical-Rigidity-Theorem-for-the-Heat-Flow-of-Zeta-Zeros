#!/usr/bin/env python3
"""Certify the exact affine-wedge remainder decomposition of the pair triangle."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
import sys
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "bi_Morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_exact_bimorse_face_fold_gate.json",
    "wedge": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_corner_kernel_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
B = 5_122_421
KEYS = ((B, 621), (B, 622), (A, 39_852), (A, 39_853), (A, 39_894), (A, 39_895))


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


def symbolic_certificate() -> dict[str, str]:
    u, v, sigma = sp.symbols("u v sigma", positive=True, real=True)
    m, endpoint, t = sp.symbols("m D t", positive=True, real=True)
    p = sp.sqrt(u) - 1 / sp.sqrt(u)
    z = 2 * m * u / endpoint
    P = sp.sqrt(sp.pi * m * endpoint) * p
    r = t / (2 * sp.pi * m**2)
    x = 1 / (1 + r * v)

    dz_dP = sp.diff(z, u) / sp.diff(P, u)
    dv_dsigma = v * sigma / (v - 1)
    minus_dx_dS = -sp.diff(x, v) * dv_dsigma * sp.sqrt(2 / t)
    transformed_jacobian = sp.simplify(
        z ** (-sp.Rational(1, 2))
        * (2 + sp.I * sp.pi * endpoint**2 * z)
        * x ** (-sp.Rational(7, 4))
        * (1 - x) ** (-sp.Rational(1, 4))
        * dz_dP
        * minus_dx_dS
    )

    K = sp.pi * endpoint * m
    W0 = 4 * (K - sp.I) * r ** sp.Rational(3, 4) / (endpoint * sp.sqrt(sp.pi * t))
    endpoint_factor = 2 * u * (K * u - sp.I) / ((u + 1) * (K - sp.I))
    outer_factor = v ** sp.Rational(3, 4) * sigma / (v - 1)
    normalized_weight = sp.simplify(transformed_jacobian / sp.I)
    require(sp.simplify(normalized_weight - W0 * endpoint_factor * outer_factor) == 0, "amplitude factorization failed")
    require(sp.simplify(endpoint_factor.subs(u, 1) - 1) == 0, "endpoint normalization failed")

    p_local = sp.symbols("p_local", real=True)
    root = (p_local + sp.sqrt(p_local**2 + 4)) / 2
    u_of_p = root**2
    endpoint_in_p = endpoint_factor.subs(u, u_of_p)
    lambda_P = (3 * K - sp.I) / (2 * (K - sp.I) * sp.sqrt(sp.pi * m * endpoint))
    require(sp.simplify(sp.diff(endpoint_in_p, p_local).subs(p_local, 0) / sp.sqrt(sp.pi * m * endpoint) - lambda_P) == 0, "P slope failed")

    v_series = 1 + sigma + sigma**2 / 3 + sigma**3 / 36 - sigma**4 / 270 + sigma**5 / 4320
    morse_defect = sp.series(2 * (v_series - 1 - sp.log(v_series)) - sigma**2, sigma, 0, 6).removeO()
    require(sp.expand(morse_defect) == 0, "logistic Morse series failed")
    outer_series = sp.series(v_series ** sp.Rational(3, 4) * sigma / (v_series - 1), sigma, 0, 3).removeO()
    require(sp.simplify(outer_series - (1 + 5 * sigma / 12 - sigma**2 / 96)) == 0, "outer amplitude series failed")
    mu_S = sp.Rational(5, 12) * sp.sqrt(2 / t)

    indicator_exact, indicator_tangent, exact_amplitude, affine_amplitude = sp.symbols(
        "chi_exact chi_tangent A_exact A_affine"
    )
    decomposition = sp.expand(
        indicator_exact * (exact_amplitude - affine_amplitude)
        + (indicator_exact - indicator_tangent) * affine_amplitude
        - (indicator_exact * exact_amplitude - indicator_tangent * affine_amplitude)
    )
    require(decomposition == 0, "signed remainder decomposition failed")

    C, C_a, C_h, lam, mu, rho = sp.symbols("C C_a C_h lambda mu rho")
    M_P = -sp.I * C_a
    M_S = sp.I * (C_h - rho * C_a)
    affine_wedge = sp.expand(C + lam * M_P + mu * M_S)
    require(
        sp.simplify(affine_wedge - (C + sp.I * mu * C_h - sp.I * (lam + mu * rho) * C_a)) == 0,
        "affine wedge moment reduction failed",
    )

    P_of_v = sp.sqrt(sp.pi * m * endpoint) * (
        sp.sqrt(endpoint * x / (2 * m)) - sp.sqrt(2 * m / (endpoint * x))
    )
    P_v = sp.diff(P_of_v, v)
    face_curvature = sp.factor(2 * (sp.diff(P_of_v, v, 2) + sp.Rational(2, 3) * P_v).subs(v, 1) / t)
    expected_curvature = -(
        8 * sp.pi**2 * endpoint * m**3
        - 5 * sp.pi * endpoint * m * t
        + 16 * sp.pi**2 * m**4
        + 10 * sp.pi * m**2 * t
        + t**2
    ) / (6 * (2 * sp.pi * m**2 + t) ** sp.Rational(5, 2))
    require(sp.simplify(face_curvature - expected_curvature) == 0, "face curvature failed")

    return {
        "exact_endpoint_contribution": "eps_D exp(i Phi_D,m) W0_D,m (2pi)^(-1) int_(S>h) int_(P<P_D(S)) A_D,m(P,S) exp(i(P^2-S^2)/2)dP dS",
        "weight_prefactor": "W0=4(pi*D*m-i)r^(3/4)/(D sqrt(pi*t)), r=t/(2pi*m^2)",
        "normalized_amplitude": "A_D,m(P,S)=E_D,m(u(P)) O(v(S))",
        "endpoint_factor": "E=2u(pi*D*m*u-i)/[(u+1)(pi*D*m-i)]",
        "outer_factor": "O=v^(3/4)s/(v-1)",
        "affine_amplitude": "A_aff=1+lambda_P P+mu_S S",
        "lambda_P": "(3pi*D*m-i)/[2(pi*D*m-i)sqrt(pi*m*D)]",
        "mu_S": "(5/12)sqrt(2/t)",
        "signed_remainder": "chi_exact*A_exact-chi_tangent*A_aff=chi_exact*(A_exact-A_aff)+(chi_exact-chi_tangent)*A_aff",
        "amplitude_remainder": "R_amp=(2pi)^(-1) int_(S>h)int_(P<P_D(S))(A_exact-A_aff)exp(i(P^2-S^2)/2)dP dS",
        "face_remainder": "R_face=(2pi)^(-1) int_(S>h)int_(a+rho*S)^(P_D(S)) A_aff exp(i(P^2-S^2)/2)dP dS, with oriented inner integral",
        "affine_wedge": "C_aff=C+i*mu_S*partial_h C-i(lambda_P+mu_S*rho)*partial_a C",
        "face_curvature": "P_D''(0)=-[8pi^2 Dm^3-5pi Dmt+16pi^2m^4+10pi m^2t+t^2]/[6(2pi m^2+t)^(5/2)]",
        "interpretation": "The phase has no remainder. The exact defect from the affine tangent wedge is exactly one second-order local amplitude term plus one oriented curved-face strip; no absolute values have yet been taken.",
    }


def geometry_certificate(wedge: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = 100
    pi, t = arb.pi(), arb(HEIGHT)
    wedge_rows = {
        (row["endpoint"], row["mode"]): row
        for row in wedge["geometry_certificate"]["rows"]
    }
    mu = arb(5) / 12 * (arb(2) / t).sqrt()
    rows: list[dict[str, Any]] = []
    for endpoint_int, mode_int in KEYS:
        endpoint, mode = arb(endpoint_int), arb(mode_int)
        denominator = 2 * pi * mode**2 + t
        alpha = 2 * mode + t / (pi * mode)
        a = (pi * mode / alpha).sqrt() * (endpoint - alpha)
        rho = -(2 * t).sqrt() * (pi * endpoint * mode + denominator) / (2 * denominator ** arb("1.5"))
        c = rho**2 - 1
        curvature = -(
            8 * pi**2 * endpoint * mode**3
            - 5 * pi * endpoint * mode * t
            + 16 * pi**2 * mode**4
            + 10 * pi * mode**2 * t
            + t**2
        ) / (6 * denominator ** arb("2.5"))
        K = pi * endpoint * mode
        scale = (pi * mode * endpoint).sqrt()
        lambda_re = (3 * K**2 + 1) / (2 * (K**2 + 1) * scale)
        lambda_im = K / ((K**2 + 1) * scale)
        r = t / (2 * pi * mode**2)
        w0_scale = 4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt())
        rows.append(
            {
                "endpoint": endpoint_int,
                "mode": mode_int,
                "a_face_detuning_ball": a.str(80, more=True),
                "rho_signed_ball": rho.str(80, more=True),
                "c_characteristic_defect_ball": c.str(80, more=True),
                "face_second_derivative_ball": curvature.str(80, more=True),
                "lambda_P_real_ball": lambda_re.str(80, more=True),
                "lambda_P_imag_ball": lambda_im.str(80, more=True),
                "mu_S_ball": mu.str(80, more=True),
                "W0_real_ball": (w0_scale * K).str(80, more=True),
                "W0_imag_ball": (-w0_scale).str(80, more=True),
                "h_half_boundary_ball": wedge_rows[(endpoint_int, mode_int)]["h_half_boundary_ball"],
            }
        )

    for row in rows:
        require(abs(arb(row["face_second_derivative_ball"])).upper() < arb("7e-6"), "local face curvature exceeds audit cap")
        require(arb(row["lambda_P_real_ball"]).upper() < arb("2e-5"), "local P slope exceeds audit cap")
        require(abs(arb(row["lambda_P_imag_ball"])).upper() < arb("1e-14"), "local P imaginary slope exceeds audit cap")
    require(mu.upper() < arb("6e-6"), "outer S slope exceeds audit cap")
    return {
        "height": HEIGHT,
        "rows": rows,
        "uniform_local_coefficient_audit": {
            "abs_face_second_derivative_below": "7e-6",
            "lambda_P_real_below": "2e-5",
            "abs_lambda_P_imag_below": "1e-14",
            "mu_S_below": "6e-6",
        },
        "warning": "These are pointwise tangent coefficients at S=P=0, not bounds for either global remainder integral.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = {(row["endpoint"], row["mode"]): row for row in artifact["geometry_certificate"]["rows"]}
    return f"""# Exact affine-wedge remainder identity

Date: 2026-08-13

Status: exact transformed-amplitude and signed-domain decomposition; not a
bound for either remainder integral

In the exact bi-Morse variables, including both Jacobians, the endpoint-`D`
pair triangle is

```text
eps_D exp(i Phi_D,m) W0_D,m (2pi)^(-1)
 int_(S>h) int_(P<P_D(S)) A_D,m(P,S)
 exp(i(P^2-S^2)/2)dP dS,                              (AR1)

W0_D,m=4(pi Dm-i)r^(3/4)/(D sqrt(pi t)),
r=t/(2pi m^2).                                        (AR2)
```

The normalized transformed amplitude factorizes exactly:

```text
A_D,m(P,S)=E_D,m(u(P)) O(v(S)),
E_D,m(u)=2u(pi Dm u-i)/[(u+1)(pi Dm-i)],
O(v)=v^(3/4)s/(v-1).                                  (AR3)
```

Both factors equal one at the joint saddle.  Their complete affine tangent is

```text
A_aff(P,S)=1+lambda_P P+mu_S S,
lambda_P=(3pi Dm-i)/[2(pi Dm-i)sqrt(pi mD)],
mu_S=(5/12)sqrt(2/t).                                 (AR4)
```

This retains the mandatory endpoint-density current rather than silently
dropping it.  If `C=C_rho(a,h)` is the wedge from Section 11.427, its affine
moment is still exact:

```text
C_aff=C+i mu_S partial_h C
        -i(lambda_P+mu_S rho)partial_a C.             (AR5)
```

The exact-minus-affine-tangent defect now splits before norms as

```text
R_exact-aff=R_amp+R_face,                              (AR6)

R_amp=(2pi)^(-1) int_(S>h) int_(P<P_D(S))
      [A_D,m(P,S)-A_aff(P,S)]exp(i(P^2-S^2)/2)dP dS,

R_face=(2pi)^(-1) int_(S>h) int_(a+rho S)^(P_D(S))
       A_aff(P,S)exp(i(P^2-S^2)/2)dP dS.              (AR7)
```

The inner integral in `R_face` is oriented, so (AR6) remains valid when the
curved face crosses its tangent.  There is no phase remainder.

At the six saved-height transition modes, interval arithmetic gives

```text
|P_D''(0)|<7e-6,   Re(lambda_P)<2e-5,
|Im(lambda_P)|<1e-14,   mu_S<6e-6.                    (AR8)
```

For orientation, the near-corner values are

```text
A,39894: P_D''(0)={rows[(A,39894)]['face_second_derivative_ball']},
A,39895: P_D''(0)={rows[(A,39895)]['face_second_derivative_ball']}.
```

These small tangent coefficients are promising but are not global integral
bounds.  The next certificate must put a finite local box around the A
corner, bound the second-order amplitude and face strip there, and dispatch
its complement by phase-adapted integration by parts.  The B half-boundary
tail can be treated separately because Section 11.427 puts it over 270598
standardized units away.

Pi provenance: every `pi` in (AR1)--(AR8) comes from equation (9), the exact
Morse scalings, and the Fourier Gaussian.  No fitted constant is introduced.

Proof boundary: exact coordinate Jacobians, amplitude factorization, affine
wedge moments, signed remainder identity, and local tangent coefficients
only.  No global `R_amp` or `R_face` bound, `R_Dir` estimate, complete
`Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["triangle"]["decision"]["pair_Volterra_integral_reduced_to_separable_triangle"] is True, "triangle dependency drift")
    require(dependencies["bi_Morse"]["decision"]["pair_triangle_phase_is_exactly_hyperbolic_quadratic"] is True, "bi-Morse dependency drift")
    require(dependencies["wedge"]["decision"]["canonical_wedge_boundary_derivative_system_proved"] is True, "wedge dependency drift")

    artifact = {
        "kind": STEM,
        "status": "exact_affine_wedge_and_signed_curved_face_amplitude_remainder_identity_certified_bounds_open",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "geometry_certificate": geometry_certificate(dependencies["wedge"]),
        "decision": {
            "exact_transformed_amplitude_factorized": True,
            "mandatory_affine_P_current_retained": True,
            "outer_affine_S_current_retained": True,
            "affine_wedge_reduced_to_C_and_boundary_derivatives": True,
            "exact_defect_split_into_signed_amplitude_and_face_remainders": True,
            "phase_remainder_identically_zero": True,
            "global_amplitude_remainder_bound_proved": False,
            "global_curved_face_remainder_bound_proved": False,
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
            "sympy_threads": 1,
            "arb_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Choose a finite A-corner box in (P,S), certify second derivatives of E(u(P))O(v(S)) and P_A(S) throughout it, and bound the two signed remainders in (AR7) without dividing by c=rho^2-1. Treat the complement by hyperbolic phase integration by parts. Reuse the remote Fresnel-tail machinery for B and sum only after the common Gamma carrier is restored.",
        "proof_boundary": "Exact transformed-amplitude factorization, affine wedge moments, signed remainder identity, and saved-height tangent coefficients only. No global amplitude or curved-face remainder bound, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact affine-wedge amplitude and curved-face remainder identity", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
