#!/usr/bin/env python3
"""Rigorously evaluate the null-coordinate affine wedge at the A corner."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_null_coordinate_affine_corner_evaluation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "null_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_null_coordinate_unfolding_gate.json",
    "affine_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
MODES = (39_894, 39_895)
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


def complex_record(value: acb, digits: int = 80) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def right_ray(A_quad: arb, B_linear: arb, r0: arb, cutoff: arb) -> tuple[acb, arb]:
    require(not A_quad.contains(0), "quadratic sign is unresolved")
    sign = 1 if A_quad.lower() > 0 else -1
    magnitude = abs(A_quad)
    sqrt_two = arb(2).sqrt()
    direction = acb(1, sign) / sqrt_two
    imaginary = acb(0, 1)

    def integrand(y: acb, analytic: bool) -> acb:
        del analytic
        R = acb(r0) + direction * y
        return direction * (imaginary * (A_quad * R**2 + B_linear * R)).exp() / R

    value = acb.integral(
        integrand,
        arb(0),
        cutoff,
        abs_tol=arb("1e-30"),
        rel_tol=arb("1e-30"),
        eval_limit=400_000,
        depth_limit=55,
    )
    linear_growth = -arb(sign) * (2 * A_quad * r0 + B_linear) / sqrt_two
    decay_rate = 2 * magnitude * cutoff - linear_growth
    require(decay_rate.lower() > 0 and (cutoff - r0).lower() > 0, "right-ray tail is not decreasing")
    tail = (-magnitude * cutoff**2 + linear_growth * cutoff).exp() / (
        decay_rate * (cutoff - r0)
    )
    return value + acb(arb(0, tail), arb(0, tail)), tail


def full_lower_indented(A_quad: arb, B_linear: arb) -> acb:
    require(not A_quad.contains(0), "quadratic sign is unresolved")
    sign = 1 if A_quad.lower() > 0 else -1
    direction = acb(1, sign) / arb(2).sqrt()
    argument = direction * B_linear / (2 * abs(A_quad).sqrt())
    return acb(0, 1) * arb.pi() * (1 + argument.erf())


def canonical_wedge(a: arb, rho: arb, h: arb, cutoff: arb) -> tuple[acb, dict[str, arb]]:
    delta = 1 + rho
    require(not delta.contains(0), "delta sign is unresolved")
    r0 = a + delta * h
    require(r0.lower() > 0, "the saved corner path requires positive r0")
    A1, B1 = arb("0.5"), -h
    A2, B2 = arb("0.5") - 1 / delta, a / delta
    ray1, tail1 = right_ray(A1, B1, r0, cutoff)
    ray2, tail2 = right_ray(A2, B2, r0, cutoff)
    full1 = full_lower_indented(A1, B1)
    if delta.upper() < 0:
        value = (full1 - ray1 - full_lower_indented(A2, B2) + ray2) / (2 * arb.pi() * acb(0, 1))
        branch = "delta_negative_two_lower_indented_integrals_half_jumps_cancel"
    else:
        value = (full1 - ray1 + ray2) / (2 * arb.pi() * acb(0, 1))
        branch = "delta_positive_lower_indentation_supplies_half_jump"
    return value, {
        "delta": delta,
        "r0": r0,
        "A2": A2,
        "B2": B2,
        "ray1_tail": tail1,
        "ray2_tail": tail2,
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


def evaluate_rows(affine: dict[str, Any], cutoff: int) -> list[dict[str, Any]]:
    ctx.dps = 100
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    cutoff_ball = arb(cutoff)
    source_rows = {
        (row["endpoint"], row["mode"]): row
        for row in affine["geometry_certificate"]["rows"]
    }
    rows: list[dict[str, Any]] = []
    for mode_int in MODES:
        source = source_rows[(A, mode_int)]
        a = arb(source["a_face_detuning_ball"])
        rho = arb(source["rho_signed_ball"])
        h = arb(source["h_half_boundary_ball"])
        C, metadata = canonical_wedge(a, rho, h, cutoff_ball)
        C_a, C_h = wedge_derivatives(a, rho, h)

        mode = arb(mode_int)
        K = pi * endpoint * mode
        scale = (pi * mode * endpoint).sqrt()
        lambda_P = acb(
            (3 * K**2 + 1) / (2 * (K**2 + 1) * scale),
            K / ((K**2 + 1) * scale),
        )
        mu_S = arb(5) / 12 * (arb(2) / t).sqrt()
        C_aff = C + acb(0, 1) * mu_S * C_h - acb(0, 1) * (lambda_P + mu_S * rho) * C_a
        r = t / (2 * pi * mode**2)
        W0 = acb(K, -1) * (4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt()))
        weighted = W0 * C_aff
        rows.append(
            {
                "endpoint": A,
                "mode": mode_int,
                "branch": metadata["branch"],
                "delta_one_plus_rho_ball": metadata["delta"].str(80, more=True),
                "r0_ball": metadata["r0"].str(80, more=True),
                "second_quadratic_A2_ball": metadata["A2"].str(80, more=True),
                "second_linear_B2_ball": metadata["B2"].str(80, more=True),
                "ray_cutoff": cutoff,
                "first_ray_tail_absolute_bound": metadata["ray1_tail"].str(80, more=True),
                "second_ray_tail_absolute_bound": metadata["ray2_tail"].str(80, more=True),
                "canonical_wedge_C": complex_record(C),
                "partial_a_C": complex_record(C_a),
                "partial_h_C": complex_record(C_h),
                "affine_wedge_C_aff": complex_record(C_aff),
                "affine_minus_scalar_wedge": complex_record(C_aff - C),
                "W0_times_affine_wedge": complex_record(weighted),
            }
        )
    require(arb(rows[0]["canonical_wedge_C"]["real_ball"]).lower() > arb("0.55"), "39894 wedge real-part drift")
    require(arb(rows[1]["canonical_wedge_C"]["real_ball"]).lower() > arb("0.22"), "39895 wedge real-part drift")
    require(abs(arb(rows[0]["affine_minus_scalar_wedge"]["real_ball"])).lower() > arb("6e-4"), "mandatory 39894 affine correction disappeared")
    require(arb(rows[0]["first_ray_tail_absolute_bound"]).upper() < arb("1e-70"), "first ray tail target failed")
    require(arb(rows[1]["first_ray_tail_absolute_bound"]).upper() < arb("1e-70"), "second-mode first ray tail target failed")
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    rows = {row["mode"]: row for row in artifact["rows"]}
    return f"""# Certified affine A-corner wedge values

Date: 2026-08-13

Status: rigorous canonical scalar and affine wedge evaluation at modes 39894
and 39895; not a curved-face, amplitude-remainder, or global residual bound

For a quadratic phase `A R^2+B R`, complete the phase before deformation.
The lower-indented full-line integral is

```text
L_full(A,B)=i pi[1+erf(e^(i sgn(A)pi/4)B/(2sqrt(|A|)))]. (CE1)
```

Subtract its right ray
`R=r0+e^(i sgn(A)pi/4)y`.  Along that ray the quadratic modulus is exactly
`exp(-|A|y^2+Ly)`, so the omitted tail has an explicit monotone Gaussian
majorant.  Combining the rays according to (NU3)--(NU4) retains the Abel
half-jump and avoids epsilon extrapolation.

The certified scalar wedges are

```text
C_39894={rows[39894]['canonical_wedge_C']['real_ball']}
       +i {rows[39894]['canonical_wedge_C']['imag_ball']},

C_39895={rows[39895]['canonical_wedge_C']['real_ball']}
       +i {rows[39895]['canonical_wedge_C']['imag_ball']}. (CE2)
```

The exact affine moments from Section 11.428 give

```text
C_aff=C+i mu_S partial_h C-i(lambda_P+mu_S rho)partial_a C. (CE3)
```

Hence

```text
C_aff,39894={rows[39894]['affine_wedge_C_aff']['real_ball']}
           +i {rows[39894]['affine_wedge_C_aff']['imag_ball']},

C_aff,39895={rows[39895]['affine_wedge_C_aff']['real_ball']}
           +i {rows[39895]['affine_wedge_C_aff']['imag_ball']}. (CE4)
```

At mode 39894 the mandatory affine correction is

```text
{rows[39894]['affine_minus_scalar_wedge']['real_ball']}
+i {rows[39894]['affine_minus_scalar_wedge']['imag_ball']}. (CE5)
```

Its real component exceeds `6e-4`, so a scalar-only corner carrier would not
be an admissible approximation at the target precision.  These are values of
the canonical tangent wedge only.  The exact paired-triangle phase carrier,
curved-face strip, transformed-amplitude remainder, all other modes, and the
global projector orientation have not been summed here.

Pi provenance: (CE1) is the exact Gaussian rotation of the equation-(9)
hyperbolic phase; (CE3) uses the Fourier Gaussian moments.  No fitted
constant or geometric surrogate is introduced.

Proof boundary: rigorous canonical scalar/affine wedge values and explicit
steepest-ray tail bounds at two saved-height modes only.  No exact curved-face
or transformed-amplitude remainder, `R_Dir` estimate, complete `Q_K-T` or
`T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["null_reduction"]["decision"]["near_null_wedge_reduced_to_one_dimensional_Abel_integrals"] is True, "null-reduction dependency drift")
    require(dependencies["affine_remainder"]["decision"]["affine_wedge_reduced_to_C_and_boundary_derivatives"] is True, "affine dependency drift")

    rows = evaluate_rows(dependencies["affine_remainder"], RAY_CUTOFF)
    artifact = {
        "kind": STEM,
        "status": "rigorous_null_coordinate_scalar_and_affine_A_corner_wedge_values_certified_exact_remainders_open",
        "passed": True,
        "rows": rows,
        "decision": {
            "epsilon_extrapolation_avoided": True,
            "Abel_half_jump_retained_by_lower_indentation": True,
            "canonical_scalar_wedge_values_certified": True,
            "canonical_affine_wedge_values_certified": True,
            "mandatory_affine_current_quantitatively_nonzero": True,
            "curved_face_remainder_bound_proved": False,
            "transformed_amplitude_remainder_bound_proved": False,
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
            "ray_cutoff": RAY_CUTOFF,
        },
        "next_obligation": "Restore eps_A exp(i Phi_A,m)W0 to the two certified affine wedges and derive the exact paired-projector orientation. Then evaluate the same characteristic carrier across the complete A transition roster with deterministic mode order. In parallel, express the curved boundary R_A(S)=P_A(S)+S as eta_A plus its cubic null unfolding and bound exact-minus-canonical only after the paired carrier is assembled.",
        "proof_boundary": "Rigorous canonical scalar and affine wedge evaluation at modes 39894 and 39895 only. No curved-face or transformed-amplitude remainder, complete A roster, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified null-coordinate affine A-corner wedge values", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
