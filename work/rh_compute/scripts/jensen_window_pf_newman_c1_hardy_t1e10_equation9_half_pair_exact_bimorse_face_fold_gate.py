#!/usr/bin/env python3
"""Certify the exact bi-Morse triangle chart and its face-fold defect."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_exact_bimorse_face_fold_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "universal_morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
    "fold_atlas": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_continuous_height_gate.json",
}

T = 10_000_000_000
A = 159_577
B = 5_122_421


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
    u, p, m, endpoint, t = sp.symbols("u p m D t", positive=True, real=True)
    z_star = 2 * m / endpoint
    z = z_star * u
    endpoint_phase = sp.pi * m**2 / z + sp.pi * endpoint**2 * z / 4
    endpoint_phase_star = sp.pi * m * endpoint
    p_squared = u + 1 / u - 2
    require(
        sp.simplify(endpoint_phase - endpoint_phase_star - sp.pi * m * endpoint * p_squared / 2) == 0,
        "endpoint global Morse identity failed",
    )
    require(sp.simplify((sp.sqrt(u) - 1 / sp.sqrt(u)) ** 2 - p_squared) == 0, "signed endpoint coordinate failed")

    inverse_root = (p + sp.sqrt(p**2 + 4)) / 2
    inverse_u = inverse_root**2
    require(sp.simplify(sp.sqrt(inverse_u) - 1 / sp.sqrt(inverse_u) - p) == 0, "endpoint inverse failed")

    alpha_m = 2 * m + t / (sp.pi * m)
    x_m = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    u_at_outer_saddle = sp.simplify(endpoint * x_m / (2 * m))
    require(sp.simplify(u_at_outer_saddle - endpoint / alpha_m) == 0, "face detuning ratio failed")
    p_at_outer_saddle = sp.sqrt(endpoint / alpha_m) - sp.sqrt(alpha_m / endpoint)
    require(
        sp.simplify(p_at_outer_saddle - (endpoint - alpha_m) / sp.sqrt(endpoint * alpha_m)) == 0,
        "face detuning identity failed",
    )

    # Differentiate the exact boundary p_D(s)=p(D*x(s)/(2m)).  The known
    # logistic derivative at the outer saddle is dx/ds=-r*x_m^2.
    x = sp.symbols("x", positive=True, real=True)
    ratio = endpoint * x / (2 * m)
    p_boundary = sp.sqrt(ratio) - 1 / sp.sqrt(ratio)
    r = t / (2 * sp.pi * m**2)
    dx_ds_at_zero = -r * x_m**2
    slope = sp.simplify(sp.diff(p_boundary, x).subs(x, x_m) * dx_ds_at_zero)
    expected_slope = -t * (endpoint + alpha_m) / (
        2 * sp.pi * m * sp.sqrt(endpoint) * alpha_m ** sp.Rational(3, 2)
    )
    require(sp.simplify(slope - expected_slope) == 0, "face slope failed")

    standardized_slope_squared = sp.simplify(slope**2 * (2 * sp.pi * m * endpoint / t))
    expected_standardized = t * (endpoint + alpha_m) ** 2 / (2 * sp.pi * m * alpha_m**3)
    require(sp.simplify(standardized_slope_squared - expected_standardized) == 0, "standardized slope failed")

    event_t = sp.pi * m * (endpoint - 2 * m)
    event_alpha = sp.simplify(alpha_m.subs(t, event_t))
    require(sp.simplify(event_alpha - endpoint) == 0, "event alpha failed")
    event_slope_squared = sp.simplify(standardized_slope_squared.subs(t, event_t))
    require(sp.simplify(event_slope_squared - (2 - 4 * m / endpoint)) == 0, "event slope defect failed")
    null_defect = sp.simplify(event_slope_squared - 1)
    require(sp.simplify(null_defect - (1 - 4 * m / endpoint)) == 0, "face-fold defect failed")

    half_corner_t = sp.solve(
        [sp.Eq(endpoint, 4 * m), sp.Eq(t, event_t)],
        (m, t),
        dict=True,
    )[0][t]
    require(sp.simplify(half_corner_t - sp.pi * endpoint**2 / 8) == 0, "half-corner height failed")

    return {
        "endpoint_ratio": "u=z/z_D, z_D=2m/D",
        "endpoint_Morse_coordinate": "p=sqrt(u)-1/sqrt(u)",
        "endpoint_inverse": "u(p)=[(p+sqrt(p^2+4))/2]^2",
        "endpoint_global_phase": "phi_(m,D)(z)=pi*m*D+(pi*m*D/2)p^2",
        "outer_global_phase": "psi_m(x)=psi_m(x_m)-(t/4)s^2",
        "standard_coordinates": "P=sqrt(pi*m*D)p, S=sqrt(t/2)s, so phase defect=(P^2-S^2)/2",
        "exact_boundary": "P<P_D(S), obtained from p_D(s)=sqrt[D*x(s)/(2m)]-sqrt[2m/(D*x(s))]",
        "saddle_face_detuning": "p_D(0)=(D-alpha_m)/sqrt(D*alpha_m)",
        "boundary_slope": "p_D'(0)=-t(D+alpha_m)/(2*pi*m*sqrt(D)*alpha_m^(3/2))",
        "standard_slope_squared": "rho_D^2=t(D+alpha_m)^2/(2*pi*m*alpha_m^3)",
        "event_face_defect": "At t=pi*m(D-2m), alpha_m=D and rho_D^2-1=1-4m/D.",
        "null_corner": "The face is null exactly at m=D/4; simultaneous event height is t=pi*D^2/8, equivalently D=sqrt(8t/pi).",
        "interpretation": "The phase is globally quadratic in both variables. Only amplitude transport and curvature of the exact boundary remain; the A fold is the null face of the same chart, while the B face is noncharacteristic.",
    }


def finite_height_ledger() -> dict[str, Any]:
    mp.mp.dps = 80
    t, pi = mp.mpf(T), mp.pi

    def alpha(mode: int) -> mp.mpf:
        value = mp.mpf(mode)
        return 2 * value + t / (pi * value)

    def row(endpoint: int, mode: int) -> dict[str, str | int]:
        D, m = mp.mpf(endpoint), mp.mpf(mode)
        al = alpha(mode)
        p0 = (D - al) / mp.sqrt(D * al)
        rho2 = t * (D + al) ** 2 / (2 * pi * m * al**3)
        event_defect = 1 - 4 * m / D
        return {
            "endpoint": endpoint,
            "mode": mode,
            "alpha_m": mp.nstr(al, 65),
            "p_face_at_outer_saddle": mp.nstr(p0, 65),
            "saved_height_standard_slope_squared": mp.nstr(rho2, 65),
            "event_face_fold_defect_1_minus_4m_over_D": mp.nstr(event_defect, 65),
        }

    rows = [
        row(B, 621),
        row(B, 622),
        row(A, 39_695),
        row(A, 39_852),
        row(A, 39_853),
        row(A, 39_894),
        row(A, 39_895),
        row(A, 40_093),
    ]
    by_key = {(item["endpoint"], item["mode"]): item for item in rows}
    b_defect = mp.mpf(by_key[(B, 621)]["event_face_fold_defect_1_minus_4m_over_D"])
    a_near = mp.mpf(by_key[(A, 39_894)]["event_face_fold_defect_1_minus_4m_over_D"])
    require(b_defect > mp.mpf("0.999"), "B face unexpectedly characteristic")
    require(abs(a_near - mp.mpf(1) / A) < mp.mpf("1e-70"), "A near-null defect drift")
    require(mp.mpf(by_key[(B, 621)]["p_face_at_outer_saddle"]) < 0, "B 621 face sign drift")
    require(mp.mpf(by_key[(B, 622)]["p_face_at_outer_saddle"]) > 0, "B 622 face sign drift")
    require(mp.mpf(by_key[(A, 39_852)]["p_face_at_outer_saddle"]) < 0, "A 39852 face sign drift")
    require(mp.mpf(by_key[(A, 39_853)]["p_face_at_outer_saddle"]) > 0, "A 39853 face sign drift")
    return {
        "height": T,
        "rows": rows,
        "B_event_noncharacteristic_lower_bound": "0.999",
        "A_mode_39894_exact_event_defect": f"1/{A}",
        "A_half_corner_mode": mp.nstr(mp.mpf(A) / 4, 65),
        "A_exact_null_height_pi_A_squared_over_8": mp.nstr(pi * mp.mpf(A) ** 2 / 8, 65),
        "saved_height_minus_A_null_height": mp.nstr(t - pi * mp.mpf(A) ** 2 / 8, 65),
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = {(row["endpoint"], row["mode"]): row for row in artifact["finite_height_ledger"]["rows"]}
    return f"""# Exact bi-Morse face and fold chart

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the uniform triangle estimate

The separable triangle phase has a second global Morse coordinate.  For the
endpoint phase put

```text
z_D=2m/D,  u=z/z_D,
p=sqrt(u)-1/sqrt(u).                                  (BM1)
```

Then, globally for `z>0`,

```text
phi_(m,D)(z)=pi*m*D+(pi*m*D/2)p^2.                   (BM2)
```

Together with the exact outer coordinate `s`,

```text
psi_m(x)=psi_m(x_m)-(t/4)s^2,                         (BM3)
```

the standardized variables

```text
P=sqrt(pi*m*D)p,   S=sqrt(t/2)s                       (BM4)
```

make the complete phase defect exactly `(P^2-S^2)/2`.  There is no phase
Taylor remainder anywhere in this bi-Morse chart.

The triangle `z<x` becomes the exact curved face

```text
P<P_D(S),
p_D(s)=sqrt[D*x(s)/(2m)]-sqrt[2m/(D*x(s))].           (BM5)
```

At the outer saddle,

```text
p_D(0)=(D-alpha_m)/sqrt(D*alpha_m),                   (BM6)

rho_D^2=[dP_D/dS(0)]^2
 =t(D+alpha_m)^2/(2*pi*m*alpha_m^3).                  (BM7)
```

At an exact endpoint event `t=pi*m(D-2m)`, `alpha_m=D`, and the quadratic
defect induced along the tangent face is

```text
rho_D^2-1=1-4m/D.                                    (BM8)
```

This explains the two qualitatively different endpoints.  At the B crossing,

```text
1-4*621/B={rows[(B,621)]['event_face_fold_defect_1_minus_4m_over_D']}>0.999,
```

so the B face is safely noncharacteristic.  For the A atlas near mode 39894,

```text
1-4*39894/A={rows[(A,39894)]['event_face_fold_defect_1_minus_4m_over_D']}
            =1/{A}.                                   (BM9)
```

The null face occurs exactly at `m=D/4`.  Requiring it to be an endpoint
event also gives

```text
t=pi*D^2/8,  D=sqrt(8t/pi).                           (BM10)
```

Thus the A characteristic fold and the `x=1/2` corner are not neighboring
accidents: they are the same null face in exact bi-Morse coordinates.  The
saved `A` differs slightly from `sqrt(8t/pi)`, producing the small but nonzero
detuning recorded by the interval ledger.

The prospective triangle theorem is now simpler.  Its phase is exact; it
must bound only the transformed amplitude, the curvature of (BM5), and the
finite face/corner truncations.  An ordinary Fresnel face theorem applies at
B.  A null-face fold theorem, compatible with the certified 399-event atlas,
is required at A.

Pi provenance: `pi` is inherited from the exact equation-(9) endpoint and
outer phases.  No fitted or geometric constant is introduced.

Proof boundary: exact bi-Morse coordinates, face detuning, tangent slope, and
fold-defect identity only.  No uniform amplitude/boundary-curvature bound,
complete B face estimate, A-fold splice, complete paired residual,
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion follows.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(dependencies["triangle"].get("passed") is True, "triangle dependency is not passed")
    require(dependencies["universal_morse"].get("passed") is True, "Morse dependency is not passed")
    require(
        dependencies["triangle"]["decision"]["pair_Volterra_integral_reduced_to_separable_triangle"] is True,
        "triangle dependency drift",
    )
    require(
        dependencies["universal_morse"]["decision"]["universal_global_morse_phase_proved"] is True,
        "outer Morse dependency drift",
    )
    fold_atlas = dependencies["fold_atlas"]
    require(
        fold_atlas.get("status") == "beta_minus_4_continuous_height_certified_on_all_399_exact_event_cells",
        "fold atlas status drift",
    )
    require(fold_atlas["certificate"]["event_count"] == 399, "fold atlas event-count drift")
    require(
        "2.152" in fold_atlas["certificate"]["maximum_normalized_absolute_upper"],
        "fold atlas maximum-bound drift",
    )

    artifact = {
        "kind": STEM,
        "status": "exact_global_bimorse_triangle_face_and_null_fold_identity_proved_uniform_estimate_open",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "finite_height_ledger": finite_height_ledger(),
        "decision": {
            "endpoint_phase_has_exact_global_Morse_coordinate": True,
            "pair_triangle_phase_is_exactly_hyperbolic_quadratic": True,
            "endpoint_event_is_exact_face_crossing": True,
            "B_face_is_noncharacteristic": True,
            "A_fold_and_half_boundary_are_same_null_face": True,
            "uniform_transformed_amplitude_bound_proved": False,
            "uniform_curved_face_remainder_proved": False,
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
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "sympy_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Prove the ordinary noncharacteristic face estimate first at B using the exact hyperbolic Gaussian phase and curved boundary P_B(S), retaining the complete transformed driver amplitude. Then express the A null-face correction in the same coordinates and compare it directly with the certified beta^-4 fold atlas.",
        "proof_boundary": "Exact bi-Morse coordinates and face-fold geometry only. No uniform transformed-amplitude or curved-face bound, complete B crossing, A-fold splice, complete paired residual, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact bi-Morse triangle face: B ordinary, A null-fold", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
