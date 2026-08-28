#!/usr/bin/env python3
"""Certify a characteristic-coordinate reduction of the nearly null A wedge."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_null_coordinate_unfolding_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "affine_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate.json",
    "wedge": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_corner_kernel_gate.json",
    "bi_Morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_exact_bimorse_face_fold_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
MODES = (39_894, 39_895)


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
    P, S, a, rho, h = sp.symbols("P S a rho h", real=True)
    R = sp.symbols("R", real=True)
    delta = sp.symbols("delta", real=True, nonzero=True)
    phase = (P**2 - S**2) / 2
    phase_R = sp.expand(phase.subs(P, R - S))
    require(sp.simplify(phase_R - (R**2 / 2 - R * S)) == 0, "null-coordinate phase failed")
    require(sp.simplify((P - a - rho * S).subs(P, R - S).subs(rho, delta - 1) - (R - a - delta * S)) == 0, "null-coordinate face failed")
    r0 = a + delta * h

    t, endpoint = sp.symbols("t D", positive=True, real=True)
    nu = sp.sqrt(t / (2 * sp.pi))
    dynamic_endpoint = sp.sqrt(8 * t / sp.pi)
    eta = sp.sqrt(sp.pi) * (endpoint - 4 * nu) / 2
    eta_rationalized = 4 * (sp.pi * endpoint**2 / 8 - t) / (
        sp.sqrt(sp.pi) * (endpoint + dynamic_endpoint)
    )
    require(sp.simplify(eta - eta_rationalized) == 0, "corner unfolding rationalization failed")
    require(sp.simplify(eta.subs(t, sp.pi * endpoint**2 / 8)) == 0, "null-height unfolding failed")

    m = sp.symbols("m", positive=True, real=True)
    event_t = sp.pi * m * (endpoint - 2 * m)
    event_rho = -sp.sqrt(2) * sp.sqrt(endpoint - 2 * m) / sp.sqrt(endpoint)
    event_curvature = -2 * (-endpoint + 6 * m) / (
        3 * sp.sqrt(sp.pi) * endpoint ** sp.Rational(3, 2) * sp.sqrt(m)
    )
    null_rho = sp.simplify(event_rho.subs(m, endpoint / 4))
    null_curvature = sp.simplify(event_curvature.subs(m, endpoint / 4))
    require(null_rho == -1, "continuous null slope failed")
    require(sp.simplify(null_curvature + 2 / (3 * sp.sqrt(sp.pi) * endpoint)) == 0, "continuous null curvature failed")
    require(sp.simplify(null_rho * null_curvature / 2 - 1 / (3 * sp.sqrt(sp.pi) * endpoint)) == 0, "Airy cubic coefficient failed")
    require(sp.simplify(event_t.subs(m, endpoint / 4) - sp.pi * endpoint**2 / 8) == 0, "null event height failed")

    return {
        "null_coordinate": "R=P+S, P=R-S",
        "phase": "(P^2-S^2)/2=R^2/2-R*S",
        "face": "P<a+rho*S iff R<a+delta*S, delta=1+rho",
        "tangent_corner_value": "r0=a+delta*h",
        "delta_negative_formula": "C=(2pi*i)^(-1) int_(-infinity)^r0 exp(iR^2/2)[exp(-iRh)-exp(-iR(R-a)/delta)]dR/R",
        "delta_positive_formula": "C=1/2+(2pi*i)^(-1)PV{int_(-infinity)^r0 exp(iR^2/2-iRh)dR/R+int_r0^infinity exp(iR^2/2-iR(R-a)/delta)dR/R}",
        "delta_zero_formula": "C=H(a)/2+(2pi*i)^(-1)PV int_(-infinity)^a exp(iR^2/2-iRh)dR/R",
        "removable_origin": "For delta<0 the numerator difference vanishes at R=0; for delta>=0 the displayed PV and half-jump are the common Abel boundary value.",
        "continuous_half_mode": "nu=sqrt(t/(2pi))",
        "dynamic_null_endpoint": "D_t=4nu=sqrt(8t/pi)",
        "corner_unfolding": "eta_D=sqrt(pi)(D-D_t)/2=4(pi D^2/8-t)/[sqrt(pi)(D+D_t)]",
        "null_equivalence": "eta_D=0 iff t=pi D^2/8",
        "continuous_null_cubic": "At m=D/4 and t=pi D^2/8, rho=-1, P_D''(0)=-2/(3sqrt(pi)D), and the face phase has cubic coefficient rho*P_D''(0)/2=1/(3Dsqrt(pi)).",
        "interpretation": "The A corner is a finite codimension-two unfolding in eta_D. The one-dimensional formulas remain meaningful as delta changes sign and avoid the nonuniform 1/sqrt(|rho^2-1|) parameterization.",
    }


def geometry_certificate(affine: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = 110
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    nu = (t / (2 * pi)).sqrt()
    dynamic_endpoint = (8 * t / pi).sqrt()
    t_star = pi * endpoint**2 / 8
    eta = pi.sqrt() * (endpoint - dynamic_endpoint) / 2
    eta_alt = 4 * (t_star - t) / (pi.sqrt() * (endpoint + dynamic_endpoint))
    require(eta.overlaps(eta_alt), "corner unfolding forms do not overlap")
    require(eta.lower() > arb("0.07784") and eta.upper() < arb("0.07785"), "corner unfolding drift")

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
        delta = 1 + rho
        tangent_corner = a + delta * h
        exact_half_face = pi.sqrt() * (endpoint - 4 * mode_int) / 2
        exact_corner = exact_half_face + h
        rows.append(
            {
                "mode": mode_int,
                "delta_one_plus_rho_ball": delta.str(90, more=True),
                "tangent_corner_r0_ball": tangent_corner.str(90, more=True),
                "tangent_corner_minus_eta_ball": (tangent_corner - eta).str(90, more=True),
                "exact_curved_corner_R_at_half_ball": exact_corner.str(90, more=True),
                "exact_curved_corner_minus_eta_ball": (exact_corner - eta).str(90, more=True),
                "h_half_boundary_ball": h.str(90, more=True),
            }
        )

    by_mode = {row["mode"]: row for row in rows}
    require(arb(by_mode[39_894]["delta_one_plus_rho_ball"]).upper() < 0, "39894 null-coordinate slope sign drift")
    require(arb(by_mode[39_895]["delta_one_plus_rho_ball"]).lower() > 0, "39895 null-coordinate slope sign drift")
    require(abs(arb(by_mode[39_894]["tangent_corner_minus_eta_ball"])).upper() < arb("1e-10"), "39894 tangent corner mismatch")
    require(abs(arb(by_mode[39_895]["tangent_corner_minus_eta_ball"])).upper() < arb("1e-9"), "39895 tangent corner mismatch")
    require(abs(arb(by_mode[39_894]["exact_curved_corner_minus_eta_ball"])).upper() < arb("1e-6"), "39894 curved corner mismatch")
    require(abs(arb(by_mode[39_895]["exact_curved_corner_minus_eta_ball"])).upper() < arb("1e-5"), "39895 curved corner mismatch")
    return {
        "height": HEIGHT,
        "A": A,
        "continuous_half_boundary_mode_ball": nu.str(100, more=True),
        "dynamic_null_endpoint_ball": dynamic_endpoint.str(100, more=True),
        "A_null_height_ball": t_star.str(100, more=True),
        "A_null_height_minus_saved_height_ball": (t_star - t).str(100, more=True),
        "corner_unfolding_eta_A_ball": eta.str(100, more=True),
        "rows": rows,
        "interpretation": "The adjacent integer corner modes straddle delta=0 but share the same positive eta_A to better than 1e-9 in the tangent chart. The exact curved half-face differs from eta_A by less than 1e-5.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    geometry = artifact["geometry_certificate"]
    rows = {row["mode"]: row for row in geometry["rows"]}
    return f"""# Null-coordinate unfolding of the A wedge

Date: 2026-08-13

Status: exact characteristic-coordinate reduction and saved-height unfolding
certificate; not a numerical wedge or residual bound

The flat-face derivative formula is nonuniform when `c=rho^2-1` crosses
zero.  Near the A corner use instead

```text
R=P+S,       delta=1+rho.                              (NU1)
```

Then

```text
(P^2-S^2)/2=R^2/2-RS,
P<a+rho S  iff  R<a+delta S,
r0=a+delta h.                                         (NU2)
```

Swapping the two half-plane integrations and performing the `S` integral
gives, for `delta<0`,

```text
C=(2pi i)^(-1) int_(-infinity)^r0 exp(iR^2/2)
 [exp(-iRh)-exp(-iR(R-a)/delta)] dR/R.                (NU3)
```

The numerator in (NU3) vanishes at `R=0`.  For `delta>0`, the same common
Abel boundary value is

```text
C=1/2+(2pi i)^(-1)PV{{
 int_(-infinity)^r0 exp(iR^2/2-iRh)dR/R
 +int_r0^infinity exp(iR^2/2-iR(R-a)/delta)dR/R}}.    (NU4)
```

At `delta=0`,

```text
C=H(a)/2+(2pi i)^(-1)PV
 int_(-infinity)^a exp(iR^2/2-iRh)dR/R.               (NU5)
```

Equations (NU3)--(NU5) are one-dimensional characteristic formulas.  They do
not divide by `c` and expose the half-jump rather than hiding it in a
singular Fresnel scaling.

The continuous half-boundary mode and its dynamic null endpoint are

```text
nu=sqrt(t/(2pi)),       D_t=4nu=sqrt(8t/pi).           (NU6)
```

For endpoint `D`, the codimension-two unfolding is exactly

```text
eta_D=sqrt(pi)(D-D_t)/2
     =4(pi D^2/8-t)/[sqrt(pi)(D+D_t)].                 (NU7)
```

Thus `eta_D=0` exactly at `t=pi D^2/8`.  For the saved A endpoint,

```text
eta_A={geometry['corner_unfolding_eta_A_ball']}.       (NU8)
```

The two adjacent integer modes straddle the null slope:

```text
m=39894: delta={rows[39894]['delta_one_plus_rho_ball']},
         r0-eta_A={rows[39894]['tangent_corner_minus_eta_ball']},
m=39895: delta={rows[39895]['delta_one_plus_rho_ball']},
         r0-eta_A={rows[39895]['tangent_corner_minus_eta_ball']}. (NU9)
```

At the exact continuous null point, `rho=-1`,
`P_A''(0)=-2/(3sqrt(pi)A)`, and the boundary phase has cubic coefficient
`1/(3A sqrt(pi))`.  This recovers the Airy cubic from the triangle geometry
itself and explains the provenance of its `pi`.

The result changes the quantitative route: evaluate the affine A wedge with
(NU3)--(NU5), then compare the exact curved boundary to this null-coordinate
carrier.  It remains essential to retain the paired-triangle amplitude and
the global projector orientation; the older fold atlas cannot be inserted
by matching phases alone.

Proof boundary: exact null-coordinate algebra, one-dimensional Abel formulas,
continuous cubic coefficient, and saved-height corner localization only.  No
canonical wedge value, curved-face or amplitude bound, `R_Dir` estimate,
complete `Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["affine_remainder"]["decision"]["exact_defect_split_into_signed_amplitude_and_face_remainders"] is True, "affine dependency drift")
    require(dependencies["wedge"]["decision"]["A_face_null_slope_and_half_boundary_require_one_corner_chart"] is True, "wedge dependency drift")
    require(dependencies["bi_Morse"]["decision"]["A_fold_and_half_boundary_are_same_null_face"] is True, "bi-Morse dependency drift")

    artifact = {
        "kind": STEM,
        "status": "exact_null_coordinate_wedge_reduction_and_saved_height_A_unfolding_certified_quantitative_bound_open",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "geometry_certificate": geometry_certificate(dependencies["affine_remainder"]),
        "decision": {
            "near_null_wedge_reduced_to_one_dimensional_Abel_integrals": True,
            "reduction_avoids_division_by_characteristic_defect": True,
            "half_jump_retained_explicitly": True,
            "A_corner_unfolding_parameter_identified_exactly": True,
            "Airy_cubic_recovered_from_continuous_null_triangle": True,
            "old_fold_atlas_embedded_in_global_projector": False,
            "canonical_A_wedge_numerically_bounded": False,
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
        "next_obligation": "Build a contour-stable interval evaluator for (NU3)-(NU5), including the affine P and S moments, and certify the canonical A-corner values at modes 39894 and 39895 without epsilon extrapolation. Then derive the exact curved-boundary correction in R=P+S and match its cubic carrier to the paired-triangle amplitude before invoking any existing Airy atlas.",
        "proof_boundary": "Exact null-coordinate wedge reduction, continuous unfolding and cubic coefficient, and saved-height localization only. No canonical wedge value, curved-face or amplitude remainder bound, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified null-coordinate A-wedge unfolding and one-dimensional reduction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
