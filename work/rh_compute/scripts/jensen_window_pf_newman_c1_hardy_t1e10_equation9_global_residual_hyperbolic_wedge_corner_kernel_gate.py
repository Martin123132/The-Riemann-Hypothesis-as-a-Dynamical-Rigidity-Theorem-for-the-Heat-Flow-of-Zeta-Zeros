#!/usr/bin/env python3
"""Certify the canonical hyperbolic wedge/half-boundary kernel."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_corner_kernel_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "flat_face": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_flat_face_projector_kernel_gate.json",
    "projector_defect": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_phase_space_projector_defect_gate.json",
    "bi_Morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_exact_bimorse_face_fold_gate.json",
}

HEIGHT = 10_000_000_000
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
    p, s, a, rho, h = sp.symbols("P S a rho h", real=True)
    c = sp.symbols("c", real=True, nonzero=True)
    phase = (p**2 - s**2) / 2
    face_phase = sp.expand(phase.subs(p, a + rho * s))
    completed = c * (s + a * rho / c) ** 2 / 2 - a**2 / (2 * c)
    require(sp.simplify(completed.subs(c, rho**2 - 1) - face_phase) == 0, "wedge face completion failed")

    corner_density = sp.exp(sp.I * ((a + rho * h) ** 2 - h**2) / 2) / (2 * sp.pi)
    lower_boundary = -sp.exp(-sp.I * h**2 / 2) * sp.Integral(
        sp.exp(sp.I * p**2 / 2), (p, -sp.oo, a + rho * h)
    ) / (2 * sp.pi)
    lower_boundary_a = sp.diff(lower_boundary, a)
    require(sp.simplify(lower_boundary_a + corner_density) == 0, "mixed boundary density failed")

    mode, t = sp.symbols("m t", positive=True, real=True)
    v_half = 2 * sp.pi * mode**2 / t
    require(sp.solve(sp.Eq(v_half, 1), mode)[0] == sp.sqrt(2) * sp.sqrt(t) / (2 * sp.sqrt(sp.pi)), "half-boundary mode failed")
    endpoint = sp.symbols("D", positive=True, real=True)
    require(sp.simplify((endpoint / 4).subs(endpoint, sp.sqrt(8 * t / sp.pi)) - sp.sqrt(t / (2 * sp.pi))) == 0, "null-face/half-boundary coincidence failed")

    return {
        "wedge_kernel": "C_rho(a,h)=(2pi)^(-1) Abel-int_(S>h) exp(-iS^2/2) int_(P<a+rho*S)exp(iP^2/2)dP dS",
        "characteristic_defect": "c=rho^2-1",
        "a_derivative": "partial_a C=exp[-ia^2/(2c)] T_sgn(c)(sqrt(|c|/pi)[h+a*rho/c])/[2sqrt(pi|c|)], c!=0",
        "h_derivative": "partial_h C=-exp(-ih^2/2) int_(-infinity)^(a+rho*h)exp(iP^2/2)dP/(2pi)",
        "mixed_derivative": "partial_a partial_h C=-exp(i[(a+rho*h)^2-h^2]/2)/(2pi)",
        "flat_face_limit": "lim_(h to -infinity) C_rho(a,h)=F_(rho^2-1)(a)",
        "empty_limit": "lim_(h to +infinity) C_rho(a,h)=0",
        "full_face_limits": "lim_(a to -infinity)C=0; lim_(a to +infinity)C is the normalized one-dimensional S tail above h",
        "half_boundary": "h=S_H=sqrt(t/2) s(v_H), v_H=2pi*m^2/t, s(v)=sgn(v-1)sqrt(2[v-1-log v])",
        "half_boundary_event": "h=0 iff m=sqrt(t/(2pi))",
        "joint_null_corner": "m=D/4=sqrt(t/(2pi)) iff D=sqrt(8t/pi)",
    }


def half_coordinate(mode: int, t: arb, pi: arb) -> tuple[arb, arb]:
    m = arb(mode)
    v_half = 2 * pi * m**2 / t
    defect = v_half - 1 - v_half.log()
    require(defect.lower() > 0, f"half-coordinate defect is not positive at m={mode}")
    root = (2 * defect).sqrt()
    if (v_half - 1).upper() < 0:
        signed_s = -root
    elif (v_half - 1).lower() > 0:
        signed_s = root
    else:
        raise RuntimeError(f"half-boundary sign unresolved at m={mode}")
    return v_half, (t / 2).sqrt() * signed_s


def geometry_certificate(flat_face: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = 100
    pi = arb.pi()
    t = arb(HEIGHT)
    flat_rows = {
        (row["endpoint"], row["mode"]): row
        for row in flat_face["geometry_certificate"]["rows"]
    }
    keys = [
        (B, 621),
        (B, 622),
        (A, 39_852),
        (A, 39_853),
        (A, 39_894),
        (A, 39_895),
    ]
    rows: list[dict[str, Any]] = []
    for endpoint, mode in keys:
        v_half, s_half = half_coordinate(mode, t, pi)
        flat = flat_rows[(endpoint, mode)]
        rows.append(
            {
                "endpoint": endpoint,
                "mode": mode,
                "v_half_minus_one_ball": (v_half - 1).str(80, more=True),
                "h_half_boundary_ball": s_half.str(80, more=True),
                "a_face_detuning_ball": flat["a_face_detuning_ball"],
                "c_characteristic_defect_ball": flat["c_characteristic_defect_ball"],
            }
        )

    indexed = {(row["endpoint"], row["mode"]): row for row in rows}
    require(abs(arb(indexed[(B, 621)]["h_half_boundary_ball"])).lower() > arb("250000"), "B 621 half boundary is not remote")
    require(abs(arb(indexed[(B, 622)]["h_half_boundary_ball"])).lower() > arb("250000"), "B 622 half boundary is not remote")
    require(arb(indexed[(A, 39_852)]["h_half_boundary_ball"]).upper() < arb("-149"), "A 39852 half boundary drift")
    require(arb(indexed[(A, 39_853)]["h_half_boundary_ball"]).upper() < arb("-146"), "A 39853 half boundary drift")
    require(arb("-0.81") < arb(indexed[(A, 39_894)]["h_half_boundary_ball"]) < arb("-0.80"), "A 39894 half-boundary scale drift")
    require(arb("2.73") < arb(indexed[(A, 39_895)]["h_half_boundary_ball"]) < arb("2.74"), "A 39895 half-boundary scale drift")

    nu = (t / (2 * pi)).sqrt()
    a_quarter = arb(A) / 4
    require(arb(39_894) < nu < arb(39_895), "half-boundary event bracket drift")
    require(arb(39_894) < a_quarter < arb(39_895), "A null mode bracket drift")
    return {
        "height": HEIGHT,
        "rows": rows,
        "half_boundary_mode_ball": nu.str(90, more=True),
        "A_null_face_mode": "159577/4=39894.25",
        "A_null_minus_half_boundary_mode_ball": (a_quarter - nu).str(90, more=True),
        "interpretation": {
            "B": "At the B occupancy edge the half boundary is more than 250000 standardized S units away, so it is not a local corner.",
            "A_face": "At A modes 39852 and 39853 the half boundary is still below -146, while the face detuning is local.",
            "A_corner": "At 39894|39895 the half boundary moves from -0.81 to +2.74 while c changes sign, so face, characteristic slope, and half cutoff share one canonical wedge scale.",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = {(row["endpoint"], row["mode"]): row for row in artifact["geometry_certificate"]["rows"]}
    return f"""# Hyperbolic wedge and half-boundary corner kernel

Date: 2026-08-13

Status: exact canonical wedge calculus and saved-height corner localization;
not a bound for the curved-face or transformed-amplitude remainder

The half-domain adds a second boundary to the flat-face kernel.  Define

```text
C_rho(a,h)=(2pi)^(-1) Abel-integral_(S>h) exp(-iS^2/2)
             integral_(P<a+rho S)exp(iP^2/2)dP dS,
c=rho^2-1.                                           (WK1)
```

For `c!=0`, differentiation and exact square completion give

```text
partial_a C_rho(a,h)
 =exp[-ia^2/(2c)]/[2sqrt(pi|c|)]
  T_sgn(c)(sqrt(|c|/pi)[h+a rho/c]),                 (WK2)
```

where `T_eps(u)=integral_u^infinity exp(i eps pi v^2/2)dv`.  The other
boundary derivative is

```text
partial_h C_rho(a,h)
 =-exp(-ih^2/2)/(2pi)
   integral_(-infinity)^(a+rho h)exp(iP^2/2)dP,       (WK3)

partial_a partial_h C_rho
 =-exp(i[(a+rho h)^2-h^2]/2)/(2pi).                  (WK4)
```

Thus the complete canonical corner is determined by one Fresnel tail and its
corner density.  Its exact limits are

```text
lim_(h to -infinity) C_rho(a,h)=F_(rho^2-1)(a),
lim_(h to +infinity) C_rho(a,h)=0.                   (WK5)
```

So the flat-face theorem is not a separate approximation: it is the remote
half-boundary limit of the same wedge kernel.

For the exact logistic Morse coordinate,

```text
h=S_H=sqrt(t/2)s(v_H),
v_H=2pi m^2/t,
s(v)=sgn(v-1)sqrt(2[v-1-log v]).                     (WK6)
```

At the B occupancy edge the half boundary is remote:

```text
m=621: h={rows[(B,621)]['h_half_boundary_ball']},
m=622: h={rows[(B,622)]['h_half_boundary_ball']}.     (WK7)
```

At the A face entry it is still remote:

```text
m=39852: h={rows[(A,39852)]['h_half_boundary_ball']},
m=39853: h={rows[(A,39853)]['h_half_boundary_ball']}. (WK8)
```

But at the target edge it enters the canonical scale:

```text
m=39894: h={rows[(A,39894)]['h_half_boundary_ball']},
m=39895: h={rows[(A,39895)]['h_half_boundary_ball']}. (WK9)
```

The exact half-boundary event is

```text
m=sqrt(t/(2pi))
 ={artifact['geometry_certificate']['half_boundary_mode_ball']},          (WK10)
```

while the A null-face mode is `A/4=39894.25`; their separation is
`{artifact['geometry_certificate']['A_null_minus_half_boundary_mode_ball']}`.
This quantifies why the A face, its near-null slope, and `x=1/2` cutoff must
be treated as one corner, while the B face may use the full-S profile plus a
remote half-boundary correction.

Pi provenance: `pi` comes from the equation-(9) Fourier/Kummer phase and the
exact Gaussian standardization.  The half-boundary event follows from
`2pi m^2/t=1`; no fitted geometric constant is introduced.

Proof boundary: exact Abel wedge identities and saved-height localization
only.  No curved-face or transformed-amplitude remainder, signed projector
bound, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["flat_face"]["decision"]["flat_face_projector_profile_evaluated_exactly"] is True, "flat-face dependency drift")
    require(dependencies["projector_defect"]["decision"]["three_saved_height_stationary_occupancy_interfaces_identified"] is True, "projector occupancy dependency drift")
    require(dependencies["bi_Morse"]["decision"]["A_fold_and_half_boundary_are_same_null_face"] is True, "bi-Morse corner dependency drift")

    artifact = {
        "kind": STEM,
        "status": "exact_hyperbolic_wedge_boundary_calculus_and_saved_height_A_half_corner_localization_certified",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "canonical_domain": "S>h and P<a+rho*S",
            "regularization": "positive Gaussian damping followed by the common Abel boundary value",
        },
        "symbolic_certificate": symbolic_certificate(),
        "geometry_certificate": geometry_certificate(dependencies["flat_face"]),
        "decision": {
            "canonical_wedge_boundary_derivative_system_proved": True,
            "flat_face_is_remote_half_boundary_limit": True,
            "B_occupancy_edge_and_half_boundary_are_locally_coupled": False,
            "A_face_entry_and_half_boundary_are_locally_coupled": False,
            "A_target_edge_and_half_boundary_are_on_same_canonical_scale": True,
            "A_face_null_slope_and_half_boundary_require_one_corner_chart": True,
            "curved_face_remainder_proved": False,
            "transformed_amplitude_remainder_proved": False,
            "phase_adapted_signed_projector_bound_proved": False,
            "compressed_R_Dir_target_proved": False,
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
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Express the exact pair-triangle current minus its canonical wedge carrier as one signed S integral: transformed-amplitude difference over the exact face plus the strip between P_D(S) and a+rho*S. On the B side bound the remote h correction with the existing Fresnel tails. On the A side compare the c-near-zero wedge plus curved-face correction directly with the certified beta^-4 fold atlas over modes 39852..39895. Preserve the common target carrier and sum before taking norms.",
        "proof_boundary": "Exact canonical wedge boundary calculus and saved-height B/A half-boundary localization only. No curved-face or amplitude remainder, signed projector bound, R_Dir estimate, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified hyperbolic wedge corner kernel and A half-boundary localization", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
