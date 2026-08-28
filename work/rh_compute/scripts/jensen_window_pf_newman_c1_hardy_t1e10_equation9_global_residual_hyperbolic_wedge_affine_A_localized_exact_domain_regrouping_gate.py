#!/usr/bin/env python3
"""Certify the localized exact-domain regrouping and exterior saddle rosters."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_localized_exact_domain_regrouping_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "affine_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate.json",
    "orientation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate.json",
    "roster": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_transition_roster_gate.json",
    "local_face": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_curved_face_strip_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
LOWER_START = 39_853
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))
Y0_TEXT = "0.0037"
PRECISION = 100
EXPECTED_TANGENT_EXTERIOR = tuple(range(39_884, 39_895))
EXPECTED_EXACT_EXTERIOR = tuple(range(39_927, 39_937))


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
    exact, tangent, local = sp.symbols("chi_E chi_T chi_L")
    exterior = 1 - local
    full_face = exact - tangent
    local_face = local * full_face
    exterior_face = exterior * full_face
    require(sp.expand(tangent + local_face + exterior_face - exact) == 0, "full face regrouping failed")
    require(
        sp.expand(tangent + local_face - (local * exact + exterior * tangent)) == 0,
        "partial-localization contamination identity failed",
    )
    require(
        sp.expand(tangent + full_face - (local * exact + exterior * exact)) == 0,
        "localized exact-domain split failed",
    )

    x, m, endpoint, t = sp.symbols("x m D t", positive=True, real=True)
    endpoint_phase = sp.pi * m**2 / x + sp.pi * endpoint**2 * x / 4
    outer_phase = -sp.pi * m**2 / x + t * sp.log((1 - x) / x) / 2
    face_phase = sp.simplify(endpoint_phase + outer_phase)
    expected_face_phase = sp.pi * endpoint**2 * x / 4 + t * sp.log((1 - x) / x) / 2
    require(sp.simplify(face_phase - expected_face_phase) == 0, "common face phase failed")

    face_derivative = sp.factor(sp.diff(face_phase, x))
    expected_derivative = sp.pi * endpoint**2 / 4 - t / (2 * x * (1 - x))
    require(sp.simplify(face_derivative - expected_derivative) == 0, "face derivative failed")

    a, rho, S = sp.symbols("a rho S", real=True)
    tangent_phase = ((a + rho * S) ** 2 - S**2) / 2
    tangent_derivative = sp.diff(tangent_phase, S)
    require(
        sp.simplify(tangent_derivative - (rho * a + (rho**2 - 1) * S)) == 0,
        "tangent phase derivative failed",
    )

    return {
        "indicator_identity": "chi_T+(chi_E-chi_T)=chi_L chi_E+chi_X chi_E, chi_X=1-chi_L",
        "partial_localization_guard": "chi_T+chi_L(chi_E-chi_T)=chi_L chi_E+chi_X chi_T",
        "interpretation": "global tangent plus only the local face strip still contains the complete tangent exterior",
        "localized_exact_affine_domain": "C_aff,exact=C_aff,exact,loc+C_aff,exact,ext",
        "common_exact_face_phase": "Phi_face,A(x)=pi A^2 x/4+(t/2)log((1-x)/x), independent of m",
        "common_exact_face_derivative": "Phi_face,A'(x)=pi A^2/4-t/[2x(1-x)]",
        "tangent_boundary_phase": "F_T(S)=((a+rho S)^2-S^2)/2",
        "tangent_boundary_derivative": "F_T'(S)=rho a+(rho^2-1)S",
    }


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    y0 = arb(Y0_TEXT)
    S0 = (t * (y0 - (1 + y0).log())).sqrt()
    kappa = (1 - 8 * t / (pi * endpoint**2)).sqrt()
    x_star = (1 - kappa) / 2
    face_hessian = t * (1 - 2 * x_star) / (2 * x_star**2 * (1 - x_star) ** 2)

    rows: list[dict[str, Any]] = []
    tangent_exterior: list[int] = []
    exact_exterior: list[int] = []
    p0_values: list[arb] = []
    for mode_int in MODES:
        mode = arb(mode_int)
        r = t / (2 * pi * mode**2)
        denominator = 2 * pi * mode**2 + t
        alpha = 2 * mode + t / (pi * mode)
        a = (pi * mode / alpha).sqrt() * (endpoint - alpha)
        rho = -(2 * t).sqrt() * (pi * endpoint * mode + denominator) / (
            2 * denominator ** arb("1.5")
        )
        c = rho**2 - 1
        d = rho * a
        tangent_root: arb | None = None
        has_tangent_exterior = False
        if c.lower() > 0:
            tangent_root = -d / c
            has_tangent_exterior = tangent_root.lower() > S0.upper()
        if has_tangent_exterior:
            tangent_exterior.append(mode_int)

        x0 = 1 / (1 + r * (1 + y0))
        exact_stationary_y = (1 / x_star - 1) / r - 1
        has_exact_exterior = exact_stationary_y.lower() > y0.upper()
        if has_exact_exterior:
            exact_exterior.append(mode_int)

        g0 = 1 + r * (1 + y0)
        p0 = (pi / 2).sqrt() * (endpoint / g0.sqrt() - 2 * mode * g0.sqrt())
        p0_values.append(p0)
        rows.append(
            {
                "mode": mode_int,
                "a_ball": a.str(65, more=True),
                "rho_ball": rho.str(65, more=True),
                "c_ball": c.str(65, more=True),
                "rho_a_ball": d.str(65, more=True),
                "tangent_stationary_S_ball": tangent_root.str(65, more=True) if tangent_root is not None else None,
                "tangent_stationary_in_exterior": has_tangent_exterior,
                "x_cutoff_ball": x0.str(65, more=True),
                "exact_face_stationary_y_ball": exact_stationary_y.str(65, more=True),
                "exact_face_stationary_in_exterior": has_exact_exterior,
                "P_exact_at_y0_ball": p0.str(65, more=True),
            }
        )

    require(tuple(tangent_exterior) == EXPECTED_TANGENT_EXTERIOR, "tangent exterior saddle roster drift")
    require(tuple(exact_exterior) == EXPECTED_EXACT_EXTERIOR, "exact exterior saddle roster drift")
    require(max(value.upper() for value in p0_values) < arb("-261.2"), "exterior Fresnel-tail margin drift")
    require(face_hessian.lower() > arb("8e7"), "exact face stationary Hessian margin drift")

    by_mode = {row["mode"]: row for row in rows}
    require(arb(by_mode[39_883]["tangent_stationary_S_ball"]).upper() < S0.lower(), "39883 tangent cutoff guard drift")
    require(arb(by_mode[39_884]["tangent_stationary_S_ball"]).lower() > S0.upper(), "39884 tangent cutoff guard drift")
    require(arb(by_mode[39_894]["c_ball"]).lower() > 0, "39894 tangent c sign drift")
    require(arb(by_mode[39_895]["c_ball"]).upper() < 0, "39895 tangent c sign drift")
    require(arb(by_mode[39_926]["exact_face_stationary_y_ball"]).upper() < y0.lower(), "39926 exact cutoff guard drift")
    require(arb(by_mode[39_927]["exact_face_stationary_y_ball"]).lower() > y0.upper(), "39927 exact cutoff guard drift")

    return {
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": len(rows),
        "y_cutoff": Y0_TEXT,
        "S_cutoff_ball": S0.str(75, more=True),
        "common_exact_face_kappa_ball": kappa.str(75, more=True),
        "common_exact_face_stationary_x_ball": x_star.str(75, more=True),
        "common_exact_face_stationary_hessian_ball": face_hessian.str(75, more=True),
        "tangent_exterior_stationary_modes": tangent_exterior,
        "exact_face_exterior_stationary_modes": exact_exterior,
        "P_exact_at_cutoff_uniform_upper": "-261.2",
        "rows": rows,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    tangent = c["tangent_exterior_stationary_modes"]
    exact = c["exact_face_exterior_stationary_modes"]
    rows = {row["mode"]: row for row in c["rows"]}
    return f"""# Localized exact affine A-domain and exterior saddle rosters

Date: 2026-08-13

Status: exact regrouping and rigorous exterior-stationarity guard; not an
exterior integral bound or complete A carrier

Let `chi_E` and `chi_T` denote the exact and tangent A domains, and split
`S>h` at `S_0` into `chi_L+chi_X=1`.  Before norms,

```text
chi_T+(chi_E-chi_T)=chi_E=chi_L chi_E+chi_X chi_E.    (LR1)
```

However, retaining only the already certified local face strip gives

```text
chi_T+chi_L(chi_E-chi_T)=chi_L chi_E+chi_X chi_T.     (LR2)
```

Thus global tangent plus local face is not the local exact carrier: it still
contains the entire tangent exterior.  The cancellation-safe next object is

```text
C_aff,exact=C_aff,exact,loc+C_aff,exact,ext,           (LR3)
```

with no tangent exterior in the second term.

This matters because the tangent boundary phase

```text
F_T(S)=((a+rho S)^2-S^2)/2,
F_T'(S)=rho a+(rho^2-1)S                              (LR4)
```

has an artificial exterior stationary point for exactly

```text
{tangent[0]}..{tangent[-1]} ({len(tangent)} modes).   (LR5)
```

The adjacent guards are

```text
mode 39883: S_T*={rows[39883]['tangent_stationary_S_ball']} < S_0,
mode 39884: S_T*={rows[39884]['tangent_stationary_S_ball']} > S_0,
mode 39894: rho^2-1={rows[39894]['c_ball']} >0,
mode 39895: rho^2-1={rows[39895]['c_ball']} <0.        (LR6)
```

These saddles belong to the unbounded linearization, not the exact triangle.
They must be removed algebraically by (LR3), not bounded as though the
exterior were nonstationary.

After the common carrier is restored, the exact boundary phase on `z=x` is

```text
Phi_face,A(x)=pi A^2 x/4+(t/2)log((1-x)/x),           (LR7)
```

which is independent of `m`.  Its unique half-domain stationary point is

```text
kappa=sqrt(1-8t/(pi A^2))={c['common_exact_face_kappa_ball']},
x_*=(1-kappa)/2={c['common_exact_face_stationary_x_ball']}. (LR8)
```

At the common cutoff `y_0={Y0_TEXT}`, this true stationary point remains in
the exterior for exactly

```text
{exact[0]}..{exact[-1]} ({len(exact)} modes).          (LR9)
```

The adjacent exact guards are

```text
mode 39926: y_*={rows[39926]['exact_face_stationary_y_ball']} < y_0,
mode 39927: y_*={rows[39927]['exact_face_stationary_y_ball']} > y_0. (LR10)
```

Finally, every exact exterior inner endpoint satisfies

```text
P_A(y_0)<-261.2,                                      (LR11)
```

so the inner Fresnel tail can be expanded with an explicit inverse-`P`
remainder.  Equations (LR7)--(LR10) show that its leading boundary currents
must be summed in the common `x` phase and supplied with a ten-mode incomplete
stationary transition; a blanket exterior integration-by-parts estimate is
not valid.

Pi provenance: every `pi` in (LR4)--(LR11) comes from the exact equation-(9)
bi-Morse phase, the triangle face `z=x`, and the saved endpoint normalization.
No fitted constant is introduced.

Proof boundary: exact localization algebra, the complete 84-mode tangent and
exact exterior stationary rosters, common exact-face phase, and cutoff
Fresnel-tail margin only.  No local exact-affine value, exterior Fresnel
current estimate, transformed-amplitude remainder, complete A endpoint
theorem, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["affine_remainder"]["decision"]["exact_defect_split_into_signed_amplitude_and_face_remainders"] is True,
        "affine-remainder dependency drift",
    )
    require(
        dependencies["orientation"]["decision"]["common_phase_and_paired_endpoint_orientation_retained"]
        if "common_phase_and_paired_endpoint_orientation_retained" in dependencies["orientation"]["decision"]
        else dependencies["orientation"]["decision"]["odd_endpoint_phase_reduced_to_common_classical_carrier"],
        "orientation dependency drift",
    )
    require(dependencies["roster"]["decision"]["complete_84_mode_A_transition_roster_evaluated"] is True, "roster dependency drift")
    require(dependencies["local_face"]["decision"]["signed_local_curved_face_strip_summed_before_norms"] is True, "local-face dependency drift")

    artifact = {
        "kind": STEM,
        "status": "localized_exact_affine_A_domain_regrouping_and_exterior_stationary_rosters_certified_quantitative_exterior_open",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "global_tangent_plus_full_face_strip_equals_exact_affine_domain": True,
            "global_tangent_plus_local_face_retains_tangent_exterior": True,
            "tangent_exterior_must_be_removed_before_estimation": True,
            "eleven_artificial_tangent_exterior_stationary_modes_identified": True,
            "exact_face_phase_is_common_after_carrier_restoration": True,
            "ten_true_exact_face_stationary_transitions_remain_in_exterior": True,
            "uniform_exterior_inner_Fresnel_endpoint_below_minus_261p2": True,
            "blanket_exterior_nonstationary_IBP_is_valid": False,
            "exterior_exact_affine_current_bound_proved": False,
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
            "sympy_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Evaluate the localized exact affine domain, not global tangent plus a truncated face correction. Expand the exact exterior inner Fresnel tail through inverse P^5 with an explicit inverse P^7 remainder, sum the leading exact-face boundary currents in the common x phase, and treat modes 39927..39936 with an incomplete stationary transition. Then certify the exact-minus-affine amplitude under the same localization.",
        "proof_boundary": "Exact localization algebra, exterior stationary rosters, common exact-face phase, and cutoff Fresnel-tail margin only. No local exact-affine value, exterior current estimate, transformed-amplitude remainder, complete A endpoint theorem, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified localized exact A domain and exterior saddle rosters", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
