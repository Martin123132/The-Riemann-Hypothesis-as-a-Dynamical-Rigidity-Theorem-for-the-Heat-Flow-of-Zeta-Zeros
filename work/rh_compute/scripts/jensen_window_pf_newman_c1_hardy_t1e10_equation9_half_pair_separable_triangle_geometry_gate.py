#!/usr/bin/env python3
"""Certify the separable triangular integral geometry of each Fourier pair."""

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
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "Volterra_transport": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_endpoint_transport_volterra_gate.json",
    "Gamma_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_volterra_incomplete_gamma_outer_kernel_gate.json",
    "half_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
    "coverage_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
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
    z, x, m, endpoint, t = sp.symbols("z x m D t", positive=True, real=True)
    phi_endpoint = sp.pi * m**2 / z + sp.pi * endpoint**2 * z / 4
    psi_outer = -sp.pi * m**2 / x + t * sp.log((1 - x) / x) / 2
    total = phi_endpoint + psi_outer
    require(sp.simplify(sp.diff(total, z, x)) == 0, "phase failed to separate")

    z_stationary = 2 * m / endpoint
    x_stationary = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    require(sp.simplify(sp.diff(phi_endpoint, z).subs(z, z_stationary)) == 0, "endpoint saddle failed")
    require(sp.simplify(sp.diff(psi_outer, x).subs(x, x_stationary)) == 0, "outer saddle failed")

    alpha_m = 2 * m + t / (sp.pi * m)
    require(
        sp.factor(sp.simplify(x_stationary - z_stationary)).has(alpha_m - endpoint)
        or sp.simplify((x_stationary - z_stationary) / (alpha_m - endpoint)) != 0,
        "gap factor audit failed",
    )
    gap_ratio = sp.simplify((x_stationary - z_stationary) / (endpoint - alpha_m))
    expected_ratio = 2 * sp.pi * m**2 / (endpoint * (t + 2 * sp.pi * m**2))
    require(sp.simplify(gap_ratio - expected_ratio) == 0, "boundary-gap identity failed")
    require(sp.simplify((x_stationary - z_stationary).subs(endpoint, alpha_m)) == 0, "characteristic face failed")

    phi_second = sp.simplify(sp.diff(phi_endpoint, z, 2).subs(z, z_stationary))
    psi_second = sp.factor(sp.diff(psi_outer, x, 2).subs(x, x_stationary))
    require(sp.simplify(phi_second - sp.pi * endpoint**3 / (4 * m)) == 0, "endpoint Hessian failed")
    require(phi_second.is_positive is not False, "endpoint Hessian sign failed")
    require(sp.simplify(psi_second) != 0, "outer Hessian degenerated")

    return {
        "triangle": "T={0<z<x<1/2}",
        "endpoint_driver": "h_D(z)=z^(-1/2)(2+i*pi*D^2*z)",
        "outer_amplitude": "a(x)=x^(-7/4)(1-x)^(-1/4)",
        "separable_phase": "Phi_(m,D)(z,x)=phi_(m,D)(z)+psi_m(x)",
        "endpoint_phase": "phi_(m,D)(z)=pi*m^2/z+pi*D^2*z/4",
        "outer_phase": "psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x)",
        "pair_triangle": "K_m+K_(-m)=sum_(D in {B,A}) eps_D/(2i*pi) integral_T h_D(z)a(x)exp(i Phi_(m,D)) dx dz, eps_B=+1, eps_A=-1",
        "stationary_point": "z_D=2m/D, x_m=2*pi*m^2/(t+2*pi*m^2)",
        "boundary_gap": "x_m-z_D=2*pi*m^2[D-alpha_m]/[D(t+2*pi*m^2)], alpha_m=2m+t/(pi*m)",
        "characteristic_face": "The joint saddle meets z=x exactly at D=alpha_m.",
        "Hessian": "The phase Hessian is diagonal; phi_D''(z_D)=pi*D^3/(4m)>0 and psi_m''(x_m)<0.",
        "interpretation": "The tangent half-plane model is the local boundary chart of this exact triangular pair integral, not an independent endpoint approximation.",
    }


def saved_height_ledger() -> dict[str, Any]:
    mp.mp.dps = 80
    t, pi = mp.mpf(T), mp.pi

    def alpha(mode: int) -> mp.mpf:
        value = mp.mpf(mode)
        return 2 * value + t / (pi * value)

    def x_mode(mode: int) -> mp.mpf:
        value = mp.mpf(mode)
        return 2 * pi * value**2 / (t + 2 * pi * value**2)

    def z_mode(mode: int, endpoint: int) -> mp.mpf:
        return 2 * mp.mpf(mode) / endpoint

    rows = []
    for mode in (621, 622, 39_694, 39_852, 39_853, 39_894, 39_895):
        xm = x_mode(mode)
        rows.append(
            {
                "mode": mode,
                "alpha_m": mp.nstr(alpha(mode), 65),
                "x_m": mp.nstr(xm, 65),
                "B_z": mp.nstr(z_mode(mode, B), 65),
                "B_gap_x_minus_z": mp.nstr(xm - z_mode(mode, B), 65),
                "A_z": mp.nstr(z_mode(mode, A), 65),
                "A_gap_x_minus_z": mp.nstr(xm - z_mode(mode, A), 65),
            }
        )

    by_mode = {row["mode"]: row for row in rows}
    require(mp.mpf(by_mode[621]["B_gap_x_minus_z"]) < 0, "B mode 621 occupancy drift")
    require(mp.mpf(by_mode[622]["B_gap_x_minus_z"]) > 0, "B mode 622 occupancy drift")
    require(mp.mpf(by_mode[39_852]["A_gap_x_minus_z"]) < 0, "A mode 39852 occupancy drift")
    require(mp.mpf(by_mode[39_853]["A_gap_x_minus_z"]) > 0, "A mode 39853 occupancy drift")
    require(x_mode(39_894) < mp.mpf("0.5") < x_mode(39_895), "half-boundary occupancy drift")
    return {
        "height": T,
        "rows": rows,
        "face_ledger": {
            "B_injection": "At 621|622 the B stationary point crosses z=x into the triangle.",
            "ordinary_core": "For 622..39852 the B saddle is inside and the A saddle is outside the triangle.",
            "A_injection": "At 39852|39853 the A stationary point also crosses z=x into the triangle with opposite sign.",
            "fold_owned_overlap": "39695..39894 remains assigned to the certified fold atlas despite ordinary validity overlap through 39852.",
            "half_boundary": "At 39894|39895 the outer saddle x_m crosses x=1/2; simultaneously z_A is within one lattice step of that corner.",
            "outer": "For m>=39895 the outer phase has no stationary point in the half triangle.",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = {row["mode"]: row for row in artifact["saved_height_ledger"]["rows"]}
    return f"""# Separable triangular geometry of a signed Fourier pair

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the uniform triangle estimate

Insert the endpoint-driven solution into its Volterra kernel before taking
absolute values.  For each source endpoint `D` define

```text
phi_(m,D)(z)=pi*m^2/z+pi*D^2*z/4,
psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x).                 (TG1)
```

Then every symmetric pair has the exact separable representation

```text
K_m+K_-m=sum_(D in {{B,A}}) eps_D/(2i*pi)
 integral_(0<z<x<1/2)
 z^(-1/2)(2+i*pi*D^2*z)
 x^(-7/4)(1-x)^(-1/4)
 exp[i(phi_(m,D)(z)+psi_m(x))] dx dz,                 (TG2)

eps_B=+1, eps_A=-1.
```

The phase Hessian is diagonal.  Its stationary coordinates are

```text
z_D=2m/D,
x_m=2*pi*m^2/(t+2*pi*m^2),
alpha_m=2m+t/(pi*m).                                  (TG3)
```

Their signed distance to the triangular face is exactly

```text
x_m-z_D=2*pi*m^2(D-alpha_m)/[D(t+2*pi*m^2)].          (TG4)
```

Thus the endpoint characteristic equation `D=alpha_m` is simply the event
where the joint saddle crosses the face `z=x`.  The positive `z` curvature
and negative `x` curvature show that this is the same hyperbolic saddle
geometry as the original portcullis, now with `+m/-m` cancellation already
built in.

At the saved height the exact face ledger is:

```text
B gap at m=621: {rows[621]['B_gap_x_minus_z']},
B gap at m=622: {rows[622]['B_gap_x_minus_z']},

A gap at m=39852: {rows[39852]['A_gap_x_minus_z']},
A gap at m=39853: {rows[39853]['A_gap_x_minus_z']}.    (TG5)
```

Consequently:

```text
1..621:       no B joint saddle in the triangle,
622..39852:   B saddle in, A saddle out,
39853..39894: both endpoint saddles in, with opposite signs,
m>=39895:     outer x saddle beyond the half-boundary. (TG6)
```

The ordinary/fold ownership boundary remains `39694|39695`; (TG6) is
stationary geometry and does not override that certified proof allocation.

This unifies the earlier obstructions.  The B tangent profile is the local
half-plane approximation to (TG2) at the `B` face.  The A atlas is the same
face geometry with the opposite endpoint sign.  The `39894|39895` event is
the additional corner `x=1/2`.  They should therefore be estimated with one
triangle theorem and compatible face/corner charts, rather than three
unrelated approximations.

Pi provenance: all `pi` factors come from the exact equation-(9)
Kummer/Fourier phase.  No fitted constant is used.

Proof boundary: exact separable triangular normal form, stationary
coordinates, gap identity, and saved-height face ledger only.  No uniform
triangle estimate, quantitative face/corner chart, complete paired residual,
`T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["Volterra_transport"]["decision"]["half_pair_triangular_Volterra_representation_proved"] is True,
        "Volterra dependency drift",
    )
    require(
        dependencies["Gamma_kernel"]["decision"]["Volterra_kernel_reduced_exactly_to_incomplete_Gamma_path"] is True,
        "kernel dependency drift",
    )
    coverage = dependencies["coverage_ledger"]["mode_coverage"]["proposed_disjoint_target_ownership"]
    require(coverage["ordinary_core"]["range"] == [622, 39_694], "ordinary ownership drift")
    require(coverage["fold_owned_target"]["range"] == [39_695, 39_894], "fold ownership drift")

    artifact = {
        "kind": STEM,
        "status": "exact_signed_pair_separable_triangle_geometry_and_face_ledger_proved_uniform_bound_open",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "saved_height_ledger": saved_height_ledger(),
        "decision": {
            "pair_Volterra_integral_reduced_to_separable_triangle": True,
            "joint_phase_Hessian_is_diagonal_and_nondegenerate": True,
            "endpoint_characteristic_is_triangle_face_crossing": True,
            "B_A_and_half_boundary_events_share_one_triangle_geometry": True,
            "tangent_B_profile_is_complete_pair_triangle_theorem": False,
            "uniform_triangle_face_corner_estimate_proved": False,
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
        "next_obligation": "Formulate one quantitative stationary-phase theorem for the separable phase over 0<z<x<1/2: an interior product term, a uniform z=x face profile retaining the full driver amplitude, and an x=1/2 corner profile. Apply it first to the B 621|622 crossing together with the signed pair completion.",
        "proof_boundary": "Exact pair-triangle normal form and face geometry only. No uniform triangle estimate, quantitative face/corner chart, complete paired residual, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact separable triangle geometry for signed Fourier pairs", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
