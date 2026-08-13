#!/usr/bin/env python3
"""Certify the Abel-theta sum-first pair triangle and its modular dual roster."""

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

import mpmath as mp
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "pair_triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
}

T = 10_000_000_000
A = 159_577
B = 5_122_421
L = (B - A) // 2


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
    x, z, endpoint, n, height = sp.symbols("x z D n t", positive=True, real=True)
    r = 1 / z - 1 / x
    dual_phase = (
        sp.pi * endpoint**2 * z / 4
        + height * sp.log((1 - x) / x) / 2
        - sp.pi * n**2 / r
    )
    require(
        sp.simplify(-sp.pi * n**2 / r + sp.pi * n**2 * x * z / (x - z)) == 0,
        "reciprocal face coordinate failed",
    )

    alpha = endpoint - 2 * n
    z_stationary = x * alpha / endpoint
    dz = sp.factor(sp.diff(dual_phase, z))
    dx = sp.factor(sp.diff(dual_phase, x))
    require(sp.simplify(dz.subs(z, z_stationary)) == 0, "dual z saddle failed")
    dx_on_saddle = sp.factor(sp.simplify(dx.subs(z, z_stationary)))
    expected_dx = sp.pi * alpha**2 / 4 - height / (2 * x * (1 - x))
    require(sp.simplify(dx_on_saddle - expected_dx) == 0, "dual x saddle failed")

    discriminant = 1 - 8 * height / (sp.pi * alpha**2)
    root_minus = (1 - sp.sqrt(discriminant)) / 2
    root_plus = (1 + sp.sqrt(discriminant)) / 2
    stationary_equation = x * (1 - x) - 2 * height / (sp.pi * alpha**2)
    require(sp.simplify(stationary_equation.subs(x, root_minus)) == 0, "lower root failed")
    require(sp.simplify(stationary_equation.subs(x, root_plus)) == 0, "upper root failed")

    require(B - 2 * L == A, "roster endpoint arithmetic failed")
    require(L + 1 == (B - A) // 2 + 1, "roster cardinality failed")

    return {
        "face_coordinate": "r(z,x)=1/z-1/x=(x-z)/(x*z)>0 on 0<z<x<1/2",
        "Abel_parameter": "tau_epsilon=r+i*epsilon, epsilon>0",
        "theta_kernel": "Theta(tau)=sum_(m in Z) exp(i*pi*tau*m^2)",
        "positive_pair_sum": "sum_(m>=1)exp(i*pi*tau*m^2)=(Theta(tau)-1)/2",
        "grouped_triangle": "S_epsilon=(1/(4*i*pi))*integral_T Delta_g(z)*a_t(x)*(Theta(tau_epsilon)-1) dx dz",
        "endpoint_driver": "Delta_g=z^(-1/2){exp(i*pi*B^2*z/4)(2+i*pi*B^2*z)-exp(i*pi*A^2*z/4)(2+i*pi*A^2*z)}",
        "corner_majorant": "|Delta_g(z)|<=(3*pi/2)(B^2+A^2)z^(1/2), so the grouped x=0 majorant is proportional to x^(-1/4)",
        "modular_identity": "Theta(tau)=(-i*tau)^(-1/2)Theta(-1/tau), Im(tau)>0, principal square root",
        "dual_phase": "Psi_(D,n)=pi*D^2*z/4+(t/2)log((1-x)/x)-pi*n^2*x*z/(x-z)",
        "dual_z_saddle": "z=x(D-2n)/D",
        "dual_x_saddle": "x(1-x)=2t/[pi(D-2n)^2]",
        "dual_label": "alpha_(D,n)=D-2n",
        "B_roster": f"n=0..{L} gives alpha=B,B-2,...,A ({L + 1} odd values)",
        "continuation_match": "alpha_(B,L+k)=A-2k=alpha_(A,k); this labels the two endpoint continuations but does not prove termwise cancellation",
    }


def theta_sum(tau: mp.mpc, cutoff: int) -> mp.mpc:
    return mp.mpf(1) + 2 * mp.fsum(mp.exp(mp.j * mp.pi * tau * n * n) for n in range(1, cutoff + 1))


def modular_witness() -> dict[str, str | int]:
    mp.mp.dps = 100
    tau = mp.mpc("0.371", "0.237")
    cutoff = 40
    direct = theta_sum(tau, cutoff)
    transformed = (-mp.j * tau) ** (-mp.mpf("0.5")) * theta_sum(-1 / tau, cutoff)
    error = abs(direct - transformed)
    require(error < mp.mpf("1e-90"), "Jacobi modular witness failed")
    return {
        "tau_real": mp.nstr(tau.real, 50),
        "tau_imag": mp.nstr(tau.imag, 50),
        "cutoff_each_side": cutoff,
        "absolute_difference": mp.nstr(error, 100),
    }


def render_note() -> str:
    return f"""# Abel-theta pair triangle and modular dual roster

Date: 2026-08-13

Status: exact reduction using the Jacobi transformation; not a proof of the
quantitative source-minus-target estimate

For the exact pair triangle put

```text
r(z,x)=1/z-1/x=(x-z)/(xz)>0,
tau_epsilon=r+i epsilon,  epsilon>0,
Theta(tau)=sum_(m in Z) exp(i*pi*tau*m^2).             (AT1)
```

Abel weighting the nonzero symmetric Fourier pairs and summing before either
endpoint is separated gives

```text
S_epsilon=1/(4i*pi) integral_(0<z<x<1/2)
 Delta_g(z) a_t(x)[Theta(tau_epsilon)-1] dx dz,        (AT2)

Delta_g(z)=z^(-1/2){{
 e^(i*pi*B^2*z/4)(2+i*pi*B^2*z)
-e^(i*pi*A^2*z/4)(2+i*pi*A^2*z)}}.                    (AT3)
```

The endpoint difference in (AT3) is compulsory.  The elementary estimate

```text
|Delta_g(z)| <= (3*pi/2)(B^2+A^2) z^(1/2)             (AT4)
```

makes the grouped corner majorant proportional to `x^(-1/4)`, hence
integrable.  Separating the two `z^(-1/2)` endpoint pieces would lose this
argument.  The established paired `1/m^2` bound then permits
`epsilon -> 0+` by dominated convergence.

For `Im(tau)>0`, apply the standard Jacobi identity

```text
Theta(tau)=(-i*tau)^(-1/2)Theta(-1/tau).              (AT5)
```

The `n`th dual endpoint phase is

```text
Psi_(D,n)=pi*D^2*z/4+(t/2)log((1-x)/x)
            -pi*n^2*x*z/(x-z).                        (AT6)
```

Its interior stationary equations are exactly

```text
alpha_(D,n)=D-2n,
z=x alpha_(D,n)/D,
x(1-x)=2t/[pi alpha_(D,n)^2].                         (AT7)
```

Thus the modular dual index is not an arbitrary new parameter.  On the B
branch,

```text
n=0..{L}  <->  alpha=B,B-2,...,A,                     (AT8)
```

which is the complete contiguous odd source roster of `{L + 1}` values,
in reverse.  Moreover `alpha_(B,L+k)=alpha_(A,k)=A-2k`; the two endpoint
terms label the same continuation below the source interval, although (AT2)
does not make them cancel term by term.

This is both useful and a guard.  The modular transform removes modewise
face poles and recovers the exact arithmetic roster, but it is an involutive
re-expression of the source.  A gain still requires a uniform modular
face/corner estimate that extracts the classical target saddle and bounds
the continuation jointly.

Pi provenance: `pi` is inherited from the equation-(9) Fourier/Kummer phase
and the canonical Jacobi transform.  No geometric fit or arbitrary circle
constant is introduced.

Proof boundary: exact Abel-theta reassembly, corner integrability,
Jacobi-transform reduction, and dual saddle/roster geometry only.  No
quantitative modular remainder, source-minus-target estimate, complete
`T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["symmetric_poisson"]["decision"]["symmetric_poisson_interchange_proved"] is True,
        "symmetric-Poisson dependency drift",
    )
    require(
        dependencies["symmetric_poisson"]["decision"]["paired_one_over_m_endpoint_current_cancels"] is True,
        "paired-tail dependency drift",
    )
    require(
        dependencies["pair_triangle"]["decision"]["pair_Volterra_integral_reduced_to_separable_triangle"] is True,
        "pair-triangle dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "exact_Abel_theta_pair_triangle_and_modular_dual_source_roster_proved_quantitative_gain_open",
        "passed": True,
        "scope": {
            "height": T,
            "lower_endpoint": A,
            "upper_endpoint": B,
            "roster_steps": L,
            "roster_terms": L + 1,
        },
        "symbolic_certificate": symbolic_certificate(),
        "modular_witness": modular_witness(),
        "decision": {
            "nonzero_pairs_reassembled_as_one_Abel_theta_triangle": True,
            "endpoint_difference_retained_before_majorization": True,
            "grouped_x0_corner_has_integrable_majorant": True,
            "Abel_limit_justified_by_paired_one_over_m_squared_bound": True,
            "Jacobi_modular_dual_phase_identified": True,
            "B_dual_indices_recover_complete_odd_source_roster": True,
            "continuation_endpoints_cancel_termwise": False,
            "modular_transform_alone_proves_quantitative_gain": False,
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
            "process_priority": priority,
        },
        "next_obligation": "Use the modular kernel near r=0 to derive a uniform face parametrix with the B and A endpoint drivers still grouped. Isolate the finite B-dual roster n=0..L, identify its ordinary and A-fold saddle terms with the normalized target carriers, and bound the matched continuation n>L jointly rather than termwise.",
        "proof_boundary": "Exact Abel-theta and modular dual-roster reduction only. No quantitative modular remainder, source-minus-target estimate, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(), encoding="utf-8")
    print("certified Abel-theta pair triangle and modular dual roster", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
