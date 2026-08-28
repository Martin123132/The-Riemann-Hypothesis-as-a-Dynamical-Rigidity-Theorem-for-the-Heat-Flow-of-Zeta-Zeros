#!/usr/bin/env python3
"""Reduce the joined finite source to an entire auxiliary-kernel difference."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_riemann_auxiliary_entire_kernel_finite_difference_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "pole_safe_join": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_PW_unowned_cell_pole_safe_two_boundary_contour_join_gate.json",
    "unequal_truncation_map": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_unequal_truncation_saddle_map_gate.json",
}

HEIGHT = 10_000_000_000
N_MINUS = 79_788
N_PLUS = 2_561_211
LABEL_COUNT = 2_481_423


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def kappa(z: mp.mpc) -> mp.mpc:
    return (mp.exp(-1j * mp.pi * z) - mp.exp(-1j * mp.pi * z * z)) / (
        2j * mp.sin(mp.pi * z)
    )


def kappa_auxiliary_form(z: mp.mpc) -> mp.mpc:
    return (
        mp.exp(-0.5j * mp.pi * (z * z + z))
        * mp.sin(0.5 * mp.pi * (z * z - z))
        / mp.sin(mp.pi * z)
    )


def prefix_kernel_direct(z: mp.mpc, terms: int) -> mp.mpc:
    return mp.exp(-1j * mp.pi * z * z) * mp.fsum(
        mp.exp(1j * mp.pi * (2 * r + 1) * z) for r in range(terms)
    )


def prefix_kernel_entire(z: mp.mpc, terms: int) -> mp.mpc:
    return kappa(z) - kappa(z - terms)


def window_kernel_direct(z: mp.mpc, n_minus: int, n_plus: int) -> mp.mpc:
    return mp.exp(-1j * mp.pi * z * z) * mp.fsum(
        mp.exp(1j * mp.pi * (2 * r + 1) * z) for r in range(n_minus, n_plus)
    )


def window_kernel_entire(z: mp.mpc, n_minus: int, n_plus: int) -> mp.mpc:
    return kappa(z - n_minus) - kappa(z - n_plus)


def fmt(value: mp.mpf | mp.mpc, digits: int = 35) -> str:
    return mp.nstr(value, digits)


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    require(CHECKER.is_file(), "missing independent checker")
    dependency_rows: dict[str, dict[str, str]] = {}
    for name, path in DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {name}")
        dependency_rows[name] = {"path": relative(path), "sha256": file_hash(path)}

    mp.mp.dps = 90
    point_rows: list[dict[str, Any]] = []
    points = (
        mp.mpc("0.23", "0.31"),
        mp.mpc("1.37", "-0.19"),
        mp.mpc("2.5", "0.27"),
        mp.mpc("-0.41", "0.22"),
    )
    for z in points:
        auxiliary_discrepancy = abs(kappa(z) - kappa_auxiliary_form(z))
        prefix_rows = []
        for terms in (1, 2, 5, 9):
            discrepancy = abs(prefix_kernel_direct(z, terms) - prefix_kernel_entire(z, terms))
            require(discrepancy < mp.mpf("1e-75"), "finite-prefix entire-kernel identity failed")
            prefix_rows.append({"terms": terms, "discrepancy_absolute": fmt(discrepancy)})
        window_discrepancy = abs(window_kernel_direct(z, 2, 7) - window_kernel_entire(z, 2, 7))
        require(auxiliary_discrepancy < mp.mpf("1e-75"), "auxiliary kernel form failed")
        require(window_discrepancy < mp.mpf("1e-74"), "window finite-difference identity failed")
        point_rows.append(
            {
                "z": fmt(z),
                "auxiliary_form_discrepancy_absolute": fmt(auxiliary_discrepancy),
                "prefix_rows": prefix_rows,
                "window_discrepancy_absolute": fmt(window_discrepancy),
            }
        )

    integer_rows = []
    for integer in (-7, -2, 0, 1, 4, 11):
        expected = mp.mpf(integer) - mp.mpf("0.5")
        eps = mp.mpf("1e-25")
        symmetric = (kappa(mp.mpf(integer) + eps) + kappa(mp.mpf(integer) - eps)) / 2
        discrepancy = abs(symmetric - expected)
        require(discrepancy < mp.mpf("1e-44"), "integer removable value failed")
        integer_rows.append(
            {
                "integer": integer,
                "removable_value": fmt(expected),
                "symmetric_probe_discrepancy_absolute": fmt(discrepancy),
            }
        )

    for integer in (-3, 0, 5):
        for terms in (1, 4, 13):
            exact_difference = (mp.mpf(integer) - mp.mpf("0.5")) - (
                mp.mpf(integer - terms) - mp.mpf("0.5")
            )
            require(exact_difference == terms, "integer prefix extension failed")

    actual_zero_value = (
        mp.mpf(-N_MINUS) - mp.mpf("0.5")
    ) - (mp.mpf(-N_PLUS) - mp.mpf("0.5"))
    require(actual_zero_value == LABEL_COUNT, "actual zero value drift")

    artifact: dict[str, Any] = {
        "kind": STEM,
        "status": "riemann_auxiliary_entire_kernel_finite_difference_certified",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "s": "1/2+i*t",
            "n_minus": N_MINUS,
            "n_plus": N_PLUS,
            "label_count": LABEL_COUNT,
        },
        "dependencies": dependency_rows,
        "entire_kernel": {
            "definition": "kappa(z)=(exp(-i*pi*z)-exp(-i*pi*z^2))/(2*i*sin(pi*z))",
            "auxiliary_form": "kappa(z)=exp(-i*pi*(z^2+z)/2)*sin(pi*(z^2-z)/2)/sin(pi*z)",
            "integer_removable_value": "kappa(k)=k-1/2 for every integer k",
            "is_entire": True,
            "primary_source": "https://arxiv.org/abs/2407.02016",
            "point_rows": point_rows,
            "integer_rows": integer_rows,
        },
        "finite_difference_identities": {
            "regularized_prefix_kernel": "exp(-i*pi*z^2)*sum_(r=0)^(N-1)exp(i*pi*(2r+1)z)=kappa(z)-kappa(z-N)",
            "regularized_prefix_packet": "R_N(z)=z^(-s)*(kappa(z)-kappa(z-N))",
            "actual_window_kernel": "exp(-i*pi*z^2)*D_W(z)=kappa(z-n_minus)-kappa(z-n_plus)",
            "actual_integrand": "g_W(z)=z^(-s)*(kappa(z-n_minus)-kappa(z-n_plus))",
            "actual_zero_removable_value": str(LABEL_COUNT),
            "all_csc_integer_poles_removed_before_integration": True,
        },
        "common_functional": {
            "definition": "L_s[h]=integral_0^infinity[Fresnel] z^(-s)h(z)dz-integral_C z^(-s)h(z)dz",
            "finite_prefix_source": "S_N=L_s[kappa(z)-kappa(z-N)]",
            "actual_source": "S_W=L_s[kappa(z-n_minus)-kappa(z-n_plus)]",
            "fully_joined_target": "L_s[kappa(z-n_minus)-kappa(z-n_plus)]-C_G*sum_(m=622)^39936m^(-s)-mathcal_A_A^nat",
            "common_base_kappa_cancels_before_norms": True,
        },
        "decision": {
            "entire_kernel_route_selected": True,
            "reason": "it is an exact global finite difference with no sine poles and no imported asymptotic constant",
            "next_action": "derive an endpoint-complete contour deformation or Mellin finite-difference formula for the single common functional, retaining Gamma and A subtraction inside the transform",
            "joined_packet_enclosed": False,
            "J_Z_enclosed": False,
            "D_K_enclosed": False,
            "rh_implication": False,
        },
        "runtime": {
            "seconds": time.perf_counter() - started,
            "active_compute_workers": 1,
            "priority": priority,
            "thread_caps": 1,
        },
        "source_hashes": {"builder": file_hash(BUILDER), "checker": file_hash(CHECKER)},
        "proof_boundary": (
            "Exact entire auxiliary-kernel identity, integer removable values, finite-prefix/window finite differences, "
            "and common-functional reassembly only. No quantitative contour deformation, Mellin remainder bound, "
            "joined-packet, J_Z, or D_K enclosure, and no non-A, all-height, Lambda<=0, PF-infinity, RH, or "
            "prize-level conclusion."
        ),
    }

    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(
        f"""# Entire Riemann-auxiliary kernel for the joined equation-(4) packet

Date: 2026-08-27

Status: exact finite-difference reduction; quantitative enclosure remains open.

Define

```text
kappa(z)=[exp(-i*pi*z)-exp(-i*pi*z^2)]/[2i sin(pi*z)]
        =exp[-i*pi(z^2+z)/2]
          sin[pi(z^2-z)/2]/sin(pi*z).
```

The numerator vanishes at every integer and l'Hopital gives

```text
kappa(k)=k-1/2,       k integer.
```

Integer translation leaves the linear quotient fixed and sends the Gaussian quotient to its `N`th csc packet.  Therefore

```text
R_N(z)=z^(-s)[kappa(z)-kappa(z-N)],

g_W(z)=z^(-s)[kappa(z-79788)-kappa(z-2561211)].
```

At `z=0`, the bracket has removable value `2561211-79788={LABEL_COUNT}`, exactly reproducing the locally integrable finite source.  No sine pole remains.

For

```text
L_s[h]=integral_0^infinity[Fresnel] z^(-s)h(z)dz
       -integral_C z^(-s)h(z)dz,
```

the complete source is now

```text
S_W=L_s[kappa(z-79788)-kappa(z-2561211)].
```

Thus the fully reassembled target is one entire-kernel finite difference minus the already-owned Gamma and A carriers.  This is the selected exact route for the next endpoint-complete deformation.

The kernel is the trigonometric kernel in Proposition 9 of Arias de Reyna's integral representation of Riemann's auxiliary function: https://arxiv.org/abs/2407.02016.  Only the kernel identity is imported; no auxiliary-function value or bound is transferred to `J_Z`.

Pi provenance: every `pi` is inherited from the Riemann-Siegel Gaussian/sine kernel.  The translation identity uses only integer parity.

## Proof boundary

{artifact['proof_boundary']}
""",
        encoding="utf-8",
    )
    print("certified entire auxiliary-kernel finite difference for the joined packet")
    print(f"result: {relative(RESULT)}")
    print(f"note: {relative(NOTE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
