#!/usr/bin/env python3
"""Derive the finite geometric RSI contour residual and its released poles."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_geometric_contour_residual_and_pole_release_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "A21_provenance": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_A21_exactness_provenance_and_contour_jump_target_gate.json",
    "normalization_repair": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_domain_Kummer_normalization_repair_gate.json",
}

SOURCE_ALPHA_MIN = 159_577
SOURCE_ALPHA_MAX = 5_122_421
SOURCE_ALPHA_STEP = 2
SURROGATE_DPS = 60
SURROGATE_S = ("0.5", "2.0")
SURROGATE_N = (1, 2, 3, 5)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {relative(path)}")
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


def kernel_residual(z: mp.mpc, terms: int) -> mp.mpc:
    prefix = mp.mpc(0)
    for k in range(terms):
        prefix += mp.e ** (1j * mp.pi * (2 * k + 1) * z)
    return 1 / mp.sin(mp.pi * z) - (-2j * prefix + mp.e ** (2j * mp.pi * terms * z) / mp.sin(mp.pi * z))


def shifted_contour_integral(center: mp.mpf, terms: int, s: mp.mpc) -> mp.mpc:
    rotation = mp.e ** (-1j * mp.pi / 4)

    def integrand(q: mp.mpf) -> mp.mpc:
        u = center + q * rotation
        return (
            mp.mpf(1) / (2j)
            * mp.e ** (-1j * mp.pi * u * u)
            * mp.power(u + terms, -s)
            / mp.sin(mp.pi * u)
            * rotation
        )

    return mp.quad(integrand, [-mp.inf, -4, -1, 0, 1, 4, mp.inf])


def surrogate_rows() -> list[dict[str, str | int]]:
    mp.mp.dps = SURROGATE_DPS
    s = mp.mpf(SURROGATE_S[0]) + 1j * mp.mpf(SURROGATE_S[1])
    rows = []
    for terms in SURROGATE_N:
        left = shifted_contour_integral(mp.mpf("0.5") - terms, terms, s)
        right = shifted_contour_integral(mp.mpf("0.5"), terms, s)
        dirichlet = mp.fsum(mp.power(n, -s) for n in range(1, terms + 1))
        discrepancy = left - right - dirichlet
        rows.append(
            {
                "N": terms,
                "left_shifted_contour": mp.nstr(left, 52),
                "right_standard_contour": mp.nstr(right, 52),
                "Dirichlet_prefix": mp.nstr(dirichlet, 52),
                "residue_identity_discrepancy_absolute": mp.nstr(abs(discrepancy), 18),
            }
        )
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    roster = artifact["actual_window_identity"]
    numerical = artifact["surrogate_contour_check"]
    return f"""# Finite geometric RSI contour residual and pole release

Date: 2026-08-27

Status: exact finite identity certified; not a proof of the continuation defect

The A21 infinite interchange is unnecessary for every finite odd-label
prefix.  For any integer `N>=0` and any non-pole `z`, finite geometric algebra
gives

```text
csc(pi z)
 =-2i sum_(k=0)^(N-1) exp(i*pi*(2k+1)*z)
  +exp(2*pi*i*N*z)csc(pi*z).                         (FR1)
```

This is an identity of meromorphic functions, not an asymptotic expansion.
Insert (FR1) directly into the exact Riemann-Siegel contour before any
infinite summation.  If `P_N(s)` is the resulting finite odd-exponential
prefix and

```text
E_N(s)=1/(2i) integral_C
       exp(-i*pi*z^2+2*pi*i*N*z) z^(-s)csc(pi*z) dz, (FR2)
```

then exactly

```text
RSI(s)=P_N(s)+E_N(s).                                (FR3)
```

Translate `u=z-N`.  Since

```text
exp(i*pi*N^2)=(-1)^N,
sin(pi(u+N))=(-1)^N sin(pi*u),                       (FR4)
```

the signs cancel and

```text
E_N(s)=1/(2i) integral_(C-N)
       exp(-i*pi*u^2)(u+N)^(-s)csc(pi*u) du.         (FR5)
```

Move `C-N` back to `C`.  The branch point `u=-N` remains to the left of the
deformation strip.  Exactly the poles `m=1-N,...,0` are crossed, and

```text
Res_(u=m) [exp(-i*pi*u^2)(u+N)^(-s)csc(pi*u)]
 =(m+N)^(-s)/pi.                                     (FR6)
```

Therefore

```text
E_N(s)=sum_(n=1)^N n^(-s)+C_N(s),

C_N(s)=1/(2i) integral_C
       exp(-i*pi*u^2)(u+N)^(-s)csc(pi*u) du.         (FR7)
```

The altered complex-height numerical contour check has maximum discrepancy
`{numerical['maximum_residue_identity_discrepancy_absolute']}` for
`N=1,2,3,5`; it confirms the residue orientation independently of the finite
geometric algebra.

For the actual odd window, put

```text
n_-=(A-1)/2={roster['prefix_before_window']},
n_+=(B+1)/2={roster['prefix_through_window']}.        (FR8)
```

There are `{roster['window_label_count']}` labels, and subtraction of two
copies of (FR3) gives the exact kernel-level common-contour identity

```text
P_[A,B]=P_(n_+)-P_(n_-)
       =E_(n_-)-E_(n_+)
       =-sum_(n=n_-+1)^(n_+) n^(-s)
        +C_(n_-)-C_(n_+).                            (FR9)
```

Thus the finite window is already an exact difference of two translated
residual contours plus an explicit pole block.  No claim that the A21
infinite remainder vanishes is needed.  The next gate must align the finite
exponential-prefix normalization and branches with the corrected physical
`Q_K`, then combine the released pole block with `G` and `A_transition`
before taking norms.

Pi provenance: every `pi` in (FR1)--(FR7) comes from the original
`sin(pi z)` denominator, Gaussian `exp(-i*pi*z^2)`, and integer contour
translation.  No geometric circle construction or fitted constant is used.

Proof boundary: exact finite meromorphic identity, contour translation,
residue inventory, actual-window index arithmetic, and a surrogate numerical
orientation check only.  The finite exponential integrals have not yet been
matched branch-for-branch to the corrected physical `Q_K`; no `C_K`, `D_K`,
`Delta_KU`, `J_Z`, or non-A enclosure, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["A21_provenance"]["decision"]["A21_permitted_as_exact_whole_contour_interchange"]
        is False,
        "A21 nonpromotion guard drift",
    )

    require(SOURCE_ALPHA_MIN % 2 == 1 and SOURCE_ALPHA_MAX % 2 == 1, "source endpoints are not odd")
    prefix_before = (SOURCE_ALPHA_MIN - 1) // 2
    prefix_through = (SOURCE_ALPHA_MAX + 1) // 2
    label_count = prefix_through - prefix_before
    require(label_count == 2_481_423, "actual source count drift")
    require(SOURCE_ALPHA_MIN == 2 * prefix_before + 1, "lower prefix map failed")
    require(SOURCE_ALPHA_MAX == 2 * prefix_through - 1, "upper prefix map failed")

    mp.mp.dps = 90
    kernel_rows = []
    for z_real, z_imag in (("0.37", "0.2"), ("0.5", "-0.31"), ("-0.2", "0.7")):
        z = mp.mpc(z_real, z_imag)
        z_text = f"{z_real}{'+' if not z_imag.startswith('-') else ''}{z_imag}j"
        for terms in (1, 2, 5, 11):
            discrepancy = kernel_residual(z, terms)
            require(abs(discrepancy) < mp.mpf("1e-80"), "finite kernel identity numerical guard failed")
            kernel_rows.append(
                {
                    "z": z_text,
                    "N": terms,
                    "discrepancy_absolute": mp.nstr(abs(discrepancy), 18),
                }
            )

    rows = surrogate_rows()
    max_discrepancy = max(mp.mpf(str(row["residue_identity_discrepancy_absolute"])) for row in rows)
    require(max_discrepancy < mp.mpf("1e-48"), "surrogate pole-release orientation check failed")

    artifact = {
        "kind": STEM,
        "status": "finite_RSI_geometric_prefix_residual_translation_and_pole_release_identity_certified",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "s": "1/2+i*t",
            "contour": "C: z=1/2+q*exp(-i*pi/4), q in R",
        },
        "finite_kernel_identity": {
            "formula": "csc(pi*z)=-2*i*sum_(k=0)^(N-1)exp(i*pi*(2*k+1)*z)+exp(2*pi*i*N*z)*csc(pi*z)",
            "valid_for": "integer N>=0 and z not an integer pole",
            "meromorphic_identity": True,
            "kernel_numerical_rows": kernel_rows,
        },
        "exact_RSI_prefix_residual": {
            "prefix": "P_N=(1/(2*i))*integral_C exp(-i*pi*z^2)z^(-s)[-2*i*sum_(k=0)^(N-1)exp(i*pi*(2*k+1)*z)]dz",
            "residual": "E_N=(1/(2*i))*integral_C exp(-i*pi*z^2+2*pi*i*N*z)z^(-s)csc(pi*z)dz",
            "identity": "RSI=P_N+E_N",
            "infinite_sum_or_interchange_used": False,
        },
        "translated_residual_and_poles": {
            "translated_residual": "E_N=(1/(2*i))*integral_(C-N)exp(-i*pi*u^2)(u+N)^(-s)csc(pi*u)du",
            "crossed_poles": "m=1-N,...,0",
            "residue_at_m": "(m+N)^(-s)/pi",
            "standard_contour_remainder": "C_N=(1/(2*i))*integral_C exp(-i*pi*u^2)(u+N)^(-s)csc(pi*u)du",
            "pole_release_identity": "E_N=sum_(n=1)^N n^(-s)+C_N",
            "branch_point_minus_N_outside_deformation_strip": True,
        },
        "actual_window_identity": {
            "alpha_min": SOURCE_ALPHA_MIN,
            "alpha_max": SOURCE_ALPHA_MAX,
            "alpha_step": SOURCE_ALPHA_STEP,
            "prefix_before_window": prefix_before,
            "prefix_through_window": prefix_through,
            "window_label_count": label_count,
            "released_Dirichlet_block": [prefix_before + 1, prefix_through],
            "identity": "P_[A,B]=E_(n_-)-E_(n_+)=-sum_(n=n_-+1)^(n_+)n^(-s)+C_(n_-)-C_(n_+)",
        },
        "surrogate_contour_check": {
            "decimal_digits": SURROGATE_DPS,
            "s": f"{SURROGATE_S[0]}+{SURROGATE_S[1]}i",
            "rows": rows,
            "maximum_residue_identity_discrepancy_absolute": mp.nstr(max_discrepancy, 18),
            "purpose": "altered complex-height orientation and residue-sign check only",
        },
        "decision": {
            "A21_infinite_interchange_required_for_finite_prefix": False,
            "finite_prefix_plus_residual_exact": True,
            "integer_translation_releases_Dirichlet_prefix_exactly": True,
            "actual_Kummer_window_has_common_contour_residual_difference": True,
            "physical_QK_prefactor_and_branch_alignment_complete": False,
            "released_pole_block_joined_to_G_and_A_transition": False,
            "D_K_enclosed": False,
            "non_A_bound_proved": False,
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
            "numerical_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Derive the finite exponential-prefix integral in the corrected equation-(4)/(9) Kummer normalization with a branch-explicit constant check. Then project the actual n_-=79788 and n_+=2561211 residual difference, join its released Dirichlet block to G and A_transition algebraically, and isolate the remaining common-contour representative of D_K before any norm.",
        "proof_boundary": "Exact finite meromorphic kernel identity, finite RSI prefix/residual split, integer contour translation, pole release, actual-window index arithmetic, and an altered-height numerical orientation check only. No physical Q_K branch/prefactor match, C_K, D_K, Delta_KU, J_Z, non-A enclosure, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified finite RSI geometric residual and exact Dirichlet pole release", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
