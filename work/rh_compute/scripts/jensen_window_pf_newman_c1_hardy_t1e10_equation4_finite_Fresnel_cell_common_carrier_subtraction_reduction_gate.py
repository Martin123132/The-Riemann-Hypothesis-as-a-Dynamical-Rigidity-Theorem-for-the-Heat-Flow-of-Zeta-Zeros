#!/usr/bin/env python3
"""Put the finite Fresnel branch and owned physical carriers on one mode lattice."""

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
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]

STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_common_carrier_subtraction_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "Fresnel_connection": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_branch_Fresnel_Dirichlet_connection_and_mode_roster_gate.json",
    "finite_residual": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_geometric_contour_residual_and_pole_release_gate.json",
    "physical_ownership": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_physical_transform_ownership_ledger_gate.json",
    "Gamma_bulk": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "A_transition": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate.json",
    "exact_H": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_continuation_defect_target_gate.json",
}

HEIGHT = 10_000_000_000
SOURCE_ALPHA_MIN = 159_577
SOURCE_ALPHA_MAX = 5_122_421
SOURCE_COUNT = 2_481_423
OWNED_FIRST = 622
ORDINARY_LAST = 39_852
A_FIRST = 39_853
TARGET_LAST = 39_894
EXTENDED_LAST = 39_936


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


def source(alpha_min: int, count: int, y: mp.mpf | mp.mpc) -> mp.mpc:
    return mp.fsum(mp.e ** (1j * mp.pi * (alpha_min + 2 * j) * y) for j in range(count))


def cell_witnesses() -> dict[str, Any]:
    mp.mp.dps = 85
    t = mp.mpf("2.75")
    s = mp.mpf("0.5") + 1j * t
    alpha_min = 1
    count = 4

    covariance_rows = []
    for mode, u_text in ((1, "-0.31"), (2, "0.17"), (5, "0.41")):
        u = mp.mpf(u_text)
        direct = mp.e ** (-1j * mp.pi * (mode + u) ** 2) * source(alpha_min, count, mode + u)
        reduced = mp.e ** (-1j * mp.pi * u * u - 2j * mp.pi * mode * u) * source(alpha_min, count, u)
        covariance_rows.append(
            {
                "mode": mode,
                "u": u_text,
                "chirp_source_covariance_discrepancy_absolute": mp.nstr(abs(direct - reduced), 25),
            }
        )

    cell_rows = []
    for mode in (1, 2, 5):
        direct = mp.quad(
            lambda y: y ** (-s) * mp.e ** (-1j * mp.pi * y * y) * source(alpha_min, count, y),
            [mp.mpf(mode) - mp.mpf("0.5"), mode, mp.mpf(mode) + mp.mpf("0.5")],
        )
        beta = mp.quad(
            lambda u: (1 + u / mode) ** (-s)
            * mp.e ** (-1j * mp.pi * u * u - 2j * mp.pi * mode * u)
            * source(alpha_min, count, u),
            [mp.mpf("-0.5"), 0, mp.mpf("0.5")],
        )
        reduced = mode ** (-s) * beta
        cell_rows.append(
            {
                "mode": mode,
                "cell_factorization_discrepancy_absolute": mp.nstr(abs(direct - reduced), 25),
            }
        )

    max_covariance = max(mp.mpf(row["chirp_source_covariance_discrepancy_absolute"]) for row in covariance_rows)
    max_cell = max(mp.mpf(row["cell_factorization_discrepancy_absolute"]) for row in cell_rows)
    require(max_covariance < mp.mpf("1e-75"), "finite-source chirp covariance failed")
    require(max_cell < mp.mpf("1e-70"), "fresh Fresnel cell factorization failed")
    return {
        "surrogate_height": mp.nstr(t, 10),
        "surrogate_odd_roster": [1, 7, 2],
        "covariance_rows": covariance_rows,
        "cell_rows": cell_rows,
        "maximum_covariance_discrepancy_absolute": mp.nstr(max_covariance, 25),
        "maximum_cell_discrepancy_absolute": mp.nstr(max_cell, 25),
    }


def gamma_lift_witnesses() -> list[dict[str, Any]]:
    mp.mp.dps = 85
    t = mp.mpf("3.25")
    s = mp.mpf("0.5") + 1j * t
    theta = mp.im(mp.loggamma(mp.mpf("0.25") + 0.5j * t)) - t * mp.log(mp.pi) / 2
    theta_zero = t * (mp.log(t / (2 * mp.pi)) - 1) / 2 - mp.pi / 8
    exponent = mp.mpf("0.75") + 0.5j * t
    i_bulk = mp.e ** (0.5j * t) * mp.gamma(exponent) / (0.5j * t) ** exponent
    i_gaussian = 2 * mp.sqrt(mp.pi / t) * mp.e ** (-1j * mp.pi / 4)
    bulk_factor = i_bulk / i_gaussian
    h_factor = (
        2 ** mp.mpf("1.25")
        * t ** mp.mpf("0.25")
        * mp.pi ** mp.mpf("1.5")
        * mp.e ** (-3 * mp.pi * t / 4)
        / (
            (1 + mp.e ** (-2 * mp.pi * t))
            * abs(mp.gamma(mp.mpf("0.25") + 0.5j * t))
            * abs(mp.gamma(mp.mpf("0.75") + 0.5j * t)) ** 2
        )
    )
    common = h_factor * mp.e ** (-1j * theta) * bulk_factor * mp.e ** (1j * theta_zero)

    rows = []
    for mode in (2, 5, 11):
        complex_lift = common * mode ** (-s)
        hardy_lift = 2 * mp.re(mp.e ** (1j * theta) * complex_lift)
        physical = h_factor * 2 * mp.re(
            bulk_factor * mp.e ** (1j * (theta_zero - t * mp.log(mode))) / mp.sqrt(mode)
        )
        rows.append(
            {
                "mode": mode,
                "Hardy_lift_discrepancy_absolute": mp.nstr(abs(hardy_lift - physical), 25),
            }
        )
    require(
        max(mp.mpf(row["Hardy_lift_discrepancy_absolute"]) for row in rows) < mp.mpf("1e-75"),
        "Gamma common-coefficient Hardy lift failed",
    )
    return rows


def algebra_and_rosters() -> dict[str, Any]:
    qk, g_target, g_extra, a_endpoint, h = sp.symbols("Q_K G_target G_extra A_endpoint H")
    a_transition = a_endpoint + g_extra
    jz = qk - g_target - a_transition
    extended_form = qk - (g_target + g_extra) - a_endpoint
    require(sp.expand(jz - extended_form) == 0, "extended physical ownership reassembly failed")
    require(ORDINARY_LAST - OWNED_FIRST + 1 == 39_231, "ordinary count drift")
    require(A_FIRST - OWNED_FIRST == 39_231, "ordinary/A boundary drift")
    require(EXTENDED_LAST - A_FIRST + 1 == 84, "A-window count drift")
    require(TARGET_LAST - A_FIRST + 1 == 42, "target-side A count drift")
    require(EXTENDED_LAST - TARGET_LAST == 42, "outer-side A count drift")
    return {
        "physical_reassembly": "J_Z=Q_K-G_target-A_transition=Q_K-G_extended-A_endpoint",
        "A_transition": "A_transition=A_endpoint+G_extra",
        "G_extended": "G_extended=sum_(m=622)^39936 G_m",
        "ordinary_owned_cells": [OWNED_FIRST, ORDINARY_LAST],
        "ordinary_owned_count": 39_231,
        "A_transition_cells": [A_FIRST, EXTENDED_LAST],
        "A_transition_count": 84,
        "A_transition_halves": [[A_FIRST, TARGET_LAST], [TARGET_LAST + 1, EXTENDED_LAST]],
        "A_transition_half_counts": [42, 42],
        "unowned_lower_cells": [1, OWNED_FIRST - 1],
        "unowned_upper_cells": [EXTENDED_LAST + 1, "infinity"],
        "symbolic_HJZ_reassembly_zero": str(sp.expand(h * jz - (h * qk - h * (g_target + g_extra) - h * a_endpoint))),
    }


def render_note(artifact: dict[str, Any]) -> str:
    return f"""# Finite Fresnel cells and common-carrier subtraction

Date: 2026-08-27

Status: exact cell-amplitude reduction certified; joined amplitude bound open

From Section 11.486, write

```text
D_W(y)=sum_(j=0)^(M-1)exp[i*pi*(A+2j)y],
B_W=integral_0^infinity[Fresnel]
    y^(-s)exp(-i*pi*y^2)D_W(y)dy.                    (CC1)
```

For every integer `m` and `-1/2<=u<=1/2`, oddness of `A` gives

```text
D_W(m+u)=(-1)^m D_W(u),
exp[-i*pi*(m+u)^2]D_W(m+u)
 =exp(-i*pi*u^2-2*pi*i*m*u)D_W(u).                  (CC2)
```

Define the central half-cell and positive mode amplitudes

```text
B_0=integral_0^(1/2)y^(-s)exp(-i*pi*y^2)D_W(y)dy,

beta_m(t)=integral_(-1/2)^(1/2)
 (1+u/m)^(-s)exp(-i*pi*u^2-2*pi*i*m*u)D_W(u)du,
                                                               (CC3)
B_m=m^(-s)beta_m(t),       m>=1.                     (CC4)
```

Contiguous partition of the Fresnel/Abel limit, without exchanging an
infinite series with an integral, proves

```text
B_W=B_0+sum_(m=1)^infinity[Fresnel]m^(-s)beta_m(t).  (CC5)
```

The altered finite-roster quadrature witnesses have maximum covariance and
cell-factorization discrepancies `{artifact['cell_witnesses']['maximum_covariance_discrepancy_absolute']}`
and `{artifact['cell_witnesses']['maximum_cell_discrepancy_absolute']}`.

The exact global-Morse Gamma carrier has

```text
G_m=2 Re[B(t)exp(i*theta_0)m^(-s)],
B(t)=I_bulk(t)/I_G(t).                               (CC6)
```

Put

```text
C_G(t)=H(t)exp[-i*theta(t)]B(t)exp(i*theta_0).       (CC7)
```

Then, exactly and mode by mode,

```text
Hardy_t[C_G(t)m^(-s)]=H(t)G_m.                       (CC8)
```

Thus the ordinary owned cells are governed by the explicit coefficient
defect `beta_m-C_G`; this is subtraction before norms, not comparison of two
separately accumulated large values.

For `39853<=m<=39936`, let `a_m` denote the already-defined real physical
projection of `A_m+A_-m` and use the canonical Hardy lift

```text
mathcal A_m=(H(t)/2)exp[-i*theta(t)]a_m,
Hardy_t[mathcal A_m]=H(t)a_m.                        (CC9)
```

Because `A_transition=A_endpoint+G_extra`, the exact joined target is

```text
H(t)J_Z=Hardy_t{{
 P_W+B_0
 +sum_(m=1)^621 m^(-s)beta_m
 +sum_(m=622)^39852 m^(-s)[beta_m-C_G]
 +sum_(m=39853)^39936
       [m^(-s)(beta_m-C_G)-mathcal A_m]
 +sum_(m=39937)^infinity[Fresnel]m^(-s)beta_m}}.     (CC10)
```

Here

```text
P_W=-sum_(n=79789)^2561211 n^(-s)
    +C_79788-C_2561211.                              (CC11)
```

Equation (CC10) is the requested amplitude-level common lattice.  It does not
say that matching saddles makes `beta_m-C_G` or the A-corrected coefficient
small.  The next obligation is twofold: derive the natural integral lift of
`mathcal A_m` on the same cells, then join `P_W` to the lower and upper
unowned-cell packages through an exact Mordell/contour transformation before
attempting an interval norm.

The mode ownership is now algebraically disjoint:

```text
1..621:         unowned lower cells,
622..39852:     39231 ordinary Gamma-subtracted cells,
39853..39936:   84 Gamma-plus-A-subtracted cells (42+42),
39937..infinity:unowned upper cells.                 (CC12)
```

Pi provenance: every `pi` in (CC1)--(CC12) descends from the finite RSI
Gaussian, odd Fourier roster, Riemann-Siegel phase, or the exact global-Morse
Gamma normalization.  No fitted constant is introduced.

Proof boundary: exact finite-source chirp covariance, contiguous Fresnel-cell
partition, Gamma common-coefficient lift, A canonical Hardy lift, and joined
`H J_Z` reduction only.  No bound for any coefficient defect, natural
cell-integral A lift, Mordell join of `P_W` with unowned cells, actual-height
`J_Z` or `D_K` enclosure, non-A, all-height, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")

    ownership = dependencies["physical_ownership"]["identity_certificate"]
    require(
        ownership["joined_source_owned_transform"]
        == "J_Z=P_t[Z]=Q_K-G-A_transition=R_KGamma-A_transition",
        "physical ownership identity drift",
    )
    fresnel = dependencies["Fresnel_connection"]["finite_Fresnel_branch_sum"]
    require(fresnel["integer_singularities_removable_by_finite_sum"] is True, "finite source removability lost")

    artifact = {
        "kind": STEM,
        "status": "exact_Fresnel_cell_common_Gamma_and_A_carrier_subtraction_reduction_certified_no_amplitude_bound",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "source_labels": [SOURCE_ALPHA_MIN, SOURCE_ALPHA_MAX, 2],
            "source_label_count": SOURCE_COUNT,
            "owned_mode_range": [OWNED_FIRST, EXTENDED_LAST],
        },
        "exact_cell_reduction": {
            "source": "D_W(y)=sum_(j=0)^(M-1)exp(i*pi*(A+2j)y)",
            "covariance": "exp(-i*pi*(m+u)^2)D_W(m+u)=exp(-i*pi*u^2-2*pi*i*m*u)D_W(u)",
            "central_half_cell": "B_0=integral_0^(1/2)y^(-s)exp(-i*pi*y^2)D_W(y)dy",
            "mode_amplitude": "beta_m=integral_(-1/2)^(1/2)(1+u/m)^(-s)exp(-i*pi*u^2-2*pi*i*m*u)D_W(u)du",
            "mode_cell": "B_m=m^(-s)beta_m",
            "partition": "B_W=B_0+sum_(m=1)^infinity[Fresnel]m^(-s)beta_m",
            "construction": "contiguous finite cell partition followed by the original Fresnel limit; no infinite sum-integral interchange",
        },
        "Gamma_common_carrier": {
            "physical_mode": "G_m=2Re[B(t)exp(i*theta_0)m^(-s)]",
            "common_coefficient": "C_G=H(t)exp(-i*theta(t))B(t)exp(i*theta_0)",
            "Hardy_lift": "Hardy_t[C_G*m^(-s)]=H(t)G_m",
            "ordinary_cell_defect": "m^(-s)[beta_m-C_G]",
        },
        "A_common_lift": {
            "physical_mode": "a_m=P_t[A_m+A_-m]",
            "canonical_lift": "mathcal_A_m=(H(t)/2)exp(-i*theta(t))*a_m",
            "Hardy_lift": "Hardy_t[mathcal_A_m]=H(t)a_m",
            "natural_same_cell_integral_lift_derived": False,
        },
        "joined_identity": {
            "physical": "J_Z=Q_K-G_extended-A_endpoint",
            "complex_Hardy_argument": "P_W+B_0+sum_(1..621)m^(-s)beta_m+sum_(622..39852)m^(-s)(beta_m-C_G)+sum_(39853..39936)[m^(-s)(beta_m-C_G)-mathcal_A_m]+sum_(39937..infinity)m^(-s)beta_m",
            "result": "H(t)J_Z=Hardy_t[complex_Hardy_argument]",
            "finite_prefix": "P_W=-sum_(n=79789)^2561211 n^(-s)+C_79788-C_2561211",
        },
        "mode_rosters": algebra_and_rosters(),
        "cell_witnesses": cell_witnesses(),
        "Gamma_lift_witnesses": gamma_lift_witnesses(),
        "decision": {
            "finite_Fresnel_cell_partition_exact": True,
            "Gamma_carrier_subtracted_on_common_m_minus_s_coefficient": True,
            "A_carrier_has_exact_canonical_Hardy_lift": True,
            "natural_A_same_cell_integral_lift_derived": False,
            "P_W_joined_to_unowned_cell_packages": False,
            "owned_cell_amplitude_defects_bounded": False,
            "actual_height_J_Z_enclosed": False,
            "actual_height_D_K_enclosed": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "references": {
            "Gabcke_parabolic_RSI": "https://arxiv.org/abs/1512.01186",
            "Kuznetsov_truncated_theta_Mordell": "https://arxiv.org/abs/1306.4081",
            "reference_use": "route context for the next exact Mordell join; neither source is used to claim the open amplitude bound",
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
        "next_obligation": "Derive the natural exact same-cell integral lift of each paired A endpoint mode, then use an exact Mordell/contour transformation to join P_W with modes 1..621 and 39937..infinity. Bound only the resulting joined coefficient packets, never beta_m, P_W, or the carriers separately.",
        "proof_boundary": "Exact finite-source chirp covariance, contiguous Fresnel-cell partition, Gamma common-coefficient lift, A canonical Hardy lift, and joined H*J_Z reduction only. No coefficient-defect bound, natural same-cell A lift, Mordell join of P_W to unowned cells, actual-height J_Z or D_K enclosure, non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact Fresnel-cell common-carrier subtraction reduction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
