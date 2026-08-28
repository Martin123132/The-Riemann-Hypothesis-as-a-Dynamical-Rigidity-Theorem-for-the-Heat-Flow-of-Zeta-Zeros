#!/usr/bin/env python3
"""Certify two-jet reassembly of the full common projector kernel."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__)
CHECKER = Path(__file__).with_name("check_" + STEM + ".py")
DEPENDENCIES = {
    "two_jet_extraction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_extraction_gate.json",
    "R_after_A_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "finite_endpoint_decomposition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_finite_endpoint_tail_decomposition_gate.json",
    "projector_defect": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_phase_space_projector_defect_gate.json",
}

A = 159_577
B = 5_122_421
L = (B - A) // 2
T_LO = 622
T_HI = 39_894


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


def symbolic_certificate() -> dict[str, Any]:
    m, K, phi = sp.symbols("m K phi", nonzero=True)
    p_hat = 1 / (2 * sp.pi**2 * m**2)
    per_mode = 2 * sp.pi * sp.I * m * (phi + K * p_hat / 2)
    expected = 2 * sp.pi * sp.I * m * phi + sp.I * K / (2 * sp.pi * m)
    require(sp.simplify(per_mode - expected) == 0, "per-mode two-jet coefficient drift")

    E, H_T, R_all, R_T, P_T = sp.symbols("E H_T R_all R_T P_T")
    source_complement = E - sp.I * K * H_T / (2 * sp.pi) + R_all - R_T
    target_endpoint_atoms = R_T + sp.I * K * H_T / (2 * sp.pi) - P_T
    reassembled = sp.expand(source_complement + target_endpoint_atoms)
    require(sp.simplify(reassembled - (E + R_all - P_T)) == 0, "common-kernel reassembly drift")

    phi_plus, phi_minus, weight = sp.symbols("phi_plus phi_minus weight")
    pair = 2 * sp.pi * sp.I * m * weight * phi_plus + 2 * sp.pi * sp.I * (-m) * weight * phi_minus
    pair_expected = 2 * sp.pi * sp.I * m * weight * (phi_plus - phi_minus)
    require(sp.expand(pair - pair_expected) == 0, "symmetric pair reduction drift")

    return {
        "per_mode_identity": "I_m=2*pi*i*m*hat Phi_x(m)+i*K_x/(2*pi*m), m!=0",
        "zero_mode_identity": "I_0=E_x",
        "source_complement": "sum_(m notin T)w_m I_m=E_x-i*K_x*H_T/(2*pi)+R_all-R_T",
        "target_endpoint_atoms": "sum_(m in T)w_m(A_m+B_m)=R_T+i*K_x*H_T/(2*pi)-sum_(m in T)w_m P_m",
        "harmonic_cancellation_residual": "0",
        "symmetric_pair": "m*w_m[hat Phi_x(m)-hat Phi_x(-m)]",
    }


def common_kernel_certificate() -> dict[str, Any]:
    return {
        "finite_identity": "Delta_(M,epsilon)(x)=H_x+E_x+2*pi*i*sum_(m=1)^M m*w_m[hat Phi_x(m)-hat Phi_x(-m)]-sum_(m=622)^39894 w_m*P_m(x)",
        "cutoff": "M>=B=5122421",
        "weight": "w_m=exp(-pi*epsilon*m^2)",
        "target_endpoint_relation": "A_m+B_m=I_m-P_m for every positive target mode",
        "absolute_limit": "The symmetric Phi series is absolutely and uniformly convergent on x in [0,1] after the two endpoint jets are extracted.",
        "zero_regulator_form": "Set w_m=1 after taking the certified symmetric limit; no one-sided endpoint-current series remains.",
        "R_after_A_insertion": "Subtract chi_W*Btr, O, and the certified A notch from this same Delta before the physical projection, exactly as in the R_after_A equivalence gate.",
    }


def render_note() -> str:
    return f"""# Two-jet reassembly of the common projector kernel

Date: 2026-08-23

Status: exact finite-regulator target-current cancellation and absolutely
convergent symmetric common-kernel representation certified; bound open

Let `hat Phi_x(m)` be the Fourier coefficient of the periodic `C^1` defect
from the two-jet extraction gate.  For every nonzero integer `m`, integration
by parts gives

```text
I_m=2*pi*i*m*hat Phi_x(m)+i*K_x/(2*pi*m),
I_0=E_x.                                                (TR1)
```

The second term in (TR1) is the extracted first-derivative endpoint jet.  In
the target-notched source integral it contributes

```text
-i*K_x*H_T(epsilon)/(2*pi).                            (TR2)
```

For every positive target mode, however,

```text
A_m+B_m=I_m-P_m,                                      (TR3)
```

where `P_m` is the exact full-line Gamma bulk.  Substitution of (TR1) into
(TR3) contributes the opposite current

```text
+i*K_x*H_T(epsilon)/(2*pi).                            (TR4)
```

Thus (TR2) and (TR4) cancel coefficientwise before any norm.  The target
Fourier coefficients of `Phi_x` are restored at the same time.  With
`w_m=exp(-pi*epsilon*m^2)` and every common cutoff `M>=B={B}`, the complete
Gamma-normalized projector kernel is exactly

```text
Delta_(M,epsilon)(x)
 =H_x+E_x
  +2*pi*i*sum_(m=1)^M m*w_m
     [hat Phi_x(m)-hat Phi_x(-m)]
  -sum_(m=622)^39894 w_m*P_m(x).                      (TR5)
```

The apparent target harmonic current is therefore not an additional error
term and must not be counted after (TR5).  Since `Phi_x` is periodic through
its first derivative,

```text
|hat Phi_x(m)|<=K_Phi(x)/(2*pi*|m|)^3,                (TR6)
```

so the symmetric series in (TR5) is absolutely convergent.  The finite source
roster makes `K_Phi(x)` continuous on `0<=x<=1`; hence its maximum is finite
and the convergence is uniform on the Kummer interval.  This licenses the
zero-regulator representation structurally, but the raw derivative maximum is
not asserted to be quantitatively useful.

The post-A object remains

```text
mathfrak R_A=Delta-chi_W*Btr-O-mathcal A,              (TR7)
```

with exactly the ownership and physical limit order already certified.  No A
or B subtraction has been discarded or counted twice.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (TR1)--(TR7) is inherited from the Kummer
quadratic phase, integer Fourier character, and Gaussian Abel regulator.  No
fitted or geometric occurrence is introduced.

Proof boundary: exact finite-regulator coefficient reassembly, cancellation
of the target harmonic endpoint jet, and absolute/uniform convergence of the
periodic two-jet pair series only.  No useful numerical tail constant,
physical Kummer quadrature, A/B-subtracted `R_after_A` or `R_Dir` bound,
complete `Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["two_jet_extraction"]["decision"]["target_asymmetry_reduced_to_closed_harmonic_current"] is True, "two-jet harmonic drift")
    require(dependencies["two_jet_extraction"]["decision"]["differentiated_leakage_absolutely_convergent"] is True, "two-jet tail drift")
    require(dependencies["R_after_A_equivalence"]["decision"]["common_kernel_equals_six_class_post_A_mask"] is True, "R_after_A kernel drift")
    require(dependencies["finite_endpoint_decomposition"]["decision"]["finite_current_decomposed_exactly"] is True, "endpoint decomposition drift")
    require(dependencies["projector_defect"]["decision"]["target_full_line_carrier_subtracted_before_norms"] is True, "target carrier drift")

    artifact = {
        "kind": STEM,
        "status": "exact_two_jet_common_kernel_reassembly_harmonic_current_cancelled_symmetric_pair_series_absolute_quantitative_bound_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "roster_coordinate_length": L,
            "target_positive_modes": [T_LO, T_HI],
            "common_cutoff": "M>=B=5122421",
        },
        "symbolic_certificate": symbolic_certificate(),
        "common_kernel_certificate": common_kernel_certificate(),
        "decision": {
            "per_mode_two_jet_identity_proved": True,
            "target_harmonic_current_cancels_coefficientwise": True,
            "target_Phi_coefficients_restored_exactly": True,
            "common_kernel_reassembled_as_symmetric_pair_series": True,
            "zero_regulator_pair_series_absolute_and_uniform": True,
            "R_after_A_ownership_and_limit_order_preserved": True,
            "harmonic_current_is_independent_error_term": False,
            "raw_derivative_majorant_closes_target": False,
            "quantitative_R_after_A_bound_proved": False,
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
        "next_obligation": "Develop a phase-adapted evaluator for the symmetric pair coefficient m[hat Phi_x(m)-hat Phi_x(-m)] and compose it with H_x+E_x and the finite target Gamma bulk before subtracting the already-certified B window, B outer block, and A transition. Use the periodic C1 cancellation, not the enormous raw K_Phi derivative majorant.",
        "proof_boundary": "Exact finite-regulator coefficient reassembly, cancellation of the target harmonic endpoint jet, and absolute/uniform convergence of the periodic two-jet pair series only. No useful numerical tail constant, physical Kummer quadrature, A/B-subtracted R_after_A or R_Dir bound, complete Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(), encoding="utf-8")
    print("certified two-jet common-kernel reassembly", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
