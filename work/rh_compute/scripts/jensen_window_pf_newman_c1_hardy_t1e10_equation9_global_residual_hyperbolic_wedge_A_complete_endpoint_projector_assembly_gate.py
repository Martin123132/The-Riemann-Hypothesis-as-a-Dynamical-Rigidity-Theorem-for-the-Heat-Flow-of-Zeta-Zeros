#!/usr/bin/env python3
"""Assemble the complete exact A endpoint and its projector jump."""

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

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "local_affine": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_exact_domain_gate.json",
    "local_nonlinear": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_local_exact_nonlinear_amplitude_gate.json",
    "exact_exterior": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_rational_common_phase_contour_gate.json",
    "gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "orientation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate.json",
    "projector": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_phase_space_projector_defect_gate.json",
}

HEIGHT = 10_000_000_000
ENDPOINT = 159_577
LOWER_START = 39_853
TARGET_END = 39_894
UPPER_END = 39_936
PRECISION = 90


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


def complex_from_record(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def symbolic_certificate() -> dict[str, str]:
    P_plus, P_minus, Q_plus, Q_minus, chi = sp.symbols("P_plus P_minus Q_plus Q_minus chi")
    finite_pair = (P_plus + Q_plus) + (P_minus + Q_minus) - chi * P_plus
    projector_form = (1 - chi) * P_plus + P_minus + Q_plus + Q_minus
    require(sp.expand(finite_pair - projector_form) == 0, "paired-projector identity failed")
    return {
        "paired_projector": "I_m+I_-m-chi_T(m)P_m=(1-chi_T(m))P_m+P_-m+Q_m+Q_-m",
        "endpoint_coefficient": "Q_m+Q_-m has coefficient +1 for every mode 39853..39936",
        "positive_bulk_coefficient": "1-chi_T(m)=0 for m<=39894 and 1 for m>=39895",
        "exact_endpoint_assembly": "A_endpoint=A_local_affine+A_local_exact_minus_affine+A_exact_exterior",
        "projector_completed_transition": "A_transition=A_endpoint+sum_(m=39895)^39936 G_m",
        "Gamma_mode": "G_m=2 Re[B(t) exp(i(theta_0-t log m))]/sqrt(m)",
    }


def assembly_certificate(dependencies: dict[str, dict[str, Any]], precision: int) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    local_affine = arb(
        dependencies["local_affine"]["certificate"]["signed_sums"]["complete_physical_local_exact_affine_ball"]
    )
    local_nonlinear = arb(
        dependencies["local_nonlinear"]["certificate"]["signed_sums"]["complete_physical_local_exact_minus_affine_ball"]
    )
    exact_exterior = arb(dependencies["exact_exterior"]["certificate"]["physical_full_exact_exterior_ball"])
    complete_endpoint = local_affine + local_nonlinear + exact_exterior

    gamma = complex_from_record(dependencies["gamma_bulk"]["certificate"]["Gamma_bulk_factor_ball"])
    theta_correction = arb(dependencies["gamma_bulk"]["certificate"]["exact_theta_minus_theta_zero_ball"])
    t, pi, imaginary = arb(HEIGHT), arb.pi(), acb(0, 1)
    theta_zero = t * ((t / (2 * pi)).log() - 1) / 2 - pi / 8
    exact_phase_factor = (imaginary * theta_correction).exp()
    gamma_outer = arb(0)
    classical_outer = arb(0)
    rows: list[dict[str, Any]] = []
    for mode_int in range(LOWER_START, UPPER_END + 1):
        target_indicator = 1 if mode_int <= TARGET_END else 0
        coefficient = 1 - target_indicator
        carrier = (imaginary * (theta_zero - t * arb(mode_int).log())).exp()
        gamma_mode = 2 * (gamma * carrier).real / arb(mode_int).sqrt()
        classical_mode = 2 * (exact_phase_factor * carrier).real / arb(mode_int).sqrt()
        if coefficient:
            gamma_outer += gamma_mode
            classical_outer += classical_mode
            rows.append(
                {
                    "mode": mode_int,
                    "target_indicator": target_indicator,
                    "positive_full_line_bulk_projector_coefficient": coefficient,
                    "exact_Gamma_full_line_physical_ball": gamma_mode.str(75, more=True),
                    "exact_classical_full_line_physical_ball": classical_mode.str(75, more=True),
                    "Gamma_minus_classical_physical_ball": (gamma_mode - classical_mode).str(75, more=True),
                }
            )
    projector_completed = complete_endpoint + gamma_outer
    classical_completed = complete_endpoint + classical_outer
    require(len(rows) == 42 and rows[0]["mode"] == 39_895 and rows[-1]["mode"] == 39_936, "outer jump roster drift")
    require(complete_endpoint.upper() < arb("-0.129") and complete_endpoint.lower() > arb("-0.130"), "endpoint scale guard lost")
    require(gamma_outer.lower() > arb("0.092") and gamma_outer.upper() < arb("0.093"), "Gamma jump scale guard lost")
    require(projector_completed.upper() < arb("-0.036") and projector_completed.lower() > arb("-0.037"), "projector-completed scale guard lost")
    return {
        "height": HEIGHT,
        "endpoint": ENDPOINT,
        "endpoint_mode_range": [LOWER_START, UPPER_END],
        "endpoint_mode_count": UPPER_END - LOWER_START + 1,
        "outer_Gamma_mode_range": [TARGET_END + 1, UPPER_END],
        "outer_Gamma_mode_count": len(rows),
        "precision_decimal_digits": precision,
        "components": {
            "complete_local_exact_affine_ball": local_affine.str(75, more=True),
            "complete_local_exact_minus_affine_ball": local_nonlinear.str(75, more=True),
            "complete_exact_exterior_ball": exact_exterior.str(75, more=True),
            "complete_exact_A_endpoint_carrier_ball": complete_endpoint.str(75, more=True),
            "outer_positive_full_line_Gamma_bulk_ball": gamma_outer.str(75, more=True),
            "outer_positive_full_line_classical_bulk_ball": classical_outer.str(75, more=True),
            "outer_Gamma_minus_classical_bulk_ball": (gamma_outer - classical_outer).str(75, more=True),
            "projector_completed_exact_A_transition_ball": projector_completed.str(75, more=True),
            "classical_projector_completed_exact_A_transition_ball": classical_completed.str(75, more=True),
        },
        "outer_Gamma_rows": rows,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    v = c["components"]
    return f"""# Complete exact A endpoint and projector transition

Date: 2026-08-14

Status: complete saved-height A endpoint carrier and adjacent positive
full-line projector jump certified; not a complete `R_Dir` bound

The endpoint coefficient does not jump at the target edge.  The exact finite
paired-projector identity is

```text
I_m+I_-m-chi_T(m)P_m
 =(1-chi_T(m))P_m+P_-m+Q_m+Q_-m.                    (AP1)
```

Therefore `Q_m+Q_-m` has coefficient `+1` throughout modes `39853..39936`,
while the positive full-line Gamma bulk has coefficient zero through mode
`39894` and one on modes `39895..39936`.

The complete exact endpoint carrier is assembled before norms:

```text
local exact affine       = {v['complete_local_exact_affine_ball']},
local exact-minus-affine = {v['complete_local_exact_minus_affine_ball']},
full exact exterior      = {v['complete_exact_exterior_ball']},
complete exact A endpoint= {v['complete_exact_A_endpoint_carrier_ball']}. (AP2)
```

For the occupancy jump, the exact global-Morse Gamma factor `B(t)` gives

```text
G_m=2 Re[B(t) exp(i(theta_0-t log m))]/sqrt(m),
sum_(39895<=m<=39936) G_m
   ={v['outer_positive_full_line_Gamma_bulk_ball']},
Gamma-minus-classical sum={v['outer_Gamma_minus_classical_bulk_ball']}. (AP3)
```

Combining (AP2) and (AP3) yields the projector-completed transition

```text
A_transition={v['projector_completed_exact_A_transition_ball']}.        (AP4)
```

The endpoint and full-line terms remain separately recorded because they
have different projector ownership.  Their signed combination in (AP4) is
the admissible A-transition object for later `R_Dir` assembly.

Pi provenance: every `pi` in (AP1)--(AP4) comes from equation (9), the exact
bi-Morse/Fresnel maps, odd-endpoint parity, the Gamma/Gaussian bulk identity,
and the paper normalization.  No fitted constant is introduced.

Proof boundary: the exact A endpoint carrier and its adjacent 42-mode
positive full-line projector jump at `t=10^10` only.  The B endpoint, zero
mode, negative and remote-positive paired currents, complete `R_Dir`,
`Q_K-T`, `T_upper`, any all-height theorem, `Lambda<=0`, PF-infinity, RH,
and any prize-level conclusion remain unproved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["local_nonlinear"]["decision"]["compact_exact_minus_affine_transformed_amplitude_certified"] is True, "nonlinear dependency drift")
    require(dependencies["exact_exterior"]["decision"]["full_exact_exterior_value_after_endpoint_replacement_certified"] is True, "exterior dependency drift")
    require(dependencies["gamma_bulk"]["decision"]["ordinary_full_line_bulk_closed_at_saved_height"] is True, "Gamma dependency drift")
    require(dependencies["orientation"]["decision"]["paired_endpoint_coefficient_continuous_across_target_edge"] is True, "endpoint orientation drift")
    require(dependencies["orientation"]["decision"]["positive_full_line_bulk_projector_jump_identified"] is True, "projector jump drift")
    require(dependencies["projector"]["decision"]["target_full_line_carrier_subtracted_before_norms"] is True, "projector dependency drift")

    artifact = {
        "kind": STEM,
        "status": "complete_exact_A_endpoint_and_positive_full_line_projector_transition_certified_R_Dir_open",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "certificate": assembly_certificate(dependencies, PRECISION),
        "decision": {
            "complete_exact_A_endpoint_carrier_certified": True,
            "positive_full_line_Gamma_jump_certified": True,
            "projector_completed_exact_A_transition_certified": True,
            "endpoint_and_bulk_projector_ownership_kept_separate": True,
            "complete_R_Dir_proved": False,
            "complete_Q_K_minus_T_proved": False,
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
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "Insert the projector-completed exact A transition into the cancellation-preserving R_Dir ledger with the independently certified B endpoint, zero mode, negative modes, remote positive modes, and outer blocks. Audit ownership before taking any modulus.",
        "proof_boundary": "Complete exact A endpoint carrier and adjacent 42-mode positive full-line projector jump at t=10^10 only. No complete R_Dir or Q_K-T estimate, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified complete exact A endpoint and projector transition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
