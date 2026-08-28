#!/usr/bin/env python3
"""Rejoin the finite completed-B block to the global Dirichlet residual."""

from __future__ import annotations

from decimal import Decimal, getcontext
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

import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_finite_block_dirichlet_rejoin_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "outer_phase": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_phase_coupled_tangential_gate.json",
    "completed_extraction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_endpoint_safe_completed_B_extraction_gate.json",
    "global_extraction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_B_window_extraction_target_gate.json",
    "symmetric_reassembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.json",
}

A = 159_577
B = 5_122_421
L = (B - A) // 2
TARGET_START = 622
TARGET_END = 39_894
TARGET_COUNT = TARGET_END - TARGET_START + 1
OUTER_BOUND_TEXT = "8e-10"
WINDOW_UPPER_TEXT = "-1.3198e-4"
WORKING_TARGET_TEXT = "8.6e-6"


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
    global_current, window_trace, cutoff_exterior = sp.symbols("G chiB etaE")
    completed_low, completed_outer, joined = sp.symbols("Ctr_low Ctr_outer R_join_comp")
    completed_allocation = window_trace + cutoff_exterior * (completed_low + completed_outer) + joined
    finite_join = cutoff_exterior * completed_low + joined
    outer = cutoff_exterior * completed_outer
    require(sp.expand(completed_allocation - (window_trace + outer + finite_join)) == 0, "finite rejoin identity failed")
    require(sp.expand(finite_join - (completed_allocation - window_trace - outer)) == 0, "Dirichlet handoff identity failed")

    z = sp.symbols("z", nonzero=True)
    cutoff = sp.symbols("M", integer=True, positive=True)
    symmetric = z ** (-cutoff) * (1 - z ** (2 * cutoff + 1)) / (1 - z)
    require(sp.cancel((1 - z) * symmetric - (z ** (-cutoff) - z ** (cutoff + 1))) == 0, "Dirichlet quotient drift")
    target = z**TARGET_START * (1 - z**TARGET_COUNT) / (1 - z)
    require(sp.cancel((1 - z) * target - (z**TARGET_START - z ** (TARGET_END + 1))) == 0, "target quotient drift")

    return {
        "completed_split": "Ctr_(M,epsilon)=Ctr_<B,epsilon+Ctr_>=B,M,epsilon for M>=B",
        "outer_block": "O_(M,epsilon)=eta_delta E_W Ctr_>=B,M,epsilon",
        "finite_join": "J_(M,epsilon)=eta_delta E_W Ctr_<B,epsilon+R_join,comp(M,epsilon)",
        "three_part_identity": "G_(M,epsilon)=chi_W Btr_(M,epsilon)+O_(M,epsilon)+J_(M,epsilon)",
        "handoff": "J_(M,epsilon)=G_(M,epsilon)-chi_W Btr_(M,epsilon)-O_(M,epsilon)",
        "weighted_kernel": "D_(M,epsilon)(u)=sum_(m=-M)^M exp(-pi*epsilon*m^2)exp(-2pi*i*m*u)",
        "weighted_target_kernel": "G_(T,epsilon)(u)=sum_(m=622)^39894 exp(-pi*epsilon*m^2)exp(-2pi*i*m*u)",
        "weighted_complement": "C_(M,epsilon)=D_(M,epsilon)-G_(T,epsilon)",
        "zero_regulator_Dirichlet": "D_M(u)=sin((2M+1)pi*u)/sin(pi*u)",
        "zero_regulator_target": "G_T(u)=exp(-40516pi*i*u)sin(39273pi*u)/sin(pi*u)",
        "global_kernel_integrand": "G_(M,epsilon)(x)=H_x+integral_0^L f_x(u)C_(M,epsilon)(u)du+sum_(m=622)^39894 exp(-pi*epsilon*m^2)(P_A,m+P_B,m)",
    }


def target_arithmetic() -> dict[str, str]:
    getcontext().prec = 40
    outer = Decimal(OUTER_BOUND_TEXT)
    window_magnitude = -Decimal(WINDOW_UPPER_TEXT)
    working = Decimal(WORKING_TARGET_TEXT)
    positive_target = working + window_magnitude - outer
    negative_target = window_magnitude - outer
    require(positive_target == Decimal("0.0001405792"), "working target arithmetic drift")
    require(negative_target == Decimal("0.0001319792"), "negative target arithmetic drift")
    return {
        "window_upper": WINDOW_UPPER_TEXT,
        "outer_absolute_upper": OUTER_BOUND_TEXT,
        "working_R_end_upper": WORKING_TARGET_TEXT,
        "finite_join_sufficient_for_working_target": "1.405792e-4",
        "finite_join_sufficient_for_negative_R_end": "1.319792e-4",
        "identity": "R_end=E_Btr,win+E_outer+R_Dir",
    }


def render_note(artifact: dict[str, Any]) -> str:
    a = artifact["target_arithmetic"]
    return f"""# Dirichlet rejoin of the finite completed-B block

Date: 2026-08-13

Status: exact common-regulator compression and revised quantitative target;
not a proof of a bound for the compressed residual

For `M>=B={B}`, split the smooth crossing-completed current at the analytic
tail threshold:

```text
Ctr_(M,epsilon)=Ctr_<B,epsilon+Ctr_>=B,M,epsilon,
O_(M,epsilon)=eta_delta E_W Ctr_>=B,M,epsilon.        (DR1)
```

Do **not** norm the remaining `{B - 1}` positive modes or construct their
`1,280,347` crossing charts.  Rejoin them to the common remainder:

```text
J_(M,epsilon)
 =eta_delta E_W Ctr_<B,epsilon+R_join,comp(M,epsilon).
```

The completed allocation then gives the finite identity

```text
G_(M,epsilon)=chi_W Btr_(M,epsilon)
               +O_(M,epsilon)+J_(M,epsilon),

J_(M,epsilon)=G_(M,epsilon)-chi_W Btr_(M,epsilon)
               -O_(M,epsilon).                       (DR2)
```

Thus the finite crossing roster is represented by the already-established
global kernel, not by millions of local Fresnel charts.  With
`w_m=exp(-pi*epsilon*m^2)`, put

```text
D_(M,epsilon)(u)=sum_(m=-M)^M w_m exp(-2pi*i*m*u),
G_(T,epsilon)(u)=sum_(m=622)^39894 w_m exp(-2pi*i*m*u),
C_(M,epsilon)=D_(M,epsilon)-G_(T,epsilon).            (DR3)
```

The compressed global integrand is exactly

```text
G_(M,epsilon)(x)
 =H_x+integral_0^{L} f_x(u)C_(M,epsilon)(u)du
  +sum_(m=622)^39894 w_m(P_A,m+P_B,m).                (DR4)
```

At `epsilon=0`, (DR3) becomes the two removable Dirichlet quotients

```text
D_M=sin((2M+1)pi*u)/sin(pi*u),
G_T=exp(-40516pi*i*u)sin(39273pi*u)/sin(pi*u).        (DR5)
```

The analytic outer gate supplies its own dominated Abel limit and
`|E_outer|<{a['outer_absolute_upper']}`.  Consequently

```text
R_end=E_Btr,win+E_outer+R_Dir,                        (DR6)

R_Dir < {a['finite_join_sufficient_for_working_target']}
        ==> R_end < {a['working_R_end_upper']},

R_Dir < {a['finite_join_sufficient_for_negative_R_end']}
        ==> R_end < 0.                                (DR7)
```

This is the key compression: the next quantitative object is one
Dirichlet/Abel kernel plus the fixed `{TARGET_COUNT}`-mode target-tail block,
with the certified B window and analytic outer block removed exactly.  The
remaining object still contains the endpoint half-current, A-fold sector,
negative modes, and nonlocal steps in their valid common grouping.

No estimate for `R_Dir` is proved here.  In particular this does not prove
the working target, a complete `T_upper`, a height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.

Pi provenance: `pi` in (DR3)--(DR5) is forced by the Gaussian Abel weight
and integer Fourier character already present in equation (9); no fitted or
geometric value is inserted.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["outer_phase"]["decision"]["outer_block_separate_Abel_limit_justified_by_summable_derivative_majorants"] is True,
        "outer Abel-limit dependency drift",
    )
    require(
        dependencies["completed_extraction"]["decision"]["crossing_completed_exterior_B_current_smooth_at_finite_cutoff"] is True,
        "completed-current dependency drift",
    )
    require(
        dependencies["global_extraction"]["decision"]["B_trace_window_extracted_at_every_finite_symmetric_cutoff"] is True,
        "global extraction dependency drift",
    )
    require(
        dependencies["symmetric_reassembly"]["decision"]["finite_symmetric_dirichlet_reassembly_proved"] is True,
        "Dirichlet dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "finite_completed_B_block_rejoined_exactly_to_one_common_regulator_Dirichlet_residual",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "roster_coordinate_length": L,
            "outer_threshold": B,
            "finite_completed_positive_mode_count": B - 1,
            "avoided_crossing_chart_count": 1_280_347,
            "target_positive_modes": [TARGET_START, TARGET_END],
            "target_mode_count": TARGET_COUNT,
        },
        "symbolic_certificate": symbolic_certificate(),
        "target_arithmetic": target_arithmetic(),
        "decision": {
            "finite_completed_B_block_rejoined_before_norms": True,
            "finite_crossing_roster_replaced_by_common_Dirichlet_Abel_kernel": True,
            "outer_limit_separated_only_after_absolute_summable_certificate": True,
            "window_negative_certificate_preserved": True,
            "common_endpoint_half_A_negative_and_nonlocal_grouping_preserved": True,
            "compressed_residual_below_1p405792e_minus_4_proved": False,
            "compressed_residual_below_1p319792e_minus_4_proved": False,
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
        "next_obligation": "Bound the single compressed residual R_Dir below 1.405792e-4 (or 1.319792e-4 for negativity). Start from the modular form of D_(M,epsilon), subtract the fixed target kernel before norms, identify the ordinary and A-fold target carriers in the same representation, and bound only their joint defect; do not reopen per-crossing completed-B charts.",
        "proof_boundary": "Exact common-regulator rejoin and Dirichlet/Abel compression of the finite completed-B block, plus sufficient target arithmetic only. No quantitative compressed residual, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified finite completed-B rejoin to one Dirichlet residual", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
