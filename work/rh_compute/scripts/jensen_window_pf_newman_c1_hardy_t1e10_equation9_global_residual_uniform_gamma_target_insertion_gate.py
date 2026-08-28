#!/usr/bin/env python3
"""Insert the exact common Gamma carrier into the global endpoint residual."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "finite_rejoin": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_finite_block_dirichlet_rejoin_gate.json",
    "Gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "paired_target": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "double_Weber": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate.json",
    "coverage": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
}

TARGET_START = 622
ORDINARY_END = 39_694
FOLD_START = 39_695
TARGET_END = 39_894
TARGET_COUNT = TARGET_END - TARGET_START + 1
GAMMA_PHYSICAL_BUDGET = Decimal("5e-20")
WORKING_DIRICHLET_TARGET = Decimal("1.405792e-4")
NEGATIVE_DIRICHLET_TARGET = Decimal("1.319792e-4")


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
    half, dirichlet, target_fourier, bulk, endpoint = sp.symbols(
        "H D I_T B_T Q_T"
    )
    rejoined = half + (dirichlet - target_fourier) + endpoint
    require(
        sp.expand(rejoined.subs(target_fourier, bulk + endpoint) - (half + dirichlet - bulk)) == 0,
        "uniform bulk insertion identity failed",
    )

    source, gamma, classical, window, outer, residual = sp.symbols(
        "Q_K G T E_B E_outer R_Dir"
    )
    gamma_residual = source - gamma
    exact_residual = source - classical
    require(
        sp.expand(exact_residual - (gamma_residual + gamma - classical)) == 0,
        "Gamma-to-classical correction identity failed",
    )
    require(
        sp.expand(gamma_residual - (window + outer + residual)).subs(
            source, gamma + window + outer + residual
        ) == 0,
        "window/outer/Dirichlet allocation failed",
    )

    beta, phase, carrier = sp.symbols("Beta P z")
    gamma_mode = beta * carrier
    classical_mode = phase * carrier
    require(
        sp.expand(gamma_mode - classical_mode - (beta - phase) * carrier) == 0,
        "modewise common-factor identity failed",
    )

    ordinary_count = ORDINARY_END - TARGET_START + 1
    fold_count = TARGET_END - FOLD_START + 1
    require(ordinary_count == 39_073, "ordinary target count drift")
    require(fold_count == 200, "fold target count drift")
    require(ordinary_count + fold_count == TARGET_COUNT == 39_273, "target partition drift")

    return {
        "weighted_source": "S_(M,epsilon)=H_x+integral_0^L f_x(u)D_(M,epsilon)(u)du",
        "weighted_target_fourier": "I_(T,epsilon)=integral_0^L f_x(u)G_(T,epsilon)(u)du=sum_(m=622)^39894 w_m I_m",
        "mode_decomposition": "I_m=P_bulk,m+Q_m, Q_m=P_A,m+P_B,m",
        "uniform_insertion": "H+integral f(D-G_T)+sum_T w_m Q_m=S_(M,epsilon)-sum_T w_m P_bulk,m",
        "Gamma_residual": "R_KGamma=Q_K-G=E_Btr,win+E_outer+R_Dir",
        "exact_target_residual": "Q_K-T=R_KGamma+(G-T)",
        "Gamma_mode": "g_m=B(t)exp(i(theta_0-t log m))/sqrt(m)",
        "classical_mode": "tau_m=exp(i(theta-t log m))/sqrt(m)",
        "mode_defect": "g_m-tau_m=[B(t)-exp(i(theta-theta_0))]exp(i(theta_0-t log m))/sqrt(m)",
    }


def correction_certificate(gamma: dict[str, Any]) -> dict[str, str]:
    ctx.dps = 100
    physical = arb(gamma["certificate"]["aggregate_physical_absolute_bound_ball"])
    budget = arb(str(GAMMA_PHYSICAL_BUDGET))
    require(physical.upper() < budget, "Gamma-to-classical physical correction exceeds 5e-20")

    getcontext().prec = 50
    working = WORKING_DIRICHLET_TARGET - GAMMA_PHYSICAL_BUDGET
    negative = NEGATIVE_DIRICHLET_TARGET - GAMMA_PHYSICAL_BUDGET
    require(working == Decimal("0.00014057919999999995"), "working corrected target drift")
    require(negative == Decimal("0.00013197919999999995"), "negative corrected target drift")
    return {
        "Gamma_to_classical_physical_absolute_bound_ball": gamma["certificate"]["aggregate_physical_absolute_bound_ball"],
        "safe_Gamma_correction_budget": str(GAMMA_PHYSICAL_BUDGET),
        "working_R_Dir_sufficient_target": str(working),
        "negative_R_Dir_sufficient_target": str(negative),
        "target_identity": "Q_K-T=E_Btr,win+E_outer+R_Dir+E_Gamma, |E_Gamma|<5e-20",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["correction_certificate"]
    return f"""# Uniform Gamma-target insertion in the global residual

Date: 2026-08-13

Status: exact common-carrier identification and corrected target arithmetic;
not a proof of the compressed endpoint-defect bound

For the finite Abel weight `w_m=exp(-pi*epsilon*m^2)`, write

```text
S_(M,epsilon)=H_x+integral_0^L f_x(u)D_(M,epsilon)(u)du,
I_(T,epsilon)=integral_0^L f_x(u)G_(T,epsilon)(u)du
             =sum_(m=622)^39894 w_m I_m.             (GI1)
```

The exact current decomposition is

```text
I_m=P_bulk,m+Q_m,       Q_m=P_A,m+P_B,m.             (GI2)
```

Substituting (GI2) before any norm gives

```text
H_x+integral f_x(D_(M,epsilon)-G_(T,epsilon))
   +sum_(m=622)^39894 w_m Q_m
 =S_(M,epsilon)-sum_(m=622)^39894 w_m P_bulk,m.      (GI3)
```

Thus the target kernel in Section 11.423 already inserts one exact full-line
Gamma bulk carrier for every target mode.  The ordinary/fold partition

```text
622..39694      39073 modes,
39695..39894      200 modes                           (GI4)
```

does not create a second bulk carrier.  Its distinction begins only in the
finite endpoint defect `Q_m`.  In particular the exceptional A-face term is
retained in `Q_m`; it is not an exception to (GI3).

After the common Abel limit and physical projection, refine the notation of
Section 11.423 by writing

```text
R_KGamma=Q_K-G=E_Btr,win+E_outer+R_Dir,              (GI5)
Q_K-T=R_KGamma+(G-T).                                (GI6)
```

Here the exact common mode factors are

```text
g_m=B(t) exp(i(theta_0-t log m))/sqrt(m),
tau_m=exp(i(theta-t log m))/sqrt(m),
g_m-tau_m=[B(t)-exp(i(theta-theta_0))]
             exp(i(theta_0-t log m))/sqrt(m).        (GI7)
```

The saved Arb certificate over all 39273 target modes gives

```text
|G-T|_physical <= {c['Gamma_to_classical_physical_absolute_bound_ball']}
                 < {c['safe_Gamma_correction_budget']}.                (GI8)
```

Therefore the fully explicit sufficient targets are

```text
R_Dir < {c['working_R_Dir_sufficient_target']}
        ==> Q_K-T < 8.6e-6,

R_Dir < {c['negative_R_Dir_sufficient_target']}
        ==> Q_K-T < 0,                               (GI9)
```

using the certified negative B window and the `8e-10` outer allowance.  The
`5e-20` correction is immaterial numerically but necessary for an exact
object identity.

The remaining theorem is now unambiguous: bound the single finite endpoint
defect `R_Dir`, with its A face, B complement, zero, negative, outer, and
half-current sectors still grouped.  No separate ordinary or fold bulk
identification remains.

Pi provenance: `pi` in (GI1)--(GI7) comes from the equation-(9) Gaussian
Abel weight, Fourier character, Kummer phase, and exact Riemann--Siegel phase.
No fitted or geometric constant is introduced.

Proof boundary: exact carrier insertion, object identification, and a
saved-height Gamma-to-classical correction bound only.  No quantitative
bound for `R_Dir`, complete endpoint residual theorem, all-corridor or
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["finite_rejoin"]["decision"]["finite_completed_B_block_rejoined_before_norms"] is True,
        "finite-rejoin dependency drift",
    )
    require(
        dependencies["Gamma_bulk"]["decision"]["ordinary_full_line_bulk_closed_at_saved_height"] is True,
        "Gamma-bulk dependency drift",
    )
    require(
        dependencies["paired_target"]["decision"]["target_carrier_subtracted_exactly_once"] is True,
        "paired-target dependency drift",
    )
    require(
        dependencies["double_Weber"]["decision"]["entire_finite_equation9_source_roster_reconstructed_exactly"] is True,
        "double-Weber source dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "uniform_Gamma_bulk_target_inserted_across_ordinary_and_fold_rosters_exact_classical_correction_below_5e_minus_20",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "ordinary_target_modes": [TARGET_START, ORDINARY_END],
            "ordinary_target_count": ORDINARY_END - TARGET_START + 1,
            "fold_target_modes": [FOLD_START, TARGET_END],
            "fold_target_count": TARGET_END - FOLD_START + 1,
            "total_target_count": TARGET_COUNT,
        },
        "symbolic_certificate": symbolic_certificate(),
        "correction_certificate": correction_certificate(dependencies["Gamma_bulk"]),
        "decision": {
            "Dirichlet_target_kernel_is_exact_Gamma_bulk_insertion": True,
            "ordinary_and_fold_bulk_carriers_are_one_common_object": True,
            "A_face_exception_remains_in_endpoint_defect_not_bulk_carrier": True,
            "Gamma_residual_identified_as_QK_minus_G": True,
            "exact_classical_residual_requires_Gamma_correction": True,
            "Gamma_to_classical_physical_correction_below_5e_minus_20": True,
            "compressed_R_Dir_target_proved": False,
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
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Bound the single endpoint defect R_Dir in (GI3) below 0.00014057919999999995, preserving the A face, B complement, zero, negative, finite outer, endpoint-half, and common Abel sectors. Use the double-Weber representation only after the uniform Gamma carrier has been removed; do not create separate ordinary and fold bulk budgets.",
        "proof_boundary": "Exact uniform Gamma-carrier insertion, Q_K-G versus Q_K-T object separation, and saved-height Gamma-to-classical correction below 5e-20 only. No quantitative R_Dir bound, complete endpoint residual theorem, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified uniform Gamma target insertion and exact classical correction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
