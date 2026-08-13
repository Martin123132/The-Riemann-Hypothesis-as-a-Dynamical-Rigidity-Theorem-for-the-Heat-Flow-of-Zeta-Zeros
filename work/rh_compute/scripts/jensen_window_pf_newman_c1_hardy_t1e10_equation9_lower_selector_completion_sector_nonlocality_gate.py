#!/usr/bin/env python3
"""Certify the large completed sector outside the selector fold packet."""

from __future__ import annotations

import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = REPO_ROOT / "work/rh_compute/scripts"
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(SCRIPT_ROOT))
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

packet_gate = importlib.import_module(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate"
)
from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_selector_completion_sector_nonlocality_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "completed_profile": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "common_profile": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "exact_profile": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.json",
    "source_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.json",
    "selector_increment": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate.json",
    "two_packet_barrier": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate.json",
}
PRECISION = 70


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
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Selector completion-sector nonlocality

Date: 2026-08-13

Status: rigorous completion-sector barrier for the selector fold packet;
not a fixed-state source-minus-target estimate

The exact one-cell Poisson completion and paired fold selector identity use
the same source normalization.  With

```text
P_C=Q_half^(C)-Q_half^(C+2),
D_fold=F_C-F_(C+2),
C_out=P_C-D_fold,
```

finite algebra gives

```text
P_C=D_fold+C_out.                                      (CN1)
```

The exact completed common profile satisfies

```text
P_C=kappa_C g_ex(lambda,0),
kappa_C=exp(i*beta^3)2sqrt(2).
```

The beta-minus-four endpoint profile is enclosed on the whole ordinary top
corridor with `lambda` kept as one interval.  Appending the certified
Kummer-ODE exact-profile error and applying odd-square reflection gives

```text
{c['physical_completed_source_lower_ball']} < Q_C
 < {c['physical_completed_source_upper_ball']}.        (CN2)
```

The independent two-packet theorem supplies

```text
{c['physical_fold_increment_lower_ball']} < Delta Q_fold
 < {c['physical_fold_increment_upper_ball']}.          (CN3)
```

Projecting (CN1) and subtracting interval endpoints therefore proves

```text
{c['physical_outside_completion_lower_ball']} < Q_out
 < {c['physical_outside_completion_upper_ball']},
Q_out=2Re[exp(-i*pi/8)C_out].                          (CN4)
```

In particular `Q_out>200`, while `Delta Q_fold<35`.  The fold packet is less
than `{c['fold_to_source_fraction_upper_ball']}` of the completed source
projection, and the outside completion exceeds it by a factor greater than
`{c['outside_to_fold_ratio_lower_ball']}`.

This is a nonlocality result, not a small-error result.  A fixed-B selector
jump is not closed inside the proposed fold ownership block: its endpoint
half-current, zero and other Fourier sectors carry a mandatory projected
contribution above `200`.  Consequently neither the 200-mode fold packet nor
the 39273-mode positive target roster may replace the completed current one
cell at a time.  The fixed-state endpoint residual must be estimated in the
global cancellation-preserving form of Section 11.350, or in another form
that proves the same completion cancellation explicitly.

Pi provenance: `beta^3=pi*C^2/8`, the Fourier phases `2*pi*n`, and the
half-Kummer phase `pi/8` all descend from the exact equation-(9) transform.
No fitted or geometric constant is introduced.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No bound for the global fixed-state endpoint residual, complete `Q_K-T` or
`T_upper`, all-corridor or height-uniform theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    ctx.dps = PRECISION
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in DEPENDENCIES.items()
    }
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["completed_profile"]["decision"]["beta4_all_mode_Fourier_completion_proved"] is True,
        "completed beta4 profile drift",
    )
    require(
        dependencies["common_profile"]["decision"]["exact_transformed_strip_common_profile_identity_proved"] is True,
        "exact common-profile identity drift",
    )
    require(
        dependencies["selector_increment"]["decision"]["complete_selector_difference_is_removed_source_term"] is True,
        "selector completion identity drift",
    )

    c_value = arb(packet_gate.C)
    beta = (arb.pi() * c_value * c_value / 8) ** (arb(1) / 3)
    lambda_max = arb.pi() / (16 * beta)
    lam = arb(lambda_max / 2, lambda_max / 2)
    y_max = arb.pi() * c_value / (2 * beta)
    g4_origin, _, _ = packet_gate.profile_derivatives(lam, arb(0), beta, y_max)
    exact_profile_error = arb(
        dependencies["exact_profile"]["interval_certificate"]["uniform_canonical_exact_minus_beta4_profile_bound"]
    ).upper()
    canonical_source_lower = (g4_origin.real.lower() - exact_profile_error).lower()
    canonical_source_upper = (g4_origin.real.upper() + exact_profile_error).upper()
    physical_factor = 4 * arb(2).sqrt()
    physical_source_lower = (physical_factor * canonical_source_lower).lower()
    physical_source_upper = (physical_factor * canonical_source_upper).upper()
    require(physical_source_lower > 230, "completed source projection lost lower barrier 230")

    packet_certificate = dependencies["two_packet_barrier"]["certificate"]
    fold_lower = arb(packet_certificate["physical_exact_fold_increment_lower_ball"]).lower()
    fold_upper = arb(packet_certificate["physical_exact_fold_increment_upper_ball"]).upper()
    require(fold_lower > 20 and fold_upper < 35, "two-packet enclosure drift")
    outside_lower = (physical_source_lower - fold_upper).lower()
    outside_upper = (physical_source_upper - fold_lower).upper()
    require(outside_lower > 200, "outside completion did not clear lower barrier 200")
    fraction_upper = (fold_upper / physical_source_lower).upper()
    ratio_lower = (outside_lower / fold_upper).lower()
    require(fraction_upper < arb("0.14"), "fold/source fraction exceeded 0.14")
    require(ratio_lower > 6, "outside/fold ratio lost factor six")

    certificate = {
        "height_interval": "t*-pi/16<=t<=t*",
        "beta_ball": beta.str(PRECISION, more=True),
        "lambda_interval": "0<=lambda<=pi/(16beta)",
        "beta4_completed_profile_origin_real_ball": g4_origin.real.str(PRECISION, more=True),
        "uniform_exact_profile_error_ball": exact_profile_error.str(PRECISION, more=True),
        "canonical_completed_source_real_lower_ball": canonical_source_lower.str(PRECISION, more=True),
        "canonical_completed_source_real_upper_ball": canonical_source_upper.str(PRECISION, more=True),
        "physical_completed_source_lower_ball": physical_source_lower.str(PRECISION, more=True),
        "physical_completed_source_upper_ball": physical_source_upper.str(PRECISION, more=True),
        "physical_fold_increment_lower_ball": fold_lower.str(PRECISION, more=True),
        "physical_fold_increment_upper_ball": fold_upper.str(PRECISION, more=True),
        "physical_outside_completion_lower_ball": outside_lower.str(PRECISION, more=True),
        "physical_outside_completion_upper_ball": outside_upper.str(PRECISION, more=True),
        "fold_to_source_fraction_upper_ball": fraction_upper.str(PRECISION, more=True),
        "outside_to_fold_ratio_lower_ball": ratio_lower.str(PRECISION, more=True),
    }
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "selector_completion_sector_nonlocality_interval_certified",
        "passed": True,
        "exact_identity": {
            "completed_selector_source": "P_C=kappa_C*g_ex(lambda,0)",
            "fold_selector_increment": "D_fold=F_C-F_(C+2)",
            "outside_completion": "C_out=P_C-D_fold",
            "completion_balance": "P_C=D_fold+C_out",
            "physical_projection": "Q=2Re[exp(-i*pi/8)*complex_half_domain_quantity]",
        },
        "certificate": certificate,
        "decision": {
            "completed_selector_source_profile_identity_used": True,
            "fold_and_outside_completion_balance_proved": True,
            "outside_completion_projection_uniformly_above_200": True,
            "fold_selector_packet_is_cancellation_closed": False,
            "one_cell_positive_target_roster_is_global_fixed_state_residual": False,
            "fixed_state_endpoint_residual_bounded": False,
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
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Return to the fixed-state global endpoint residual R_end=lim[H+integral f C_M+sum_target(P_A+P_B)]. "
            "Derive its lower-A characteristic splice without assigning a local target block to one selector cell."
        ),
        "proof_boundary": (
            "Uniform completed-source, fold-packet, and outside-completion projection enclosures on one ordinary top "
            "corridor only. No global fixed-state endpoint residual, complete Q_K-T or T_upper, all-corridor or "
            "height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print(
        f"certified outside selector completion > {outside_lower}; local fold ownership is not cancellation-closed",
        flush=True,
    )


if __name__ == "__main__":
    main()
