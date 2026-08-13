#!/usr/bin/env python3
"""Certify the weighted opposite beta^-4 branch by whole-kernel L1 bounds."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
    "opposite_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_opposite_branch_roster_reflection_pairing_gate.json",
    "weighted_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.json",
    "completed_cutoffs": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_completed_cutoff_family_uniform_dirichlet_gate.json",
}

L = 39_696
Q = 39_894
U = 40_094
PANELS = 16_384
PRECISION = 80


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


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def reconstruct_weights(pairing: dict[str, Any]) -> dict[int, acb]:
    certificate = pairing["leading_multiplier_certificate"]
    weights: dict[int, acb] = {}
    for row in certificate["pair_rows"]:
        coherent = parse_complex(row["coherent_coefficient_ball"])
        leakage = parse_complex(row["conjugacy_leakage_coefficient_ball"])
        weights[int(row["minus_mode"])] = coherent + leakage
        weights[int(row["plus_mode"])] = (coherent - leakage).conjugate()
    for row in certificate["edge_rows"]:
        weights[int(row["mode"])] = parse_complex(row["leading_defect_ball"])
    require(set(weights) == set(range(L, Q)) | set(range(Q + 1, U + 1)), "weight roster drift")
    return weights


def polynomial_l1(coefficients: list[acb]) -> dict[str, Any]:
    pi = arb.pi()
    derivative = 2 * pi * sum((arb(index) * abs(value).upper() for index, value in enumerate(coefficients)), arb(0))
    midpoint_integral = arb(0)
    i = acb(0, 1)
    for panel in range(PANELS):
        s = arb(2 * panel + 1) / (2 * PANELS)
        z = (-2 * pi * i * s).exp()
        # Direct powers prevent rectangular-ball inflation from repeated rotations.
        value = sum((coefficient * z**index for index, coefficient in enumerate(coefficients)), acb(0))
        midpoint_integral += abs(value).upper()
    midpoint_integral /= PANELS
    lipschitz_padding = derivative / (2 * PANELS)
    total = midpoint_integral + lipschitz_padding
    return {
        "coefficient_count": len(coefficients),
        "panels": PANELS,
        "midpoint_absolute_integral_ball": midpoint_integral.str(PRECISION, more=True),
        "derivative_supremum_ball": derivative.str(PRECISION, more=True),
        "panel_Lipschitz_padding_ball": lipschitz_padding.str(PRECISION, more=True),
        "L1_bound_ball": total.str(PRECISION, more=True),
    }


def certificate(pairing: dict[str, Any], cutoffs: dict[str, Any]) -> dict[str, Any]:
    weights = reconstruct_weights(pairing)
    minus = polynomial_l1([weights[mode] for mode in range(L, Q)])
    plus = polynomial_l1([weights[mode] for mode in range(Q + 1, U + 1)])
    minus_l1 = arb(minus["L1_bound_ball"]).upper()
    plus_l1 = arb(plus["L1_bound_ball"]).upper()
    total_l1 = minus_l1 + plus_l1

    completed = cutoffs["certificate"]
    y_max = arb(completed["selector_y_max_ball"]).upper()
    profile = arb(completed["g4_profile_supremum_bound_ball"]).upper()
    branch_amplitude = profile / (2 * y_max)
    opposite_bound = y_max * branch_amplitude * total_l1
    full_bound = arb(completed["weighted_completed_roster_bound_ball"]).upper()
    selected_bound = full_bound + opposite_bound
    require(minus_l1 < arb("1.6e-5"), "negative weighted kernel L1 exceeded 1.6e-5")
    require(plus_l1 < arb("2.4e-5"), "positive weighted kernel L1 exceeded 2.4e-5")
    require(total_l1 < arb("4.0e-5"), "two-branch weighted kernel L1 exceeded 4.0e-5")
    require(opposite_bound < arb("0.00166"), "weighted opposite branch exceeded 0.00166")
    require(selected_bound < arb("0.02324"), "selected weighted branch exceeded 0.02324")
    return {
        "height_interval": "t*-pi/16<=t<=t*",
        "minus_branch_kernel": minus,
        "plus_branch_kernel": plus,
        "two_branch_kernel_L1_bound_ball": total_l1.str(PRECISION, more=True),
        "full_strip_single_branch_amplitude_bound_ball": branch_amplitude.str(PRECISION, more=True),
        "weighted_opposite_branch_bound_ball": opposite_bound.str(PRECISION, more=True),
        "auxiliary_weighted_full_roster_bound_ball": full_bound.str(PRECISION, more=True),
        "selected_weighted_branch_bound_ball": selected_bound.str(PRECISION, more=True),
        "selected_decomposition": "W_sel=W_full-W_opp",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    minus = c["minus_branch_kernel"]
    plus = c["plus_branch_kernel"]
    return f"""# Weighted opposite-branch Dirichlet bound

Date: 2026-08-13

Status: rigorous ordinary-top-corridor leading-defect selected/full/opposite
branch decomposition; exact finite-integral remainder remains open

For each beta-minus-four mode write

```text
G4_m=G_sel,m+G_opp,m.
```

With the extended weights from Section 11.393,

```text
W_sel=W_full-W_opp,
W_opp=sum_(m=L)^U w_mG_opp,m.                       (OD1)
```

The opposite branch is `C_+` on modes `39696..39893` and `C_-` on modes
`39895..40094`; event zero has weight zero.  With `s=y/Y` and `hY=2pi`,
factoring a unit-modulus initial frequency leaves two ordinary weighted
Fourier polynomials

```text
P_-(s)=sum_(j=0)^197 c^-_j e^(-2pi i j s),
P_+(s)=sum_(j=0)^199 c^+_j e^(-2pi i j s).          (OD2)
```

For a polynomial `P`, midpoint enclosure on `N_p={PANELS}` equal panels and

```text
sup_s |P'(s)|<=2pi sum_j j|c_j|                     (OD3)
```

give a rigorous whole-kernel `L1` bound.  The saved coefficient balls are
uniform over the complete top corridor.  The independently checked results
are

```text
Integral_0^1|P_-(s)|ds
 <{minus['L1_bound_ball']},
Integral_0^1|P_+(s)|ds
 <{plus['L1_bound_ball']},
sum <{c['two_branch_kernel_L1_bound_ball']}<4.0e-5.                  (OD4)
```

The full-strip branch envelope from Section 11.394 gives

```text
sup_(lambda,y)|C_+|=sup_(lambda,y)|C_-|
 <{c['full_strip_single_branch_amplitude_bound_ball']}.              (OD5)
```

Applying (OD4)--(OD5) once to each whole branch proves

```text
|W_opp|<{c['weighted_opposite_branch_bound_ball']}<0.00166.          (OD6)
```

Together with the auxiliary completed full-roster bound,

```text
|W_sel|<=|W_full|+|W_opp|
 <{c['selected_weighted_branch_bound_ball']}<0.02324.                (OD7)
```

No modewise finite-integral norm is used.  The estimate is conservative and
does not exploit cancellation between `W_full` and `W_opp`; its role is to
make the completion-to-selected-branch passage mathematically valid.  It
still concerns the certified leading stationary-defect weights only, not the
exact finite-integral minus beta-minus-four amplitude remainder.

Pi provenance: the Fourier phases use the exact period `hY=2pi`, while the
branch amplitude and height interval inherit `beta^3=pi C^2/8`.  No new
geometric normalization is introduced.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No exact finite-integral amplitude remainder, completed branch source/initial
data splice, complete `Q_K-T` or `T_upper`, all-corridor theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")
    ctx.dps = PRECISION
    ctx.threads = 1
    cert = certificate(dependencies["event_pairing"], dependencies["completed_cutoffs"])
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "weighted_opposite_and_selected_beta4_branches_bounded",
        "passed": True,
        "certificate": cert,
        "exact_algebra": {
            "branch_split": "G4_m=G_sel,m+G_opp,m",
            "weighted_split": "W_sel=W_full-W_opp",
            "negative_opposite_branch": "C_+ on modes 39696..39893",
            "positive_opposite_branch": "C_- on modes 39895..40094",
            "whole_kernel_norms_only": True,
        },
        "decision": {
            "weighted_opposite_branch_bound_proved": True,
            "selected_branch_leading_defect_weighted_bound_proved": True,
            "exact_finite_integral_amplitude_remainder_proved": False,
            "completed_branch_initial_data_source_splice_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "panels_per_branch": PANELS,
        },
        "next_action": (
            "Sharpen W_full and W_opp jointly, or derive the exact finite-integral-minus-beta^-4 operator remainder "
            "in a completion-preserving form. The present selected bound is valid but too coarse for the local headroom."
        ),
        "proof_boundary": (
            "Uniform leading-defect weighted opposite and selected beta^-4 branch bounds on the ordinary top corridor only. "
            "No exact finite-integral amplitude remainder, completed branch source/initial-data splice, complete Q_K-T or "
            "T_upper, all-corridor theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified weighted opposite-branch Dirichlet bound", flush=True)


if __name__ == "__main__":
    main()
