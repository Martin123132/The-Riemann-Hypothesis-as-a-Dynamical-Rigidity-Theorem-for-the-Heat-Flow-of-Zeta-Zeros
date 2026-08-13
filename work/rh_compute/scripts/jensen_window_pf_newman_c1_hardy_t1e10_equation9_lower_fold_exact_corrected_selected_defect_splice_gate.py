#!/usr/bin/env python3
"""Splice the exact finite-t correction into the selected leading defect."""

from __future__ import annotations

from fractions import Fraction
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

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "selected_beta4": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate.json",
    "common_profile_operator": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "exact_profile_closure": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.json",
    "branch_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
}

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


def exact_finite_roster_check() -> dict[str, str]:
    modes = range(-7, 9)
    weights = {m: Fraction(5 * m - 2, 11 * abs(m) + 17) for m in modes}
    selected = {m: Fraction(3 * m + 1, 13 * abs(m) + 19) for m in modes}
    opposite = {m: Fraction(7 - 2 * m, 17 * abs(m) + 23) for m in modes}
    correction = {m: Fraction(4 * m - 3, 19 * abs(m) + 29) for m in modes}
    beta4 = {m: selected[m] + opposite[m] for m in modes}
    exact = {m: beta4[m] + correction[m] for m in modes}
    left = sum((weights[m] * (exact[m] - opposite[m]) for m in modes), Fraction())
    right = sum((weights[m] * selected[m] for m in modes), Fraction()) + sum(
        (weights[m] * correction[m] for m in modes), Fraction()
    )
    require(left == right, "corrected selected-roster splice failed")
    return {
        "beta4_split": "G4_m=G_sel,m+G_opp,m",
        "exact_correction": "Delta G_m=G_ex,m-G4_m",
        "corrected_selected_coefficient": "G_ex,m-G_opp,m=G_sel,m+Delta G_m",
        "weighted_splice": "W_corr=sum_m w_m(G_ex,m-G_opp,m)=W_sel,4+Delta W_finite-t",
        "independent_exact_test": "16 unrelated rational modes",
    }


def numerical_certificate(dependencies: dict[str, Any]) -> dict[str, str]:
    selected = dependencies["selected_beta4"]["certificate"]
    profile = dependencies["exact_profile_closure"]["interval_certificate"]
    selected_bound = arb(selected["uniform_selected_weighted_bound_ball"]).upper()
    finite_t_bound = arb(profile["uniform_weighted_finite_t_profile_correction_bound"]).upper()
    corrected_modulus = (selected_bound + finite_t_bound).upper()

    node_real_upper = max(
        arb(row["value_ball"]["real_ball"]).upper() for row in selected["node_rows"]
    )
    quadrature_error = arb(selected["y_midpoint_error_ball"]).upper()
    height_error = arb(selected["height_transport_error_ball"]).upper()
    corrected_real_upper = (node_real_upper + quadrature_error + height_error + finite_t_bound).upper()
    relative_correction = (finite_t_bound / selected_bound).upper()
    require(corrected_modulus < arb("4.088814e-5"), "corrected selected modulus exceeds 4.088814e-5")
    require(corrected_real_upper < -arb("2.233e-5"), "corrected selected real part is not uniformly negative")
    require(relative_correction < arb("2.935e-7"), "finite-t correction ratio exceeds 2.935e-7")
    return {
        "height_interval": "t*-pi/16<=t<=t*",
        "selected_beta4_uniform_modulus_bound": selected_bound.str(PRECISION, more=True),
        "exact_finite_t_weighted_correction_bound": finite_t_bound.str(PRECISION, more=True),
        "finite_t_to_selected_bound_ratio": relative_correction.str(PRECISION, more=True),
        "corrected_selected_uniform_modulus_bound": corrected_modulus.str(PRECISION, more=True),
        "maximum_node_real_part_upper": node_real_upper.str(PRECISION, more=True),
        "quadrature_error_bound": quadrature_error.str(PRECISION, more=True),
        "height_transport_error_bound": height_error.str(PRECISION, more=True),
        "corrected_selected_uniform_real_part_upper": corrected_real_upper.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Exact corrected selected-defect splice

Date: 2026-08-13

Status: exact cancellation-preserving splice and a rigorous signed
top-corridor bound; the final source/initial-data orientation remains open

For every mode in the augmented roster, the beta-minus-four branch split and
the exact finite-height correction are

```text
G4_m=G_sel,m+G_opp,m,
Delta G_m=G_ex,m-G4_m.                               (CS1)
```

Therefore exact algebra, without assigning an exact branch decomposition to
`G_ex,m`, gives

```text
G_ex,m-G_opp,m=G_sel,m+Delta G_m,
W_corr=sum_m w_m(G_ex,m-G_opp,m)
      =W_sel,4+Delta W_finite-t.                     (CS2)
```

The opposite beta-minus-four branch remains inside the completed remainder;
the new Kummer-ODE correction is inserted as one whole Fourier profile.  No
modewise finite-integral bound is used.

Combining the two certified estimates gives

```text
sup|W_corr|<{c['corrected_selected_uniform_modulus_bound']}<4.088814e-5. (CS3)
```

The exact correction is less than
`{c['finite_t_to_selected_bound_ratio']}` times the previous selected bound.
More importantly, all five grouped real parts are negative.  Appending the
quadrature, height-transport, and exact finite-height errors in the adverse
direction proves throughout the complete ordinary top corridor

```text
Re W_corr<{c['corrected_selected_uniform_real_part_upper']}<-2.233e-5. (CS4)
```

Thus the old scalar-headroom comparison by modulus is not the decisive test:
the corrected term has a certified sign.  Whether that sign is favorable
depends on the still-unproved common source/initial-data orientation in the
final `Q_K-T` or `T_upper` identity.  This gate does not assume that
orientation.

Pi provenance is inherited unchanged from `beta^3=pi*C^2/8`, `hY=2pi`, the
inverse Airy Fourier normalization, and the Kummer-ODE profile theorem.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No source/initial-data sign orientation, local-headroom closure, complete
`Q_K-T` or `T_upper`, all-corridor or height-uniform theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")
    require(dependencies["common_profile_operator"]["decision"]["exact_transformed_strip_common_profile_identity_proved"] is True, "common-profile decision drift")
    require(dependencies["branch_projection"]["decision"]["selected_plus_opposite_branch_equals_boundary_block"] is True, "branch split decision drift")

    ctx.dps = PRECISION
    algebra = exact_finite_roster_check()
    certificate = numerical_certificate(dependencies)
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_corrected_selected_defect_splice_and_negative_real_part_certified",
        "passed": True,
        "exact_algebra": algebra,
        "certificate": certificate,
        "decision": {
            "exact_finite_t_correction_spliced_into_selected_defect": True,
            "opposite_beta4_branch_retained_in_completed_remainder": True,
            "corrected_selected_real_part_uniformly_negative": True,
            "source_initial_data_sign_orientation_proved": False,
            "local_normalized_headroom_closed": False,
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
        },
        "next_action": (
            "Trace the corrected selected object through the exact lower-A source and initial-data coordinates into the common "
            "Q_K-T or T_upper scalar, fixing its real-part sign and normalization before any headroom comparison."
        ),
        "proof_boundary": (
            "An exact corrected-selected algebraic splice, modulus bound, and negative real-part theorem on one top corridor only. "
            "No source/initial-data sign orientation, local-headroom closure, complete Q_K-T or T_upper, all-corridor or "
            "height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified exact corrected selected defect: Re W_corr < -2.233e-5", flush=True)


if __name__ == "__main__":
    main()
