#!/usr/bin/env python3
"""Certify a uniform completed cutoff-family bound by one Dirichlet kernel."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_completed_cutoff_family_uniform_dirichlet_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "branch_envelope": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate.json",
    "weighted_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.json",
}

A = 159_577
MODE_COUNT = 399
PRECISION = 100


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


def upper_ball(text: str) -> arb:
    return arb(text).upper()


def profile_certificate(branch: dict[str, Any], weighted: dict[str, Any]) -> dict[str, Any]:
    pi = arb.pi()
    beta = (pi * arb(A) ** 2 / 8) ** (arb(1) / 3)
    epsilon = beta**-2
    y_max = pi * A / (2 * beta)
    lambda_max = pi / (16 * beta)
    x_max = y_max + lambda_max

    envelope = branch["certificate"]
    certified_x_min = upper_ball(envelope["X_min_ball"])
    certified_x_max = arb(envelope["X_max_ball"]).lower()
    require(lambda_max.upper() >= certified_x_min, "near-zero Airy panel does not meet the saved atlas")
    require(x_max.upper() < certified_x_max, "full selector strip exceeds the saved Airy atlas")

    # Directly close the tiny interval omitted by the positive-X event atlas.
    near_zero_x = arb(lambda_max.upper() / 2, lambda_max.upper() / 2)
    ai, ai_prime, _, _ = acb(-near_zero_x).airy()
    airy_bound = max(upper_ball(envelope["W_modulus_bound_ball"]), abs(ai).upper())
    airy_derivative_bound = max(upper_ball(envelope["W_X_modulus_bound_ball"]), abs(ai_prime).upper())

    lam = lambda_max.upper()
    y = y_max.upper()
    u_correction = epsilon.upper() * (13 * lam + 3 * y) / 60
    u_correction += epsilon.upper() ** 2 * (
        448 * lam**5
        + 280 * lam**2 * y**3
        + 4565 * lam**2
        + 105 * lam * y**4
        + 30 * lam * y
        + 63 * y**5
        + 405 * y**2
    ) / 50400
    u_absolute = 1 + u_correction
    v_absolute = epsilon.upper() * (8 * lam**2 + 4 * lam * y + 3 * y**2) / 60
    v_absolute += epsilon.upper() ** 2 * (
        40 * lam**3 + 20 * lam**2 * y + lam * y**2 + 9 * y**3 + 27
    ) / 1680
    profile_bound = y * (u_absolute * airy_bound + v_absolute * airy_derivative_bound)

    dirichlet_l1 = 1 + (pi * MODE_COUNT / 2).log()
    completed_family_bound = profile_bound * (1 + dirichlet_l1)
    variation = upper_ball(weighted["coefficient_certificate"]["endpoint_augmented_total_variation_ball"])
    weighted_roster_bound = completed_family_bound * variation
    require(profile_bound < arb("82.66"), "profile supremum exceeded 82.66")
    require(dirichlet_l1 < arb("7.441"), "Dirichlet L1 majorant exceeded 7.441")
    require(completed_family_bound < arb("698"), "completed cutoff family exceeded 698")
    require(weighted_roster_bound < arb("0.02159"), "weighted completed roster exceeded 0.02159")

    return {
        "beta_ball": beta.str(PRECISION, more=True),
        "ordinary_top_lambda_interval": ["0", lambda_max.str(PRECISION, more=True)],
        "selector_y_max_ball": y_max.str(PRECISION, more=True),
        "profile_X_interval": ["0", x_max.str(PRECISION, more=True)],
        "near_zero_Ai_absolute_ball": abs(ai).str(PRECISION, more=True),
        "near_zero_Ai_prime_absolute_ball": abs(ai_prime).str(PRECISION, more=True),
        "Ai_absolute_bound_ball": airy_bound.str(PRECISION, more=True),
        "Ai_prime_absolute_bound_ball": airy_derivative_bound.str(PRECISION, more=True),
        "U_absolute_bound_ball": u_absolute.str(PRECISION, more=True),
        "V_absolute_bound_ball": v_absolute.str(PRECISION, more=True),
        "g4_profile_supremum_bound_ball": profile_bound.str(PRECISION, more=True),
        "maximum_one_sided_Dirichlet_L1_bound_ball": dirichlet_l1.str(PRECISION, more=True),
        "uniform_completed_cutoff_family_bound_ball": completed_family_bound.str(PRECISION, more=True),
        "weighted_completed_roster_bound_ball": weighted_roster_bound.str(PRECISION, more=True),
        "Dirichlet_L1_derivation": (
            "For D_N(s)=sum_(j=0)^(N-1)exp(-2pi i j s), symmetry, sin(pi s)>=2s on [0,1/2], "
            "and |sin(pi N s)|<=min(pi N s,1) give Integral_0^1|D_N|<=1+log(pi N/2)."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Uniform completed cutoff-family Dirichlet bound

Date: 2026-08-13

Status: rigorous ordinary-top-corridor bound for every beta-minus-four
completed cutoff; exact finite-integral remainder and final splice remain open

For `L=39696`, `U=40094`, and `L<=k<=U`, Section 11.393 defines

```text
S_k=sum_(m=L)^k G4_m,
R_k=H4+C4+sum_(m=k+1)^U G4_m
   =g4_lambda(0)-S_k.                                (CD1)
```

The coefficient formula from the completed projection gives

```text
G4_(39895+n)=Integral_0^1 g4_lambda(s)e^(-2pi i n s)ds.
```

Hence a prefix of length `N=k-L+1` is one contiguous Fourier packet,

```text
S_k=Integral_0^1 g4_lambda(s)e^(2pi i 199s)D_N(s)ds,
D_N(s)=sum_(j=0)^(N-1)e^(-2pi i j s).                (CD2)
```

No coefficient is bounded separately.  For `0<=s<=1/2`,
`sin(pi s)>=2s` and `|sin(pi Ns)|<=min(pi Ns,1)`.  Splitting at
`s=1/(pi N)` and using symmetry proves

```text
Integral_0^1 |D_N(s)|ds <=1+log(pi N/2)
 <=1+log(399pi/2)
 ={c['maximum_one_sided_Dirichlet_L1_bound_ball']}<7.441.              (CD3)
```

On the ordinary half of the selector top cell,

```text
0<=lambda<=pi/(16beta),
0<=y<=Y=pi A/(2beta),
0<=X=lambda+y<117.                                   (CD4)
```

The saved Airy modulus atlas covers the range above its first positive
point; one direct interval closes the remaining interval down to `X=0`.
Triangle evaluation of the exact beta-minus-four polynomials gives

```text
|U|<{c['U_absolute_bound_ball']},
|V|<{c['V_absolute_bound_ball']},
sup_s |g4_lambda(s)|
 <{c['g4_profile_supremum_bound_ball']}<82.66.                         (CD5)
```

Combining (CD1)--(CD5) proves simultaneously for all 399 cutoffs

```text
max(|g4_lambda(0)|,max_k|R_k|)
 <{c['uniform_completed_cutoff_family_bound_ball']}<698.              (CD6)
```

Inserting the certified endpoint-augmented weight variation from Section
11.393 gives a cancellation-preserving uniform bound for the auxiliary
leading-defect-weighted full beta-minus-four roster:

```text
|W_full|=|sum_(m=L)^U w_m G4_m|
 <{c['weighted_completed_roster_bound_ball']}<0.02159.                 (CD7)
```

This closes the full-roster partial-completion information deficit at a
conservative scale.  It is not yet a selected-branch estimate: exact branch
algebra gives `W_sel=W_full-W_opp`, and `W_opp=sum w_mG_opp,m` remains to be
bounded.  It also does not identify the exact finite-integral amplitude
remainder with `w_m` or splice (CD7) into the branch initial-data and source
coordinates of Section 11.392.  Sharpening (CD3) or exploiting the profile's
variation may reduce the constant, but is not needed for validity.

Pi provenance: `Y`, the height interval, and the Fourier kernel inherit the
exact identities `beta^3=pi A^2/8` and `hY=2pi`; (CD3) uses the standard
Fourier exponential `e^(-2pi i n s)`.  No fitted circle constant enters.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No weighted opposite-branch estimate, selected-branch weighted bound, exact
finite-integral amplitude remainder, completed branch initial-data or source
splice, complete `Q_K-T` or `T_upper`, all-corridor theorem,
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
    certificate = profile_certificate(dependencies["branch_envelope"], dependencies["weighted_completion"])
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "uniform_completed_cutoff_family_bounded_on_ordinary_top_corridor",
        "passed": True,
        "certificate": certificate,
        "exact_algebra": {
            "prefix_kernel": "S_k=Integral_0^1 g4_lambda(s) exp(2pi i 199s) D_(k-L+1)(s) ds",
            "completed_cutoff": "R_k=g4_lambda(0)-S_k=H4+C4+sum_(m=k+1)^U G4_m",
            "cutoff_count": MODE_COUNT,
            "no_termwise_coefficient_norm": True,
        },
        "decision": {
            "uniform_completed_partial_remainder_bound_proved": True,
            "auxiliary_full_beta4_weighted_roster_bound_proved": True,
            "weighted_opposite_branch_bound_proved": False,
            "selected_branch_weighted_roster_bound_proved": False,
            "exact_finite_integral_amplitude_remainder_proved": False,
            "completed_branch_initial_data_source_splice_proved": False,
            "all_height_corridors_proved": False,
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
            "Bound W_opp=sum w_m G_opp,m by the event-ordered opposite-branch reflection/Dirichlet structure. "
            "Then combine W_sel=W_full-W_opp before addressing the exact finite-integral amplitude remainder."
        ),
        "proof_boundary": (
            "Uniform ordinary-top-corridor completed beta^-4 cutoff-family and auxiliary full-roster weighted bounds only. "
            "No weighted opposite-branch or selected-branch bound, exact finite-integral amplitude remainder, completed branch initial-data/source splice, complete Q_K-T or "
            "T_upper, all-corridor theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified uniform completed cutoff-family Dirichlet bound", flush=True)


if __name__ == "__main__":
    main()
