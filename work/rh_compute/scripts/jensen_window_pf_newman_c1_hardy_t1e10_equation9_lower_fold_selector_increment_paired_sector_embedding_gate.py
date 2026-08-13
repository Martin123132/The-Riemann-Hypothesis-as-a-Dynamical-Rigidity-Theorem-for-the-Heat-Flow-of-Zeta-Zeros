#!/usr/bin/env python3
"""Certify the exact selector increment in paired fold coordinates."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "selector_strip": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate.json",
    "paired_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "source_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.json",
    "target_allocation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_residual_target_allocation_embedding_obligation_gate.json",
}

C = 159_577
C_NEXT = 159_579
B = 5_122_423
FOLD_LO = 39_695
FOLD_HI = 39_894
AUGMENTED_LO = 39_696
AUGMENTED_HI = 40_094
OUTER_POSITIVE_LO = FOLD_HI + 1
OUTER_POSITIVE_HI = AUGMENTED_HI


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


def coefficient_balance() -> dict[str, Any]:
    fold_signed = Counter(range(FOLD_LO, FOLD_HI + 1))
    fold_signed.update(-mode for mode in range(FOLD_LO, FOLD_HI + 1))
    augmented = Counter(range(AUGMENTED_LO, AUGMENTED_HI + 1))
    defect = fold_signed.copy()
    defect.subtract(augmented)
    defect = Counter({mode: value for mode, value in defect.items() if value})

    expected = Counter({FOLD_LO: 1})
    expected.update({-mode: 1 for mode in range(FOLD_LO, FOLD_HI + 1)})
    expected.update({mode: -1 for mode in range(OUTER_POSITIVE_LO, OUTER_POSITIVE_HI + 1)})
    require(defect == expected, "raw selector/fold coefficient balance failed")
    require(sum(value for value in defect.values() if value > 0) == 201, "positive defect count drifted")
    require(-sum(value for value in defect.values() if value < 0) == 200, "negative defect count drifted")
    return {
        "fold_signed_mode_count": sum(fold_signed.values()),
        "augmented_positive_mode_count": sum(augmented.values()),
        "cancelled_positive_intersection": [AUGMENTED_LO, FOLD_HI],
        "cancelled_positive_intersection_count": FOLD_HI - AUGMENTED_LO + 1,
        "surviving_lower_edge": FOLD_LO,
        "surviving_negative_target_block": [-FOLD_HI, -FOLD_LO],
        "subtracted_outer_positive_block": [OUTER_POSITIVE_LO, OUTER_POSITIVE_HI],
        "positive_terms_after_cancellation": 201,
        "negative_terms_after_cancellation": 200,
    }


def exact_rational_check() -> dict[str, str]:
    # Independent finite data exercise the state difference, fixed-target
    # cancellation, and both completed-remainder forms without numerics.
    cutoff = 40_101
    strip = {mode: Fraction(7 * mode - 3, 113) for mode in range(-cutoff, cutoff + 1)}
    old = {mode: Fraction(11 * mode + 5, 127) for mode in range(-cutoff, cutoff + 1)}
    new = {mode: old[mode] - strip[mode] for mode in old}
    target = {mode: Fraction(5 * mode + 1, 131) for mode in range(FOLD_LO, FOLD_HI + 1)}

    fold_old = sum((old[r] + old[-r] - target[r] for r in target), Fraction())
    fold_new = sum((new[r] + new[-r] - target[r] for r in target), Fraction())
    fold_increment = fold_old - fold_new
    signed_strip = sum((strip[r] + strip[-r] for r in target), Fraction())
    require(fold_increment == signed_strip, "fixed-target fold increment failed")

    augmented = sum((strip[m] for m in range(AUGMENTED_LO, AUGMENTED_HI + 1)), Fraction())
    raw_defect = (
        strip[FOLD_LO]
        + sum((strip[-r] for r in range(FOLD_LO, FOLD_HI + 1)), Fraction())
        - sum((strip[m] for m in range(OUTER_POSITIVE_LO, OUTER_POSITIVE_HI + 1)), Fraction())
    )
    require(fold_increment == augmented + raw_defect, "raw boundary-block decomposition failed")

    endpoint_half = Fraction(19, 151)
    removed_source = endpoint_half + sum(strip.values(), Fraction())
    outside_fold = removed_source - fold_increment
    raw_roster_remainder = removed_source - augmented
    require(
        fold_increment - augmented == raw_roster_remainder - outside_fold,
        "raw completed-remainder difference identity failed",
    )
    return {
        "cutoff": str(cutoff),
        "fold_increment": str(fold_increment),
        "raw_augmented_sum": str(augmented),
        "raw_boundary_defect": str(raw_defect),
        "raw_roster_remainder": str(raw_roster_remainder),
        "completed_remainder_difference": str(raw_roster_remainder - outside_fold),
    }


def exact_identities() -> dict[str, str]:
    return {
        "per_mode_strip": "I_m^(C)(x)-I_m^(C+2)(x)=S_m(x)",
        "integrated_per_mode_strip": "K_m^(C)(t)-K_m^(C+2)(t)=Z_m(t), Z_m=int_0^(1/2)W_t S_m",
        "endpoint_half_strip": "K_H^(C)-K_H^(C+2)=Z_H",
        "complete_selector_cocycle": "P_C=Q_half^(C)-Q_half^(C+2)=Z_H+lim_sym sum_m Z_m",
        "fold_state": "F_C=sum_(r=39695)^39894[K_r^(C)+K_-r^(C)-tauhat_r]",
        "fold_increment": "D_fold=F_C-F_(C+2)=sum_(r=39695)^39894[Z_r+Z_-r]",
        "raw_augmented_roster": "Z_R=sum_(m=39696)^40094 Z_m",
        "raw_embedding_defect": "E_raw=D_fold-Z_R=Z_39695+sum_(r=39695)^39894 Z_-r-sum_(m=39895)^40094 Z_m",
        "outside_fold_completion": "C_out=P_C-D_fold",
        "raw_roster_completion": "R_raw=P_C-Z_R",
        "cancellation_safe_raw_remainder": "D_fold-Z_R=R_raw-C_out=E_raw",
        "weighted_defect_guard": "C_corr^K=e^(i*beta^3)2sqrt(2)W_corr retains the common source carrier, is weighted by a_m, and is not identified with the unweighted raw roster Z_R",
        "weighted_defect_projection": "2Re[e^(-i*pi/8)C_corr^K]=4sqrt(2)Re(W_corr) because exp(i*pi*(C^2-1)/8)=1",
        "physical_increment": "Delta Q_fold=2Re[exp(-i*pi/8)D_fold]",
    }


def render_note(artifact: dict[str, Any]) -> str:
    balance = artifact["coefficient_balance"]
    return f"""# Selector increment in paired fold coordinates

Date: 2026-08-13

Status: exact selector-cocycle and paired-sector embedding identity proved;
the increment remainder is not bounded

Fix the upper odd endpoint `B={B}` and let the lower selector advance from
`C={C}` to `C+2={C_NEXT}`.  If `I_m^(C)` is the finite-Poisson coefficient
before the advance, the integer character is unchanged by the unit shift, so

```text
I_m^(C)(x)-I_m^(C+2)(x)=S_m(x),
K_m^(C)(t)-K_m^(C+2)(t)=Z_m(t),
Z_m(t)=integral_0^(1/2) W_t(x)S_m(x)dx.              (SI1)
```

The endpoint half-current difference is `Z_H`.  The completed symmetric
Poisson theorem therefore gives the exact selector cocycle

```text
P_C=Q_half^(C)-Q_half^(C+2)
   =Z_H+lim_sym sum_m Z_m,                            (SI2)
```

where `P_C` is the single removed source term after the same half-domain
weight and integration.  This identity keeps the endpoint, zero, negative,
and outer-positive sectors attached.

For the fold-owned target block define

```text
F_C=sum_(r={FOLD_LO})^{FOLD_HI}
 [K_r^(C)+K_-r^(C)-tauhat_r].                         (SI3)
```

The target carrier is independent of the selector and cancels exactly in the
state difference.  Hence

```text
D_fold=F_C-F_(C+2)
      =sum_(r={FOLD_LO})^{FOLD_HI}(Z_r+Z_-r).         (SI4)
```

This is the natural location of the one-cell strip: it is an increment of
the paired residual, not a summand of one fixed selector state.

Let `Z_R=sum_(m={AUGMENTED_LO})^{AUGMENTED_HI}Z_m` be the raw positive
399-mode roster.  Exact coefficient cancellation leaves

```text
E_raw=D_fold-Z_R
 =Z_{FOLD_LO}+sum_(r={FOLD_LO})^{FOLD_HI}Z_-r
  -sum_(m={OUTER_POSITIVE_LO})^{OUTER_POSITIVE_HI}Z_m. (SI5)
```

The common positive intersection `{balance['cancelled_positive_intersection'][0]}..{balance['cancelled_positive_intersection'][1]}`
contains exactly {balance['cancelled_positive_intersection_count']} modes.  Thus the true raw mismatch is one lower edge,
the 200 negative target modes, and the 200 outer-positive modes.  The earlier
399-to-200 event allocation remains useful internal branch bookkeeping, but
it is not the Fourier-sector embedding (SI5).

There is also a cancellation-safe completed form of the raw mismatch:

```text
C_out=P_C-D_fold,       R_raw=P_C-Z_R,
D_fold-Z_R=R_raw-C_out=E_raw.                         (SI6)
```

Equation (SI6) is cancellation-safe: `R_raw` and `C_out` are completed
objects before their difference is estimated.  It does **not** insert the
corrected selected quantity from Section 11.403.  That quantity is

```text
C_corr^K=e^(i beta^3)2sqrt(2)W_corr,
W_corr=sum_m a_m(G_ex,m-G_opp,m),                     (SI7)
```

The factor `e^(i beta^3)` is mandatory: it is the common source carrier
suppressed in the common-profile coefficients.  Since
`exp(i*pi*(C^2-1)/8)=1`, its paired physical projection is exactly the
quantity in Section 11.403.  The term still carries the nonconstant
leading-defect weights `a_m`, whereas `Z_R` is the unweighted raw strip
roster.  Equal physical units do not make these objects equal.  The next
exact task is to decompose `Z_R` into (SI7) plus its compulsory weighted
companion before attempting a fixed-state sign theorem.

The physical increment is

```text
Delta Q_fold=2 Re[exp(-i*pi/8)D_fold].                (SI8)
```

The `pi/8` is inherited from odd-square Kummer reflection.  The selector
cocycle and coefficient cancellation introduce no new `pi` and no fitted
constant.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No weighted-defect-to-fixed-state embedding, no bound for `E_raw`, no signed
fold-increment theorem, no fixed-state fold residual bound, no all-corridor telescope, complete `Q_K-T` or
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")

    strip = dependencies["selector_strip"]
    paired = dependencies["paired_residual"]
    require(strip["decision"]["complete_symmetric_poisson_strip_reassembles_removed_term"] is True, "selector cocycle drift")
    require(strip["decision"]["positive_399_mode_roster_alone_sufficient"] is False, "399-mode sufficiency guard drift")
    require(paired["scope"]["fold_owned_modes"] == [FOLD_LO, FOLD_HI], "fold roster drift")
    require(paired["decision"]["target_carrier_subtracted_exactly_once"] is True, "target allocation drift")
    require(dependencies["completed_projection"]["decision"]["zero_negative_outer_positive_complement_retained_as_one_object"] is True, "completed complement drift")
    require(dependencies["source_projection"]["decision"]["corrected_selected_component_embedded_in_fold_owned_paired_residual"] is False, "embedding guard drift")
    require(dependencies["source_projection"]["decision"]["source_carrier_times_half_projection_phase_is_exactly_one"] is True, "source carrier phase drift")
    require(dependencies["target_allocation"]["decision"]["exact_399_to_200_group_allocation_proved"] is True, "target allocation dependency drift")
    require(dependencies["target_allocation"]["decision"]["selector_to_paired_residual_embedding_proved"] is False, "target allocation overpromotion")

    balance = coefficient_balance()
    rational = exact_rational_check()
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_selector_cocycle_and_paired_sector_embedding_identity_proved_remainder_open",
        "passed": True,
        "scope": {
            "fixed_upper_endpoint_B": B,
            "old_selector_C": C,
            "new_selector_C": C_NEXT,
            "fold_target_modes": [FOLD_LO, FOLD_HI],
            "augmented_positive_modes": [AUGMENTED_LO, AUGMENTED_HI],
        },
        "exact_identities": exact_identities(),
        "coefficient_balance": balance,
        "exact_rational_check": rational,
        "decision": {
            "per_mode_selector_difference_equals_strip_coefficient": True,
            "complete_selector_difference_is_removed_source_term": True,
            "fixed_target_carrier_cancels_in_selector_increment": True,
            "fold_increment_equals_signed_strip_pair_block": True,
            "raw_positive_399_to_signed_fold_defect_decomposed_exactly": True,
            "cancellation_safe_completed_remainder_identity_proved": True,
            "event_allocation_is_fourier_sector_embedding": False,
            "state_level_embedding_remainder_must_vanish_for_increment_route": False,
            "corrected_selected_weighted_defect_equals_raw_positive_roster": False,
            "common_source_carrier_retained_in_weighted_defect_coordinate": True,
            "corrected_selected_weighted_defect_embedded_in_fixed_state_residual": False,
            "increment_embedding_remainder_bounded": False,
            "signed_fold_increment_proved": False,
            "all_corridor_telescope_proved": False,
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
            "Decompose the unweighted raw strip Z_R exactly as the carrier-restored weighted corrected defect "
            "C_corr^K=e^(i beta^3)2sqrt(2)W_corr plus its compulsory coefficient companion, then combine that "
            "identity with E_raw before any sign estimate."
        ),
        "proof_boundary": (
            "Exact fixed-B selector cocycle, paired fold-sector increment, and raw roster remainder decomposition only. "
            "No embedding of the weighted corrected defect in the fixed-state residual, no bound for the raw increment "
            "remainder, signed fold-increment theorem, fixed-state fold residual, all-corridor telescope, "
            "complete Q_K-T or T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified exact selector cocycle and paired fold-sector increment; remainder remains open", flush=True)


if __name__ == "__main__":
    main()
