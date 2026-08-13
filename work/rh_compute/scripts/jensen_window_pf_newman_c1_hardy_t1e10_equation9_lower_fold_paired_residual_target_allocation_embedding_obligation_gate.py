#!/usr/bin/env python3
"""Certify the 399-to-200 target allocation and expose the open embedding."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_residual_target_allocation_embedding_obligation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "paired_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "ownership": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "common_profile": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "source_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.json",
    "selector_strip": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate.json",
    "event_zero": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_selector_center_recentered_transport_gate.json",
}

Q = 39_894
AUGMENTED_LO = 39_696
AUGMENTED_HI = 40_094
FOLD_LO = 39_695
FOLD_HI = 39_894


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


def allocation_groups() -> list[dict[str, Any]]:
    groups: list[dict[str, Any]] = []
    groups.append(
        {
            "target_mode": FOLD_LO,
            "kind": "selector_exchange_edge",
            "augmented_modes": [Q + 199, Q + 200],
        }
    )
    for target in range(Q - 198, Q):
        groups.append(
            {
                "target_mode": target,
                "kind": "reflected_event_pair",
                "augmented_modes": [target, 2 * Q - target],
            }
        )
    groups.append({"target_mode": Q, "kind": "event_zero", "augmented_modes": [Q]})

    target_modes = [row["target_mode"] for row in groups]
    augmented_modes = [mode for row in groups for mode in row["augmented_modes"]]
    require(target_modes == list(range(FOLD_LO, FOLD_HI + 1)), "target allocation is not ordered and exhaustive")
    require(len(augmented_modes) == len(set(augmented_modes)), "allocation groups overlap")
    require(set(augmented_modes) == set(range(AUGMENTED_LO, AUGMENTED_HI + 1)), "augmented roster is not exhaustive")
    require(len(groups) == 200 and len(augmented_modes) == 399, "allocation counts drifted")
    return groups


def exact_identity_check(groups: list[dict[str, Any]]) -> dict[str, str]:
    # Unrelated exact rational data check the four-term assembly without
    # relying on the symbolic spelling used in the theorem note.
    source = {
        mode: Fraction(7 * mode - 19, 11 * (mode - AUGMENTED_LO + 3))
        for mode in range(AUGMENTED_LO, AUGMENTED_HI + 1)
    }
    correction = {
        mode: Fraction(5 * mode + 3, 13 * (mode - AUGMENTED_LO + 5))
        for mode in range(AUGMENTED_LO, AUGMENTED_HI + 1)
        if mode != Q
    }
    target = {
        mode: Fraction(3 * mode - 2, 17 * (mode - FOLD_LO + 7))
        for mode in range(FOLD_LO, FOLD_HI + 1)
    }
    physical_pairs = {
        mode: Fraction(2 * mode + 1, 19 * (mode - FOLD_LO + 11))
        for mode in range(FOLD_LO, FOLD_HI + 1)
    }

    grouped_source = {
        row["target_mode"]: sum((source[m] for m in row["augmented_modes"]), Fraction())
        for row in groups
    }
    grouped_correction = {
        row["target_mode"]: sum((correction.get(m, Fraction()) for m in row["augmented_modes"]), Fraction())
        for row in groups
    }
    embedding = sum(physical_pairs.values(), Fraction()) - sum(source.values(), Fraction())
    fold_half = sum((physical_pairs[r] - target[r] for r in target), Fraction())
    corrected = sum(correction.values(), Fraction())
    event_zero = grouped_source[Q] - target[Q]
    target_allocation = sum(
        (
            grouped_source[r] - grouped_correction[r] - target[r]
            for r in range(FOLD_LO, Q)
        ),
        Fraction(),
    )
    require(fold_half == embedding + corrected + event_zero + target_allocation, "four-term fold assembly identity failed")
    return {
        "fold_half_block": "F_fold=sum_(r=39695)^39894 (K_r+K_-r-tauhat_r)",
        "embedding_remainder": "E_emb=sum_r(K_r+K_-r)-sum_(m=39696)^40094 Z_m",
        "corrected_selected_term": "c_m=e^(i*beta^3)2sqrt(2)a_m(G_ex,m-G_opp,m), C_corr=sum_(m!=39894)c_m=e^(i*beta^3)2sqrt(2)W_corr",
        "event_zero_term": "E_0=Z_39894-tauhat_39894",
        "target_allocation_term": "A_tar=sum_(r=39695)^39893[sum_(m in A_r)(Z_m-c_m)-tauhat_r]",
        "assembly": "F_fold=E_emb+C_corr+E_0+A_tar",
        "physical_projection": "Q_fold=2Re[e^(-i*pi/8)F_fold]",
        "identity_check": "independent exact rational values on 399 augmented and 200 target modes",
    }


def render_note(artifact: dict[str, Any]) -> str:
    return f"""# Fold-owned target allocation and embedding obligation

Date: 2026-08-13

Status: exact 399-to-200 allocation and four-term assembly identity proved;
the analytic selector-to-paired-residual embedding remains open

The exact paired residual contains the half-domain fold block

```text
F_fold=sum_(r={FOLD_LO})^{FOLD_HI}(K_r+K_-r-tauhat_r).       (FA1)
```

Its target roster has 200 modes.  The corrected selector calculation instead
uses the 399-mode augmented roster `{AUGMENTED_LO}..{AUGMENTED_HI}`.  Exact
event ordering gives a disjoint allocation of that roster to the target:

```text
A_39894={{39894}},
A_(39894-j)={{39894-j,39894+j}},       1<=j<=198,
A_39695={{40093,40094}}.                              (FA2)
```

Thus (FA2) consists of one event-zero group, 198 reflected pairs, and one
selector-exchange edge.  Its 200 target labels are exactly
`39695..39894`, and its 399 augmented modes are exactly
`39696..40094`, with no overlap or omission.  This is an allocation theorem,
not an analytic equality between a group and its target mode.

Let `Z_m` denote the augmented selector contribution after an exact
embedding into the paired-residual coordinates has been supplied, and put

```text
E_emb=sum_(r=39695)^39894(K_r+K_-r)
      -sum_(m=39696)^40094 Z_m.                       (FA3)
```

For `m!=39894`, restore the suppressed common source carrier and define

```text
c_m=e^(i beta^3)2sqrt(2)a_m(G_ex,m-G_opp,m),
C_corr=sum_(m!=39894)c_m=e^(i beta^3)2sqrt(2)W_corr.
```

Its projected value is exactly the quantity certified in Section 11.403.
Pure finite algebra then gives

```text
C_corr=sum_(m!=39894)c_m,
E_0=Z_39894-tauhat_39894,
A_tar=sum_(r=39695)^39893
 [sum_(m in A_r)(Z_m-c_m)-tauhat_r],

F_fold=E_emb+C_corr+E_0+A_tar.                        (FA4)
```

An independent exact-rational test verifies (FA4) on all 399 augmented and
200 target slots.  Applying the common half-domain projection gives

```text
Q_fold=2 Re[e^(-i*pi/8)F_fold].                       (FA5)
```

Section 11.403 fixes the scale and negative sign of the projection of
`C_corr`; it does not show `E_emb=0`.  In fact the existing completed
selector theorem retains the zero, negative, and outer-positive complement
and endpoint half-current as one nonlocal object.  The paired-target theorem
also states explicitly that `tauhat_r` is an algebraic allocation, not a
modewise half-Gamma identity.  Therefore neither completion may be localized
to (FA1) by index counting alone.

The event-zero gate controls only the exact-minus-beta-minus-four transport
at mode 39894; it does not bound `E_0`.  No current gate bounds `A_tar` or
gives an exact source formula for `E_emb`.  The first decisive obligation is
to derive (FA3) from the original paired finite-Poisson current and the
completed one-cell selector identity, preserving all complement sectors.
Only then should `E_emb+E_0+A_tar` be bounded jointly with the signed
`C_corr` projection.

Pi provenance: `pi/8` in (FA5) is the exact odd-square half-Kummer
reflection phase.  The allocation (FA2)--(FA4) introduces no new `pi` and no
fitted constant.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No selector-to-paired-residual embedding, event-zero target residual,
target-allocation bound, complete fold-owned residual, complete `Q_K-T` or
`T_upper`, all-corridor or height-uniform theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")

    paired = dependencies["paired_residual"]
    require(paired["scope"]["fold_owned_modes"] == [FOLD_LO, FOLD_HI], "paired fold roster drift")
    require(paired["decision"]["target_allocation_is_modewise_half_Gamma_bulk_claim"] is False, "target allocation guard drift")
    require(dependencies["ownership"]["decision"]["mode_ownership_disjoint_and_exhaustive"] is True, "ownership drift")
    require(dependencies["event_pairing"]["decision"]["remaining_398_terms_equal_198_event_pairs_plus_one_selector_edge_pair"] is True, "event pairing drift")
    require(dependencies["completed_projection"]["decision"]["zero_negative_outer_positive_complement_retained_as_one_object"] is True, "completion locality guard drift")
    require(dependencies["common_profile"]["decision"]["complete_physical_source_initial_data_splice_proved"] is False, "physical splice unexpectedly changed")
    require(dependencies["source_projection"]["decision"]["corrected_selected_component_embedded_in_fold_owned_paired_residual"] is False, "source projection embedding guard drift")
    require(dependencies["selector_strip"]["decision"]["positive_399_mode_roster_alone_sufficient"] is False, "selector completion guard drift")

    groups = allocation_groups()
    event_pairs = dependencies["event_pairing"]["roster"]["consecutive_event_mode_pairs"]
    require(
        {tuple(row["augmented_modes"]) for row in groups if row["kind"] == "reflected_event_pair"}
        == {tuple(pair) for pair in event_pairs},
        "allocation does not match certified event pairs",
    )
    exact = exact_identity_check(groups)
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_399_to_200_target_allocation_and_four_term_assembly_proved_embedding_open",
        "passed": True,
        "scope": {
            "augmented_selector_modes": [AUGMENTED_LO, AUGMENTED_HI],
            "augmented_mode_count": 399,
            "fold_target_modes": [FOLD_LO, FOLD_HI],
            "fold_target_count": 200,
            "event_zero_mode": Q,
        },
        "allocation_groups": groups,
        "exact_assembly": exact,
        "decision": {
            "exact_399_to_200_group_allocation_proved": True,
            "allocation_groups_disjoint_and_exhaustive": True,
            "target_carrier_subtracted_exactly_once_per_fold_mode": True,
            "four_term_fold_assembly_identity_proved": True,
            "corrected_selected_projection_scale_and_sign_available": True,
            "selector_to_paired_residual_embedding_proved": False,
            "embedding_remainder_zero_proved": False,
            "event_zero_full_target_residual_bounded": False,
            "target_allocation_remainder_bounded": False,
            "complete_fold_owned_residual_bounded": False,
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
            "Derive E_emb directly from the original paired finite-Poisson current and the completed exact one-cell selector "
            "identity. Keep the zero, negative, outer-positive, and endpoint-half sectors attached, then combine E_emb, "
            "event zero, and target allocation before using the signed corrected projection."
        ),
        "proof_boundary": (
            "Exact combinatorial target allocation and formal four-term fold assembly only. No selector-to-paired-residual "
            "embedding, event-zero target residual, target-allocation bound, complete fold-owned residual, complete Q_K-T or "
            "T_upper, all-corridor or height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified exact 399-to-200 target allocation; analytic embedding remains open", flush=True)


if __name__ == "__main__":
    main()
