#!/usr/bin/env python3
"""Insert the weighted corrected defect into the raw selector increment exactly."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "common_profile": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "weighted_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.json",
    "corrected_selected": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate.json",
    "source_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.json",
    "selector_increment": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate.json",
}

C = 159_577
Q = 39_894
FOLD_LO = 39_695
FOLD_HI = 39_894
AUGMENTED_LO = 39_696
AUGMENTED_HI = 40_094


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


def exact_rational_check() -> dict[str, str]:
    """Exercise both coefficient and Fourier-sector identities exactly."""

    kappa = Fraction(17, 19)
    exact_profile = {
        mode: Fraction(5 * (mode - Q) + 7, 211 + abs(mode - Q))
        for mode in range(AUGMENTED_LO, AUGMENTED_HI + 1)
    }
    opposite = {
        mode: Fraction(3 * (mode - Q) - 2, 223 + abs(mode - Q))
        for mode in exact_profile
    }
    weights = {
        mode: (Fraction() if mode == Q else Fraction(2 * (mode - Q) + 1, 100_003 + abs(mode - Q)))
        for mode in exact_profile
    }

    raw = {mode: kappa * exact_profile[mode] for mode in exact_profile}
    corrected = {
        mode: kappa * weights[mode] * (exact_profile[mode] - opposite[mode])
        for mode in exact_profile
    }
    companion = {
        mode: kappa
        * ((1 - weights[mode]) * exact_profile[mode] + weights[mode] * opposite[mode])
        for mode in exact_profile
    }
    require(all(raw[m] == corrected[m] + companion[m] for m in raw), "modewise companion identity failed")
    require(corrected[Q] == 0 and companion[Q] == raw[Q], "event-zero companion identity failed")

    raw_sum = sum(raw.values(), Fraction())
    corrected_sum = sum(corrected.values(), Fraction())
    companion_sum = sum(companion.values(), Fraction())
    require(raw_sum == corrected_sum + companion_sum, "roster companion identity failed")

    corrected_inner = sum((corrected[m] for m in range(AUGMENTED_LO, Q + 1)), Fraction())
    corrected_outer = sum((corrected[m] for m in range(Q + 1, AUGMENTED_HI + 1)), Fraction())
    require(corrected_sum == corrected_inner + corrected_outer, "corrected branch split failed")

    strip = {
        mode: Fraction(7 * mode - 11, 307 + abs(mode))
        for mode in range(-AUGMENTED_HI, AUGMENTED_HI + 1)
    }
    strip.update(raw)
    fold_increment = sum(
        (strip[mode] + strip[-mode] for mode in range(FOLD_LO, FOLD_HI + 1)),
        Fraction(),
    )
    raw_sector_remainder = (
        strip[FOLD_LO]
        + sum((strip[-mode] for mode in range(FOLD_LO, FOLD_HI + 1)), Fraction())
        - sum((strip[mode] for mode in range(FOLD_HI + 1, AUGMENTED_HI + 1)), Fraction())
    )
    require(fold_increment == raw_sum + raw_sector_remainder, "paired-sector identity failed")
    require(
        fold_increment == corrected_sum + companion_sum + raw_sector_remainder,
        "weighted increment decomposition failed",
    )
    paired_companion = (
        strip[FOLD_LO]
        + strip[-FOLD_LO]
        + sum(
            (companion[mode] + strip[-mode] for mode in range(AUGMENTED_LO, Q + 1)),
            Fraction(),
        )
    )
    require(fold_increment == corrected_inner + paired_companion, "inner-pair collapse failed")
    require(
        companion_sum + raw_sector_remainder == paired_companion - corrected_outer,
        "outer corrected cancellation failed",
    )

    outside = Fraction(23, 317)
    removed_source = fold_increment + outside
    require(
        removed_source == corrected_sum + companion_sum + raw_sector_remainder + outside,
        "completed source decomposition failed",
    )
    return {
        "carrier_fixture": str(kappa),
        "raw_roster_sum": str(raw_sum),
        "corrected_weighted_sum": str(corrected_sum),
        "coefficient_companion_sum": str(companion_sum),
        "corrected_inner_sum": str(corrected_inner),
        "corrected_outer_sum": str(corrected_outer),
        "paired_companion_sum": str(paired_companion),
        "raw_sector_remainder": str(raw_sector_remainder),
        "fold_increment": str(fold_increment),
        "completed_outside_remainder": str(outside),
    }


def exact_identities() -> dict[str, str]:
    return {
        "source_carrier": "kappa_C=e^(i*beta^3)2sqrt(2), beta^3=pi*C^2/8",
        "raw_strip_mode": "Z_m=kappa_C*G_ex,m",
        "weight_convention": "a_39894=0; other a_m are the certified leading stationary-defect weights",
        "corrected_mode": "c_m=kappa_C*a_m*(G_ex,m-G_opp,m)",
        "companion_mode": "b_m=kappa_C*((1-a_m)G_ex,m+a_m*G_opp,m)",
        "modewise_split": "Z_m=c_m+b_m",
        "roster_split": "Z_R=C_corr+B_comp",
        "corrected_roster": "C_corr=kappa_C*W_corr",
        "fold_increment": "D_fold=C_corr+B_comp+E_raw",
        "completed_source": "P_C=C_corr+B_comp+E_raw+C_out",
        "physical_split": "Delta Q_fold=R_corr+R_joint",
        "signed_component": "R_corr=2Re[e^(-i*pi/8)C_corr]=4sqrt(2)Re(W_corr)",
        "joint_remainder": "R_joint=2Re[e^(-i*pi/8)(B_comp+E_raw)]",
        "corrected_branch_split": "C_corr=C_inner+C_outer",
        "paired_companion": "J_pair=Z_39695+Z_-39695+sum_(m=39696)^39894(b_m+Z_-m)",
        "outer_cancellation": "B_comp+E_raw=J_pair-C_outer",
        "collapsed_increment": "D_fold=C_inner+J_pair",
        "collapsed_completed_source": "P_C=C_inner+J_pair+C_out",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["projection_certificate"]
    return f"""# Weighted defect and raw-strip companion decomposition

Date: 2026-08-13

Status: exact carrier-restored insertion into the selector increment proved;
the compulsory joint companion is not bounded

Let

```text
kappa_C=e^(i beta^3)2sqrt(2),       beta^3=pi*C^2/8,
Z_m=kappa_C G_ex,m.                                  (WD1)
```

The factor `kappa_C` restores both the exact equation-(9) half-mode
normalization and the common source carrier suppressed in the common-profile
Fourier coefficients.  Extend the certified leading-defect weights by the
existing event-zero convention `a_39894=0`, and put

```text
c_m=kappa_C a_m(G_ex,m-G_opp,m),
b_m=kappa_C[(1-a_m)G_ex,m+a_m G_opp,m].              (WD2)
```

Direct coefficient algebra, with no estimate, gives

```text
Z_m=c_m+b_m.                                         (WD3)
```

Thus on the augmented roster `R=39696..40094`,

```text
Z_R=sum_R Z_m=C_corr+B_comp,
C_corr=sum_R c_m=kappa_C W_corr,
B_comp=sum_R b_m.                                    (WD4)
```

The event-zero mode has `c_39894=0` and
`b_39894=Z_39894`; it has not disappeared.  Combining (WD4) with the exact
paired-sector identity of Section 11.405 proves

```text
D_fold=C_corr+B_comp+E_raw,                           (WD5)
P_C=C_corr+B_comp+E_raw+C_out.                       (WD6)
```

Equation (WD6) is the completed one-cell source identity.  It retains the
coefficient companion, lower edge, negative target modes, outer-positive
modes, and completed outside sector before any absolute value is taken.

The paired physical projection splits exactly as

```text
Delta Q_fold=R_corr+R_joint,
R_corr =2Re[e^(-i*pi/8)C_corr]=4sqrt(2)Re W_corr,
R_joint=2Re[e^(-i*pi/8)(B_comp+E_raw)].              (WD7)
```

The second equality uses
`exp(i*pi*(C^2-1)/8)=1`, with `(C^2-1)/8` even.  The existing interval
certificate therefore gives

```text
R_corr<{c['source_projected_corrected_real_part_upper']}<-1.263e-4,
|R_corr|<{c['source_projected_corrected_modulus_bound']}<2.314e-4. (WD8)
```

There is a sharper cancellation.  Split `C_corr=C_inner+C_outer` at event
zero and define

```text
J_pair=Z_39695+Z_-39695
       +sum_(m=39696)^39894(b_m+Z_-m).               (WD9)
```

Using `Z_m=c_m+b_m` on the 200 outer-positive modes gives

```text
B_comp+E_raw=J_pair-C_outer,
D_fold=C_inner+J_pair,
P_C=C_inner+J_pair+C_out.                            (WD10)
```

Thus the outer weighted correction cancels exactly from the operative fold
increment.  The signed theorem (WD8) concerns `C_inner+C_outer`; it does not
by itself determine the sign of `C_inner`.  The next quantitative step is to
certify the inner corrected projection and then estimate `J_pair` as one
paired object.  Bounding `B_comp` and `E_raw` separately is not licensed.

Pi provenance: `beta^3=pi*C^2/8` and `pi/8` are inherited from the exact
Kummer phase and odd-square half-domain reflection.  No fitted scale or
geometric insertion is used.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No bound for `R_joint`, signed fold increment, fixed-state fold residual,
all-corridor telescope, complete `Q_K-T` or `T_upper`, height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in DEPENDENCIES.items()
    }
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["common_profile"]["decision"]["exact_transformed_strip_common_profile_identity_proved"] is True,
        "exact common-profile identity drift",
    )
    require(
        dependencies["weighted_completion"]["coefficient_certificate"]["extended_weight_definition"]
        == "w_m=a_m for m!=q and w_q=0",
        "event-zero weight convention drift",
    )
    require(
        dependencies["source_projection"]["decision"]["source_carrier_times_half_projection_phase_is_exactly_one"] is True,
        "source carrier phase drift",
    )
    require(
        dependencies["source_projection"]["decision"]["corrected_selected_statistic_equals_unweighted_selector_strip"] is False,
        "weighted statistic was overidentified",
    )
    require(
        dependencies["selector_increment"]["decision"]["fold_increment_equals_signed_strip_pair_block"] is True,
        "selector increment identity drift",
    )
    require(
        dependencies["selector_increment"]["decision"]["common_source_carrier_retained_in_weighted_defect_coordinate"] is True,
        "selector carrier repair drift",
    )

    parity_integer = (C * C - 1) // 8
    require((C * C - 1) % 8 == 0 and parity_integer % 2 == 0, "odd-square carrier parity failed")
    exact_check = exact_rational_check()
    projection = dependencies["source_projection"]["certificate"]
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_weighted_defect_raw_strip_companion_decomposition_proved_joint_remainder_open",
        "passed": True,
        "scope": {
            "lower_selector_C": C,
            "event_zero_mode": Q,
            "fold_target_modes": [FOLD_LO, FOLD_HI],
            "augmented_positive_modes": [AUGMENTED_LO, AUGMENTED_HI],
            "augmented_mode_count": AUGMENTED_HI - AUGMENTED_LO + 1,
        },
        "carrier_parity": {
            "integer": parity_integer,
            "is_even": True,
            "identity": "exp(i*pi*(C^2-1)/8)=1",
        },
        "exact_identities": exact_identities(),
        "exact_rational_check": exact_check,
        "projection_certificate": {
            "source_projected_corrected_real_part_upper": projection["source_projected_corrected_real_part_upper"],
            "source_projected_corrected_modulus_bound": projection["source_projected_corrected_modulus_bound"],
            "sufficient_joint_remainder_upper_target": "1.263e-4",
        },
        "decision": {
            "raw_strip_mode_is_carrier_restored_exact_common_profile_coefficient": True,
            "modewise_weighted_defect_plus_companion_equals_raw_strip": True,
            "weighted_corrected_defect_inserted_into_selector_increment": True,
            "common_source_carrier_retained": True,
            "event_zero_retained_in_coefficient_companion": True,
            "completed_source_decomposition_proved": True,
            "corrected_component_projection_uniformly_negative": True,
            "outer_corrected_block_cancels_from_collapsed_fold_increment": True,
            "collapsed_fold_increment_equals_inner_corrected_plus_paired_companion": True,
            "total_corrected_projection_determines_inner_corrected_sign": False,
            "inner_corrected_projection_bounded": False,
            "paired_companion_projection_bounded": False,
            "coefficient_companion_plus_raw_sector_remainder_bounded": False,
            "signed_fold_increment_proved": False,
            "fixed_state_fold_residual_bounded": False,
            "all_corridor_telescope_proved": False,
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
            "Split the grouped corrected-defect quadrature into its inner modes 39696..39894 and outer modes "
            "39895..40094, including the split exact finite-t profile correction, to certify the projection of "
            "C_inner. Then bound J_pair in the collapsed identity D_fold=C_inner+J_pair without separating its pairs."
        ),
        "proof_boundary": (
            "Exact carrier-restored weighted-defect/raw-strip companion and completed selector-increment identities only. "
            "No bound for the joint companion-sector remainder, signed fold increment, fixed-state fold residual, "
            "all-corridor telescope, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, "
            "or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified exact weighted-defect/raw-strip companion decomposition; joint remainder remains open", flush=True)


if __name__ == "__main__":
    main()
