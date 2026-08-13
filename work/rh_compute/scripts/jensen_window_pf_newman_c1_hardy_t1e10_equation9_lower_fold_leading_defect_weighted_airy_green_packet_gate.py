#!/usr/bin/env python3
"""Certify leading-defect weighted Airy Green packets by finite Abel variation."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
    "grouped_green_packet": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate.json",
}

Q = 39_894
MINUS_MODES = tuple(range(Q - 198, Q))
PLUS_MODES = tuple(range(Q + 1, Q + 201))
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


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def reconstruct_defects(pairing: dict[str, Any]) -> dict[int, acb]:
    certificate = pairing["leading_multiplier_certificate"]
    defects: dict[int, acb] = {}
    for row in certificate["pair_rows"]:
        coherent = parse_complex(row["coherent_coefficient_ball"])
        leakage = parse_complex(row["conjugacy_leakage_coefficient_ball"])
        defects[int(row["minus_mode"])] = coherent + leakage
        defects[int(row["plus_mode"])] = (coherent - leakage).conjugate()
    for row in certificate["edge_rows"]:
        defects[int(row["mode"])] = parse_complex(row["leading_defect_ball"])
    require(set(defects) == set(MINUS_MODES) | set(PLUS_MODES), "leading-defect roster reconstruction failed")
    return defects


def branch_certificate(modes: tuple[int, ...], defects: dict[int, acb]) -> dict[str, Any]:
    values = [defects[mode] for mode in modes]
    magnitudes = [abs(value).upper() for value in values]
    differences = [abs(values[index + 1] - values[index]).upper() for index in range(len(values) - 1)]
    difference_prefix = [arb(0)]
    for value in differences:
        difference_prefix.append(difference_prefix[-1] + value)

    maximum_coefficient_variation = arb(0)
    maximum_effective_variation = arb(0)
    maximum_indices = (modes[0], modes[0])
    for start in range(len(modes)):
        block_maximum = arb(0)
        for stop in range(start, len(modes)):
            block_maximum = max(block_maximum, magnitudes[stop])
            total_variation = difference_prefix[stop] - difference_prefix[start]
            coefficient_variation = min(
                magnitudes[start] + total_variation,
                magnitudes[stop] + total_variation,
            )
            effective_variation = coefficient_variation + block_maximum
            maximum_coefficient_variation = max(maximum_coefficient_variation, coefficient_variation)
            if effective_variation > maximum_effective_variation:
                maximum_effective_variation = effective_variation
                maximum_indices = (modes[start], modes[stop])

    return {
        "mode_interval": [modes[0], modes[-1]],
        "mode_count": len(modes),
        "maximum_coefficient_absolute_ball": max(magnitudes).str(PRECISION, more=True),
        "maximum_any_contiguous_coefficient_Abel_variation_ball": maximum_coefficient_variation.str(PRECISION, more=True),
        "maximum_any_contiguous_modulus_adjusted_variation_ball": maximum_effective_variation.str(PRECISION, more=True),
        "maximum_modulus_adjusted_variation_candidate_modes": list(maximum_indices),
    }


def weighted_packet_certificate(pairing: dict[str, Any], grouped: dict[str, Any]) -> dict[str, Any]:
    defects = reconstruct_defects(pairing)
    minus = branch_certificate(MINUS_MODES, defects)
    plus = branch_certificate(PLUS_MODES, defects)
    minus_effective = arb(minus["maximum_any_contiguous_modulus_adjusted_variation_ball"]).upper()
    plus_effective = arb(plus["maximum_any_contiguous_modulus_adjusted_variation_ball"]).upper()
    total_effective = minus_effective + plus_effective

    packet = grouped["certificate"]
    unweighted_kernel = arb(packet["grouped_any_contiguous_Green_kernel_bound_ball"]).upper()
    unweighted_derivative = arb(packet["grouped_any_contiguous_derivative_Green_kernel_bound_ball"]).upper()
    weighted_kernel = (unweighted_kernel * total_effective).upper()
    weighted_derivative = (unweighted_derivative * total_effective).upper()
    require(weighted_kernel < arb("7.58e-7"), "weighted two-branch K packet exceeded 7.58e-7")
    require(weighted_derivative < arb("0.001634"), "weighted two-branch partial_s K packet exceeded 0.001634")

    return {
        "coefficient_definition": "a_m=[R_A(s_m,mu) exp(i beta^3 F(s_m,mu))-1]/sqrt(m)",
        "height_interval": "t*-pi/16<=t<=t*",
        "minus_branch": minus,
        "plus_branch": plus,
        "two_branch_modulus_adjusted_variation_sum_ball": total_effective.str(PRECISION, more=True),
        "weighted_two_branch_K_packet_bound_ball": weighted_kernel.str(PRECISION, more=True),
        "weighted_two_branch_partial_s_K_packet_bound_ball": weighted_derivative.str(PRECISION, more=True),
        "weighted_Abel_identity": (
            "For z partial sums bounded by B and monotone positive M_j, |sum a_j M_j z_j| is at most "
            "B M_max [min(|a_left|,|a_right|)+sum|Delta a|+max|a|] on each contiguous block."
        ),
        "modulus_guard": (
            "The added max|a| term pays for the full monotone variation of M_j without assuming a lower modulus bound."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    minus = c["minus_branch"]
    plus = c["plus_branch"]
    return f"""# Leading-defect weighted Airy Green packets

Date: 2026-08-13

Status: rigorous top-corridor weighted packet interface for the certified
leading stationary defects; no finite-integral amplitude remainder or
completed endpoint splice

Retain the exact leading-defect coefficients from Section 11.385,

```text
a_m=[R_A(s_m,mu)exp(i beta^3F(s_m,mu))-1]/sqrt(m).    (WG1)
```

The saved coherent/leakage balls reconstruct every coefficient on the two
selected branches, with no new midpoint evaluation.  For a contiguous block
`I=[a,b]`, define

```text
V_I=min(|a_a|,|a_b|)+sum_(j=a)^(b-1)|a_(j+1)-a_j|.  (WG2)
```

Finite interval enumeration of every branch block proves

```text
minus branch: V_I <{minus['maximum_any_contiguous_coefficient_Abel_variation_ball']},
plus branch:  V_I <{plus['maximum_any_contiguous_coefficient_Abel_variation_ball']}. (WG3)
```

If `M_j=M(-r_j)`, its positive monotonicity gives, in either orientation,

```text
V_I(aM)/M_max <= V_I(a)+max_(j in I)|a_j|.            (WG4)
```

The final term pays for the entire monotone modulus variation; no lower
bound for `M` and no constant-amplitude replacement is used.  Exhausting all
blocks yields

```text
minus effective variation
 <{minus['maximum_any_contiguous_modulus_adjusted_variation_ball']},
plus effective variation
 <{plus['maximum_any_contiguous_modulus_adjusted_variation_ball']},
two-branch sum
 <{c['two_branch_modulus_adjusted_variation_sum_ball']}.               (WG5)
```

Applying (WG4) to the exact phase packets in Section 11.390 proves, uniformly
over every active contiguous block and the complete top corridor,

```text
sum of the two branchwise weighted K-packet bounds
 <{c['weighted_two_branch_K_packet_bound_ball']}<7.58e-7,

sum of the two branchwise weighted partial_s K-packet bounds
 <{c['weighted_two_branch_partial_s_K_packet_bound_ball']}<0.001634.   (WG6)
```

No termwise Green-kernel absolute sum occurs in (WG6).  The first number is
small, but it is not yet an additive residual budget: it must be integrated
against the correctly normalized completed source.  The second number must
still be multiplied by the composed initial-data channel and is not itself
small enough to discard.  Moreover, `a_m` is only the certified leading
stationary defect; the exact finite-integral amplitude and endpoint/Fresnel
remainder is not represented by (WG1).

Pi provenance is inherited from `beta^3=pi C^2/8`, the exact stationary
action bridge, and the Airy Wronskian packet in Section 11.390.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No completed source integral, initial-data cancellation, exact finite-
integral selected/logistic match, grouped 398-mode splice, complete
`Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["event_pairing"].get("passed") is True, "event-pairing dependency failed")
    require(dependencies["grouped_green_packet"].get("passed") is True, "grouped Green dependency failed")

    ctx.dps = PRECISION
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "leading_defect_weighted_Airy_Green_packets_certified_on_top_corridor",
        "passed": True,
        "certificate": weighted_packet_certificate(
            dependencies["event_pairing"], dependencies["grouped_green_packet"]
        ),
        "decision": {
            "all_saved_leading_defects_reconstructed_from_pair_balls": True,
            "all_contiguous_branch_coefficient_variations_enumerated": True,
            "Airy_modulus_variation_retained": True,
            "weighted_two_branch_K_packet_below_7_58e_minus_7": True,
            "weighted_two_branch_derivative_packet_below_0_001634": True,
            "finite_integral_amplitude_remainder_proved": False,
            "completed_endpoint_cancellation_proved": False,
            "grouped_398_mode_splice_proved": False,
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
            "Insert the weighted K and partial_s K packets into the exact branchwise variation-of-constants reassembly. "
            "Identify and compose the initial data and source integrals with the beta^-4 Poisson completion before norms; "
            "derive a separate coefficient-variation certificate for the finite-integral amplitude remainder."
        ),
        "proof_boundary": (
            "Leading stationary-defect weighted packets only. No exact finite-integral amplitude remainder, completed "
            "source or initial-data cancellation, grouped 398-mode splice, complete Q_K-T or T_upper, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified leading-defect weighted Airy Green packets", flush=True)


if __name__ == "__main__":
    main()
