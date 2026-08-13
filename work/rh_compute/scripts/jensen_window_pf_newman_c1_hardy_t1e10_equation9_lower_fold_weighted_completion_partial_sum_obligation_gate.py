#!/usr/bin/env python3
"""Reduce weighted completion to a cancellation-preserving partial-sum family."""

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

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
    "weighted_packets": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate.json",
    "variation_reassembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_variation_of_constants_roster_reassembly_gate.json",
}

L = 39_696
Q = 39_894
U = 40_094
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
    require(set(weights) == set(range(L, Q)) | set(range(Q + 1, U + 1)), "398-mode weight roster drift")
    weights[Q] = acb(0)
    return weights


def exact_identity_check() -> None:
    weights = {m: Fraction((m - L + 3) ** 2 + 1, 7 * (m - L + 2) + 5) for m in range(L, U + 1)}
    weights[Q] = Fraction(0)
    terms = {m: Fraction(5 * (m - L) - 11, 13 * (m - L + 1) + 17) for m in range(L, U + 1)}
    h4 = Fraction(7, 19)
    c4 = Fraction(-11, 23)
    g4 = h4 + c4 + sum(terms.values(), Fraction())
    prefixes = {k: sum((terms[m] for m in range(L, k + 1)), Fraction()) for k in range(L, U + 1)}
    remainders = {
        k: h4 + c4 + sum((terms[m] for m in range(k + 1, U + 1)), Fraction())
        for k in range(L, U + 1)
    }
    direct = sum((weights[m] * terms[m] for m in range(L, U + 1)), Fraction())
    abel = weights[U] * prefixes[U] + sum(
        ((weights[k] - weights[k + 1]) * prefixes[k] for k in range(L, U)), Fraction()
    )
    completed = weights[L] * g4 - weights[U] * remainders[U] - sum(
        ((weights[k] - weights[k + 1]) * remainders[k] for k in range(L, U)), Fraction()
    )
    require(direct == abel == completed, "399-mode completed Abel identity failed")


def coefficient_certificate(pairing: dict[str, Any]) -> dict[str, Any]:
    weights = reconstruct_weights(pairing)
    total_variation = sum((abs(weights[k + 1] - weights[k]).upper() for k in range(L, U)), arb(0))
    left_endpoint = abs(weights[L]).upper()
    right_endpoint = abs(weights[U]).upper()
    augmented_variation = left_endpoint + total_variation + right_endpoint
    edge_difference = abs(weights[U - 1] - weights[U])
    require(edge_difference.lower() > arb(0), "last two selected weights were not certified distinct")
    require(total_variation < arb("1.547052e-5"), "full weight variation exceeded certificate")
    require(augmented_variation < arb("3.093180e-5"), "augmented weight variation exceeded certificate")
    return {
        "extended_weight_definition": "w_m=a_m for m!=q and w_q=0",
        "mode_interval": [L, U],
        "mode_count": U - L + 1,
        "removed_event_mode": Q,
        "left_endpoint_absolute_ball": left_endpoint.str(PRECISION, more=True),
        "right_endpoint_absolute_ball": right_endpoint.str(PRECISION, more=True),
        "internal_total_variation_ball": total_variation.str(PRECISION, more=True),
        "endpoint_augmented_total_variation_ball": augmented_variation.str(PRECISION, more=True),
        "fixed_event_zero_countermodel_modes": [U - 1, U],
        "last_two_weight_difference_absolute_ball": edge_difference.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["coefficient_certificate"]
    return f"""# Weighted completion partial-sum obligation

Date: 2026-08-13

Status: exact cancellation-preserving Abel reduction and non-identifiability
guard; uniform completed partial-remainder bound remains open

Let `L={L}`, `q={Q}`, and `U={U}`.  Extend the 398 selected
leading-defect weights to an auxiliary complete 399-mode beta-minus-four
roster by

```text
w_m=a_m (m!=q),             w_q=0.                  (WC1)
```

For the unweighted beta-minus-four terms define

```text
S_k=sum_(m=L)^k G4_m,
R_k=H4+C4+sum_(m=k+1)^U G4_m.                       (WC2)
```

The completed projection `H4+sum G4+C4=g4_lambda(0)` gives, without
separating `H4` or `C4`,

```text
S_k=g4_lambda(0)-R_k.                               (WC3)
```

Finite Abel summation is exactly

```text
W_full=sum_(m=L)^U w_m G4_m
 =w_U S_U+sum_(k=L)^(U-1)(w_k-w_(k+1))S_k           (WC4)

 =w_L g4_lambda(0)-w_U R_U
  -sum_(k=L)^(U-1)(w_k-w_(k+1))R_k.                 (WC5)
```

An exact rational 399-mode calculation checks both equalities.  Each `R_k`
keeps the endpoint half-current, the zero/negative/outer-positive complement,
and the untouched upper suffix together.  Thus (WC5) does not split bare
endpoint exponentials or take termwise absolute values of the `G4_m`.

The certified weight budget is

```text
sum_(k=L)^(U-1)|w_(k+1)-w_k|
 <{c['internal_total_variation_ball']},

|w_L|+sum|Delta w|+|w_U|
 <{c['endpoint_augmented_total_variation_ball']}<3.093180e-5.          (WC6)
```

Consequently a uniform completed-family estimate

```text
max(|g4_lambda(0)|, max_(L<=k<=U)|R_k|) <= B_comp
```

would imply the cancellation-preserving bound

```text
|W_full| < 3.093180e-5 B_comp.                      (WC7)
```

This is not yet the selected-branch weighted sum.  Since
`G4_m=G_sel,m+G_opp,m`, exact algebra requires

```text
W_sel=W_full-W_opp,
W_opp=sum_(m=L)^U w_m G_opp,m.                      (WC7a)
```

The completed identity controls `W_full`; a separate cancellation-preserving
bound for `W_opp` is required before (WC7) can control `W_sel`.

The one global completion is not enough by itself.  The last two selected
weights satisfy

```text
|w_40093-w_40094|
 ={c['last_two_weight_difference_absolute_ball']}>0.                 (WC8)
```

Holding event zero fixed, replace `G4_40093` by `G4_40093+z` and
`G4_40094` by `G4_40094-z`.  The completed total and every other term are
unchanged, while `W_full` changes by `(w_40093-w_40094)z`.  Therefore the global
identity alone cannot identify or cancel the auxiliary weighted full roster.  This is a
logical non-identifiability result for that identity, not a claim that the
actual Fourier coefficients can be varied arbitrarily.

The next theorem targets are a uniform representation or bound for the
completed family `R_k` and a separate weighted opposite-branch estimate;
only their combination through (WC7a) reaches the selected branch.

Pi provenance: this gate introduces no new `pi`; all weights and completed
terms inherit `beta^3=pi C^2/8` and the prior beta-minus-four normalization.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No uniform completed partial-remainder or weighted opposite-branch bound,
selected-branch initial-data/source cancellation, exact finite-integral amplitude remainder, grouped 398-mode
splice, complete `Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
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
    exact_identity_check()
    certificate = coefficient_certificate(dependencies["event_pairing"])
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "weighted_completion_reduced_to_completed_partial_sum_family",
        "passed": True,
        "exact_algebra": {
            "weighted_full_sum": "W_full=sum_(m=L)^U w_m G4_m with w_q=0",
            "selected_full_opposite_decomposition": "W_sel=W_full-W_opp, where W_opp=sum w_m G_opp,m",
            "prefix": "S_k=sum_(m=L)^k G4_m",
            "completed_partial_remainder": "R_k=H4+C4+sum_(m=k+1)^U G4_m",
            "completion_relation": "S_k=g4_lambda(0)-R_k",
            "Abel_identity": "W=w_U S_U+sum_(k=L)^(U-1)(w_k-w_(k+1))S_k",
            "completed_Abel_identity": "W=w_L g4_lambda(0)-w_U R_U-sum_(k=L)^(U-1)(w_k-w_(k+1))R_k",
            "uniform_bound": "If max(|g4_lambda(0)|,max_k|R_k|)<=B_comp, then |W|<=B_comp(|w_L|+sum|Delta w|+|w_U|)",
        },
        "coefficient_certificate": certificate,
        "countermodel": {
            "event_zero_held_fixed": True,
            "perturbation": "G4_40093 -> G4_40093+z; G4_40094 -> G4_40094-z",
            "completed_total_change": "0",
            "weighted_sum_change": "(w_40093-w_40094)z",
            "scope": "Shows non-identifiability from the global completed identity alone; does not vary the actual Fourier model.",
        },
        "decision": {
            "global_unweighted_completion_alone_identifies_weighted_roster": False,
            "event_zero_term_needed_for_countermodel": False,
            "completed_partial_sum_family_exactly_reduces_weighted_roster": True,
            "uniform_completed_partial_remainder_bound_proved": False,
            "weighted_opposite_branch_bound_proved": False,
            "selected_branch_weighted_roster_bound_proved": False,
            "completed_initial_data_cancellation_proved": False,
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
            "Construct cancellation-preserving bounds for the completed family "
            "R_k=H4+C4+sum_(m=k+1)^U G4_m and for W_opp=sum w_m G_opp,m. "
            "Only then infer W_sel=W_full-W_opp; do not split H4 or C4."
        ),
        "proof_boundary": (
            "Exact auxiliary full-roster weighted-completion Abel reduction, certified coefficient variation, and a fixed-event-zero "
            "non-identifiability countermodel for the global identity alone. No uniform completed partial-remainder or "
            "weighted opposite-branch bound, selected-branch initial-data/source cancellation, exact finite-integral amplitude remainder, grouped 398-mode "
            "splice, complete Q_K-T or T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified weighted-completion partial-sum obligation", flush=True)


if __name__ == "__main__":
    main()
