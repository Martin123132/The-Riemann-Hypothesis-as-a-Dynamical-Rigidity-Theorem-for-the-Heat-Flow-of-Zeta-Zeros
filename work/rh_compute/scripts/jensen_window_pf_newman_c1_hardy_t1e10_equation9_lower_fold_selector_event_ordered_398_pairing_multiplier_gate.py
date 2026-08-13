#!/usr/bin/env python3
"""Remove event zero and pair the remaining selector carriers in event order."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "event0_extension": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate.json",
    "stationary_bridge": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate.json",
    "event_handoff": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.json",
}

C = 159_577
Q = (C - 1) // 4
PAIR_COUNT = 198
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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def action_defect(s: arb, mu: arb) -> arb:
    root = (2 * s + mu - 1).sqrt()
    return (
        3 * s * s
        - 8 * s * root
        + 12 * s
        - 4 * mu * root
        + 6 * mu * s.log()
        - 3 * mu * (1 - mu).log()
        + 9 * mu
        + 4 * root
        - 6 * s.log()
        + 3 * (1 - mu).log()
        - 11
    ) / 6


def center_action_defect(s: arb) -> arb:
    root = (2 * s - 1).sqrt()
    return (3 * s * s - 8 * s * root + 12 * s + 4 * root - 6 * s.log() - 11) / 6


def height_derivative_defect(s: arb, mu: arb) -> arb:
    root = (2 * s + mu - 1).sqrt()
    numerator = 4 * s + 2 * mu - 2 * root * s.log() + root * (1 - mu).log() - 2 * root - 2
    return -numerator / (2 * root)


def leading_defect(mode: int, c: arb, beta: arb, mu: arb) -> acb:
    s = 4 * arb(mode) / c
    ratio = ((1 - mu) * (2 * s + mu - 1)) ** (arb(1) / 4) / s.sqrt()
    center_phase = beta**3 * center_action_defect(s)
    phase_radius = beta**3 * mu.upper() * abs(height_derivative_defect(s, mu)).upper()
    phase = arb(center_phase, phase_radius)
    multiplier = ratio * acb(0, phase).exp()
    return (multiplier - 1) / arb(mode).sqrt()


def point_leading_defect(mode: int, c: arb, beta: arb, mu: arb) -> acb:
    s = 4 * arb(mode) / c
    ratio = ((1 - mu) * (2 * s + mu - 1)) ** (arb(1) / 4) / s.sqrt()
    phase = beta**3 * action_defect(s, mu)
    return (ratio * acb(0, phase).exp() - 1) / arb(mode).sqrt()


def panel_leading_defect(mode: int, c: arb, beta: arb, mu_center: arb, mu_ball: arb, mu_half_width: arb) -> acb:
    s = 4 * arb(mode) / c
    ratio = ((1 - mu_ball) * (2 * s + mu_ball - 1)) ** (arb(1) / 4) / s.sqrt()
    phase_center = beta**3 * action_defect(s, mu_center)
    phase_radius = beta**3 * mu_half_width * abs(height_derivative_defect(s, mu_ball)).upper()
    phase = arb(phase_center, phase_radius)
    return (ratio * acb(0, phase).exp() - 1) / arb(mode).sqrt()


def roster_certificate() -> dict[str, Any]:
    original = {
        *(('C_minus', 4 * j + 1) for j in range(199)),
        *(('C_plus', -(4 * j + 3)) for j in range(200)),
    }
    event_zero = ("C_minus", 1)
    remaining = original - {event_zero}
    paired: set[tuple[str, int]] = set()
    mode_pairs = []
    for j in range(1, PAIR_COUNT + 1):
        minus_term = ("C_minus", 4 * j + 1)
        plus_term = ("C_plus", -(4 * j - 1))
        paired.update((minus_term, plus_term))
        mode_pairs.append([Q - j, Q + j])
    edge = {("C_plus", -795), ("C_plus", -799)}
    require(paired.isdisjoint(edge), "event pairs overlap selector edge")
    require(remaining == paired | edge, "398-term event-order pairing does not reconstruct the roster")
    modes = {mode for pair in mode_pairs for mode in pair} | {Q + 199, Q + 200}
    require(modes == set(range(39_696, 40_095)) - {Q}, "mode projection drift")
    return {
        "original_selected_term_count": len(original),
        "removed_event_zero": {"mode": Q, "detuning_label": -1, "branch": "C_minus", "exponent": 1},
        "remaining_term_count": len(remaining),
        "consecutive_event_pair_count": PAIR_COUNT,
        "consecutive_event_mode_pairs": mode_pairs,
        "selector_exchange_edge_modes": [Q + 199, Q + 200],
        "selector_exchange_edge_labels": [795, 799],
        "pairwise_identity": (
            "C_-*exp(i*(4j+1)x)+C_+*exp(-i*(4j-1)x)="
            "2*exp(i*x)*Re(C_-*exp(i*4j*x)), 1<=j<=198"
        ),
        "edge_identity": (
            "C_+*(exp(-i*795x)+exp(-i*799x))="
            "2*C_+*exp(-i*797x)*cos(2x)"
        ),
        "closed_remaining_projection": (
            "B_rem=2*exp(i*x)*sin(396x)/sin(2x)*Re(C_-*exp(i*398x))"
            "+2*C_+*exp(-i*797x)*cos(2x)"
        ),
        "removable_quotient_rule": "sin(396x)/sin(2x) is its original 198-term finite exponential sum",
    }


def multiplier_certificate() -> dict[str, Any]:
    pi = arb.pi()
    c = arb(C)
    tstar = pi * c * c / 8
    beta = tstar ** (arb(1) / 3)
    mu_radius = pi / (16 * tstar)
    mu = arb(mu_radius / 2, mu_radius / 2)
    defects = {
        mode: leading_defect(mode, c, beta, mu)
        for mode in range(39_696, 40_095)
        if mode != Q
    }
    rows = []
    coherent_values: list[acb] = []
    leakage_values: list[acb] = []
    for j in range(1, PAIR_COUNT + 1):
        minus_mode = Q - j
        plus_mode = Q + j
        minus = defects[minus_mode]
        plus = defects[plus_mode]
        coherent = (minus + plus.conjugate()) / 2
        leakage = (minus - plus.conjugate()) / 2
        coherent_values.append(coherent)
        leakage_values.append(leakage)
        rows.append({
            "j": j,
            "minus_mode": minus_mode,
            "plus_mode": plus_mode,
            "coherent_coefficient_ball": complex_record(coherent),
            "conjugacy_leakage_coefficient_ball": complex_record(leakage),
            "conjugacy_defect_absolute_ball": (2 * abs(leakage)).str(PRECISION, more=True),
        })

    maximum_defect = max(abs(value).upper() for value in defects.values())
    maximum_pair_conjugacy_defect = max(2 * abs(value).upper() for value in leakage_values)
    leakage_l1 = sum((abs(value).upper() for value in leakage_values), arb(0))
    coherent_variation = abs(coherent_values[-1]).upper() + sum(
        (abs(coherent_values[index] - coherent_values[index + 1]).upper() for index in range(PAIR_COUNT - 1)),
        arb(0),
    )
    edge_modes = [Q + 199, Q + 200]
    edge_l1 = sum((abs(defects[mode]).upper() for mode in edge_modes), arb(0))
    sample_rows = []
    maximum_sample_aggregate = arb(0)
    for index in range(33):
        point_mu = mu_radius * index / 32
        point_t = tstar * (1 - point_mu)
        aggregate = acb(0)
        for mode in sorted(defects):
            coefficient = point_leading_defect(mode, c, beta, point_mu)
            relative_carrier = acb(0, -point_t * (arb(mode) / Q).log()).exp()
            aggregate += coefficient * relative_carrier
        maximum_sample_aggregate = max(maximum_sample_aggregate, abs(aggregate).upper())
        sample_rows.append({
            "height_fraction_index_over_32": index,
            "mu_ball": point_mu.str(PRECISION, more=True),
            "relative_carrier_aggregate_ball": complex_record(aggregate),
            "absolute_ball": abs(aggregate).str(PRECISION, more=True),
        })
    aggregate_panels = 256
    continuous_aggregate_bound = arb(0)
    continuous_maximum_panel = -1
    continuous_maximum_ball = acb(0)
    for index in range(aggregate_panels):
        lower = mu_radius * index / aggregate_panels
        upper = mu_radius * (index + 1) / aggregate_panels
        center = (lower + upper) / 2
        half_width = (upper - lower) / 2
        mu_ball = arb(center, half_width)
        t_ball = tstar * (1 - mu_ball)
        aggregate = acb(0)
        for mode in sorted(defects):
            coefficient = panel_leading_defect(mode, c, beta, center, mu_ball, half_width)
            relative_carrier = acb(0, -t_ball * (arb(mode) / Q).log()).exp()
            aggregate += coefficient * relative_carrier
        if abs(aggregate).upper() > continuous_aggregate_bound:
            continuous_aggregate_bound = abs(aggregate).upper()
            continuous_maximum_panel = index
            continuous_maximum_ball = aggregate
    require(maximum_defect < arb("7.9e-6"), "leading coefficient defect exceeds prior bridge")
    require(maximum_pair_conjugacy_defect < arb("1.6e-5"), "pair conjugacy diagnostic escaped fail-safe cap")
    require(coherent_variation < arb("0.001"), "coherent Abel variation escaped fail-safe cap")
    require(leakage_l1 < arb("0.002"), "leakage l1 diagnostic escaped fail-safe cap")
    require(continuous_aggregate_bound < arb("1.64e-5"), "continuous grouped leading aggregate exceeds 1.64e-5")
    return {
        "height_interval": "t*-pi/16<=t<=t*",
        "mu_interval": "0<=mu<=pi/(16t*)",
        "pair_decomposition": (
            "a_-D+a_+conj(D)=2Re(rD)+2iIm(lD), "
            "r=(a_-+conj(a_+))/2, l=(a_--conj(a_+))/2"
        ),
        "abel_identity": (
            "sum_(j=1)^N r_j z^j=r_N D_N+sum_(j=1)^(N-1)(r_j-r_(j+1))D_j, "
            "D_k=sum_(j=1)^k z^j"
        ),
        "maximum_normalized_leading_defect_ball": maximum_defect.str(PRECISION, more=True),
        "maximum_pair_conjugacy_defect_ball": maximum_pair_conjugacy_defect.str(PRECISION, more=True),
        "pair_conjugacy_leakage_l1_bound_ball": leakage_l1.str(PRECISION, more=True),
        "coherent_Abel_variation_bound_ball": coherent_variation.str(PRECISION, more=True),
        "selector_edge_leading_defect_l1_bound_ball": edge_l1.str(PRECISION, more=True),
        "sampled_relative_carrier_aggregate_maximum_ball": maximum_sample_aggregate.str(PRECISION, more=True),
        "sampled_relative_carrier_aggregate_rows": sample_rows,
        "sampled_aggregate_is_a_proof": False,
        "continuous_relative_carrier_aggregate_panel_count": aggregate_panels,
        "continuous_relative_carrier_aggregate_bound_ball": continuous_aggregate_bound.str(PRECISION, more=True),
        "continuous_relative_carrier_aggregate_maximum_panel": continuous_maximum_panel,
        "continuous_relative_carrier_aggregate_maximum_complex_ball": complex_record(continuous_maximum_ball),
        "continuous_relative_carrier_aggregate_certified": True,
        "pair_rows": rows,
        "edge_rows": [
            {"mode": mode, "leading_defect_ball": complex_record(defects[mode])}
            for mode in edge_modes
        ],
    }


def render_note(artifact: dict[str, Any]) -> str:
    r = artifact["roster"]
    c = artifact["leading_multiplier_certificate"]
    return f"""# Event-ordered pairing after removing event zero

Date: 2026-08-13

Status: exact 398-term selector regrouping and leading-multiplier diagnostic;
not a uniform finite-integral splice

At the completed selector boundary, remove only the already transported
selected summand `m={Q}`, whose detuning label is `-1`.  The remaining
`{r['remaining_term_count']}` selected terms split exactly into
`{r['consecutive_event_pair_count']}` consecutive event pairs

```text
(m_-,m_+)=(q-j,q+j),       q={Q}, 1<=j<=198,          (EP1)
```

and the selector-exchange edge `(40093,40094)`.  For each event pair,

```text
C_- exp(i(4j+1)x)+C_+ exp(-i(4j-1)x)
 =2 exp(ix) Re[C_- exp(i4jx)].                        (EP2)
```

The two edge terms obey

```text
C_+[exp(-i795x)+exp(-i799x)]
 =2 C_+ exp(-i797x) cos(2x).                          (EP3)
```

Therefore the complete remainder of the selected projection is

```text
{r['closed_remaining_projection']}.                  (EP4)
```

The sine quotient in (EP4) always means its original finite 198-term sum.
This derivation removes no complement mode and does not split the endpoint
half-current from the completed remainder.

For the exact leading stationary multiplier, define

```text
a_m=[R_A(s_m,mu) exp(i beta^3 F(s_m,mu))-1]/sqrt(m).
```

Across the top corridor, the measured interval bounds are

```text
max |a_m| < {c['maximum_normalized_leading_defect_ball']},
max |a_(q+j)-conj(a_(q-j))|
 < {c['maximum_pair_conjugacy_defect_ball']},
sum pair leakage < {c['pair_conjugacy_leakage_l1_bound_ball']},
coherent Abel variation < {c['coherent_Abel_variation_bound_ball']},
two selector-edge defects < {c['selector_edge_leading_defect_l1_bound_ball']}.
                                                               (EP5)
```

Factoring out the common classical carrier leaves

```text
sup_(t*-pi/16<=t<=t*)
 |sum_(m!=q) a_m(t) exp[-it log(m/q)]|
 < {c['continuous_relative_carrier_aggregate_bound_ball']} <1.64e-5. (EP6)
```

This is a rigorous 256-panel height enclosure.  The factor `pi/16` is the
exact selector-cell width inherited from the Kummer/Fourier event geometry.

These are leading-carrier coefficients only.  Their exact coherent/leakage
decomposition identifies the finite differences needed by Abel summation,
but (EP5) is not yet multiplied by a rigorously controlled finite-integral
kernel.  It therefore cannot be added to the one-mode transport bound as a
complete splice estimate.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No uniform finite-integral selected-branch theorem, grouped 398-mode
remainder bound, all-corridor continuation, complete `Q_K-T` or `T_upper`,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["completed_projection"].get("passed") is True, "completed projection dependency failed")
    require(dependencies["event0_extension"]["decision"]["event0_fold_chart_covers_complete_top_corridor"] is True, "event-zero dependency drift")
    require(dependencies["stationary_bridge"]["decision"]["exact_leading_stationary_amplitude_ratio_proved"] is True, "stationary bridge dependency drift")

    ctx.dps = PRECISION
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "event_zero_removed_and_remaining_398_selected_terms_paired_in_exact_event_order",
        "passed": True,
        "roster": roster_certificate(),
        "leading_multiplier_certificate": multiplier_certificate(),
        "decision": {
            "event_zero_removed_by_exact_finite_algebra": True,
            "remaining_398_terms_equal_198_event_pairs_plus_one_selector_edge_pair": True,
            "completed_remainder_and_half_current_preserved": True,
            "leading_multiplier_coherent_leakage_decomposition_proved": True,
            "continuous_grouped_leading_multiplier_aggregate_below_1_64e_minus_5": True,
            "finite_integral_kernel_bound_proved": False,
            "grouped_398_mode_splice_proved": False,
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
            "Derive a uniform finite-integral kernel for the coherent and leakage coefficients in (EP5). "
            "Use the exact Abel identity before absolute values and keep the two selector-edge modes inside the completed remainder."
        ),
        "proof_boundary": (
            "Exact finite-roster/event-order algebra and leading stationary-multiplier finite differences on the top corridor only. "
            "No uniform finite-integral kernel bound, grouped 398-mode splice, all-corridor continuation, complete Q_K-T or T_upper theorem, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("paired 398 selected terms as 198 event pairs plus one selector edge pair", flush=True)


if __name__ == "__main__":
    main()
