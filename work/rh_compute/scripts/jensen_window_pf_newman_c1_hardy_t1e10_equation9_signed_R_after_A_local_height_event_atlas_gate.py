#!/usr/bin/env python3
"""Certify the maximal same-roster height cell around t=1e10."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_local_height_event_atlas_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)


def result_path(stem: str) -> Path:
    return ROOT / "work" / "rh_compute" / "results" / f"{stem}.json"


DEPENDENCIES = {
    "saved_closure": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_one_sided_closure_gate"
    ),
    "finite_mode_roster": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_branch_Fresnel_Dirichlet_connection_and_mode_roster_gate"
    ),
    "turning_event_atlas": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate"
    ),
    "ordinary_complement": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate"
    ),
    "closed_752_layer": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_label_exterior_transition_layer_gate"
    ),
    "source_schedule": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate"
    ),
}

PRECISION_BITS = 384
HEIGHT = 10_000_000_000
A = 159_577
B = 5_122_421
SOURCE_COUNT = 2_481_423
N_MINUS = 79_788
N_PLUS = 2_561_211

ORDINARY_START = 622
ORDINARY_END = 39_852
TRANSITION_START = 39_853
TRANSITION_END = 39_936
CLASSICAL_END = 39_894

# Exact rational boundaries.  The two order-16 points are deliberately distinct.
L_NUM, L_DEN = 1_243, 2
LOWER_A_CELL_NUM, LOWER_A_CELL_DEN = 79_705, 2
UPPER_A_CELL_NUM, UPPER_A_CELL_DEN = 79_873, 2
CORE_SPLIT_NUM, CORE_SPLIT_DEN = 9_945, 16   # 621 + 9/16
FAR_SPLIT_NUM, FAR_SPLIT_DEN = 9_947, 16    # 621 + 11/16

Q_MINUS_NUM, Q_MINUS_DEN = 159_577, 2
Q_PLUS_NUM, Q_PLUS_DEN = 5_122_423, 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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


def ball_record(value: arb, digits: int = 90) -> dict[str, str]:
    return {"ball": value.str(digits, more=True)}


def rational(num: int, den: int) -> arb:
    return arb(num) / den


def event_height(alpha: int, mode: int) -> arb:
    return arb.pi() * mode * (alpha - 2 * mode)


def saddle_roots(alpha: int, height: arb) -> tuple[arb, arb]:
    radical = (arb(alpha) ** 2 - 8 * height / arb.pi()).sqrt()
    return (arb(alpha) - radical) / 4, (arb(alpha) + radical) / 4


def q_value(height: arb, y: arb) -> arb:
    return y + height / (2 * arb.pi() * y)


def q_prime(height: arb, y: arb) -> arb:
    return 1 - height / (2 * arb.pi() * y * y)


def eight_width_margin(height: arb, y: arb, q: arb) -> arb:
    detuning = q - q_value(height, y)
    return detuning * detuning - 64 * abs(q_prime(height, y))


def closed_ball(lower: arb, upper: arb) -> arb:
    require(lower < upper, "invalid interval endpoints")
    midpoint = (lower + upper) / 2
    radius = ((upper - lower) / 2).upper()
    return arb(midpoint, radius)


def dependency_guards(deps: dict[str, dict[str, Any]], lower: arb, upper: arb) -> dict[str, Any]:
    require(deps["saved_closure"]["passed"], "saved closure dependency failed")

    roster = deps["finite_mode_roster"]["saddle_mode_roster"]
    require(roster["alpha_A"] == A and roster["alpha_B"] == B, "finite roster endpoint drift")
    require(roster["ordinary_integer_roster"] == [ORDINARY_START, ORDINARY_END], "ordinary roster drift")
    require(roster["transition_integer_roster"] == [TRANSITION_START, TRANSITION_END], "transition roster drift")
    require(roster["transition_lower_half"] == [TRANSITION_START, CLASSICAL_END], "lower half drift")
    require(roster["transition_upper_half"] == [CLASSICAL_END + 1, TRANSITION_END], "upper half drift")

    atlas = deps["turning_event_atlas"]["certified_event_atlas"]
    prior_lower = arb(atlas["next_downward_event_height_ball"])
    prior_upper = arb(atlas["next_upward_event_height_ball"])
    require(lower.overlaps(prior_lower), "lower event misses prior turning atlas")
    require(upper.overlaps(prior_upper), "upper event misses prior turning atlas")
    require(atlas["saved_transition_roster"] == [TRANSITION_START, TRANSITION_END], "prior atlas roster drift")

    ordinary = deps["ordinary_complement"]["certificate"]
    require(ordinary["upper_stationary_label_count"] == 230, "stationary complement count drift")
    require(
        ordinary["upper_stationary_q_range"] == ["2561211.5", "2561440.5"],
        "stationary complement labels drift",
    )

    layer = deps["closed_752_layer"]["certificate"]
    require(layer["stationary_label_count"] == 230, "752 layer stationary count drift")
    require(layer["transition_label_count"] == 752, "752 layer count drift")
    require(layer["root_of_unity_order_at_split"] == 16, "752 layer root order drift")

    schedule = deps["source_schedule"]["certificate"]
    require(schedule["pre_effective_start"] == A, "source pre-start drift")
    require(schedule["pre_endpoint"] == B, "source pre-endpoint drift")
    require(schedule["pre_roster_count"] == SOURCE_COUNT, "source pre-count drift")
    source_transition = arb(schedule["transition_entry_height_ball"])
    require(source_transition > upper, "audited source scheduling transition entered local cell")

    return {
        "prior_turning_atlas_lower_overlap": True,
        "prior_turning_atlas_upper_overlap": True,
        "saved_mode_roster_reproduced": True,
        "saved_complement_census_reproduced": True,
        "saved_752_layer_combinatorics_reproduced": True,
        "audited_next_source_transition_ball": ball_record(source_transition),
        "audited_next_source_transition_margin_above_cell_ball": ball_record(source_transition - upper),
        "source_scope": (
            "The theorem chart fixes A and B. The audited next implementation scheduling transition lies "
            "above this cell; this gate does not promote the compiled source approximation to an exact theorem."
        ),
    }


def build_certificate(deps: dict[str, dict[str, Any]]) -> dict[str, Any]:
    ctx.prec = PRECISION_BITS
    pi = arb.pi()
    t0 = arb(HEIGHT)
    lower = event_height(A, ORDINARY_END)
    upper = event_height(A, TRANSITION_END)
    require(lower < t0 < upper, "saved height is outside primary cell")
    t_cell = closed_ball(lower, upper)

    L = rational(L_NUM, L_DEN)
    a = rational(LOWER_A_CELL_NUM, LOWER_A_CELL_DEN)
    U = rational(UPPER_A_CELL_NUM, UPPER_A_CELL_DEN)
    core_split = rational(CORE_SPLIT_NUM, CORE_SPLIT_DEN)
    far_split = rational(FAR_SPLIT_NUM, FAR_SPLIT_DEN)
    q_minus = rational(Q_MINUS_NUM, Q_MINUS_DEN)
    q_plus = rational(Q_PLUS_NUM, Q_PLUS_DEN)

    require(A == 2 * (ORDINARY_END + U), "lower-wall complementary root identity drift")
    require(A == 2 * (TRANSITION_END + a), "upper-wall complementary root identity drift")
    require((lower - 2 * pi * U * ORDINARY_END).contains(0), "lower-wall Q(U) formula drift")
    require((upper - 2 * pi * a * TRANSITION_END).contains(0), "upper-wall Q(a) formula drift")

    lower_A_at_lower, upper_A_at_lower = saddle_roots(A, lower)
    lower_A_at_upper, upper_A_at_upper = saddle_roots(A, upper)
    require((lower_A_at_lower - ORDINARY_END).contains(0), "lower A root misses lower wall")
    require((upper_A_at_lower - U).contains(0), "upper A root misses lower-wall half cell")
    require((lower_A_at_upper - a).contains(0), "lower A root misses upper-wall half cell")
    require((upper_A_at_upper - TRANSITION_END).contains(0), "upper A root misses upper wall")

    transition_margin = arb(A) - (8 * upper / pi).sqrt()
    require(transition_margin > 0, "A ceased to be above the saddle transition")

    center_lower = (lower / (2 * pi)).sqrt()
    center_upper = (upper / (2 * pi)).sqrt()
    require(center_lower > CLASSICAL_END, "classical floor crossed lower integer")
    require(center_upper < CLASSICAL_END + 1, "classical floor crossed upper integer")

    B_lower_at_lower, B_upper_at_lower = saddle_roots(B, lower)
    B_lower_at_upper, B_upper_at_upper = saddle_roots(B, upper)
    require(arb(621) < B_lower_at_lower < B_lower_at_upper < arb(622), "B lower saddle floor drift")
    require(
        arb(2_560_588) < B_upper_at_upper < B_upper_at_lower < arb(2_560_589),
        "B upper saddle floor drift",
    )

    qL_lower = q_value(lower, L)
    qL_upper = q_value(upper, L)
    qa_lower = q_value(lower, a)
    qa_upper = q_value(upper, a)
    qU_lower = q_value(lower, U)
    qU_upper = q_value(upper, U)

    existing_upper_margin = qL_lower - (q_plus + 229)
    next_upper_gap = q_plus + 230 - qL_upper
    lower_complement_gap = qa_lower - (q_minus - 1)
    core_saddle_exit_gap = q_plus - q_value(upper, core_split)
    require(existing_upper_margin > 0, "upper stationary label 230 left the roster")
    require(next_upper_gap > 0, "upper stationary label 231 entered the roster")
    require(lower_complement_gap > 0, "lower complementary label became stationary")
    require(core_saddle_exit_gap > 0, "first upper saddle escaped the rational core")
    require((qU_lower - q_minus).contains(0), "lower wall Q(U)=q_minus identity failed")
    require((q_minus - qa_upper).contains(0), "upper wall Q(a)=q_minus identity failed")
    require(qU_upper > q_minus and qa_lower < q_minus, "A saddle orientation drift")

    below_center_margin = lower / (2 * pi * a * a) - 1
    require(below_center_margin > 0, "ordinary interval reached the Q turning point")

    finite_margin = eight_width_margin(t_cell, far_split, q_plus)
    remote_margin = eight_width_margin(t_cell, L, q_plus + 752)
    previous_grid_margin = eight_width_margin(t_cell, L + rational(2, 16), q_plus)
    previous_count_margin = eight_width_margin(t_cell, L, q_plus + 736)
    require(finite_margin.lower() > 0, "finite 752 launch lost eight-width separation")
    require(remote_margin.lower() > 0, "remote 752 launch lost eight-width separation")
    require(previous_grid_margin.upper() < 0, "previous 1/16 grid point became admissible")
    require(previous_count_margin.upper() < 0, "previous multiple of 16 became admissible")

    require(math.gcd(9, 16) == 1 and math.gcd(11, 16) == 1, "root order arithmetic drift")
    require(240 == 15 * 16 and 752 == 47 * 16, "closed roster factorization drift")
    require(230 <= 240 < 752, "closed roster ordering drift")

    # Section 11.492 is a diagnostic donor branch, but its integer route split is stable too.
    sixth_root = (t_cell.log() / 6).exp()
    lambda_star = (sixth_root + (sixth_root * sixth_root - 4).sqrt()) / 2
    y_star = (t_cell / (2 * pi)).sqrt() / lambda_star
    require(y_star > rational(1719, 2), "diagnostic route split crossed 859.5")
    require(y_star < rational(1721, 2), "diagnostic route split crossed 860.5")

    dependency_checks = dependency_guards(deps, lower, upper)

    return {
        "saved_height": HEIGHT,
        "precision_bits": PRECISION_BITS,
        "primary_open_cell": {
            "lower_exact": "pi*39852*(159577-2*39852)",
            "upper_exact": "pi*39936*(159577-2*39936)",
            "lower_ball": ball_record(lower),
            "upper_ball": ball_record(upper),
            "saved_height_distance_above_lower_ball": ball_record(t0 - lower),
            "saved_height_distance_below_upper_ball": ball_record(upper - t0),
            "width_ball": ball_record(upper - lower),
            "maximality": (
                "Maximal open interval containing 1e10 on which the fixed-A transition roster "
                "39853..39936 and its ordinary/transition cell ownership remain unchanged."
            ),
            "lower_wall": {
                "crossing_mode": ORDINARY_END,
                "A_lower_root_ball": ball_record(lower_A_at_lower),
                "A_upper_root_ball": ball_record(upper_A_at_lower),
                "simultaneous_exact_relations": ["y_-(A)=39852", "y_+(A)=39936.5=U"],
            },
            "upper_wall": {
                "crossing_mode": TRANSITION_END,
                "A_lower_root_ball": ball_record(lower_A_at_upper),
                "A_upper_root_ball": ball_record(upper_A_at_upper),
                "simultaneous_exact_relations": ["y_-(A)=39852.5=a", "y_+(A)=39936"],
            },
        },
        "height_independent_inventory": {
            "fixed_odd_window": {
                "A": A,
                "B": B,
                "M": SOURCE_COUNT,
                "n_minus": N_MINUS,
                "n_plus": N_PLUS,
                "q_minus": "159577/2",
                "q_plus": "5122423/2",
            },
            "fixed_cell_boundaries": {
                "L": "1243/2=621.5",
                "a": "79705/2=39852.5",
                "U": "79873/2=39936.5",
                "core_split": "9945/16=621+9/16",
                "far_split": "9947/16=621+11/16",
            },
            "root_of_unity_closures": {
                "alternating_order": 2,
                "core_order": 16,
                "core_closed_count": 240,
                "far_closed_count": 752,
                "proof": (
                    "240=15*16 and 752=47*16; gcd(9,16)=gcd(11,16)=1. "
                    "These endpoint cancellations are exact and independent of height."
                ),
            },
        },
        "uniform_combinatorial_guards": {
            "A_remains_above_transition_margin_at_upper_ball": ball_record(transition_margin),
            "classical_center_at_lower_ball": ball_record(center_lower),
            "classical_center_at_upper_ball": ball_record(center_upper),
            "classical_floor": CLASSICAL_END,
            "B_lower_saddle_at_lower_ball": ball_record(B_lower_at_lower),
            "B_lower_saddle_at_upper_ball": ball_record(B_lower_at_upper),
            "B_upper_saddle_at_lower_ball": ball_record(B_upper_at_lower),
            "B_upper_saddle_at_upper_ball": ball_record(B_upper_at_upper),
            "B_lower_floor": 621,
            "B_upper_floor": 2_560_588,
            "ordinary_roster": [ORDINARY_START, ORDINARY_END],
            "ordinary_count": ORDINARY_END - ORDINARY_START + 1,
            "transition_roster": [TRANSITION_START, TRANSITION_END],
            "transition_count": TRANSITION_END - TRANSITION_START + 1,
            "transition_halves": [[TRANSITION_START, CLASSICAL_END], [CLASSICAL_END + 1, TRANSITION_END]],
            "transition_half_counts": [42, 42],
            "upper_complement_stationary_count": 230,
            "upper_complement_q_range": ["2561211.5", "2561440.5"],
            "existing_230th_label_margin_at_lower_ball": ball_record(existing_upper_margin),
            "next_231st_label_gap_at_upper_ball": ball_record(next_upper_gap),
            "lower_complement_nonstationary_gap_at_lower_ball": ball_record(lower_complement_gap),
            "first_upper_saddle_before_core_split_gap_at_upper_ball": ball_record(core_saddle_exit_gap),
            "Q_at_L_lower_ball": ball_record(qL_lower),
            "Q_at_L_upper_ball": ball_record(qL_upper),
            "Q_at_a_lower_ball": ball_record(qa_lower),
            "Q_at_a_upper_ball": ball_record(qa_upper),
            "Q_at_U_lower_ball": ball_record(qU_lower),
            "Q_at_U_upper_ball": ball_record(qU_upper),
            "physical_interval_below_sqrt_p_margin_ball": ball_record(below_center_margin),
        },
        "uniform_route_guards": {
            "closed_cell_ball": ball_record(t_cell),
            "finite_752_eight_width_margin_ball": ball_record(finite_margin),
            "remote_752_eight_width_margin_ball": ball_record(remote_margin),
            "previous_grid_eight_width_margin_ball": ball_record(previous_grid_margin),
            "previous_count_eight_width_margin_ball": ball_record(previous_count_margin),
            "selected_launches_remain_admissible": True,
            "immediate_predecessors_remain_inadmissible": True,
            "diagnostic_unequal_truncation_y_star_ball": ball_record(y_star),
            "diagnostic_route_split": {"endpoint_cells": [622, 860], "contracting_core": [861, 39936]},
            "diagnostic_scope": "The Section 11.492 donor split is stable but is not used by the final signed closure route.",
        },
        "dependency_checks": dependency_checks,
        "identity_scope": {
            "uniform_exact_algebra": (
                "With A and B fixed, the finite RSI prefix, branch correction, finite Dirichlet source, "
                "cell covariance, pole-safe join, entire-kernel finite difference, and complementary-lattice "
                "identities remain exact throughout the open cell."
            ),
            "uniform_numerical_enclosures": False,
            "saved_height_sign_transported": False,
            "crossing_handoff_needed_at_lower_wall": True,
            "crossing_handoff_needed_at_upper_wall": True,
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    cell = c["primary_open_cell"]
    guards = c["uniform_combinatorial_guards"]
    route = c["uniform_route_guards"]
    return f"""# Signed-closure local height event atlas

Date: 2026-08-28

Status: maximal same-roster cell and all active discrete guards certified;
height-uniform packet bounds and sign transport remain open

Keep the exact odd window fixed at

```text
A=159577, B=5122421, M=2481423.
```

For a fixed odd label `alpha`, a saddle reaches mode `m` at

```text
tau_m(alpha)=pi*m*(alpha-2m).
```

The nearest two active events around `t0=10^10` are therefore

```text
t_-=pi*39852*(159577-2*39852)
   ={cell['lower_ball']['ball']},

t_+=pi*39936*(159577-2*39936)
   ={cell['upper_ball']['ball']}.
```

Hence the maximal open interval containing the saved height on which the
same `A` transition roster is valid is

```text
t_- < t < t_+,
t0-t_-={cell['saved_height_distance_above_lower_ball']['ball']},
t_+-t0={cell['saved_height_distance_below_upper_ball']['ball']}.
```

The walls are not decimal guesses.  Exact complementary-root identities give

```text
t=t_-: y_-(A)=39852,   y_+(A)=39936.5=U,
t=t_+: y_-(A)=39852.5=a, y_+(A)=39936.
```

Inside the open cell the complete active ownership data stay fixed:

```text
ordinary modes:       622..39852 ({guards['ordinary_count']} modes),
transition modes:     39853..39936 ({guards['transition_count']} modes),
classical split:      39853..39894 | 39895..39936,
B lower/upper floors: 621 and 2560588,
upper complement:     230 stationary labels,
                      2561211.5..2561440.5.
```

The rational boundaries are proof choices, not height events:

```text
L=621.5,
B_core=621+9/16=621.5625,
C_far=621+11/16=621.6875.
```

Their exact Fourier orders remain 2 and 16 independently of `t`.
`240=15*16` and `752=47*16`, so both grouped endpoint cancellations remain
exact.  On the whole closed event cell, interval arithmetic also proves

```text
finite 752-launch margin  > 0: {route['finite_752_eight_width_margin_ball']['ball']},
remote 752-launch margin  > 0: {route['remote_752_eight_width_margin_ball']['ball']},
previous grid margin      < 0: {route['previous_grid_eight_width_margin_ball']['ball']},
previous count margin     < 0: {route['previous_count_eight_width_margin_ball']['ball']}.
```

Thus the chosen far launches remain admissible and their immediate tested
predecessors remain inadmissible everywhere in the cell.  The dormant
unequal-truncation donor split also stays `622..860 | 861..39936`, but it is
not part of the final signed route.

What this certifies is the discrete and exact-algebra layer.  With `A,B`
fixed, the finite RSI identity, branch correction, Dirichlet-cell partition,
pole-safe contour join, entire-kernel finite difference, and complementary
half-lattice identities all retain the same form in this cell.

What it does **not** certify is just as important: the complex balls for
`V_L`, `C_U`, `T_upper`, `T_lower`, `O_join`, `T_A^join`, `J_Z`, the signed B
trace, or the outer allowance have not yet been transported in height.  The
saved inequality `Q_K-T<0` therefore remains a theorem only at `t=10^10`.
Crossing either wall requires the already-known exact one-mode ownership
handoff plus matched interval estimates; this atlas does not silently cross
it.

Next obligation: derive interval-valued transport bounds for the joined
complex packets, signed B trace, and outer term on a first certified subcell,
retaining the exact cancellations before every norm.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing required path: {path}")
    deps = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    certificate = build_certificate(deps)
    artifact = {
        "kind": STEM,
        "status": "maximal_same_roster_local_height_event_cell_and_active_discrete_guards_certified",
        "passed": True,
        "certificate": certificate,
        "decision": {
            "maximal_same_roster_cell_certified": True,
            "all_active_discrete_guards_stable_on_cell": True,
            "root_of_unity_closures_height_independent": True,
            "uniform_packet_enclosures_certified": False,
            "saved_sign_transported_off_height": False,
            "all_height_theorem": False,
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
        "resource_mode": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_obligation": (
            "Transport the joined complex packets, signed B trace, and outer allowance on a first interval "
            "subcell while retaining exact cancellation; then add exact crossing handoffs at both event walls."
        ),
        "proof_boundary": (
            "Exact local event atlas and uniform combinatorial guards only. No height-uniform packet balls, "
            "off-height sign theorem, all-height theorem, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print(
        "built local height event atlas: same-roster cell "
        "(-206.835884..., +57.057899...), all active discrete guards stable"
    )


if __name__ == "__main__":
    main()
