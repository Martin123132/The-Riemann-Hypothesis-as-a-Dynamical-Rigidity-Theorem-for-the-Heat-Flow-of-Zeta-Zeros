#!/usr/bin/env python3
"""Audit the exact mode and height ownership for the lower-fold/ordinary join."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

from flint import arb, ctx


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "portcullis": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json",
    "fresnel_partition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json",
    "turning_handoff": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.json",
    "continuous_height_atlas": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_continuous_height_gate.json",
    "classical_upper": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate.json",
}

PRECISION = 90
C = 159_577
Q = (C - 1) // 4
T0 = 10_000_000_000
EVENT_COUNT = 399
EVENT_RADIUS_DENOMINATOR = 16

SOURCE_START = 622
LOWER_INTERIOR_END = 39_852
SOURCE_END = 39_894
REFLECTED_TRANSITION_END = 39_936
EVENT_SUPPORT_START = 39_695
EVENT_SUPPORT_END = 40_093

TARGET_NORMALIZED = arb("0.000019")
TARGET_PHYSICAL = arb("0.0000086")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


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


def count(lo: int, hi: int) -> int:
    require(lo <= hi, f"empty integer range {lo}..{hi}")
    return hi - lo + 1


def event_mode(index: int) -> int:
    require(0 <= index < EVENT_COUNT, "event index outside atlas")
    return Q - index // 2 if index % 2 == 0 else Q + (index + 1) // 2


def triangular(index: int) -> int:
    return index * (index + 1) // 2


def event_height(index: int, pi: arb) -> arb:
    return pi * C**2 / 8 - pi * (triangular(index) + arb(1) / 8)


def interval_record(value: arb) -> str:
    return value.str(PRECISION, more=True)


def mode_coverage() -> dict[str, Any]:
    source_target = set(range(SOURCE_START, SOURCE_END + 1))
    lower_interior = set(range(SOURCE_START, LOWER_INTERIOR_END + 1))
    lower_transition = set(range(LOWER_INTERIOR_END + 1, SOURCE_END + 1))
    reflected_transition = set(range(SOURCE_END + 1, REFLECTED_TRANSITION_END + 1))
    event_support = set(range(EVENT_SUPPORT_START, EVENT_SUPPORT_END + 1))
    overlap = lower_interior & event_support
    ordinary_core = lower_interior - event_support
    fold_owned_target = source_target & event_support
    event_extra = event_support - source_target

    require(source_target == lower_interior | lower_transition, "source partition is not exhaustive")
    require(lower_interior.isdisjoint(lower_transition), "source partition overlaps")
    require(overlap == set(range(EVENT_SUPPORT_START, LOWER_INTERIOR_END + 1)), "overlap roster drift")
    require(ordinary_core == set(range(SOURCE_START, EVENT_SUPPORT_START)), "ordinary-core roster drift")
    require(fold_owned_target == set(range(EVENT_SUPPORT_START, SOURCE_END + 1)), "fold-owned target drift")
    require(source_target == ordinary_core | fold_owned_target, "proposed target ownership is not exhaustive")
    require(ordinary_core.isdisjoint(fold_owned_target), "proposed target ownership overlaps")
    require(event_extra == set(range(SOURCE_END + 1, EVENT_SUPPORT_END + 1)), "event extra roster drift")
    require(len(reflected_transition) == 42, "reflected transition count drift")

    return {
        "source_classical_target": {"range": [SOURCE_START, SOURCE_END], "count": len(source_target)},
        "fixed_height_lower_interior": {"range": [SOURCE_START, LOWER_INTERIOR_END], "count": len(lower_interior)},
        "fixed_height_lower_transition": {"range": [LOWER_INTERIOR_END + 1, SOURCE_END], "count": len(lower_transition)},
        "fixed_height_reflected_transition": {"range": [SOURCE_END + 1, REFLECTED_TRANSITION_END], "count": len(reflected_transition)},
        "all_event_atlas_support": {"range": [EVENT_SUPPORT_START, EVENT_SUPPORT_END], "count": len(event_support)},
        "ordinary_atlas_validity_overlap": {"range": [EVENT_SUPPORT_START, LOWER_INTERIOR_END], "count": len(overlap)},
        "proposed_disjoint_target_ownership": {
            "ordinary_core": {"range": [SOURCE_START, EVENT_SUPPORT_START - 1], "count": len(ordinary_core)},
            "fold_owned_target": {"range": [EVENT_SUPPORT_START, SOURCE_END], "count": len(fold_owned_target)},
            "disjoint": True,
            "exhaustive_for_source_target": True,
        },
        "atlas_modes_excluded_from_source_target": {"range": [SOURCE_END + 1, EVENT_SUPPORT_END], "count": len(event_extra)},
        "boundary_ledger": {
            "B_tail_to_lower_interior": {"left_mode": 621, "right_mode": 622},
            "ordinary_core_to_fold_owned_overlap": {"left_mode": 39_694, "right_mode": 39_695},
            "lower_interior_to_endpoint_transition": {"left_mode": 39_852, "right_mode": 39_853},
            "lower_to_reflected_transition": {"left_mode": 39_894, "right_mode": 39_895},
            "reflected_transition_to_upper_interior": {"left_mode": 39_936, "right_mode": 39_937},
        },
        "count_identity": "39073 ordinary-core + 158 overlap + 42 lower-transition = 39273 source-target modes",
        "target_projection_guard": "Modes 39895..40093 in the event atlas are exact-Poisson bookkeeping, not terms of the source classical T_upper target.",
    }


def height_coverage(pi: arb) -> dict[str, Any]:
    top = pi * C**2 / 8
    lower = pi * (C - 2) ** 2 / 8
    radius = pi / EVENT_RADIUS_DENOMINATOR
    events = [event_height(index, pi) for index in range(EVENT_COUNT)]
    modes = [event_mode(index) for index in range(EVENT_COUNT)]

    require(set(modes) == set(range(EVENT_SUPPORT_START, EVENT_SUPPORT_END + 1)), "event support is not contiguous")
    require(all(lower < tau < top for tau in events), "event escaped selector cell")
    require(all(events[index] - radius > events[index + 1] + radius for index in range(EVENT_COUNT - 1)), "event buffers overlap")

    event_rows: list[dict[str, Any]] = []
    for index, (mode, tau) in enumerate(zip(modes, events)):
        event_rows.append(
            {
                "event_index": index,
                "mode": mode,
                "center_ball": interval_record(tau),
                "lower_face_ball": interval_record(tau - radius),
                "upper_face_ball": interval_record(tau + radius),
                "ordinary_stationary_side": "t>tau_m",
                "lower_endpoint_transition_side": "t<tau_m",
            }
        )

    corridor_rows: list[dict[str, Any]] = []
    for index in range(EVENT_COUNT + 1):
        upper = top if index == 0 else events[index - 1] - radius
        lower_edge = lower if index == EVENT_COUNT else events[index] + radius
        require(upper > lower_edge, f"empty ordinary-height corridor {index}")
        crossed = modes[:index]
        roster = None if not crossed else [min(crossed), max(crossed)]
        if crossed:
            require(set(crossed) == set(range(roster[0], roster[1] + 1)), f"corridor {index} transition roster is not contiguous")
        corridor_rows.append(
            {
                "corridor_index": index,
                "lower_open_ball": interval_record(lower_edge),
                "upper_open_ball": interval_record(upper),
                "width_ball": interval_record(upper - lower_edge),
                "transition_count": index,
                "transition_roster": roster,
                "ownership": "open corridor; both boundary points belong to adjacent event buffers when present",
            }
        )

    total_event_width = EVENT_COUNT * 2 * radius
    total_corridor_width = sum((arb(row["width_ball"]) for row in corridor_rows), arb(0))
    selector_width = top - lower
    require((total_event_width + total_corridor_width - selector_width).contains(0), "height partition width identity failed")

    spacings = [events[index] - events[index + 1] for index in range(EVENT_COUNT - 1)]
    interior_corridor_widths = [spacing - 2 * radius for spacing in spacings]
    minimum = min(interior_corridor_widths, key=lambda value: value.lower())
    require((minimum - 7 * pi / 8).contains(0), "minimum corridor width identity failed")

    crossed_at_t0 = sum(tau > arb(T0) for tau in events)
    t0_roster = modes[:crossed_at_t0]
    require(crossed_at_t0 == 84, "saved-height event count drift")
    require(set(t0_roster) == set(range(39_853, 39_937)), "saved-height transition roster drift")

    return {
        "selector_cell": {"lower_ball": interval_record(lower), "upper_ball": interval_record(top), "width_ball": interval_record(selector_width)},
        "event_buffer_radius_ball": interval_record(radius),
        "event_buffer_count": EVENT_COUNT,
        "ordinary_height_corridor_count": EVENT_COUNT + 1,
        "minimum_interevent_corridor_width_ball": interval_record(minimum),
        "top_corridor_width_ball": corridor_rows[0]["width_ball"],
        "bottom_corridor_width_ball": corridor_rows[-1]["width_ball"],
        "partition_width_residual_ball": interval_record(total_event_width + total_corridor_width - selector_width),
        "event_buffers_pairwise_disjoint": True,
        "event_buffers_alone_exhaust_selector_height": False,
        "event_buffers_plus_open_corridors_exhaust_selector_height": True,
        "boundary_assignment": "Event buffers own their closed faces; corridor interiors are open, so the ownership partition is disjoint.",
        "saved_height_snapshot": {"height": T0, "crossed_event_count": crossed_at_t0, "transition_roster": [min(t0_roster), max(t0_roster)]},
        "event_rows": event_rows,
        "corridor_rows": corridor_rows,
    }


def render_note(artifact: dict[str, Any]) -> str:
    modes = artifact["mode_coverage"]
    heights = artifact["height_coverage"]
    overlap = modes["ordinary_atlas_validity_overlap"]
    proposed = modes["proposed_disjoint_target_ownership"]
    obligations = artifact["quantitative_frontier"]
    return f"""# Lower-fold/ordinary mode-and-height coverage ledger

Date: 2026-08-13

Status: exact ownership ledger complete; not a proof of the quantitative ordinary-Morse splice

At the saved height, the source classical target is exactly `622..39894`.
The finite-Poisson geometry partitions it as

```text
622..39852       lower ordinary interior (39231 modes),
39853..39894     lower endpoint transition (42 modes).
```

The 399-event atlas supports `39695..40093`.  Its exact overlap with the
lower ordinary interior is `{overlap['range'][0]}..{overlap['range'][1]}`
({overlap['count']} modes).  A disjoint target ownership suitable for a
matched proof is therefore

```text
622..39694       ordinary core ({proposed['ordinary_core']['count']} modes),
39695..39894     fold-owned target ({proposed['fold_owned_target']['count']} modes).
```

The overlap `39695..39852` is retained as a validation region, not counted
twice.  Modes `39895..40093` belong to exact Poisson/selector bookkeeping but
not to the source `T_upper` target.

The 399 closed event buffers have radius `pi/16` and are pairwise disjoint.
They do **not** cover the selector cell by themselves.  Their complement has
{heights['ordinary_height_corridor_count']} open ordinary-height corridors;
the smallest inter-event corridor has exact width `7*pi/8`, enclosed by
`{heights['minimum_interevent_corridor_width_ball']}`.  Assigning every event
face to its event buffer makes the 399 buffers plus 400 corridor interiors a
disjoint and exhaustive height ownership table.

The first missing normalized inequality is

```text
sup_n |K_m(tau_n+pi/16)-S_m(tau_n+pi/16)| <= eps_Morse,
```

where `K_m` is the exact grouped finite-Poisson mode and `S_m` is the ordinary
Morse approximation with its classical carrier.  The corresponding physical
unit is `sqrt(2)/pi` times the normalized value.  The already certified fold
error leaves local prototype headroom
`{obligations['local_normalized_headroom_ball']}` normalized and
`{obligations['local_physical_headroom_ball']}` physical.  These are local
splice diagnostics, not a complete aggregate budget.

The next aggregate obligation is to bound the ordinary remainder on
`622..39852`, preserve the 158-mode overlap before taking absolute values,
and then prove that the disjoint ordinary-core plus fold-owned target equals
the source `T_upper` in the exact phase normalization.

Proof boundary: exact integer ownership, event ordering, and height coverage
only.  No ordinary stationary-phase remainder, quantitative Airy/Morse
overlap, complete `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1

    loaded = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(
        all(item.get("passed") is True for name, item in loaded.items() if name != "continuous_height_atlas"),
        "a dependency gate is not passed",
    )
    continuous_claims = loaded["continuous_height_atlas"].get("claims", {})
    require(
        continuous_claims.get("all_399_exact_event_cells_uniformly_certified") is True
        and continuous_claims.get("normalized_target_proved_on_all_cells") is True
        and continuous_claims.get("physical_target_proved_on_all_cells") is True,
        "continuous-height atlas certification flags are not all true",
    )
    require(loaded["portcullis"]["stationary_classification"]["classical_source_upper_term_count"] == 39_273, "portcullis source count drift")
    require(loaded["turning_handoff"]["certified_event_atlas"]["event_count_inside_selector_cell"] == EVENT_COUNT, "event atlas count drift")
    require(loaded["continuous_height_atlas"]["certificate"]["event_count"] == EVENT_COUNT, "continuous atlas count drift")
    require(loaded["classical_upper"]["scope"]["classical_upper_index_start"] == SOURCE_START, "classical start drift")
    require(loaded["classical_upper"]["scope"]["classical_upper_index_end"] == SOURCE_END, "classical end drift")

    pi = arb.pi()
    modes = mode_coverage()
    heights = height_coverage(pi)

    fold_max_normalized = arb(loaded["continuous_height_atlas"]["certificate"]["maximum_normalized_absolute_upper"])
    fold_max_physical = arb(loaded["continuous_height_atlas"]["certificate"]["maximum_physical_absolute_upper"])
    normalized_headroom = TARGET_NORMALIZED - fold_max_normalized
    physical_headroom = TARGET_PHYSICAL - fold_max_physical
    require(normalized_headroom > 0 and physical_headroom > 0, "existing fold result has no local target headroom")

    artifact = {
        "kind": STEM,
        "status": "exact_lower_fold_ordinary_mode_and_height_ownership_complete_quantitative_splice_open",
        "passed": True,
        "mode_coverage": modes,
        "height_coverage": heights,
        "quantitative_frontier": {
            "exact_mode_notation": "K_m(t)=paper-normalized exact grouped endpoint/Fresnel finite-Poisson mode",
            "fold_notation": "A_m^(4)(t)=beta^-4 corrected lower-fold chart",
            "ordinary_notation": "S_m(t)=ordinary Morse chart with classical m^(-1/2) carrier",
            "first_unproved_normalized_inequality": "sup_n |K_(m_n)(tau_n+pi/16)-S_(m_n)(tau_n+pi/16)| <= epsilon_Morse",
            "physical_conversion": "epsilon_physical=(sqrt(2)/pi)*epsilon_normalized",
            "existing_maximum_fold_normalized_ball": fold_max_normalized.str(PRECISION, more=True),
            "existing_maximum_fold_physical_ball": fold_max_physical.str(PRECISION, more=True),
            "local_normalized_headroom_ball": normalized_headroom.str(PRECISION, more=True),
            "local_physical_headroom_ball": physical_headroom.str(PRECISION, more=True),
            "aggregate_ordinary_remainder": "R_ord(t)=sum_(m=622)^39852 [K_m(t)-m^(-1/2)exp(i(theta(t)-t log m))]",
            "aggregate_target_identity_to_prove": "ordinary-core 622..39694 plus fold-owned 39695..39894 equals source T_upper 622..39894 after exact carrier/conjugate assembly",
            "local_headroom_is_complete_T_upper_budget": False,
        },
        "decision": {
            "mode_ownership_disjoint_and_exhaustive": True,
            "height_ownership_disjoint_and_exhaustive": True,
            "event_cells_alone_cover_selector_height": False,
            "ordinary_fold_validity_overlap_nonempty": True,
            "ordinary_morse_remainder_proved": False,
            "quantitative_splice_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"elapsed_seconds": round(time.time() - started, 3), "workers": 1, "flint_threads": 1, "process_priority": priority},
        "next_obligation": "Derive an explicit ordinary-Morse mode remainder in the exact grouped normalization, certify it first at the stationary faces tau_n+pi/16 and on the 158-mode saved-height overlap, then sum it on 622..39852 without termwise loss that exceeds the source budget.",
        "proof_boundary": "Exact mode-index and selector-height ownership only. No ordinary stationary-phase remainder, quantitative Airy/Morse splice, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }

    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built lower-fold/ordinary coverage ledger: "
        f"target={modes['source_classical_target']['count']}, overlap={modes['ordinary_atlas_validity_overlap']['count']}, "
        f"event-buffers={heights['event_buffer_count']}, corridors={heights['ordinary_height_corridor_count']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
