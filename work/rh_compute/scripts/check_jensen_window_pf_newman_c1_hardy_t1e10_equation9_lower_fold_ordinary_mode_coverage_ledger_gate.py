#!/usr/bin/env python3
"""Independent checks for the lower-fold/ordinary ownership ledger."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from flint import arb, ctx


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()

C = 159_577
Q = (C - 1) // 4
EVENT_COUNT = 399


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def event_mode(index: int) -> int:
    return Q - index // 2 if index % 2 == 0 else Q + (index + 1) // 2


def event_height(index: int, pi: arb) -> arb:
    triangular = index * (index + 1) // 2
    return pi * C**2 / 8 - pi * (triangular + arb(1) / 8)


def main() -> int:
    ctx.dps = 120
    ctx.threads = 1
    require(RESULT.is_file(), "coverage result is missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "coverage result is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency {path}")
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    modes = artifact["mode_coverage"]
    require(modes["source_classical_target"] == {"range": [622, 39894], "count": 39273}, "source target drift")
    require(modes["fixed_height_lower_interior"] == {"range": [622, 39852], "count": 39231}, "lower interior drift")
    require(modes["ordinary_atlas_validity_overlap"] == {"range": [39695, 39852], "count": 158}, "overlap drift")
    require(modes["proposed_disjoint_target_ownership"]["ordinary_core"]["count"] == 39073, "ordinary core count drift")
    require(modes["proposed_disjoint_target_ownership"]["fold_owned_target"]["count"] == 200, "fold target count drift")

    height = artifact["height_coverage"]
    require(len(height["event_rows"]) == EVENT_COUNT, "event row count drift")
    require(len(height["corridor_rows"]) == EVENT_COUNT + 1, "corridor row count drift")
    require(height["event_buffers_alone_exhaust_selector_height"] is False, "event-only coverage overclaim")
    require(height["event_buffers_plus_open_corridors_exhaust_selector_height"] is True, "height exhaustion missing")

    pi = arb.pi()
    radius = pi / 16
    events = [event_height(index, pi) for index in range(EVENT_COUNT)]
    modes_independent = [event_mode(index) for index in range(EVENT_COUNT)]
    require(set(modes_independent) == set(range(39695, 40094)), "independent event support drift")
    require(all(events[index] - radius > events[index + 1] + radius for index in range(EVENT_COUNT - 1)), "independent event overlap")
    minimum_corridor = min(events[index] - events[index + 1] - 2 * radius for index in range(EVENT_COUNT - 1))
    require((minimum_corridor - 7 * pi / 8).contains(0), "independent minimum corridor identity failed")
    require(arb(height["minimum_interevent_corridor_width_ball"]).overlaps(minimum_corridor), "stored minimum corridor does not overlap")

    top = pi * C**2 / 8
    lower = pi * (C - 2) ** 2 / 8
    widths = sum((arb(row["width_ball"]) for row in height["corridor_rows"]), arb(0))
    residual = EVENT_COUNT * 2 * radius + widths - (top - lower)
    require(residual.contains(0), "independent height width identity failed")
    require(arb(height["partition_width_residual_ball"]).contains(0), "stored width residual excludes zero")

    crossed = [mode for mode, tau in zip(modes_independent, events) if tau > arb(10_000_000_000)]
    require(len(crossed) == 84 and set(crossed) == set(range(39853, 39937)), "saved-height roster drift")

    frontier = artifact["quantitative_frontier"]
    factor = arb(2).sqrt() / pi
    fold_normalized = arb(frontier["existing_maximum_fold_normalized_ball"])
    fold_physical = arb(frontier["existing_maximum_fold_physical_ball"])
    norm_headroom = arb(frontier["local_normalized_headroom_ball"])
    phys_headroom = arb(frontier["local_physical_headroom_ball"])
    require((factor * fold_normalized).overlaps(fold_physical), "normalized/physical fold conversion drift")
    require(norm_headroom > 0 and phys_headroom > 0, "local splice headroom is not positive")
    require(phys_headroom > factor * norm_headroom, "independent physical target is unexpectedly tighter")
    require(artifact["decision"]["ordinary_morse_remainder_proved"] is False, "ordinary remainder overclaim")
    require(artifact["decision"]["complete_T_upper_proved"] is False, "T_upper overclaim")

    print("checked lower-fold/ordinary coverage ledger: 39273 target modes, 158-mode overlap, 399+400 height partition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
