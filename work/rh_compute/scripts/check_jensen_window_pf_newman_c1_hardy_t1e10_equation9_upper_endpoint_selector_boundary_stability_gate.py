#!/usr/bin/env python3
"""Independently check upper-endpoint stability at the selector boundary."""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, ctx
import mpmath as mp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_upper_endpoint_selector_boundary_stability_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def point_recurrence(height: mp.mpf, effective_start: int) -> tuple[int, list[tuple[int, int, int]]]:
    p = mp.pi
    et = mp.mpf("0.005")
    a = mp.sqrt(8 * height / p)
    fifth = mp.mpf(float(0.2).as_integer_ratio()[0]) / float(0.2).as_integer_ratio()[1]
    ib = 2 * int(mp.floor(height**fifth))
    x = mp.exp((fifth - mp.log(et) / mp.log(height)) / 4)
    endpoint = effective_start + ib
    odd_length = 0
    gsum_count = 0
    rows: list[tuple[int, int, int]] = []
    for block in range(1, 36):
        if endpoint < mp.mpf("1.11") * a:
            gsum_count = int(mp.floor(x**block * ib / 2))
            raw = mp.sqrt(mp.sqrt(et / p) * a) * (endpoint / a - 1) ** mp.mpf("0.625")
        else:
            previous = odd_length
            b = a / (2 * endpoint)
            raw = et**mp.mpf("0.25") * height**mp.mpf("0.25") * (1 / b - b) / (mp.sqrt(p) * 2**mp.mpf("0.875"))
        candidate = int(mp.floor(raw))
        if block > 7 and candidate > mp.mpf("1.5") * previous:
            candidate = int(mp.floor(mp.mpf("1.5") * previous))
        if candidate % 2 == 0:
            candidate -= 1
        odd_length = candidate
        endpoint += 2 * (odd_length + 1) * gsum_count
        rows.append((gsum_count, odd_length, endpoint))
    return endpoint, rows


def main() -> None:
    priority = gate.selector_gate.event_gate.window_gate.cell_gate.ode_gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing upper-endpoint artifact")
    ctx.dps = 130
    mp.mp.dps = 110
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "upper-endpoint artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["pinned_source"]["sha256"] == gate.file_hash(gate.SOURCE), "source hash mismatch")
    require(artifact["dependencies"]["saved_fixture"]["sha256"] == gate.file_hash(gate.FIXTURE), "fixture hash mismatch")
    require(artifact["dependencies"]["selector_strip_gate"]["sha256"] == gate.file_hash(gate.SELECTOR_GATE), "selector hash mismatch")

    saved = artifact["certified_upper_endpoint"]
    fresh = gate.certificate()
    require(saved == fresh, "upper-endpoint certificate mismatch")

    fixture_text = gate.FIXTURE.read_text(encoding="utf-8")
    fixture_endpoints = [int(value) for value in re.findall(r"alpha value at end of this block=\s+([0-9]+)(?:\.0)?", fixture_text)]
    require(len(fixture_endpoints) == 36 and fixture_endpoints[-1] == 5_122_421, "independent fixture parse failed")

    t_saved = mp.mpf(10_000_000_000)
    saved_endpoint, saved_rows = point_recurrence(t_saved, 159577)
    require(saved_endpoint == 5_122_421, "independent saved recurrence failed")
    require([159777] + [row[2] for row in saved_rows] == fixture_endpoints, "independent fixture recurrence mismatch")

    selector_height = mp.pi * mp.mpf(159577) ** 2 / 8
    expected_rows = [(row["gsum_count"], row["odd_length"], row["endpoint"]) for row in saved["boundary_rows"]]
    for offset in (-8000, -1, 0, 1, 200000):
        endpoint, rows = point_recurrence(selector_height + offset, 159579)
        require(endpoint == 5_122_423, f"independent boundary endpoint failed at offset {offset}")
        require(rows == expected_rows, f"independent boundary roster failed at offset {offset}")

    threshold = gate.source_default_real(3.2)
    radius = arb(gate.ROOT_BRACKET_RADIUS)
    for target, center_text in ((159577, gate.LOWER_TRANSITION_ROOT_CENTER), (159579, gate.UPPER_TRANSITION_ROOT_CENTER)):
        center = arb(center_text)
        left = gate.transition_residual(center - radius, target, threshold)
        right = gate.transition_residual(center + radius, target, threshold)
        require(left.lower() > 0 and right.upper() < 0, f"independent transition bracket failed for {target}")

    require(saved["boundary_effective_start_from_left"] == saved["boundary_effective_start_at_and_from_right"] == 159579, "effective-start continuity failed")
    require(saved["boundary_upper_endpoint"] - saved["saved_upper_endpoint"] == 2, "endpoint shift size failed")
    require(saved["boundary_fixed_B_old_selector_term_count"] - saved["boundary_fixed_B_new_selector_term_count"] == 1, "fixed-B term count failed")
    decisions = artifact["decision"]
    require(decisions["upper_endpoint_changes_at_selector_boundary"] is False, "selector-boundary decision corrupted")
    require(decisions["earlier_coupled_start_transition_reassembled"] is False, "coupled-transition boundary corrupted")
    require(decisions["complete_source_selector_splice_proved"] is False, "source-splice boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated upper-endpoint selector stability independently: saved B=5122421, boundary B=5122423, no B jump at t*; priority={priority}")


if __name__ == "__main__":
    main()
