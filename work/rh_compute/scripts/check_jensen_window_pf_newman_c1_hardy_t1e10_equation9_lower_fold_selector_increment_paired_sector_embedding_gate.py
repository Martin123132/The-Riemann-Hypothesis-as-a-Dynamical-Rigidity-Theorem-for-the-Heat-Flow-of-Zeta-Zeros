#!/usr/bin/env python3
"""Independently check the selector-increment paired-sector identity."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
FOLD_LO = 39_695
FOLD_HI = 39_894
AUGMENTED_LO = 39_696
AUGMENTED_HI = 40_094


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    require(RESULT.is_file(), "missing result artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    # Rebuild the signed coefficient identity directly from indicator
    # functions, independently of the builder's state dictionaries.
    actual = Counter()
    for mode in range(FOLD_LO, FOLD_HI + 1):
        actual[mode] += 1
        actual[-mode] += 1
    for mode in range(AUGMENTED_LO, AUGMENTED_HI + 1):
        actual[mode] -= 1
    actual = Counter({mode: value for mode, value in actual.items() if value})

    expected = Counter({FOLD_LO: 1})
    for mode in range(FOLD_LO, FOLD_HI + 1):
        expected[-mode] += 1
    for mode in range(FOLD_HI + 1, AUGMENTED_HI + 1):
        expected[mode] -= 1
    require(actual == expected, "independent Fourier-sector coefficient balance failed")

    # A second exact fixture starts from a new selector state and reconstructs
    # the old state, opposite to the builder's old-minus-strip construction.
    cutoff = 40_103
    strip = {mode: Fraction(5 * mode + 17, 97) for mode in range(-cutoff, cutoff + 1)}
    new = {mode: Fraction(3 * mode - 4, 101) for mode in strip}
    old = {mode: new[mode] + strip[mode] for mode in strip}
    tau = {mode: Fraction(2 * mode + 9, 103) for mode in range(FOLD_LO, FOLD_HI + 1)}

    old_fold = sum((old[r] + old[-r] - tau[r] for r in tau), Fraction())
    new_fold = sum((new[r] + new[-r] - tau[r] for r in tau), Fraction())
    increment = old_fold - new_fold
    signed_strip = sum((strip[r] + strip[-r] for r in tau), Fraction())
    require(increment == signed_strip, "independent target-cancellation check failed")

    raw = sum((strip[m] for m in range(AUGMENTED_LO, AUGMENTED_HI + 1)), Fraction())
    boundary = (
        strip[FOLD_LO]
        + sum((strip[-m] for m in range(FOLD_LO, FOLD_HI + 1)), Fraction())
        - sum((strip[m] for m in range(FOLD_HI + 1, AUGMENTED_HI + 1)), Fraction())
    )
    require(increment == raw + boundary, "independent boundary-block identity failed")

    endpoint = Fraction(29, 109)
    source = endpoint + sum(strip.values(), Fraction())
    outside = source - increment
    raw_remainder = source - raw
    require(increment - raw == raw_remainder - outside, "independent raw completed remainder identity failed")

    decision = artifact["decision"]
    require(decision["fold_increment_equals_signed_strip_pair_block"] is True, "fold increment decision drift")
    require(decision["event_allocation_is_fourier_sector_embedding"] is False, "event allocation was overpromoted")
    require(decision["corrected_selected_weighted_defect_equals_raw_positive_roster"] is False, "weighted defect was identified with raw roster")
    require(decision["common_source_carrier_retained_in_weighted_defect_coordinate"] is True, "common source carrier was dropped")
    require(decision["corrected_selected_weighted_defect_embedded_in_fixed_state_residual"] is False, "weighted defect embedding was overpromoted")
    require(decision["increment_embedding_remainder_bounded"] is False, "unproved remainder bound was promoted")
    require(decision["complete_T_upper_proved"] is False, "unproved endpoint theorem was promoted")
    print("PASS: independent selector-cocycle and paired-sector increment check; remainder remains open")


if __name__ == "__main__":
    main()
