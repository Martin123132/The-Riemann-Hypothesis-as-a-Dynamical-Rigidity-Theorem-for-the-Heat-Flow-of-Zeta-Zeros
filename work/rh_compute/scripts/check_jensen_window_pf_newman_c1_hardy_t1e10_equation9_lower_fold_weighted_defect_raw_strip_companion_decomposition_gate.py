#!/usr/bin/env python3
"""Independently check the weighted-defect/raw-strip companion identity."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
Q = 39_894
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

    parity_integer = (C * C - 1) // 8
    require((C * C - 1) % 8 == 0 and parity_integer % 2 == 0, "independent carrier parity failed")
    require(artifact["carrier_parity"]["integer"] == parity_integer, "saved carrier parity drift")

    # Rebuild the identity with data unrelated to the builder fixture and in
    # the opposite direction: construct corrected and companion coefficients,
    # recover every raw strip coefficient, then assemble the paired increment.
    kappa = Fraction(-29, 31)
    exact_profile = {
        mode: Fraction(11 * (mode - Q) - 5, 401 + 2 * abs(mode - Q))
        for mode in range(AUGMENTED_LO, AUGMENTED_HI + 1)
    }
    opposite = {
        mode: Fraction(13 * (mode - Q) + 3, 419 + 3 * abs(mode - Q))
        for mode in exact_profile
    }
    weights = {
        mode: (Fraction() if mode == Q else Fraction(5 * (mode - Q) - 1, 200_003 + abs(mode - Q)))
        for mode in exact_profile
    }
    corrected = {
        mode: kappa * weights[mode] * (exact_profile[mode] - opposite[mode])
        for mode in exact_profile
    }
    companion = {
        mode: kappa * ((1 - weights[mode]) * exact_profile[mode] + weights[mode] * opposite[mode])
        for mode in exact_profile
    }
    recovered = {mode: corrected[mode] + companion[mode] for mode in exact_profile}
    expected_raw = {mode: kappa * exact_profile[mode] for mode in exact_profile}
    require(recovered == expected_raw, "independent modewise recovery failed")
    require(corrected[Q] == 0 and companion[Q] == expected_raw[Q], "independent event-zero recovery failed")

    strip = {
        mode: Fraction(17 * mode + 7, 431 + abs(mode))
        for mode in range(-AUGMENTED_HI, AUGMENTED_HI + 1)
    }
    strip.update(recovered)
    raw_sum = sum(recovered.values(), Fraction())
    corrected_sum = sum(corrected.values(), Fraction())
    companion_sum = sum(companion.values(), Fraction())
    corrected_inner = sum((corrected[m] for m in range(AUGMENTED_LO, Q + 1)), Fraction())
    corrected_outer = sum((corrected[m] for m in range(Q + 1, AUGMENTED_HI + 1)), Fraction())
    sector = (
        strip[FOLD_LO]
        + sum((strip[-mode] for mode in range(FOLD_LO, FOLD_HI + 1)), Fraction())
        - sum((strip[mode] for mode in range(FOLD_HI + 1, AUGMENTED_HI + 1)), Fraction())
    )
    fold_increment = sum(
        (strip[mode] + strip[-mode] for mode in range(FOLD_LO, FOLD_HI + 1)),
        Fraction(),
    )
    require(raw_sum == corrected_sum + companion_sum, "independent roster split failed")
    require(fold_increment == corrected_sum + companion_sum + sector, "independent increment split failed")
    paired_companion = (
        strip[FOLD_LO]
        + strip[-FOLD_LO]
        + sum(
            (companion[mode] + strip[-mode] for mode in range(AUGMENTED_LO, Q + 1)),
            Fraction(),
        )
    )
    require(corrected_sum == corrected_inner + corrected_outer, "independent corrected branch split failed")
    require(fold_increment == corrected_inner + paired_companion, "independent collapsed increment failed")
    require(companion_sum + sector == paired_companion - corrected_outer, "independent outer cancellation failed")

    outside = Fraction(-37, 443)
    removed_source = fold_increment + outside
    require(removed_source == corrected_sum + companion_sum + sector + outside, "independent completion failed")

    decision = artifact["decision"]
    require(decision["raw_strip_mode_is_carrier_restored_exact_common_profile_coefficient"] is True, "raw mode identity drift")
    require(decision["modewise_weighted_defect_plus_companion_equals_raw_strip"] is True, "modewise split drift")
    require(decision["weighted_corrected_defect_inserted_into_selector_increment"] is True, "increment insertion drift")
    require(decision["common_source_carrier_retained"] is True, "source carrier was dropped")
    require(decision["event_zero_retained_in_coefficient_companion"] is True, "event zero was dropped")
    require(decision["outer_corrected_block_cancels_from_collapsed_fold_increment"] is True, "outer corrected block did not cancel")
    require(decision["collapsed_fold_increment_equals_inner_corrected_plus_paired_companion"] is True, "collapsed increment drift")
    require(decision["total_corrected_projection_determines_inner_corrected_sign"] is False, "total sign was transferred to inner block")
    require(decision["inner_corrected_projection_bounded"] is False, "inner projection was overpromoted")
    require(decision["coefficient_companion_plus_raw_sector_remainder_bounded"] is False, "joint remainder was overpromoted")
    require(decision["signed_fold_increment_proved"] is False, "signed increment was overpromoted")
    require(decision["complete_T_upper_proved"] is False, "T_upper was overpromoted")
    print("PASS: independent weighted-defect/raw-strip companion and selector-increment check; joint remainder remains open")


if __name__ == "__main__":
    main()
