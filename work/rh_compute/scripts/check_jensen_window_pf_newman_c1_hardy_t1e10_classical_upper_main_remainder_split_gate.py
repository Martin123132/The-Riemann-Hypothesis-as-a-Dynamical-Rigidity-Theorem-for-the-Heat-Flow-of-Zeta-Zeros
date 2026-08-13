#!/usr/bin/env python3
"""Independently validate the classical upper-main residual split."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
SCRIPT_ROOT = Path(__file__).resolve().parent
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate as gate
import jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate as component_split


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing classical upper-main result")
    require(gate.NOTE.is_file(), "missing classical upper-main note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate", "kind drift")
    require(artifact["status"] == "rigorous_fifteen_point_source_zp_vs_classical_upper_main_and_exact_remaining_correction_split", "status drift")
    require(bool(artifact["passed"]), "classical upper-main split failed")
    require(artifact["scope"]["precisions_decimal_digits"] == list(gate.PRECISIONS), "precision roster drift")
    require(int(artifact["scope"]["classical_upper_summands"]) == gate.UPPER_END - gate.LOWER_START + 1, "summand roster drift")

    rows = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_OUTPUTS, "output roster drift")
    tolerance = arb("0.005")
    fresh_remainders: list[arb] = []
    fresh_gaps: list[arb] = []
    ctx.dps = CHECK_PRECISION
    for index, row in enumerate(rows, start=1):
        offset = Decimal(index - 8) / Decimal(100)
        require(int(row["output_index"]) == index, "output index drift")
        require(row["output_label"] == f"t{offset:+.2f}", f"output label drift at {index}")
        require(Decimal(row["target_t"]) == Decimal("1e10") + offset, f"target height drift at {index}")
        fresh_main = gate.classical_upper_main(row["target_t"], CHECK_PRECISION)
        fresh_components = component_split.exact_components(row["target_t"], CHECK_PRECISION)
        fresh_upper = arb(fresh_components["exact_complementary_upper_ball"])
        require(fresh_main.overlaps(arb(row["high_precision"]["classical_upper_main_ball"])), f"fresh upper main misses saved ball at {index}")
        require(fresh_upper.overlaps(arb(row["high_precision"]["exact_complementary_upper_ball"])), f"fresh upper target misses saved ball at {index}")
        remainder = fresh_upper - fresh_main
        source = arb(row["source_upper_zp_binary128_decimal"])
        hybrid_main = source - fresh_main
        hybrid_target = source - fresh_upper
        identity = hybrid_main - remainder - hybrid_target
        require(identity.contains(0), f"fresh upper identity excludes zero at {index}")
        require(remainder.overlaps(arb(row["exact_remaining_correction_ball"])), f"fresh remainder misses saved ball at {index}")
        require(hybrid_main.overlaps(arb(row["hybrid_minus_classical_main_ball"])), f"fresh hybrid-main gap misses saved ball at {index}")
        require(abs(hybrid_main).lower() > abs(remainder).upper(), f"fresh remainder subordination drift at {index}")
        require(abs(hybrid_main).lower() > tolerance, f"fresh tolerance rejection drift at {index}")
        fresh_remainders.append(abs(remainder))
        fresh_gaps.append(abs(hybrid_main))

    aggregate = artifact["aggregate"]
    require(int(aggregate["output_count"]) == gate.EXPECTED_OUTPUTS, "output aggregate drift")
    require(int(aggregate["precision_overlap_count"]) == gate.EXPECTED_OUTPUTS, "precision aggregate drift")
    require(int(aggregate["identity_contains_zero_count"]) == gate.EXPECTED_OUTPUTS, "identity aggregate drift")
    require(int(aggregate["remainder_subordinate_count"]) == gate.EXPECTED_OUTPUTS, "remainder aggregate drift")
    require(int(aggregate["hybrid_minus_classical_main_rejects_0p005_count"]) == gate.EXPECTED_OUTPUTS, "tolerance aggregate drift")
    require(max(fresh_remainders, key=lambda ball: ball.upper()).overlaps(arb(aggregate["maximum_exact_remaining_correction_absolute_ball"])), "maximum remainder aggregate misses fresh result")
    require(min(fresh_gaps, key=lambda ball: ball.lower()).overlaps(arb(aggregate["minimum_hybrid_minus_classical_main_absolute_ball"])), "minimum hybrid gap aggregate misses fresh result")

    decision = artifact["decision"]
    require(not bool(decision["exact_remaining_correction_explains_saved_zp_error"]), "remainder was overpromoted")
    require(bool(decision["saved_error_is_inside_hybrid_approximation_to_upper_main"]), "hybrid localization decision lost")
    require(not bool(decision["finite_split_is_height_uniform_theorem"]), "height-uniform overpromotion")
    require(not bool(decision["finite_split_has_rh_implication"]), "RH overpromotion")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "exact identity", "cannot explain", "inside the source", "no RH implication"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated classical upper-main split: "
        f"15 identities, remainder-subordinate={aggregate['remainder_subordinate_count']}, "
        f"hybrid-main-rejects-0.005={aggregate['hybrid_minus_classical_main_rejects_0p005_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
