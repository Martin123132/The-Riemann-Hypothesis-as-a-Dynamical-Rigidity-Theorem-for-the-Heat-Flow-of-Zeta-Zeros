#!/usr/bin/env python3
"""Independently validate the exact Hardy-Z Arb calibration gate."""

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

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate as gate


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_hardy(target_t: str) -> acb:
    ctx.dps = CHECK_PRECISION
    ctx.threads = 1
    t = arb(target_t)
    theta = acb(arb("0.25"), t / 2).lgamma().imag - t * arb.pi().log() / 2
    return acb(0, theta).exp() * acb(arb("0.5"), t).zeta()


def main() -> int:
    require(gate.RESULT.is_file(), "missing exact Hardy calibration result")
    require(gate.NOTE.is_file(), "missing exact Hardy calibration note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate", "kind drift")
    require(artifact["status"] == "rigorous_fifteen_point_arb_calibration_rejects_absolute_0p005_for_saved_source_run", "status drift")
    require(bool(artifact["passed"]), "calibration gate failed")
    require(artifact["scope"]["precisions_decimal_digits"] == list(gate.PRECISIONS), "precision roster drift")

    rows = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_OUTPUTS, "output roster drift")
    tolerance = arb("0.005")
    independent_absolute: list[arb] = []
    for index, row in enumerate(rows, start=1):
        offset = Decimal(index - 8) / Decimal(100)
        require(int(row["output_index"]) == index, "output index drift")
        require(row["output_label"] == f"t{offset:+.2f}", f"output label drift at {index}")
        require(Decimal(row["target_t"]) == Decimal("1e10") + offset, f"target height drift at {index}")
        low = arb(row["low_precision"]["hardy_real_ball"])
        high = arb(row["high_precision"]["hardy_real_ball"])
        require(low.overlaps(high), f"saved precision overlap drift at {index}")
        require(arb(row["low_precision"]["hardy_imag_ball"]).contains(0), f"low imaginary drift at {index}")
        require(arb(row["high_precision"]["hardy_imag_ball"]).contains(0), f"high imaginary drift at {index}")

        fresh = independent_hardy(row["target_t"])
        require(fresh.imag.contains(0), f"fresh Hardy imaginary excludes zero at {index}")
        require(fresh.real.overlaps(high), f"fresh Hardy enclosure misses saved high precision at {index}")
        source = arb(row["source_binary128_decimal"])
        signed = source - fresh.real
        absolute = abs(signed)
        internal = arb(row["current_internal_column_upper"])
        residual = absolute - internal
        require(signed.upper() < 0, f"fresh source-minus-Hardy sign drift at {index}")
        require(absolute.lower() > tolerance, f"fresh tolerance rejection drift at {index}")
        require(residual.lower() > tolerance, f"fresh internal-column guard drift at {index}")
        require(bool(row["requested_tolerance_rejected"]), f"saved decision drift at {index}")
        require(bool(row["internal_column_alone_cannot_restore_tolerance"]), f"saved internal decision drift at {index}")
        independent_absolute.append(absolute)

    aggregate = artifact["aggregate"]
    require(int(aggregate["output_count"]) == gate.EXPECTED_OUTPUTS, "output aggregate drift")
    require(int(aggregate["precision_overlap_count"]) == gate.EXPECTED_OUTPUTS, "precision aggregate drift")
    require(int(aggregate["strict_negative_source_minus_hardy_count"]) == gate.EXPECTED_OUTPUTS, "sign aggregate drift")
    require(int(aggregate["requested_tolerance_rejected_count"]) == gate.EXPECTED_OUTPUTS, "tolerance aggregate drift")
    require(int(aggregate["internal_column_insufficient_count"]) == gate.EXPECTED_OUTPUTS, "internal aggregate drift")
    require(not bool(aggregate["absolute_0p005_target_holds_for_saved_source_run"]), "saved source tolerance was overpromoted")
    require(not bool(aggregate["rh_implication"]), "finite calibration was overpromoted to RH")

    decision = artifact["decision"]
    require(not bool(decision["saved_source_outputs_meet_absolute_0p005"]), "decision tolerance drift")
    require(not bool(decision["current_internal_column_alone_can_explain_gap"]), "decision internal-column drift")
    require(bool(decision["outer_remainder_still_required"]), "outer-remainder obligation lost")
    require(not bool(decision["finite_calibration_is_evidence_for_rh"]), "RH nonpromotion guard lost")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "All 15", "farther than `0.005`", "falsifies", "does not falsify RH", "outer column"):
        require(token in note, f"note boundary token missing: {token}")

    print(
        "validated exact Hardy Arb calibration: "
        f"15/15 reject 0.005, min={aggregate['minimum_absolute_error_ball']}, "
        f"max={aggregate['maximum_absolute_error_ball']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
