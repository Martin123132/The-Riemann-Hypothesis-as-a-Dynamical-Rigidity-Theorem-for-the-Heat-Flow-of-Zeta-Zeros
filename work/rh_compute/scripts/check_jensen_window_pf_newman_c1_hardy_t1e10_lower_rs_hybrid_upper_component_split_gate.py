#!/usr/bin/env python3
"""Independently check the lower-RS versus hybrid-upper residual split."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate as gate


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing component-split result")
    require(gate.NOTE.is_file(), "missing component-split note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate", "kind drift")
    require(artifact["status"] == "rigorous_fifteen_point_lower_rs_hybrid_upper_and_final_addition_residual_split", "status drift")
    require(bool(artifact["passed"]), "component split failed")
    require(artifact["scope"]["precisions_decimal_digits"] == list(gate.PRECISIONS), "precision roster drift")
    require(int(artifact["scope"]["lower_sum_cutoff"]) == gate.EXPECTED_NC, "NC drift")
    require(int(artifact["scope"]["riemann_siegel_index"]) == gate.EXPECTED_RSN, "RSN drift")

    rows = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_OUTPUTS, "output roster drift")
    tolerance = arb("0.005")
    fresh_dominance: list[str] = []
    fresh_upper_rejections = 0
    ctx.dps = CHECK_PRECISION
    for index, row in enumerate(rows, start=1):
        offset = Decimal(index - 8) / Decimal(100)
        require(int(row["output_index"]) == index, "output index drift")
        require(row["output_label"] == f"t{offset:+.2f}", f"output label drift at {index}")
        require(Decimal(row["target_t"]) == Decimal("1e10") + offset, f"target height drift at {index}")
        fresh = gate.exact_components(row["target_t"], CHECK_PRECISION)
        for field in ("hardy_real_ball", "exact_lower_component_ball", "exact_complementary_upper_ball"):
            require(arb(fresh[field]).overlaps(arb(row["high_precision"][field])), f"fresh component misses saved enclosure at output {index}: {field}")

        source_upper = arb(row["source_upper_zp_binary128_decimal"])
        source_lower = arb(row["source_lower_rszsum_binary128_decimal"])
        source_final = arb(row["source_final_binary128_decimal"])
        exact_hardy = arb(fresh["hardy_real_ball"])
        exact_lower = arb(fresh["exact_lower_component_ball"])
        exact_upper = arb(fresh["exact_complementary_upper_ball"])
        lower_gap = source_lower - exact_lower
        upper_gap = source_upper - exact_upper
        addition_gap = source_final - source_upper - source_lower
        direct_total = source_final - exact_hardy
        identity = lower_gap + upper_gap + addition_gap - direct_total
        require(identity.contains(0), f"fresh identity excludes zero at output {index}")
        for field, fresh_ball in (
            ("lower_source_minus_exact_ball", lower_gap),
            ("hybrid_upper_source_minus_exact_ball", upper_gap),
            ("final_addition_rounding_ball", addition_gap),
            ("direct_source_minus_hardy_ball", direct_total),
        ):
            require(fresh_ball.overlaps(arb(row[field])), f"fresh residual misses saved ball at output {index}: {field}")

        magnitudes = {"lower": abs(lower_gap), "upper": abs(upper_gap), "addition": abs(addition_gap)}
        dominant = max(magnitudes, key=lambda key: magnitudes[key].upper())
        require(magnitudes[dominant].lower() > max(magnitudes[key].upper() for key in magnitudes if key != dominant), f"fresh dominance is not rigorous at output {index}")
        require(row["dominant_component"] == dominant, f"saved dominance drift at output {index}")
        fresh_dominance.append(dominant)
        upper_rejects = magnitudes["upper"].lower() > tolerance
        require(bool(row["hybrid_upper_alone_rejects_0p005"]) == upper_rejects, f"upper tolerance decision drift at output {index}")
        fresh_upper_rejections += int(upper_rejects)

    aggregate = artifact["aggregate"]
    require(int(aggregate["output_count"]) == gate.EXPECTED_OUTPUTS, "output aggregate drift")
    require(int(aggregate["precision_overlap_count"]) == gate.EXPECTED_OUTPUTS, "precision aggregate drift")
    require(int(aggregate["identity_contains_zero_count"]) == gate.EXPECTED_OUTPUTS, "identity aggregate drift")
    require(int(aggregate["lower_component_dominant_count"]) == fresh_dominance.count("lower"), "lower dominance aggregate drift")
    require(int(aggregate["hybrid_upper_component_dominant_count"]) == fresh_dominance.count("upper"), "upper dominance aggregate drift")
    require(int(aggregate["final_addition_component_dominant_count"]) == fresh_dominance.count("addition"), "addition dominance aggregate drift")
    require(int(aggregate["hybrid_upper_alone_rejects_0p005_count"]) == fresh_upper_rejections, "upper tolerance aggregate drift")

    decision = artifact["decision"]
    require(bool(decision["component_split_is_exact_identity"]), "identity decision drift")
    require(not bool(decision["finite_split_is_height_uniform_theorem"]), "height-uniform overpromotion")
    require(not bool(decision["finite_split_has_rh_implication"]), "RH overpromotion")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "splits identically", "Lower-component dominance", "Hybrid-upper dominance", "does not", "exact identity"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated lower-RS/hybrid-upper split: "
        f"15 identities, lower-dominant={aggregate['lower_component_dominant_count']}, "
        f"upper-dominant={aggregate['hybrid_upper_component_dominant_count']}, "
        f"upper-rejects-0.005={aggregate['hybrid_upper_alone_rejects_0p005_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
