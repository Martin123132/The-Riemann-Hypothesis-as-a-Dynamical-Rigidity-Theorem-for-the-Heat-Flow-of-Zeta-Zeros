#!/usr/bin/env python3
"""Audit the cancellation retained by pairing endpoint remainders with W2--W4."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate as log_full
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells


SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.json"
)
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_paired_exceptional_cancellation_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_paired_exceptional_cancellation_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_paired_exceptional_cancellation_gate.py"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper_abs(value: Any) -> Fraction:
    return log_full.upper(abs(value))


def lower_abs(value: Any) -> Fraction:
    return log_full.lower(abs(value))


def main() -> int:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    require(source["status"] == "rigorous_all_374_endpoint_linearization_decomposition_validated", "paired source drift")
    rows: list[dict[str, Any]] = []
    max_ratio = (Fraction(0), 0)
    max_gap = (Fraction(0), 0)
    for row in source["rows"]:
        chain = int(row["chain"])
        upper_pair = native_q.complex_from_record(row["b_remainder"]) - native_q.complex_from_record(row["paper_w2"])
        lower_pair = (
            native_q.complex_from_record(row["c_remainder"])
            + native_q.complex_from_record(row["zero_mode"])
            - native_q.complex_from_record(row["paper_w3_or_w4"])
        )
        correction = native_q.complex_from_record(row["exact_correction"])
        gap = upper_pair + lower_pair - correction
        require(gap.contains(0), f"chain {chain} paired correction identity excludes zero")
        correction_lower = lower_abs(correction)
        require(correction_lower > 0, f"chain {chain} paired correction contains zero")
        triangle = upper_abs(upper_pair) + upper_abs(lower_pair)
        ratio = triangle / correction_lower
        gap_upper = upper_abs(gap)
        max_ratio = max(max_ratio, (ratio, chain))
        max_gap = max(max_gap, (gap_upper, chain))
        rows.append(
            {
                "chain": chain,
                "upper_endpoint_paired_magnitude_upper": decimal(upper_abs(upper_pair)),
                "lower_endpoint_paired_magnitude_upper": decimal(upper_abs(lower_pair)),
                "paired_triangle_sum_upper": decimal(triangle),
                "exact_correction_magnitude_lower": decimal(correction_lower),
                "paired_triangle_to_correction_ratio_upper": decimal(ratio),
                "paired_identity_gap_upper": decimal(gap_upper),
            }
        )
    raw_ratio = Fraction(source["aggregate"]["maximum_component_triangle_to_correction_ratio_upper"])
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_paired_exceptional_cancellation_gate",
        "status": "rigorous_all_374_paired_exceptional_cancellation_audited",
        "identity": "Q_exact-Q_paper=(R_b-W2)+(R_c+zero_mode-(W3_or_W4))",
        "rows": rows,
        "aggregate": {
            "call_count": len(rows),
            "paired_identity_count": sum(float(row["paired_identity_gap_upper"]) < 1e-10 for row in rows),
            "maximum_paired_identity_gap_upper": decimal(max_gap[0]),
            "maximum_paired_identity_gap_witness": max_gap[1],
            "maximum_raw_component_ratio_upper": decimal(raw_ratio),
            "maximum_paired_triangle_to_correction_ratio_upper": decimal(max_ratio[0]),
            "maximum_paired_triangle_to_correction_ratio_witness": max_ratio[1],
            "raw_to_paired_worst_ratio_reduction_lower": decimal(raw_ratio / max_ratio[0]),
        },
        "dependencies": {
            "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Rigorous finite-roster cancellation audit only. The ratios diagnose absolute-value loss and are not "
            "height-uniform bounds, fitted proof constants, recurrence budgets, or proofs of Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    aggregate = artifact["aggregate"]
    NOTE.write_text(
        f"""# Hardy block-20 paired exceptional-cancellation gate

Date: 2026-08-07

Status: rigorous all-374 finite-roster cancellation audit; not a height-uniform proof and not a proof of RH

The exact correction is grouped according to the paper approximation it
replaces:

```text
Q_exact-Q_paper=(R_b-W2)+(R_c+zero_mode-(W3 or W4)).
```

Both grouped endpoint families sum to the exact correction on all 374 calls.
The maximum identity gap is
`{aggregate['maximum_paired_identity_gap_upper']}`.  Pairing reduces the worst
component-triangle/correction ratio from
`{aggregate['maximum_raw_component_ratio_upper']}` to
`{aggregate['maximum_paired_triangle_to_correction_ratio_upper']}`, a factor
of at least `{aggregate['raw_to_paired_worst_ratio_reduction_lower']}`.  The
remaining worst witness is chain
`{aggregate['maximum_paired_triangle_to_correction_ratio_witness']}`.

Thus W2--W4 must be bounded jointly with the exact endpoint remainder they
replace.  Even after that pairing, separate absolute values for the upper and
lower endpoint families can lose a factor above 210, so a competitive uniform
theorem must preserve their common Euler--Maclaurin structure or cancellation.
These finite ratios are diagnostics, not fitted constants, recurrence budgets,
or proofs of `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
""",
        encoding="utf-8",
    )
    print("built paired exceptional-cancellation gate: 374/374 identities, worst paired ratio below 211")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
