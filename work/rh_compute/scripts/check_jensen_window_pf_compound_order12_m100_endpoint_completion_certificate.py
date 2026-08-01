#!/usr/bin/env python3
"""Validate the finite lambda=-100 order-twelve endpoint completion."""

from __future__ import annotations

from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order12_m100_endpoint_completion_certificate as generator  # noqa: E402


def main() -> int:
    artifact = json.loads(generator.DEFAULT_OUT.read_text(encoding="utf-8"))
    issues = []
    if artifact.get("kind") != "jensen_window_pf_compound_order12_m100_endpoint_completion_certificate":
        issues.append("bad artifact kind")
    if artifact.get("status") != "rigorous lambda=-100 order-twelve finite endpoint sign chart through n=1492":
        issues.append("bad artifact status")

    rebuilt = generator.build_artifact()
    if artifact != rebuilt:
        issues.append("rebuilt artifact differs")

    expected_exact = {
        "condensation": "Q_(12,n)*Q_(10,n+2)=Q_(11,n+1)^2-Q_(11,n)*Q_(11,n+2)",
        "relative_coordinate": "R_n=Q_(11,n+1)^2/(Q_(11,n)*Q_(11,n+2))-1",
        "inherited_signs": "Q_(12,n)(-100)<0 for n=0,1,2,3 and Q_(12,n)(-100)>0 for every 4<=n<=1240",
        "collar_theorem": "Q_(12,n)(-100)>0 for every 1241<=n<=1492",
        "finite_sign_chart": "Q_(12,n)(-100)<0 for n=0,1,2,3 and Q_(12,n)(-100)>0 for every 4<=n<=1492",
        "remaining_tail_target": "Q_(12,n)(-100)>0 for every n>=1493",
    }
    if artifact.get("exact") != expected_exact:
        issues.append("exact theorem block changed")

    collar = artifact.get("order11_collar", {})
    collar_rows = collar.get("rows", [])
    if [row.get("n") for row in collar_rows] != list(range(1241, 1495)):
        issues.append("Q11 collar rows are not contiguous n=1241..1494")
    if collar.get("inherited_overlap_rows") != 1243 or collar.get("inherited_overlap_failures") != []:
        issues.append("inherited Q11 overlap gate changed")
    if any(row.get("Q11_sign") != "positive" for row in collar_rows):
        issues.append("nonpositive Q11 collar row")

    finite = artifact.get("finite", {})
    finite_rows = finite.get("rows", [])
    if [row.get("n") for row in finite_rows] != list(range(1241, 1493)):
        issues.append("Q12 endpoint rows are not contiguous n=1241..1492")
    for row in finite_rows:
        n = row.get("n")
        if row.get("Q12_sign") != "positive" or not all(row.get("checks", {}).values()):
            issues.append(f"failed Q12 row at n={n}")
        for key in ("relative_Q11_margin_lower", "Q12_lower"):
            try:
                lower = Decimal(row[key])
            except Exception as exc:
                issues.append(f"unparseable {key} at n={n}: {exc}")
            else:
                if lower <= 0:
                    issues.append(f"nonpositive {key} at n={n}")

    expected_summary = {
        "rows": 7,
        "ready_rows": 6,
        "open_rows": 1,
        "coefficient_rows": 1515,
        "new_coefficient_rows": 252,
        "coefficient_chunks": 4,
        "rebuilt_Q11_rows": 1495,
        "inherited_Q11_overlap_rows": 1243,
        "Q11_collar_rows": 254,
        "endpoint_collar_rows": 252,
        "positive_endpoint_collar_rows": 252,
        "negative_endpoint_collar_rows": 0,
        "inconclusive_endpoint_collar_rows": 0,
        "finite_sign_chart_theorems": 1,
        "open_analytic_tail_targets": 1,
        "all_shift_endpoint_positivity_theorems": 0,
        "rh_claims": 0,
    }
    if artifact.get("summary") != expected_summary:
        issues.append("summary changed")

    note = generator.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "exactly 252 coefficients",
        expected_exact["collar_theorem"],
        expected_exact["finite_sign_chart"],
        "The first four endpoint rows remain rigorously negative",
        "remaining endpoint obligation",
    ):
        if marker not in note:
            issues.append(f"note misses marker: {marker}")

    if issues:
        for issue in issues:
            print(f"issue: {issue}")
        return 1
    print(
        "validated order-twelve endpoint completion: "
        "1515 coefficients, 1243 inherited overlaps, "
        "252 positive collar rows, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
