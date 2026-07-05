#!/usr/bin/env python3
"""Finite reshaped-Hankel sign-consistency point audit.

This is a theorem-search diagnostic for the Grussler-Damm style
sign-consistency reduction.  It tests all k-column minors of the k-row
reshaped Hankel matrix built from cached A_k(lambda) values.

The input coefficients are decimal cache values rationalized exactly.  This
is not an interval certificate for the transcendental coefficients.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
import json
import math
from math import gcd
from pathlib import Path
from itertools import combinations


DEFAULT_LAMBDAS = ("0", "1e-6", "1e-4", "1e-2", "1e-1")
REPO_ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class AuditRow:
    lam: str
    order_k: int
    n_cols: int
    max_coefficient_index: int
    expected_sign: int
    tests: int
    positive: int
    negative: int
    zero: int
    first_bad_columns: tuple[int, ...] | None
    minimum_signed_value: str | None
    coefficient_common_denominator_digits: int


def dec(value: object) -> Decimal:
    return Decimal(str(value))


def short_fraction(value: Fraction, digits: int = 18) -> str:
    as_decimal = Decimal(value.numerator) / Decimal(value.denominator)
    return f"{as_decimal:.{digits}E}"


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def det_bareiss_int(matrix: list[list[int]]) -> int:
    n = len(matrix)
    work = [list(row) for row in matrix]
    sign = 1
    previous_pivot = 1
    for col in range(n - 1):
        pivot = None
        for row in range(col, n):
            if work[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return 0
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            sign *= -1
        pivot_value = work[col][col]
        for row in range(col + 1, n):
            row_factor = work[row][col]
            for j in range(col + 1, n):
                work[row][j] = (
                    work[row][j] * pivot_value - row_factor * work[col][j]
                ) // previous_pivot
            work[row][col] = 0
        previous_pivot = pivot_value
    return sign * work[n - 1][n - 1]


def integer_coefficients(values: list[str]) -> tuple[list[int], int]:
    fractions = [Fraction(value) for value in values]
    common_denominator = 1
    for value in fractions:
        common_denominator = lcm(common_denominator, value.denominator)
    return (
        [
            value.numerator * (common_denominator // value.denominator)
            for value in fractions
        ],
        common_denominator,
    )


def signed_hankel_sign(order_k: int) -> int:
    return -1 if (order_k * (order_k - 1) // 2) % 2 else 1


def load_coeff_rows(path: Path) -> dict[Decimal, dict]:
    rows: dict[Decimal, dict] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") != "coefficients":
                continue
            rows[dec(row["lam"])] = row
    return rows


def audit_row(row: dict, order_k: int, n_cols: int) -> AuditRow:
    max_index = n_cols + order_k - 2
    if max_index > int(row["max_k"]):
        raise ValueError(
            f"need coefficients through {max_index}, but row has max_k={row['max_k']}"
        )

    coefficients, common_denominator = integer_coefficients(row["A"])
    expected = signed_hankel_sign(order_k)
    positive = negative = zero = 0
    first_bad_columns: tuple[int, ...] | None = None
    minimum_signed_int: int | None = None

    for columns in combinations(range(n_cols), order_k):
        determinant = det_bareiss_int(
            [[coefficients[i + column] for column in columns] for i in range(order_k)]
        )
        signed_value = expected * determinant
        if signed_value > 0:
            positive += 1
            if minimum_signed_int is None or signed_value < minimum_signed_int:
                minimum_signed_int = signed_value
        elif signed_value < 0:
            negative += 1
            if first_bad_columns is None:
                first_bad_columns = tuple(columns)
        else:
            zero += 1
            if first_bad_columns is None:
                first_bad_columns = tuple(columns)

    return AuditRow(
        lam=str(row["lam"]),
        order_k=order_k,
        n_cols=n_cols,
        max_coefficient_index=max_index,
        expected_sign=expected,
        tests=math.comb(n_cols, order_k),
        positive=positive,
        negative=negative,
        zero=zero,
        first_bad_columns=first_bad_columns,
        minimum_signed_value=(
            None
            if minimum_signed_int is None
            else short_fraction(Fraction(minimum_signed_int, common_denominator ** order_k))
        ),
        coefficient_common_denominator_digits=len(str(common_denominator)),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coeff-cache",
        type=Path,
        default=REPO_ROOT / "work/rh_compute/results/repro_hankel_15c_coefficients.jsonl",
    )
    parser.add_argument(
        "--lambdas",
        default=",".join(DEFAULT_LAMBDAS),
        help="Comma-separated lambda list to load from the coefficient cache.",
    )
    parser.add_argument("--orders", default="2,3,4,5", help="Comma-separated k orders.")
    parser.add_argument("--n-cols", type=int, default=18, help="Number of reshaped columns.")
    parser.add_argument("--output", type=Path, default=None, help="Optional JSON output path.")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    getcontext().prec = 120
    args = build_parser().parse_args()
    rows = load_coeff_rows(args.coeff_cache)
    targets = [dec(part.strip()) for part in args.lambdas.split(",") if part.strip()]
    orders = [int(part.strip()) for part in args.orders.split(",") if part.strip()]

    audit_rows: list[AuditRow] = []
    for lam in targets:
        if lam not in rows:
            raise KeyError(f"lambda {lam} not found in {args.coeff_cache}")
        for order_k in orders:
            audit_rows.append(audit_row(rows[lam], order_k, args.n_cols))

    ok = all(row.negative == 0 and row.zero == 0 and row.positive == row.tests for row in audit_rows)
    payload = {
        "kind": "hankel_sign_consistency_reduction_point_audit",
        "date": "2026-07-05",
        "proof_boundary": (
            "Exact rational point audit of cached decimal coefficients only; "
            "not an interval certificate, not all-order sign consistency, "
            "and not a proof of Jensen hyperbolicity, RH, or Lambda <= 0."
        ),
        "coeff_cache": str(args.coeff_cache.relative_to(REPO_ROOT) if args.coeff_cache.is_absolute() else args.coeff_cache),
        "orders": orders,
        "n_cols": args.n_cols,
        "rows": [asdict(row) for row in audit_rows],
        "ok": ok,
    }

    if args.output is not None:
        output = args.output
        if not output.is_absolute():
            output = REPO_ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        for row in audit_rows:
            print(
                "OK" if row.negative == 0 and row.zero == 0 else "FAIL",
                "reshaped Hankel sign-consistency point audit:",
                f"lambda={row.lam}",
                f"k={row.order_k}",
                f"N={row.n_cols}",
                f"{row.positive}/{row.tests} positive",
            )
        print(f"validated {sum(1 for row in audit_rows if row.negative == 0 and row.zero == 0)}/{len(audit_rows)} reshaped Hankel sign-consistency point audits")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
