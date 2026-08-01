#!/usr/bin/env python3
"""Certify the four lambda-zero complement rows for order twelve."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order10_lambda0_prefix_certificate as source  # noqa: E402


ORDER = 12
N_MIN = 0
N_MAX = 3
MAX_K = N_MAX + 2 * (ORDER - 1)
PRECISION_DPS = source.PRECISION_DPS
COEFFICIENT_SOURCES = source.COEFFICIENT_SOURCES
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order12_lambda0_prefix_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order12_lambda0_prefix_certificate.md"


@dataclass(frozen=True)
class PrefixRow:
    n: int
    coefficient_range: list[int]
    matrix_dimension: int
    epsilon_12: int
    raw_H12_ball: str
    Q12_ball: str
    classification: str


def load_coefficients() -> dict[int, flint.arb]:
    coefficients: dict[int, flint.arb] = {}
    for path in COEFFICIENT_SOURCES:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                if row.get("kind") != "acb_coefficient_enclosure":
                    continue
                if not source.decimal_equal(row.get("lam"), "0"):
                    continue
                k = int(row["k"])
                if k > MAX_K:
                    continue
                if k in coefficients:
                    raise RuntimeError(f"duplicate lambda-zero coefficient A_{k}")
                full_mu = flint.arb(row["full_mu_ball"])
                coefficient = flint.arb(row["A_ball"])
                normalized = full_mu * math.factorial(k) / math.factorial(2 * k)
                if not coefficient.overlaps(normalized):
                    raise RuntimeError(f"A_{k} misses mu_{{2k}} k!/(2k)!")
                if not coefficient > 0:
                    raise RuntimeError(f"A_{k}(0) is not strictly positive")
                coefficients[k] = coefficient
    missing = [k for k in range(MAX_K + 1) if k not in coefficients]
    if missing:
        raise RuntimeError(f"missing lambda-zero coefficient balls: {missing}")
    return coefficients


def build_artifact() -> dict:
    flint.ctx.dps = PRECISION_DPS
    contracts = source.validate_source_contracts()
    coefficients = load_coefficients()
    epsilon = source.signed_hankel_sign(ORDER)
    if epsilon != 1:
        raise RuntimeError("order-twelve orientation mismatch")
    rows = []
    for n in range(N_MIN, N_MAX + 1):
        raw = flint.arb_mat(
            [[coefficients[n + i + j] for j in range(ORDER)] for i in range(ORDER)]
        ).det()
        signed = epsilon * raw
        if not signed > 0:
            raise RuntimeError(f"Q12 lambda-zero row failed at n={n}: {signed}")
        rows.append(
            PrefixRow(
                n=n,
                coefficient_range=[n, n + 2 * (ORDER - 1)],
                matrix_dimension=ORDER,
                epsilon_12=epsilon,
                raw_H12_ball=source.compact_ball(raw),
                Q12_ball=source.compact_ball(signed),
                classification="positive",
            )
        )
    source_paths = (
        *COEFFICIENT_SOURCES,
        source.HANKEL_SUMMARY_SOURCE,
        source.DEFINITION_SOURCE,
    )
    return {
        "kind": "jensen_window_pf_compound_order12_lambda0_prefix_certificate",
        "date": "2026-07-22",
        "status": "rigorous lambda-zero order-twelve prefix certificate for n=0,1,2,3",
        "proof_boundary": "This proves four signed order-twelve Hankel rows at lambda zero only. It does not prove the delayed heat ray, remaining endpoint block, continuum theorem, PF-infinity, RH, or Lambda<=0.",
        "sources": [source.source_record(path) for path in source_paths],
        "normalization": {
            "coefficient": "A_k(lambda)=mu_(2k)(lambda)*k!/(2k)!",
            "definition": contracts["definition"],
            "order": ORDER,
            "epsilon_12": epsilon,
            "orientation": "Q_(12,n)=H_(12,n)",
        },
        "finite": {
            "lambda": "0",
            "n_range": [N_MIN, N_MAX],
            "coefficient_range": [0, MAX_K],
            "precision_dps": PRECISION_DPS,
            "all_Q12_positive": True,
            "theorem": "Q_(12,n)(0)>0 for every integer 0<=n<=3",
            "rows": [asdict(row) for row in rows],
        },
        "summary": {
            "coefficient_rows": MAX_K + 1,
            "prefix_rows": 4,
            "positive_Q12_rows": 4,
            "inconclusive_rows": 0,
            "orders_above_12": 0,
            "rh_claims": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order12_lambda0_prefix_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order12_lambda0_prefix_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Order-Twelve Lambda-Zero Prefix Certificate",
        "",
        "Date: 2026-07-22",
        "",
        "Status: rigorous certificate for four lambda-zero complement rows.",
        "This is not all-shift order twelve or PF-infinity, and it is not a proof of RH or `Lambda<=0`.",
        "",
        "Fresh 520-digit outward-rounded determinant rebuilds prove",
        "",
        "```text",
        "Q_(12,n)(0)>0 for every integer 0<=n<=3.",
        "```",
        "",
        "| n | rigorous Q_(12,n)(0) enclosure |",
        "|---:|---|",
    ]
    for row in artifact["finite"]["rows"]:
        lines.append(f"| {row['n']} | `{row['Q12_ball']}` |")
    lines.extend(
        [
            "",
            "These rows are ready to complement a future delayed heat ray on",
            "`n>=4`. They do not prove that ray or any order above twelve.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print("certified lambda-zero order-twelve prefix: 26 coefficients, 4 positive Q12 rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
