#!/usr/bin/env python3
"""Certify the four shifted-flow complement rows for order eleven at lambda zero."""

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

import jensen_window_pf_compound_order10_lambda0_prefix_certificate as order10  # noqa: E402


ORDER = 11
PROBE_M = ORDER - 1
N_MIN = 0
N_MAX = 3
MAX_K = N_MAX + 2 * (ORDER - 1)
PRECISION_DPS = order10.PRECISION_DPS
COEFFICIENT_SOURCES = order10.COEFFICIENT_SOURCES
HANKEL_ROWS_SOURCE = order10.HANKEL_ROWS_SOURCE
HANKEL_SUMMARY_SOURCE = order10.HANKEL_SUMMARY_SOURCE
DEFINITION_SOURCE = order10.DEFINITION_SOURCE
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_lambda0_prefix_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order11_lambda0_prefix_certificate.md"
)


@dataclass(frozen=True)
class PrefixRow:
    n: int
    coefficient_range: list[int]
    probe_m: int
    matrix_dimension: int
    probe_sigma: int
    epsilon_11: int
    raw_H11_ball: str
    Q11_ball: str
    source_Q11_ball: str
    source_overlap: bool
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
                if not order10.decimal_equal(row.get("lam"), "0"):
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


def load_probe_rows() -> dict[int, dict]:
    selected: dict[int, dict] = {}
    with HANKEL_ROWS_SOURCE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if (
                row.get("kind") == "arb_hankel_enclosure_sign_probe"
                and order10.decimal_equal(row.get("lam"), "0")
                and int(row.get("m", -1)) == PROBE_M
                and N_MIN <= int(row.get("shift", -1)) <= N_MAX
            ):
                n = int(row["shift"])
                if n in selected:
                    raise RuntimeError(f"duplicate source probe row for n={n}")
                selected[n] = row
    missing = [n for n in range(N_MIN, N_MAX + 1) if n not in selected]
    if missing:
        raise RuntimeError(f"missing source probe rows: {missing}")
    return selected


def build_artifact() -> dict:
    flint.ctx.dps = PRECISION_DPS
    contracts = order10.validate_source_contracts()
    coefficients = load_coefficients()
    source_rows = load_probe_rows()
    epsilon = order10.signed_hankel_sign(ORDER)
    sigma = order10.probe_sign(PROBE_M)
    if epsilon != -1 or sigma != epsilon:
        raise RuntimeError("order-eleven probe orientation mismatch")

    rows = []
    for n in range(N_MIN, N_MAX + 1):
        raw = flint.arb_mat(
            [
                [coefficients[n + i + j] for j in range(ORDER)]
                for i in range(ORDER)
            ]
        ).det()
        signed = epsilon * raw
        source = source_rows[n]
        source_signed = flint.arb(source["signed_det"])
        checks = {
            "source_sigma": int(source.get("sigma", 0)) == sigma,
            "source_positive": source.get("classification") == "positive",
            "source_ok": source.get("ok") is True,
            "source_excludes_zero": (
                source.get("contains_zero") is False
                and source.get("sign_separated") is True
            ),
            "rebuild_positive": bool(signed > 0),
            "source_positive_ball": bool(source_signed > 0),
            "source_overlap": bool(signed.overlaps(source_signed)),
        }
        if not all(checks.values()):
            raise RuntimeError(f"order-eleven lambda-zero row failed at n={n}: {checks}")
        rows.append(
            PrefixRow(
                n=n,
                coefficient_range=[n, n + 2 * (ORDER - 1)],
                probe_m=PROBE_M,
                matrix_dimension=ORDER,
                probe_sigma=sigma,
                epsilon_11=epsilon,
                raw_H11_ball=order10.compact_ball(raw),
                Q11_ball=order10.compact_ball(signed),
                source_Q11_ball=source["signed_det"],
                source_overlap=True,
                classification="positive",
            )
        )

    source_paths = (
        *COEFFICIENT_SOURCES,
        HANKEL_ROWS_SOURCE,
        HANKEL_SUMMARY_SOURCE,
        DEFINITION_SOURCE,
    )
    return {
        "kind": "jensen_window_pf_compound_order11_lambda0_prefix_certificate",
        "date": "2026-07-16",
        "status": "rigorous lambda-zero order-eleven prefix certificate for n=0,1,2,3",
        "proof_boundary": (
            "This proves four signed order-eleven Hankel rows at lambda zero only. "
            "It does not prove the shifted heat ray, the analytic endpoint tail, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [order10.source_record(path) for path in source_paths],
        "normalization": {
            "coefficient": "A_k(lambda)=mu_(2k)(lambda)*k!/(2k)!",
            "definition": contracts["definition"],
            "order": ORDER,
            "probe_parameter": PROBE_M,
            "probe_matrix_dimension": PROBE_M + 1,
            "epsilon_11": epsilon,
            "probe_sigma": sigma,
            "index_identification": "probe m=10 and shift=n equals Q_(11,n)",
        },
        "finite": {
            "lambda": "0",
            "n_range": [N_MIN, N_MAX],
            "coefficient_range": [0, MAX_K],
            "precision_dps": PRECISION_DPS,
            "all_Q11_positive": True,
            "theorem": "Q_(11,n)(0)>0 for every integer 0<=n<=3",
            "rows": [asdict(row) for row in rows],
        },
        "summary": {
            "coefficient_rows": MAX_K + 1,
            "prefix_rows": N_MAX - N_MIN + 1,
            "positive_Q11_rows": N_MAX - N_MIN + 1,
            "source_overlaps": N_MAX - N_MIN + 1,
            "inconclusive_rows": 0,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_lambda0_prefix_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_lambda0_prefix_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Order-Eleven Lambda-Zero Prefix Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous certificate for the four shifted-flow complement rows. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order11_lambda0_prefix_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_lambda0_prefix_certificate.py",
        "```",
        "",
        "The source probe parameter `m=10` builds an `11x11` determinant and",
        "uses `sigma=(-1)^(10*11/2)=-1=epsilon_11`. Therefore its shift `n`",
        "is exactly the project coordinate `Q_(11,n)`.",
        "",
        "Fresh determinant rebuilds from the rigorous 520-digit coefficient",
        "balls independently overlap all four source rows and prove",
        "",
        "```text",
        "Q_(11,n)(0)>0 for every integer 0<=n<=3.",
        "```",
        "",
        "| n | rigorous Q_(11,n)(0) enclosure |",
        "|---:|---|",
    ]
    for row in artifact["finite"]["rows"]:
        lines.append(f"| {row['n']} | `{row['Q11_ball']}` |")
    lines.extend(
        [
            "",
            "This prefix is ready to compose with a future shifted heat ray on",
            "`n>=4`. It does not itself prove that ray, order twelve,",
            "PF-infinity, RH, or `Lambda<=0`.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "certified lambda-zero order-eleven prefix: "
        "24 coefficients, 4 positive Q11 rows, 4 source overlaps"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
