#!/usr/bin/env python3
"""Certify the four missing order-ten Hankel rows at lambda zero."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


ORDER = 10
PROBE_M = ORDER - 1
N_MIN = 0
N_MAX = 3
MAX_K = N_MAX + 2 * (ORDER - 1)
PRECISION_DPS = 520

RESULTS = REPO_ROOT / "work/rh_compute/results"
COEFFICIENT_SOURCES = (
    RESULTS / "acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl",
    RESULTS / "acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl",
    RESULTS / "acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl",
)
HANKEL_ROWS_SOURCE = (
    RESULTS / "arb_hankel_enclosure_lam0_m0_m19_s0_s24_dps520.jsonl"
)
HANKEL_SUMMARY_SOURCE = (
    RESULTS / "arb_hankel_enclosure_lam0_m0_m19_s0_s24_dps520_summary.json"
)
DEFINITION_SOURCE = (
    RESULTS / "jensen_window_pf_all_order_endpoint_heat_reduction.json"
)
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order10_lambda0_prefix_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order10_lambda0_prefix_certificate.md"


@dataclass(frozen=True)
class PrefixRow:
    n: int
    coefficient_range: list[int]
    probe_m: int
    matrix_dimension: int
    probe_sigma: int
    epsilon_10: int
    raw_H10_ball: str
    Q10_ball: str
    source_Q10_ball: str
    source_overlap: bool
    classification: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_path(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT)).replace("\\", "/")


def source_record(path: Path) -> dict:
    return {"path": relative_path(path), "sha256": sha256(path)}


def decimal_equal(left: object, right: object) -> bool:
    return Decimal(str(left)) == Decimal(str(right))


def signed_hankel_sign(order: int) -> int:
    return -1 if ((order * (order - 1) // 2) % 2) else 1


def probe_sign(probe_m: int) -> int:
    return -1 if ((probe_m * (probe_m + 1) // 2) % 2) else 1


def compact_ball(value: flint.arb, digits: int = 90) -> str:
    text = value.str(digits)
    if not flint.arb(text).contains(value):
        raise RuntimeError("serialized Arb ball does not contain the computed ball")
    return text


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
                if not decimal_equal(row.get("lam"), "0"):
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
                    raise RuntimeError(f"A_{k} does not overlap mu_{{2k}} k!/(2k)!")
                if not coefficient > 0:
                    raise RuntimeError(f"A_{k}(0) is not sign-separated positive")
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
                and decimal_equal(row.get("lam"), "0")
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


def validate_source_contracts() -> dict:
    summary = json.loads(HANKEL_SUMMARY_SOURCE.read_text(encoding="utf-8"))
    expected_summary = {
        "kind": "arb_hankel_enclosure_sign_probe_summary",
        "lam": "0",
        "m_min": 0,
        "m_max": 19,
        "shift_min": 0,
        "shift_max": 24,
        "dps": PRECISION_DPS,
        "needed_max_k": 62,
        "rows": 500,
        "ok": 500,
        "failed_or_inconclusive": 0,
        "all_ok": True,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            raise RuntimeError(f"Hankel summary contract changed at {key}")
    if summary.get("counts") != {"positive": 500}:
        raise RuntimeError("Hankel summary sign counts changed")

    definitions = json.loads(DEFINITION_SOURCE.read_text(encoding="utf-8"))
    expected_definition = (
        "H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m), H_(0,n)=1, "
        "epsilon_m=(-1)^binom(m,2), Q_(m,n)=epsilon_m*H_(m,n)"
    )
    if definitions.get("exact", {}).get("definitions") != expected_definition:
        raise RuntimeError("signed Hankel definition contract changed")
    return {"summary": summary, "definition": expected_definition}


def build_artifact() -> dict:
    flint.ctx.dps = PRECISION_DPS
    contracts = validate_source_contracts()
    coefficients = load_coefficients()
    source_rows = load_probe_rows()

    epsilon = signed_hankel_sign(ORDER)
    sigma = probe_sign(PROBE_M)
    if epsilon != -1 or sigma != epsilon:
        raise RuntimeError("order/probe sign identification failed")

    rows: list[PrefixRow] = []
    for n in range(N_MIN, N_MAX + 1):
        matrix = flint.arb_mat(
            [[coefficients[n + i + j] for j in range(ORDER)] for i in range(ORDER)]
        )
        raw = matrix.det()
        signed = raw if epsilon > 0 else -raw
        source = source_rows[n]
        source_signed = flint.arb(source["signed_det"])
        if int(source.get("sigma", 0)) != sigma:
            raise RuntimeError(f"source probe sign changed at n={n}")
        if source.get("classification") != "positive" or not source.get("ok"):
            raise RuntimeError(f"source probe is not certified positive at n={n}")
        if source.get("contains_zero") or not source.get("sign_separated"):
            raise RuntimeError(f"source probe is not sign-separated at n={n}")
        if not signed > 0 or not source_signed > 0:
            raise RuntimeError(f"Q_(10,{n})(0) is not sign-separated positive")
        if not signed.overlaps(source_signed):
            raise RuntimeError(f"independent determinant misses source row at n={n}")
        rows.append(
            PrefixRow(
                n=n,
                coefficient_range=[n, n + 2 * (ORDER - 1)],
                probe_m=PROBE_M,
                matrix_dimension=ORDER,
                probe_sigma=sigma,
                epsilon_10=epsilon,
                raw_H10_ball=compact_ball(raw),
                Q10_ball=compact_ball(signed),
                source_Q10_ball=source["signed_det"],
                source_overlap=True,
                classification="positive",
            )
        )

    sources = [source_record(path) for path in COEFFICIENT_SOURCES]
    sources.extend(
        source_record(path)
        for path in (HANKEL_ROWS_SOURCE, HANKEL_SUMMARY_SOURCE, DEFINITION_SOURCE)
    )
    return {
        "kind": "jensen_window_pf_compound_order10_lambda0_prefix_certificate",
        "date": "2026-07-16",
        "status": "rigorous lambda-zero order-ten prefix certificate for n=0,1,2,3",
        "proof_boundary": (
            "This proves four contiguous signed Hankel rows at lambda=0 only. "
            "It does not prove positivity at negative lambda, the order-ten all-shift "
            "theorem until the delayed heat tail is composed, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "normalization": {
            "coefficient": "A_k(lambda)=mu_(2k)(lambda)*k!/(2k)!",
            "definition": contracts["definition"],
            "order": ORDER,
            "probe_parameter": PROBE_M,
            "probe_matrix_dimension": PROBE_M + 1,
            "epsilon_10": epsilon,
            "probe_sigma": sigma,
            "index_identification": "probe m=9 and shift=n equals Q_(10,n)",
        },
        "finite": {
            "lambda": "0",
            "n_range": [N_MIN, N_MAX],
            "coefficient_range": [0, MAX_K],
            "precision_dps": PRECISION_DPS,
            "all_Q10_positive": True,
            "theorem": "Q_(10,n)(0)>0 for every integer 0<=n<=3",
            "rows": [asdict(row) for row in rows],
        },
        "summary": {
            "coefficient_rows": MAX_K + 1,
            "prefix_rows": N_MAX - N_MIN + 1,
            "positive_Q10_rows": N_MAX - N_MIN + 1,
            "source_overlaps": N_MAX - N_MIN + 1,
            "inconclusive_rows": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order10_lambda0_prefix_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order10_lambda0_prefix_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Order-Ten Lambda-Zero Prefix Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous finite certificate for the four lambda-zero prefix rows. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order10_lambda0_prefix_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_lambda0_prefix_certificate.py",
        "```",
        "",
        "## Index And Normalization",
        "",
        "The coefficient enclosure source uses",
        "`A_k(lambda)=mu_(2k)(lambda) k!/(2k)!`. The project coordinate is",
        "",
        "```text",
        "H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m)",
        "epsilon_m=(-1)^binom(m,2)",
        "Q_(m,n)=epsilon_m H_(m,n).",
        "```",
        "",
        "The source probe parameter `m=9` builds a `10x10` matrix and uses",
        "`sigma=(-1)^(9*10/2)=-1=epsilon_10`. Thus probe shift `n` is exactly",
        "`Q_(10,n)`, with no rescaling or sign conversion left over.",
        "",
        "## Certified Rows",
        "",
        "Rebuilding each determinant from the rigorous lambda-zero `A_ball` inputs",
        "at 520 decimal digits, and independently overlapping the pre-existing",
        "signed-Hankel row, proves",
        "",
        "```text",
        "Q_(10,n)(0)>0 for every integer 0<=n<=3.",
        "```",
        "",
        "| n | rigorous Q_(10,n)(0) enclosure |",
        "|---:|---|",
    ]
    for row in artifact["finite"]["rows"]:
        lines.append(f"| {row['n']} | `{row['Q10_ball']}` |")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            artifact["proof_boundary"],
            "In particular, this finite certificate alone is not an all-shift",
            "order-ten theorem and is not a proof of PF-infinity, RH, or `Lambda<=0`.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "certified lambda-zero order-ten prefix: "
        "4 positive Q10 rows, 4 source overlaps, 0 inconclusive"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
