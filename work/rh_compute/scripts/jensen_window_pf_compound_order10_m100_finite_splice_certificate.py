#!/usr/bin/env python3
"""Certify the two lambda=-100 order-ten rows after the finite scan."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order9_m100_prefix_certificate as prefix  # noqa: E402
import jensen_window_pf_endpoint_order10_counterexample as endpoint  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_m100_finite_splice_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order10_m100_finite_splice_certificate.md"
)
ENDPOINT_SOURCE = endpoint.DEFAULT_OUT
NEW_COEFFICIENT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order10_k1259_k1260_dps220.jsonl"
)
NEW_COEFFICIENT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order10_k1259_k1260_dps220_summary.json"
)
SOURCE_PATHS = (*endpoint.SOURCE_PATHS, NEW_COEFFICIENT_SOURCE)
SPLICE_INDICES = (1241, 1242)
MAX_COEFFICIENT_INDEX = 1260
PRECISION_BITS = 4096


@dataclass(frozen=True)
class SpliceRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_sources() -> dict:
    old_endpoint = load_json(ENDPOINT_SOURCE)
    finite = old_endpoint.get("finite_scan", {})
    coefficient_summary = load_json(NEW_COEFFICIENT_SUMMARY)
    if finite.get("negative_indices") != [0, 1, 2, 3]:
        raise RuntimeError("order-ten negative endpoint prefix changed")
    if finite.get("positive_range") != [4, 1240]:
        raise RuntimeError("order-ten finite positive range changed")
    if finite.get("inconclusive_rows") != 0:
        raise RuntimeError("order-ten finite scan has inconclusive rows")
    expected_summary = {
        "kind": "acb_coefficient_enclosure_summary",
        "rows": 2,
        "k_min": 1259,
        "k_max": 1260,
        "lambdas": ["-100.0"],
        "n_sum": 70,
        "cutoff": "7",
        "dps": 220,
        "abs_tol": "1e-150",
    }
    for key, expected in expected_summary.items():
        if coefficient_summary.get(key) != expected:
            raise RuntimeError(
                f"two-row coefficient source changed at {key}: "
                f"{coefficient_summary.get(key)!r}"
            )
    return {
        "old_endpoint": (
            "Q_(10,n)(-100)<0 for 0<=n<=3 and "
            "Q_(10,n)(-100)>0 for 4<=n<=1240"
        ),
        "old_endpoint_sha256": sha256(ENDPOINT_SOURCE),
        "new_coefficients": "A_1259(-100)>0 and A_1260(-100)>0",
        "new_coefficient_sha256": sha256(NEW_COEFFICIENT_SOURCE),
        "new_coefficient_summary_sha256": sha256(NEW_COEFFICIENT_SUMMARY),
        "new_coefficient_summary": coefficient_summary,
    }


def direct_splice_rows(values: dict[int, flint.arb]) -> list[dict]:
    rows = []
    for n in SPLICE_INDICES:
        q8 = endpoint.direct_signed_hankel(values, 8, n + 2)
        q9_left = endpoint.direct_signed_hankel(values, 9, n)
        q9_center = endpoint.direct_signed_hankel(values, 9, n + 1)
        q9_right = endpoint.direct_signed_hankel(values, 9, n + 2)
        q10 = endpoint.direct_signed_hankel(values, 10, n)
        schur = endpoint.direct_deep_schur(values, 10, n + 9)
        relative = q9_center**2 / (q9_left * q9_right) - 1
        toda_rhs = q9_center**2 - q9_left * q9_right
        toda_lhs = q10 * q8
        toda_residual = toda_lhs - toda_rhs
        schur_residual = schur - q10 / values[0] ** 10
        checks = {
            "q8_positive": endpoint.sign_class(q8) == "positive",
            "q9_left_positive": endpoint.sign_class(q9_left) == "positive",
            "q9_center_positive": endpoint.sign_class(q9_center) == "positive",
            "q9_right_positive": endpoint.sign_class(q9_right) == "positive",
            "q10_positive": endpoint.sign_class(q10) == "positive",
            "deep_schur_positive": endpoint.sign_class(schur) == "positive",
            "relative_margin_positive": endpoint.sign_class(relative) == "positive",
            "toda_lhs_positive": endpoint.sign_class(toda_lhs) == "positive",
            "toda_rhs_positive": endpoint.sign_class(toda_rhs) == "positive",
            "toda_residual_contains_zero": bool(toda_residual.contains(0)),
            "schur_residual_contains_zero": bool(schur_residual.contains(0)),
        }
        if not all(checks.values()):
            raise RuntimeError(f"direct order-ten splice failed at n={n}: {checks}")
        rows.append(
            {
                "n": n,
                "N": n + 9,
                "shape": f"(({n + 9}^10))",
                "Q8_ball": endpoint.ball_text(q8),
                "Q9_left_ball": endpoint.ball_text(q9_left),
                "Q9_center_ball": endpoint.ball_text(q9_center),
                "Q9_right_ball": endpoint.ball_text(q9_right),
                "Q10_ball": endpoint.ball_text(q10),
                "deep_schur_ball": endpoint.ball_text(schur),
                "relative_Q9_margin_ball": endpoint.ball_text(relative),
                "relative_Q9_margin_lower": prefix.arb_lower_text(relative),
                "toda_lhs_ball": endpoint.ball_text(toda_lhs),
                "toda_rhs_ball": endpoint.ball_text(toda_rhs),
                "toda_residual_ball": endpoint.ball_text(toda_residual),
                "schur_residual_ball": endpoint.ball_text(schur_residual),
                "checks": checks,
            }
        )
    return rows


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    sources = validate_sources()
    values, source_diagnostics = prefix.merged_coefficients(
        source_paths=SOURCE_PATHS,
        max_coefficient_index=MAX_COEFFICIENT_INDEX,
    )
    splice_rows = direct_splice_rows(values)
    exact = {
        "orientation": "epsilon_10=(-1)^45=-1, Q_(10,n)=-H_(10,n)",
        "toda_identity": (
            "Q_(10,n)*Q_(8,n+2)=Q_(9,n+1)^2-"
            "Q_(9,n)*Q_(9,n+2)"
        ),
        "relative_margin": (
            "L_n=Q_(9,n+1)^2/(Q_(9,n)*Q_(9,n+2))-1"
        ),
        "splice": "Q_(10,n)(-100)>0 for n=1241,1242",
        "combined_positive_block": (
            "Q_(10,n)(-100)>0 for every 4<=n<=1242"
        ),
        "preserved_negative_prefix": (
            "Q_(10,n)(-100)<0 for n=0,1,2,3"
        ),
        "next_analytic_index": "n>=1243, equivalently k=n+9>=1252",
    }
    rows = [
        SpliceRow(
            "co10m100fsc_01_old_endpoint",
            "theorem_input",
            "ready_to_apply",
            "The established scan supplies the four negative shifts and the contiguous positive block through n=1240.",
            sources["old_endpoint"],
            "Previously certified finite endpoint rows only.",
        ),
        SpliceRow(
            "co10m100fsc_02_new_coefficients",
            "interval_input",
            "ready_to_apply",
            "Retained-integral Arb quadrature encloses the two coefficients needed by the next two order-ten determinants.",
            sources["new_coefficients"],
            "Two lambda=-100 coefficients only.",
            sources["new_coefficient_summary"],
        ),
        SpliceRow(
            "co10m100fsc_03_direct_determinants",
            "interval_certificate",
            "ready_to_apply",
            "Independent 4096-bit signed Hankel and Jacobi-Trudi determinants are strictly positive at both splice rows.",
            exact["splice"],
            "Direct finite Arb determinants only.",
            {"rows": splice_rows},
        ),
        SpliceRow(
            "co10m100fsc_04_toda_consistency",
            "interval_certificate",
            "ready_to_apply",
            "The positive Toda numerators and relative Q9 margins agree with the direct order-ten signs at both rows.",
            exact["toda_identity"] + "; " + exact["relative_margin"],
            "Interval consistency check for the exact condensation identity.",
        ),
        SpliceRow(
            "co10m100fsc_05_combined_block",
            "analytic_composition",
            "ready_to_apply",
            "The old scan and the two direct rows form one contiguous positive order-ten block.",
            exact["combined_positive_block"],
            "The four negative endpoint shifts remain excluded; no analytic tail is claimed.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_m100_finite_splice_certificate",
        "date": "2026-07-16",
        "status": "rigorous two-index order-ten finite splice at lambda=-100",
        "proof_boundary": (
            "This artifact proves the two new order-ten endpoint signs and the "
            "combined positive block 4<=n<=1242. It preserves the negative "
            "signs at n=0,1,2,3 and does not prove the analytic tail, delayed "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "source_diagnostics": source_diagnostics,
        "finite": {
            "lambda": "-100",
            "coefficient_range": [0, MAX_COEFFICIENT_INDEX],
            "precision_bits": PRECISION_BITS,
            "splice_indices": list(SPLICE_INDICES),
            "splice_rows": splice_rows,
            "combined_positive_range": [4, 1242],
            "preserved_negative_indices": [0, 1, 2, 3],
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "coefficients": len(values),
            "new_coefficient_rows": 2,
            "direct_Q10_checks": len(splice_rows),
            "direct_schur_checks": len(splice_rows),
            "toda_consistency_checks": len(splice_rows),
            "finite_splice_rows": len(splice_rows),
            "combined_positive_Q10_rows": 1239,
            "preserved_negative_endpoint_shifts": 4,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_m100_finite_splice_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_m100_finite_splice_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = artifact["finite"]
    lines = [
        "# Order-Ten Lambda=-100 Finite Splice Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous two-index endpoint splice. This is not a proof of",
        "the analytic order-ten tail, delayed entry, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order10_m100_finite_splice_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order10_m100_finite_splice_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_m100_finite_splice_certificate.py",
        "```",
        "",
        "## Added Coefficients",
        "",
        "Retained-integral quadrature with `n_sum=70`, cutoff 7, 220 decimal",
        "digits, and absolute tolerance `1e-150` gives rigorous balls for",
        "`A_1259(-100)` and `A_1260(-100)`.",
        "",
        "## Direct Signs",
        "",
        "At 4096-bit precision, signed Hankel determinants, Jacobi-Trudi",
        "determinants, Toda numerators, and relative `Q_9` margins all give",
        "the same strict signs:",
        "",
        "```text",
    ]
    for row in finite["splice_rows"]:
        lines.append(
            f"n={row['n']}: L_n={row['relative_Q9_margin_ball']}; Q_(10,n)>0"
        )
    lines.extend(
        [
            "```",
            "",
            "## Splice",
            "",
            "```text",
            exact["splice"],
            exact["combined_positive_block"],
            exact["preserved_negative_prefix"],
            exact["next_analytic_index"],
            "```",
            "",
            "The finite theorem now reaches two rows beyond the old scan. The",
            "remaining endpoint task is a full-kernel continuum curvature bound",
            "starting at `n=1243`; no sign is inferred there by this artifact.",
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
    summary = artifact["summary"]
    print(
        "wrote order-ten finite splice: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['finite_splice_rows']} direct splice rows, "
        f"{summary['combined_positive_Q10_rows']} combined positive signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
