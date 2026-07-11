#!/usr/bin/env python3
"""Build an exact log-concave Mellin countermodel to adjacent-k monotonicity."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.md"
DEFAULT_PRECISION_BITS = 512


@dataclass(frozen=True)
class CountermodelRow:
    id: str
    role: str
    claim: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 70) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def arb_upper_text(value: flint.arb, digits: int = 70) -> str:
    rounded = value.upper().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_plus(), "E")


def normalized_mellin(p: flint.arb) -> flint.arb:
    five = flint.arb(5)
    return five.gamma_lower(p) / (five**p * p.gamma())


def contraction(k: int) -> flint.arb:
    p = flint.arb(k) + flint.arb("0.5")
    return normalized_mellin(p + 1) * normalized_mellin(p - 1) / normalized_mellin(p) ** 2


def build_artifact() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    x1 = contraction(1)
    x2 = contraction(2)
    gap = x2 - x1
    if not bool(x1 > 0 and x1 < 1 and x2 > 0 and x2 < 1 and gap < 0):
        raise RuntimeError("log-concave Mellin countermodel did not separate the required signs")
    rows = [
        CountermodelRow(
            id="lcmmc_01_log_concave_density",
            role="exact_input",
            claim="f(y)=exp(-5*y)*1_[0,1](y) is an integrable log-concave density on [0,infinity).",
            proof_boundary="Abstract density only; not the Riemann/Newman kernel.",
        ),
        CountermodelRow(
            id="lcmmc_02_normalized_mellin_closed_form",
            role="exact_reduction",
            claim="H(p)=integral_0^1 y^(p-1)*exp(-5*y)dy/Gamma(p)=gamma_lower(p,5)/(5^p*Gamma(p)).",
            proof_boundary="Exact incomplete-Gamma evaluation for p>0.",
        ),
        CountermodelRow(
            id="lcmmc_03_first_contraction",
            role="interval_evaluation",
            claim="The first Gamma-normalized Mellin contraction lies strictly between zero and one.",
            proof_boundary="One Arb-enclosed abstract contraction.",
            diagnostics={"x_1_ball": arb_text(x1), "x_1_upper": arb_upper_text(x1)},
        ),
        CountermodelRow(
            id="lcmmc_04_second_contraction",
            role="interval_evaluation",
            claim="The second Gamma-normalized Mellin contraction lies strictly between zero and one.",
            proof_boundary="One Arb-enclosed abstract contraction.",
            diagnostics={"x_2_ball": arb_text(x2), "x_2_upper": arb_upper_text(x2)},
        ),
        CountermodelRow(
            id="lcmmc_05_monotone_wall_violation",
            role="countermodel_gate",
            claim="The strict Arb inequality x_2<x_1 disproves any generic implication from log-concavity to adjacent-k contraction monotonicity.",
            proof_boundary="Exact abstract countermodel gate; not a failure of the actual zeta sequence.",
            diagnostics={
                "x_2_minus_x_1_ball": arb_text(gap),
                "x_2_minus_x_1_upper": arb_upper_text(gap),
            },
        ),
        CountermodelRow(
            id="lcmmc_06_zeta_specific_handoff",
            role="non_promotion_gate",
            claim="A proof of x_(k+1)>=x_k must use structure beyond kernel log-concavity and Berwald-Borell.",
            proof_boundary="No claim against a zeta-specific recurrence, determinant identity, or saddle theorem.",
        ),
    ]
    summary = {
        "countermodel_rows": len(rows),
        "log_concave_density_rows": 1,
        "upper_wall_rows": 2,
        "monotone_wall_violations": 1,
        "x_1_lower": arb_lower_text(x1),
        "x_1_upper": arb_upper_text(x1),
        "x_2_lower": arb_lower_text(x2),
        "x_2_upper": arb_upper_text(x2),
        "x_2_minus_x_1_upper": arb_upper_text(gap),
        "main_finding": (
            "The compactly supported log-concave density exp(-5y) on [0,1] has "
            "Gamma-normalized Mellin contractions x_1 and x_2 inside the Berwald-Borell "
            "upper wall but x_2<x_1. Therefore the remaining adjacent-k wall is not a "
            "generic consequence of log-concavity and needs zeta-specific structure."
        ),
    }
    return {
        "kind": "jensen_window_pf_log_concave_mellin_monotone_wall_countermodel",
        "status": "exact interval countermodel gate",
        "date": "2026-07-10",
        "proof_boundary": (
            "This exact Arb countermodel blocks promotion from generic log-concavity and "
            "Berwald-Borell to adjacent-k monotonicity. It is not the zeta kernel, does not "
            "disprove zeta monotone contractions, and does not prove or disprove RH or Lambda <= 0."
        ),
        "source_upper_wall_certificate": "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
        "source_monotone_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py",
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Log-Concave Mellin Monotone-Wall Countermodel",
        "",
        "Date: 2026-07-10",
        "",
        "Status: exact interval countermodel gate. This is not a proof or disproof",
        "of zeta monotone contractions, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_log_concave_mellin_monotone_wall_countermodel`.",
        "",
        "Machine-readable result:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF log-concave Mellin monotone-wall countermodel: 6 rows, 0 issues, 2 upper-wall contractions, 1 monotone-wall violation",
        "```",
        "",
        "## Exact Witness",
        "",
        "Take the log-concave density",
        "",
        "```text",
        "f(y)=exp(-5*y)*1_[0,1](y).",
        "```",
        "",
        "Its Gamma-normalized Mellin transform is",
        "",
        "```text",
        "H(p)=integral_0^1 y^(p-1)*exp(-5*y)dy/Gamma(p)",
        "    =gamma_lower(p,5)/(5^p*Gamma(p)).",
        "x_k=H(k+3/2)*H(k-1/2)/H(k+1/2)^2.",
        "```",
        "",
        "Arb certifies:",
        "",
        f"- `x_1` lies between `{summary['x_1_lower']}` and `{summary['x_1_upper']}`.",
        f"- `x_2` lies between `{summary['x_2_lower']}` and `{summary['x_2_upper']}`.",
        f"- `x_2-x_1 < {summary['x_2_minus_x_1_upper']}`.",
        "",
        "Both contractions satisfy the Berwald-Borell upper wall `0<x_k<1`,",
        "but `x_2<x_1`. Thus log-concavity of the measure, even with compact",
        "support and a log-affine density, does not imply `x_(k+1)>=x_k`.",
        "",
        "## Proof Boundary",
        "",
        "This witness blocks only a generic theorem. It does not touch the strict",
        "zeta-kernel log-concavity certificate and does not rule out a zeta-specific",
        "ratio recurrence, determinant identity, or uniform saddle theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
        "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF log-concave Mellin monotone-wall countermodel: "
        "6 rows, 0 issues, 2 upper-wall contractions, 1 monotone-wall violation"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
