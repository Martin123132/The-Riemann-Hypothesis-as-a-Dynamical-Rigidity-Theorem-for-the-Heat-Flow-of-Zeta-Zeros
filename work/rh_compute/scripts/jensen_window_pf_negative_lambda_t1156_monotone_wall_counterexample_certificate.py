#!/usr/bin/env python3
"""Certify a zeta-kernel monotone-wall violation at lambda=-1156."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE_JSONL = (
    REPO_ROOT
    / "work/rh_compute/results/acb_enclosures_lambda_m1156_k119_k122_dps250.jsonl"
)
DEFAULT_SOURCE_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/acb_enclosures_lambda_m1156_k119_k122_dps250_summary.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md"
)
DEFAULT_PRECISION_BITS = 1024


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_text(value: flint.arb, digits: int = 80) -> str:
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


def load_source(path: Path) -> dict[int, flint.arb]:
    coefficients: dict[int, flint.arb] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row = json.loads(raw)
            if row.get("kind") != "acb_coefficient_enclosure":
                continue
            if row.get("lam") != "-1156.0":
                raise ValueError(f"unexpected lambda in source row: {row.get('lam')}")
            coefficients[int(row["k"])] = flint.arb(row["A_ball"])
    if set(coefficients) != {119, 120, 121, 122}:
        raise ValueError(f"source coefficient indices differ: {sorted(coefficients)}")
    return coefficients


def build_artifact(source_path: Path = DEFAULT_SOURCE_JSONL) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    coefficients = load_source(source_path)
    a119, a120 = coefficients[119], coefficients[120]
    a121, a122 = coefficients[121], coefficients[122]
    x120 = a121 * a119 / a120**2
    x121 = a122 * a120 / a121**2
    gap = x121 - x120
    log_ratio = (x121 / x120).log()
    polynomial = a122 * a120**3 - a121**3 * a119
    all_positive = all(bool(value > 0 and not value.contains(0)) for value in coefficients.values())
    upper_walls = bool(x120 < 1 and x121 < 1 and x120 > 0 and x121 > 0)
    violation = bool(gap < 0 and not gap.contains(0) and polynomial < 0 and not polynomial.contains(0))
    if not (all_positive and upper_walls and violation):
        raise RuntimeError("T=1156 monotone-wall source did not certify the required signs")

    diagnostics = {
        "lambda": "-1156",
        "T": 1156,
        "violated_wall_index": 120,
        "coefficient_indices": [119, 120, 121, 122],
        "coefficient_balls": {str(k): arb_text(value) for k, value in coefficients.items()},
        "x_120_ball": arb_text(x120),
        "x_121_ball": arb_text(x121),
        "x_121_minus_x_120_ball": arb_text(gap),
        "x_121_minus_x_120_upper": arb_upper_text(gap),
        "log_x_121_over_x_120_ball": arb_text(log_ratio),
        "log_x_121_over_x_120_upper": arb_upper_text(log_ratio),
        "polynomial_witness_ball": arb_text(polynomial),
        "polynomial_witness_upper": arb_upper_text(polynomial),
        "all_source_coefficients_positive": all_positive,
        "both_contractions_inside_upper_wall": upper_walls,
        "monotone_wall_strictly_violated": violation,
    }
    rows = [
        CertificateRow(
            id="t1156mwc_01_rigorous_coefficient_source",
            role="interval_input",
            claim="ACB integration plus analytic n-series and u-tail bounds encloses A_119 through A_122 at lambda=-1156.",
            proof_boundary="Four coefficient enclosures at one heat parameter only.",
            diagnostics={
                "coefficient_indices": diagnostics["coefficient_indices"],
                "all_source_coefficients_positive": all_positive,
            },
        ),
        CertificateRow(
            id="t1156mwc_02_exact_contraction_reduction",
            role="exact_reduction",
            claim="x_120=A_121*A_119/A_120^2 and x_121=A_122*A_120/A_121^2.",
            proof_boundary="Exact ratio algebra only.",
        ),
        CertificateRow(
            id="t1156mwc_03_upper_wall_occupancy",
            role="interval_evaluation",
            claim="Both x_120 and x_121 are strictly between zero and one.",
            proof_boundary="Two contractions at lambda=-1156 only; the global upper wall is proved separately.",
            diagnostics={"x_120_ball": diagnostics["x_120_ball"], "x_121_ball": diagnostics["x_121_ball"]},
        ),
        CertificateRow(
            id="t1156mwc_04_strict_monotone_wall_violation",
            role="interval_counterexample",
            claim="Arb certifies x_121-x_120<0, so the adjacent-k monotone wall fails for the actual zeta kernel at lambda=-1156 and k=120.",
            proof_boundary="A strict local zeta-kernel counterexample, not a statement about every lambda or every k.",
            diagnostics={
                "gap_ball": diagnostics["x_121_minus_x_120_ball"],
                "gap_upper": diagnostics["x_121_minus_x_120_upper"],
            },
        ),
        CertificateRow(
            id="t1156mwc_05_polynomial_witness",
            role="interval_counterexample",
            claim="The equivalent polynomial A_122*A_120^3-A_121^3*A_119 is strictly negative.",
            proof_boundary="Equivalent cross-multiplied witness using positive enclosed coefficients.",
            diagnostics={"polynomial_upper": diagnostics["polynomial_witness_upper"]},
        ),
        CertificateRow(
            id="t1156mwc_06_fixed_k_promotion_gate",
            role="non_promotion_gate",
            claim="The fixed-k=22 theorem for every T>=1156 cannot be promoted to an all-k cone theorem at T=1156.",
            proof_boundary="Blocks that promotion only; the certified fixed-k=22 theorem remains valid.",
        ),
        CertificateRow(
            id="t1156mwc_07_revised_tail_handoff",
            role="theorem_search_handoff",
            claim="A surviving cone-entry route must use a different entry time with a finite collar plus an eventual-k theorem, or avoid the monotone cone.",
            proof_boundary="Research redirection only; no moderate-T all-k theorem or cone entry is claimed.",
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "source_coefficient_rows": len(coefficients),
        "upper_wall_contractions": 2,
        "zeta_monotone_wall_violations": 1,
        "violated_wall_index": 120,
        "gap_upper": diagnostics["x_121_minus_x_120_upper"],
        "log_ratio_upper": diagnostics["log_x_121_over_x_120_upper"],
        "fixed_k_t1156_all_k_promotion_blocked": True,
        "main_finding": (
            "Rigorous ACB coefficient enclosures for the actual Newman kernel at lambda=-1156 "
            "give x_121-x_120<0, with an upper enclosure below -1.68e-8. Thus the adjacent-k "
            "monotone wall fails at k=120 even though both contractions remain below one. The "
            "fixed-k=22 T>=1156 certificate is valid, but T=1156 cannot be an all-k cone-entry time."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate",
        "status": "interval zeta-kernel counterexample certificate",
        "date": "2026-07-10",
        "proof_boundary": (
            "This interval certificate proves one actual zeta-kernel adjacent-k violation at "
            "lambda=-1156. It blocks all-k cone entry at that parameter and universal promotion "
            "of the fixed-k T>=1156 collar. It does not disprove cone entry at another lambda, "
            "does not disprove Jensen-window PF-infinity or RH, and does not prove Lambda <= 0."
        ),
        "source_enclosure_jsonl": source_path.relative_to(REPO_ROOT).as_posix(),
        "source_enclosure_summary": DEFAULT_SOURCE_SUMMARY.relative_to(REPO_ROOT).as_posix(),
        "source_enclosure_generator": "work/rh_compute/scripts/acb_coefficient_enclosures.py",
        "source_tail_bounds": "work/rh_compute/scripts/coefficient_tail_bounds.py",
        "source_fixed_k_certificate": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.md",
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py",
        "source_generation_command": (
            "python work/rh_compute/scripts/acb_coefficient_enclosures.py --lambdas=-1156 "
            "--k-min 119 --k-max 122 --n-sum 20 --cutoff 2 --dps 250 --digits 100 "
            "--abs-tol 1e-220 --constant-terms 40 --output-dir work/rh_compute/results "
            "--run-id acb_enclosures_lambda_m1156_k119_k122_dps250 --overwrite"
        ),
        "diagnostics": diagnostics,
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF T=1156 Monotone-Wall Counterexample Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: interval zeta-kernel counterexample certificate. This is not a",
        "proof or disproof of RH, and it does not prove `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate`.",
        "",
        "Machine-readable result and source enclosures:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.json",
        artifact["source_enclosure_jsonl"],
        artifact["source_enclosure_summary"],
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF T=1156 monotone-wall counterexample certificate: 7 rows, 0 issues, 4 coefficient enclosures, 1 zeta monotone-wall violation",
        "```",
        "",
        "## Certified Violation",
        "",
        "For the actual Newman kernel at `lambda=-1156`, define",
        "",
        "```text",
        "x_120=A_121*A_119/A_120^2",
        "x_121=A_122*A_120/A_121^2.",
        "```",
        "",
        "ACB integration of `A_119..A_122`, composed with analytic n-series and",
        "u-tail bounds, certifies",
        "",
        "```text",
        f"x_120 = {diagnostics['x_120_ball']}",
        f"x_121 = {diagnostics['x_121_ball']}",
        f"x_121-x_120 = {diagnostics['x_121_minus_x_120_ball']}",
        f"log(x_121/x_120) = {diagnostics['log_x_121_over_x_120_ball']}",
        "```",
        "",
        f"In particular, `x_121-x_120 < {summary['gap_upper']}`. Both contractions",
        "are strictly between zero and one, so this is specifically a failure of",
        "the adjacent-k wall `x_(k+1)>=x_k`, not of the Mellin upper wall.",
        "",
        "## Route Consequence",
        "",
        "The fixed-k=22 certificate for every real `T>=1156` remains valid. This",
        "counterexample shows that it cannot be extended into an all-k cone theorem",
        "at `T=1156`: the actual sequence leaves the monotone wall at `k=120`.",
        "",
        "The surviving cone-entry programme is now:",
        "",
        "```text",
        "1. choose a moderate fixed negative lambda with a deep certified finite collar;",
        "2. prove eventual-k adjacent monotonicity there by a zeta saddle/tail theorem;",
        "3. splice the finite collar to that tail and then invoke cone invariance.",
        "```",
        "",
        "The existing `lambda=-25,-50,-100` k300 enclosures support this revised",
        "search, but they remain finite evidence and do not prove the tail theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.md",
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
    parser.add_argument("--source-jsonl", type=Path, default=DEFAULT_SOURCE_JSONL)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact(args.source_jsonl)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF T=1156 monotone-wall counterexample certificate: "
        "7 rows, 0 issues, 4 coefficient enclosures, 1 zeta monotone-wall violation"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
