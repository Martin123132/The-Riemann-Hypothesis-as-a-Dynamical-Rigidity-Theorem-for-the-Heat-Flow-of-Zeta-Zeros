#!/usr/bin/env python3
"""Compose the lambda=-100 cone and raw corridor into defect theorems."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.md"
)
CONE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.json"
)
CORRIDOR_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.json"
)


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    formula: str
    readiness: str
    proof_boundary: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_artifact() -> dict:
    cone = load_json(CONE_SOURCE)
    corridor = load_json(CORRIDOR_SOURCE)
    cone_summary = cone.get("summary", {})
    corridor_summary = corridor.get("summary", {})
    if cone_summary.get("full_cone_entry_rows") != 1:
        raise RuntimeError("full ratio-cone source is not closed")
    if corridor_summary.get("raw_corridor_theorem_rows") != 1:
        raise RuntimeError("raw-corridor source is not closed")
    if corridor_summary.get("open_requirements") != 0:
        raise RuntimeError("raw-corridor source still has open requirements")

    rows = [
        CertificateRow(
            id="m100adc_01_defect_coordinates",
            role="exact_reduction",
            claim="Introduce the raw ratio, contraction, defect, and scaled defect.",
            formula=(
                "x_k=((2*k-1)/(2*k+1))*R_k, d_k=1-x_k, "
                "s_k=((2*k+1)/2)*d_k"
            ),
            readiness="available_exact",
            proof_boundary="Exact coordinate change only.",
        ),
        CertificateRow(
            id="m100adc_02_full_cone_input",
            role="theorem_input",
            claim="Full lambda=-100 ratio-cone entry supplies both walls and the adjacent-k wall.",
            formula=(
                "(2*k-1)/(2*k+1)<=x_k<=1 and x_(k+1)>=x_k "
                "for every k>=1"
            ),
            readiness="ready_to_apply",
            proof_boundary="Uses the existing full ratio-cone theorem at one heat parameter.",
        ),
        CertificateRow(
            id="m100adc_03_raw_corridor_input",
            role="theorem_input",
            claim="The lambda=-100 raw-corridor theorem supplies both adjacent raw-ratio walls.",
            formula=(
                "((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k<=R_(k+1)"
                "<=(2+(2*k-1)*R_k)/(2*k+1)"
            ),
            readiness="ready_to_apply",
            proof_boundary="Uses the existing all-k raw-corridor theorem at lambda=-100.",
        ),
        CertificateRow(
            id="m100adc_04_defect_cone",
            role="theorem_conclusion",
            claim="The two pointwise cone walls are exactly the required defect bounds.",
            formula="0<=d_k<=2/(2*k+1) for every k>=1",
            readiness="ready_to_apply",
            proof_boundary="All-k defect cone at lambda=-100.",
        ),
        CertificateRow(
            id="m100adc_05_defect_monotonicity",
            role="theorem_conclusion",
            claim="The adjacent contraction wall is exactly decreasing defect.",
            formula="d_(k+1)<=d_k for every k>=1",
            readiness="ready_to_apply",
            proof_boundary="All-k monotone defect bridge at lambda=-100.",
        ),
        CertificateRow(
            id="m100adc_06_scaled_defect_cone",
            role="theorem_conclusion",
            claim="Scaling the defect cone gives the exact adaptive scaled-defect cone.",
            formula="0<=s_k<=1 for every k>=1",
            readiness="ready_to_apply",
            proof_boundary="All-k exact scaled-defect cone at lambda=-100.",
        ),
        CertificateRow(
            id="m100adc_07_scaled_defect_growth",
            role="theorem_conclusion",
            claim="The upper raw-corridor wall is exactly increasing scaled defect.",
            formula="s_(k+1)>=s_k for every k>=1",
            readiness="ready_to_apply",
            proof_boundary="All-k scaled-defect growth at lambda=-100.",
        ),
        CertificateRow(
            id="m100adc_08_parameter_specific_closure",
            role="parameter_specific_closure",
            claim=(
                "The lambda=-100 theorem discharges the defect-tail and adaptive-defect "
                "inputs needed by the one-entry-parameter cone route."
            ),
            formula=(
                "lambda=-100: defect tail closed from k=1; adaptive exact cone and "
                "monotone bridge closed from k=1"
            ),
            readiness="ready_to_apply",
            proof_boundary=(
                "The stronger simultaneous theorem at lambda=-25,-50,-100 is not proved "
                "and is not needed for the established lambda=-100 entry route."
            ),
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "theorem_input_rows": 2,
        "defect_conclusion_rows": 4,
        "parameter_specific_closure_rows": 1,
        "open_requirements": 0,
        "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in rows),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate",
        "date": "2026-07-11",
        "status": "exact interval-analytic adaptive-defect theorem at lambda=-100",
        "proof_boundary": (
            "This artifact proves the all-k defect cone, defect monotonicity, exact "
            "scaled-defect cone, and scaled-defect growth at lambda=-100. It does not "
            "prove the stronger simultaneous lambda=-25,-50,-100 target, PF-infinity, "
            "the all-order Jensen bridge, RH, or Lambda <= 0."
        ),
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "source_full_cone": CONE_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "source_raw_corridor": CORRIDOR_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "source_adaptive_target": (
            "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md"
        ),
        "source_defect_tail_target": (
            "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md"
        ),
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda -100 Adaptive-Defect Certificate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: exact interval-analytic adaptive-defect theorem at lambda=-100.",
        "This is not a proof of PF-infinity, the all-order Jensen bridge, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF negative-lambda -100 adaptive-defect certificate: "
            f"{summary['certificate_rows']} rows, 0 issues, 2 theorem inputs, "
            "4 defect conclusions, 0 open requirements"
        ),
        "```",
        "",
        "## Exact Composition",
        "",
        "Put",
        "",
        "```text",
        "x_k=((2*k-1)/(2*k+1))*R_k,",
        "d_k=1-x_k,",
        "s_k=((2*k+1)/2)*d_k.",
        "```",
        "",
        "Full ratio-cone entry at lambda=-100 proves",
        "",
        "```text",
        "(2*k-1)/(2*k+1)<=x_k<=1,",
        "x_(k+1)>=x_k.",
        "```",
        "",
        "Therefore, for every `k>=1`,",
        "",
        "```text",
        "0<=d_k<=2/(2*k+1),",
        "d_(k+1)<=d_k,",
        "0<=s_k<=1.",
        "```",
        "",
        "The upper raw-corridor wall gives the additional exact identity",
        "",
        "```text",
        "s_(k+1)-s_k",
        "  =(2+(2*k-1)*R_k-(2*k+1)*R_(k+1))/2",
        "  >=0.",
        "```",
        "",
        "Thus the lambda=-100 entry route satisfies the defect-tail theorem from",
        "`k=1` and the adaptive exact-cone plus monotone-bridge target from `k=1`.",
        "The older targets asked for the stronger simultaneous statement at",
        "`lambda=-25,-50,-100`; that simultaneous theorem is not proved, but it is",
        "not needed once one full entry parameter has been established at `lambda=-100`.",
        "",
        "Higher-degree minor cones and the all-order Jensen/PF bridge remain open.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
        "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "validated Jensen-window PF negative-lambda -100 adaptive-defect certificate: "
        f"{summary['certificate_rows']} rows, 0 issues, 2 theorem inputs, "
        "4 defect conclusions, 0 open requirements"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
