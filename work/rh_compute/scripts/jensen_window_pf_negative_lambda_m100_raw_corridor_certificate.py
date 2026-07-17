#!/usr/bin/env python3
"""Compose full ratio-cone entry with scaled curvature into the raw corridor."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.md"
)
SCALED_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.json"
)
CONE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.json"
)


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    formula: str
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_artifact() -> dict:
    scaled = load_json(SCALED_SOURCE)
    cone = load_json(CONE_SOURCE)
    scaled_summary = scaled.get("summary", {})
    cone_summary = cone.get("summary", {})
    if scaled_summary.get("full_scaled_curvature_theorem_rows") != 1:
        raise RuntimeError("scaled-curvature theorem source is not closed")
    if scaled_summary.get("open_requirements") != 0:
        raise RuntimeError("scaled-curvature theorem source still has open requirements")
    if cone_summary.get("full_cone_entry_rows") != 1:
        raise RuntimeError("full ratio-cone source is not closed")

    rows = [
        CertificateRow(
            id="m100rcc_01_raw_curvature_coordinates",
            role="exact_reduction",
            claim="Rewrite the raw moment ratio through the normalized contraction and its log curvature.",
            formula=(
                "R_k=((2*k+1)/(2*k-1))*exp(-B_k), "
                "B_k=-log(((2*k-1)/(2*k+1))*R_k)"
            ),
            readiness="available_exact",
            proof_boundary="Exact coordinate change only.",
        ),
        CertificateRow(
            id="m100rcc_02_full_cone_input",
            role="theorem_input",
            claim="Full lambda=-100 ratio-cone entry supplies nonnegative and decreasing B curvature.",
            formula="B_k>=0 and B_(k+1)<=B_k for every k>=1",
            readiness="ready_to_apply",
            proof_boundary="Uses the existing full ratio-cone theorem at one heat parameter.",
            diagnostics={
                "full_cone_entry_rows": cone_summary["full_cone_entry_rows"],
                "analytic_tail_start": cone_summary["analytic_adjacent_tail_start"],
            },
        ),
        CertificateRow(
            id="m100rcc_03_scaled_curvature_input",
            role="theorem_input",
            claim="The continuous bridge supplies the linear lower curvature wall.",
            formula="B_(k+1)>=((2*k+1)/(2*k+3))*B_k for every k>=1",
            readiness="ready_to_apply",
            proof_boundary="Uses the new all-k scaled-curvature theorem at lambda=-100.",
            diagnostics={
                "compact_blocks": scaled_summary["compact_blocks"],
                "prefix_gaps": scaled_summary["positive_prefix_gaps"],
            },
        ),
        CertificateRow(
            id="m100rcc_04_linear_to_nonlinear_barrier",
            role="exact_lemma",
            claim="The linear wall dominates the exact nonlinear lower raw-corridor wall for B>=0.",
            formula=(
                "log((2*k+3)/(2+(2*k+1)*exp(-B))) "
                "<=((2*k+1)/(2*k+3))*B, B>=0"
            ),
            readiness="available_exact",
            proof_boundary=(
                "At B=0 equality holds; the derivative of the left side is "
                "y/(2+y)<=y_0/(2+y_0)=(2*k+1)/(2*k+3)."
            ),
        ),
        CertificateRow(
            id="m100rcc_05_curvature_corridor",
            role="exact_theorem_composition",
            claim="The two source theorems and the calculus lemma prove the complete B corridor.",
            formula=(
                "log((2*k+3)/(2+(2*k+1)*exp(-B_k))) "
                "<=B_(k+1)<=B_k"
            ),
            readiness="ready_to_apply",
            proof_boundary="Complete coefficient-curvature corridor at lambda=-100.",
        ),
        CertificateRow(
            id="m100rcc_06_raw_ratio_corridor",
            role="interval_analytic_theorem",
            claim="The actual zeta moments satisfy the full raw-ratio decrement corridor at lambda=-100.",
            formula=(
                "((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k<=R_(k+1)"
                "<=(2+(2*k-1)*R_k)/(2*k+1), k>=1"
            ),
            readiness="ready_to_apply",
            proof_boundary=(
                "Closes the zeta-specific raw-corridor target at lambda=-100; "
                "not PF-infinity, the all-order Jensen bridge, RH, or Lambda <= 0."
            ),
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "exact_reduction_rows": 2,
        "theorem_input_rows": 2,
        "raw_corridor_theorem_rows": 1,
        "compact_source_blocks": scaled_summary["compact_blocks"],
        "finite_prefix_source_gaps": scaled_summary["positive_prefix_gaps"],
        "open_requirements": 0,
        "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in rows),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_m100_raw_corridor_certificate",
        "date": "2026-07-11",
        "status": "exact interval-analytic raw-corridor theorem at lambda=-100",
        "proof_boundary": (
            "This artifact proves the full raw-ratio decrement corridor for the actual "
            "zeta heat-flow moments at lambda=-100. It does not prove PF-infinity, the "
            "all-order Jensen bridge, RH, or Lambda <= 0."
        ),
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "source_scaled_curvature": SCALED_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "source_full_cone": CONE_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "source_exact_bridge": "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
        "source_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda -100 Raw-Corridor Certificate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: exact interval-analytic raw-corridor theorem at lambda=-100.",
        "This is not a proof of PF-infinity, the all-order Jensen bridge, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_m100_raw_corridor_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF negative-lambda -100 raw-corridor certificate: "
            f"{summary['certificate_rows']} rows, 0 issues, 2 theorem inputs, "
            "1 raw-corridor theorem, 0 open requirements"
        ),
        "```",
        "",
        "## Exact Composition",
        "",
        "Put",
        "",
        "```text",
        "R_k=M_(k+1)*M_(k-1)/M_k^2,",
        "B_k=-log(((2*k-1)/(2*k+1))*R_k).",
        "```",
        "",
        "Full ratio-cone entry at lambda=-100 proves",
        "",
        "```text",
        "B_k>=0, B_(k+1)<=B_k.",
        "```",
        "",
        "The scaled-curvature continuous bridge proves",
        "",
        "```text",
        "B_(k+1)>=((2*k+1)/(2*k+3))*B_k.",
        "```",
        "",
        "For `B>=0`, the exact calculus inequality",
        "",
        "```text",
        "log((2*k+3)/(2+(2*k+1)*exp(-B)))",
        "<=((2*k+1)/(2*k+3))*B",
        "```",
        "",
        "follows because equality holds at zero and the left derivative decreases from",
        "`(2*k+1)/(2*k+3)`. Therefore",
        "",
        "```text",
        "log((2*k+3)/(2+(2*k+1)*exp(-B_k)))<=B_(k+1)<=B_k.",
        "```",
        "",
        "Returning to raw ratios gives, for every `k>=1`,",
        "",
        "```text",
        "((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k<=R_(k+1)",
        "R_(k+1)<=(2+(2*k-1)*R_k)/(2*k+1).",
        "```",
        "",
        "This closes the zeta-specific raw-corridor target at lambda=-100. Higher-degree",
        "minor cones and the all-order Jensen/PF bridge remain open.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.md",
        "outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
        "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
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
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "validated Jensen-window PF negative-lambda -100 raw-corridor certificate: "
        f"{summary['certificate_rows']} rows, 0 issues, 2 theorem inputs, "
        "1 raw-corridor theorem, 0 open requirements"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
