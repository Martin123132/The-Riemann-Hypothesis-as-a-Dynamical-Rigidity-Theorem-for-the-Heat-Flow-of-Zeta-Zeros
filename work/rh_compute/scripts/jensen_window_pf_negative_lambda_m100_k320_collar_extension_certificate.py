#!/usr/bin/env python3
"""Promote the repaired lambda=-100 A_245..A_320 source into a cone collar."""

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
    / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl"
)
DEFAULT_SOURCE_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220_summary.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md"
)
DEFAULT_PRECISION_BITS = 1024


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_text(value: flint.arb, digits: int = 60) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def load_coefficients(path: Path) -> dict[int, flint.arb]:
    values: dict[int, flint.arb] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row = json.loads(raw)
            if row.get("kind") != "acb_coefficient_enclosure":
                continue
            if row.get("lam") != "-100.0":
                raise ValueError(f"unexpected lambda: {row.get('lam')}")
            values[int(row["k"])] = flint.arb(row["A_ball"])
    if set(values) != set(range(245, 321)):
        raise ValueError(f"source index range differs: {min(values)}..{max(values)}, count={len(values)}")
    return values


def build_artifact(source_path: Path = DEFAULT_SOURCE_JSONL) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    values = load_coefficients(source_path)
    if not all(bool(value > 0 and not value.contains(0)) for value in values.values()):
        raise RuntimeError("source coefficient positivity failed")

    contractions: dict[int, flint.arb] = {}
    cone_rows = []
    for k in range(246, 320):
        x = values[k + 1] * values[k - 1] / values[k] ** 2
        lower_margin = x - flint.arb(2 * k - 1) / (2 * k + 1)
        upper_margin = 1 - x
        if not bool(lower_margin > 0 and upper_margin > 0):
            raise RuntimeError(f"cone wall failed at k={k}")
        contractions[k] = x
        cone_rows.append(
            {
                "k": k,
                "x_ball": arb_text(x),
                "lower_wall_margin_lower": arb_lower_text(lower_margin),
                "upper_wall_margin_lower": arb_lower_text(upper_margin),
            }
        )

    wall_rows = []
    for k in range(246, 319):
        gap = contractions[k + 1] - contractions[k]
        log_gap = (contractions[k + 1] / contractions[k]).log()
        if not bool(gap > 0 and log_gap > 0):
            raise RuntimeError(f"adjacent wall failed at k={k}")
        wall_rows.append(
            {
                "k": k,
                "gap_ball": arb_text(gap),
                "gap_lower": arb_lower_text(gap),
                "log_gap_ball": arb_text(log_gap),
                "log_gap_lower": arb_lower_text(log_gap),
                "new_extension_row": k >= 300,
            }
        )

    extension_rows = [row for row in wall_rows if row["new_extension_row"]]
    min_extension = min(extension_rows, key=lambda row: Decimal(row["log_gap_lower"]))
    min_all = min(wall_rows, key=lambda row: Decimal(row["log_gap_lower"]))
    diagnostics = {
        "lambda": "-100",
        "coefficient_range": [245, 320],
        "coefficient_rows": len(values),
        "cone_k_range": [246, 319],
        "cone_rows": len(cone_rows),
        "adjacent_wall_k_range": [246, 318],
        "adjacent_wall_rows": len(wall_rows),
        "new_extension_k_range": [300, 318],
        "new_extension_rows": len(extension_rows),
        "minimum_extension_log_gap_lower": min_extension["log_gap_lower"],
        "minimum_extension_log_gap_at_k": min_extension["k"],
        "minimum_all_log_gap_lower": min_all["log_gap_lower"],
        "minimum_all_log_gap_at_k": min_all["k"],
        "selected_extension_rows": [
            row for row in extension_rows if row["k"] in {300, 301, 310, 317, 318}
        ],
        "cone_detail_rows": cone_rows,
        "wall_detail_rows": wall_rows,
        "all_source_coefficients_positive": True,
        "all_cone_rows_certified": True,
        "all_adjacent_wall_rows_certified": True,
    }
    rows = [
        CertificateRow(
            id="m100k320_01_repaired_source",
            role="interval_input",
            readiness="finite_validated",
            claim="Import the repaired dps220 A_245..A_320 coefficient enclosures at lambda=-100.",
            proof_boundary="One finite coefficient source at one heat parameter.",
            diagnostics={"coefficient_rows": len(values), "coefficient_range": [245, 320]},
        ),
        CertificateRow(
            id="m100k320_02_coefficient_positivity",
            role="interval_certificate",
            readiness="finite_validated",
            claim="All 76 source coefficient balls are strictly positive.",
            proof_boundary="Finite positivity only.",
        ),
        CertificateRow(
            id="m100k320_03_cone_walls",
            role="interval_certificate",
            readiness="finite_validated",
            claim="The lower and upper ratio-cone walls hold for x_246 through x_319.",
            proof_boundary="74 finite contraction rows only; the all-k upper wall is a separate theorem.",
            diagnostics={"cone_rows": len(cone_rows), "cone_k_range": [246, 319]},
        ),
        CertificateRow(
            id="m100k320_04_adjacent_wall",
            role="interval_certificate",
            readiness="finite_validated",
            claim="Arb certifies x_(k+1)>x_k for every k=246..318.",
            proof_boundary="73 finite adjacent rows only; no eventual-k theorem.",
            diagnostics={"adjacent_wall_rows": len(wall_rows), "adjacent_wall_k_range": [246, 318]},
        ),
        CertificateRow(
            id="m100k320_05_new_collar_extension",
            role="finite_handoff",
            readiness="finite_validated",
            claim="The previously unused source rows extend the certified lambda=-100 adjacent-wall collar through k=318.",
            proof_boundary="Finite extension k=300..318 only.",
            diagnostics={
                "new_extension_rows": len(extension_rows),
                "minimum_extension_log_gap_lower": min_extension["log_gap_lower"],
                "minimum_extension_log_gap_at_k": min_extension["k"],
            },
        ),
        CertificateRow(
            id="m100k320_06_tail_handoff",
            role="non_promotion_gate",
            readiness="not_ready_to_apply",
            claim="Cone entry still needs an analytic adjacent-wall theorem from k=319 onward, together with the first-summand perturbation transfer.",
            proof_boundary="No promotion from a finite collar to cone entry, RH, or Lambda <= 0.",
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "coefficient_rows": len(values),
        "cone_rows": len(cone_rows),
        "adjacent_wall_rows": len(wall_rows),
        "new_extension_rows": len(extension_rows),
        "new_collar_end_k": 318,
        "minimum_extension_log_gap_lower": min_extension["log_gap_lower"],
        "minimum_extension_log_gap_at_k": min_extension["k"],
        "ready_to_apply_rows": 0,
        "main_finding": (
            "The existing repaired dps220 lambda=-100 source contains more rigorous information "
            "than the k300 audit exposed. Arb certifies 76 positive coefficients, both cone walls "
            "for x_246..x_319, and the adjacent wall for k=246..318. The 19 newly promoted rows "
            "k=300..318 have minimum adjacent log gap above 3.709e-6 at k=318, extending the "
            "finite handoff while leaving the analytic k>=319 tail open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate",
        "status": "finite Arb cone-collar extension certificate",
        "date": "2026-07-10",
        "proof_boundary": (
            "This finite Arb certificate extends the lambda=-100 adjacent-wall collar through "
            "k=318. It does not prove an eventual-k wall, full cone entry, RH, or Lambda <= 0."
        ),
        "source_enclosure_jsonl": source_path.relative_to(REPO_ROOT).as_posix(),
        "source_enclosure_summary": DEFAULT_SOURCE_SUMMARY.relative_to(REPO_ROOT).as_posix(),
        "source_k300_audit": "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "source_dominance_certificate": "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py",
        "diagnostics": diagnostics,
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Negative-Lambda -100 k320 Collar Extension Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: finite Arb cone-collar extension certificate. This is not a proof",
        "of an eventual-k theorem, cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate`.",
        "",
        "Machine-readable result and source:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.json",
        artifact["source_enclosure_jsonl"],
        artifact["source_enclosure_summary"],
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF negative-lambda -100 k320 collar extension certificate: 6 rows, 0 issues, 76 positive coefficients, 74 cone rows, 73 adjacent-wall rows, 19 new extension rows, 0 ready-to-apply rows",
        "```",
        "",
        "## Certified Extension",
        "",
        "The repaired dps220 source already covered `A_245..A_320`. Direct Arb",
        "ratio arithmetic certifies",
        "",
        "```text",
        "(2*k-1)/(2*k+1) < x_k < 1,       k=246..319",
        "x_(k+1)>x_k,                      k=246..318.",
        "```",
        "",
        "The second line promotes 19 previously unused adjacent rows",
        "`k=300..318`. The weakest new logarithmic gap is",
        "",
        "```text",
        f"L_318 > {diagnostics['minimum_extension_log_gap_lower']}",
        "```",
        "",
        "where `L_k=log(x_(k+1)/x_k)`.",
        "",
        "Selected extension rows:",
        "",
        "| k | Arb log gap |",
        "|---:|---:|",
    ]
    for row in diagnostics["selected_extension_rows"]:
        lines.append(f"| `{row['k']}` | `{row['log_gap_ball']}` |")
    lines.extend(
        [
            "",
            "## Handoff",
            "",
            "The finite lambda=-100 collar now reaches `k=318`. The analytic tail",
            "must start at `k=319`; the all-k first-summand dominance certificate",
            "already controls the full-kernel perturbation there, but the dominant",
            "`n=1` adjacent-wall lower bound remains open.",
            "",
            "```text",
            "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
            "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
            "outputs/signed_hankel_jensen_dependency_graph.md",
            "```",
            "",
            "Summary:",
            "",
            summary["main_finding"],
        ]
    )
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
        "validated Jensen-window PF negative-lambda -100 k320 collar extension certificate: "
        "6 rows, 0 issues, 76 positive coefficients, 74 cone rows, 73 adjacent-wall rows, "
        "19 new extension rows, 0 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
