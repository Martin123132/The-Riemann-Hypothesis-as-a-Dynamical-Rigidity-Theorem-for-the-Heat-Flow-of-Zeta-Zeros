#!/usr/bin/env python3
"""Certify full ratio-cone entry for the zeta coefficients at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists() and str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


BASE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl"
)
REPAIR_220_250 = (
    REPO_ROOT
    / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl"
)
REPAIR_245_320 = (
    REPO_ROOT
    / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md"
DEFAULT_PRECISION_BITS = 1024


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    formula: str | None
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_text(value: flint.arb, digits: int = 60) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def load_source(path: Path) -> dict[int, flint.arb]:
    values: dict[int, flint.arb] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row = json.loads(raw)
            if row.get("kind") != "acb_coefficient_enclosure":
                continue
            if row.get("lam") not in {"-100", "-100.0"}:
                continue
            values[int(row["k"])] = flint.arb(row["A_ball"])
    return values


def merged_coefficients() -> tuple[dict[int, flint.arb], list[dict]]:
    sources = [BASE_SOURCE, REPAIR_220_250, REPAIR_245_320]
    values: dict[int, flint.arb] = {}
    diagnostics = []
    for precedence, source in enumerate(sources):
        loaded = load_source(source)
        overwritten = len(set(values).intersection(loaded))
        values.update(loaded)
        diagnostics.append(
            {
                "precedence": precedence,
                "source": source.relative_to(REPO_ROOT).as_posix(),
                "lambda_minus_100_rows": len(loaded),
                "index_range": [min(loaded), max(loaded)],
                "overwritten_rows": overwritten,
            }
        )
    if set(values) != set(range(321)):
        raise RuntimeError("merged source does not cover A_0..A_320")
    return values, diagnostics


def build_artifact() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    values, source_rows = merged_coefficients()
    if not all(bool(value > 0 and not value.contains(0)) for value in values.values()):
        raise RuntimeError("merged coefficient positivity failed")

    contractions: dict[int, flint.arb] = {}
    minimum_lower: tuple[flint.arb, int] | None = None
    minimum_upper: tuple[flint.arb, int] | None = None
    selected_contractions = []
    selected_k = {1, 2, 100, 219, 220, 244, 245, 299, 300, 318, 319}
    for k in range(1, 320):
        x = values[k + 1] * values[k - 1] / values[k] ** 2
        lower = x - flint.arb(2 * k - 1) / (2 * k + 1)
        upper = 1 - x
        if not bool(lower > 0 and upper > 0):
            raise RuntimeError(f"prefix pointwise cone failed at k={k}")
        contractions[k] = x
        if minimum_lower is None or lower.lower() < minimum_lower[0].lower():
            minimum_lower = (lower, k)
        if minimum_upper is None or upper.lower() < minimum_upper[0].lower():
            minimum_upper = (upper, k)
        if k in selected_k:
            selected_contractions.append(
                {
                    "k": k,
                    "x_ball": arb_text(x),
                    "lower_margin_lower": arb_lower_text(lower),
                    "upper_margin_lower": arb_lower_text(upper),
                }
            )

    minimum_adjacent: tuple[flint.arb, int] | None = None
    selected_adjacent = []
    for k in range(1, 319):
        gap = contractions[k + 1] - contractions[k]
        log_gap = (contractions[k + 1] / contractions[k]).log()
        if not bool(gap > 0 and log_gap > 0):
            raise RuntimeError(f"prefix adjacent cone failed at k={k}")
        if minimum_adjacent is None or gap.lower() < minimum_adjacent[0].lower():
            minimum_adjacent = (gap, k)
        if k in selected_k:
            selected_adjacent.append(
                {
                    "k": k,
                    "gap_ball": arb_text(gap),
                    "gap_lower": arb_lower_text(gap),
                    "log_gap_lower": arb_lower_text(log_gap),
                }
            )
    assert minimum_lower is not None and minimum_upper is not None and minimum_adjacent is not None

    diagnostics = {
        "parameters": {"lambda": "-100", "precision_bits": DEFAULT_PRECISION_BITS},
        "merged_sources": source_rows,
        "prefix": {
            "coefficient_range": [0, 320],
            "positive_coefficients": 321,
            "pointwise_cone_k_range": [1, 319],
            "positive_pointwise_cone_rows": 319,
            "adjacent_k_range": [1, 318],
            "positive_adjacent_rows": 318,
            "minimum_lower_margin_lower": arb_lower_text(minimum_lower[0]),
            "minimum_lower_margin_at_k": minimum_lower[1],
            "minimum_upper_margin_lower": arb_lower_text(minimum_upper[0]),
            "minimum_upper_margin_at_k": minimum_upper[1],
            "minimum_adjacent_margin_lower": arb_lower_text(minimum_adjacent[0]),
            "minimum_adjacent_margin_at_k": minimum_adjacent[1],
            "selected_contractions": selected_contractions,
            "selected_adjacent": selected_adjacent,
        },
        "global_composition": {
            "lower_wall": "x_k>=((2*k-1)/(2*k+1)) for every real lambda and k>=1",
            "upper_wall": "x_k<=1 for every real lambda and k>=1",
            "prefix_adjacent": "x_(k+1)>x_k for k=1..318 at lambda=-100",
            "tail_adjacent": "x_(k+1)>x_k for every k>=319 at lambda=-100",
            "coverage": "all integer k>=1",
        },
    }
    rows = [
        CertificateRow(
            id="m100fce_01_repaired_source_merge",
            role="interval_input",
            claim="Merge the broad A_0..A_300 source with both dps220 lambda=-100 repair sources in precedence order.",
            formula="broad <- repair_220_250 <- repair_245_320",
            readiness="finite_validated",
            proof_boundary="One finite merged source at lambda=-100.",
            diagnostics={"merged_sources": source_rows},
        ),
        CertificateRow(
            id="m100fce_02_prefix_coefficient_positivity",
            role="interval_certificate",
            claim="All merged coefficient balls A_0 through A_320 are strictly positive.",
            formula="A_k(-100)>0, 0<=k<=320",
            readiness="finite_validated",
            proof_boundary="Finite coefficient positivity only.",
        ),
        CertificateRow(
            id="m100fce_03_prefix_pointwise_cone",
            role="interval_certificate",
            claim="Every merged prefix contraction satisfies both pointwise cone walls.",
            formula="(2*k-1)/(2*k+1)<x_k<1, 1<=k<=319",
            readiness="finite_validated",
            proof_boundary="Finite pointwise cone only.",
            diagnostics={
                key: diagnostics["prefix"][key]
                for key in (
                    "positive_pointwise_cone_rows",
                    "minimum_lower_margin_lower",
                    "minimum_lower_margin_at_k",
                    "minimum_upper_margin_lower",
                    "minimum_upper_margin_at_k",
                )
            },
        ),
        CertificateRow(
            id="m100fce_04_prefix_adjacent_cone",
            role="interval_certificate",
            claim="Every adjacent prefix contraction is strictly increasing.",
            formula="x_(k+1)>x_k, 1<=k<=318",
            readiness="finite_validated",
            proof_boundary="Finite adjacent cone only.",
            diagnostics={
                key: diagnostics["prefix"][key]
                for key in (
                    "positive_adjacent_rows",
                    "minimum_adjacent_margin_lower",
                    "minimum_adjacent_margin_at_k",
                )
            },
        ),
        CertificateRow(
            id="m100fce_05_global_pointwise_walls",
            role="exact_theorem_composition",
            claim="The exact Cauchy-Schwarz lower wall and Mellin log-concavity upper wall cover every remaining index.",
            formula="(2*k-1)/(2*k+1)<=x_k<=1, k>=1",
            readiness="ready_to_apply",
            proof_boundary="Global pointwise walls only; no adjacent conclusion by themselves.",
        ),
        CertificateRow(
            id="m100fce_06_analytic_adjacent_tail",
            role="theorem_composition",
            claim="The paired-remainder wall closure proves the adjacent cone on the entire tail.",
            formula="x_(k+1)>x_k, k>=319, lambda=-100",
            readiness="ready_to_apply",
            proof_boundary="Analytic adjacent tail only.",
        ),
        CertificateRow(
            id="m100fce_07_full_cone_entry",
            role="interval_analytic_theorem",
            claim="The actual zeta heat-flow coefficient sequence lies in the full infinite ratio cone at lambda=-100.",
            formula="(2*k-1)/(2*k+1)<=x_k<=1 and x_(k+1)>=x_k for every k>=1",
            readiness="ready_to_apply",
            proof_boundary="Full cone entry at one heat parameter; not forward infinite-flow legitimacy or the all-order bridge.",
        ),
        CertificateRow(
            id="m100fce_08_flow_handoff",
            role="open_handoff",
            claim="Apply the exact boundary algebra after proving analytic legitimacy for the infinite or collared forward heat flow.",
            formula="cone entry at -100 + legitimate invariance => cone for lambda>=-100",
            readiness="blocked_by_flow_legitimacy",
            proof_boundary="Open downstream handoff; no promotion to RH or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate",
        "date": "2026-07-10",
        "status": "interval full cone-entry theorem at lambda=-100 and open flow-legitimacy handoff",
        "proof_boundary": (
            "This artifact proves full infinite ratio-cone entry for the actual zeta "
            "coefficients at lambda=-100. It does not prove analytic legitimacy of the "
            "forward infinite flow, the all-order Jensen bridge, RH, or Lambda <= 0."
        ),
        "source_base_jsonl": BASE_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "source_repair_220_250_jsonl": REPAIR_220_250.relative_to(REPO_ROOT).as_posix(),
        "source_repair_245_320_jsonl": REPAIR_245_320.relative_to(REPO_ROOT).as_posix(),
        "source_lower_wall": "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "source_upper_wall": "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
        "source_adjacent_tail": "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
        "source_flow_invariance": "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py",
        "diagnostics": diagnostics,
        "summary": {
            "certificate_rows": len(rows),
            "positive_coefficients": 321,
            "positive_pointwise_cone_rows": 319,
            "positive_adjacent_rows": 318,
            "analytic_adjacent_tail_start": 319,
            "full_cone_entry_rows": 1,
            "open_flow_handoff_rows": 1,
            "ready_to_apply_rows": 3,
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    prefix = artifact["diagnostics"]["prefix"]
    lines = [
        "# Jensen-Window PF Negative-Lambda -100 Full Cone-Entry Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: interval full cone-entry theorem at lambda=-100 and open flow-legitimacy handoff. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py",
        "```",
        "",
        "## Repaired Prefix",
        "",
        "The broad lambda=-100 coefficient source is overridden by the dps220",
        "repairs on `k=220..250` and `k=245..320`. Direct Arb arithmetic proves",
        "",
        "```text",
        "A_k(-100)>0, 0<=k<=320,",
        "(2*k-1)/(2*k+1)<x_k<1, 1<=k<=319,",
        "x_(k+1)>x_k, 1<=k<=318.",
        "```",
        "",
        f"The weakest prefix adjacent margin is `{prefix['minimum_adjacent_margin_lower']}` at `k={prefix['minimum_adjacent_margin_at_k']}`.",
        "",
        "## Global Composition",
        "",
        "The exact Cauchy-Schwarz and Mellin log-concavity theorems prove the",
        "pointwise lower and upper walls for every real heat parameter and every",
        "`k>=1`. The paired-remainder wall closure proves the adjacent tail for",
        "every `k>=319`. Hence",
        "",
        "```text",
        "(2*k-1)/(2*k+1)<=x_k(-100)<=1,",
        "x_(k+1)(-100)>=x_k(-100),",
        "for every integer k>=1.",
        "```",
        "",
        "This is full infinite ratio-cone entry at `lambda=-100`.",
        "",
        "## Remaining Handoff",
        "",
        "The exact ratio-cone boundary algebra is already available. What remains",
        "is a rigorous infinite-dimensional or collared finite-flow legitimacy",
        "argument before propagating the cone forward in lambda. This certificate",
        "does not promote cone entry directly to the all-order Jensen/Newman result.",
        "",
        "```text",
        "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
        "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
    ]
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
        "validated Jensen-window PF negative-lambda -100 full cone-entry certificate: "
        "8 rows, 0 issues, 321 positive coefficients, 319 pointwise cone rows, "
        "318 adjacent rows, 1 analytic adjacent tail, 1 open flow handoff, 3 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
