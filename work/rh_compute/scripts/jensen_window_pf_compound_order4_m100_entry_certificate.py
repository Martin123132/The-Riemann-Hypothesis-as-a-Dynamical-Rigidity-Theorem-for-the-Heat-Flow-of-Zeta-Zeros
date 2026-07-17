#!/usr/bin/env python3
"""Compose the complete contiguous order-four entry theorem at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md"
)
SOURCE_CONDENSATION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_condensation_gate.json"
)
SOURCE_COMPACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json"
)
SOURCE_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.json"
)
SOURCE_BRIDGE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_first_summand_curvature_bridge.json"
)

PREFIX_LAST_N = 316
TAIL_FIRST_N = 317
TAIL_FIRST_K = 320


@dataclass(frozen=True)
class EntryRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def source_row(artifact: dict, row_id: str) -> dict:
    matches = [row for row in artifact.get("rows", []) if row.get("id") == row_id]
    if len(matches) != 1:
        raise RuntimeError(f"expected one source row {row_id}, found {len(matches)}")
    return matches[0]


def exact_tail_arithmetic() -> dict:
    k0 = TAIL_FIRST_K
    first_summand_margin = Fraction(k0**2 - 36, 10 * k0**2 * (k0**2 - 1))
    if first_summand_margin <= 0:
        raise RuntimeError("first-summand rational margin failed")

    perturbation_coefficients = (
        1,
        1902,
        1_482_055,
        604_691_300,
        135_889_991_935,
        15_877_422_036_942,
        748_002_501_678_169,
    )
    if any(coefficient <= 0 for coefficient in perturbation_coefficients):
        raise RuntimeError("full-kernel perturbation polynomial failed")

    buffer_numerator_at_k0 = 753 * k0**2 - 1000 * (2 * k0 + 1)
    if buffer_numerator_at_k0 != 76_466_200:
        raise RuntimeError("tail buffer endpoint numerator changed")
    shifted_buffer_coefficients = (
        753,
        2 * 753 * k0 - 2000,
        buffer_numerator_at_k0,
    )
    if any(coefficient <= 0 for coefficient in shifted_buffer_coefficients):
        raise RuntimeError("tail buffer shifted polynomial failed")

    first_cap = Fraction(18, 5)
    perturbation_cap = Fraction(2, 5)
    total_cap = first_cap + perturbation_cap
    if total_cap != 4:
        raise RuntimeError("tail budget sum failed")
    return {
        "indices": "n>=317, k=n+3>=320, |s|<=1 => k+s>=319",
        "tent_identity": (
            "integral_[-1,1](1-|s|)/(k+s)^2 ds"
            "=-log(1-1/k^2)"
        ),
        "log_inequality": "-log(1-z)<=z/(1-z), 0<=z<1",
        "first_summand_bound": (
            "P_n^(1)<=7/(2*(k^2-1))<=18/(5*k^2)"
        ),
        "first_summand_margin_at_k320": fraction_text(first_summand_margin),
        "first_summand_margin_formula": (
            "18/(5*k^2)-7/(2*(k^2-1))"
            "=(k^2-36)/(10*k^2*(k^2-1))"
        ),
        "perturbation_bound": (
            "|P_n-P_n^(1)|<=10112*(k+1)^2/(k-3)^6<=2/(5*k^2)"
        ),
        "perturbation_shifted_coefficients_descending": list(
            perturbation_coefficients
        ),
        "total_penalty_bound": "P_n<=4/k^2",
        "defect_buffer": (
            "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))"
        ),
        "buffer_comparison": "4/k^2<753/(250*(2*k+1))",
        "buffer_shifted_coefficients_descending": list(
            shifted_buffer_coefficients
        ),
        "buffer_endpoint_numerator": buffer_numerator_at_k0,
        "strict_tail_consequence": "P_n<-3*log(x_k), n>=317",
    }


def build_artifact() -> dict:
    condensation = json.loads(SOURCE_CONDENSATION.read_text(encoding="utf-8"))
    compact = json.loads(SOURCE_COMPACT.read_text(encoding="utf-8"))
    ray = json.loads(SOURCE_RAY.read_text(encoding="utf-8"))
    bridge = json.loads(SOURCE_BRIDGE.read_text(encoding="utf-8"))

    condensation_summary = condensation.get("summary", {})
    if condensation_summary.get("prefix_order_four_margins") != 317:
        raise RuntimeError("order-four prefix source is incomplete")
    if condensation_summary.get("prefix_scaled_penalty_caps") != 317:
        raise RuntimeError("prefix penalty source is incomplete")
    compact_summary = compact.get("summary", {})
    if compact_summary.get("all_blocks_passed") is not True:
        raise RuntimeError("compact curvature source is not closed")
    ray_summary = ray.get("summary", {})
    if ray_summary.get("global_corridor_to_curvature_closed") is not True:
        raise RuntimeError("global corridor-to-curvature source is not closed")
    perturbation_row = source_row(bridge, "co4fcb_07_full_kernel_transfer")
    if perturbation_row.get("readiness") != "ready_to_apply":
        raise RuntimeError("full-kernel perturbation source is not ready")
    localization_row = source_row(bridge, "co4fcb_04b_localized_curvature")
    if localization_row.get("readiness") != "ready_to_apply":
        raise RuntimeError("localized curvature reduction is not ready")

    arithmetic = exact_tail_arithmetic()
    prefix = condensation.get("finite", {})
    rows = [
        EntryRow(
            id="co4m100ec_01_condensation_coordinate",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Desnanot-Jacobi condensation identifies the positive contiguous order-four sign with strict log-concavity of the positive order-three magnitudes.",
            formula="H_(4,n)>0 iff T_(n+1)^2>T_n*T_(n+2)",
            proof_boundary="Contiguous order-four sign coordinate only.",
        ),
        EntryRow(
            id="co4m100ec_02_gap_penalty_coordinate",
            role="exact_sign_equivalence",
            readiness="ready_to_apply",
            claim="The defect-gap factorization converts the order-four sign into one logarithmic penalty inequality.",
            formula="H_(4,n)>0 iff P_n<-3*log(x_(n+3))",
            proof_boundary="Exact lambda-independent coordinate reduction only.",
        ),
        EntryRow(
            id="co4m100ec_03_repaired_prefix",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="The repaired Arb coefficient source proves every prefix order-four margin strictly positive at lambda=-100.",
            formula="H_(4,n)>0, 0<=n<=316, lambda=-100",
            proof_boundary="Finite prefix only.",
            diagnostics={
                "prefix_rows": 317,
                "n_range": [0, PREFIX_LAST_N],
                "minimum_margin_lower": prefix.get("minimum_margin_lower"),
                "minimum_margin_n": prefix.get("minimum_margin_n"),
            },
        ),
        EntryRow(
            id="co4m100ec_04_global_first_curvature",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The compact interval, exact-corridor finite cover, analytic ray, and stable localization prove the first-summand curvature ceiling on the complete real tail.",
            formula="K_1(t)<=7/(2*t^2), every real t>=319",
            proof_boundary="First Newman summand at lambda=-100 only.",
            diagnostics={
                "compact_range": "319<=t<=V'(2)",
                "mode_tail": "u>=2",
                "compact_blocks": compact_summary.get("accepted_central_blocks"),
                "ray_scaled_U_upper": ray.get("scalar_comparison", {}).get(
                    "scaled_U_upper"
                ),
            },
        ),
        EntryRow(
            id="co4m100ec_05_tent_transfer",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The exact tent identity transfers the continuous curvature ceiling to the first-summand discrete penalty.",
            formula="P_n^(1)<=18/(5*k^2), k=n+3>=320",
            proof_boundary="First-summand discrete penalty only.",
            diagnostics=arithmetic,
        ),
        EntryRow(
            id="co4m100ec_06_full_kernel_transfer",
            role="exact_perturbation_theorem",
            readiness="ready_to_apply",
            claim="The proved all-summand dominance estimate changes the logarithmic order-four penalty by at most two fifths inverse square.",
            formula="|P_n-P_n^(1)|<=2/(5*k^2), k>=320",
            proof_boundary="Full Newman kernel at lambda=-100 on the analytic tail.",
            diagnostics={
                "source_formula": perturbation_row.get("formula"),
                "shifted_polynomial": arithmetic[
                    "perturbation_shifted_coefficients_descending"
                ],
            },
        ),
        EntryRow(
            id="co4m100ec_07_complete_penalty_tail",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The first-summand and full-kernel budgets prove the sufficient complete-kernel penalty ceiling.",
            formula="P_n<=4/(n+3)^2, n>=317",
            proof_boundary="Complete-kernel logarithmic penalty only.",
        ),
        EntryRow(
            id="co4m100ec_08_strict_buffer",
            role="exact_rational_theorem",
            readiness="ready_to_apply",
            claim="The scaled-defect anchor leaves a strict buffer above the complete penalty ceiling at every tail index.",
            formula="4/k^2<753/(250*(2*k+1))<=-3*log(x_k), k>=320",
            proof_boundary="Tail conversion from the penalty ceiling to the order-four sign.",
            diagnostics={
                "endpoint_numerator": arithmetic["buffer_endpoint_numerator"],
                "shifted_coefficients_descending": arithmetic[
                    "buffer_shifted_coefficients_descending"
                ],
            },
        ),
        EntryRow(
            id="co4m100ec_09_all_shift_entry",
            role="exact_entry_theorem",
            readiness="ready_to_apply",
            claim="The 317-row interval prefix and strict analytic tail prove every contiguous order-four Hankel minor has the required positive sign at lambda=-100.",
            formula="H_(4,n)(-100)>0 for every integer n>=0",
            proof_boundary="Contiguous order-four entry at one heat parameter only.",
        ),
        EntryRow(
            id="co4m100ec_10_forward_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Derive and close a legitimate infinite-index order-four compound-flow propagation theorem from lambda=-100 through lambda=0.",
            formula="H_(4,n)(-100)>0 for all n => H_(4,n)(lambda)>0 for -100<=lambda<=0",
            proof_boundary="Open; no order-four propagation, arbitrary-column order four, PF-infinity, RH, or Lambda<=0 is asserted.",
        ),
    ]
    source_paths = (
        SOURCE_CONDENSATION,
        SOURCE_COMPACT,
        SOURCE_RAY,
        SOURCE_BRIDGE,
    )
    return {
        "kind": "jensen_window_pf_compound_order4_m100_entry_certificate",
        "date": "2026-07-13",
        "status": "all-shift contiguous order-four entry theorem at lambda=-100",
        "proof_boundary": (
            "This artifact proves H_(4,n)(-100)>0 for every n>=0. It does not "
            "prove forward order-four invariance, arbitrary-column order four, "
            "PF-infinity, the Riemann hypothesis, or Lambda<=0."
        ),
        "exact": {
            "condensation_sign": "H_(4,n)>0 iff T_(n+1)^2>T_n*T_(n+2)",
            "gap_penalty_sign": "H_(4,n)>0 iff P_n<-3*log(x_(n+3))",
            "global_first_summand_curvature": "K_1(t)<=7/(2*t^2), t>=319",
            "complete_tail_penalty": "P_n<=4/(n+3)^2, n>=317",
            "all_shift_entry": "H_(4,n)(-100)>0, every n>=0",
        },
        "tail_arithmetic": arithmetic,
        "prefix": {
            "n_range": [0, PREFIX_LAST_N],
            "rows": 317,
            "tail_n_range": f"n>={TAIL_FIRST_N}",
            "tail_k_range": f"k>={TAIL_FIRST_K}",
            "minimum_margin_lower": prefix.get("minimum_margin_lower"),
            "minimum_margin_n": prefix.get("minimum_margin_n"),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": len(rows) - 1,
            "open_analytic_rows": 1,
            "prefix_order_four_margins": 317,
            "analytic_order_four_tails": 1,
            "global_curvature_theorems": 1,
            "full_kernel_transfers": 1,
            "all_shift_order_four_entry_theorems": 1,
            "open_forward_handoffs": 1,
        },
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in source_paths
        ],
        "source_hashes": {
            str(path.relative_to(REPO_ROOT)).replace("\\", "/"): sha256(path)
            for path in source_paths
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_m100_entry_certificate.py"
        ),
        "remaining_target": (
            "Derive the order-four compound heat-flow system, establish a "
            "weighted infinite maximum principle, and propagate the strict "
            "lambda=-100 entry through lambda=0."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    arithmetic = artifact["tail_arithmetic"]
    prefix = artifact["prefix"]
    lines = [
        "# Jensen-Window PF Compound Order-Four Lambda=-100 Entry Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: all-shift contiguous order-four entry theorem at `lambda=-100`.",
        "This proves one complete compound layer at one heat parameter. It is not",
        "a proof of forward order-four invariance, arbitrary-column order four,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_m100_entry_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_m100_entry_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_m100_entry_certificate.py",
        "```",
        "",
        "## Sign Coordinate",
        "",
        "Desnanot-Jacobi condensation and the completed order-two and order-three",
        "signs give",
        "",
        "```text",
        "H_(4,n)>0 iff T_(n+1)^2>T_n*T_(n+2).",
        "```",
        "",
        "The exact order-three gap factorization then gives the equivalent stable",
        "penalty condition",
        "",
        "```text",
        "P_n=log(G_n*G_(n+2)/G_(n+1)^2)",
        "H_(4,n)>0 iff P_n<-3*log(x_(n+3)).",
        "```",
        "",
        "## Repaired Prefix",
        "",
        "The 1024-bit repaired Arb source proves",
        "",
        "```text",
        "H_(4,n)(-100)>0, 0<=n<=316.",
        "```",
        "",
        f"All `{prefix['rows']}` margins are strict. The minimum lower enclosure is",
        f"`{prefix['minimum_margin_lower']}` at `n={prefix['minimum_margin_n']}`.",
        "",
        "## Global Curvature",
        "",
        "The compact interval certificate covers `319<=t<=V'(2)`. The complete",
        "finite exact-corridor cover and the analytic `u>=20` ray prove the",
        "localized ceiling on every mode `u>=2`; the stable localization converts",
        "that ceiling to the first-summand curvature. Therefore",
        "",
        "```text",
        "K_1(t)<=7/(2*t^2) for every real t>=319.",
        "```",
        "",
        "## Analytic Tail",
        "",
        "For `n>=317`, put `k=n+3>=320`. Then `k+s>=319` throughout the",
        "tent interval, and",
        "",
        "```text",
        "P_n^(1)=integral_[-1,1](1-|s|)*K_1(k+s) ds",
        "       <=(7/2)*[-log(1-1/k^2)]",
        "       <=7/(2*(k^2-1))",
        "       <=18/(5*k^2).",
        "```",
        "",
        "The already-proved complete-kernel perturbation gives",
        "",
        "```text",
        "|P_n-P_n^(1)|<=2/(5*k^2),",
        "P_n<=4/k^2.",
        "```",
        "",
        "The scaled-defect anchor supplies the strict remaining buffer:",
        "",
        "```text",
        "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))>4/k^2.",
        "```",
        "",
        "After `k=320+m`, the numerator of the last strict comparison has",
        "positive coefficients",
        "",
        "```text",
        f"{arithmetic['buffer_shifted_coefficients_descending']}",
        "```",
        "",
        f"and endpoint value `{arithmetic['buffer_endpoint_numerator']}`. Hence",
        "`P_n<-3*log(x_k)` and therefore `H_(4,n)>0` for every `n>=317`.",
        "",
        "## Entry Theorem",
        "",
        "Combining the finite prefix and analytic tail proves",
        "",
        "```text",
        "H_(4,n)(-100)>0 for every integer n>=0.",
        "```",
        "",
        "This closes the contiguous order-four entry problem at `lambda=-100`.",
        "The next live theorem is forward propagation: derive the order-four",
        "compound flow and prevent a first loss of positivity from escaping to",
        "infinite index before carrying the layer through `lambda=0`.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_condensation_gate.md",
        "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "certified Jensen-window PF compound order-four lambda=-100 entry: "
        f"{summary['rows']} rows, {summary['exact_rows']} exact rows, "
        f"{summary['prefix_order_four_margins']} prefix margins, "
        f"{summary['analytic_order_four_tails']} analytic tail, "
        f"{summary['all_shift_order_four_entry_theorems']} all-shift entry theorem, "
        f"{summary['open_forward_handoffs']} open forward handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
