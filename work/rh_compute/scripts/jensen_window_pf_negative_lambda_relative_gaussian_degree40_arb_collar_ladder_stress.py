#!/usr/bin/env python3
"""Build an Arb degree-ladder stress for the relative-Gaussian real-T collar."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import REPO_ROOT  # noqa: E402
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_diagnostics as build_arb_collar_diagnostics,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md"


@dataclass(frozen=True)
class Extremum:
    degree: int
    name: str
    sample: str


@dataclass(frozen=True)
class DegreeLadderRow:
    max_taylor_degree: int
    continuation_M: int
    ratio_ball_rows: int
    highest_ratio_degree: int
    highest_ratio_to_c0: str
    highest_ratio_sign: str
    positive_normalizer_rows: int
    certified_stencil_rows: int
    failed_bernstein_rows: int
    finite_degree_arb_collar_certified: bool
    normalizer_min_lowers: dict[str, str]
    stencil_min_lowers: dict[str, str]
    stencil_endpoints: dict[str, str]
    weakest_stencil_name: str
    weakest_stencil_min_lower: str
    proof_boundary: str


def arb_lower(value_text: str) -> flint.arb:
    return flint.arb(value_text).lower()


def weaker(left: tuple[flint.arb, Extremum] | None, degree: int, name: str, sample: str) -> tuple[flint.arb, Extremum]:
    candidate = (arb_lower(sample), Extremum(degree=degree, name=name, sample=sample))
    if left is None or candidate[0] < left[0]:
        return candidate
    return left


def parse_degree_values(text: str) -> list[int]:
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("at least one degree is required")
    if any(value < 0 or value % 2 for value in values):
        raise argparse.ArgumentTypeError("degrees must be non-negative even integers")
    if values != sorted(set(values)):
        raise argparse.ArgumentTypeError("degrees must be strictly increasing")
    return values


def build_ladder_rows(
    degrees: list[int],
    cutoff_n: int,
    precision_bits: int,
    k: int,
    collar_start_T: int,
) -> tuple[list[DegreeLadderRow], Extremum, Extremum]:
    rows: list[DegreeLadderRow] = []
    weakest_stencil: tuple[flint.arb, Extremum] | None = None
    weakest_normalizer: tuple[flint.arb, Extremum] | None = None

    for degree in degrees:
        diagnostics = build_arb_collar_diagnostics(
            degree,
            cutoff_n,
            precision_bits,
            k,
            degree // 2,
            collar_start_T,
        )
        ratio_rows = diagnostics["ratio_ball_rows"]
        normalizer_min_lowers = {
            row["name"]: row["min_bernstein_lower"] for row in diagnostics["normalizer_rows"]
        }
        stencil_min_lowers = {
            row["name"]: row["min_bernstein_lower"] for row in diagnostics["stencil_rows"]
        }
        stencil_endpoints = {
            row["name"]: row["endpoint_value"] for row in diagnostics["stencil_rows"]
        }
        row_weakest: tuple[flint.arb, Extremum] | None = None
        for name, sample in normalizer_min_lowers.items():
            weakest_normalizer = weaker(weakest_normalizer, degree, name, sample)
        for name, sample in stencil_min_lowers.items():
            weakest_stencil = weaker(weakest_stencil, degree, name, sample)
            row_weakest = weaker(row_weakest, degree, name, sample)
        if row_weakest is None:
            raise RuntimeError(f"no stencil rows for degree {degree}")

        highest_ratio = ratio_rows[-1]
        rows.append(
            DegreeLadderRow(
                max_taylor_degree=degree,
                continuation_M=degree // 2,
                ratio_ball_rows=len(ratio_rows),
                highest_ratio_degree=int(highest_ratio["degree"]),
                highest_ratio_to_c0=str(highest_ratio["ratio_to_c0"]),
                highest_ratio_sign=str(highest_ratio["sign"]),
                positive_normalizer_rows=int(diagnostics["positive_normalizer_rows"]),
                certified_stencil_rows=int(diagnostics["certified_stencil_rows"]),
                failed_bernstein_rows=int(diagnostics["failed_bernstein_rows"]),
                finite_degree_arb_collar_certified=bool(diagnostics["finite_degree_arb_collar_certified"]),
                normalizer_min_lowers=normalizer_min_lowers,
                stencil_min_lowers=stencil_min_lowers,
                stencil_endpoints=stencil_endpoints,
                weakest_stencil_name=row_weakest[1].name,
                weakest_stencil_min_lower=row_weakest[1].sample,
                proof_boundary="Finite degree-ladder stress only; not an infinite Taylor-tail theorem.",
            )
        )

    if weakest_stencil is None or weakest_normalizer is None:
        raise RuntimeError("degree ladder produced no extrema")
    return rows, weakest_stencil[1], weakest_normalizer[1]


def build_diagnostics(
    degrees: list[int],
    cutoff_n: int,
    precision_bits: int,
    k: int,
    collar_start_T: int,
) -> dict:
    rows, weakest_stencil, weakest_normalizer = build_ladder_rows(
        degrees,
        cutoff_n,
        precision_bits,
        k,
        collar_start_T,
    )
    return {
        "parameters": {
            "min_taylor_degree": min(degrees),
            "max_taylor_degree": max(degrees),
            "degree_step": min(b - a for a, b in zip(degrees, degrees[1:])) if len(degrees) > 1 else 0,
            "degree_values": degrees,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "collar_start_T": collar_start_T,
            "real_interval_u": f"[0, 1/{collar_start_T}]",
            "real_interval_T": f"[{collar_start_T}, infinity)",
        },
        "degree_ladder_rows": [asdict(row) for row in rows],
        "degree_ladder_row_count": len(rows),
        "all_degree_ladder_rows_certified": all(row.finite_degree_arb_collar_certified for row in rows),
        "total_failed_bernstein_rows": sum(row.failed_bernstein_rows for row in rows),
        "degree40_row": asdict(rows[-1]),
        "weakest_stencil_across_ladder": asdict(weakest_stencil),
        "weakest_normalizer_across_ladder": asdict(weakest_normalizer),
        "proof_boundary_note": (
            "This is an Arb finite-surrogate stress ladder through degree 40 on the same real-T collar. "
            "It does not estimate or bound the infinite residual Taylor tail beyond the last tested degree."
        ),
    }


def build_artifact(
    degrees: list[int],
    cutoff_n: int,
    precision_bits: int,
    k: int,
    collar_start_T: int,
) -> dict:
    diagnostics = build_diagnostics(degrees, cutoff_n, precision_bits, k, collar_start_T)
    rows = [
        {
            "id": "nlrgd40acl_01_degree_ladder_setup",
            "role": "finite_interval_ladder",
            "readiness": "not_ready_to_apply",
            "claim": "The Arb real-T collar test is repeated for every even finite surrogate degree from 16 through 40 at fixed k=22.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md",
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
            ],
            "proof_boundary": "Finite degree-ladder setup only; not an infinite Taylor-tail theorem.",
        },
        {
            "id": "nlrgd40acl_02_all_ladder_levels_bernstein_positive",
            "role": "finite_interval_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "For each tested degree 16,18,...,40, Bernstein coefficients certify all four normalizers and all three structured stencil polynomials on 0<=u<=1/1156.",
            "proof_boundary": "Finite finite-degree certificate ladder only; not an all-degree or all-k theorem.",
        },
        {
            "id": "nlrgd40acl_03_first_omitted_terms_survive",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The first omitted even Taylor levels after degree 16 do not introduce a Bernstein sign failure on the collar.",
            "proof_boundary": "Finite omitted-term stress only; not a bound for the remaining infinite tail.",
        },
        {
            "id": "nlrgd40acl_04_degree40_collar_endpoint_stability",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-40 endpoint stencil values remain positive on the T>=1156 collar, giving a stronger finite signal than the degree-16 surrogate alone.",
            "proof_boundary": "Finite endpoint-stability diagnostic only; not an analytic residual estimate.",
        },
        {
            "id": "nlrgd40acl_05_companion_stencil_bottleneck",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Across the tested ladder, the weakest stencil lower bound remains the companion stencil at the original degree-16 level.",
            "proof_boundary": "Finite bottleneck diagnostic only; not proof that all higher degrees or the infinite tail preserve the same bottleneck.",
        },
        {
            "id": "nlrgd40acl_06_live_residual_tail_upgrade",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "The next analytic upgrade is a residual-tail theorem beyond degree 40 for the B, companion, and weighted-gap stencils on the same collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md",
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md",
            ],
            "proof_boundary": "Live theorem-search route only; the residual-tail theorem is not proved.",
        },
        {
            "id": "nlrgd40acl_07_finite_ladder_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "Passing the finite Arb ladder through degree 40 proves the signed infinite-tail stencil theorem.",
            "gap": "The finite ladder stops at degree 40 and gives no analytic estimate for the residual series.",
            "proof_boundary": "Rejected finite-promotion template only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
        {
            "id": "nlrgd40acl_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Any promoted proof must state a residual-tail bound beyond a stated degree, the real-T collar, the fixed or uniform k range, and interval-safe comparisons against all three stencil margins.",
            "source_artifacts": [
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of the required residual-tail bound.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "degree_ladder_rows": diagnostics["degree_ladder_row_count"],
        "min_taylor_degree": diagnostics["parameters"]["min_taylor_degree"],
        "max_taylor_degree": diagnostics["parameters"]["max_taylor_degree"],
        "collar_start_T": str(collar_start_T),
        "degree40_positive_normalizer_rows": diagnostics["degree40_row"]["positive_normalizer_rows"],
        "degree40_certified_stencil_rows": diagnostics["degree40_row"]["certified_stencil_rows"],
        "total_failed_bernstein_rows": diagnostics["total_failed_bernstein_rows"],
        "all_degree_ladder_rows_certified": diagnostics["all_degree_ladder_rows_certified"],
        "weakest_stencil_across_ladder": diagnostics["weakest_stencil_across_ladder"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The same real-T collar T>=1156 survives Arb/Bernstein finite-surrogate stress for every even "
            "Taylor degree 16 through 40 at fixed k=22: all four normalizers and the B, companion, and "
            "weighted-gap stencils certify positive at every tested degree, with zero Bernstein failures. "
            "This makes the remaining gap sharper, but still finite: the proof still needs a residual-tail "
            "bound beyond the last tested degree and then removal of the fixed-k limitation."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress",
        "date": "2026-07-07",
        "status": "finite theorem-search diagnostic",
        "source_degree16_arb_collar_certificate": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md"
        ),
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "source_stencil_remainder_obligations": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md"
        ),
        "source_dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py"
        ),
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It stress-tests finite Arb coefficient-ball surrogates "
            "through degree 40 at fixed k=22 on T>=1156, but it does not bound the infinite residual Taylor tail, "
            "does not prove scaled-curvature monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The ladder is finite-degree through degree 40 only.",
            "The infinite residual Taylor tail beyond degree 40 remains unbounded.",
            "The result is fixed at k=22 and is not an all-k theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][0]["diagnostics"]
    degree40 = diagnostics["degree40_row"]
    weakest_stencil = diagnostics["weakest_stencil_across_ladder"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['degree_ladder_rows']} degree levels, "
        f"max degree {summary['max_taylor_degree']}, "
        f"{summary['total_failed_bernstein_rows']} failed Bernstein rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-40 Arb Collar Ladder Stress",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof",
        "of an infinite Taylor-tail theorem, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress`.",
        "",
        "Proof boundary: this artifact stress-tests finite Arb coefficient-ball",
        "surrogates through degree 40 on the same real-`T` collar. It does not",
        "bound the infinite residual Taylor tail.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Degree Ladder",
        "",
        "```text",
        f"real T interval: {diagnostics['parameters']['real_interval_T']}",
        f"u interval: {diagnostics['parameters']['real_interval_u']}",
        f"degree values: {diagnostics['parameters']['degree_values']}",
        f"degree levels: {summary['degree_ladder_rows']}",
        f"all degree levels certified: {summary['all_degree_ladder_rows_certified']}",
        f"total failed Bernstein rows: {summary['total_failed_bernstein_rows']}",
        f"weakest stencil lower: degree {weakest_stencil['degree']} {weakest_stencil['name']} = {weakest_stencil['sample']}",
        "```",
        "",
        "Per-degree stress rows:",
        "",
        "```text",
    ]
    for row in diagnostics["degree_ladder_rows"]:
        lines.append(
            f"degree={row['max_taylor_degree']}, M={row['continuation_M']}: "
            f"normalizers={row['positive_normalizer_rows']}/4, "
            f"stencils={row['certified_stencil_rows']}/3, "
            f"failed Bernstein={row['failed_bernstein_rows']}, "
            f"weakest stencil={row['weakest_stencil_name']} {row['weakest_stencil_min_lower']}"
        )
    lines.extend(
        [
            "```",
            "",
            "## Degree-40 Endpoint Snapshot",
            "",
            "```text",
            f"highest ratio c{degree40['highest_ratio_degree']}/c0: {degree40['highest_ratio_to_c0']} ({degree40['highest_ratio_sign']})",
            f"B_product endpoint: {degree40['stencil_endpoints']['B_product']}",
            f"companion_product endpoint: {degree40['stencil_endpoints']['companion_product']}",
            f"weighted_gap_derivative endpoint: {degree40['stencil_endpoints']['weighted_gap_derivative']}",
            "```",
            "",
            "Interpretation:",
            "",
            "The degree-16 Arb collar was not an isolated finite-truncation fluke",
            "through the next twelve even truncation levels. The remaining task is",
            "now cleaner but still open: replace finite degree-by-degree stress with",
            "a signed residual-tail bound beyond degree 40 on the same collar, then",
            "remove the fixed-`k` limitation.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree16_arb_collar_certificate"],
            artifact["source_uniform_remainder_target"],
            artifact["source_stencil_remainder_obligations"],
            artifact["source_dependency_graph"],
            "```",
            "",
            "Summary:",
            "",
            summary["main_finding"],
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degrees", type=parse_degree_values, default=parse_degree_values("16,18,20,22,24,26,28,30,32,34,36,38,40"))
    parser.add_argument("--cutoff-n", type=int, default=80)
    parser.add_argument("--precision-bits", type=int, default=384)
    parser.add_argument("--tail-start-k", type=int, default=22)
    parser.add_argument("--collar-start-T", type=int, default=1156)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(args.degrees, args.cutoff_n, args.precision_bits, args.tail_start_k, args.collar_start_T)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
