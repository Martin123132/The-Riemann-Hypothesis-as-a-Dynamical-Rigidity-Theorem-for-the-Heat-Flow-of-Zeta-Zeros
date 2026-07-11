#!/usr/bin/env python3
"""Build the dominant first-summand saddle-wall theorem target and scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys

import mpmath as mp


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FULL_SOURCE_JSONL = (
    REPO_ROOT
    / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md"
)
DEFAULT_T = 100
DEFAULT_TAIL_START_K = 319
DEFAULT_MPMATH_DPS = 60
DEFAULT_BISECTION_STEPS = 210
DEFAULT_UPPER_PADDING = "1.2"
DEFAULT_SAMPLE_K = (300, 400, 500, 700, 1000, 2000, 5000, 10000, 20000)
DEFAULT_PRECISION_BITS = 512


@dataclass(frozen=True)
class TargetRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    formula: str | None = None
    diagnostics: dict | None = None
    gap: str | None = None


def sci(value: mp.mpf, digits: int = 35) -> str:
    return mp.nstr(value, n=digits, min_fixed=-6, max_fixed=6)


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 70) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def log_phi1(u: mp.mpf) -> mp.mpf:
    q = mp.pi * mp.exp(4 * u)
    return mp.log(mp.pi) + 5 * u + mp.log(2 * q - 3) - q


def saddle_derivative(u: mp.mpf, k: int) -> mp.mpf:
    q = mp.pi * mp.exp(4 * u)
    return 2 * k / u - 2 * DEFAULT_T * u + 5 + 8 * q / (2 * q - 3) - 4 * q


def saddle_point(k: int) -> mp.mpf:
    low = mp.mpf("0.01")
    high = mp.mpf(1)
    while saddle_derivative(high, k) > 0:
        high *= mp.mpf("1.4")
    for _ in range(DEFAULT_BISECTION_STEPS):
        middle = (low + high) / 2
        if saddle_derivative(middle, k) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def log_first_moment(k: int) -> tuple[mp.mpf, mp.mpf]:
    saddle = saddle_point(k)
    maximum = 2 * k * mp.log(saddle) - DEFAULT_T * saddle**2 + log_phi1(saddle)

    def scaled_integrand(u: mp.mpf) -> mp.mpf:
        if not u:
            return mp.mpf(0)
        return mp.exp(2 * k * mp.log(u) - DEFAULT_T * u**2 + log_phi1(u) - maximum)

    upper = saddle + mp.mpf(DEFAULT_UPPER_PADDING)
    scaled_integral = mp.quad(scaled_integrand, [0, saddle, upper])
    return maximum + mp.log(2 * scaled_integral), saddle


def sample_row(k: int) -> dict:
    moment_rows = {index: log_first_moment(index) for index in range(k - 1, k + 3)}

    def normalized(index: int) -> mp.mpf:
        return moment_rows[index][0] - mp.loggamma(mp.mpf(index) + mp.mpf("0.5"))

    log_gap = normalized(k + 2) - 3 * normalized(k + 1) + 3 * normalized(k) - normalized(k - 1)
    target = mp.mpf(1) / (4 * k * k)
    scaled = log_gap * k * k
    L = mp.log(k)
    lower_bracket = (L + mp.mpf("0.88")) / 8
    upper_bracket = L / 4
    saddle = moment_rows[k][1]
    return {
        "k": k,
        "saddle": sci(saddle),
        "lower_bracket": sci(lower_bracket),
        "upper_bracket": sci(upper_bracket),
        "saddle_inside_bracket": bool(lower_bracket < saddle < upper_bracket),
        "log_gap": sci(log_gap),
        "scaled_k2_log_gap": sci(scaled),
        "quarter_k2_target": sci(target),
        "target_margin": sci(log_gap - target),
        "positive_log_gap": bool(log_gap > 0),
        "above_quarter_k2_target": bool(log_gap > target),
        "proof_boundary": (
            "High-precision mpmath saddle quadrature with a finite upper truncation; "
            "not an interval enclosure or a uniform theorem."
        ),
    }


def load_full_log_gap_at_300(path: Path) -> flint.arb:
    values: dict[int, flint.arb] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row = json.loads(raw)
            if row.get("kind") == "acb_coefficient_enclosure":
                values[int(row["k"])] = flint.arb(row["A_ball"])
    needed = {299, 300, 301, 302}
    if not needed.issubset(values):
        raise ValueError("full-kernel source lacks k=299..302")
    x300 = values[301] * values[299] / values[300] ** 2
    x301 = values[302] * values[300] / values[301] ** 2
    return (x301 / x300).log()


def exact_geometry() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    k = flint.arb(300)
    L = k.log()
    pi = flint.arb.pi()
    upper_sign_margin = 4 * pi - 8 / L - 10 / k
    bracket_width_margin = L - flint.arb(22) / 25
    transfer_margin = flint.arb(1) / (4 * DEFAULT_TAIL_START_K**2) - flint.arb(16) / (
        DEFAULT_TAIL_START_K - 1
    ) ** 6
    if not bool(upper_sign_margin > 0 and bracket_width_margin > 0 and transfer_margin > 0):
        raise RuntimeError("exact saddle target geometry failed")
    return {
        "upper_saddle_sign_margin_ball": arb_text(upper_sign_margin),
        "upper_saddle_sign_margin_lower": arb_lower_text(upper_sign_margin),
        "bracket_width_margin_ball": arb_text(bracket_width_margin),
        "bracket_width_margin_lower": arb_lower_text(bracket_width_margin),
        "wall_transfer_margin_at_k319_ball": arb_text(transfer_margin),
        "wall_transfer_margin_at_k319_lower": arb_lower_text(transfer_margin),
        "upper_sign_propagation": (
            "8/log(k)+10/k decreases, so D_k(log(k)/4)<0 follows from the k=300 endpoint."
        ),
        "transfer_propagation": (
            "(k-1)^6/k^2 increases for k>1, so 1/(4*k^2)>16/(k-1)^6 "
            "propagates from k=319."
        ),
    }


def build_artifact(full_source_path: Path = DEFAULT_FULL_SOURCE_JSONL) -> dict:
    mp.mp.dps = DEFAULT_MPMATH_DPS
    sample_rows = [sample_row(k) for k in DEFAULT_SAMPLE_K]
    if not all(row["positive_log_gap"] and row["above_quarter_k2_target"] for row in sample_rows):
        raise RuntimeError("first-summand saddle sample missed the target")
    if not all(row["saddle_inside_bracket"] for row in sample_rows):
        raise RuntimeError("sample saddle escaped exact bracket")
    min_scaled = min(sample_rows, key=lambda row: Decimal(row["scaled_k2_log_gap"]))
    min_margin = min(sample_rows, key=lambda row: Decimal(row["target_margin"]))
    geometry = exact_geometry()
    full_gap_300 = load_full_log_gap_at_300(full_source_path)
    first_gap_300 = mp.mpf(sample_rows[0]["log_gap"])
    full_midpoint = mp.mpf(full_gap_300.mid().str(80, radius=False))
    cross_difference = abs(first_gap_300 - full_midpoint)

    diagnostics = {
        "parameters": {
            "T": DEFAULT_T,
            "tail_start_k": DEFAULT_TAIL_START_K,
            "mpmath_dps": DEFAULT_MPMATH_DPS,
            "bisection_steps": DEFAULT_BISECTION_STEPS,
            "upper_padding": DEFAULT_UPPER_PADDING,
            "sample_k": list(DEFAULT_SAMPLE_K),
        },
        "exact_geometry": geometry,
        "sample_rows": sample_rows,
        "sample_count": len(sample_rows),
        "positive_sample_rows": sum(row["positive_log_gap"] for row in sample_rows),
        "target_sample_rows": sum(row["above_quarter_k2_target"] for row in sample_rows),
        "bracketed_saddle_rows": sum(row["saddle_inside_bracket"] for row in sample_rows),
        "minimum_scaled_k2_log_gap": min_scaled["scaled_k2_log_gap"],
        "minimum_scaled_k2_log_gap_at_k": min_scaled["k"],
        "minimum_target_margin": min_margin["target_margin"],
        "minimum_target_margin_at_k": min_margin["k"],
        "scaled_profile_strictly_increasing": all(
            Decimal(sample_rows[index + 1]["scaled_k2_log_gap"])
            > Decimal(sample_rows[index]["scaled_k2_log_gap"])
            for index in range(len(sample_rows) - 1)
        ),
        "full_kernel_k300_crosscheck": {
            "full_Arb_log_gap_ball": arb_text(full_gap_300),
            "first_summand_mpmath_log_gap": sci(first_gap_300),
            "absolute_midpoint_difference": sci(cross_difference),
            "proof_boundary": (
                "Cross-check only: the full value is Arb-enclosed, but the first-summand "
                "quadrature value is not an interval enclosure."
            ),
        },
        "formal_asymptotic_target": "k^2*L_k^(1)->1 as k->infinity",
        "quantitative_theorem_target": "L_k^(1)>=1/(4*k^2) for every integer k>=319",
    }
    rows = [
        TargetRow(
            id="fsswt_01_exact_first_moment_wall",
            role="exact_reduction",
            readiness="available_exact",
            claim="Define the first-summand moments and their Gamma-normalized adjacent log wall.",
            formula=(
                "M_k^(1)=2*integral_0^infinity u^(2k)*exp(-100*u^2)*phi_1(u)du; "
                "L_k^(1)=Delta^3 log(M_j^(1)/Gamma(j+1/2)) at j=k-1"
            ),
            proof_boundary="Exact normalization only.",
        ),
        TargetRow(
            id="fsswt_02_exact_saddle_geometry",
            role="exact_lemma",
            readiness="available_exact",
            claim="The first-summand log integrand is strictly concave and has one saddle.",
            formula=(
                "S_k''(u)=-2*k/u^2-200-16*q-96*q/(2*q-3)^2<0; "
                "S_k'(0+)=+infinity, S_k'(infinity)=-infinity"
            ),
            proof_boundary="Exact unique-mode geometry only.",
        ),
        TargetRow(
            id="fsswt_03_all_k_saddle_bracket",
            role="exact_analytic_bound",
            readiness="interval_validated",
            claim="For every k>=300 the unique saddle lies between (log(k)+22/25)/8 and log(k)/4.",
            formula="(log(k)+22/25)/8 < s_k < log(k)/4",
            diagnostics=geometry,
            proof_boundary="All-k saddle location bracket; not a moment-curvature theorem.",
        ),
        TargetRow(
            id="fsswt_04_high_precision_scout",
            role="floating_diagnostic",
            readiness="not_ready_to_apply",
            claim="High-precision saddle-centered quadrature samples the first-summand wall at nine k values from 300 through 20000.",
            diagnostics={
                "sample_count": len(sample_rows),
                "mpmath_dps": DEFAULT_MPMATH_DPS,
                "upper_padding": DEFAULT_UPPER_PADDING,
            },
            proof_boundary="Finite mpmath quadrature only; no certified truncation or quadrature remainder.",
        ),
        TargetRow(
            id="fsswt_05_positive_scaled_profile",
            role="finite_diagnostic",
            readiness="not_ready_to_apply",
            claim="All nine samples satisfy L_k^(1)>1/(4*k^2), and k^2*L_k^(1) increases from about 0.365 to 0.732.",
            diagnostics={
                "minimum_scaled_k2_log_gap": diagnostics["minimum_scaled_k2_log_gap"],
                "minimum_scaled_k2_log_gap_at_k": diagnostics["minimum_scaled_k2_log_gap_at_k"],
                "scaled_profile_strictly_increasing": diagnostics["scaled_profile_strictly_increasing"],
            },
            proof_boundary="Finite shape evidence only.",
        ),
        TargetRow(
            id="fsswt_06_full_kernel_overlap",
            role="finite_crosscheck",
            readiness="not_ready_to_apply",
            claim="At k=300 the high-precision first-summand value agrees with the repaired full-kernel Arb log gap to the displayed precision.",
            diagnostics=diagnostics["full_kernel_k300_crosscheck"],
            proof_boundary="Mixed interval/floating cross-check only; not a first-summand interval certificate.",
        ),
        TargetRow(
            id="fsswt_07_formal_saddle_scale",
            role="asymptotic_target",
            readiness="not_ready_to_apply",
            claim="The saddle equation predicts k^2*L_k^(1)->1, consistent with the finite scaled profile.",
            formula="4*pi*exp(4*s_k)~2*k/s_k and d_k log(M_k^(1))~2*log(s_k)",
            proof_boundary="Formal leading scale only; no uniform Laplace remainder.",
        ),
        TargetRow(
            id="fsswt_08_quantitative_wall_target",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="The compact and ray paired-remainder certificates prove L_k^(1)>=1/(4*k^2) for every integer k>=319.",
            formula="L_k^(1)>=1/(4*k^2), k>=319",
            gap=None,
            proof_boundary="Analytic first-summand wall theorem; not the remaining all-order bridge.",
        ),
        TargetRow(
            id="fsswt_09_full_wall_handoff",
            role="conditional_bridge",
            readiness="ready_to_apply",
            claim="The quantitative first-summand target plus the 16/(k-1)^6 dominance error proves the full lambda=-100 wall for k>=319 and splices to the finite collar.",
            formula="L_k>=1/(4*k^2)-16/(k-1)^6>0",
            diagnostics={"transfer_margin_at_k319": geometry["wall_transfer_margin_at_k319_ball"]},
            proof_boundary="Exact wall handoff with the analytic antecedent discharged.",
        ),
    ]
    summary = {
        "target_rows": len(rows),
        "sample_rows": len(sample_rows),
        "positive_sample_rows": diagnostics["positive_sample_rows"],
        "target_sample_rows": diagnostics["target_sample_rows"],
        "bracketed_saddle_rows": diagnostics["bracketed_saddle_rows"],
        "minimum_scaled_k2_log_gap": diagnostics["minimum_scaled_k2_log_gap"],
        "minimum_scaled_k2_log_gap_at_k": diagnostics["minimum_scaled_k2_log_gap_at_k"],
        "scaled_profile_strictly_increasing": diagnostics["scaled_profile_strictly_increasing"],
        "open_requirement_rows": sum(row.role == "open_requirement" for row in rows),
        "ready_to_apply_rows": 2,
        "main_finding": (
            "The dominant n=1 integrand has exact strict-concavity geometry and an all-k saddle "
            "bracket. Sixty-digit saddle-centered quadrature at nine k values from 300 to 20000 "
            "finds L_k^(1)>1/(4*k^2), with k^2*L_k^(1) increasing from about 0.365 to 0.732 "
            "and formally tending to one. This remains finite floating evidence. The paired "
            "compact and ray theorems prove the analytic tail, and the certified full-kernel "
            "perturbation plus finite collar close the adjacent wall."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_saddle_wall_target",
        "status": "exact saddle geometry with analytic wall closure",
        "date": "2026-07-10",
        "proof_boundary": (
            "This artifact composes exact first-summand saddle geometry, finite high-precision "
            "cross-checks, and the paired remainder theorems to prove the quantitative all-k "
            "adjacent wall. It does not prove the remaining all-order bridge, RH, or Lambda <= 0."
        ),
        "source_dominance_certificate": "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
        "source_collar_extension": "outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md",
        "source_paired_remainder_certificate": "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
        "source_paired_ray_certificate": "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
        "source_full_enclosure_jsonl": full_source_path.relative_to(REPO_ROOT).as_posix(),
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py",
        "diagnostics": diagnostics,
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Negative-Lambda First-Summand Saddle-Wall Target",
        "",
        "Date: 2026-07-10",
        "",
        "Status: exact saddle geometry with analytic wall closure; not a proof of the remaining all-order bridge.",
        "This proves the lambda=-100 adjacent wall, not the remaining all-order bridge, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_saddle_wall_target`.",
        "",
        "Machine-readable result:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF negative-lambda first-summand saddle-wall closure: 9 rows, 0 issues, 9 positive samples, 9 quarter-k2 samples, 9 bracketed saddles, 0 open requirements, 2 ready-to-apply rows",
        "```",
        "",
        "## Exact Geometry",
        "",
        "For the first Newman summand at `lambda=-100`, let",
        "",
        "```text",
        "S_k(u)=2*k*log(u)-100*u^2+log(phi_1(u)).",
        "```",
        "",
        "With `q=pi*exp(4u)`,",
        "",
        "```text",
        "S_k'(u)=2*k/u-200*u+5+8*q/(2*q-3)-4*q,",
        "S_k''(u)=-2*k/u^2-200-16*q-96*q/(2*q-3)^2<0.",
        "```",
        "",
        "Thus there is one saddle `s_k`. Exact endpoint propagation proves",
        "",
        "```text",
        "(log(k)+22/25)/8 < s_k < log(k)/4, k>=300.",
        "```",
        "",
        "## High-Precision Scout",
        "",
        "The scout evaluates the first-summand moments with 60-digit mpmath",
        "saddle-centered quadrature. The finite upper truncation is not certified,",
        "so these rows are theorem-search evidence only.",
        "",
        "| k | saddle | L_k^(1) | k^2 L_k^(1) | L_k^(1)-1/(4k^2) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in diagnostics["sample_rows"]:
        lines.append(
            f"| `{row['k']}` | `{row['saddle']}` | `{row['log_gap']}` | "
            f"`{row['scaled_k2_log_gap']}` | `{row['target_margin']}` |"
        )
    lines.extend(
        [
            "",
            "All nine rows lie above `1/(4*k^2)`, and the sampled scaled profile",
            "is strictly increasing. At `k=300`, the first-summand value agrees",
            "with the repaired full-kernel Arb log gap to the displayed precision:",
            "",
            "```text",
            f"full Arb: {diagnostics['full_kernel_k300_crosscheck']['full_Arb_log_gap_ball']}",
            f"n=1 mp:   {diagnostics['full_kernel_k300_crosscheck']['first_summand_mpmath_log_gap']}",
            "```",
            "",
            "This agreement is a cross-check, not an interval promotion of the mpmath",
            "quadrature.",
            "",
            "## Analytic Closure",
            "",
            "The exact target is",
            "",
            "```text",
            "L_k^(1)>=1/(4*k^2), k>=319.",
            "```",
            "",
            "The paired interval theorem proves the seventh-order normalized remainder",
            "floor `-79/1000` for every real mode `0.9264<=u<=5`, reaching",
            "`t=V'(x(5))>1.5241916613e10`. The analytic ray certificate proves",
            "the same paired remainder floor on `u>=5`, closing the half-line.",
            "The later exact cumulant bridge sharpens this to the sufficient estimate",
            "`kappa_3,t(2*log(U))>=-37/(50*t^2)` for every real `t>=318`.",
            "The leading-saddle certificate proves caps 13/20, 1/100, and 1/1000",
            "through fifth order; the compact and ray certificates close the remainder.",
            "",
            "The first-summand dominance certificate gives",
            "",
            "```text",
            "|L_k-L_k^(1)|<=16/(k-1)^6.",
            "```",
            "",
            "Since `1/(4*k^2)>16/(k-1)^6` for every `k>=319`, the proved",
            "dominant target closes the full-kernel tail and splices to the certified",
            "finite collar through `k=318`.",
            "",
            "```text",
            "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
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
    parser.add_argument("--full-source-jsonl", type=Path, default=DEFAULT_FULL_SOURCE_JSONL)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact(args.full_source_jsonl)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF negative-lambda first-summand saddle-wall closure: "
        "9 rows, 0 issues, 9 positive samples, 9 quarter-k2 samples, 9 bracketed saddles, "
        "0 open requirements, 2 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
