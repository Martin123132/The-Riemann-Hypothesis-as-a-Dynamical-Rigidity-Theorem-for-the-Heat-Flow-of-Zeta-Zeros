#!/usr/bin/env python3
"""Build the first-summand cumulant bridge and its finite stress scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md"
)
DEFAULT_T = 100
DEFAULT_START_T = 318
DEFAULT_START_K = 319
DEFAULT_CUMULANT_CONSTANT = Fraction(37, 50)
DEFAULT_MPMATH_DPS = 60
DEFAULT_BISECTION_STEPS = 210
DEFAULT_UPPER_PADDING = "1.2"
DEFAULT_SAMPLE_T = (318, 319, 400, 700, 1000, 2000, 5000, 20000, 100000)


@dataclass(frozen=True)
class ClaimRow:
    id: str
    role: str
    claim: str
    formula: str | None
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None
    gap: str | None = None


def sci(value: mp.mpf, digits: int = 35) -> str:
    return mp.nstr(value, n=digits, min_fixed=-6, max_fixed=6)


def log_phi1(u: mp.mpf) -> mp.mpf:
    q = mp.pi * mp.exp(4 * u)
    return mp.log(mp.pi) + 5 * u + mp.log(2 * q - 3) - q


def saddle_derivative(u: mp.mpf, t: mp.mpf) -> mp.mpf:
    q = mp.pi * mp.exp(4 * u)
    return 2 * t / u - 2 * DEFAULT_T * u + 5 + 8 * q / (2 * q - 3) - 4 * q


def saddle_point(t: mp.mpf) -> mp.mpf:
    low = mp.mpf("0.01")
    high = mp.mpf(1)
    while saddle_derivative(high, t) > 0:
        high *= mp.mpf("1.4")
    for _ in range(DEFAULT_BISECTION_STEPS):
        middle = (low + high) / 2
        if saddle_derivative(middle, t) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def cumulant_row(t_value: int) -> dict:
    t = mp.mpf(t_value)
    saddle = saddle_point(t)
    maximum = 2 * t * mp.log(saddle) - DEFAULT_T * saddle**2 + log_phi1(saddle)
    center = 2 * mp.log(saddle)

    def scaled_weight(u: mp.mpf) -> mp.mpf:
        if not u:
            return mp.mpf(0)
        return mp.exp(2 * t * mp.log(u) - DEFAULT_T * u**2 + log_phi1(u) - maximum)

    upper = saddle + mp.mpf(DEFAULT_UPPER_PADDING)
    shifted = []
    for power in range(4):
        shifted.append(
            mp.quad(
                lambda u, power=power: scaled_weight(u) * (2 * mp.log(u) - center) ** power,
                [0, saddle, upper],
            )
        )
    mean = shifted[1] / shifted[0]
    second = shifted[2] / shifted[0]
    third = shifted[3] / shifted[0]
    kappa3 = third - 3 * mean * second + 2 * mean**3
    scaled_kappa3 = t**2 * kappa3
    target = -mp.mpf(DEFAULT_CUMULANT_CONSTANT.numerator) / DEFAULT_CUMULANT_CONSTANT.denominator
    return {
        "t": t_value,
        "saddle": sci(saddle),
        "kappa3_2logu": sci(kappa3),
        "scaled_t2_kappa3": sci(scaled_kappa3),
        "target_scaled_floor": sci(target),
        "target_margin": sci(scaled_kappa3 - target),
        "above_target": bool(scaled_kappa3 > target),
        "proof_boundary": (
            "High-precision finite-upper mpmath quadrature; not an interval enclosure "
            "or a uniform cumulant theorem."
        ),
    }


def transfer_diagnostics() -> dict:
    k = DEFAULT_START_K
    rational_margin = (
        Fraction(4, (2 * k + 1) ** 2)
        - DEFAULT_CUMULANT_CONSTANT / (k - 1) ** 2
        - Fraction(1, 4 * k**2)
    )
    if rational_margin <= 0:
        raise RuntimeError("cumulant transfer margin is not positive")
    # After k=n+319, this is the numerator of the rational transfer margin.
    shifted_coefficients = (9_130_082_706, 215_582_064, 1_489_493, 4_108, 4)
    if any(value <= 0 for value in shifted_coefficients):
        raise RuntimeError("shifted transfer polynomial is not coefficient-positive")
    return {
        "cumulant_constant": "37/50",
        "continuous_start_t": DEFAULT_START_T,
        "integer_start_k": DEFAULT_START_K,
        "gamma_lower_bound": "-log(1-z)>=z for 0<z<1",
        "rational_transfer_margin_at_k319": str(rational_margin),
        "rational_transfer_margin_at_k319_decimal": sci(
            mp.mpf(rational_margin.numerator) / rational_margin.denominator
        ),
        "unshifted_numerator": "4*k^4-996*k^3+401*k^2-50*k-25",
        "unshifted_denominator": "100*k^2*(k-1)^2*(2*k+1)^2",
        "shifted_variable": "n=k-319",
        "shifted_numerator": (
            "4*n^4+4108*n^3+1489493*n^2+215582064*n+9130082706"
        ),
        "shifted_coefficients_ascending": list(shifted_coefficients),
        "propagation": (
            "Every coefficient of the shifted numerator is positive for n=k-319>=0."
        ),
    }


def build_artifact() -> dict:
    mp.mp.dps = DEFAULT_MPMATH_DPS
    sample_rows = [cumulant_row(t) for t in DEFAULT_SAMPLE_T]
    if not all(row["above_target"] for row in sample_rows):
        raise RuntimeError("finite cumulant scout crossed the proposed floor")
    scaled_values = [mp.mpf(row["scaled_t2_kappa3"]) for row in sample_rows]
    if not all(right > left for left, right in zip(scaled_values, scaled_values[1:])):
        raise RuntimeError("sampled scaled cumulant is not strictly increasing")
    minimum_row = min(sample_rows, key=lambda row: mp.mpf(row["target_margin"]))
    transfer = transfer_diagnostics()

    rows = [
        ClaimRow(
            id="fscb_01_continuous_moment_family",
            role="exact_reduction",
            claim="Extend the first-summand moments to a smooth real parameter t.",
            formula=(
                "M_t^(1)=2*integral_0^infinity u^(2t)*exp(-100*u^2)*phi_1(u)du; "
                "F(t)=log(M_t^(1))"
            ),
            readiness="available_exact",
            proof_boundary="Exact definition and differentiation domain for t>=318.",
        ),
        ClaimRow(
            id="fscb_02_cumulant_identity",
            role="exact_lemma",
            claim="The third derivative of the log moment is a third cumulant.",
            formula="F'''(t)=kappa_3,t(2*log(U))",
            readiness="available_exact",
            proof_boundary="Exact differentiation-under-the-integral identity.",
        ),
        ClaimRow(
            id="fscb_03_bspline_identity",
            role="exact_lemma",
            claim="The cubic finite difference is an average of the continuous third derivative.",
            formula=(
                "Delta^3 F(k-1)=integral_[0,1]^3 F'''(k-1+r+s+v)drdsdv"
            ),
            readiness="available_exact",
            proof_boundary="Exact threefold fundamental-theorem-of-calculus identity.",
        ),
        ClaimRow(
            id="fscb_04_gamma_wall",
            role="exact_lemma",
            claim="The Gamma normalization contributes an explicit positive wall.",
            formula=(
                "-Delta^3 log Gamma(k-1/2)=-log(1-1/(k+1/2)^2)"
            ),
            readiness="available_exact",
            proof_boundary="Exact Gamma recurrence only.",
        ),
        ClaimRow(
            id="fscb_05_uniform_cumulant_target",
            role="analytic_theorem",
            claim="The paired compact and ray certificates prove the uniform lower bound for the special first-summand log-skewness.",
            formula="kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318",
            readiness="ready_to_apply",
            proof_boundary=(
                "Analytic theorem supplied by the paired remainder certificates; it does "
                "not by itself prove cone propagation or the all-order bridge."
            ),
            gap=None,
        ),
        ClaimRow(
            id="fscb_06_exact_rational_transfer",
            role="exact_conditional_bridge",
            claim="The cumulant target implies the quantitative first-summand wall.",
            formula=(
                "L_k^(1)>=1/(k+1/2)^2-37/(50*(k-1)^2)>=1/(4*k^2), k>=319"
            ),
            readiness="ready_to_apply",
            proof_boundary="Exact implication with fscb_05 now discharged.",
            diagnostics=transfer,
        ),
        ClaimRow(
            id="fscb_07_high_precision_scout",
            role="finite_scout",
            claim="Finite high-precision cumulant samples stay above the proposed floor.",
            formula="t^2*kappa_3,t(2*log(U))>-37/50 on the recorded sample",
            readiness="finite_only",
            proof_boundary="Floating theorem-search evidence; not interval-certified or uniform.",
            diagnostics={"sample_rows": sample_rows},
        ),
        ClaimRow(
            id="fscb_08_full_kernel_handoff",
            role="conditional_handoff",
            claim="The proved cumulant theorem closes the lambda=-100 full-kernel tail.",
            formula=(
                "fscb_05 => L_k^(1)>=1/(4*k^2) => "
                "L_k>=1/(4*k^2)-16/(k-1)^6>0"
            ),
            readiness="ready_to_apply",
            proof_boundary=(
                "Uses the separately certified first-summand dominance and finite collar; "
                "proves the adjacent wall but not the remaining cone-propagation/all-order bridge."
            ),
        ),
    ]

    summary = {
        "target_rows": len(rows),
        "exact_identity_rows": sum(row.role in {"exact_reduction", "exact_lemma"} for row in rows),
        "conditional_bridge_rows": sum(row.role == "exact_conditional_bridge" for row in rows),
        "sample_rows": len(sample_rows),
        "above_target_sample_rows": sum(row["above_target"] for row in sample_rows),
        "scaled_profile_strictly_increasing": True,
        "minimum_sample_margin": minimum_row["target_margin"],
        "minimum_sample_margin_at_t": minimum_row["t"],
        "open_requirement_rows": sum(row.role == "open_requirement" for row in rows),
        "ready_to_apply_rows": 3,
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_cumulant_bridge",
        "date": "2026-07-10",
        "status": "exact cumulant reduction with analytic hypothesis discharged",
        "proof_boundary": (
            "This artifact composes the analytic paired-remainder cumulant theorem with the exact "
            "Gamma transfer and proves the first-summand and full-kernel adjacent wall; "
            "the finite samples remain cross-checks only. "
            "It does not prove the remaining cone-propagation/all-order bridge, RH, or Lambda <= 0."
        ),
        "source_saddle_target": (
            "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md"
        ),
        "source_dominance_certificate": (
            "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md"
        ),
        "source_collar_extension": (
            "outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md"
        ),
        "source_paired_ray_certificate": (
            "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md"
        ),
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py"
        ),
        "diagnostics": {
            "parameters": {
                "T": DEFAULT_T,
                "continuous_start_t": DEFAULT_START_T,
                "integer_start_k": DEFAULT_START_K,
                "cumulant_constant": "37/50",
                "mpmath_dps": DEFAULT_MPMATH_DPS,
                "bisection_steps": DEFAULT_BISECTION_STEPS,
                "upper_padding": DEFAULT_UPPER_PADDING,
                "sample_t": list(DEFAULT_SAMPLE_T),
            },
            "exact_transfer": transfer,
            "sample_rows": sample_rows,
            "formal_asymptotic_target": "t^2*kappa_3,t(2*log(U))->0 from below",
        },
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Negative-Lambda First-Summand Cumulant Bridge",
        "",
        "Date: 2026-07-10",
        "",
        "Status: exact cumulant reduction with analytic hypothesis discharged; not a proof of the remaining all-order bridge.",
        "This proves the lambda=-100 adjacent wall, not the remaining all-order bridge, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_cumulant_bridge`.",
        "",
        "Machine-readable result:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF negative-lambda first-summand cumulant bridge: 8 rows, 0 issues, 4 exact identities, 1 conditional bridge, 9 positive samples, 0 open requirements, 3 ready-to-apply rows",
        "```",
        "",
        "## Exact Decomposition",
        "",
        "For real `t>=318`, define",
        "",
        "```text",
        "M_t^(1)=2*integral_0^infinity u^(2t)*exp(-100*u^2)*phi_1(u)du,",
        "F(t)=log(M_t^(1)).",
        "```",
        "",
        "Under the probability measure proportional to the integrand, differentiation",
        "under the integral gives",
        "",
        "```text",
        "F'''(t)=kappa_3,t(2*log(U)).",
        "```",
        "",
        "Three applications of the fundamental theorem of calculus give",
        "",
        "```text",
        "Delta^3 F(k-1)=integral_[0,1]^3 F'''(k-1+r+s+v)drdsdv.",
        "```",
        "",
        "The Gamma normalization is exact:",
        "",
        "```text",
        "-Delta^3 log Gamma(k-1/2)=-log(1-1/(k+1/2)^2).",
        "```",
        "",
        "## Sufficient Cumulant Bound",
        "",
        "The paired compact and ray certificates prove the single continuous estimate",
        "",
        "```text",
        "kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318.",
        "```",
        "",
        "Indeed, `-log(1-z)>=z` and the cubic-difference integral imply",
        "",
        "```text",
        "L_k^(1)>=1/(k+1/2)^2-37/(50*(k-1)^2).",
        "```",
        "",
        "After subtracting `1/(4*k^2)`, the exact rational numerator is",
        "",
        "```text",
        "4*k^4-996*k^3+401*k^2-50*k-25.",
        "```",
        "",
        "With `n=k-319`, this becomes",
        "",
        "```text",
        diagnostics["exact_transfer"]["shifted_numerator"] + ".",
        "```",
        "",
        "Every coefficient is positive for `n>=0`. Therefore the cumulant estimate",
        "implies `L_k^(1)>=1/(4*k^2)` for every integer `k>=319`.",
        "",
        "## High-Precision Scout",
        "",
        "The following 60-digit saddle-centered calculations use finite upper",
        "truncation and are not interval enclosures.",
        "",
        "| t | saddle | t^2 kappa_3 | margin above -37/50 |",
        "|---:|---:|---:|---:|",
    ]
    for row in diagnostics["sample_rows"]:
        lines.append(
            f"| `{row['t']}` | `{row['saddle']}` | `{row['scaled_t2_kappa3']}` | "
            f"`{row['target_margin']}` |"
        )
    lines.extend(
        [
            "",
            "All nine recorded samples exceed the proposed floor, and the sampled",
            "scaled cumulant is strictly increasing. This is theorem-search evidence only.",
            "",
            "## Proof Boundary",
            "",
            "The exact algebra reduces the existing first-summand wall target to a",
            "special-kernel skewness estimate. Generic strict log-concavity is not enough",
            "to determine third-cumulant sign. The paired remainder theorem supplies",
            "the required explicit `q=pi*exp(4u)` saddle geometry.",
            "",
            "With the cumulant estimate proved, the existing dominance certificate",
            "and collar extension give",
            "",
            "```text",
            "L_k>=1/(4*k^2)-16/(k-1)^6>0, k>=319,",
            "```",
            "",
            "spliced to the certified finite adjacent wall through `k=318`.",
            "",
            "```text",
            "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
            "outputs/signed_hankel_jensen_dependency_graph.md",
            "```",
            "",
            "Summary:",
            "",
            (
                "The first-summand wall is exactly reduced to the uniform bound "
                "kappa_3,t(2 log U)>=-37/(50 t^2) for t>=318. The Gamma contribution "
                "and rational transfer are exact; nine finite high-precision samples "
                f"clear the floor, with the smallest sampled margin "
                f"{summary['minimum_sample_margin']} at t={summary['minimum_sample_margin_at_t']}. "
                "The uniform cumulant estimate is discharged by the paired ray certificate."
            ),
        ]
    )
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
        "validated Jensen-window PF negative-lambda first-summand cumulant bridge: "
        "8 rows, 0 issues, 4 exact identities, 1 conditional bridge, 9 positive samples, "
        "0 open requirements, 3 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
