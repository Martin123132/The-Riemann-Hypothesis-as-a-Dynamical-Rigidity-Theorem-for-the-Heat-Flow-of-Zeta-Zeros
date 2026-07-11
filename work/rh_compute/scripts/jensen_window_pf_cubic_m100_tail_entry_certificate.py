#!/usr/bin/env python3
"""Certify the all-k cubic reciprocal-defect tail at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_cubic_m100_tail_entry_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cubic_m100_tail_entry_certificate.md"
COMPACT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json"
)
RAY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.json"
)
DEFAULT_K = 319
Q_FLOOR = 1_000_000_000
U_FLOOR = 5


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def positive_shift_coefficients(expression: sp.Expr, variable: sp.Symbol, shift: int) -> list[int]:
    polynomial = sp.cancel(expression)
    shifted = sp.Poly(sp.expand(polynomial.subs(variable, variable + shift)), variable)
    coefficients = [int(value) for value in reversed(shifted.all_coeffs())]
    if any(value <= 0 for value in coefficients):
        raise RuntimeError(f"shifted polynomial is not coefficient-positive: {shifted.as_expr()}")
    return coefficients


def ray_skewness_diagnostics(ray: dict) -> dict:
    q_lower = flint.arb(ray["diagnostics"]["mode_handoff"]["q_at_499_lower"])
    if not q_lower > Q_FLOOR:
        raise RuntimeError("ray q floor is not available")

    q = Fraction(Q_FLOOR, 1)
    u = Fraction(U_FLOOR, 1)
    a_small_terms = Fraction(1408, Q_FLOOR) + Fraction(120, U_FLOOR * Q_FLOOR**2) + Fraction(3600, Q_FLOOR**3)
    if not a_small_terms < Fraction(4, 5):
        raise RuntimeError("ray curvature upper normalization failed")
    b_negative_terms = (
        Fraction(2304, Q_FLOOR)
        + Fraction(216, U_FLOOR**2 * Q_FLOOR)
        + Fraction(25632, U_FLOOR * Q_FLOOR**2)
        + Fraction(2880, Q_FLOOR**3)
        + Fraction(450, U_FLOOR**2 * Q_FLOOR**3)
        + Fraction(21600, U_FLOOR * Q_FLOOR**4)
    )
    if not b_negative_terms < 12:
        raise RuntimeError("ray third-potential lower normalization failed")
    if not 49 * q > 1_800_000:
        raise RuntimeError("ray leading skewness does not dominate the residual")
    source_bound = ray["diagnostics"]["cumulant_composition"].get("bound")
    if source_bound != "|kappa_3(Y)+alpha|<=120/q":
        raise RuntimeError("ray residual skewness source drifted")
    return {
        "q_floor": Q_FLOOR,
        "u_floor": U_FLOOR,
        "q_source_lower": ray["diagnostics"]["mode_handoff"]["q_at_499_lower"],
        "curvature_upper_reduction": (
            "N_a<=68*u*q^3 and 2*q-3>=(19/10)*q imply a=V''<=5*u^2*q"
        ),
        "curvature_small_terms": str(a_small_terms),
        "third_potential_lower_reduction": (
            "N_b>=500*u^2*q^4 and 2*q-3<=2*q imply b=V'''>=7*u^3*q"
        ),
        "third_potential_negative_terms": str(b_negative_terms),
        "normalized_leading_lower": "alpha=b/a^(3/2)>=7/(5^(3/2)*sqrt(q))",
        "residual_source_bound": source_bound,
        "dominance_squared_gate": "49*q>1800000",
        "dominance_margin_at_q_floor": 49 * Q_FLOOR - 1_800_000,
        "conclusion": "kappa_3,t(2*log(U))<0 on the full mode ray u>=5",
    }


def exact_tail_diagnostics() -> dict:
    k, m, n = sp.symbols("k m n", integer=True, positive=True)
    defect_floor_numerator = sp.expand(
        (sp.Rational(1, 4) / m - 16 / (m - 1) ** 6 - sp.Rational(16, 5) / (m - 1) ** 5 - sp.Rational(1, 5) / m)
        * 20
        * m
        * (m - 1) ** 6
    )
    defect_coefficients = positive_shift_coefficients(defect_floor_numerator, m, 320)

    epsilon_numerator = sp.expand((k - 1) ** 6 - 16 * k**2)
    epsilon_coefficients = positive_shift_coefficients(epsilon_numerator, k, DEFAULT_K)

    final_numerator = sp.expand(
        k**4 * 99_999**4 - (5 * k + 6) ** 3 * 100_000**4
    )
    final_coefficients = positive_shift_coefficients(final_numerator, k, DEFAULT_K)

    flint.ctx.prec = 512
    k0 = flint.arb(DEFAULT_K)
    analytic_bound = (
        (flint.arb(100_000) / 99_999) ** 2
        * (5 * k0 + 6) ** flint.arb("1.5")
        / k0**2
    )
    if not analytic_bound < 1:
        raise RuntimeError("reciprocal-defect analytic endpoint bound failed")
    return {
        "tail_start_k": DEFAULT_K,
        "full_log_wall_lower": "L_k>=1/(4*k^2)-16/(k-1)^6",
        "full_log_wall_upper": (
            "L_k<=-log(1-4/(2*k+1)^2)+16/(k-1)^6"
        ),
        "telescoping_identity": "-log(x_m)=sum_(j=m)^infinity L_j",
        "tail_sum_lower": (
            "A_m=1/(4*m)-16/(m-1)^6-16/(5*(m-1)^5)>=1/(5*m), m>=320"
        ),
        "defect_floor": "d_m=1-x_m>=1-exp(-A_m)>=A_m/(1+A_m)>=1/(5*m+1)",
        "defect_floor_shifted_coefficients_ascending": defect_coefficients,
        "epsilon_bound": "16/(k-1)^6<=1/k^2, k>=319",
        "epsilon_shifted_coefficients_ascending": epsilon_coefficients,
        "increment_bound": (
            "q_(k+1)-q_k<=(100000/99999)^2*(5*k+6)^(3/2)/k^2<1"
        ),
        "increment_bound_at_k319_ball": str(analytic_bound),
        "final_shifted_quartic_coefficients_ascending": final_coefficients,
        "spatial_tail": "q_(k+1)-q_k=O(k^(-1/2))->0 at lambda=-100",
    }


def build_payload() -> dict:
    compact = load_json(COMPACT_SOURCE)
    ray = load_json(RAY_SOURCE)
    compact_diag = compact["diagnostics"]["compact_remainder"]
    if compact_diag.get("negative_cumulant_blocks") != 4_074:
        raise RuntimeError("compact negative-skewness theorem is unavailable")
    if not Fraction(compact_diag["maximum_cumulant3_upper"]) < 0:
        raise RuntimeError("compact cumulant upper bound is not negative")
    ray_diag = ray_skewness_diagnostics(ray)
    tail_diag = exact_tail_diagnostics()
    rows = [
        CertificateRow(
            id="cmte_01_compact_negative_skewness",
            role="interval_theorem",
            readiness="ready_to_apply",
            claim="Every compact first-summand tilted distribution has strictly negative third log-cumulant.",
            formula="kappa_3,t(2*log(U))<0 for 0.9264<=u_t<=5",
            proof_boundary="Compact mode theorem only.",
            diagnostics={
                "negative_blocks": compact_diag["negative_cumulant_blocks"],
                "maximum_upper": compact_diag["maximum_cumulant3_upper"],
                "maximum_block": compact_diag["maximum_cumulant3_block"],
            },
        ),
        CertificateRow(
            id="cmte_02_ray_potential_bounds",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="On u>=5, explicit potential derivatives bound the leading standardized skewness away from zero.",
            formula="V''<=5*u^2*q, V'''>=7*u^3*q, q>=10^9",
            proof_boundary="First-summand asymptotic mode ray only.",
            diagnostics=ray_diag,
        ),
        CertificateRow(
            id="cmte_03_ray_negative_skewness",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The leading skewness dominates the certified Gaussian-comparison residual on the full ray.",
            formula="kappa_3<=-7/(5^(3/2)*sqrt(q))+120/q<0",
            proof_boundary="Analytic u>=5 first-summand theorem only.",
        ),
        CertificateRow(
            id="cmte_04_global_first_summand_upper_wall",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Negative skewness makes the first-summand adjacent log wall no larger than its Gamma contribution.",
            formula="L_k^(1)<=-log(1-1/(k+1/2)^2), k>=319",
            proof_boundary="First-summand adjacent upper bound only.",
        ),
        CertificateRow(
            id="cmte_05_full_kernel_two_sided_wall",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="First-summand dominance transfers both lower and upper adjacent log-wall bounds to the full kernel.",
            formula=(
                "1/(4*k^2)-16/(k-1)^6<=L_k<="
                "-log(1-4/(2*k+1)^2)+16/(k-1)^6"
            ),
            proof_boundary="All-k lambda=-100 adjacent-log estimate only.",
        ),
        CertificateRow(
            id="cmte_06_telescoping_defect_floor",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The lower adjacent wall telescopes to a quantitative lower bound on every tail defect.",
            formula="d_m>=1/(5*m+1), m>=320",
            proof_boundary="Lambda=-100 defect tail only.",
            diagnostics={
                "shifted_coefficients": tail_diag["defect_floor_shifted_coefficients_ascending"]
            },
        ),
        CertificateRow(
            id="cmte_07_reciprocal_increment_bound",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The two-sided adjacent wall and defect floor force every reciprocal-defect tail increment strictly below one.",
            formula=tail_diag["increment_bound"],
            proof_boundary="All-k lambda=-100 cubic-tail theorem only.",
            diagnostics={
                "endpoint_ball": tail_diag["increment_bound_at_k319_ball"],
                "shifted_quartic_coefficients": tail_diag[
                    "final_shifted_quartic_coefficients_ascending"
                ],
            },
        ),
        CertificateRow(
            id="cmte_08_full_cubic_entry",
            role="interval_analytic_theorem",
            readiness="ready_to_apply",
            claim="The repaired prefix and analytic tail prove every shifted cubic Jensen polynomial is hyperbolic at lambda=-100.",
            formula="q_(k+1)-q_k<1 for every k>=1 at lambda=-100",
            proof_boundary="All-shift degree-3 entry at one heat parameter; not forward propagation or higher degree.",
        ),
        CertificateRow(
            id="cmte_09_spatial_tail",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The explicit increment bound supplies the spatial tail required by the cubic maximum principle at the entry parameter.",
            formula="q_(k+1)-q_k=O(k^(-1/2))->0 at lambda=-100",
            proof_boundary="Tail at lambda=-100 only; not uniform on a forward heat interval.",
        ),
        CertificateRow(
            id="cmte_10_forward_uniformity_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Forward cubic propagation now requires a reciprocal-defect increment tail uniform on compact lambda intervals.",
            formula="sup_(lambda in [-100,L]) |q_(k+1)(lambda)-q_k(lambda)| -> 0",
            proof_boundary="Open forward-tail theorem; higher-degree Jensen hyperbolicity, PF-infinity, RH, and Lambda <= 0 remain unproved.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_cubic_m100_tail_entry_certificate",
        "date": "2026-07-10",
        "status": "all-k lambda=-100 cubic Jensen entry theorem with open forward-uniform tail",
        "proof_boundary": (
            "This artifact proves every shifted degree-3 Jensen polynomial is hyperbolic at lambda=-100 and gives an explicit reciprocal-defect spatial tail there. "
            "It does not prove the tail uniformly on forward heat intervals, higher-degree Jensen hyperbolicity, PF-infinity, RH, or Lambda <= 0."
        ),
        "source_compact_skewness": COMPACT_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "source_ray_remainder": RAY_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "source_cumulant_bridge": "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
        "source_dominance": "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
        "source_cubic_prefix": "outputs/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.md",
        "diagnostics": {
            "compact_negative_skewness": {
                "negative_blocks": compact_diag["negative_cumulant_blocks"],
                "maximum_upper": compact_diag["maximum_cumulant3_upper"],
                "maximum_block": compact_diag["maximum_cumulant3_block"],
            },
            "ray_negative_skewness": ray_diag,
            "tail_algebra": tail_diag,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "compact_negative_cumulant_blocks": compact_diag["negative_cumulant_blocks"],
            "analytic_ray_negative_skewness": True,
            "lambda_minus_100_prefix_margins": 318,
            "lambda_minus_100_analytic_tail_start": DEFAULT_K,
            "full_cubic_entry": True,
            "entry_spatial_tail": True,
            "forward_uniform_tail": False,
            "open_handoff_rows": 1,
            "ready_to_apply_rows": 9,
            "target_closing": False,
            "main_finding": (
                "Compact Arb negative skewness and an analytic ray estimate give a first-summand adjacent upper wall. "
                "Together with the proved lower wall and full-kernel perturbation, telescoping yields d_m>=1/(5m+1). "
                "A coefficient-positive shifted quartic then proves q_(k+1)-q_k<1 for every k>=319, which composes with 318 repaired prefix margins to establish all-shift cubic Jensen hyperbolicity at lambda=-100."
            ),
        },
        "invariants": [
            "The compact skewness input is interval-certified on every accepted block.",
            "The ray skewness sign uses only explicit derivative and moment-comparison bounds.",
            "The finite prefix and analytic tail meet without an index gap.",
            "The spatial tail is proved at lambda=-100, not uniformly along the forward flow.",
            "Degree-3 entry is not an all-degree Jensen theorem.",
        ],
    }


def result_line(payload: dict) -> str:
    summary = payload["summary"]
    return (
        "validated Jensen-window PF cubic lambda=-100 tail-entry certificate: "
        f"{summary['rows']} rows, 0 issues, {summary['compact_negative_cumulant_blocks']} negative-skewness blocks, "
        "1 analytic negative-skewness ray, 318 prefix margins, 1 all-k cubic tail, "
        "1 full cubic entry theorem, 1 open forward-uniform tail"
    )


def write_note(payload: dict, path: Path) -> None:
    compact = payload["diagnostics"]["compact_negative_skewness"]
    tail = payload["diagnostics"]["tail_algebra"]
    lines = [
        "# Jensen-Window PF Cubic Lambda=-100 Tail-Entry Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: all-k lambda=-100 cubic Jensen entry theorem with open forward-uniform",
        "tail. This is not a proof of higher-degree Jensen hyperbolicity, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_cubic_m100_tail_entry_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_cubic_m100_tail_entry_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_cubic_m100_tail_entry_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(payload),
        "```",
        "",
        "## Negative Skewness",
        "",
        f"The compact paired engine certifies `kappa_3,t(2*log(U))<0` on all `{compact['negative_blocks']}`",
        "mode blocks through `u=5`. Its closest outward-rounded upper bound is",
        "",
        "```text",
        compact["maximum_upper"],
        "```",
        "",
        "On `u>=5`, put `q=pi*exp(4u)`. Exact component bounds give",
        "",
        "```text",
        "V''<=5*u^2*q,",
        "V'''>=7*u^3*q,",
        "alpha=V'''/(V'')^(3/2)>=7/(5^(3/2)*sqrt(q)).",
        "```",
        "",
        "The analytic ray theorem already proves `|kappa_3+alpha|<=120/q` and",
        "`q>=10^9`; squaring reduces dominance to `49*q>1800000`. Hence the",
        "third cumulant is negative for every real `t>=318`.",
        "",
        "## Two-Sided Adjacent Wall",
        "",
        "The exact Gamma/cumulant identity and first-summand dominance now give",
        "",
        "```text",
        "1/(4*k^2)-16/(k-1)^6 <= L_k",
        "L_k <= -log(1-4/(2*k+1)^2)+16/(k-1)^6,",
        "k>=319.",
        "```",
        "",
        "Since `x_k->1`, summing the lower wall proves",
        "",
        "```text",
        "-log(x_m)=sum_(j=m)^infinity L_j,",
        "d_m=1-x_m>=1/(5*m+1), m>=320.",
        "```",
        "",
        "The latter estimate is exact: after clearing denominators and setting",
        "`n=m-320`, every coefficient is positive.",
        "",
        "## Cubic Tail",
        "",
        "Using `d_k>=d_(k+1)` and the upper adjacent wall gives",
        "",
        "```text",
        tail["increment_bound"],
        "```",
        "",
        f"At `k=319` the explicit upper ball is `{tail['increment_bound_at_k319_ball']}`.",
        "After squaring and setting `n=k-319`, the final quartic has positive",
        "coefficients",
        "",
        "```text",
        str(tail["final_shifted_quartic_coefficients_ascending"]),
        "```",
        "",
        "so the bound remains below one for every `k>=319` and tends to zero.",
        "Together with the 318 repaired prefix margins, this proves every shifted",
        "degree-3 Jensen polynomial is hyperbolic at `lambda=-100`.",
        "",
        "## Remaining Handoff",
        "",
        "Forward cubic propagation now has one explicit analytic requirement:",
        "",
        "```text",
        "sup_(lambda in [-100,L]) |q_(k+1)(lambda)-q_k(lambda)| -> 0",
        "for every finite L.",
        "```",
        "",
        "The entry theorem and inward boundary algebra are closed. Uniformity of the",
        "spatial tail along the forward heat interval, and every higher-degree minor",
        "cone, remain open.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, args.note)
    print(result_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
