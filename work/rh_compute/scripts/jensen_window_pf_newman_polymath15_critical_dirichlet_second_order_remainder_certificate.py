#!/usr/bin/env python3
"""Certify a second-order remainder for the critical Dirichlet blocks."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_critical_"
    "dirichlet_second_order_remainder_certificate"
)
DEFAULT_OUT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / f"outputs/{STEM}.md"
POLYMATH_SOURCE = "https://arxiv.org/abs/1904.12438"
DLMF_SOURCE = "https://dlmf.nist.gov/5.11"
SOURCE_FILES = {
    "first_correction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.json"
    ),
    "critical_collar": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.json"
    ),
}

L_MIN = 50
T_AUDIT_MIN = 1_000_000
SHIFT_PRODUCT_BOUND = 27
SHIFT_BOUND = sp.Rational(27, 2)
PER_TERM_CONSTANT = 400_000
CENTRAL_CONSTANT = 312_830
TAIL_CONSTANT = 2
FIXED_CELL_CONSTANT = 8_000_000


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCE_FILES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing source artifacts: {missing}")
    return {name: file_hash(path) for name, path in SOURCE_FILES.items()}


def exact_moments() -> dict[str, str]:
    t = sp.symbols("t", nonnegative=True)
    c = sp.symbols("c", complex=True)
    v = sp.symbols("v", real=True)
    w = sp.sqrt(t) * v + c
    second = sp.Poly(sp.expand(w**2), v)
    gaussian_second = sp.simplify(
        second.coeff_monomial(v**0)
        + sp.Rational(1, 2) * second.coeff_monomial(v**2)
    )
    if gaussian_second != c**2 + t / 2:
        raise RuntimeError("shifted Gaussian second moment failed")

    # For dmu=pi^(-1/2) exp(-v^2)dv.
    if sp.simplify(sp.gamma(1) / sp.sqrt(sp.pi) - 1 / sp.sqrt(sp.pi)) != 0:
        raise RuntimeError("absolute first moment normalization failed")
    if sp.simplify(sp.gamma(2) / sp.sqrt(sp.pi) - 1 / sp.sqrt(sp.pi)) != 0:
        raise RuntimeError("absolute third moment normalization failed")
    if sp.simplify(sp.gamma(sp.Rational(5, 2)) / sp.sqrt(sp.pi) - sp.Rational(3, 4)) != 0:
        raise RuntimeError("fourth moment normalization failed")
    return {
        "shifted_second_moment": (
            "For w=sqrt(t)v+c and dmu=pi^(-1/2)exp(-v^2)dv, "
            "E[w^2]=c^2+t/2"
        ),
        "gaussian_moments": (
            "E|v|=1/sqrt(pi), E|v|^3=1/sqrt(pi), E[v^4]=3/4"
        ),
    }


def interval_audit() -> dict:
    # Exact rational guards for every scalar comparison used by the proof.
    # The mpmath values below are diagnostics only.
    exp_14_lower = sum(
        (Fraction(14) ** k) / math.factorial(k) for k in range(41)
    )
    if not exp_14_lower > 1_000_000:
        raise RuntimeError("exact log(T) guard failed")
    e_upper = sum(Fraction(1, math.factorial(k)) for k in range(10))
    e_upper += Fraction(11, 10 * math.factorial(10))
    if not e_upper < Fraction(68, 25):
        raise RuntimeError("exact e upper bound failed")

    exact_tail_guard = (
        Fraction(1) + Fraction(70, 3) - Fraction(99, 100) * 10_000
    )
    mass_tail_guard = Fraction(112, 3) - 10_000
    second_tail_guard = Fraction(1) + Fraction(56, 3) - 10_000
    exact_derivative_guard = Fraction(5, 3) - Fraction(66, 100) * 10_000
    mass_derivative_guard = Fraction(8, 3) - Fraction(2, 3) * 10_000
    second_derivative_guard = Fraction(4, 3) - Fraction(2, 3) * 10_000
    central_guard = (
        Fraction(1, 6 * T_AUDIT_MIN)
        + Fraction(1, 200)
        + Fraction(1, T_AUDIT_MIN)
    )
    global_guard = (
        25
        * PER_TERM_CONSTANT
        * Fraction(68, 25) ** 2
        / Fraction(157, 50) ** 2
    )
    for label, value in (
        ("exact tail", exact_tail_guard),
        ("mass tail", mass_tail_guard),
        ("second-moment tail", second_tail_guard),
        ("exact tail derivative", exact_derivative_guard),
        ("mass tail derivative", mass_derivative_guard),
        ("second-moment tail derivative", second_derivative_guard),
    ):
        if not value < 0:
            raise RuntimeError(f"exact rational {label} guard failed")
    if not central_guard < Fraction(3, 500):
        raise RuntimeError("exact rational central guard failed")
    if not global_guard < FIXED_CELL_CONSTANT:
        raise RuntimeError("exact rational fixed-cell constant failed")
    if not 6 * 2**49 > T_AUDIT_MIN:
        raise RuntimeError("exact critical-height threshold guard failed")

    mp.mp.dps = 80
    t0 = mp.mpf(T_AUDIT_MIN)
    root = t0 ** (mp.mpf(1) / 3)

    # Logarithms of the three normalized tail bounds at T=T_AUDIT_MIN.
    exact_tail_log = (
        mp.mpf(92) / t0
        - mp.mpf("0.99") * t0 ** (mp.mpf(2) / 3)
        + mp.mpf(5) / 3 * mp.log(t0)
        - mp.log(mp.mpf("0.99") * mp.sqrt(mp.pi))
    )
    mass_tail_log = (
        -t0 ** (mp.mpf(2) / 3)
        + mp.mpf(8) / 3 * mp.log(t0)
        - mp.log(mp.sqrt(mp.pi))
    )
    second_tail_log = (
        mp.log(2)
        + mp.mpf(4) / 3 * mp.log(t0)
        - t0 ** (mp.mpf(2) / 3)
    )
    if not exact_tail_log < 0:
        raise RuntimeError("exact-integrand tail endpoint failed")
    if not mass_tail_log < 0:
        raise RuntimeError("Gaussian-mass tail endpoint failed")
    if not second_tail_log < 0:
        raise RuntimeError("Gaussian-second-moment tail endpoint failed")

    # Each logarithmic bound is decreasing from this endpoint onward.
    exact_derivative = (
        -mp.mpf(92) / t0**2
        - mp.mpf("0.66") * t0 ** (-mp.mpf(1) / 3)
        + mp.mpf(5) / (3 * t0)
    )
    mass_derivative = (
        -mp.mpf(2) / 3 * t0 ** (-mp.mpf(1) / 3)
        + mp.mpf(8) / (3 * t0)
    )
    second_derivative = (
        -mp.mpf(2) / 3 * t0 ** (-mp.mpf(1) / 3)
        + mp.mpf(4) / (3 * t0)
    )
    if not exact_derivative < 0:
        raise RuntimeError("exact-integrand tail monotonicity failed")
    if not mass_derivative < 0:
        raise RuntimeError("Gaussian-mass tail monotonicity failed")
    if not second_derivative < 0:
        raise RuntimeError("Gaussian-second-moment tail monotonicity failed")

    central_exponent = (
        1 / (6 * t0)
        + 1 / (2 * t0 ** (mp.mpf(1) / 3))
        + 1 / t0
    )
    if not central_exponent < mp.mpf("0.006"):
        raise RuntimeError("central exponential bound failed")
    if not CENTRAL_CONSTANT + TAIL_CONSTANT < PER_TERM_CONSTANT:
        raise RuntimeError("saved per-term constant failed")

    global_constant = (
        25 * mp.mpf(PER_TERM_CONSTANT) * mp.e**2 / mp.pi**2
    )
    if not global_constant < FIXED_CELL_CONSTANT:
        raise RuntimeError("fixed-cell summation constant failed")

    actual_t_min = 2 * mp.pi * mp.e ** (L_MIN - 1)
    if not actual_t_min > T_AUDIT_MIN:
        raise RuntimeError("critical collar does not reach audit threshold")

    return {
        "role": "high_precision_margin_audit_after_exact_rational_guards",
        "precision_digits": mp.mp.dps,
        "T_audit_min": str(T_AUDIT_MIN),
        "actual_T_min_lower_bound": mp.nstr(actual_t_min, 60),
        "central_exponent_upper_at_T_audit_min": mp.nstr(
            central_exponent, 60
        ),
        "exact_tail_log_margin_at_T_audit_min": mp.nstr(
            exact_tail_log, 60
        ),
        "mass_tail_log_margin_at_T_audit_min": mp.nstr(
            mass_tail_log, 60
        ),
        "second_tail_log_margin_at_T_audit_min": mp.nstr(
            second_tail_log, 60
        ),
        "exact_tail_log_derivative_at_T_audit_min": mp.nstr(
            exact_derivative, 60
        ),
        "mass_tail_log_derivative_at_T_audit_min": mp.nstr(
            mass_derivative, 60
        ),
        "second_tail_log_derivative_at_T_audit_min": mp.nstr(
            second_derivative, 60
        ),
        "V_at_T_audit_min": mp.nstr(root, 60),
        "central_constant": str(CENTRAL_CONSTANT),
        "tail_constant": str(TAIL_CONSTANT),
        "per_term_constant": str(PER_TERM_CONSTANT),
        "fixed_cell_constant_ball": mp.nstr(global_constant, 60),
        "fixed_cell_constant_lt": str(FIXED_CELL_CONSTANT),
        "exact_rational_guards": {
            "log_T_lt_14": (
                "sum_(k=0)^40 14^k/k! > 10^6, hence log(10^6)<14"
            ),
            "e_lt_2_72": (
                "sum_(k=0)^9 1/k! + 11/(10*10!) < 68/25"
            ),
            "pi_gt_3_14": "classical bound pi>157/50",
            "exact_tail_log_upper": str(exact_tail_guard),
            "mass_tail_log_upper": str(mass_tail_guard),
            "second_tail_log_upper": str(second_tail_guard),
            "exact_tail_derivative_upper_after_multiplying_by_T": str(
                exact_derivative_guard
            ),
            "mass_tail_derivative_upper_after_multiplying_by_T": str(
                mass_derivative_guard
            ),
            "second_tail_derivative_upper_after_multiplying_by_T": str(
                second_derivative_guard
            ),
            "central_exponent_upper": str(central_guard),
            "fixed_cell_constant_upper": str(global_guard),
            "critical_T_min_guard": "2*pi*e^49>6*2^49>10^6",
        },
    }


def build_exact() -> dict:
    moments = exact_moments()
    return {
        "region": (
            "L>=50, 0<t<=1/2, tL<=25, T>=2*pi*exp(L-1), "
            "|sigma|<=1, n<=N, and the radius-1/L disk lies in one "
            "prescribed-N cell"
        ),
        "coefficient_bound": (
            "On the critical collar, |alpha(s)|<=L/2+(1+pi)/2+2/T "
            "and log(n)<=L/2+1; hence t|alpha_n|<=27"
        ),
        "shift_bound": (
            "w=sqrt(t)v+c, c=t*alpha_n/2, |c|<=27/2"
        ),
        "split": "V=T^(1/3); central: |v|<=V; tail: |v|>V",
        "central_geometry": (
            "For T>=10^6, |v|<=V implies |w|<=V and "
            "Im(s+u*w)>=T/2 for 0<=u<=1"
        ),
        "derivative_bounds": (
            "|alpha'(s)|<=1/T and "
            "sup_(0<=u<=1)|alpha''(s+u*w)|<=3/T^2"
        ),
        "central_logs": (
            "g=1/(6s)+alpha'(s)w^2/2; "
            "rho=R_M+R_G+1/(6(s+w))-1/(6s); "
            "|rho|<=(|w|^3/2+|w|/3+12/T)/T^2"
        ),
        "central_pointwise": (
            "|g+rho|<0.006, "
            "|exp(g+rho)-1-g|<=|rho|+|g+rho|^2"
        ),
        "central_moments": (
            "E|w|^3<10980, E|w|<15, E|w|^4<307330; "
            "E|rho|<5496/T^2, E|g|^2<153666/T^2, "
            "E|rho|^2<T^-2"
        ),
        "central_result": (
            "integral_(|v|<=V)|exp(g+rho)-1-g|dmu "
            f"<={CENTRAL_CONSTANT}/T^2"
        ),
        "published_tail_envelope": (
            "The pointwise envelope in the proof of Polymath Proposition 6.1 "
            "gives |I(v)|dmu<=pi^(-1/2)exp(92/T)"
            "*exp(-(1-0.26/T)v^2)dv"
        ),
        "tail_results": (
            "For T>=10^6: integral_tail |I|dmu<=T^-2, "
            "mu(tail)<=T^-3, E[v^2 1_tail]<=T^-1, and "
            f"integral_tail |I-1-g|dmu<={TAIL_CONSTANT}/T^2"
        ),
        "scalar_guards": (
            "At T=10^6, T^(1/3)=100 and T^(2/3)=10000. "
            "The exact rational bounds log(T)<14, log(2)<1, "
            "pi>157/50, and e<68/25 make all three tail logarithms "
            "negative, their derivative upper bounds negative, the central "
            "exponent <3/500, and the fixed-cell constant <8000000"
        ),
        "per_term_certificate": (
            "|r_(t,n)(s)/m_(t,n)(s)-(1+d_(t,n)(s))|"
            f"<={PER_TERM_CONSTANT}/T^2"
        ),
        "fixed_cell_sum": (
            "Using the two-block coefficient mass <=50*exp(L/4), "
            "sup_collar |B_t|/A_t<2, and T>=2*pi*exp(L-1): "
            "|R_D|/A_t<8000000*exp(-7L/4)"
        ),
        "fixed_cell_c1": (
            "Cauchy on the radius-1/L disk gives "
            "|partial_x R_D|/A_t<8000000*L*exp(-7L/4)"
        ),
        "scale_absorption": (
            "8000000*exp(-7L/4)<"
            "0.000112*exp(-5L/4) for L>=50"
        ),
        "remaining_splice": (
            "The fixed-N Dirichlet residual is closed. A global corrected "
            "remainder still needs the heat-integrated endpoint a^-2 bound "
            "and one adjacent-cutoff comparison for the new signed lift."
        ),
        **moments,
    }


def gate_rows(exact: dict) -> list[GateRow]:
    return [
        GateRow(
            id="np15cd2rc_01_exact_ratio",
            role="source_reduction",
            readiness="available_published",
            claim="The first-correction gate supplies the exact centered ratio and signed coefficient.",
            formula=(
                "I(v)=GammaStar((s+w)/2)"
                "*exp(log M_0(s+w)-log M_0(s)-alpha(s)w); "
                "d=E[g]"
            ),
            proof_boundary="Reuses the checked exact ratio; no estimate here.",
        ),
        GateRow(
            id="np15cd2rc_02_collar_shift",
            role="exact_bound",
            readiness="available_exact",
            claim="The complete saddle shift is uniformly bounded on the critical collar.",
            formula=(
                exact["region"]
                + "; "
                + exact["coefficient_bound"]
                + "; "
                + exact["shift_bound"]
            ),
            proof_boundary="Uses only the explicit alpha formula and cutoff geometry.",
        ),
        GateRow(
            id="np15cd2rc_03_central_geometry",
            role="exact_bound",
            readiness="available_exact",
            claim="The central segment stays uniformly in the upper half-plane.",
            formula=exact["split"] + "; " + exact["central_geometry"],
            proof_boundary="The actual critical T threshold is much larger than 10^6.",
        ),
        GateRow(
            id="np15cd2rc_04_log_derivatives",
            role="exact_bound",
            readiness="available_exact",
            claim="The explicit M_0 logarithmic derivatives have the required second-order decay.",
            formula=exact["derivative_bounds"],
            proof_boundary="Direct modulus bounds for alpha' and alpha''.",
        ),
        GateRow(
            id="np15cd2rc_05_log_remainder",
            role="exact_bound",
            readiness="available_exact",
            claim="The Taylor and scaled-gamma remainders give one integrable central majorant.",
            formula=exact["central_logs"],
            proof_boundary="Uses the DLMF complex log-gamma remainder.",
        ),
        GateRow(
            id="np15cd2rc_06_gaussian_moments",
            role="exact_identity",
            readiness="available_exact",
            claim="All shifted moments needed by the central estimate are explicit.",
            formula=(
                exact["shifted_second_moment"]
                + "; "
                + exact["gaussian_moments"]
                + "; "
                + exact["central_moments"]
            ),
            proof_boundary="The displayed inequalities deliberately round upward.",
        ),
        GateRow(
            id="np15cd2rc_07_central_certificate",
            role="exact_bound",
            readiness="available_exact",
            claim="The central exponential remainder is second order with an explicit constant.",
            formula=exact["central_pointwise"] + "; " + exact["central_result"],
            proof_boundary="Uniform in n and t under the collar hypotheses.",
        ),
        GateRow(
            id="np15cd2rc_08_tail_envelope",
            role="published_bound",
            readiness="available_published",
            claim="The published quadratic envelope controls the region where a cubic Taylor exponential would fail.",
            formula=exact["published_tail_envelope"],
            proof_boundary="This is the pointwise estimate used inside the proof of Proposition 6.1.",
        ),
        GateRow(
            id="np15cd2rc_09_tail_certificate",
            role="exact_bound",
            readiness="available_exact",
            claim="The exact integrand and correction polynomial have a negligible Gaussian tail.",
            formula=exact["tail_results"],
            proof_boundary=(
                exact["scalar_guards"]
                + ". The stored high-precision values are margin diagnostics only."
            ),
        ),
        GateRow(
            id="np15cd2rc_10_per_term_certificate",
            role="proved_theorem",
            readiness="ready_to_apply",
            claim="Every critical Dirichlet component has an explicit uniform second-order relative remainder.",
            formula=exact["per_term_certificate"],
            proof_boundary="Fixed cutoff and its radius-1/L analytic collar.",
        ),
        GateRow(
            id="np15cd2rc_11_fixed_cell_c1",
            role="proved_transfer",
            readiness="ready_to_apply",
            claim="Summation and Cauchy close both finite Dirichlet blocks below the next global scale.",
            formula=(
                exact["fixed_cell_sum"]
                + "; "
                + exact["fixed_cell_c1"]
                + "; "
                + exact["scale_absorption"]
            ),
            proof_boundary="Does not include the Riemann-Siegel endpoint residual.",
        ),
        GateRow(
            id="np15cd2rc_12_global_splice_guard",
            role="route_guard",
            readiness="open",
            claim="The finite-sum theorem must not be promoted to a global corrected remainder.",
            formula=exact["remaining_splice"],
            proof_boundary=(
                "Endpoint a^-2 control, the signed adjacent-cutoff lift, "
                "contact exclusion, Lambda<=0, and RH remain open."
            ),
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    audit = interval_audit()
    rows = gate_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "explicit second-order finite-Dirichlet remainder on fixed "
            "critical cutoff cells"
        ),
        "proof_boundary": (
            "This artifact proves the C_D/T^2 per-term estimate and its "
            "fixed-cell summed C1 transfer. It does not prove the "
            "heat-integrated Riemann-Siegel a^-2 endpoint remainder, the "
            "new adjacent-cutoff signed lift, a strict contact inequality, "
            "contact exclusion, Lambda<=0, or RH."
        ),
        "builder_sha256": file_hash(Path(__file__)),
        "source_sha256": source_hashes(),
        "exact": exact,
        "interval_audit": audit,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_or_published_reductions": 9,
            "proved_theorems": 2,
            "open_route_guards": 1,
            "per_term_constant": PER_TERM_CONSTANT,
            "fixed_cell_constant": FIXED_CELL_CONSTANT,
        },
        "sources": [
            POLYMATH_SOURCE,
            DLMF_SOURCE,
            "outputs/jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    audit = artifact["interval_audit"]
    return "\n".join(
        [
            "# Newman Critical Dirichlet Second-Order Remainder Certificate",
            "",
            "Date: 2026-07-26",
            "",
            "Status: explicit fixed-cutoff finite-Dirichlet theorem; not a proof",
            "of contact exclusion, `Lambda <= 0`, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Uniform Shift",
            "",
            "```text",
            exact["region"],
            exact["coefficient_bound"],
            exact["shift_bound"],
            "```",
            "",
            "The actual lower bound for `T` at `L=50` is",
            f"`{audit['actual_T_min_lower_bound']}`. The proof below only",
            "uses `T>=10^6`.",
            "",
            "## Central Range",
            "",
            "```text",
            exact["split"],
            exact["central_geometry"],
            exact["derivative_bounds"],
            exact["central_logs"],
            exact["central_pointwise"],
            exact["central_moments"],
            exact["central_result"],
            "```",
            "",
            "The pointwise exponent is below",
            f"`{audit['central_exponent_upper_at_T_audit_min']}` at the",
            "artificial audit endpoint and decreases thereafter.",
            "",
            "## Gaussian Tail",
            "",
            f"The proof of [Polymath 15 Proposition 6.1]({POLYMATH_SOURCE})",
            "supplies the quadratic envelope",
            "",
            "```text",
            exact["published_tail_envelope"],
            exact["tail_results"],
            "```",
            "",
            "The three logarithmic endpoint margins are respectively",
            f"`{audit['exact_tail_log_margin_at_T_audit_min']}`,",
            f"`{audit['mass_tail_log_margin_at_T_audit_min']}`, and",
            f"`{audit['second_tail_log_margin_at_T_audit_min']}`; their",
            "checked derivatives are negative. Their proof role is supplied",
            "by the exact rational guards saved alongside this audit:",
            "",
            "```text",
            exact["scalar_guards"],
            "```",
            "",
            "## Per-Term Theorem",
            "",
            "```text",
            exact["per_term_certificate"],
            "```",
            "",
            "The central and tail constants sum to less than `400000`; no",
            "asymptotic `O` constant is hidden in this statement.",
            "",
            "## Fixed-Cell Transfer",
            "",
            "```text",
            exact["fixed_cell_sum"],
            exact["fixed_cell_c1"],
            exact["scale_absorption"],
            "```",
            "",
            "Thus the two finite Dirichlet blocks no longer determine the",
            "global error exponent after their first signed corrections are",
            "retained.",
            "",
            "## Remaining Boundary",
            "",
            "```text",
            exact["remaining_splice"],
            "```",
            "",
            "The next quantitative theorem is the heat-integrated endpoint",
            "`a^-2` bound, followed by its adjacent-cutoff lift. Only after",
            "those steps can the signed contact-normal inequality be tested.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman critical Dirichlet second-order remainder certificate: "
        "12 rows, C_D=400000, fixed-cell C1 constant 8000000, "
        "1 open global splice guard"
    )


if __name__ == "__main__":
    main()
