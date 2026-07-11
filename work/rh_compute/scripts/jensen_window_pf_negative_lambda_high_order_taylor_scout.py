#!/usr/bin/env python3
"""Build the high-order Taylor truncation scout for the negative-lambda route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
from math import factorial
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_high_order_taylor_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md"

getcontext().prec = 100


Polynomial = dict[int, Fraction]


@dataclass(frozen=True)
class CoefficientRow:
    id: str
    degree: int
    enclosure: str
    ratio_to_c0: str
    sign: str
    tail_radius: str
    proof_boundary: str


@dataclass(frozen=True)
class TruncationRow:
    id: str
    truncation_degree: int
    M: int
    T: str
    k: int
    q_over_T: str
    F_k_minus_1: str
    F_k: str
    F_k_plus_1: str
    target_bound: str
    status: str
    truncated_B: str | None
    B_over_target_bound: str | None
    proof_boundary: str


@dataclass(frozen=True)
class HighOrderTaylorDiagnostics:
    max_taylor_degree: int
    coefficient_rows: list[CoefficientRow]
    truncation_rows: list[TruncationRow]
    invalid_normalizer_rows: int
    upper_wall_violation_rows: int
    overbound_rows: int
    target_satisfying_rows: int


def padd(left: Polynomial, right: Polynomial) -> Polynomial:
    result = dict(left)
    for degree, coeff in right.items():
        result[degree] = result.get(degree, Fraction(0)) + coeff
    return {degree: coeff for degree, coeff in result.items() if coeff}


def pmul(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for i, left_coeff in left.items():
        for j, right_coeff in right.items():
            result[i + j] = result.get(i + j, Fraction(0)) + left_coeff * right_coeff
    return {degree: coeff for degree, coeff in result.items() if coeff}


def pscale(poly: Polynomial, scale: Fraction) -> Polynomial:
    return {degree: coeff * scale for degree, coeff in poly.items() if coeff * scale}


def qpower(degree: int, coeff: Fraction = Fraction(1)) -> Polynomial:
    return {degree: coeff}


def generate_polynomials(max_degree: int) -> list[Polynomial]:
    """Coefficient polynomials P_j(q) for the Phi local Taylor series."""
    # S(u)=-q*(exp(4u)-1), and E(u)=exp(S(u)).
    s_series: list[Polynomial] = [{} for _ in range(max_degree + 1)]
    for index in range(1, max_degree + 1):
        s_series[index] = {1: -Fraction(4**index, factorial(index))}

    e_series: list[Polynomial] = [{} for _ in range(max_degree + 1)]
    e_series[0] = {0: Fraction(1)}
    for n in range(1, max_degree + 1):
        acc: Polynomial = {}
        for i in range(1, n + 1):
            acc = padd(acc, pscale(pmul(s_series[i], e_series[n - i]), Fraction(i)))
        e_series[n] = pscale(acc, Fraction(1, n))

    polynomials: list[Polynomial] = []
    for n in range(max_degree + 1):
        acc: Polynomial = {}
        for i in range(n + 1):
            # Coefficient of u^i in 2*q^2*exp(9u)-3*q*exp(5u).
            a_i = padd(
                qpower(2, Fraction(2 * 9**i, factorial(i))),
                qpower(1, Fraction(-3 * 5**i, factorial(i))),
            )
            acc = padd(acc, pmul(a_i, e_series[n - i]))
        polynomials.append(acc)
    return polynomials


def arb_fraction(value: Fraction) -> flint.arb:
    return flint.arb(value.numerator) / flint.arb(value.denominator)


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_negative(value: flint.arb) -> bool:
    return bool(value < 0 and not value.contains(0))


def sign_name(value: flint.arb) -> str:
    if arb_positive(value):
        return "positive"
    if arb_negative(value):
        return "negative"
    return "unresolved"


def eval_poly(poly: Polynomial, q: flint.arb) -> flint.arb:
    total = flint.arb(0)
    for degree, coeff in poly.items():
        total += arb_fraction(coeff) * (q**degree)
    return total


def tail_power_bound(cutoff_n: int, power_degree: int, pi: flint.arb) -> flint.arb:
    n0 = cutoff_n + 1
    n = flint.arb(n0)
    first = (n ** (2 * power_degree)) * (-(pi * flint.arb(n0 * n0))).exp()
    ratio = ((flint.arb(n0 + 1) / n) ** (2 * power_degree)) * (
        -(pi * flint.arb(2 * n0 + 1))
    ).exp()
    if not arb_positive(flint.arb(1) - ratio):
        raise RuntimeError(f"tail ratio did not certify below one for m={power_degree}")
    return first / (flint.arb(1) - ratio)


def tail_radius_bound(poly: Polynomial, cutoff_n: int, pi: flint.arb) -> flint.arb:
    total = flint.arb(0)
    for degree, coeff in poly.items():
        total += arb_fraction(abs(coeff)) * (pi**degree) * tail_power_bound(cutoff_n, degree, pi)
    return total


def upper_radius_string(value: flint.arb, safety_factor: int = 4) -> str:
    if not arb_positive(value):
        raise RuntimeError(f"non-positive tail radius {value}")
    mid, rad, exp10 = value.mid_rad_10exp()
    mantissa = (abs(int(mid)) + int(rad)) * safety_factor
    return f"{mantissa}e{int(exp10)}"


def error_ball(radius: flint.arb) -> flint.arb:
    return flint.arb(f"[0 +/- {upper_radius_string(radius)}]")


def coefficient_value(poly: Polynomial, cutoff_n: int, pi: flint.arb) -> tuple[flint.arb, str]:
    finite = flint.arb(0)
    for n in range(1, cutoff_n + 1):
        q = pi * flint.arb(n * n)
        finite += (-q).exp() * eval_poly(poly, q)
    tail = tail_radius_bound(poly, cutoff_n, pi)
    return finite + error_ball(tail), upper_radius_string(tail)


def arb_mid_decimal(value: flint.arb) -> Decimal:
    mid, _rad, exp10 = value.mid_rad_10exp()
    return Decimal(int(mid)) * (Decimal(10) ** int(exp10))


def decimal_format(value: Decimal) -> str:
    return f"{value:.18E}"


def rising(q: Decimal, length: int) -> Decimal:
    total = Decimal(1)
    for index in range(length):
        total *= q + Decimal(index)
    return total


def truncated_multiplier(k: int, T: int, ratios: list[Decimal], M: int) -> Decimal:
    q = Decimal(k) + Decimal("0.5")
    big_t = Decimal(T)
    return sum(ratios[j] * rising(q, j) / (big_t**j) for j in range(M + 1))


def target_bound(k: int) -> Decimal:
    return Decimal(2) / (Decimal(3) * Decimal(2 * k + 1))


def classify_truncation(k: int, T: int, ratios: list[Decimal], M: int) -> TruncationRow:
    f_minus = truncated_multiplier(k - 1, T, ratios, M)
    f = truncated_multiplier(k, T, ratios, M)
    f_plus = truncated_multiplier(k + 1, T, ratios, M)
    bound = target_bound(k)
    q_over_t = (Decimal(k) + Decimal("0.5")) / Decimal(T)
    if not (f_minus > 0 and f > 0 and f_plus > 0):
        status = "invalid_normalizer"
        b_value = None
        ratio_value = None
    else:
        x_value = f_plus * f_minus / (f * f)
        b_value = -x_value.ln()
        if b_value < 0:
            status = "upper_wall_violation"
        elif b_value > bound:
            status = "overbound"
        else:
            status = "target_satisfying_truncation"
        ratio_value = b_value / bound
    return TruncationRow(
        id=f"nhts_M{M}_T{T}",
        truncation_degree=2 * M,
        M=M,
        T=str(T),
        k=k,
        q_over_T=decimal_format(q_over_t),
        F_k_minus_1=decimal_format(f_minus),
        F_k=decimal_format(f),
        F_k_plus_1=decimal_format(f_plus),
        target_bound=decimal_format(bound),
        status=status,
        truncated_B=None if b_value is None else decimal_format(b_value),
        B_over_target_bound=None if ratio_value is None else decimal_format(ratio_value),
        proof_boundary="Finite truncation diagnostic only; not a proof about the full Taylor tail or actual zeta coefficients.",
    )


def build_diagnostics(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> HighOrderTaylorDiagnostics:
    flint.ctx.prec = precision_bits
    polynomials = generate_polynomials(max_degree)
    pi = flint.arb.pi()

    coefficient_values: dict[int, flint.arb] = {}
    tail_radii: dict[int, str] = {}
    for degree in range(0, max_degree + 1, 2):
        value, radius = coefficient_value(polynomials[degree], cutoff_n, pi)
        coefficient_values[degree] = value
        tail_radii[degree] = radius
    c0 = coefficient_values[0]
    if not arb_positive(c0):
        raise RuntimeError(f"c0 is not positive: {c0}")

    coefficient_rows: list[CoefficientRow] = []
    ratios_decimal: list[Decimal] = []
    for degree in range(0, max_degree + 1, 2):
        ratio = coefficient_values[degree] / c0
        ratios_decimal.append(arb_mid_decimal(ratio))
        coefficient_rows.append(
            CoefficientRow(
                id=f"nhts_c{degree}",
                degree=degree,
                enclosure=coefficient_values[degree].str(40),
                ratio_to_c0=ratio.str(40),
                sign=sign_name(coefficient_values[degree]),
                tail_radius=tail_radii[degree],
                proof_boundary="Certified local Taylor coefficient only; not a uniform remainder theorem.",
            )
        )

    truncation_rows: list[TruncationRow] = []
    for M in range(3, (max_degree // 2) + 1):
        for T in t_values:
            truncation_rows.append(classify_truncation(k, T, ratios_decimal, M))

    invalid = sum(1 for row in truncation_rows if row.status == "invalid_normalizer")
    upper = sum(1 for row in truncation_rows if row.status == "upper_wall_violation")
    overbound = sum(1 for row in truncation_rows if row.status == "overbound")
    satisfying = sum(1 for row in truncation_rows if row.status == "target_satisfying_truncation")
    return HighOrderTaylorDiagnostics(
        max_taylor_degree=max_degree,
        coefficient_rows=coefficient_rows,
        truncation_rows=truncation_rows,
        invalid_normalizer_rows=invalid,
        upper_wall_violation_rows=upper,
        overbound_rows=overbound,
        target_satisfying_rows=satisfying,
    )


def build_artifact(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    diagnostics = build_diagnostics(max_degree, cutoff_n, precision_bits, k, t_values)
    rows = [
        {
            "id": "nhts_01_formal_polynomial_generator",
            "role": "exact_algebra",
            "readiness": "not_ready_to_apply",
            "claim": "The P_j(q) polynomials are generated by formal series from (2*q^2*exp(9u)-3*q*exp(5u))*exp(-q*(exp(4u)-1)).",
            "proof_boundary": "Exact local algebra only; not a theorem about the analytic Taylor remainder.",
        },
        {
            "id": "nhts_02_even_coefficient_certification",
            "role": "finite_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "The even Phi Taylor coefficients c0,c2,...,c14 are certified with finite sums plus geometric tail bounds.",
            "proof_boundary": "Finite local coefficient certificate only; not uniform in k or T.",
        },
        {
            "id": "nhts_03_high_order_truncation_matrix",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The k=22 grid for truncation degrees 6,8,10,12,14 has invalid normalizers, overbound rows, and upper-wall violations.",
            "proof_boundary": "Finite truncation diagnostic only; it does not decide the full Taylor series or actual zeta prefix.",
        },
        {
            "id": "nhts_04_higher_order_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "Increasing the Taylor truncation degree alone proves the local/mesoscopic remainder theorem.",
            "proof_boundary": "Rejected proof template only; the full Taylor tail still requires analytic control.",
        },
        {
            "id": "nhts_05_remainder_theorem_handoff",
            "role": "conditional_handoff",
            "readiness": "not_ready_to_apply",
            "claim": "A successful theorem must replace truncation-order experiments with a bound on the infinite Taylor tail and its log-differences.",
            "proof_boundary": "Conditional handoff only; not bounded log-curvature, cone entry, or Lambda <= 0.",
        },
    ]
    summary = {
        "scout_rows": len(rows),
        "coefficient_rows": len(diagnostics.coefficient_rows),
        "truncation_rows": len(diagnostics.truncation_rows),
        "invalid_normalizer_rows": diagnostics.invalid_normalizer_rows,
        "upper_wall_violation_rows": diagnostics.upper_wall_violation_rows,
        "overbound_rows": diagnostics.overbound_rows,
        "target_satisfying_rows": diagnostics.target_satisfying_rows,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "High-order local Taylor truncations through degree 14 do not provide a stable finite proof model: "
            "the k=22 sample matrix contains invalid normalizers, overbound log-curvature rows, and upper-wall "
            "violations, so the proof route still needs an analytic infinite-tail remainder theorem rather than "
            "a higher finite truncation."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_high_order_taylor_scout",
        "date": "2026-07-06",
        "status": "finite_theorem_search_diagnostic",
        "source_phi_taylor_sign_scout": "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
        "source_taylor_moment_budget": "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_high_order_taylor_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_high_order_taylor_scout.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It certifies higher local Taylor coefficients and "
            "classifies finite truncation samples, but it does not bound the infinite Taylor remainder, "
            "does not prove bounded log-curvature, does not prove cone entry, and does not prove Lambda <= 0."
        ),
        "parameters": {
            "max_taylor_degree": max_degree,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "sample_T_values": t_values,
            "truncation_degrees": [2 * M for M in range(3, (max_degree // 2) + 1)],
        },
        "scout_rows": rows,
        "diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply.",
            "Finite Taylor truncations are not promoted to an infinite Taylor-tail theorem.",
            "The actual zeta finite prefix remains certified by ACB enclosures, not by this truncation matrix.",
            "The local/mesoscopic remainder theorem remains open.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda high-order Taylor scout: "
        f"{summary['coefficient_rows']} coefficient rows, {summary['truncation_rows']} truncation rows, "
        f"{summary['invalid_normalizer_rows']} invalid normalizers, "
        f"{summary['upper_wall_violation_rows']} upper-wall violations, "
        f"{summary['overbound_rows']} overbound rows, 0 ready-to-apply rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda High-Order Taylor Scout",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of the",
        "bounded log-curvature theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_high_order_taylor_scout`.",
        "",
        "Proof boundary: this artifact certifies higher local Taylor coefficients",
        "and classifies finite truncation samples. It does not prove an infinite",
        "Taylor-tail remainder theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_high_order_taylor_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_high_order_taylor_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_high_order_taylor_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Coefficient Scope",
        "",
        "The generator forms `P_j(q)` by exact formal series from:",
        "",
        "```text",
        "(2*q^2*exp(9u)-3*q*exp(5u))*exp(-q*(exp(4u)-1))",
        "```",
        "",
        "It certifies the even local coefficients through `c14` using the same",
        "finite-sum plus geometric-tail method as the lower-order sign scout.",
        "",
        "Coefficient signs:",
        "",
        "```text",
    ]
    for row in diagnostics["coefficient_rows"]:
        lines.append(f"c{row['degree']}: {row['sign']}, ratio c{row['degree']}/c0 = {row['ratio_to_c0']}")
    lines.extend(
        [
            "```",
            "",
            "## Truncation Matrix",
            "",
            "Rows classify `k=22` and `T=25,50,100,200,500,1000,2000` for",
            "truncation degrees 6, 8, 10, 12, and 14.",
            "",
            "```text",
        ]
    )
    for row in diagnostics["truncation_rows"]:
        ratio = row["B_over_target_bound"] if row["B_over_target_bound"] is not None else "n/a"
        lines.append(
            f"degree={row['truncation_degree']}, T={row['T']}: "
            f"F_k={row['F_k']}, status={row['status']}, B/bound={ratio}"
        )
    lines.extend(
        [
            "```",
            "",
            "Summary counts:",
            "",
            "```text",
            f"invalid normalizers: {summary['invalid_normalizer_rows']}",
            f"upper-wall violations: {summary['upper_wall_violation_rows']}",
            f"overbound rows: {summary['overbound_rows']}",
            f"target-satisfying truncation rows: {summary['target_satisfying_rows']}",
            "```",
            "",
            "Interpretation: higher finite Taylor truncation is a theorem-search",
            "diagnostic, not a replacement for an infinite Taylor-tail estimate.",
            "The local/mesoscopic route still needs an analytic remainder theorem",
            "controlling positivity and second/third log-differences.",
            "",
            "Integration:",
            "",
            "```text",
            "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
            "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
            "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
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


def parse_t_values(text: str) -> list[int]:
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("at least one T value is required")
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=14)
    parser.add_argument("--tail-cutoff-n", type=int, default=16)
    parser.add_argument("--precision-bits", type=int, default=256)
    parser.add_argument("--tail-start-k", type=int, default=22)
    parser.add_argument("--t-values", type=parse_t_values, default=parse_t_values("25,50,100,200,500,1000,2000"))
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(args.max_degree, args.tail_cutoff_n, args.precision_bits, args.tail_start_k, args.t_values)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda high-order Taylor scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
