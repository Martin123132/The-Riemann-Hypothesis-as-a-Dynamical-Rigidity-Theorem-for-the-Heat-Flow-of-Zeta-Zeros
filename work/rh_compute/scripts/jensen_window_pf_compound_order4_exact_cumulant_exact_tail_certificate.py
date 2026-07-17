#!/usr/bin/env python3
"""Certify the two exact-density tails in the epsilon-ten disk contract."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402
import sympy as sp  # noqa: E402

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    sha256,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.md"
)
SOURCE_FORMAL_TAIL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.json"
)
SOURCE_GLOBAL_GEOMETRY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json"
)
FINITE_Q_FLOOR = 9_000
SHIFTED_U_FLOOR = Fraction(39, 20)
SHIFTED_Q_FLOOR = 7_200
RAY_START = 20
LOCAL_CURVATURE_RATIO = Fraction(59_319, 100_000)


@dataclass(frozen=True)
class TailRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def potential_curvature_polynomials() -> dict:
    u, q, v, capital_q = sp.symbols("u q v capital_q", nonnegative=True)
    numerator = (
        64 * q**3 * u
        + 16 * q**3
        + 1408 * q**2 * u
        - 84 * q**2
        - 4560 * q * u
        + 120 * q
        + 3600 * u
        - 45
    )
    curvature = u * numerator / (4 * (2 * q - 3) ** 2)

    upper_margin = sp.expand(20 * u * q * (2 * q - 3) ** 2 - numerator)
    upper_shifted = sp.Poly(
        sp.expand(
            upper_margin.subs({u: 2 + v, q: FINITE_Q_FLOOR + capital_q})
        ),
        v,
        capital_q,
    )
    lower_margin = sp.expand(5 * numerator - 78 * u * q * (2 * q - 3) ** 2)
    lower_shifted = sp.Poly(
        sp.expand(
            lower_margin.subs(
                {
                    u: sp.Rational(
                        SHIFTED_U_FLOOR.numerator, SHIFTED_U_FLOOR.denominator
                    )
                    + v,
                    q: SHIFTED_Q_FLOOR + capital_q,
                }
            )
        ),
        v,
        capital_q,
    )
    if not all(coefficient > 0 for coefficient in upper_shifted.coeffs()):
        raise RuntimeError("upper curvature polynomial positivity failed")
    if not all(coefficient > 0 for coefficient in lower_shifted.coeffs()):
        raise RuntimeError("lower curvature polynomial positivity failed")
    return {
        "potential": (
            "V=100*u^2+q-5*u-log(2*q-3)-log(u), q=pi*exp(4*u)"
        ),
        "curvature_formula": f"V''={sp.sstr(curvature)}",
        "central_upper_margin_after_substitution": sp.sstr(upper_shifted.as_expr()),
        "central_upper_margin_terms": len(upper_shifted.terms()),
        "central_upper_minimum_coefficient": str(min(upper_shifted.coeffs())),
        "shifted_lower_margin_after_substitution": sp.sstr(lower_shifted.as_expr()),
        "shifted_lower_margin_terms": len(lower_shifted.terms()),
        "shifted_lower_minimum_coefficient": str(min(lower_shifted.coeffs())),
        "central_bound": "V''(x_u)<=5*u^2*q for u>=2, q>=9000",
        "shifted_bound": (
            "V''(x_v)>=(39/10)*v^2*q_v for v>=39/20, q_v>=7200"
        ),
    }


def adaptive_collar_geometry() -> dict:
    flint.ctx.prec = 256
    q = flint.arb(FINITE_Q_FLOOR)
    s = (32 * q.log()).sqrt()
    y = 1 + s
    drift_ceiling = (flint.arb(39) * q / 10).sqrt() / 10
    drift_margin = drift_ceiling - y
    monotonicity_margin = (flint.arb(39) * q / 10).sqrt() * s - 320
    actual_q_at_two = flint.arb.pi() * flint.arb(8).exp()
    if not bool(
        actual_q_at_two > FINITE_Q_FLOOR
        and y > 13
        and drift_margin > 0
        and monotonicity_margin > 0
    ):
        raise RuntimeError("adaptive exact-tail collar gate failed")

    curvature_ratio = (
        Fraction(39, 50) * Fraction(39, 40) ** 2 * Fraction(4, 5)
    )
    if curvature_ratio != LOCAL_CURVATURE_RATIO or curvature_ratio <= Fraction(1, 2):
        raise RuntimeError("local strong-convexity ratio failed")
    return {
        "window": "Y(q)=1+sqrt(32*log(q))",
        "q_at_u_two_lower": arb_lower_text(actual_q_at_two),
        "q_floor": FINITE_Q_FLOOR,
        "Y_at_q_floor_lower": arb_lower_text(y),
        "drift_ceiling_at_q_floor_lower": arb_lower_text(drift_ceiling),
        "drift_margin_at_q_floor_lower": arb_lower_text(drift_margin),
        "drift_monotonicity_margin_lower": arb_lower_text(monotonicity_margin),
        "drift_inequality": "Y<sqrt((39/10)*q)/10 for q>=9000",
        "standardized_shift": (
            "delta=Y/(2*sqrt(a)); u*delta<1/20 and delta<1/40"
        ),
        "left_u_ratio": "u_-/u=exp(-delta)>39/40",
        "left_q_ratio": "q_-/q=exp(-4*(u-u_-))>4/5",
        "left_floors": "u_->39/20 and q_->7200",
        "curvature_ratio": str(curvature_ratio),
        "strong_convexity": "W_u''(y)>59319/100000>1/2 for |y|<=Y",
    }


def exact_tail_budgets() -> dict:
    flint.ctx.prec = 256
    q = flint.arb(FINITE_Q_FLOOR)
    s = (32 * q.log()).sqrt()
    finite_log_margin = (
        5 * q.log() - flint.arb(500_000).log() - flint.arb(3) / 4 - s / 2
    )
    finite_derivative_margin = 5 - 8 / s

    u = flint.arb(RAY_START)
    ray_q = flint.arb.pi() * (4 * u).exp()
    ray_s = (32 * ray_q.log()).sqrt()
    ray_log_margin = (
        5 * ray_q.log()
        - (300_000 * u).log()
        - flint.arb(3) / 4
        - ray_s / 2
    )
    ray_derivative_margin = 20 - flint.arb(1) / u - 32 / ray_s
    if not bool(
        finite_log_margin > 0
        and finite_derivative_margin > 0
        and ray_log_margin > 0
        and ray_derivative_margin > 0
    ):
        raise RuntimeError("exact-tail scalar budget failed")
    return {
        "unit_disk_prefactor": (
            "|exp(-z^2/2)|/sqrt(2*pi)<=sqrt(e/(2*pi))<1 for |z|<=1"
        ),
        "endpoint_geometry": "W_u(+-Y)>=Y^2/4; outward slope>=Y/2",
        "effective_tail_slope": "Y/2-1>1",
        "common_one_sided_bound": (
            "exact tail <exp(Y-Y^2/4)=exp(3/4+sqrt(32*log(q))/2)*q^-8"
        ),
        "finite": {
            "target": "each exact tail <1/(500000*q^3)",
            "endpoint_log_margin_lower": arb_lower_text(finite_log_margin),
            "log_q_derivative_margin_lower": arb_lower_text(
                finite_derivative_margin
            ),
            "propagation": (
                "5*log(q)-log(500000)-3/4-sqrt(32*log(q))/2 is increasing"
            ),
        },
        "ray": {
            "target": "each exact tail <1/(300000*u*q^3)",
            "q_at_u_twenty_lower": arb_lower_text(ray_q),
            "endpoint_log_margin_lower": arb_lower_text(ray_log_margin),
            "u_derivative_margin_lower": arb_lower_text(ray_derivative_margin),
            "propagation": (
                "5*log(q)-log(300000*u)-3/4-sqrt(32*log(q))/2 is increasing"
            ),
        },
    }


def inherited_global_geometry(source: dict) -> dict:
    matching = [
        row
        for row in source.get("rows", [])
        if row.get("id") == "fsprc_02_left_tail_monotonicity"
    ]
    if len(matching) != 1:
        raise RuntimeError("global left-tail geometry row is missing")
    row = matching[0]
    if row.get("readiness") != "interval_validated":
        raise RuntimeError("global left-tail geometry row is not validated")
    diagnostics = row.get("diagnostics", {})
    if diagnostics.get("positive_curvature_intervals") != 9164:
        raise RuntimeError("global curvature cover changed")
    return {
        "source_row": row["id"],
        "curvature_cover": "V''>0 for u>=1/100",
        "tiny_u_slope": "V'<0 for 0<u<=1/100",
        "positive_curvature_intervals": diagnostics["positive_curvature_intervals"],
        "tiny_u_slope_upper": diagnostics["tiny_u_slope_upper"],
        "left_reference_slope_lower": diagnostics["left_reference_slope_lower"],
        "consequence": (
            "both outward slopes stay above their values at y=+-Y on the full tails"
        ),
    }


def build_artifact() -> dict:
    formal_tail = load_json(SOURCE_FORMAL_TAIL)
    if formal_tail.get("summary", {}).get("formal_tails_closed") != 2:
        raise RuntimeError("formal-tail source is not closed")
    global_source = load_json(SOURCE_GLOBAL_GEOMETRY)
    global_geometry = inherited_global_geometry(global_source)
    curvature = potential_curvature_polynomials()
    collar = adaptive_collar_geometry()
    budgets = exact_tail_budgets()
    rows = [
        TailRow(
            id="co4ecetc_01_curvature_identity",
            role="exact_symbolic_identity",
            readiness="ready_to_apply",
            claim="The first-summand potential has an explicit rational-exponential curvature formula.",
            formula=curvature["curvature_formula"],
            proof_boundary="Exact potential algebra only.",
        ),
        TailRow(
            id="co4ecetc_02_curvature_bounds",
            role="exact_polynomial_inequality",
            readiness="ready_to_apply",
            claim="Positive-coefficient substitutions prove the required central upper and shifted lower curvature bounds.",
            formula="(39/10)*v^2*q_v<=V''(x_v); V''(x_u)<=5*u^2*q",
            proof_boundary="Curvature bounds on the displayed u and q domains only.",
            diagnostics=curvature,
        ),
        TailRow(
            id="co4ecetc_03_adaptive_shift",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The adaptive standardized collar changes the mode coordinate and q by controlled fixed ratios.",
            formula="u_-/u>39/40; q_-/q>4/5",
            proof_boundary="Adaptive collar geometry only.",
            diagnostics=collar,
        ),
        TailRow(
            id="co4ecetc_04_local_strong_convexity",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The standardized exact potential is uniformly strongly convex throughout the adaptive collar.",
            formula=collar["strong_convexity"],
            proof_boundary="Strong convexity only for |y|<=Y.",
        ),
        TailRow(
            id="co4ecetc_05_global_outward_monotonicity",
            role="exact_interval_analytic_lemma",
            readiness="ready_to_apply",
            claim="The inherited global potential geometry propagates both endpoint slopes through the complete exact tails.",
            formula="V''>0 for u>=1/100; V'<0 for 0<u<=1/100",
            proof_boundary="Potential monotonicity only; no partition estimate by itself.",
            diagnostics=global_geometry,
        ),
        TailRow(
            id="co4ecetc_06_common_exact_tail_bound",
            role="exact_complex_analytic_theorem",
            readiness="ready_to_apply",
            claim="Strong convexity, outward slopes, and the unit-disk tilt give one common bound for either exact-density tail.",
            formula=budgets["common_one_sided_bound"],
            proof_boundary="Each exact-density tail separately on |z|<=1.",
            diagnostics={
                "prefactor": budgets["unit_disk_prefactor"],
                "endpoint": budgets["endpoint_geometry"],
                "effective_slope": budgets["effective_tail_slope"],
            },
        ),
        TailRow(
            id="co4ecetc_07_finite_exact_tails",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Each exact-density tail fits one fifth of the finite partition residual budget.",
            formula=budgets["finite"]["target"],
            proof_boundary="Exact tails on 2<=u<=20 only.",
            diagnostics=budgets["finite"],
        ),
        TailRow(
            id="co4ecetc_08_ray_exact_tails",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Each exact-density tail fits its assigned asymptotic partition residual budget.",
            formula=budgets["ray"]["target"],
            proof_boundary="Exact tails on u>=20 only.",
            diagnostics=budgets["ray"],
        ),
        TailRow(
            id="co4ecetc_09_central_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Only the central exact-minus-epsilon-ten formal partition residual remains open.",
            formula=(
                "finite <1/(500000*q^3); ray <1/(300000*u*q^3)"
            ),
            proof_boundary="Open central exact-density theorem; no cumulant corridor is asserted.",
        ),
    ]
    source_paths = (SOURCE_FORMAL_TAIL, SOURCE_GLOBAL_GEOMETRY)
    return {
        "kind": "jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate",
        "date": "2026-07-13",
        "status": "exact-density two-tail theorem with one open central residual",
        "proof_boundary": (
            "This artifact proves both exact-density tails fit the epsilon-ten "
            "complex-disk partition budgets. Together with the inherited formal-tail "
            "certificate, only the central exact-minus-formal residual remains. It "
            "does not prove that residual, the exact cumulant corridors, curvature "
            "ray, order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "curvature_certificate": curvature,
        "adaptive_collar": collar,
        "global_tail_geometry": global_geometry,
        "tail_budgets": budgets,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 8,
            "open_analytic_rows": 1,
            "positive_coefficient_polynomials": 2,
            "local_curvature_ratio": str(LOCAL_CURVATURE_RATIO),
            "formal_tails_inherited": 2,
            "exact_tails_closed": 2,
            "open_exact_components": 1,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            path.relative_to(REPO_ROOT).as_posix(): sha256(path)
            for path in source_paths
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.py"
        ),
        "remaining_target": (
            "On |y|<=Y(q), prove the exact Gaussian-factored density differs from "
            "its epsilon-ten formal density by less than 1/(500000*q^3) on "
            "2<=u<=20 and 1/(300000*u*q^3) on u>=20, uniformly for |z|<=1."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Exact Cumulant Exact-Tail Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact-density two-tail theorem with one open central residual.",
        "This is not a proof of the exact cumulant corridors, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.py",
        "```",
        "",
        "## Curvature Geometry",
        "",
        "The exact curvature formula and two positive-coefficient substitutions prove",
        "",
        "```text",
        "V''(x_u)<=5*u^2*q                  (u>=2, q>=9000),",
        "V''(x_v)>=(39/10)*v^2*q_v          (v>=39/20, q_v>=7200).",
        "```",
        "",
        "For `Y(q)=1+sqrt(32*log(q))`, the adaptive shift obeys",
        "",
        "```text",
        "u_-/u>39/40,",
        "q_-/q>4/5,",
        "W_u''(y)>59319/100000>1/2          (|y|<=Y).",
        "```",
        "",
        "Hence `W_u(+-Y)>=Y^2/4`, with outward slopes at least `Y/2`.",
        "The inherited global potential geometry propagates those slopes through",
        "both complete tails.",
        "",
        "## Exact Tails",
        "",
        "Uniformly for `|z|<=1`, either exact-density tail satisfies",
        "",
        "```text",
        artifact["tail_budgets"]["common_one_sided_bound"],
        artifact["tail_budgets"]["finite"]["target"],
        artifact["tail_budgets"]["ray"]["target"],
        "```",
        "",
        "Together with the formal-tail theorem, all four tails in the complex-disk",
        "partition decomposition are closed.",
        "",
        "## Remaining Boundary",
        "",
        "Only the central exact-minus-formal residual remains:",
        "",
        "```text",
        "finite: <1/(500000*q^3),             2<=u<=20,",
        "ray:    <1/(300000*u*q^3),           u>=20.",
        "```",
        "",
        "This central residual is an open target. No exact cumulant corridor is",
        "promoted by the two-tail theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.md",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
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
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four exact cumulant exact tails: "
        "9 rows, 8 exact rows, 2 positive-coefficient polynomials, "
        "2 closed exact tails, 1 open central residual"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
