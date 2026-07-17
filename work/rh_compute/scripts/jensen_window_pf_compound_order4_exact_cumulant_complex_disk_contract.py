#!/usr/bin/env python3
"""Build the cancellation-preserving complex-disk exact-cumulant contract."""

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

import sympy as sp  # noqa: E402

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    sha256,
)
from jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate import (  # noqa: E402
    formal_cumulant_expansion,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md"
)
SOURCE_SECOND_FORMAL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.json"
)
SOURCE_SECOND_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.json"
)
SOURCE_SECOND_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.json"
)
SOURCE_BUDGET = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.json"
)
MAX_EPSILON_ORDER = 10
MAX_POTENTIAL_ORDER = 12
FINITE_Q_FLOOR = 9_000
RAY_Q_FLOOR = 10**35
RAY_START = 20


@dataclass(frozen=True)
class ContractRow:
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


def formal_partition_logarithm() -> tuple[
    tuple[sp.Symbol, ...], sp.Symbol, list[sp.Expr], list[sp.Expr]
]:
    z, y = sp.symbols("z y")
    symbols = tuple(sp.symbols("L_3:13"))
    perturbation = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for order, symbol in zip(range(3, MAX_POTENTIAL_ORDER + 1), symbols):
        perturbation[order - 2] = symbol * y**order / sp.factorial(order)

    exponential = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    exponential[0] = sp.Integer(1)
    for degree in range(1, MAX_EPSILON_ORDER + 1):
        exponential[degree] = sp.expand(
            -sum(
                index * perturbation[index] * exponential[degree - index]
                for index in range(1, degree + 1)
            )
            / degree
        )

    maximum_y_degree = max(sp.Poly(value, y).degree() for value in exponential)
    tilted_moments = [sp.Integer(1)]
    for _degree in range(maximum_y_degree):
        previous = tilted_moments[-1]
        tilted_moments.append(sp.expand(z * previous + sp.diff(previous, z)))

    def tilted_expectation(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient * tilted_moments[monomial[0]]
                for monomial, coefficient in sp.Poly(polynomial, y).terms()
            )
        )

    partition = [tilted_expectation(value) for value in exponential]
    logarithm = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for degree in range(1, MAX_EPSILON_ORDER + 1):
        logarithm[degree] = sp.expand(
            partition[degree]
            - sum(
                index * logarithm[index] * partition[degree - index]
                for index in range(1, degree)
            )
            / degree
        )
    return symbols, z, partition, logarithm


def audit_cumulants(logarithm: list[sp.Expr], z: sp.Symbol) -> dict:
    expected = formal_cumulant_expansion()
    comparisons = 0
    for order in range(2, 9):
        for degree in range(1, MAX_EPSILON_ORDER + 1):
            coefficient = sp.factor(
                sp.diff(logarithm[degree], z, order).subs(z, 0)
            )
            if sp.expand(coefficient - expected[order][degree]) != 0:
                raise RuntimeError(
                    f"partition-log audit failed at kappa_{order}, epsilon^{degree}"
                )
            comparisons += 1
    return {
        "coefficient_comparisons": comparisons,
        "result": "all cumulant coefficients through epsilon^10 agree exactly",
    }


def coefficient_norm(
    expression: sp.Expr,
    symbols: tuple[sp.Symbol, ...],
    z: sp.Symbol,
    caps: dict[int, Fraction],
) -> sp.Rational:
    total = sp.Rational(0)
    for powers, coefficient in sp.Poly(expression, *symbols, z).terms():
        term = abs(sp.Rational(coefficient))
        for order, power in zip(range(3, 13), powers[:-1]):
            cap = caps[order]
            term *= sp.Rational(cap.numerator, cap.denominator) ** power
        total += term
    return sp.factor(total)


def log_majorant_coefficients(partition_norms: list[sp.Rational], order: int) -> list[sp.Rational]:
    coefficients = [sp.Rational(0) for _ in range(order + 1)]
    for degree in range(1, order + 1):
        source = partition_norms[degree] if degree < len(partition_norms) else 0
        coefficients[degree] = sp.factor(
            (
                degree * source
                + sum(
                    index
                    * coefficients[index]
                    * (
                        partition_norms[degree - index]
                        if degree - index < len(partition_norms)
                        else 0
                    )
                    for index in range(1, degree)
                )
            )
            / degree
        )
    return coefficients


def formal_majorants(
    symbols: tuple[sp.Symbol, ...],
    z: sp.Symbol,
    partition: list[sp.Expr],
) -> dict:
    finite_source = load_json(SOURCE_SECOND_FINITE)
    finite_caps = {
        int(order): Fraction(value)
        for order, value in finite_source["normalized_jet_caps"].items()
    }
    expected_caps = {
        3: Fraction(6, 5),
        4: Fraction(3, 2),
        5: Fraction(2),
        6: Fraction(3),
        7: Fraction(9, 2),
        8: Fraction(7),
        9: Fraction(12),
        10: Fraction(21),
        11: Fraction(38),
        12: Fraction(71),
    }
    if finite_caps != expected_caps:
        raise RuntimeError("finite normalized-jet caps changed")
    ray_caps = {order: Fraction(2) for order in range(3, 13)}
    finite_norms = [sp.Rational(0)] + [
        coefficient_norm(partition[degree], symbols, z, finite_caps)
        for degree in range(1, 11)
    ]
    ray_norms = [sp.Rational(0)] + [
        coefficient_norm(partition[degree], symbols, z, ray_caps)
        for degree in range(1, 11)
    ]

    epsilon_cap = sp.Rational(1, 94)
    finite_partition_defect = sum(
        finite_norms[degree] * epsilon_cap**degree
        for degree in range(1, 11)
    )
    if not finite_partition_defect < sp.Rational(1, 100):
        raise RuntimeError("finite formal partition nonvanishing gate failed")

    finite_log = log_majorant_coefficients(finite_norms, 20)
    partial_ratio = sum(
        finite_log[degree] / sp.Integer(94) ** (degree - 6)
        for degree in range(11, 21)
    )
    scaled_partition = sum(
        finite_norms[degree] * sp.Rational(20, 94) ** degree
        for degree in range(1, 11)
    )
    if not scaled_partition < sp.Rational(71, 100):
        raise RuntimeError("finite formal log-tail scaling gate failed")
    tail_ratio = sp.Rational(94**6 * 71, 29 * 20**21)
    finite_log_ratio = sp.factor(partial_ratio + tail_ratio)
    if not finite_log_ratio < sp.Rational(1, 7500):
        raise RuntimeError("finite formal log recomposition gate failed")

    ray_radius = sp.Rational(1, 100)
    ray_partition_at_radius = sum(
        ray_norms[degree] * ray_radius**degree for degree in range(1, 11)
    )
    if not ray_partition_at_radius < sp.Rational(1, 50):
        raise RuntimeError("ray formal partition radius gate failed")
    ray_tail_endpoint_left = (100**11 * 100_000 * RAY_START) ** 2
    ray_tail_endpoint_right = 49**2 * RAY_Q_FLOOR**5
    if ray_tail_endpoint_left >= ray_tail_endpoint_right:
        raise RuntimeError("ray formal log recomposition endpoint failed")
    ray_monotonicity = Fraction(10) - Fraction(1, RAY_START)
    if ray_monotonicity <= 0:
        raise RuntimeError("q^(5/2)/u monotonicity failed")

    rows = {}
    for degree in range(1, 11):
        rows[str(degree)] = {
            "partition_formula": sp.sstr(partition[degree]),
            "partition_terms": len(sp.Poly(partition[degree], *symbols, z).terms()),
            "z_degree": sp.Poly(partition[degree], z).degree(),
            "finite_coefficient_norm": str(finite_norms[degree]),
            "ray_coefficient_norm": str(ray_norms[degree]),
        }
    return {
        "partition_rows": rows,
        "finite": {
            "normalized_jet_caps": {
                str(order): str(cap) for order, cap in finite_caps.items()
            },
            "epsilon_cap": "1/94",
            "partition_defect_bound": str(finite_partition_defect),
            "partition_nonvanishing": "|P^[10](z)-1|<1/100 on |z|<=1",
            "majorant_log_partial_ratio": str(partial_ratio),
            "majorant_log_tail_ratio": str(tail_ratio),
            "majorant_log_total_ratio": str(finite_log_ratio),
            "formal_log_recomposition": (
                "|log P^[10]-S^[10]|<1/(7500*q^3)"
            ),
        },
        "ray": {
            "normalized_jet_cap": "|L_r|<2, r=3,...,12",
            "majorant_radius": "1/100",
            "partition_majorant_at_radius": str(ray_partition_at_radius),
            "partition_nonvanishing": "|P^[10](z)-1|<1/100 on |z|<=1",
            "formal_log_tail_majorant": (
                "100^11/(49*q^(11/2))<1/(100000*u*q^3)"
            ),
            "endpoint_left_squared": str(ray_tail_endpoint_left),
            "endpoint_right_squared": str(ray_tail_endpoint_right),
            "q_power_log_derivative_lower": str(ray_monotonicity),
            "formal_log_recomposition": (
                "|log P^[10]-S^[10]|<1/(100000*u*q^3)"
            ),
        },
    }


def scalar_cauchy_budgets() -> dict:
    finite_log_constant = Fraction(1, 7500) + Fraction(200, 99 * 100_000)
    finite_scaled_cap = 56 * finite_log_constant
    finite_target = Fraction(9, 1000)
    if not finite_scaled_cap < finite_target:
        raise RuntimeError("finite Cauchy cumulant budget failed")

    ray_log_constant = Fraction(1, 100_000) + Fraction(200, 99 * 20_000)
    ray_scaled_cap = 56 * ray_log_constant
    ray_target = Fraction(1, 100)
    if not ray_scaled_cap < ray_target:
        raise RuntimeError("ray Cauchy cumulant budget failed")
    return {
        "cauchy_formula": (
            "sup_|z|<=1 |H_u(z)|<=delta => |H_u^(r)(0)|<=r!*delta"
        ),
        "scaled_derivative_factor": (
            "r!/(r-2)!=r(r-1)<=56, r=2,...,8"
        ),
        "finite": {
            "required_partition_residual": (
                "sup_|z|<=1 |A_u(z)-P_u^[10](z)|<1/(100000*q^3)"
            ),
            "partition_to_log_factor": "200/99",
            "formal_log_constant": "1/7500",
            "total_log_constant": str(finite_log_constant),
            "scaled_cumulant_cap": str(finite_scaled_cap),
            "required_scaled_budget": str(finite_target),
            "margin": str(finite_target - finite_scaled_cap),
        },
        "ray": {
            "required_partition_residual": (
                "sup_|z|<=1 |A_u(z)-P_u^[10](z)|<1/(20000*u*q^3)"
            ),
            "partition_to_log_factor": "200/99",
            "formal_log_constant": "1/100000",
            "total_log_constant": str(ray_log_constant),
            "scaled_cumulant_cap": f"{ray_scaled_cap}/u",
            "required_scaled_budget": "1/(100u)",
            "margin": f"{ray_target-ray_scaled_cap}/u",
        },
    }


def build_artifact() -> dict:
    symbols, z, partition, logarithm = formal_partition_logarithm()
    audit = audit_cumulants(logarithm, z)
    majorants = formal_majorants(symbols, z, partition)
    cauchy = scalar_cauchy_budgets()
    budget = load_json(SOURCE_BUDGET)
    if budget.get("summary", {}).get("global_epsilon_ten_layer_closed") is not True:
        raise RuntimeError("epsilon-ten remainder budget source is not closed")

    logarithm_rows = {
        str(degree): {
            "formula": sp.sstr(logarithm[degree]),
            "terms": len(sp.Poly(logarithm[degree], *symbols, z).terms()),
        }
        for degree in range(1, 11)
    }
    rows = [
        ContractRow(
            id="co4eccdc_01_exact_partition_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The exact cumulant generating function factors into the Gaussian term and one Gaussian-factored partition A_u(z).",
            formula="K_u(z)=z^2/2+log A_u(z)-log A_u(0)",
            proof_boundary="Exact standardized-coordinate identity only.",
        ),
        ContractRow(
            id="co4eccdc_02_epsilon_ten_partition_algebra",
            role="exact_formal_algebra",
            readiness="ready_to_apply",
            claim="Exact graded recurrences give the Gaussian-factored partition P^[10] and logarithm S^[10] through epsilon ten.",
            formula="P^[10]=1+sum_(n=1)^10 epsilon^n Z_n; S^[10]=sum_(n=1)^10 epsilon^n S_n",
            proof_boundary="Exact finite formal algebra only.",
            diagnostics=audit,
        ),
        ContractRow(
            id="co4eccdc_03_formal_partition_nonvanishing",
            role="exact_interval_composition",
            readiness="ready_to_apply",
            claim="Global normalized-jet caps and exact coefficient norms keep the epsilon-ten formal partition uniformly away from zero on the unit disk.",
            formula="|P_u^[10](z)-1|<1/100 for |z|<=1, u>=2",
            proof_boundary="Formal partition only; not the exact density partition.",
            diagnostics={
                "finite": majorants["finite"]["partition_nonvanishing"],
                "ray": majorants["ray"]["partition_nonvanishing"],
            },
        ),
        ContractRow(
            id="co4eccdc_04_finite_formal_log_recomposition",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="An exact positive majorant recurrence controls the formal log recomposition tail on the finite interval.",
            formula=majorants["finite"]["formal_log_recomposition"],
            proof_boundary="Formal logarithm recomposition only.",
            diagnostics=majorants["finite"],
        ),
        ContractRow(
            id="co4eccdc_05_ray_formal_log_recomposition",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="A fixed-radius majorant and q-growth control the formal log recomposition tail on the asymptotic ray.",
            formula=majorants["ray"]["formal_log_recomposition"],
            proof_boundary="Formal logarithm recomposition only.",
            diagnostics=majorants["ray"],
        ),
        ContractRow(
            id="co4eccdc_06_relative_partition_log_lemma",
            role="exact_complex_analytic_lemma",
            readiness="ready_to_apply",
            claim="A small exact-minus-formal partition residual implies nonvanishing and a controlled logarithmic residual on the unit disk.",
            formula="|P-1|<1/100 and |A-P|=rho => |log A-log P|<(200/99)*rho",
            proof_boundary="Exact complex inequality under the displayed residual hypothesis.",
        ),
        ContractRow(
            id="co4eccdc_07_cauchy_cumulant_lemma",
            role="exact_complex_analytic_lemma",
            readiness="ready_to_apply",
            claim="Cauchy's estimate converts one unit-disk logarithmic residual bound into simultaneous cumulant bounds through order eight.",
            formula=cauchy["cauchy_formula"],
            proof_boundary="Exact complex-analytic reduction only.",
        ),
        ContractRow(
            id="co4eccdc_08_finite_partition_target",
            role="exact_theorem_reduction",
            readiness="ready_to_apply",
            claim="The displayed finite exact partition residual is sufficient for every epsilon-ten cumulant budget on 2<=u<=20.",
            formula=cauchy["finite"]["required_partition_residual"],
            proof_boundary="Sufficient target only; the partition residual is not proved here.",
            diagnostics=cauchy["finite"],
        ),
        ContractRow(
            id="co4eccdc_09_ray_partition_target",
            role="exact_theorem_reduction",
            readiness="ready_to_apply",
            claim="The displayed asymptotic exact partition residual is sufficient for every epsilon-ten cumulant budget on u>=20.",
            formula=cauchy["ray"]["required_partition_residual"],
            proof_boundary="Sufficient target only; the partition residual is not proved here.",
            diagnostics=cauchy["ray"],
        ),
        ContractRow(
            id="co4eccdc_10_central_tail_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove the exact partition residual by subtracting P^[10] inside the central integral and bounding the left and right tails separately.",
            formula="central partition residual + exact left tail + exact right tail + formal Gaussian tails",
            proof_boundary="Open exact-density theorem; no cumulant corridor is asserted.",
        ),
    ]
    source_paths = (
        SOURCE_SECOND_FORMAL,
        SOURCE_SECOND_FINITE,
        SOURCE_SECOND_RAY,
        SOURCE_BUDGET,
    )
    return {
        "kind": "jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract",
        "date": "2026-07-13",
        "status": "exact complex-disk reduction with open partition residual",
        "proof_boundary": (
            "This artifact proves the epsilon-ten partition algebra, formal partition "
            "nonvanishing, formal logarithm recomposition bounds, and exact Cauchy "
            "reductions to two explicit unit-disk partition residual targets. It does "
            "not prove those exact central or tail residuals, the exact cumulant "
            "corridors, curvature ray, order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "exact_partition_definition": {
            "mode": "x_u=2*log(u), t_u=V'(x_u), a_u=V''(x_u)",
            "standardized_potential": (
                "W_u(y)=V(x_u+y/sqrt(a_u))-V(x_u)-t_u*y/sqrt(a_u)"
            ),
            "gaussian_factored_partition": (
                "A_u(z)=exp(-z^2/2)/sqrt(2*pi)*integral_R exp(z*y-W_u(y))dy"
            ),
            "cumulant_generator": "K_u(z)=z^2/2+log A_u(z)-log A_u(0)",
        },
        "formal_partition": {
            "grading": "lambda_r=L_r*epsilon^(r-2), epsilon=q^(-1/2)",
            "partition_rows": majorants["partition_rows"],
            "logarithm_rows": logarithm_rows,
            "cumulant_audit": audit,
        },
        "majorants": {
            "finite": majorants["finite"],
            "ray": majorants["ray"],
        },
        "cauchy_budgets": cauchy,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 9,
            "open_analytic_rows": 1,
            "partition_degrees": 10,
            "cumulant_coefficient_audits": 70,
            "finite_partition_target": "1/(100000*q^3)",
            "ray_partition_target": "1/(20000*u*q^3)",
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            path.relative_to(REPO_ROOT).as_posix(): sha256(path)
            for path in source_paths
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.py"
        ),
        "remaining_target": (
            "Prove the unit-disk exact-minus-formal partition residual below "
            "1/(100000*q^3) on 2<=u<=20 and below 1/(20000*u*q^3) on u>=20, "
            "with central and two-tail terms separated."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["cauchy_budgets"]["finite"]
    ray = artifact["cauchy_budgets"]["ray"]
    lines = [
        "# Jensen-Window PF Order-Four Exact Cumulant Complex-Disk Contract",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact complex-disk reduction with open partition residual.",
        "This is not a proof of the exact cumulant corridors, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.py",
        "```",
        "",
        "## Exact Factorization",
        "",
        "At the mode, put",
        "",
        "```text",
        artifact["exact_partition_definition"]["standardized_potential"],
        artifact["exact_partition_definition"]["gaussian_factored_partition"],
        artifact["exact_partition_definition"]["cumulant_generator"],
        "```",
        "",
        "Exact graded Gaussian algebra constructs `P_u^[10](z)` and",
        "`S_u^[10](z)` and audits all 70 cumulant coefficients through epsilon ten.",
        "The global normalized-jet caps give",
        "",
        "```text",
        "|P_u^[10](z)-1|<1/100, |z|<=1, u>=2.",
        "```",
        "",
        "Thus the formal partition has no unit-disk zero. Exact positive-majorant",
        "recurrences prove",
        "",
        "```text",
        artifact["majorants"]["finite"]["formal_log_recomposition"],
        artifact["majorants"]["ray"]["formal_log_recomposition"],
        "```",
        "",
        "## Cauchy Reduction",
        "",
        "The relative partition logarithm lemma and Cauchy's estimate reduce the",
        "simultaneous cumulant problem to the following sufficient targets:",
        "",
        "```text",
        finite["required_partition_residual"] + ", 2<=u<=20,",
        ray["required_partition_residual"] + ", u>=20.",
        "```",
        "",
        f"The resulting finite scaled cumulant cap is `{finite['scaled_cumulant_cap']}`",
        f"against `9/1000`; the remaining margin is `{finite['margin']}`.",
        f"The ray cap is `{ray['scaled_cumulant_cap']}` against `1/(100u)`,",
        f"leaving `{ray['margin']}`.",
        "",
        "## Remaining Boundary",
        "",
        "The targets above are not yet proved. Their advantage is structural: the",
        "epsilon-ten subtraction occurs at partition level before differentiation,",
        "so the order-seven/eight cancellations are preserved automatically. The",
        "next proof must split the exact central residual, exact left tail, exact",
        "right tail, and formal Gaussian tails without returning to independent raw",
        "moment boxes.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
        "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.md",
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
        "certified order-four exact cumulant complex-disk contract: "
        "10 rows, 9 exact rows, 70 cumulant audits, "
        "2 sufficient partition targets, 1 open central-tail theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
