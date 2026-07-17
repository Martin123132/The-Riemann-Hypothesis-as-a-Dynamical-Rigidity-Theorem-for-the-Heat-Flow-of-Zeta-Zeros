#!/usr/bin/env python3
"""Prove the universal first term of every fixed graded Hankel determinant."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import itertools
import json
from math import factorial
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md"
)
HEAT_TILT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.json"
)
HIGHER_THETA_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_superpolynomial_first_summand_dominance.json"
)
XI_RATIO_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.json"
)
MAX_RECORDED_ORDER = 12
MAX_STRESS_ORDER = 8
MAX_COEFFICIENT_DEGREE = 12


@dataclass(frozen=True)
class LemmaRow:
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


def validate_sources() -> dict:
    heat_tilt = load_json(HEAT_TILT_SOURCE)
    higher_theta = load_json(HIGHER_THETA_SOURCE)
    xi_ratio = load_json(XI_RATIO_SOURCE)
    if heat_tilt.get("summary", {}).get("uniform_heat_tilt_theorems") != 1:
        raise RuntimeError("uniform heat-tilt source is not closed")
    if higher_theta.get("summary", {}).get("open_analytic_rows") != 0:
        raise RuntimeError("uniform higher-theta source is not closed")
    if xi_ratio.get("eventual_theorem", {}).get("threshold_effective_here") is not False:
        raise RuntimeError("published Xi ratio source contract changed")
    return {
        "heat_tilt_status": heat_tilt.get("status"),
        "higher_theta_status": higher_theta.get("status"),
        "xi_ratio_source": xi_ratio.get("published_ratio_input", {}).get("source"),
    }


def signed_leading_constant(order: int) -> int:
    degree = order * (order - 1) // 2
    orientation = (-1) ** degree
    return orientation * 2**degree * prod_factorials(order)


def prod_factorials(order: int) -> int:
    value = 1
    for index in range(1, order):
        value *= factorial(index)
    return value


def order_rows() -> list[dict]:
    rows = []
    for order in range(1, MAX_RECORDED_ORDER + 1):
        degree = order * (order - 1) // 2
        raw = signed_leading_constant(order)
        rows.append(
            {
                "order": order,
                "leading_degree": degree,
                "epsilon": (-1) ** degree,
                "positive_constant": abs(raw),
                "raw_leading_constant": raw,
                "raw_first_term": f"{raw}*G2**{degree}*h**{degree}",
                "signed_first_term": f"{abs(raw)}*G2**{degree}*h**{degree}",
            }
        )
    return rows


def coefficient_valuation_table() -> list[dict]:
    rows = []
    for p in range(MAX_COEFFICIENT_DEGREE + 1):
        for q in range(MAX_COEFFICIENT_DEGREE + 1):
            if p == q == 0:
                classification = "constant"
                lower_bound: int | str = 0
                diagonal = "1"
            elif p == 0 or q == 0:
                classification = "zero_pure_axis"
                lower_bound = "infinity"
                diagonal = None
            elif p == q:
                classification = "diagonal"
                lower_bound = p
                diagonal = f"(2*G2)**{p}/{factorial(p)}"
            else:
                classification = "off_diagonal"
                lower_bound = max(p, q)
                diagonal = None
            rows.append(
                {
                    "p": p,
                    "q": q,
                    "classification": classification,
                    "h_valuation_lower_bound": lower_bound,
                    "diagonal_leading_coefficient": diagonal,
                }
            )
    return rows


def permutation_penalty_stress() -> dict:
    rows = []
    total = 0
    for order in range(1, MAX_STRESS_ORDER + 1):
        degree = order * (order - 1) // 2
        zero_penalties = 0
        minimum_positive = None
        cases = 0
        for permutation in itertools.permutations(range(order)):
            cases += 1
            penalty = sum(
                max(index, permutation[index]) for index in range(order)
            ) - degree
            absolute_formula = sum(
                abs(index - permutation[index]) for index in range(order)
            ) // 2
            if penalty != absolute_formula:
                raise RuntimeError("canonical permutation penalty identity failed")
            if penalty == 0:
                zero_penalties += 1
            elif minimum_positive is None or penalty < minimum_positive:
                minimum_positive = penalty
        if zero_penalties != 1:
            raise RuntimeError("identity is not the unique zero-penalty permutation")
        if order > 1 and minimum_positive != 1:
            raise RuntimeError("unexpected minimum positive permutation penalty")
        rows.append(
            {
                "order": order,
                "permutations": cases,
                "zero_penalties": zero_penalties,
                "minimum_positive_penalty": minimum_positive,
            }
        )
        total += cases
    return {"rows": rows, "permutations_checked": total}


def proof_rows(stress: dict, constants: list[dict]) -> list[LemmaRow]:
    order7 = constants[6]
    return [
        LemmaRow(
            "gkvaol_01_graded_kernel",
            "exact_definition",
            "ready_to_apply",
            "The normalized fixed-order Hankel determinant has one graded translation-kernel form.",
            "K_(i,j)=exp(-sum_(r>=2)G_r*h^(r-1)*(L-i-j)^r), L=2*(m-1)",
            "Formal coefficient identity at one fixed order m.",
        ),
        LemmaRow(
            "gkvaol_02_row_column_gauge",
            "exact_algebraic_lemma",
            "ready_to_apply",
            "Pure powers of z and y factor into row and column gauges with constant term one.",
            "(z-y)^r=z^r+(-y)^r+sum_(a,b>=1,a+b=r) binom(r,a)z^a(-y)^b",
            "Gauge factors cannot change the first nonzero determinant coefficient.",
        ),
        LemmaRow(
            "gkvaol_03_mixed_valuation",
            "exact_combinatorial_lemma",
            "ready_to_apply",
            "Every mixed coefficient c_(p,q) has h-valuation at least max(p,q).",
            "ell mixed factors give h-degree p+q-ell>=max(p,q) because ell<=min(p,q)",
            "Coefficientwise formal-series statement; no asymptotic sign is used.",
        ),
        LemmaRow(
            "gkvaol_04_diagonal_leader",
            "exact_combinatorial_lemma",
            "ready_to_apply",
            "At the valuation floor of c_(k,k), only k quadratic mixed interactions occur.",
            "[h^k]c_(k,k)=(2*G_2)^k/k!",
            "Higher G_r cannot enter this coefficient.",
        ),
        LemmaRow(
            "gkvaol_05_cauchy_binet_floor",
            "exact_determinant_lemma",
            "ready_to_apply",
            "Cauchy-Binet forces every m-by-m determinant term to h-degree at least D=binom(m,2).",
            "sum_i max(p_i,q_(pi(i)))>=max(sum_i p_i,sum_i q_i)>=D",
            "Strictly increasing nonnegative exponent tuples only.",
        ),
        LemmaRow(
            "gkvaol_06_unique_leading_block",
            "exact_determinant_lemma",
            "ready_to_apply",
            "Equality at degree D uses P=Q=(0,...,m-1) and the identity coefficient permutation only.",
            "sum_i max(i,pi(i))=D+(1/2)*sum_i|i-pi(i)|",
            "The finite permutation stress is diagnostic support for the exact identity.",
            stress,
        ),
        LemmaRow(
            "gkvaol_07_vandermonde_constant",
            "exact_all_fixed_order_theorem",
            "ready_to_apply",
            "The universal raw leading coefficient has the signed Vandermonde orientation.",
            "[h^D]det K=epsilon_m*2^D*prod_(j=1)^(m-1)j!*G_2^D",
            "Exact for every fixed integer m>=1.",
        ),
        LemmaRow(
            "gkvaol_08_signed_first_term",
            "exact_all_fixed_order_theorem",
            "ready_to_apply",
            "The signed coordinate Q_(m,n)=epsilon_m H_(m,n) always has a positive universal first term.",
            "Q_m=positive_scale*(2^D*prod_(j=1)^(m-1)j!*G_2^D*h^D+o(h^D))",
            "Fixed m; the asymptotic threshold may depend on m.",
        ),
        LemmaRow(
            "gkvaol_09_all_fixed_heat_differences",
            "published_theorem_composition",
            "ready_to_apply",
            "The compact heat tilt obeys the required local difference estimate at every fixed order.",
            "Delta^s log R_T^(1)(k)=O_s(log(k)/k^s) for every fixed s>=2, uniformly for 0<=T<=100",
            "For fixed s choose the published saddle truncation R>s; no uniformity in s is asserted.",
        ),
        LemmaRow(
            "gkvaol_10_uniform_eventual_tail",
            "published_theorem_composition",
            "ready_to_apply",
            "For every fixed order there is a compact-heat uniform eventual positive signed tail.",
            "for every fixed m, exists N_m: Q_(m,n)(lambda)>0 for n>=N_m and -100<=lambda<=0",
            "No effective N_m, no common threshold, and no finite-prefix theorem.",
        ),
        LemmaRow(
            "gkvaol_11_order7_specialization",
            "exact_specialization",
            "ready_to_apply",
            "Order seven vanishes through h^20 and its signed first coefficient is positive.",
            order7["raw_first_term"] + "; Q_7=-H_7 has first coefficient " + str(order7["positive_constant"]),
            "Eventual signed order-seven tail only.",
            order7,
        ),
        LemmaRow(
            "gkvaol_12_nonpromotion_boundary",
            "non_promotion_guard",
            "ready_to_apply",
            "Eventual positivity at every fixed order is not all-shift sign regularity or PF-infinity.",
            "forall m exists N_m is not exists N forall m and is not forall m forall n",
            "Blocks quantifier reversal and finite-prefix promotion.",
        ),
    ]


def build_artifact() -> dict:
    sources = validate_sources()
    constants = order_rows()
    valuations = coefficient_valuation_table()
    stress = permutation_penalty_stress()
    rows = proof_rows(stress, constants)
    source_paths = (HEAT_TILT_SOURCE, HIGHER_THETA_SOURCE, XI_RATIO_SOURCE)
    return {
        "kind": "jensen_window_pf_graded_kernel_vandermonde_all_order_lemma",
        "date": "2026-07-13",
        "status": (
            "exact all-fixed-order graded-kernel leading determinant lemma and "
            "compact-uniform eventual signed tails; no all-shift theorem"
        ),
        "proof_boundary": (
            "This artifact proves the universal first coefficient and, for every "
            "fixed order m, an eventual signed tail uniform on -100<=lambda<=0. "
            "It supplies no effective threshold, all-shift sign regularity, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [path.relative_to(REPO_ROOT).as_posix() for path in source_paths],
        "source_contract": sources,
        "kernel": {
            "entry": "K_(i,j)=exp(-sum_(r>=2)G_r*h^(r-1)*(L-i-j)^r)",
            "L": "2*(m-1)",
            "mixed_entry": (
                "M_h(z,y)=exp(sum_(a,b>=1)c_(a,b)*h^(a+b-1)*z^a*y^b)"
            ),
            "quadratic_mixed_coefficient": "c_(1,1)=2*G_2",
        },
        "valuation_lemma": {
            "coefficient_bound": "ord_h c_(p,q)>=max(p,q)",
            "diagonal_leader": "[h^k]c_(k,k)=(2*G_2)^k/k!",
            "tuple_floor": "sum_i p_i>=binom(m,2) for increasing nonnegative P",
            "permutation_identity": (
                "sum_i max(i,pi(i))=binom(m,2)+(1/2)*sum_i|i-pi(i)|"
            ),
            "coefficient_table": valuations,
        },
        "heat_tilt_all_fixed_order_extension": {
            "primary_source": "https://arxiv.org/abs/2007.13582",
            "published_input": (
                "O'Sullivan Theorem 5.2 gives arbitrary-order saddle expansions "
                "for suitable multipliers"
            ),
            "uniform_suitability": (
                "Cauchy estimates for the compact Gaussian-log family give every "
                "fixed suitability coefficient and remainder uniformly in T"
            ),
            "lambert_recurrence": (
                "F_0(w)=w^2; F_(s+1)(w)=-s*F_s(w)+"
                "w/(1+w)*F_s'(w)"
            ),
            "lambert_derivative": (
                "d^s(w^2)/dk^s=k^(-s)*F_s(w), F_s(w)=O_s(w)"
            ),
            "finite_difference_identity": (
                "Delta^s f(k)=integral_[0,1]^s f^(s)(k+t_1+...+t_s)dt"
            ),
            "remainder_choice": (
                "for fixed s choose R>s; 2^s O_R(w^(3R)/k^R)="
                "o(w/k^s)"
            ),
            "theorem": (
                "Delta^s log R_T^(1)(k)=O_s(log(k)/k^s) for every "
                "fixed s>=2, uniformly for 0<=T<=100"
            ),
            "quantifier_boundary": (
                "the implied constant and chosen truncation may depend on s"
            ),
        },
        "vandermonde_theorem": {
            "degree": "D=binom(m,2)",
            "raw_first_coefficient": (
                "epsilon_m*2^D*prod_(j=1)^(m-1)j!*G_2^D"
            ),
            "signed_first_coefficient": (
                "2^D*prod_(j=1)^(m-1)j!*G_2^D>0"
            ),
            "eventual_tail": (
                "for every fixed m, exists N_m such that Q_(m,n)(lambda)>0 "
                "for n>=N_m and -100<=lambda<=0"
            ),
            "quantifier_boundary": (
                "N_m may depend on m; this does not prove all-shift signs"
            ),
        },
        "order_rows": constants,
        "permutation_stress": stress,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": len(rows),
            "recorded_orders": len(constants),
            "coefficient_valuation_cells": len(valuations),
            "permutation_stress_cases": stress["permutations_checked"],
            "all_fixed_order_leading_theorems": 1,
            "all_fixed_order_eventual_tail_theorems": 1,
            "effective_thresholds": 0,
            "all_shift_theorems": 0,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    theorem = artifact["vandermonde_theorem"]
    order7 = artifact["order_rows"][6]
    stress = artifact["permutation_stress"]
    lines = [
        "# Graded-Kernel Vandermonde Lemma At Every Fixed Order",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact all-fixed-order leading-determinant theorem and",
        "compact-heat uniform eventual signed-tail theorem. This is not a proof of",
        "all-shift sign regularity, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.json",
        "python work/rh_compute/scripts/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.py",
        "```",
        "",
        "## Kernel And Gauge",
        "",
        "Fix an integer `m>=1`, put `D=binom(m,2)`, `L=2(m-1)`,",
        "`z_i=L-i`, and `y_j=j`. The normalized determinant kernel is",
        "",
        "```text",
        artifact["kernel"]["entry"],
        "```",
        "",
        "The binomial identity",
        "",
        "```text",
        "(z-y)^r=z^r+(-y)^r",
        "        +sum_(a,b>=1,a+b=r) binom(r,a)z^a(-y)^b",
        "```",
        "",
        "splits off row and column factors with constant term one. They cannot",
        "alter the first nonzero determinant coefficient. The remaining mixed",
        "kernel is",
        "",
        "```text",
        artifact["kernel"]["mixed_entry"],
        artifact["kernel"]["quadratic_mixed_coefficient"],
        "```",
        "",
        "## Valuation Lemma",
        "",
        "Write `M_h(z,y)=sum c_(p,q)(h)z^p y^q`. A product of `ell` mixed",
        "factors contributing bidegree `(p,q)` has exact `h`-degree",
        "",
        "```text",
        "p+q-ell>=max(p,q), because ell<=min(p,q).",
        "```",
        "",
        "Therefore",
        "",
        "```text",
        "ord_h c_(p,q)>=max(p,q),",
        "[h^k]c_(k,k)=(2*G_2)^k/k!.",
        "```",
        "",
        "Equality in the diagonal formula requires `ell=k` and hence `k`",
        "copies of the sole `(a,b)=(1,1)` quadratic interaction. No `G_r`",
        "with `r>=3` can enter.",
        "",
        "## Cauchy-Binet Floor",
        "",
        "Expand the mixed kernel as a monomial matrix product. For increasing",
        "nonnegative exponent tuples `P,Q` and a coefficient permutation `pi`,",
        "",
        "```text",
        "sum_i max(p_i,q_(pi(i)))",
        " >=max(sum_i p_i,sum_i q_i)>=D.",
        "```",
        "",
        "Equality forces `P=Q=(0,1,...,m-1)`. In that block",
        "",
        "```text",
        "sum_i max(i,pi(i))=D+(1/2)*sum_i|i-pi(i)|,",
        "```",
        "",
        "so only the identity permutation remains at degree `D`. The checker",
        f"stress-tests this identity on all `{stress['permutations_checked']}`",
        "permutations through order eight.",
        "",
        "## Universal First Term",
        "",
        "The two monomial alternants are Vandermonde determinants. Since",
        "`z_i` is decreasing and `y_j` increasing, their product leaves",
        "",
        "```text",
        theorem["raw_first_coefficient"],
        "```",
        "",
        "as the first raw coefficient. Thus `Q_(m,n)=epsilon_m H_(m,n)` has",
        "the positive first coefficient",
        "",
        "```text",
        theorem["signed_first_coefficient"],
        "```",
        "",
        "## Every Fixed Heat-Difference Order",
        "",
        "O'Sullivan's Theorem 5.2 supplies the saddle expansion at arbitrary",
        "fixed truncation order for suitable multipliers. Cauchy estimates in",
        "the published complex disk make the compact Gaussian-log family",
        "uniformly suitable at every fixed order. For `w=W(2k/pi)`, put",
        "",
        "```text",
        "F_0(w)=w^2,",
        "F_(s+1)(w)=-s*F_s(w)+w/(1+w)*F_s'(w).",
        "```",
        "",
        "Induction using `w'=w/(k(1+w))` gives",
        "",
        "```text",
        "d^s(w^2)/dk^s=k^(-s)F_s(w), F_s(w)=O_s(w).",
        "```",
        "",
        "For one fixed difference order `s`, choose the published saddle",
        "truncation `R>s`. The explicit correction terms gain `s` inverse",
        "powers under the integral finite-difference formula. The difference",
        "of the remainder is bounded directly, without differentiating it:",
        "",
        "```text",
        "2^s O_R(w^(3R)/k^R)=o(w/k^s).",
        "```",
        "",
        "Consequently",
        "",
        "```text",
        "Delta^s log R_T^(1)(k)=O_s(log(k)/k^s)",
        "for every fixed s>=2, uniformly for 0<=T<=100.",
        "```",
        "",
        "The constants and the truncation may depend on `s`; this is not a",
        "uniform-in-order estimate.",
        "",
        "The all-order Xi ratio theorem, compact heat-tilt theorem, and",
        "superpolynomial higher-theta suppression now imply",
        "",
        "```text",
        theorem["eventual_tail"],
        "```",
        "",
        "This is an `all fixed m` statement: the finite threshold may depend on",
        "`m`. It neither reverses the quantifiers nor fills any finite prefix.",
        "",
        "## Order Seven",
        "",
        "At the new frontier `D=21` and `epsilon_7=-1`, so",
        "",
        "```text",
        "[h^0,...,h^21]det K=[0,...,0,-52183852646400*G_2^21],",
        order7["signed_first_term"] + " for Q_7=-H_7.",
        "```",
        "",
        "This closes the eventual order-seven tail. Endpoint entry at",
        "`lambda=-100`, its analytic finite-prefix splice, and the all-order",
        "sign-regularity bridge remain separate obligations.",
        "",
        "```text",
        "outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md",
        "outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md",
        "outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote graded-kernel all-order Vandermonde lemma: "
        f"{summary['recorded_orders']} order rows, "
        f"{summary['permutation_stress_cases']} permutation stress cases, "
        f"{summary['coefficient_valuation_cells']} valuation cells"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
