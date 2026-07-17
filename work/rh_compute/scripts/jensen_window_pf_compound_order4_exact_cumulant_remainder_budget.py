#!/usr/bin/env python3
"""Compose the proved formal layers into explicit exact-cumulant remainder budgets."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
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

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    COEFFICIENT_BOUNDS,
    sha256,
)
from jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate import (  # noqa: E402
    COEFFICIENT_BOUNDS as SECOND_COEFFICIENT_BOUNDS,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md"
)
SOURCE_FORMAL_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.json"
)
SOURCE_FORMAL_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.json"
)
SOURCE_NEXT_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.json"
)
SOURCE_NEXT_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.json"
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
FINITE_Q_FLOOR = 9_000
RAY_Q_FLOOR = 10**35
FINITE_FORMAL_SIX_MARGIN_FLOOR = Fraction(18, 1000)
FINITE_EIGHTH_CORRECTION_CAP = Fraction(1, 1000)
FINITE_EXACT_REMAINDER_BUDGET = Fraction(1, 100)
RAY_FORMAL_SIX_MARGIN = Fraction(1, 20)
RAY_EIGHTH_CORRECTION_CAP = Fraction(1, 100)
RAY_EXACT_REMAINDER_BUDGET = Fraction(1, 50)
FINITE_TEN_CORRECTION_CAP = Fraction(1, 10**7)
FINITE_EXACT_TEN_REMAINDER_BUDGET = Fraction(9, 1000)
RAY_TEN_CORRECTION_CAP = Fraction(1, 1000)
RAY_EXACT_TEN_REMAINDER_BUDGET = Fraction(1, 100)


@dataclass(frozen=True)
class BudgetRow:
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


def correction_power(order: int) -> int:
    return 3 if order <= 4 else 2 if order <= 6 else 1


def absolute_coefficient_cap(order: int) -> Fraction:
    floor, ceiling = COEFFICIENT_BOUNDS[order]
    return max(abs(floor), abs(ceiling))


def second_correction_power(order: int) -> int:
    return 4 if order <= 4 else 3 if order <= 6 else 2


def second_absolute_coefficient_cap(order: int) -> Fraction:
    floor, ceiling = SECOND_COEFFICIENT_BOUNDS[order]
    return max(abs(floor), abs(ceiling))


def finite_budget() -> dict:
    formal = load_json(SOURCE_FORMAL_FINITE)
    weakest = Decimal(formal["weakest_margin"]["value"])
    if weakest <= Decimal(FINITE_FORMAL_SIX_MARGIN_FLOOR.numerator) / Decimal(
        FINITE_FORMAL_SIX_MARGIN_FLOOR.denominator
    ):
        raise RuntimeError("finite epsilon-six margin floor failed")
    corrections = {}
    for order in range(2, 9):
        cap = absolute_coefficient_cap(order) / FINITE_Q_FLOOR ** correction_power(order)
        if cap >= FINITE_EIGHTH_CORRECTION_CAP:
            raise RuntimeError(f"finite epsilon-eight correction cap failed at order {order}")
        corrections[str(order)] = {
            "coefficient_cap": str(absolute_coefficient_cap(order)),
            "q_power": correction_power(order),
            "correction_cap_at_q_floor": str(cap),
        }
    post_eight_margin = (
        FINITE_FORMAL_SIX_MARGIN_FLOOR - FINITE_EIGHTH_CORRECTION_CAP
    )
    final_margin = post_eight_margin - FINITE_EXACT_REMAINDER_BUDGET
    if final_margin <= 0:
        raise RuntimeError("finite exact remainder reserve failed")
    return {
        "mode_interval": "2<=u<=20",
        "q_floor": FINITE_Q_FLOOR,
        "q_at_2_lower": None,
        "formal_six_actual_weakest_margin": str(weakest),
        "formal_six_margin_floor": str(FINITE_FORMAL_SIX_MARGIN_FLOOR),
        "epsilon_eight_corrections": corrections,
        "uniform_epsilon_eight_correction_cap": str(FINITE_EIGHTH_CORRECTION_CAP),
        "post_epsilon_eight_margin_floor": str(post_eight_margin),
        "sufficient_exact_remainder_budget": str(FINITE_EXACT_REMAINDER_BUDGET),
        "reserved_final_corridor_margin": str(final_margin),
    }


def ray_budget() -> dict:
    formal = load_json(SOURCE_FORMAL_RAY)
    if formal.get("scalar_geometry", {}).get("formal_transfer") != (
        "|R_r^[6]-F_r|<=22000000*u^6/q<1/(20u), u>=20"
    ):
        raise RuntimeError("formal asymptotic margin source changed")
    largest_coefficient_cap = max(
        absolute_coefficient_cap(order) for order in range(2, 9)
    )
    if largest_coefficient_cap != Fraction(61, 10):
        raise RuntimeError("unexpected largest next-parity coefficient cap")
    endpoint_left = 100 * RAY_START * largest_coefficient_cap
    if endpoint_left >= RAY_Q_FLOOR:
        raise RuntimeError("asymptotic epsilon-eight correction endpoint failed")
    monotonicity_margin = Fraction(4) - Fraction(1, RAY_START)
    if monotonicity_margin <= 0:
        raise RuntimeError("q/u monotonicity failed")
    post_eight_margin = RAY_FORMAL_SIX_MARGIN - RAY_EIGHTH_CORRECTION_CAP
    final_margin = post_eight_margin - RAY_EXACT_REMAINDER_BUDGET
    if not (post_eight_margin == Fraction(1, 25) and final_margin == Fraction(1, 50)):
        raise RuntimeError("asymptotic remainder arithmetic failed")
    return {
        "mode_interval": "u>=20",
        "q_floor": str(RAY_Q_FLOOR),
        "largest_coefficient_cap": str(largest_coefficient_cap),
        "epsilon_eight_correction": "scaled |kappa_r^[8]-kappa_r^[6]|<1/(100u)",
        "endpoint_comparison_left": str(endpoint_left),
        "endpoint_comparison_right": str(RAY_Q_FLOOR),
        "q_over_u_log_derivative_lower": str(monotonicity_margin),
        "formal_six_margin": "1/(20u)",
        "post_epsilon_eight_margin": "1/(25u)",
        "sufficient_exact_remainder_budget": "1/(50u)",
        "reserved_final_corridor_margin": "1/(50u)",
    }


def epsilon_ten_budget() -> dict:
    finite_source = load_json(SOURCE_SECOND_FINITE)
    ray_source = load_json(SOURCE_SECOND_RAY)
    if finite_source.get("parameters", {}).get("block_count") != 3600:
        raise RuntimeError("finite second-next coefficient source is not closed")
    if ray_source.get("summary", {}).get("global_second_next_layer_closed") is not True:
        raise RuntimeError("global second-next coefficient layer is not closed")

    finite_corrections = {}
    for order in range(2, 9):
        cap = (
            second_absolute_coefficient_cap(order)
            / FINITE_Q_FLOOR ** second_correction_power(order)
        )
        if cap >= FINITE_TEN_CORRECTION_CAP:
            raise RuntimeError(f"finite epsilon-ten correction cap failed at order {order}")
        finite_corrections[str(order)] = {
            "coefficient_cap": str(second_absolute_coefficient_cap(order)),
            "q_power": second_correction_power(order),
            "correction_cap_at_q_floor": str(cap),
        }
    finite_total_error = (
        FINITE_TEN_CORRECTION_CAP + FINITE_EXACT_TEN_REMAINDER_BUDGET
    )
    finite_budget_margin = FINITE_EXACT_REMAINDER_BUDGET - finite_total_error
    finite_corridor_margin = (
        FINITE_FORMAL_SIX_MARGIN_FLOOR
        - FINITE_EIGHTH_CORRECTION_CAP
        - finite_total_error
    )
    if not (
        finite_budget_margin > 0
        and finite_corridor_margin > Fraction(7, 1000)
    ):
        raise RuntimeError("finite epsilon-ten remainder reserve failed")

    largest_cap = max(
        second_absolute_coefficient_cap(order) for order in range(2, 9)
    )
    if largest_cap != Fraction(25, 4):
        raise RuntimeError("unexpected largest second-next coefficient cap")
    endpoint_left = 1000 * RAY_START * largest_cap
    endpoint_right = RAY_Q_FLOOR**2
    if endpoint_left >= endpoint_right:
        raise RuntimeError("ray epsilon-ten correction endpoint failed")
    monotonicity_margin = Fraction(8) - Fraction(1, RAY_START)
    if monotonicity_margin <= 0:
        raise RuntimeError("q^2/u monotonicity failed")
    ray_total_numerator = RAY_TEN_CORRECTION_CAP + RAY_EXACT_TEN_REMAINDER_BUDGET
    ray_budget_margin = RAY_EXACT_REMAINDER_BUDGET - ray_total_numerator
    ray_corridor_margin = Fraction(1, 25) - ray_total_numerator
    if not (
        ray_total_numerator == Fraction(11, 1000)
        and ray_budget_margin == Fraction(9, 1000)
        and ray_corridor_margin == Fraction(29, 1000)
    ):
        raise RuntimeError("ray epsilon-ten remainder arithmetic failed")
    return {
        "finite": {
            "mode_interval": "2<=u<=20",
            "corrections": finite_corrections,
            "uniform_epsilon_ten_correction_cap": str(FINITE_TEN_CORRECTION_CAP),
            "sufficient_exact_minus_epsilon_ten_budget": str(
                FINITE_EXACT_TEN_REMAINDER_BUDGET
            ),
            "total_exact_minus_epsilon_eight_cap": str(finite_total_error),
            "margin_inside_previous_exact_remainder_budget": str(
                finite_budget_margin
            ),
            "reserved_final_corridor_margin": str(finite_corridor_margin),
        },
        "ray": {
            "mode_interval": "u>=20",
            "largest_coefficient_cap": str(largest_cap),
            "epsilon_ten_correction": (
                "scaled |kappa_r^[10]-kappa_r^[8]|<1/(1000u)"
            ),
            "endpoint_comparison_left": str(endpoint_left),
            "endpoint_comparison_right": str(endpoint_right),
            "q_squared_over_u_log_derivative_lower": str(monotonicity_margin),
            "sufficient_exact_minus_epsilon_ten_budget": "1/(100u)",
            "total_exact_minus_epsilon_eight_cap": "11/(1000u)",
            "margin_inside_previous_exact_remainder_budget": "9/(1000u)",
            "reserved_final_corridor_margin": "29/(1000u)",
        },
    }


RAY_START = 20


def build_artifact() -> dict:
    flint.ctx.prec = 256
    q2 = flint.arb.pi() * flint.arb(8).exp()
    if not bool(q2 > FINITE_Q_FLOOR):
        raise RuntimeError("finite q floor failed")
    finite = finite_budget()
    finite["q_at_2_lower"] = arb_lower_text(q2)
    ray = ray_budget()
    epsilon_ten = epsilon_ten_budget()
    next_ray = load_json(SOURCE_NEXT_RAY)
    if next_ray.get("summary", {}).get("global_coefficient_layer_closed") is not True:
        raise RuntimeError("global next-parity coefficient layer is not closed")
    next_finite = load_json(SOURCE_NEXT_FINITE)
    if next_finite.get("parameters", {}).get("block_count") != 1800:
        raise RuntimeError("finite next-parity coefficient source is not closed")

    rows = [
        BudgetRow(
            id="co4ecrb_01_global_epsilon_eight_layer",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The exact formal epsilon-eight correction has explicit signed coefficient bounds for every u>=2.",
            formula="scaled(kappa_r^[8]-kappa_r^[6])=q^-p_r*C_r(u)",
            proof_boundary="Exact formal coefficient layer only.",
            diagnostics={"powers": {"2-4": 3, "5-6": 2, "7-8": 1}},
        ),
        BudgetRow(
            id="co4ecrb_02_finite_correction_cap",
            role="exact_inequality",
            readiness="ready_to_apply",
            claim="On the finite interval, every epsilon-eight scaled correction is smaller than 1/1000.",
            formula="scaled |kappa_r^[8]-kappa_r^[6]|<1/1000, 2<=u<=20",
            proof_boundary="Formal epsilon-eight correction only.",
            diagnostics=finite["epsilon_eight_corrections"],
        ),
        BudgetRow(
            id="co4ecrb_03_finite_exact_reserve",
            role="exact_theorem_reduction",
            readiness="ready_to_apply",
            claim="A uniform scaled exact-minus-epsilon-eight error below 1/100 is sufficient for every finite corridor.",
            formula="scaled |kappa_r-kappa_r^[8]|<1/100 on 2<=u<=20",
            proof_boundary="Sufficient theorem target, not a proved exact-density estimate.",
            diagnostics=finite,
        ),
        BudgetRow(
            id="co4ecrb_04_ray_correction_cap",
            role="exact_analytic_inequality",
            readiness="ready_to_apply",
            claim="On the asymptotic ray, exponential q growth makes every epsilon-eight correction smaller than 1/(100u).",
            formula=ray["epsilon_eight_correction"],
            proof_boundary="Formal epsilon-eight correction only.",
            diagnostics=ray,
        ),
        BudgetRow(
            id="co4ecrb_05_ray_exact_reserve",
            role="exact_theorem_reduction",
            readiness="ready_to_apply",
            claim="A scaled exact-minus-epsilon-eight error below 1/(50u) is sufficient for every asymptotic corridor.",
            formula="scaled |kappa_r-kappa_r^[8]|<1/(50u), u>=20",
            proof_boundary="Sufficient theorem target, not a proved exact-density estimate.",
            diagnostics=ray,
        ),
        BudgetRow(
            id="co4ecrb_06_global_epsilon_ten_layer",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The exact formal epsilon-nine/ten correction has explicit signed coefficient bounds for every u>=2.",
            formula="scaled(kappa_r^[10]-kappa_r^[8])=q^-s_r*D_r(u)",
            proof_boundary="Exact formal coefficient layer only.",
            diagnostics={"powers": {"2-4": 4, "5-6": 3, "7-8": 2}},
        ),
        BudgetRow(
            id="co4ecrb_07_epsilon_ten_correction_caps",
            role="exact_analytic_inequality",
            readiness="ready_to_apply",
            claim="The complete epsilon-ten correction is below 10^-7 on the finite interval and below 1/(1000u) on the ray.",
            formula="finite <10^-7; ray <1/(1000u)",
            proof_boundary="Formal epsilon-nine/ten correction only.",
            diagnostics=epsilon_ten,
        ),
        BudgetRow(
            id="co4ecrb_08_epsilon_ten_exact_reserve",
            role="exact_theorem_reduction",
            readiness="ready_to_apply",
            claim="It is sufficient to bound the scaled exact-minus-epsilon-ten residual by 9/1000 on the finite interval and 1/(100u) on the ray.",
            formula="finite <9/1000; ray <1/(100u)",
            proof_boundary="Sufficient theorem target, not a proved exact-density estimate.",
            diagnostics=epsilon_ten,
        ),
        BudgetRow(
            id="co4ecrb_09_central_tail_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove the sharpened exact-minus-epsilon-ten budgets by a cancellation-preserving central remainder and two adaptive tails.",
            formula="finite <9/1000; ray <1/(100u), simultaneously for r=2,...,8",
            proof_boundary="Open exact-density theorem; no exact cumulant corridor is asserted.",
        ),
    ]
    source_paths = (
        SOURCE_FORMAL_FINITE,
        SOURCE_FORMAL_RAY,
        SOURCE_NEXT_FINITE,
        SOURCE_NEXT_RAY,
        SOURCE_SECOND_FINITE,
        SOURCE_SECOND_RAY,
    )
    return {
        "kind": "jensen_window_pf_compound_order4_exact_cumulant_remainder_budget",
        "date": "2026-07-13",
        "status": "exact formal composition and sharpened exact-density remainder target",
        "proof_boundary": (
            "This artifact composes the epsilon-six corridor margins with the proved "
            "global epsilon-eight and epsilon-ten coefficient layers and derives sufficient exact-density "
            "remainder budgets. It does not prove those central or tail remainders, the "
            "exact cumulant corridors, curvature ray, order-four entry, PF-infinity, RH, "
            "or Lambda<=0."
        ),
        "finite_budget": finite,
        "asymptotic_budget": ray,
        "epsilon_ten_budget": epsilon_ten,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 8,
            "open_analytic_rows": 1,
            "finite_exact_remainder_budget": "1/100",
            "ray_exact_remainder_budget": "1/(50u)",
            "global_epsilon_eight_layer_closed": True,
            "global_epsilon_ten_layer_closed": True,
            "finite_exact_minus_epsilon_ten_budget": "9/1000",
            "ray_exact_minus_epsilon_ten_budget": "1/(100u)",
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            path.relative_to(REPO_ROOT).as_posix(): sha256(path)
            for path in source_paths
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.py"
        ),
        "remaining_target": (
            "Prove scaled exact-minus-epsilon-ten cumulant errors below 9/1000 on "
            "2<=u<=20 and below 1/(100u) on u>=20, including both tails."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite_budget"]
    ray = artifact["asymptotic_budget"]
    epsilon_ten = artifact["epsilon_ten_budget"]
    lines = [
        "# Jensen-Window PF Order-Four Exact Cumulant Remainder Budget",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact formal composition and sharpened exact-density remainder target.",
        "This is not a proof of the exact cumulant corridors, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_remainder_budget`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.py",
        "```",
        "",
        "## Finite Budget",
        "",
        f"Arb gives `q(2)>{finite['q_floor']}`. The epsilon-six formal corridor",
        f"margin is greater than `{finite['formal_six_margin_floor']}` throughout",
        "`2<=u<=20`, while the global next-parity coefficient theorem gives",
        "",
        "```text",
        "scaled |kappa_r^[8]-kappa_r^[6]|<1/1000, r=2,...,8.",
        "```",
        "",
        "Thus the epsilon-eight formal model remains at least `17/1000` inside",
        "every corridor. It is sufficient to prove",
        "",
        "```text",
        "scaled |kappa_r-kappa_r^[8]|<1/100, 2<=u<=20.",
        "```",
        "",
        f"This still reserves a final corridor margin of `{finite['reserved_final_corridor_margin']}`.",
        "",
        "## Ray Budget",
        "",
        "The epsilon-six formal ray lies `1/(20u)` inside every corridor. Since",
        "the largest next-parity coefficient cap is `61/10`, exponential q growth",
        "proves",
        "",
        "```text",
        ray["epsilon_eight_correction"],
        "```",
        "",
        "The epsilon-eight formal model is therefore at least `1/(25u)` inside",
        "every corridor. It is sufficient to prove",
        "",
        "```text",
        "scaled |kappa_r-kappa_r^[8]|<1/(50u), u>=20.",
        "```",
        "",
        "## Epsilon-Ten Sharpening",
        "",
        "The globally certified second-next coefficients give",
        "",
        "```text",
        "scaled |kappa_r^[10]-kappa_r^[8]|<1/10000000, 2<=u<=20,",
        epsilon_ten["ray"]["epsilon_ten_correction"] + ", u>=20.",
        "```",
        "",
        "It is therefore sufficient to prove the cancellation-preserving residual",
        "bounds",
        "",
        "```text",
        "scaled |kappa_r-kappa_r^[10]|<9/1000, 2<=u<=20,",
        "scaled |kappa_r-kappa_r^[10]|<1/(100u), u>=20.",
        "```",
        "",
        "These leave final corridor margins greater than `7/1000` on the finite",
        "interval and `29/(1000u)` on the asymptotic ray.",
        "",
        "## Remaining Boundary",
        "",
        "These are sufficient budgets, not exact-density estimates. The remaining",
        "work is a cancellation-preserving central residual after epsilon ten",
        "plus rigorous left and right tails. Raw high-moment interval boxes are not",
        "used because their dependency widths destroy the q-scaled cancellations.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
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
        "certified order-four exact cumulant remainder budget: "
        "9 rows, 8 exact rows, finite epsilon-ten budget 9/1000, "
        "ray epsilon-ten budget 1/(100u), 1 open central-tail theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
