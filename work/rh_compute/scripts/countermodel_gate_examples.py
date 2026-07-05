#!/usr/bin/env python3
"""Executable countermodel gates for the RH proof programme.

These examples are deliberately small.  They do not model the zeta kernel.
They test whether a proposed proof step is relying on information that is too
local or too finite to imply the Newman/RH conclusion.
"""

from __future__ import annotations

import argparse
import json
import math
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path


DEFAULT_LAMBDAS = ("0", "1e-6", "1e-4", "1e-2", "1e-1")


def dec(value: object) -> Decimal:
    return Decimal(str(value))


def short(value: Decimal, digits: int = 18) -> str:
    return f"{value:.{digits}E}"


def short_fraction(value: Fraction, digits: int = 18) -> str:
    as_decimal = Decimal(value.numerator) / Decimal(value.denominator)
    return short(as_decimal, digits)


def exact_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def load_coeff_rows(path: Path) -> dict[Decimal, dict]:
    rows: dict[Decimal, dict] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") != "coefficients":
                continue
            rows[dec(row["lam"])] = row
    return rows


def a_values(row: dict) -> list[Decimal]:
    return [dec(value) for value in row["A"]]


def c_values(a: list[Decimal]) -> list[Decimal]:
    return [value / Decimal(math.factorial(k)) for k, value in enumerate(a)]


def heat_birth_gate(tau: Decimal) -> dict:
    """Check the exact local square-root birth model.

    F_tau(t) = t^2 - 2 tau satisfies F_tau = -F_tt.  For tau > 0 its two real
    zeros have gap g(tau) = 2 sqrt(2 tau), and g' = 4/g.  Thus the local
    repulsion law is compatible with a birth from a complex pair.
    """
    if tau <= 0:
        raise ValueError("tau must be positive")
    sqrt_2tau = (Decimal(2) * tau).sqrt()
    gap = Decimal(2) * sqrt_2tau
    gap_derivative = Decimal(2) / sqrt_2tau
    repulsion_rhs = Decimal(4) / gap
    return {
        "name": "local_heat_birth",
        "model": "F_tau(t)=t^2-2*tau",
        "pde_check": "partial_tau F = -2 = -partial_tt F",
        "tau": short(tau),
        "gap": short(gap),
        "gap_derivative": short(gap_derivative),
        "four_over_gap": short(repulsion_rhs),
        "repulsion_identity_holds": gap_derivative == repulsion_rhs,
        "gate": "Any proof using only g'=4/g on the real side cannot exclude positive birth.",
    }


def solve_fraction(A: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    n = len(A)
    matrix = [list(row) + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if matrix[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            raise ArithmeticError("singular rational system")
        if pivot != col:
            matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        pivot_value = matrix[col][col]
        for j in range(col, n + 1):
            matrix[col][j] /= pivot_value
        for row in range(n):
            if row == col:
                continue
            factor = matrix[row][col]
            if factor == 0:
                continue
            for j in range(col, n + 1):
                matrix[row][j] -= factor * matrix[col][j]
    return [matrix[i][n] for i in range(n)]


def det_fraction(A: list[list[Fraction]]) -> Fraction:
    n = len(A)
    matrix = [list(row) for row in A]
    sign = 1
    det = Fraction(1)
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if matrix[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
            sign *= -1
        pivot_value = matrix[col][col]
        det *= pivot_value
        for row in range(col + 1, n):
            factor = matrix[row][col] / pivot_value
            if factor == 0:
                continue
            for j in range(col, n):
                matrix[row][j] -= factor * matrix[col][j]
    return det * sign


def hankel_matrix(moments: list[Fraction], size: int) -> list[list[Fraction]]:
    return [[moments[i + j] for j in range(size)] for i in range(size)]


def signed_hankel_sigma(m: int) -> int:
    return -1 if (m * (m + 1) // 2) % 2 else 1


def moment_recurrence_prefix_gate(order: int) -> dict:
    """Preserve a positive recurrence prefix, then break the next moment gate.

    The factorial moments m_n = n! are Stieltjes moments of exp(-x) dx on
    [0, infinity), so every finite recurrence datum before the edited moment is
    genuinely positive.  The order-N Arb scout uses shifted moments only through
    index 2N-1.  We preserve that whole prefix and choose the next even moment
    below the exact Schur-complement threshold for the next Hankel matrix.
    """
    if order < 1:
        raise ValueError("moment recurrence order must be positive")

    moments = [Fraction(math.factorial(n), 1) for n in range(2 * order + 2)]
    previous_size = order
    next_size = order + 1

    previous_hankel = hankel_matrix(moments, previous_size)
    leading_dets = [
        det_fraction(hankel_matrix(moments, size))
        for size in range(1, previous_size + 1)
    ]
    previous_det = det_fraction(previous_hankel)
    if not all(det > 0 for det in leading_dets):
        raise AssertionError("factorial moment prefix should have positive Hankel determinants")

    vector = [moments[order + i] for i in range(order)]
    solution = solve_fraction(previous_hankel, vector)
    threshold = sum(vector[i] * solution[i] for i in range(order))
    factorial_next = moments[2 * order]
    adversarial_next = threshold / 2
    moments[2 * order] = adversarial_next

    next_hankel = hankel_matrix(moments, next_size)
    next_det = det_fraction(next_hankel)
    return {
        "name": "finite_moment_recurrence_prefix_extension",
        "preserved_recurrence_order": order,
        "preserved_moment_indices": f"0..{2 * order - 1}",
        "edited_moment_index": 2 * order,
        "previous_hankel_det": short_fraction(previous_det),
        "all_preserved_leading_hankel_dets_positive": all(det > 0 for det in leading_dets),
        "previous_hankel_positive": previous_det > 0,
        "schur_threshold_for_next_moment": short_fraction(threshold),
        "factorial_next_moment": short_fraction(factorial_next),
        "adversarial_positive_next_moment": short_fraction(adversarial_next),
        "adversarial_next_moment_is_positive": adversarial_next > 0,
        "next_hankel_det_after_extension": short_fraction(next_det),
        "next_moment_gate_breaks": next_det < 0,
        "gate": (
            "Finite positive moment/recurrence data, even through the current "
            "certified order, cannot imply an all-order Edrei moment "
            "representation without an additional infinite theorem."
        ),
    }


def stieltjes_multiplier_trap_gate() -> dict:
    """Show that Stieltjes/Hankel moment positivity is not coefficient PF.

    The measure 10*delta_0 + delta_1 + delta_2 is positive on [0, infinity).
    Its moments therefore form a Stieltjes moment sequence.  But after the
    coefficient-route normalization c_k = mu_k/(2k)!, the first nontrivial
    upper-Toeplitz 2x2 minor is negative.
    """
    support_weights = [(Fraction(0), Fraction(10)), (Fraction(1), Fraction(1)), (Fraction(2), Fraction(1))]
    moments = [
        sum(weight * (point ** k) for point, weight in support_weights)
        for k in range(7)
    ]
    leading_hankel_dets = [
        det_fraction(hankel_matrix(moments, size))
        for size in range(1, 5)
    ]
    c = [moments[k] / Fraction(math.factorial(2 * k), 1) for k in range(3)]
    toeplitz_det = c[1] * c[1] - c[0] * c[2]
    return {
        "name": "stieltjes_moment_multiplier_trap",
        "positive_measure": "10*delta_0 + delta_1 + delta_2 on [0, infinity)",
        "moments_mu_0_to_6": [exact_fraction(value) for value in moments],
        "leading_hankel_dets_size_1_to_4": [exact_fraction(value) for value in leading_hankel_dets],
        "stieltjes_hankel_nonnegative": all(value >= 0 for value in leading_hankel_dets),
        "stieltjes_hankel_strict_through_rank": all(value > 0 for value in leading_hankel_dets[:3]),
        "normalized_c_0_to_2": [exact_fraction(value) for value in c],
        "toeplitz_order2_det_c1_squared_minus_c0_c2": exact_fraction(toeplitz_det),
        "toeplitz_pf_breaks_after_factorial_normalization": toeplitz_det < 0,
        "gate": (
            "A proof of coefficient PF-infinity cannot use only Stieltjes "
            "moment/Hankel positivity of mu_k plus the factor 1/(2k)!; it "
            "needs a theorem whose hypotheses imply Toeplitz total positivity."
        ),
    }


def finite_consecutive_hankel_grid_extension_gate(max_m: int = 4, max_shift: int = 8) -> dict:
    """Preserve a finite shifted-principal signed-Hankel grid, then break next shift.

    The exact sequence a_k = 1/k! has the same alternating shifted-principal
    Hankel determinant sign pattern on the small grid below.  After that finite
    grid is fixed, one positive coefficient beyond the grid can break the next
    m=1 signed-Hankel determinant and the degree-2 Jensen discriminant.
    """
    if max_m < 1:
        raise ValueError("max_m must be at least 1")
    if max_shift < 0:
        raise ValueError("max_shift must be nonnegative")

    preserved_max_index = max_shift + 2 * max_m
    base = [Fraction(1, math.factorial(k)) for k in range(preserved_max_index + 2)]

    signed_values: list[Fraction] = []
    for m in range(max_m + 1):
        size = m + 1
        sigma = signed_hankel_sigma(m)
        for shift in range(max_shift + 1):
            determinant = det_fraction(
                [[base[shift + i + j] for j in range(size)] for i in range(size)]
            )
            signed_value = sigma * determinant
            signed_values.append(signed_value)
            if signed_value <= 0:
                raise AssertionError(
                    f"base sequence failed grid at m={m}, shift={shift}: {signed_value}"
                )

    edited = list(base)
    break_shift = preserved_max_index - 1
    edited_index = preserved_max_index + 1
    adversarial_next = Fraction(2) * base[preserved_max_index] ** 2 / base[preserved_max_index - 1]
    edited[edited_index] = adversarial_next

    break_det = det_fraction(
        [[edited[break_shift + i + j] for j in range(2)] for i in range(2)]
    )
    signed_break = signed_hankel_sigma(1) * break_det
    jensen_degree2_discriminant = Fraction(4) * (
        edited[break_shift + 1] ** 2
        - edited[break_shift] * edited[break_shift + 2]
    )

    return {
        "name": "finite_consecutive_hankel_grid_extension",
        "base_sequence": "a_k=1/k!",
        "max_m_preserved": max_m,
        "max_shift_preserved": max_shift,
        "preserved_coefficient_indices": f"0..{preserved_max_index}",
        "preserved_grid_count": len(signed_values),
        "all_preserved_grid_values_positive": all(value > 0 for value in signed_values),
        "minimum_preserved_signed_grid_value": short_fraction(min(signed_values)),
        "edited_coefficient_index": edited_index,
        "adversarial_positive_next_coefficient": exact_fraction(adversarial_next),
        "adversarial_next_coefficient_is_positive": adversarial_next > 0,
        "broken_next_shift": break_shift,
        "signed_hankel_m1_value_after_extension": short_fraction(signed_break),
        "next_shift_signed_hankel_breaks": signed_break < 0,
        "jensen_degree2_discriminant_after_extension": short_fraction(jensen_degree2_discriminant),
        "jensen_hyperbolicity_breaks": jensen_degree2_discriminant < 0,
        "gate": (
            "A finite shifted-principal signed-Hankel grid, even when exact and "
            "strictly positive, cannot be promoted to all shifts, all-order "
            "sign-regularity, or Jensen hyperbolicity without an additional theorem."
        ),
    }


def finite_prefix_gates(row: dict, prefix_k: int | None = None) -> dict:
    """Build positive adversarial extensions preserving the current prefix."""
    a = a_values(row)
    c = c_values(a)
    max_k = int(row["max_k"])
    k = max_k if prefix_k is None else min(prefix_k, max_k)
    if k < 1:
        raise ValueError("prefix must contain at least k=0,1")

    # Coefficient-PF / Toeplitz trap for c_k:
    # det [[c_k, c_{k+1}], [c_{k-1}, c_k]] = c_k^2 - c_{k-1} c_{k+1}.
    # Choosing c_{k+1} larger than c_k^2/c_{k-1} breaks order-2 PF while
    # preserving all coefficients through k.
    adversarial_c_next = Decimal(2) * c[k] * c[k] / c[k - 1]
    toeplitz_det = c[k] * c[k] - c[k - 1] * adversarial_c_next

    # Signed-Hankel / Jensen trap for A_k:
    # For m=1, the observed signed-Hankel convention requires
    # -(A_{s} A_{s+2} - A_{s+1}^2) > 0.  Set s=k-1 and choose A_{k+1}
    # larger than A_k^2/A_{k-1}; this also makes the degree-2 Jensen
    # discriminant at shift k-1 negative.
    adversarial_a_next = Decimal(2) * a[k] * a[k] / a[k - 1]
    hankel_det = a[k - 1] * adversarial_a_next - a[k] * a[k]
    signed_hankel_value = -hankel_det
    jensen_degree2_discriminant = Decimal(4) * (
        a[k] * a[k] - a[k - 1] * adversarial_a_next
    )

    return {
        "name": "finite_prefix_extension",
        "lambda": str(row["lam"]),
        "preserved_prefix_through_k": k,
        "toeplitz_order2_det_after_positive_c_extension": short(toeplitz_det),
        "toeplitz_pf_breaks": toeplitz_det < 0,
        "signed_hankel_m1_value_after_positive_A_extension": short(signed_hankel_value),
        "signed_hankel_breaks": signed_hankel_value < 0,
        "jensen_degree2_discriminant_after_same_A_extension": short(jensen_degree2_discriminant),
        "jensen_hyperbolicity_breaks": jensen_degree2_discriminant < 0,
        "adversarial_c_next": short(adversarial_c_next),
        "adversarial_A_next": short(adversarial_a_next),
        "gate": (
            "Any proof step using only a finite coefficient prefix is invalid "
            "unless it adds a structural all-order hypothesis."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coeff-cache",
        type=Path,
        default=Path("work/rh_compute/results/repro_hankel_15c_coefficients.jsonl"),
    )
    parser.add_argument(
        "--lambdas",
        default=",".join(DEFAULT_LAMBDAS),
        help="Comma-separated lambda list to load from the coefficient cache.",
    )
    parser.add_argument("--prefix-k", type=int, default=None)
    parser.add_argument(
        "--moment-order",
        type=int,
        default=12,
        help="Positive recurrence-prefix order for the exact rational moment gate.",
    )
    parser.add_argument("--tau", default="0.125", help="Positive tau for the birth model.")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    args = parser.parse_args()

    getcontext().prec = 120
    rows = load_coeff_rows(args.coeff_cache)
    targets = [dec(part.strip()) for part in args.lambdas.split(",") if part.strip()]

    results: list[dict] = [
        heat_birth_gate(dec(args.tau)),
        moment_recurrence_prefix_gate(args.moment_order),
        stieltjes_multiplier_trap_gate(),
        finite_consecutive_hankel_grid_extension_gate(),
    ]
    for lam in targets:
        if lam not in rows:
            raise KeyError(f"lambda {lam} not found in {args.coeff_cache}")
        results.append(finite_prefix_gates(rows[lam], args.prefix_k))

    for result in results:
        if result["name"] == "local_heat_birth":
            if not result["repulsion_identity_holds"]:
                raise AssertionError(result)
        elif result["name"] == "finite_prefix_extension":
            required = (
                result["toeplitz_pf_breaks"],
                result["signed_hankel_breaks"],
                result["jensen_hyperbolicity_breaks"],
            )
            if not all(required):
                raise AssertionError(result)
        elif result["name"] == "finite_moment_recurrence_prefix_extension":
            required = (
                result["all_preserved_leading_hankel_dets_positive"],
                result["previous_hankel_positive"],
                result["adversarial_next_moment_is_positive"],
                result["next_moment_gate_breaks"],
            )
            if not all(required):
                raise AssertionError(result)
        elif result["name"] == "finite_consecutive_hankel_grid_extension":
            required = (
                result["all_preserved_grid_values_positive"],
                result["adversarial_next_coefficient_is_positive"],
                result["next_shift_signed_hankel_breaks"],
                result["jensen_hyperbolicity_breaks"],
            )
            if not all(required):
                raise AssertionError(result)
        else:
            required = (
                result["stieltjes_hankel_nonnegative"],
                result["stieltjes_hankel_strict_through_rank"],
                result["toeplitz_pf_breaks_after_factorial_normalization"],
            )
            if not all(required):
                raise AssertionError(result)

    if args.json:
        print(json.dumps({"ok": True, "results": results}, indent=2, sort_keys=True))
    else:
        for result in results:
            print(f"OK countermodel gate: {result['name']}")
            if result["name"] == "finite_prefix_extension":
                print(
                    "  "
                    f"lambda={result['lambda']}, "
                    f"prefix k<={result['preserved_prefix_through_k']}, "
                    "Toeplitz/Hankel/Jensen all break after positive one-term extension"
                )
            if result["name"] == "finite_moment_recurrence_prefix_extension":
                print(
                    "  "
                    f"recurrence order<={result['preserved_recurrence_order']}, "
                    f"moments {result['preserved_moment_indices']} preserved, "
                    "next Hankel/moment gate breaks after positive one-term extension"
                )
            if result["name"] == "finite_consecutive_hankel_grid_extension":
                print(
                    "  "
                    f"base {result['base_sequence']}, "
                    f"m<={result['max_m_preserved']}, "
                    f"shifts<={result['max_shift_preserved']} preserved, "
                    "next shifted m=1 signed-Hankel/Jensen gate breaks"
                )
            if result["name"] == "stieltjes_moment_multiplier_trap":
                print(
                    "  "
                    "positive Stieltjes moments stay Hankel-nonnegative, "
                    "but c_k=mu_k/(2k)! fails an order-2 Toeplitz/PF minor"
                )
        print(f"validated {len(results)} countermodel gate examples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
