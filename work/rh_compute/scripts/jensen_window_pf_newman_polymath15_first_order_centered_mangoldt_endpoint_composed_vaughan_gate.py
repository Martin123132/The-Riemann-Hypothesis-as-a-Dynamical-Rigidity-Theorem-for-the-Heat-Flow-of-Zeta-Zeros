#!/usr/bin/env python3
"""Build the centered Mangoldt endpoint-composed Vaughan gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "mangoldt_endpoint_composed_vaughan_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "contact_centering": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_contact_centering_gate.json"
    ),
    "cutoff_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_adjacent_cutoff_transport_gate.json"
    ),
    "corrected_mangoldt": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_abel_contact_gate.json"
    ),
}

ComplexQ = tuple[Fraction, Fraction]
LogCoeff = dict[int, int]
LogMoment = dict[int, ComplexQ]


@dataclass(frozen=True)
class GateRow:
    id: str
    claim: str
    readiness: str
    exact_statement: str
    implication: str
    boundary: str


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def complex_text(value: ComplexQ) -> dict[str, str]:
    return {
        "real": fraction_text(value[0]),
        "imag": fraction_text(value[1]),
    }


def moment_text(value: LogMoment) -> dict[str, dict[str, str]]:
    return {
        str(prime): complex_text(coefficient)
        for prime, coefficient in sorted(value.items())
    }


def real_log_text(value: dict[int, Fraction]) -> dict[str, str]:
    return {
        str(prime): fraction_text(coefficient)
        for prime, coefficient in sorted(value.items())
        if coefficient
    }


def add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def negate(value: ComplexQ) -> ComplexQ:
    return -value[0], -value[1]


def scale(value: ComplexQ, scalar: Fraction | int) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def moment_add(left: LogMoment, right: LogMoment) -> LogMoment:
    result = dict(left)
    for prime, value in right.items():
        result[prime] = add(
            result.get(prime, (Fraction(0), Fraction(0))),
            value,
        )
        if result[prime] == (Fraction(0), Fraction(0)):
            del result[prime]
    return result


def moment_negate(value: LogMoment) -> LogMoment:
    return {prime: negate(coefficient) for prime, coefficient in value.items()}


def moment_subtract(left: LogMoment, right: LogMoment) -> LogMoment:
    return moment_add(left, moment_negate(right))


def moment_multiply(value: LogMoment, scalar: ComplexQ) -> LogMoment:
    return {
        prime: multiply(scalar, coefficient)
        for prime, coefficient in value.items()
        if multiply(scalar, coefficient) != (Fraction(0), Fraction(0))
    }


def real_projection(value: LogMoment) -> dict[int, Fraction]:
    return {
        prime: coefficient[0]
        for prime, coefficient in value.items()
        if coefficient[0]
    }


def coeff_add(
    left: LogCoeff,
    right: LogCoeff,
    sign: int = 1,
) -> LogCoeff:
    result = dict(left)
    for prime, value in right.items():
        result[prime] = result.get(prime, 0) + sign * value
        if result[prime] == 0:
            del result[prime]
    return result


def factor_vector(number: int) -> LogCoeff:
    return {
        int(prime): int(exponent)
        for prime, exponent in sp.factorint(number).items()
    }


def lambda_vector(number: int) -> LogCoeff:
    factors = factor_vector(number)
    if number <= 1 or len(factors) != 1:
        return {}
    return {next(iter(factors)): 1}


def mobius_value(number: int) -> int:
    return int(sp.mobius(number))


def convolve_vector_scalar(
    vector_sequence: list[LogCoeff],
    scalar_sequence: list[int],
    limit: int,
) -> list[LogCoeff]:
    result: list[LogCoeff] = [{} for _ in range(limit + 1)]
    for left in range(1, limit + 1):
        coefficient = vector_sequence[left]
        if not coefficient:
            continue
        for right in range(1, limit // left + 1):
            scalar = scalar_sequence[right]
            if not scalar:
                continue
            index = left * right
            for prime, exponent in coefficient.items():
                result[index][prime] = (
                    result[index].get(prime, 0) + scalar * exponent
                )
                if result[index][prime] == 0:
                    del result[index][prime]
    return result


def vaughan_coefficients(
    limit: int,
    cutoff_lambda: int,
    cutoff_mu: int,
) -> dict[str, list[LogCoeff]]:
    one = [0] + [1] * limit
    mu = [0] + [mobius_value(number) for number in range(1, limit + 1)]
    mu_low = [
        mu[number] if number <= cutoff_mu else 0
        for number in range(limit + 1)
    ]
    mu_high = [
        mu[number] if number > cutoff_mu else 0
        for number in range(limit + 1)
    ]
    logarithm = [{}] + [
        factor_vector(number) for number in range(1, limit + 1)
    ]
    lambda_low = [{}] + [
        lambda_vector(number) if number <= cutoff_lambda else {}
        for number in range(1, limit + 1)
    ]
    lambda_high = [{}] + [
        lambda_vector(number) if number > cutoff_lambda else {}
        for number in range(1, limit + 1)
    ]

    low = convolve_vector_scalar(lambda_low, one, limit)
    type_i_main = convolve_vector_scalar(
        convolve_vector_scalar(logarithm, mu_low, limit),
        one,
        limit,
    )
    type_i_correction = convolve_vector_scalar(
        convolve_vector_scalar(
            convolve_vector_scalar(lambda_low, mu_low, limit),
            one,
            limit,
        ),
        one,
        limit,
    )
    type_ii = convolve_vector_scalar(
        convolve_vector_scalar(
            convolve_vector_scalar(lambda_high, mu_high, limit),
            one,
            limit,
        ),
        one,
        limit,
    )
    type_i = [
        coeff_add(type_i_main[number], type_i_correction[number], -1)
        for number in range(limit + 1)
    ]
    for number in range(1, limit + 1):
        rebuilt = coeff_add(
            coeff_add(low[number], type_i[number]),
            type_ii[number],
        )
        if rebuilt != logarithm[number]:
            raise RuntimeError(
                "Vaughan coefficient identity failed at "
                f"n={number}, U={cutoff_lambda}, V={cutoff_mu}"
            )
    return {
        "log": logarithm,
        "low": low,
        "type_i": type_i,
        "type_ii": type_ii,
    }


def deterministic_carrier(number: int) -> ComplexQ:
    return (
        Fraction((13 * number) % 29 - 14, number + 9),
        Fraction((17 * number) % 31 - 15, number + 11),
    )


def weighted_moment(
    coefficients: list[LogCoeff],
    cutoff: int,
    carriers: dict[int, ComplexQ] | None = None,
) -> LogMoment:
    result: LogMoment = {}
    for number in range(1, cutoff + 1):
        value = (
            carriers.get(number, (Fraction(0), Fraction(0)))
            if carriers is not None
            else deterministic_carrier(number)
        )
        if value == (Fraction(0), Fraction(0)):
            continue
        for prime, exponent in coefficients[number].items():
            contribution = scale(value, exponent)
            result[prime] = add(
                result.get(prime, (Fraction(0), Fraction(0))),
                contribution,
            )
            if result[prime] == (Fraction(0), Fraction(0)):
                del result[prime]
    return result


def coefficient_entry_moment(
    coefficient: LogCoeff,
    value: ComplexQ,
) -> LogMoment:
    return {
        prime: scale(value, exponent)
        for prime, exponent in coefficient.items()
        if scale(value, exponent) != (Fraction(0), Fraction(0))
    }


def integer_cube_root(number: int) -> int:
    root = max(0, int(round(number ** (1 / 3))))
    while (root + 1) ** 3 <= number:
        root += 1
    while root**3 > number:
        root -= 1
    return root


def transition_audit(
    first_cutoff: int = 2,
    last_cutoff: int = 215,
) -> list[dict]:
    maximum = last_cutoff + 1
    roots = {
        max(1, integer_cube_root(number))
        for number in range(first_cutoff, maximum + 1)
    }
    caches = {
        root: vaughan_coefficients(maximum, root, root)
        for root in sorted(roots)
    }
    rows: list[dict] = []
    components = ("low", "type_i", "type_ii")

    for cutoff in range(first_cutoff, last_cutoff + 1):
        next_cutoff = cutoff + 1
        root_before = max(1, integer_cube_root(cutoff))
        root_after = max(1, integer_cube_root(next_cutoff))
        old_coefficients = caches[root_before]
        new_coefficients = caches[root_after]

        old_moments = {
            component: weighted_moment(
                old_coefficients[component],
                cutoff,
            )
            for component in components
        }
        new_moments = {
            component: weighted_moment(
                new_coefficients[component],
                next_cutoff,
            )
            for component in components
        }
        deltas = {
            component: moment_subtract(
                new_moments[component],
                old_moments[component],
            )
            for component in components
        }
        total_delta: LogMoment = {}
        for component in components:
            total_delta = moment_add(total_delta, deltas[component])
        expected = coefficient_entry_moment(
            factor_vector(next_cutoff),
            deterministic_carrier(next_cutoff),
        )
        if total_delta != expected:
            raise RuntimeError(
                f"total Vaughan transition failed at N={cutoff}"
            )

        cube_transition = root_after != root_before
        transfers: dict[str, LogMoment] = {}
        entries: dict[str, LogMoment] = {}
        for component in components:
            transfer_coefficients = [{} for _ in range(maximum + 1)]
            if cube_transition:
                for number in range(1, cutoff + 1):
                    transfer_coefficients[number] = coeff_add(
                        new_coefficients[component][number],
                        old_coefficients[component][number],
                        -1,
                    )
            transfers[component] = weighted_moment(
                transfer_coefficients,
                cutoff,
            )
            entries[component] = coefficient_entry_moment(
                new_coefficients[component][next_cutoff],
                deterministic_carrier(next_cutoff),
            )
            if deltas[component] != moment_add(
                transfers[component],
                entries[component],
            ):
                raise RuntimeError(
                    f"component transition failed at N={cutoff}: "
                    f"{component}"
                )

        total_transfer: LogMoment = {}
        for component in components:
            total_transfer = moment_add(
                total_transfer,
                transfers[component],
            )
        if total_transfer:
            raise RuntimeError(
                f"internal cube transfer failed to cancel at N={cutoff}"
            )

        if not cube_transition and any(transfers.values()):
            raise RuntimeError(
                f"ordinary transition has an internal transfer at N={cutoff}"
            )

        square_transition = math.isqrt(next_cutoff) ** 2 == next_cutoff
        rows.append(
            {
                "N": cutoff,
                "n": next_cutoff,
                "K_before": root_before,
                "K_after": root_after,
                "transition": (
                    "perfect_cube" if cube_transition else "ordinary"
                ),
                "square_transition": square_transition,
                "sixth_power_overlap": cube_transition
                and square_transition,
                "delta_low_coordinates": len(deltas["low"]),
                "delta_type_i_coordinates": len(deltas["type_i"]),
                "delta_type_ii_coordinates": len(deltas["type_ii"]),
                "transfer_low_coordinates": len(transfers["low"]),
                "transfer_type_i_coordinates": len(
                    transfers["type_i"]
                ),
                "transfer_type_ii_coordinates": len(
                    transfers["type_ii"]
                ),
                "boundary_coordinates": len(expected),
                "mismatches": 0,
            }
        )
    return rows


def sum_carriers(carriers: dict[int, ComplexQ]) -> ComplexQ:
    total = (Fraction(0), Fraction(0))
    for value in carriers.values():
        total = add(total, value)
    return total


def prime_power_chain_countermodel() -> dict:
    cutoff = 64
    split = 4
    coefficients = vaughan_coefficients(cutoff, split, split)
    carriers = {
        1: (Fraction(1), Fraction(0)),
        5: (Fraction(-2), Fraction(0)),
        25: (Fraction(1), Fraction(0)),
    }
    low = weighted_moment(coefficients["low"], cutoff, carriers)
    type_i = weighted_moment(
        coefficients["type_i"],
        cutoff,
        carriers,
    )
    type_ii = weighted_moment(
        coefficients["type_ii"],
        cutoff,
        carriers,
    )
    total = moment_add(moment_add(low, type_i), type_ii)
    expected_type_i = {5: (Fraction(1), Fraction(0))}
    expected_type_ii = {5: (Fraction(-1), Fraction(0))}
    if (
        sum_carriers(carriers) != (Fraction(0), Fraction(0))
        or low
        or type_i != expected_type_i
        or type_ii != expected_type_ii
        or total
    ):
        raise RuntimeError("prime-power Vaughan countermodel failed")

    s_prime = (Fraction(1, 16), Fraction(-1, 2))
    projected_i = real_projection(moment_multiply(type_i, s_prime))
    projected_ii = real_projection(moment_multiply(type_ii, s_prime))
    penalty = {
        prime: abs(projected_i.get(prime, Fraction(0)))
        + abs(projected_ii.get(prime, Fraction(0)))
        for prime in set(projected_i) | set(projected_ii)
    }
    if penalty != {5: Fraction(1, 8)}:
        raise RuntimeError("prime-power absolute-loss witness failed")

    return {
        "N": cutoff,
        "U": split,
        "V": split,
        "carriers": {
            str(number): complex_text(value)
            for number, value in sorted(carriers.items())
        },
        "s_prime": complex_text(s_prime),
        "endpoint_e": complex_text((Fraction(0), Fraction(0))),
        "endpoint_g": complex_text((Fraction(0), Fraction(0))),
        "W_0": complex_text(sum_carriers(carriers)),
        "Z_low": moment_text(low),
        "Z_type_i": moment_text(type_i),
        "Z_type_ii": moment_text(type_ii),
        "Z_Lambda": moment_text(total),
        "Re_sZ_type_i": real_log_text(projected_i),
        "Re_sZ_type_ii": real_log_text(projected_ii),
        "separate_absolute_penalty": real_log_text(penalty),
        "Delta_endpoint_dominance": {
            "5": fraction_text(Fraction(-1, 8))
        },
        "interpretation": (
            "At W_0=0 with zero endpoint, the nonzero Type-I and Type-II "
            "pieces are +log(5) and -log(5). Their signed sum vanishes, "
            "while separate absolute estimates pay log(5)/8 after real "
            "projection by s_*'. This is an exact carrier guard, not an "
            "actual Xi state."
        ),
    }


def prime_edge_endpoint_countermodel() -> dict:
    cutoff = 11
    split = 2
    coefficients = vaughan_coefficients(cutoff, split, split)
    carriers = {
        1: (Fraction(-1), Fraction(0)),
        11: (Fraction(1), Fraction(0)),
    }
    low = weighted_moment(coefficients["low"], cutoff, carriers)
    type_i = weighted_moment(
        coefficients["type_i"],
        cutoff,
        carriers,
    )
    type_ii = weighted_moment(
        coefficients["type_ii"],
        cutoff,
        carriers,
    )
    total = moment_add(moment_add(low, type_i), type_ii)
    if (
        sum_carriers(carriers) != (Fraction(0), Fraction(0))
        or low
        or type_i != {11: (Fraction(1), Fraction(0))}
        or type_ii
        or total != type_i
        or not (11 > math.isqrt(cutoff))
    ):
        raise RuntimeError("prime-edge endpoint countermodel failed")

    s_prime = (Fraction(1, 16), Fraction(-1, 2))
    endpoint_g = moment_multiply(total, s_prime)
    mathfrak_b = moment_subtract(
        endpoint_g,
        moment_multiply(total, s_prime),
    )
    if mathfrak_b:
        raise RuntimeError("prime-edge endpoint cancellation failed")

    return {
        "N": cutoff,
        "U": split,
        "V": split,
        "prime_edge": 11,
        "prime_edge_strict": True,
        "carriers": {
            str(number): complex_text(value)
            for number, value in sorted(carriers.items())
        },
        "s_prime": complex_text(s_prime),
        "endpoint_e": complex_text((Fraction(0), Fraction(0))),
        "endpoint_g_log_coordinates": moment_text(endpoint_g),
        "W_0": complex_text(sum_carriers(carriers)),
        "Z_low": moment_text(low),
        "Z_type_i": moment_text(type_i),
        "Z_type_ii": moment_text(type_ii),
        "Z_Lambda": moment_text(total),
        "mathfrak_B_log_coordinates": moment_text(mathfrak_b),
        "interpretation": (
            "At W_0=0 and H_a=0, the formal endpoint slope "
            "g=s_*'log(11) cancels the prime-edge Type-I contribution "
            "exactly. This proves that endpoint phase correlation cannot "
            "be omitted. The freely assigned endpoint slope makes this a "
            "route guard, not an actual Xi counterexample."
        ),
    }


def source_audit() -> dict:
    hashes: dict[str, str] = {}
    payloads: dict[str, str] = {}
    for key, path in SOURCES.items():
        if not path.is_file():
            raise RuntimeError(f"missing source {key}: {path}")
        hashes[key] = hashlib.sha256(path.read_bytes()).hexdigest()
        payloads[key] = json.dumps(
            json.loads(path.read_text(encoding="utf-8")),
            sort_keys=True,
        )
    markers = {
        "contact_centering": (
            "mathfrak B_N",
            "Z_Lambda",
            "balanced_hyperbola",
        ),
        "cutoff_transport": (
            "endpoint-value jump cancels exactly",
            "perfect_square_transitions",
        ),
        "corrected_mangoldt": (
            "symmetric Mangoldt",
            "indefinite_prime_edge_minors",
        ),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in payloads[key]:
                raise RuntimeError(
                    f"source marker missing: {key}: {marker}"
                )
    return {"source_sha256": hashes}


def exact_statements() -> dict[str, str]:
    return {
        "convolution_identity": (
            "Let lambda_0=Lambda*1_(n<=U), lambda_1=Lambda-lambda_0, "
            "mu_0=mu*1_(n<=V), and mu_1=mu-mu_0. Since "
            "Lambda=mu*log, log=Lambda*1, and mu*1=delta, "
            "Lambda=lambda_0+mu_0*log-lambda_0*mu_0*1"
            "+mu_1*lambda_1*1. Therefore "
            "log=(lambda_0*1)+(mu_0*log*1"
            "-lambda_0*mu_0*1*1)+(mu_1*lambda_1*1*1)."
        ),
        "weighted_vaughan": (
            "Define C_0=lambda_0*1, "
            "C_I=mu_0*log*1-lambda_0*mu_0*1*1, and "
            "C_II=mu_1*lambda_1*1*1. Then C_0+C_I+C_II=log "
            "coefficientwise and "
            "Z_Lambda=Z_0+Z_I+Z_II, where "
            "Z_j=sum_(n<=N)C_j(n)z_n."
        ),
        "endpoint_composition": (
            "With R_(0;U,V)=E_N+s_*'log(a)W_0-s_*'Z_0, "
            "mathfrak B_N=R_(0;U,V)-s_*'Z_I-s_*'Z_II exactly. "
            "The endpoint and W_0 term are composed before any real "
            "projection or absolute value."
        ),
        "contact_fibres": (
            "At W_0=0, R_(0;U,V)=E_N-s_*'Z_0. At mathsf_X=0, "
            "mathcal C_N=Re[R_(0;U,V)-s_*'Z_I-s_*'Z_II]. "
            "No division by W_0, H_a, or a Vaughan component is used."
        ),
        "dominance_criterion": (
            "Put R_0=Re R_(0;U,V), R_I=Re(s_*'Z_I), and "
            "R_II=Re(s_*'Z_II). The reverse triangle inequality gives "
            "|Re mathfrak B_N|>=|R_0|-|R_I|-|R_II|. Hence "
            "|R_0|-|R_I|-|R_II|>A_L+2epsilon_term is sufficient for "
            "the Abel gap, but no such uniform inequality is proved."
        ),
        "joint_signed_requirement": (
            "Squaring the exact real scalar retains "
            "(R_0-R_I-R_II)^2=R_0^2+R_I^2+R_II^2"
            "-2R_0R_I-2R_0R_II+2R_IR_II. The three cross terms "
            "are part of the theorem. Separate Type-I/II norms or "
            "post-decomposition absolute values do not control them."
        ),
        "canonical_split": (
            "The canonical finite choice U=V=K_N=floor(N^(1/3)) "
            "is an exact cubic-balanced coordinate. It creates no "
            "estimate and no sign by itself."
        ),
        "ordinary_transport": (
            "If N+1 is not a perfect cube, K_(N+1)=K_N and each "
            "Vaughan component gains only C_j(N+1)z_(N+1). Their "
            "sum is log(N+1)z_(N+1)."
        ),
        "cube_transport": (
            "If N+1=r^3, K changes from r-1 to r. For m<=N put "
            "tau_j(m)=C_j^(r)(m)-C_j^(r-1)(m). Then "
            "sum_j tau_j(m)=0 coefficientwise, so all historical "
            "parameter transfers cancel. The entering coefficients "
            "still sum to log(r^3)."
        ),
        "overlap_guard": (
            "At n=64 the balanced square/wing boundary and the cubic "
            "Vaughan parameter both change. The square-row transfer and "
            "the Vaughan coefficient transfer are separate internal "
            "reindexings; both recombine to the same physical "
            "log(64)z_64 jump."
        ),
        "countermodel_consequence": (
            "The exact 1,5,25 chain has W_0=0, zero endpoint, "
            "Z_I=log(5), Z_II=-log(5), and Z_Lambda=0. Thus signed "
            "Type-I/II cancellation can be complete while separate "
            "absolute projections pay log(5)/8."
        ),
        "prime_edge_consequence": (
            "The exact 1,11 prime-edge model has W_0=0 and "
            "Z_I=log(11). Choosing the allowed formal endpoint slope "
            "g=s_*'log(11) gives mathfrak B_N=0. A proof must use the "
            "actual Xi coupling of the endpoint to the carriers."
        ),
        "theorem_fit": (
            "The new decomposition identifies a precise missing input: "
            "an Xi-specific pointwise signed correlation theorem for the "
            "endpoint core and the joint Type-I/II pair on the contact "
            "band. Mean-square Vaughan estimates, componentwise upper "
            "bounds, and generic amplitude ordering do not imply it."
        ),
        "route_decision": (
            "Retain the cubic-balanced Vaughan split as a proof-facing "
            "coordinate and falsification harness. The next candidate "
            "must state a joint signed endpoint/Type-I/Type-II inequality "
            "for the actual logarithmic carriers and survive W_0=0, "
            "mathsf_X=0, H_a=0, q=1, prime edges, cube transfers, the "
            "n=64 square/cube overlap, and adjacent real projection."
        ),
    }


def build_rows(exact: dict[str, str]) -> list[GateRow]:
    return [
        GateRow(
            "mecvg_01_convolution",
            "The physical logarithm has an exact finite Vaughan decomposition.",
            "proved",
            exact["convolution_identity"],
            "The outer divisor convolution is retained.",
            "This is arithmetic algebra, not cancellation.",
        ),
        GateRow(
            "mecvg_02_weighted",
            "The physical Mangoldt moment splits into low, Type-I, and Type-II pieces.",
            "proved",
            exact["weighted_vaughan"],
            "Every corrected carrier remains inside the weighted sum.",
            "No component has a sign.",
        ),
        GateRow(
            "mecvg_03_endpoint",
            "The Vaughan pieces compose exactly with the centered endpoint.",
            "proved",
            exact["endpoint_composition"],
            "The proof-facing object remains mathfrak B_N.",
            "Z_Lambda is not estimated in isolation.",
        ),
        GateRow(
            "mecvg_04_fibres",
            "The endpoint-composed split remains division-free on both exceptional fibres.",
            "proved",
            exact["contact_fibres"],
            "W_0=0 and mathsf_X=0 remain explicit tests.",
            "No transversality is inferred.",
        ),
        GateRow(
            "mecvg_05_criterion",
            "A componentwise reverse-triangle criterion is exactly sufficient.",
            "proved_sufficient",
            exact["dominance_criterion"],
            "It is a falsifiable stronger theorem.",
            "Its uniform Xi validity remains open.",
        ),
        GateRow(
            "mecvg_06_joint",
            "Any quadratic proof must retain all endpoint/Type-I/Type-II cross terms.",
            "proved",
            exact["joint_signed_requirement"],
            "The missing theorem is a signed joint estimate.",
            "Separate norm bounds are not equivalent.",
        ),
        GateRow(
            "mecvg_07_cubic",
            "The floor-cube-root choice is a canonical finite coordinate.",
            "proved",
            exact["canonical_split"],
            "It gives a deterministic Type-I/II registry.",
            "The choice supplies no power gain.",
        ),
        GateRow(
            "mecvg_08_ordinary",
            "Ordinary cutoff transport is componentwise exact.",
            "proved",
            exact["ordinary_transport"],
            "The physical entering logarithm is recovered.",
            "Only fixed-parameter insertion occurs.",
        ),
        GateRow(
            "mecvg_09_cube",
            "Perfect-cube parameter transfers cancel coefficientwise.",
            "proved",
            exact["cube_transport"],
            "No cube-root floor residual survives.",
            "Individual component transfers can be large.",
        ),
        GateRow(
            "mecvg_10_overlap",
            "The square/cube overlap at 64 has compatible internal transfers.",
            "finite_certificate",
            exact["overlap_guard"],
            "Both decompositions represent one physical jump.",
            "This is not a signed estimate.",
        ),
        GateRow(
            "mecvg_11_chain_guard",
            "Separate Type-I and Type-II absolute estimates destroy an exact chain cancellation.",
            "guard_validated",
            exact["countermodel_consequence"],
            "Signed joint control is necessary.",
            "The carrier model is not an actual Xi state.",
        ),
        GateRow(
            "mecvg_12_edge_guard",
            "Prime-edge control requires actual endpoint/carrier correlation.",
            "guard_validated",
            exact["prime_edge_consequence"],
            "Endpoint phase cannot be discarded.",
            "The freely assigned endpoint is not an Xi counterexample.",
        ),
        GateRow(
            "mecvg_13_fit",
            "The surviving analytic input is an Xi-specific signed correlation theorem.",
            "open",
            exact["theorem_fit"],
            "This is narrower than generic Vaughan bounds.",
            "No such pointwise theorem is currently proved.",
        ),
        GateRow(
            "mecvg_14_route",
            "The next candidate must preserve the complete signed combination.",
            "open",
            exact["route_decision"],
            "The gate supplies its exact test matrix.",
            (
                "No signed lower bound, Abel gap, winding cap, contact "
                "exclusion, Lambda<=0, PF-infinity, RH, or prize-level "
                "conclusion is proved."
            ),
        ),
    ]


def build_note(payload: dict) -> str:
    exact = payload["exact"]
    summary = payload["summary"]
    return f"""# Centered Mangoldt Endpoint-Composed Vaughan Gate

Date: 2026-07-31

Status: exact Vaughan handoff and cancellation-guard artifact; not a proof
of a signed centered-jet lower bound, Abel gap, `Lambda<=0`, PF-infinity,
RH, or a Clay-prize result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Why The Extra Convolution Matters

The physical observable is `sum_(n<=N)log(n)z_n`, not
`sum_(n<=N)Lambda(n)z_n`.  Therefore Vaughan is applied to the Mangoldt
divisor in `log=Lambda*1`, and the final convolution by `1` must remain.

```text
{exact["convolution_identity"]}
```

Put

```text
{exact["weighted_vaughan"]}
```

This is an exact coefficient identity for every finite `N,U,V`.

## Endpoint Composition

```text
{exact["endpoint_composition"]}

{exact["contact_fibres"]}
```

The decomposition has not changed the live scalar or deleted either
exceptional fibre.

## First Stronger Candidate

```text
{exact["dominance_criterion"]}
```

This criterion is rigorously sufficient, but it is deliberately stronger
than the target.  Squaring does not permit the pieces to be separated:

```text
{exact["joint_signed_requirement"]}
```

## Canonical Cubic Split

```text
{exact["canonical_split"]}

{exact["ordinary_transport"]}

{exact["cube_transport"]}
```

The exact audit checks `{summary["cutoff_transitions"]}` transitions:
`{summary["ordinary_transitions"]}` ordinary and
`{summary["perfect_cube_transitions"]}` perfect-cube transitions.  It
also contains `{summary["sixth_power_overlaps"]}` simultaneous
square/cube overlap, at `n=64`, and has
`{summary["transition_mismatches"]}` coefficient mismatches.

## Exact Cancellation Guards

```text
{exact["countermodel_consequence"]}
```

The countermodel uses `s_*'=1/16-i/2`.  The two real projections are
`+log(5)/16` and `-log(5)/16`; the exact signed sum is zero, while separate
absolute values cost `log(5)/8`.

```text
{exact["prime_edge_consequence"]}
```

These are route countermodels.  They prove no physical Xi failure.

## Surviving Theorem

```text
{exact["theorem_fit"]}

{exact["route_decision"]}
```

The earlier Mobius Vaughan gates concern a different bridge problem and
do not supply this pointwise signed theorem.

## Pi Provenance

No new `pi` appears in the Vaughan identity.  The only `pi` retained by
the endpoint-complete Xi coordinate remains inherited from the
completed-zeta and Riemann-Siegel normalization.

## Boundary

This gate proves the finite Vaughan convolution identity for the physical
logarithmic moment, endpoint composition, exceptional-fibre formulas,
reverse-triangle sufficient criterion, cubic-root cutoff transport, and
two exact cancellation guards.  It proves no pointwise Xi-specific signed
correlation estimate, signed centered-jet lower bound, Abel-scalar gap,
horizontal successor winding cap, contact exclusion, Q209, cofinal
descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level result.
"""


def build_payload() -> dict:
    exact = exact_statements()
    transitions = transition_audit()
    rows = build_rows(exact)
    return {
        "kind": "mangoldt_endpoint_composed_vaughan_gate",
        "date": "2026-07-31",
        "status": (
            "exact endpoint-composed Vaughan handoff and cancellation "
            "guards; signed lower bound remains open"
        ),
        "source_audit": source_audit(),
        "exact": exact,
        "witnesses": {
            "prime_power_chain": prime_power_chain_countermodel(),
            "prime_edge_endpoint": prime_edge_endpoint_countermodel(),
            "cutoff_transitions": transitions,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "convolution_identities": 1,
            "vaughan_components": 3,
            "cutoff_transitions": len(transitions),
            "ordinary_transitions": sum(
                row["transition"] == "ordinary" for row in transitions
            ),
            "perfect_cube_transitions": sum(
                row["transition"] == "perfect_cube"
                for row in transitions
            ),
            "square_transitions": sum(
                bool(row["square_transition"]) for row in transitions
            ),
            "sixth_power_overlaps": sum(
                bool(row["sixth_power_overlap"]) for row in transitions
            ),
            "transition_mismatches": sum(
                int(row["mismatches"]) for row in transitions
            ),
            "exact_countermodels": 2,
            "proved_componentwise_dominance_theorems": 0,
            "proved_signed_lower_bounds": 0,
            "proved_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
        "proof_boundary": (
            "This gate proves an exact finite Vaughan decomposition of "
            "the physical logarithmic moment, its endpoint composition, "
            "ordinary and perfect-cube transport, a reverse-triangle "
            "sufficient criterion, and two exact route guards. It proves "
            "no Xi-specific joint signed correlation estimate, signed "
            "lower bound, Abel gap, winding cap, contact exclusion, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    args = parser.parse_args()

    payload = build_payload()
    atomic_write(
        args.result,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    atomic_write(args.note, build_note(payload))
    summary = payload["summary"]
    print(
        "built centered Mangoldt endpoint-composed Vaughan gate: "
        f"{summary['rows']} rows, "
        f"{summary['cutoff_transitions']} transitions, "
        f"{summary['perfect_cube_transitions']} cube transfers, "
        f"{summary['exact_countermodels']} exact countermodels"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
