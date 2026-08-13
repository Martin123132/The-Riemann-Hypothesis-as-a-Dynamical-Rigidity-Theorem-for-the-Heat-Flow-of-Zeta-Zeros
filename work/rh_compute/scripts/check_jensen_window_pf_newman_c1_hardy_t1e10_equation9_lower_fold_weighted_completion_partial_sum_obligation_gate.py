#!/usr/bin/env python3
"""Independently check the weighted-completion partial-sum obligation."""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
L = 39_696
Q = 39_894
U = 40_094


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_ball(text: str) -> tuple[Decimal, Decimal]:
    body = text.strip().removeprefix("[").removesuffix("]").strip()
    if "+/-" not in body:
        return Decimal(body), Decimal(0)
    midpoint, radius = body.split("+/-", 1)
    return Decimal(midpoint.strip()), Decimal(radius.strip())


def add_ball(left: tuple[Decimal, Decimal], right: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return left[0] + right[0], left[1] + right[1]


def sub_ball(left: tuple[Decimal, Decimal], right: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return left[0] - right[0], left[1] + right[1]


def conjugate_ball(value: tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]) -> tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]:
    return value[0], (-value[1][0], value[1][1])


def add_complex(left: tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]], right: tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]) -> tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]:
    return add_ball(left[0], right[0]), add_ball(left[1], right[1])


def sub_complex(left: tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]], right: tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]) -> tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]:
    return sub_ball(left[0], right[0]), sub_ball(left[1], right[1])


def parse_complex(record: dict[str, str]) -> tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]:
    return parse_ball(record["real_ball"]), parse_ball(record["imag_ball"])


def absolute_upper(value: tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]) -> Decimal:
    real = abs(value[0][0]) + value[0][1]
    imag = abs(value[1][0]) + value[1][1]
    return (real * real + imag * imag).sqrt()


def reconstruct_weights(pairing: dict) -> dict[int, tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]]:
    weights = {}
    certificate = pairing["leading_multiplier_certificate"]
    for row in certificate["pair_rows"]:
        coherent = parse_complex(row["coherent_coefficient_ball"])
        leakage = parse_complex(row["conjugacy_leakage_coefficient_ball"])
        weights[int(row["minus_mode"])] = add_complex(coherent, leakage)
        weights[int(row["plus_mode"])] = conjugate_ball(sub_complex(coherent, leakage))
    for row in certificate["edge_rows"]:
        weights[int(row["mode"])] = parse_complex(row["leading_defect_ball"])
    weights[Q] = ((Decimal(0), Decimal(0)), (Decimal(0), Decimal(0)))
    require(set(weights) == set(range(L, U + 1)), "independent roster reconstruction failed")
    return weights


def rational_identity_check() -> None:
    weights = {m: Fraction(2 * (m - L) + 3, 5 * (m - L) + 11) for m in range(L, U + 1)}
    weights[Q] = Fraction(0)
    terms = {m: Fraction(7 * (m - L) - 13, 17 * (m - L + 1) + 19) for m in range(L, U + 1)}
    h4, c4 = Fraction(5, 29), Fraction(-3, 31)
    g4 = h4 + c4 + sum(terms.values(), Fraction())
    direct = sum((weights[m] * terms[m] for m in range(L, U + 1)), Fraction())
    prefixes = {k: sum((terms[m] for m in range(L, k + 1)), Fraction()) for k in range(L, U + 1)}
    remainders = {
        k: h4 + c4 + sum((terms[m] for m in range(k + 1, U + 1)), Fraction())
        for k in range(L, U + 1)
    }
    abel = weights[U] * prefixes[U] + sum(
        ((weights[k] - weights[k + 1]) * prefixes[k] for k in range(L, U)), Fraction()
    )
    completed = weights[L] * g4 - weights[U] * remainders[U] - sum(
        ((weights[k] - weights[k + 1]) * remainders[k] for k in range(L, U)), Fraction()
    )
    require(direct == abel == completed, "independent completed Abel identity failed")

    z = Fraction(37, 41)
    perturbed = dict(terms)
    perturbed[U - 1] += z
    perturbed[U] -= z
    require(sum(perturbed.values(), Fraction()) == sum(terms.values(), Fraction()), "countermodel changed total")
    require(perturbed[Q] == terms[Q], "countermodel changed event zero")
    shifted = sum((weights[m] * perturbed[m] for m in range(L, U + 1)), Fraction()) - direct
    require(shifted == (weights[U - 1] - weights[U]) * z, "countermodel weighted shift failed")


def main() -> None:
    getcontext().prec = 80
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")

    rational_identity_check()
    pairing_path = REPO_ROOT / artifact["dependencies"]["event_pairing"]["path"]
    pairing = json.loads(pairing_path.read_text(encoding="utf-8"))
    weights = reconstruct_weights(pairing)
    edge_difference = sub_complex(weights[U - 1], weights[U])
    require(abs(edge_difference[1][0]) > edge_difference[1][1], "last-edge imaginary intervals overlap zero")
    total_variation = sum((absolute_upper(sub_complex(weights[k + 1], weights[k])) for k in range(L, U)), Decimal(0))
    augmented = absolute_upper(weights[L]) + total_variation + absolute_upper(weights[U])
    require(total_variation < Decimal("1.547052e-5"), "independent total variation failed")
    require(augmented < Decimal("3.093180e-5"), "independent augmented variation failed")

    certificate = artifact["coefficient_certificate"]
    stored_augmented = parse_ball(certificate["endpoint_augmented_total_variation_ball"])
    require(stored_augmented[0] + stored_augmented[1] < Decimal("3.093180e-5"), "stored variation bound failed")
    decision = artifact["decision"]
    require(decision["global_unweighted_completion_alone_identifies_weighted_roster"] is False, "identifiability guard drift")
    require(decision["completed_partial_sum_family_exactly_reduces_weighted_roster"] is True, "partial-family reduction drift")
    require(decision["uniform_completed_partial_remainder_bound_proved"] is False, "proof boundary drift")
    print("independently checked weighted-completion partial-sum obligation", flush=True)


if __name__ == "__main__":
    main()
