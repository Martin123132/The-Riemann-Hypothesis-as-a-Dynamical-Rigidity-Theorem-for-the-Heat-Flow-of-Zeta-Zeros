#!/usr/bin/env python3
"""Independently check the complementary root-of-unity endpoint sums."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
    "nonstationary_complement_root_unity_hurwitz_endpoint_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

MAX_POWER = 16
CHECK_POWER = 18


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def root_of_unity(order: int, numerator: int) -> acb:
    return (acb(0, 1) * (2 * arb.pi() * numerator / order)).exp()


def phase(q0: arb, y: arb) -> acb:
    return (acb(0, 2 * arb.pi() * q0 * y)).exp()


def saved_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def contains_zero(value: acb) -> bool:
    return value.real.contains(0) and value.imag.contains(0)


def direct_finite(
    delta: arb,
    order: int,
    numerator: int,
    count: int,
    power: int,
) -> acb:
    omega = root_of_unity(order, numerator)
    return sum(
        (omega**index / (delta + index) ** power for index in range(count)),
        acb(0),
    )


def residue_finite(
    delta: arb,
    order: int,
    numerator: int,
    cycles: int,
    power: int,
) -> acb:
    omega = root_of_unity(order, numerator)
    h = arb(order)
    blocks: list[acb] = []
    for residue in reversed(range(order)):
        argument = (delta + residue) / h
        if power == 1:
            block = (argument + cycles).digamma() - argument.digamma()
        else:
            exponent = arb(power)
            block = exponent.zeta(argument) - exponent.zeta(argument + cycles)
        blocks.append(omega**residue * block)
    return sum(blocks, acb(0)) / h**power


def eta(delta: arb, power: int) -> arb:
    if power == 1:
        return (((delta + 1) / 2).digamma() - (delta / 2).digamma()) / 2
    exponent = arb(power)
    return (
        exponent.zeta(delta / 2) - exponent.zeta((delta + 1) / 2)
    ) / arb(2) ** power


def alternating_direct_ball(delta: arb, power: int, count: int) -> arb:
    require(count % 2 == 0, "even alternating witness length required")
    partial = arb(0)
    for index in range(count):
        term = (delta + index) ** (-power)
        partial += term if index % 2 == 0 else -term
    next_term = (delta + count) ** (-power)
    return arb(partial, next_term)


def verify_saved_row(
    row: dict[str, Any],
    formula: acb,
    weighted: acb,
) -> None:
    saved_formula = saved_complex(row["unphased_formula_ball"])
    saved_weighted = saved_complex(row["phase_weighted_sum_ball"])
    require(formula.overlaps(saved_formula), f"saved formula drift at k={row['power']}")
    require(weighted.overlaps(saved_weighted), f"saved weighted sum drift at k={row['power']}")


def check_finite_endpoint(endpoint: dict[str, Any]) -> None:
    y = arb(endpoint["endpoint_y"])
    q0 = arb(endpoint["base_label_q0"])
    delta = arb(endpoint["positive_shift_delta_ball"]["ball"])
    order = endpoint["root_order"]
    numerator = endpoint["root_numerator"]
    cycles = endpoint["cycles"]
    count = order * cycles
    require(len(endpoint["rows"]) == MAX_POWER, "finite endpoint row count drift")
    for row in endpoint["rows"]:
        power = row["power"]
        formula = residue_finite(delta, order, numerator, cycles, power)
        direct = direct_finite(delta, order, numerator, count, power)
        require(formula.overlaps(direct), f"independent finite mismatch at k={power}")
        require(contains_zero(formula - direct), f"independent finite zero miss at k={power}")
        verify_saved_row(row, formula, phase(q0, y) * formula)


def check_infinite_endpoint(endpoint: dict[str, Any], lower_orientation: bool) -> None:
    y = arb(endpoint["endpoint_y"])
    q0 = arb(endpoint["base_label_q0"])
    delta = arb(endpoint["positive_shift_delta_ball"]["ball"])
    require(endpoint["root_order"] == 2, "production infinite root order drift")
    require(len(endpoint["rows"]) == MAX_POWER, "infinite endpoint row count drift")
    for row in endpoint["rows"]:
        power = row["power"]
        formula = acb(eta(delta, power))
        recurrence = eta(delta, power) + eta(delta + 1, power) - delta ** (-power)
        require(recurrence.contains(0), f"independent alternating recurrence miss at k={power}")
        orientation = -1 if lower_orientation and power % 2 else 1
        verify_saved_row(row, formula, phase(q0, y) * orientation * formula)


def altered_finite_witness() -> acb:
    delta = arb("2.375")
    formula = residue_finite(delta, 5, 2, 5, CHECK_POWER)
    direct = direct_finite(delta, 5, 2, 25, CHECK_POWER)
    discrepancy = formula - direct
    require(formula.overlaps(direct), "altered order-5 finite witness mismatch")
    require(contains_zero(discrepancy), "altered order-5 discrepancy excludes zero")
    return discrepancy


def altered_alternating_witness() -> arb:
    delta = arb("3.125")
    widest = arb(0)
    for power in range(1, CHECK_POWER + 1):
        formula = eta(delta, power)
        direct_ball = alternating_direct_ball(delta, power, 8192)
        require(formula.overlaps(direct_ball), f"altered alternating witness miss at k={power}")
        difference = formula - direct_ball
        require(difference.contains(0), f"altered alternating difference excludes zero at k={power}")
        widest = max(widest, abs(difference).upper())
    return widest


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 448
    ctx.threads = 1
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing artifact: {path}")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production gate not passed")
    require(artifact["certificate"]["maximum_denominator_power"] == MAX_POWER, "power drift")
    require(artifact["decision"]["integrated_remainders_bounded"] is False, "remainder overpromotion")
    require(artifact["decision"]["complete_joined_packet_enclosed"] is False, "packet overpromotion")

    for source in artifact["sources"].values():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], "source hash drift")
    for dependency in artifact["dependencies"].values():
        path = ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")

    families = artifact["certificate"]["families"]
    finite = families["R_B_finite_240"]
    check_finite_endpoint(finite["endpoint_B"])
    check_finite_endpoint(finite["endpoint_a"])
    upper = families["R_L_infinite_upper"]
    check_infinite_endpoint(upper["endpoint_L"], False)
    check_infinite_endpoint(upper["endpoint_a"], False)
    lower = families["R_lower_infinite"]
    check_infinite_endpoint(lower["endpoint_L"], True)
    check_infinite_endpoint(lower["endpoint_a"], True)

    finite_discrepancy = altered_finite_witness()
    alternating_width = altered_alternating_witness()

    text = " ".join(NOTE.read_text(encoding="utf-8").split())
    for fragment in (
        "sum_(r=0)^15 exp(2*pi*i*9*r/16)=0",
        "F_inf(k,delta,-1)+F_inf(k,delta+1,-1)=delta^(-k)",
        "does not bound the integral containing `h_8`",
        "it is not a fitted circle constant",
    ):
        require(fragment in text, f"missing note fragment: {fragment}")

    print(
        "independently checked root-of-unity Hurwitz endpoint sums; "
        f"altered finite discrepancy {abs(finite_discrepancy).str(8, more=True)}, "
        f"alternating witness width {alternating_width.str(8, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
