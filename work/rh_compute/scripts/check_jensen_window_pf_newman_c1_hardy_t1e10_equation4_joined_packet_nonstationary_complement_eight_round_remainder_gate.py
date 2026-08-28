#!/usr/bin/env python3
"""Independently check the eight-round complementary-tail enclosures."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[3]
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
    "nonstationary_complement_eight_round_remainder_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

CHECK_ORDER = 9
CHECK_SLABS = 10_000
CHECK_CLUSTER_POWER = 5
HEIGHT = arb("10000000000")
L = arb("621.5")
B = arb("621.5625")
A = arb("39852.5")
Q_PLUS = arb("2561211.5")
Q_REMOTE = Q_PLUS + 240
Q_LOW = arb("79787.5")


CoefficientTable = list[dict[int, dict[tuple[int, int], Fraction]]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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


def saved_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def add_error(value: acb, error: arb) -> acb:
    return acb(arb(value.real, error), arb(value.imag, error))


def coefficient_table(order: int) -> CoefficientTable:
    table: CoefficientTable = [{0: {(0, -1): Fraction(1)}}]
    for _ in range(order):
        next_level: dict[int, dict[tuple[int, int], Fraction]] = {}
        for power, terms in table[-1].items():
            for (p_power, exponent2), coefficient in terms.items():
                for new_power, monomial, value in (
                    (
                        power + 1,
                        (p_power, exponent2 - 2),
                        coefficient * Fraction(exponent2, 2),
                    ),
                    (power + 2, (p_power, exponent2), coefficient * (power + 1)),
                    (
                        power + 2,
                        (p_power + 1, exponent2 - 4),
                        -coefficient * (power + 1),
                    ),
                ):
                    target = next_level.setdefault(new_power, {})
                    target[monomial] = target.get(monomial, Fraction(0)) + value
        table.append(
            {
                power: {key: value for key, value in terms.items() if value}
                for power, terms in next_level.items()
            }
        )
    return table


def half_power(y: arb, exponent2: int) -> arb:
    return y.sqrt() ** exponent2


def coefficient_value(
    terms: dict[tuple[int, int], Fraction], p: arb, y: arb
) -> arb:
    return sum(
        (
            arb(coefficient.numerator)
            / coefficient.denominator
            * p**p_power
            * half_power(y, exponent2)
            for (p_power, exponent2), coefficient in terms.items()
        ),
        arb(0),
    )


def coefficient_majorant(
    terms: dict[tuple[int, int], Fraction], p: arb, y_left: arb
) -> arb:
    return sum(
        (
            abs(arb(coefficient.numerator) / coefficient.denominator)
            * p**p_power
            * half_power(y_left, exponent2)
            for (p_power, exponent2), coefficient in terms.items()
        ),
        arb(0),
    ).upper()


def q_curve(y: arb, p: arb) -> arb:
    return y + p / y


def phase(q0: arb, y: arb) -> acb:
    return (acb(0, 2 * arb.pi() * q0 * y)).exp()


def carrier(t: arb, y: arb) -> acb:
    return (acb(0, -t * y.log() - arb.pi() * y * y)).exp()


def root(order: int, numerator: int) -> acb:
    return (acb(0, 2 * arb.pi() * numerator / order)).exp()


def finite_sum(
    delta: arb,
    order: int,
    numerator: int,
    count: int,
    power: int,
) -> acb:
    omega = root(order, numerator)
    return sum(
        (omega**index / (delta + index) ** power for index in range(count)),
        acb(0),
    )


def eta(delta: arb, power: int) -> arb:
    if power == 1:
        return (((delta + 1) / 2).digamma() - (delta / 2).digamma()) / 2
    exponent = arb(power)
    return (
        exponent.zeta(delta / 2) - exponent.zeta((delta + 1) / 2)
    ) / arb(2) ** power


def production_endpoint_sum(
    family: str, y: arb, p: arb, power: int
) -> acb:
    qy = q_curve(y, p)
    if family == "R_B":
        delta = Q_PLUS - qy
        if y == B:
            return phase(Q_PLUS, y) * finite_sum(delta, 16, 9, 240, power)
        return phase(Q_PLUS, y) * finite_sum(delta, 2, 1, 240, power)
    if family == "R_L":
        delta = Q_REMOTE - qy
        return phase(Q_REMOTE, y) * eta(delta, power)
    if family == "R_lower":
        delta = qy - Q_LOW
        orientation = -1 if power % 2 else 1
        return phase(Q_LOW, y) * orientation * eta(delta, power)
    raise RuntimeError(f"unknown production family: {family}")


def endpoint_value(
    table_level: dict[int, dict[tuple[int, int], Fraction]],
    family: str,
    p: arb,
    y: arb,
) -> acb:
    value = acb(0)
    for denominator_power, terms in table_level.items():
        value += coefficient_value(terms, p, y) * production_endpoint_sum(
            family, y, p, denominator_power + 1
        )
    return carrier(HEIGHT, y) * value


def endpoint_partial(
    table: CoefficientTable,
    family: str,
    p: arb,
    left: arb,
    right: arb,
    order: int,
) -> acb:
    value = acb(0)
    scale = acb(0, 2 * arb.pi())
    for level in range(order):
        boundary = (
            endpoint_value(table[level], family, p, right)
            - endpoint_value(table[level], family, p, left)
        )
        value += (-1) ** level * scale ** (-(level + 1)) * boundary
    return value


def clustered_edge(left: arb, right: arb, index: int, side: str) -> arb:
    width = right - left
    denominator = CHECK_SLABS**CHECK_CLUSTER_POWER
    if side == "left":
        return left + width * arb(index**CHECK_CLUSTER_POWER) / denominator
    return right - width * arb((CHECK_SLABS - index) ** CHECK_CLUSTER_POWER) / denominator


def series_bound(delta: arb, power: int, finite_count: int | None) -> arb:
    infinite = (delta ** (-power) + delta ** (1 - power) / (power - 1)).upper()
    if finite_count is None:
        return infinite
    finite = (finite_count * delta ** (-power)).upper()
    return finite if finite < infinite else infinite


def production_remainder(
    level: dict[int, dict[tuple[int, int], Fraction]],
    p: arb,
    left: arb,
    right: arb,
    q0: arb,
    orientation: str,
    finite_count: int | None,
    order: int,
) -> arb:
    side = "left" if orientation == "upper" else "right"
    total = arb(0)
    for index in range(CHECK_SLABS):
        slab_left = clustered_edge(left, right, index, side)
        slab_right = clustered_edge(left, right, index + 1, side)
        if orientation == "upper":
            delta = (q0 - q_curve(slab_left, p)).lower()
        else:
            delta = (q_curve(slab_right, p) - q0).lower()
        require(delta > 0, f"changed-order production gap lost at slab {index}")
        integrand = arb(0)
        for power, terms in level.items():
            integrand += coefficient_majorant(terms, p, slab_left) * series_bound(
                delta, power, finite_count
            )
        total += (slab_right - slab_left) * integrand
    return (total / (2 * arb.pi()) ** order).upper()


def changed_production_ball(
    table: CoefficientTable,
    family: str,
    p: arb,
) -> acb:
    if family == "R_B":
        left, q0, orientation, count = B, Q_PLUS, "upper", 240
    elif family == "R_L":
        left, q0, orientation, count = L, Q_REMOTE, "upper", None
    else:
        left, q0, orientation, count = L, Q_LOW, "lower", None
    partial = endpoint_partial(table, family, p, left, A, CHECK_ORDER)
    remainder = production_remainder(
        table[CHECK_ORDER], p, left, A, q0, orientation, count, CHECK_ORDER
    )
    return add_error(partial, remainder)


def toy_endpoint_sum(labels: list[arb], p: arb, y: arb, power: int) -> acb:
    qy = q_curve(y, p)
    return sum(
        (
            (acb(0, 2 * arb.pi() * q * y)).exp() / (q - qy) ** power
            for q in labels
        ),
        acb(0),
    )


def toy_endpoint_value(
    level: dict[int, dict[tuple[int, int], Fraction]],
    labels: list[arb],
    p: arb,
    t: arb,
    y: arb,
) -> acb:
    value = acb(0)
    for power, terms in level.items():
        value += coefficient_value(terms, p, y) * toy_endpoint_sum(
            labels, p, y, power + 1
        )
    return (acb(0, -t * y.log() - arb.pi() * y * y)).exp() * value


def toy_remainder(
    level: dict[int, dict[tuple[int, int], Fraction]],
    labels: list[arb],
    p: arb,
    left: arb,
    right: arb,
    order: int,
) -> arb:
    slabs = 1024
    total = arb(0)
    for index in range(slabs):
        slab_left = left + (right - left) * index / slabs
        slab_right = left + (right - left) * (index + 1) / slabs
        y_ball = arb((slab_left + slab_right) / 2, (slab_right - slab_left) / 2)
        integrand = arb(0)
        for power, terms in level.items():
            denominator_sum = arb(0)
            for q in labels:
                gap = abs(q - y_ball - p / y_ball).lower()
                require(gap > 0, "altered witness acquired a stationary point")
                denominator_sum += gap ** (-power)
            integrand += coefficient_majorant(terms, p, slab_left) * denominator_sum
        total += (slab_right - slab_left) * integrand
    return (total / (2 * arb.pi()) ** order).upper()


def integrate(
    integrand: Callable[[acb, bool], acb], width: arb, panels: int
) -> acb:
    value = acb(0)
    for index in range(panels):
        left = width * index / panels
        right = width * (index + 1) / panels
        value += acb.integral(
            integrand,
            left,
            right,
            abs_tol=arb("1e-55"),
            rel_tol=arb("1e-55"),
            eval_limit=250_000,
            depth_limit=40,
        )
    return value


def direct_toy_integral(labels: list[arb], p: arb, left: arb, right: arb) -> acb:
    t = 2 * arb.pi() * p
    imaginary = acb(0, 1)
    total = acb(0)
    for q in labels:
        def integrand(x: acb, analytic: bool, q_value: arb = q) -> acb:
            y = left + x
            amplitude = 1 / y.sqrt(analytic=analytic)
            phase_value = (
                -t * y.log(analytic=analytic)
                - arb.pi() * y * y
                + 2 * arb.pi() * q_value * y
            )
            return amplitude * (imaginary * phase_value).exp()

        total += integrate(integrand, right - left, 12)
    return total


def altered_sign_witness(labels: list[arb]) -> tuple[acb, acb]:
    order = 5
    table = coefficient_table(order)
    p = arb("4.75")
    t = 2 * arb.pi() * p
    left = arb("1.5")
    right = arb("1.75")
    partial = acb(0)
    scale = acb(0, 2 * arb.pi())
    for level in range(order):
        boundary = (
            toy_endpoint_value(table[level], labels, p, t, right)
            - toy_endpoint_value(table[level], labels, p, t, left)
        )
        partial += (-1) ** level * scale ** (-(level + 1)) * boundary
    error = toy_remainder(table[order], labels, p, left, right, order)
    enclosure = add_error(partial, error)
    direct = direct_toy_integral(labels, p, left, right)
    require(enclosure.overlaps(direct), "altered signed-D witness misses direct quadrature")
    return enclosure, direct


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 448
    ctx.threads = 1
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production gate not passed")
    require(artifact["decision"]["complete_joined_packet_enclosed"] is False, "join overpromotion")
    require(
        artifact["decision"]["unsplit_eight_round_absolute_route_rejected"] is True,
        "failed route was not rejected",
    )
    require(
        artifact["decision"]["eight_round_absolute_bounds_useful_at_D_K_scale"] is False,
        "useless bounds were overpromoted",
    )
    for source in artifact["sources"].values():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], "source hash drift")
    for dependency in artifact["dependencies"].values():
        path = ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")

    p = HEIGHT / (2 * arb.pi())
    table = coefficient_table(CHECK_ORDER)
    for family in ("R_B", "R_L", "R_lower"):
        replay = changed_production_ball(table, family, p)
        saved = saved_complex(artifact["certificate"]["families"][family]["tail_integral_ball"])
        require(replay.overlaps(saved), f"changed-order {family} replay does not overlap")

    positive_labels = [arb("7.5") + index for index in range(4)]
    negative_labels = [arb("2.5") - index for index in range(4)]
    positive_enclosure, positive_direct = altered_sign_witness(positive_labels)
    negative_enclosure, negative_direct = altered_sign_witness(negative_labels)

    text = " ".join(NOTE.read_text(encoding="utf-8").split())
    for fragment in (
        "h_(n+1)=d/dy[h_n/D_q]",
        "sum_(m>=0)(delta+m)^(-k)",
        "They have not yet been joined with `P_240`",
        "No fitted circle constant is used",
    ):
        require(fragment in text, f"missing note fragment: {fragment}")

    print(
        "independently checked eight-round complementary tails; "
        f"positive-D overlap width {abs(positive_enclosure-positive_direct).str(8, more=True)}, "
        f"negative-D overlap width {abs(negative_enclosure-negative_direct).str(8, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
