#!/usr/bin/env python3
"""Certify the production lower finite cell by a reversed-roster endpoint expansion."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

from flint import acb, arb, ctx


ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_cell_reversed_roster_endpoint_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

QUARTER_DISK = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_quarter_disk_vertical_arc_tail_gate.json"
COMMON_CARRIER = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_common_carrier_subtraction_reduction_gate.json"
TRANSITION_PACKET = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_alternating_boundary_jet_packet_gate.json"
SADDLE_MAP = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_unequal_truncation_saddle_map_gate.json"

HEIGHT = 10_000_000_000
LOWER_BOUNDARY = Fraction(1243, 2)
Q_PLUS = Fraction(5_122_423, 2)
SOURCE_COUNT = 2_481_423
EXPANSION_ORDER = 18
SMALL_Y_SPLIT = 300
UPPER_SLABS = 4096
PRECISION_BITS = 384


CoefficientTable = list[dict[int, dict[tuple[int, int], Fraction]]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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


def arb_from_fraction(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def arb_record(value: arb, digits: int = 80) -> dict[str, Any]:
    return {
        "ball": value.str(digits, more=True),
        "lower_float": float(value.lower()),
        "mid_float": float(value.mid()),
        "upper_float": float(value.upper()),
    }


def complex_record(value: acb, digits: int = 80) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def complex_from_record(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def add_complex_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def coefficient_table(order: int) -> CoefficientTable:
    """Return exact monomial tables for h_(n+1)=(h_n/d)' and h_0=y^-1/2.

    A monomial key (j,e2) denotes p^j y^(e2/2).  The surrounding integer
    key is the power of d in c_(n,k)d^-k.
    """

    require(order >= 1, "positive expansion order required")
    table: CoefficientTable = [{0: {(0, -1): Fraction(1)}}]
    for _ in range(order):
        next_level: dict[int, dict[tuple[int, int], Fraction]] = {}
        for power, terms in table[-1].items():
            for (p_power, exponent2), coefficient in terms.items():
                updates = (
                    (
                        power + 1,
                        (p_power, exponent2 - 2),
                        coefficient * Fraction(exponent2, 2),
                    ),
                    (power + 2, (p_power, exponent2), -coefficient * (power + 1)),
                    (
                        power + 2,
                        (p_power + 1, exponent2 - 4),
                        coefficient * (power + 1),
                    ),
                )
                for new_power, monomial, value in updates:
                    destination = next_level.setdefault(new_power, {})
                    destination[monomial] = destination.get(monomial, Fraction(0)) + value
        table.append(
            {
                power: {monomial: value for monomial, value in terms.items() if value}
                for power, terms in next_level.items()
            }
        )
    return table


def half_power(base: arb, exponent2: int) -> arb:
    return base.sqrt() ** exponent2


def evaluate_coefficient(terms: dict[tuple[int, int], Fraction], p: arb, y: arb) -> arb:
    value = arb(0)
    for (p_power, exponent2), coefficient in terms.items():
        value += arb_from_fraction(coefficient) * p**p_power * half_power(y, exponent2)
    return value


def coefficient_absolute_upper(
    terms: dict[tuple[int, int], Fraction], p: arb, y_left: arb
) -> arb:
    value = arb(0)
    for (p_power, exponent2), coefficient in terms.items():
        require(exponent2 < 0, "coefficient monotonicity guard failed")
        value += (
            abs(arb_from_fraction(coefficient))
            * p**p_power
            * half_power(y_left, exponent2)
        )
    return value.upper()


def alternating_eta(z: arb, power: int) -> arb:
    require(power >= 1 and z.lower() > 0, "invalid shifted eta parameters")
    if power == 1:
        return (((z + 1) / 2).digamma() - (z / 2).digamma()) / 2
    exponent = arb(power)
    return (exponent.zeta(z / 2) - exponent.zeta((z + 1) / 2)) / (arb(2) ** power)


def finite_alternating_shifted_sum(a: arb, count: int, power: int) -> arb:
    """Return sum_(r=1)^count (-1)^r/(r+a)^power for odd count."""

    require(count > 0 and count % 2 == 1, "odd positive roster required")
    z = a + 1
    return -alternating_eta(z, power) - alternating_eta(z + count, power)


def small_y_remainder_bound(
    level: dict[int, dict[tuple[int, int], Fraction]],
    p: arb,
    count: int,
    split: arb,
) -> arb:
    """Bound 0<y<=split using d_r>=d_1>=p/(2y)."""

    value = arb(0)
    for power, terms in level.items():
        for (p_power, exponent2), coefficient in terms.items():
            denominator2 = exponent2 + 2 * power + 2
            require(denominator2 > 0, "small-y monomial is not integrable")
            integral_power = half_power(split, denominator2)
            value += (
                abs(arb_from_fraction(coefficient))
                * count
                * arb(2) ** power
                * p ** (p_power - power)
                * integral_power
                * 2
                / denominator2
            )
    return value


def upper_remainder_bound(
    level: dict[int, dict[tuple[int, int], Fraction]],
    p: arb,
    q_plus: arb,
    count: int,
    split: arb,
    boundary: arb,
    slabs: int,
) -> arb:
    """Bound split<=y<=boundary by monotone slab majorants."""

    require(slabs >= 32, "too few upper remainder slabs")
    width = (boundary - split) / slabs
    value = arb(0)
    for index in range(slabs):
        left = split + width * index
        right = split + width * (index + 1)
        d_lower = (p / right + right - q_plus + 1).lower()
        require(d_lower > 0, "nonstationary endpoint gap lost")
        for power, terms in level.items():
            coefficient_bound = coefficient_absolute_upper(terms, p, left)
            finite_bound = (arb(count) * d_lower ** (-power)).upper()
            infinite_bound = (
                d_lower ** (-power) + d_lower ** (1 - power) / (power - 1)
            ).upper()
            roster_bound = min(finite_bound, infinite_bound)
            value += width * coefficient_bound * roster_bound
    return value


def lower_cell_certificate(
    t: arb,
    boundary: arb,
    q_plus: arb,
    count: int,
    order: int,
    split: arb,
    upper_slabs: int,
    resource_mode: str,
) -> dict[str, Any]:
    pi = arb.pi()
    p = t / (2 * pi)
    q_minus = q_plus - count
    a = p / boundary + boundary - q_plus
    nearest_gap = a + 1
    require(count % 2 == 1, "production source count must be odd")
    require(boundary < p.sqrt(), "lower branch must lie before the central saddle")
    require(nearest_gap.lower() > 0, "a source phase is stationary in the lower cell")
    require(split > 0 and split < boundary, "invalid small-y split")
    split_guard = p / (2 * split) + split - q_plus + 1
    require(split_guard.lower() > 0, "small-y d>=p/(2y) guard failed")

    table = coefficient_table(order)
    partial = acb(0)
    terms: list[dict[str, Any]] = []
    common_phase = acb(
        0,
        -t * boundary.log() - pi * boundary**2 + 2 * pi * q_plus * boundary,
    ).exp()
    for n in range(order):
        endpoint_sum = arb(0)
        for power, coefficient_terms in table[n].items():
            coefficient = evaluate_coefficient(coefficient_terms, p, boundary)
            endpoint_sum += coefficient * finite_alternating_shifted_sum(
                a, count, power + 1
            )
        factor = acb(0, 1) / (2 * pi) * (-acb(0, 1) / (2 * pi)) ** n
        term = common_phase * factor * endpoint_sum
        partial += term
        terms.append({"order": n, "term_ball": complex_record(term, 60)})

    low_raw = small_y_remainder_bound(table[order], p, count, split)
    high_raw = upper_remainder_bound(
        table[order], p, q_plus, count, split, boundary, upper_slabs
    )
    scale = (2 * pi) ** (-order)
    low = (scale * low_raw).upper()
    high = (scale * high_raw).upper()
    remainder = (low + high).upper()
    enclosure = add_complex_error(partial, remainder)
    endpoint_slope = 2 * pi * nearest_gap

    return {
        "precision_bits": ctx.prec,
        "resource_mode": resource_mode,
        "t": t.str(50, more=True),
        "p_equals_t_over_2pi_ball": p.str(80, more=True),
        "boundary_L": boundary.str(50, more=True),
        "q_minus": q_minus.str(50, more=True),
        "q_plus": q_plus.str(50, more=True),
        "source_count": count,
        "expansion_order": order,
        "small_y_split": split.str(50, more=True),
        "upper_remainder_slabs": upper_slabs,
        "endpoint_offset_a_ball": arb_record(a),
        "nearest_phase_gap_a_plus_1_ball": arb_record(nearest_gap),
        "nearest_endpoint_slope_ball": arb_record(endpoint_slope),
        "common_endpoint_phase_ball": complex_record(common_phase),
        "endpoint_terms": terms,
        "endpoint_partial_sum_ball": complex_record(partial),
        "small_y_remainder_bound_ball": arb_record(low),
        "upper_interval_remainder_bound_ball": arb_record(high),
        "total_remainder_bound_ball": arb_record(remainder),
        "lower_cell_integral_ball": complex_record(enclosure),
        "identity": "V_L=integral_0^L y^(-s)exp(-i*pi*y^2)D_W(y)dy",
        "passed": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Reversed-roster endpoint certificate for the lower finite cell

Date: 2026-08-27

Status: rigorous production lower-cell complex enclosure certified; ordinary
Gamma-subtracted join remains open.

Put `s=1/2+it`, `p=t/(2*pi)`, `L=621.5`, and reverse the odd finite roster as

```text
D_W(y)=sum_(r=1)^M exp[2*pi*i(q_+-r)y],
q_+=2561211.5,       M=2481423.                     (LC1)
```

For the `r`th label absorb `y^(-it)` into the phase.  Its derivative is

```text
phi_r'(y)=-2*pi*d_r(y),
d_r(y)=p/y+y-q_++r.                                 (LC2)
```

On `0<y<=L`, `d_r` decreases and the nearest label satisfies

```text
d_1(L)={c['nearest_phase_gap_a_plus_1_ball']['ball']},
|phi_1'(L)|={c['nearest_endpoint_slope_ball']['ball']}.
```

Thus no lower-cell label has a stationary point.  At the half-integer endpoint,

```text
exp[2*pi*i(q_+-r)L]=exp(2*pi*i*q_+*L)(-1)^r.        (LC3)
```

Starting with `h_0(y)=y^(-1/2)`, define

```text
h_(n+1,r)=d/dy[h_(n,r)/d_r].                        (LC4)
```

Repeated integration by parts is exact for the finite roster.  Every zero-end
boundary term vanishes, while every `L`-end term is a finite linear combination
of

```text
sum_(r=1)^M (-1)^r/(r+a)^k,
a=p/L+L-q_+={c['endpoint_offset_a_ball']['ball']}.   (LC5)
```

These finite alternating sums are evaluated without iterating the roster by
the shifted eta/Hurwitz identity

```text
sum_(r=1)^M (-1)^r/(r+a)^k
 =-eta_k(a+1)-eta_k(a+1+M),                         (LC6)

eta_k(z)=2^(-k)[zeta(k,z/2)-zeta(k,(z+1)/2)],       k>1,
eta_1(z)=[psi((z+1)/2)-psi(z/2)]/2.
```

The order-{c['expansion_order']} endpoint sum is

```text
{c['endpoint_partial_sum_ball']['real_ball']}
+ i {c['endpoint_partial_sum_ball']['imag_ball']}.
```

For the discarded integral, split at `y={SMALL_Y_SPLIT}`.  Below the split,
`d_r>=d_1>=p/(2y)` gives a termwise exact-monomial integral.  Above it,
{c['upper_remainder_slabs']} monotone slabs use

```text
sum_(r=1)^M d_r^(-k)
 <=min[M d_1^(-k), d_1^(-k)+d_1^(1-k)/(k-1)].       (LC7)
```

The certified remainder is

```text
small-y: {c['small_y_remainder_bound_ball']['ball']}
upper:   {c['upper_interval_remainder_bound_ball']['ball']}
total:   {c['total_remainder_bound_ball']['ball']}.
```

Therefore the complete lower finite cell is

```text
V_L = {c['lower_cell_integral_ball']['real_ball']}
    + i {c['lower_cell_integral_ball']['imag_ball']}.           (LC8)
```

This is a complex ball retained for later addition.  No absolute value of
`V_L` is substituted into the joined packet.

Pi provenance: every `pi` in (LC1)--(LC8) comes from the inherited quadratic
Riemann-Siegel/Fresnel phase, the exact Fourier spacing, or `p=t/(2*pi)`.
No geometric fit or inserted circle constant is used.

Proof boundary: exact reversed-roster identity, nonstationary lower-cell phase,
finite alternating-Hurwitz endpoint expansion, and rigorous production complex
lower-cell enclosure only.  This gate does not yet enclose the ordinary
Gamma-subtracted packet, the complete joined packet, `J_Z`, or `D_K`, and it
does not prove a non-A, all-height, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion.
"""


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    dependencies = {
        "quarter_disk": QUARTER_DISK,
        "common_carrier": COMMON_CARRIER,
        "transition_packet": TRANSITION_PACKET,
        "saddle_map": SADDLE_MAP,
    }
    for path in (*dependencies.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")

    quarter = load_json(QUARTER_DISK)
    require(quarter.get("passed") is True, "quarter-disk gate is not certified")
    require(quarter["production"]["t"] == str(HEIGHT), "height drift")
    require(quarter["production"]["label_count"] == SOURCE_COUNT, "source count drift")
    require(quarter["production"]["q_plus"] == "2561211.5", "q_plus drift")

    certificate = lower_cell_certificate(
        arb(HEIGHT),
        arb_from_fraction(LOWER_BOUNDARY),
        arb_from_fraction(Q_PLUS),
        SOURCE_COUNT,
        EXPANSION_ORDER,
        arb(SMALL_Y_SPLIT),
        UPPER_SLABS,
        resource_mode,
    )
    require(
        arb(certificate["total_remainder_bound_ball"]["ball"]).upper() < arb("3e-16"),
        "production lower-cell remainder misses target",
    )

    artifact: dict[str, Any] = {
        "kind": "rh_c1_hardy_equation4_joined_packet_lower_cell_reversed_roster_endpoint_gate",
        "date": "2026-08-27",
        "status": "production_lower_cell_reversed_roster_endpoint_complex_ball_certified",
        "passed": True,
        "production": {
            "t": str(HEIGHT),
            "L": "621.5",
            "q_minus": "79788.5",
            "q_plus": "2561211.5",
            "source_count": SOURCE_COUNT,
        },
        "exact_identities": {
            "reversed_roster": "D_W(y)=sum_(r=1)^M exp(2*pi*i*(q_+-r)*y)",
            "phase_gap": "d_r(y)=p/y+y-q_++r and phi_r'(y)=-2*pi*d_r(y)",
            "endpoint_parity": "exp(2*pi*i*(q_+-r)*L)=exp(2*pi*i*q_+*L)*(-1)^r",
            "recurrence": "h_0=y^(-1/2), h_(n+1,r)=d/dy[h_(n,r)/d_r]",
            "finite_alternating_sum": "sum_(r=1)^M(-1)^r/(r+a)^k=-eta_k(a+1)-eta_k(a+1+M) for odd M",
        },
        "certificate": certificate,
        "decision": {
            "all_lower_source_labels_nonstationary": True,
            "full_finite_roster_collapsed_without_iteration": True,
            "lower_cell_complex_ball_certified": True,
            "lower_cell_norm_taken_before_join": False,
            "ordinary_Gamma_subtracted_packet_bounded": False,
            "complete_joined_packet_bounded": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in dependencies.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Join this retained lower-cell complex ball to the ordinary 622..39852 Gamma-subtracted "
            "current and the certified transition packet before taking a norm. Use the completed-series "
            "endpoint block 622..860 and the contracting 861..39936 core; retain the certified upper arc."
        ),
        "proof_boundary": (
            "Rigorous production lower finite-cell complex enclosure only. No ordinary Gamma-subtracted "
            "packet, complete joined packet, J_Z or D_K enclosure, non-A, all-height, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified reversed-roster lower finite cell", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
