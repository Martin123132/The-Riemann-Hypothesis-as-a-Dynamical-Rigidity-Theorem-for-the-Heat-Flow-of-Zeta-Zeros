#!/usr/bin/env python3
"""Build the Phi Taylor sign scout for the heat-flow cone-entry route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md"


# Coefficient of u^j in
# (2*q^2*exp(9*u)-3*q*exp(5*u))*exp(-q*(exp(4*u)-1)).
# The full Phi coefficient is sum_n exp(-pi*n^2) * P_j(pi*n^2).
POLYS: dict[int, dict[int, Fraction]] = {
    0: {1: Fraction(-3), 2: Fraction(2)},
    1: {1: Fraction(-15), 2: Fraction(30), 3: Fraction(-8)},
    2: {1: Fraction(-75, 2), 2: Fraction(165), 3: Fraction(-112), 4: Fraction(16)},
    3: {
        1: Fraction(-375, 6),
        2: Fraction(3270, 6),
        3: Fraction(-4232, 6),
        4: Fraction(1440, 6),
        5: Fraction(-128, 6),
    },
    4: {
        1: Fraction(-1875, 24),
        2: Fraction(30930, 24),
        3: Fraction(-68096, 24),
        4: Fraction(41408, 24),
        5: Fraction(-8448, 24),
        6: Fraction(512, 24),
    },
    5: {
        1: Fraction(-9375, 120),
        2: Fraction(285870, 120),
        3: Fraction(-1008968, 120),
        4: Fraction(976320, 120),
        5: Fraction(-343040, 120),
        6: Fraction(46592, 120),
        7: Fraction(-2048, 120),
    },
    6: {
        1: Fraction(-46875, 720),
        2: Fraction(2610330, 720),
        3: Fraction(-14260064, 720),
        4: Fraction(20633312, 720),
        5: Fraction(-11109120, 720),
        6: Fraction(2536960, 720),
        7: Fraction(-245760, 720),
        8: Fraction(8192, 720),
    },
}

CERTIFIED_DEGREES = (0, 2, 4, 6)


@dataclass(frozen=True)
class CoefficientBall:
    id: str
    coefficient: str
    taylor_degree: int
    finite_sum: str
    tail_radius: str
    enclosure: str
    sign: str


@dataclass(frozen=True)
class SignCombination:
    id: str
    role: str
    formula: str
    enclosure: str
    certified_sign: str
    consequence: str
    proof_boundary: str


def arb_fraction(value: Fraction) -> flint.arb:
    return flint.arb(value.numerator) / flint.arb(value.denominator)


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_negative(value: flint.arb) -> bool:
    return bool(value < 0 and not value.contains(0))


def sign_name(value: flint.arb) -> str:
    if arb_positive(value):
        return "positive"
    if arb_negative(value):
        return "negative"
    return "unresolved"


def eval_poly(j: int, q: flint.arb) -> flint.arb:
    total = flint.arb(0)
    for degree, coeff in POLYS[j].items():
        total += arb_fraction(coeff) * (q**degree)
    return total


def tail_power_bound(cutoff_n: int, power_degree: int, pi: flint.arb) -> flint.arb:
    """Bound sum_{n>cutoff_n} n^(2m) exp(-pi*n^2) by a geometric tail."""
    n0 = cutoff_n + 1
    n = flint.arb(n0)
    first = (n ** (2 * power_degree)) * (-(pi * flint.arb(n0 * n0))).exp()
    ratio = ((flint.arb(n0 + 1) / n) ** (2 * power_degree)) * (
        -(pi * flint.arb(2 * n0 + 1))
    ).exp()
    if not arb_positive(flint.arb(1) - ratio):
        raise RuntimeError(f"tail ratio did not certify below one for m={power_degree}")
    return first / (flint.arb(1) - ratio)


def tail_radius_bound(j: int, cutoff_n: int, pi: flint.arb) -> flint.arb:
    total = flint.arb(0)
    for degree, coeff in POLYS[j].items():
        total += arb_fraction(abs(coeff)) * (pi**degree) * tail_power_bound(cutoff_n, degree, pi)
    return total


def upper_radius_string(value: flint.arb, safety_factor: int = 4) -> str:
    """Return a decimal radius safely above a positive Arb ball."""
    if not arb_positive(value):
        raise RuntimeError(f"non-positive tail radius {value}")
    mid, rad, exp10 = value.mid_rad_10exp()
    mantissa = (abs(int(mid)) + int(rad)) * safety_factor
    return f"{mantissa}e{int(exp10)}"


def error_ball(radius: flint.arb) -> flint.arb:
    return flint.arb(f"[0 +/- {upper_radius_string(radius)}]")


def coefficient_ball(j: int, cutoff_n: int, pi: flint.arb) -> tuple[flint.arb, CoefficientBall]:
    finite = flint.arb(0)
    for n in range(1, cutoff_n + 1):
        q = pi * flint.arb(n * n)
        finite += (-q).exp() * eval_poly(j, q)
    tail = tail_radius_bound(j, cutoff_n, pi)
    enclosure = finite + error_ball(tail)
    return enclosure, CoefficientBall(
        id=f"phi_taylor_c{j}",
        coefficient=f"c{j}",
        taylor_degree=j,
        finite_sum=finite.str(40),
        tail_radius=upper_radius_string(tail),
        enclosure=enclosure.str(40),
        sign=sign_name(enclosure),
    )


def build_scout(cutoff_n: int, precision_bits: int) -> dict:
    flint.ctx.prec = precision_bits
    pi = flint.arb.pi()

    coefficient_values: dict[int, flint.arb] = {}
    coefficient_rows: list[CoefficientBall] = []
    for degree in CERTIFIED_DEGREES:
        value, row = coefficient_ball(degree, cutoff_n, pi)
        coefficient_values[degree] = value
        coefficient_rows.append(row)

    c0 = coefficient_values[0]
    c2 = coefficient_values[2]
    c4 = coefficient_values[4]
    c6 = coefficient_values[6]
    if not arb_positive(c0):
        raise RuntimeError("c0 did not certify positive")

    a = c2 / c0
    b = c4 / c0
    c = c6 / c0
    upper_wall = 2 * b - a * a
    monotone_wall = 2 * (a**3 - 3 * a * b + 3 * c)

    if not arb_negative(upper_wall):
        raise RuntimeError("upper-wall Taylor combination did not certify negative")
    if not arb_positive(monotone_wall):
        raise RuntimeError("monotone-wall Taylor combination did not certify positive")

    combinations = [
        SignCombination(
            id="ptces_01_upper_wall_sign",
            role="fixed_k_upper_wall_subleading_sign",
            formula="2*b-a^2=(2*c4*c0-c2^2)/c0^2",
            enclosure=upper_wall.str(40),
            certified_sign="negative",
            consequence="For fixed k, the first nonzero correction to log x_k at order T^-2 points strictly below the upper wall x_k=1.",
            proof_boundary="Local Taylor sign only; it does not give a uniform-in-k cone-entry theorem or a controlled asymptotic remainder.",
        ),
        SignCombination(
            id="ptces_02_monotone_wall_sign",
            role="fixed_k_monotone_wall_subleading_sign",
            formula="2*(a^3-3*a*b+3*c)",
            enclosure=monotone_wall.str(40),
            certified_sign="positive",
            consequence="For fixed k, the first nonzero correction to log(x_{k+1}/x_k) at order T^-3 points toward increasing contractions.",
            proof_boundary="Local Taylor sign only; it does not give a uniform-in-k cone-entry theorem or a controlled asymptotic remainder.",
        ),
    ]

    return {
        "kind": "jensen_window_pf_phi_taylor_cone_entry_sign_scout",
        "date": "2026-07-06",
        "status": "finite_taylor_sign_certificate",
        "source_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_phi_taylor_cone_entry_sign_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_phi_taylor_cone_entry_sign_scout.py",
        "proof_boundary": (
            "Certified local Taylor-combination signs only. This confirms the "
            "fixed-k sign inputs isolated by the cone-entry asymptotic target, "
            "but it does not prove zeta cone entry, does not control the "
            "asymptotic remainder uniformly in k, does not prove the full "
            "infinite or collared finite ratio cone, does not prove jwpf_06, "
            "and leaves RH and Lambda <= 0 unsettled."
        ),
        "parameters": {
            "precision_bits": precision_bits,
            "tail_cutoff_n": cutoff_n,
            "certified_taylor_degrees": list(CERTIFIED_DEGREES),
            "tail_safety_factor": 4,
        },
        "taylor_model": {
            "phi_term": "(2*q^2*exp(9*u)-3*q*exp(5*u))*exp(-q*(exp(4*u)-1)), q=pi*n^2",
            "coefficient_formula": "c_j=sum_{n>=1} exp(-pi*n^2)*P_j(pi*n^2)",
            "normalizations": {
                "a": "c2/c0",
                "b": "c4/c0",
                "c": "c6/c0",
            },
            "tail_bound": (
                "For each polynomial degree m, bound sum_{n>N} n^(2*m)*exp(-pi*n^2) "
                "by first/(1-r), where first=(N+1)^(2*m)*exp(-pi*(N+1)^2) and "
                "r=((N+2)/(N+1))^(2*m)*exp(-pi*(2*N+3))."
            ),
        },
        "polynomial_rows": {
            f"P{j}": {str(degree): f"{coeff.numerator}/{coeff.denominator}" for degree, coeff in POLYS[j].items()}
            for j in CERTIFIED_DEGREES
        },
        "coefficient_enclosures": [asdict(row) for row in coefficient_rows],
        "ratio_enclosures": {
            "a=c2/c0": a.str(40),
            "b=c4/c0": b.str(40),
            "c=c6/c0": c.str(40),
        },
        "sign_combinations": [asdict(row) for row in combinations],
        "summary": {
            "coefficient_balls": len(coefficient_rows),
            "certified_signs": len(combinations),
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The local Phi Taylor combinations required by the fixed-k "
                "large-negative-lambda route have the desired certified signs. "
                "The remaining bridge is not these local signs, but controlled "
                "asymptotic remainders and uniform-in-k or collared finite cone coverage."
            ),
        },
        "invariants": [
            "The certificate uses the Phi series at u=0 only.",
            "The tail bound is attached as an explicit Arb error ball.",
            "The signs support only the fixed-k asymptotic route.",
            "The cone-entry theorem target remains open.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(scout: dict, path: Path) -> None:
    coeffs = scout["coefficient_enclosures"]
    combos = scout["sign_combinations"]
    summary = scout["summary"]

    coefficient_lines = "\n".join(
        f"{row['coefficient']}: {row['enclosure']} ({row['sign']})" for row in coeffs
    )
    combo_lines = "\n\n".join(
        f"{row['id']}:\n  {row['formula']}\n  enclosure: {row['enclosure']}\n  certified sign: {row['certified_sign']}\n  boundary: {row['proof_boundary']}"
        for row in combos
    )
    text = f"""# Jensen-Window PF Phi Taylor Cone-Entry Sign Scout

Date: 2026-07-06

Status: finite Taylor sign certificate. This is not a proof of zeta cone entry,
monotone contractions, Jensen-window PF-infinity, RH, or the Newman-direction
goal.

Artifact kind: `jensen_window_pf_phi_taylor_cone_entry_sign_scout`.

Proof boundary: this artifact certifies only the local Taylor-combination
signs used by the fixed-k large-negative-lambda asymptotic route. It does not
control the asymptotic remainder uniformly in `k`, does not prove the full
infinite or collared finite ratio cone, does not prove `jwpf_06`, and does not
settle `Lambda <= 0`.

Machine-readable certificate:

```text
work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_phi_taylor_cone_entry_sign_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_phi_taylor_cone_entry_sign_scout.py
```

Current result:

```text
validated Jensen-window PF Phi Taylor cone-entry sign scout: {summary['coefficient_balls']} coefficient balls, {summary['certified_signs']} certified signs, {summary['ready_to_apply_rows']} ready-to-apply rows, 0 issues
```

## Role

The cone-entry asymptotic target isolated the fixed-k sign requirements:

```text
Phi(u)=c0+c2*u^2+c4*u^4+c6*u^6+...
a=c2/c0
b=c4/c0
c=c6/c0
2*b-a^2 < 0
2*(a^3-3*a*b+3*c) > 0
```

This certificate supplies those local signs with an explicit finite sum plus
tail bound. It does not supply the missing uniform-in-`k` or collared finite
cone coverage.

Source target:

```text
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```

## Coefficient Enclosures

```text
{coefficient_lines}
```

Parameters:

```text
precision_bits = {scout['parameters']['precision_bits']}
tail_cutoff_n = {scout['parameters']['tail_cutoff_n']}
tail_safety_factor = {scout['parameters']['tail_safety_factor']}
```

Tail method:

```text
{scout['taylor_model']['tail_bound']}
```

## Certified Sign Combinations

```text
{combo_lines}
```

## Consequence

The local Phi Taylor signs are no longer the immediate obstruction in the
fixed-k cone-entry route. The remaining hard work is to convert the fixed-k
asymptotic indication into a theorem with controlled remainders and either
uniform-in-`k` coverage or a rigorous finite-prefix collar plus tail argument.
"""
    path.write_text(text, encoding="utf-8")


def write_json(scout: dict, path: Path) -> None:
    path.write_text(json.dumps(scout, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--tail-cutoff-n", type=int, default=12)
    parser.add_argument("--precision-bits", type=int, default=256)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    scout = build_scout(args.tail_cutoff_n, args.precision_bits)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    write_json(scout, args.out_json)
    write_note(scout, args.note)
    print(
        "wrote Jensen-window PF Phi Taylor cone-entry sign scout: "
        f"{args.out_json.relative_to(REPO_ROOT).as_posix()} and {args.note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
