#!/usr/bin/env python3
"""Build the exact coefficient-to-physical-RAE roster inversion gate."""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
FLINT_ROOT = Path(
    os.environ.get(
        "RH_PYTHON_FLINT_ROOT",
        r"C:\Users\ollet\Documents\Codex\third_party\python_flint_0_8_0",
    )
)
sys.path.insert(0, str(FLINT_ROOT))

try:
    from flint import arb, ctx
except ImportError as exc:  # pragma: no cover
    raise RuntimeError(f"python-flint runtime unavailable at {FLINT_ROOT}") from exc


SOURCE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
)
FIXTURE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10/enabled"
)
TELEMETRY = FIXTURE / "chain_telemetry.jsonl"
CHECKPOINT = FIXTURE / "checkpoint.jsonl"
INPUT = FIXTURE / "inputs3.nml"
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.md"
)
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.py"
)
SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"
PRECISION_LADDER = (192, 320)
COEFFICIENT_REPLAY_TOLERANCE = Fraction(1, 2**88)
RAE_RECOVERY_TOLERANCE = Fraction(1, 2**55)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def fraction_record(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def decimal_text(value: Fraction, digits: int = 70) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".55E")


def binary128_fraction(payload: str) -> Fraction:
    require(isinstance(payload, str) and len(payload) == 32, "invalid binary128 payload")
    bits = int(payload, 16)
    sign = -1 if bits >> 127 else 1
    exponent = (bits >> 112) & 0x7FFF
    fraction = bits & ((1 << 112) - 1)
    require(exponent != 0x7FFF, "non-finite binary128 payload")
    if exponent == 0:
        if fraction == 0:
            return Fraction(0)
        mantissa = fraction
        binary_exponent = 1 - 16383 - 112
    else:
        mantissa = (1 << 112) + fraction
        binary_exponent = exponent - 16383 - 112
    if binary_exponent >= 0:
        return Fraction(sign * mantissa * (2**binary_exponent), 1)
    return Fraction(sign * mantissa, 2 ** (-binary_exponent))


def fraction_ball(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def dyadic(value: arb) -> list[int]:
    mantissa, exponent = value.man_exp()
    return [int(mantissa), int(exponent)]


def arb_record(value: arb, digits: int = 54) -> dict[str, Any]:
    lower = value.lower()
    upper = value.upper()
    return {
        "display": value.str(digits),
        "lower_dyadic": dyadic(lower),
        "upper_dyadic": dyadic(upper),
        "lower_decimal": lower.str(digits, radius=False),
        "upper_decimal": upper.str(digits, radius=False),
    }


def nearest_integer_from_ball(value: arb) -> int:
    midpoint = float(value.mid())
    if midpoint >= 0:
        candidate = math.floor(midpoint + 0.5)
    else:
        candidate = math.ceil(midpoint - 0.5)
    margin = fraction_ball(Fraction(1, 2)) - abs(value - candidate)
    require(margin.lower() > 0, "nearest-integer selector is not separated")
    return candidate


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "pi_from_arctangent": "p=4*ATAN(C)",
        "pari_two_pi": "v           = Pi2n(1_8, prec) ! 2pi",
        "pari_pi": "y           = gdiv(v, stor(2_8,prec)) ! pi",
        "physical_scale": "a_pari      = gsqrt(v ,prec)                                ! sqrt(8*t/pi)",
        "block_length": "MT=floor(((((et/p)**rmc5)*a)**rmc2)*((e-1.0)**rmc3)/(2.0**rmc4),dp1)",
        "even_length_adjustment": "if (mod(MT,2).lt.0.25) then",
        "first_pivot": "rae=RN1+2.0+real(MT,dp)",
        "roster_increment": "RN4=2.0*(real(MT,dp)+1)*(real(jsum-1,dp))",
        "roster_pivot": "rae=raestartofblock  + RN4",
        "central_ratio": "e=(rae*ain)**2",
        "hyperbolic_parameter": "pc=2*(e+sqrt(e)*s)-1",
        "xx_definition": "xx=1/(pc-1)",
        "cubic_scale": "con2=((sqrt(pc)*xx)**3)/(aar(numbercalc))",
        "linear_plus": "initialcoeff(1)=xx/2.0+s1/4.0+con2",
        "linear_minus": "initialcoeffcc(1)=initialcoeff(1)-s1/2.0-2*con2",
        "quadratic_plus": "initialcoeff(2)=xx/2.0+2*con2",
        "quadratic_minus": "initialcoeffcc(2)=initialcoeff(2)-4*con2",
        "cubic_plus": "initialcoeff(i)=4*con2/3.0",
        "cubic_minus": "initialcoeffcc(i)=-initialcoeff(i)",
    }
    locations: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [index for index, line in enumerate(lines, 1) if token in line]
        require(matches, f"source RAE-map token missing: {token}")
        locations[name] = matches[0]
    return locations


def load_chain_headers() -> dict[tuple[int, int], dict[str, Any]]:
    headers: dict[tuple[int, int], dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        if record.get("type") != "chain":
            continue
        key = (int(record["sum_index"]), int(record["branch"]))
        require(key not in headers, f"duplicate telemetry header {key}")
        headers[key] = record
    require(len(headers) == 64, "telemetry header count drift")
    require(
        sorted(headers) == [(sum_index, branch) for sum_index in range(1, 33) for branch in (1, 2)],
        "telemetry physical roster gap",
    )
    return headers


def load_checkpoint_rows() -> dict[tuple[int, int], dict[str, Any]]:
    rows: dict[tuple[int, int], dict[str, Any]] = {}
    for line in CHECKPOINT.read_text(encoding="utf-8").splitlines():
        record = json.loads(line, parse_float=Decimal)
        key = (int(record["stage"]), int(record["unit"]))
        rows[key] = record
    return rows


def exact_symbolic_fixtures() -> list[dict[str, Any]]:
    fixtures: list[dict[str, Any]] = []
    for index, w in enumerate((Fraction(3, 2), Fraction(5, 3), Fraction(7, 4), Fraction(9, 5)), 1):
        pc = w**2
        r = (w + 1 / w) / 2
        s = (w - 1 / w) / 2
        e = r**2
        xx = 1 / (pc - 1)
        require(2 * (e + r * s) - 1 == pc, "exact pc inversion fixture failed")
        require(r**2 - s**2 == 1, "exact hyperbola fixture failed")
        require((w + 1 / w) / 2 == r, "exact r inversion fixture failed")
        for branch, con2 in ((1, Fraction(index, 10_000)), (2, Fraction(-index, 10_000))):
            a2 = xx / 2 + 2 * con2
            a3 = 4 * con2 / 3
            require(2 * a2 - 3 * a3 == xx, "exact coefficient cancellation failed")
            fixtures.append(
                {
                    "index": len(fixtures) + 1,
                    "branch": branch,
                    "sqrt_pc": fraction_record(w),
                    "r": fraction_record(r),
                    "xx": fraction_record(xx),
                    "con2": fraction_record(con2),
                }
            )
    return fixtures


def ideal_coefficients(t: Fraction, rae: int) -> dict[str, Any]:
    pi = arb.pi()
    a = (8 * fraction_ball(t) / pi).sqrt()
    r = arb(rae) / a
    e = r**2
    s = (e - 1).sqrt()
    pc = 2 * (e + e.sqrt() * s) - 1
    xx = 1 / (pc - 1)
    con2 = (pc.sqrt() * xx) ** 3 / a
    s1 = a / pc.sqrt()
    raw_plus = xx / 2 + s1 / 4 + con2
    raw_minus = raw_plus - s1 / 2 - 2 * con2
    lift_plus = nearest_integer_from_ball(raw_plus)
    lift_minus = nearest_integer_from_ball(raw_minus)
    return {
        "pi": pi,
        "a": a,
        "r": r,
        "pc": pc,
        "xx": xx,
        "con2": con2,
        "branches": {
            1: {
                "lift": lift_plus,
                "raw_a1": raw_plus,
                "coefficients": [raw_plus - lift_plus, xx / 2 + 2 * con2, 4 * con2 / 3],
            },
            2: {
                "lift": lift_minus,
                "raw_a1": raw_minus,
                "coefficients": [raw_minus - lift_minus, xx / 2 - 2 * con2, -4 * con2 / 3],
            },
        },
    }


def evaluate_precision(
    precision: int,
    t: Fraction,
    rn1_start: int,
    mt: int,
    headers: dict[tuple[int, int], dict[str, Any]],
) -> dict[str, Any]:
    ctx.prec = precision
    coefficient_tolerance = fraction_ball(COEFFICIENT_REPLAY_TOLERANCE)
    rae_tolerance = fraction_ball(RAE_RECOVERY_TOLERANCE)
    rows: list[dict[str, Any]] = []
    all_replay_gaps: list[tuple[float, arb, int, int, int]] = []
    all_recovery_errors: list[tuple[float, arb, int, int]] = []
    modulo_margins: list[tuple[float, arb, int, int]] = []
    pair_x_gaps: list[tuple[Fraction, int]] = []
    recovered_roster: list[int] = []

    for sum_index in range(1, 33):
        target_rae = rn1_start + mt + 2 + 2 * (mt + 1) * (sum_index - 1)
        ideal = ideal_coefficients(t, target_rae)
        branch_rows: list[dict[str, Any]] = []
        x_values: dict[int, Fraction] = {}
        recovered_values: dict[int, arb] = {}

        for branch in (1, 2):
            header = headers[(sum_index, branch)]
            require(int(header["block"]) == 20, "telemetry block drift")
            logged = [binary128_fraction(item) for item in header["initial_coefficients_hex"]]
            require(len(logged) == 3, "initial coefficient degree drift")
            x_value = 2 * logged[1] - 3 * logged[2]
            require(x_value > 0, "recovered xx is nonpositive")
            x_values[branch] = x_value

            x_ball = fraction_ball(x_value)
            pc_hat = 1 + 1 / x_ball
            root_pc = pc_hat.sqrt()
            r_hat = (root_pc + 1 / root_pc) / 2
            recovered_rae = r_hat * ideal["a"]
            recovery_error = abs(recovered_rae - target_rae)
            require(
                recovery_error.upper() < rae_tolerance.lower(),
                f"chain {header['chain']} RAE recovery tolerance failed",
            )
            require(
                recovery_error.upper() < fraction_ball(Fraction(1, 2)).lower(),
                f"chain {header['chain']} nearest physical pivot not unique",
            )
            recovered_values[branch] = recovered_rae
            all_recovery_errors.append(
                (float(recovery_error.upper()), recovery_error, int(header["chain"]), branch)
            )

            ideal_branch = ideal["branches"][branch]
            replay_gaps: list[arb] = []
            for coefficient_index, (saved, expected) in enumerate(
                zip(logged, ideal_branch["coefficients"]), 1
            ):
                gap = abs(fraction_ball(saved) - expected)
                require(
                    gap.upper() < coefficient_tolerance.lower(),
                    f"chain {header['chain']} coefficient {coefficient_index} replay gap",
                )
                replay_gaps.append(gap)
                all_replay_gaps.append(
                    (
                        float(gap.upper()),
                        gap,
                        int(header["chain"]),
                        branch,
                        coefficient_index,
                    )
                )

            modulo_margin = fraction_ball(Fraction(1, 2)) - abs(
                ideal_branch["raw_a1"] - ideal_branch["lift"]
            )
            require(modulo_margin.lower() > 0, "linear modulo selector touches a half-integer")
            modulo_margins.append(
                (float(modulo_margin.lower()), modulo_margin, int(header["chain"]), branch)
            )
            branch_rows.append(
                {
                    "chain": int(header["chain"]),
                    "branch": branch,
                    "initial_coefficients_hex": header["initial_coefficients_hex"],
                    "linear_modulo_lift": int(ideal_branch["lift"]),
                    "linear_modulo_margin": arb_record(modulo_margin),
                    "xx_exact": fraction_record(x_value),
                    "xx_decimal": decimal_text(x_value),
                    "pc_from_saved_coefficients": arb_record(pc_hat),
                    "r_from_saved_coefficients": arb_record(r_hat),
                    "recovered_rae": arb_record(recovered_rae),
                    "recovered_rae_error": arb_record(recovery_error),
                    "nearest_physical_pivot": target_rae,
                    "ideal_formula_replay_gaps": {
                        f"a{index}": arb_record(gap) for index, gap in enumerate(replay_gaps, 1)
                    },
                }
            )

        pair_x_gap = abs(x_values[1] - x_values[2])
        pair_x_gaps.append((pair_x_gap, sum_index))
        branch_recovery_gap = abs(recovered_values[1] - recovered_values[2])
        require(
            branch_recovery_gap.upper() < 2 * rae_tolerance.lower(),
            f"sum index {sum_index} branch recovery disagreement",
        )
        recovered_roster.append(target_rae)
        rows.append(
            {
                "block": 20,
                "sum_index": sum_index,
                "target_rae": target_rae,
                "roster_formula": f"{rn1_start}+{mt}+2+2*({mt}+1)*({sum_index}-1)",
                "ideal_r_equals_rae_over_a": arb_record(ideal["r"]),
                "ideal_xx": arb_record(ideal["xx"]),
                "pair_xx_exact_gap": fraction_record(pair_x_gap),
                "pair_xx_gap_decimal": decimal_text(pair_x_gap),
                "branch_recovered_rae_gap": arb_record(branch_recovery_gap),
                "branches": branch_rows,
            }
        )

    max_replay = max(all_replay_gaps, key=lambda item: item[0])
    max_recovery = max(all_recovery_errors, key=lambda item: item[0])
    min_modulo = min(modulo_margins, key=lambda item: item[0])
    max_pair_x = max(pair_x_gaps, key=lambda item: item[0])
    scale = ideal_coefficients(t, recovered_roster[0])
    return {
        "precision_bits": precision,
        "pi": arb_record(scale["pi"]),
        "a": arb_record(scale["a"]),
        "rows": rows,
        "selectors": [
            [
                row["branches"][0]["linear_modulo_lift"],
                row["branches"][1]["linear_modulo_lift"],
            ]
            for row in rows
        ],
        "recovered_roster": recovered_roster,
        "maximum_coefficient_replay_gap": {
            "chain": max_replay[2],
            "branch": max_replay[3],
            "coefficient": max_replay[4],
            "ball": arb_record(max_replay[1]),
        },
        "maximum_recovered_rae_error": {
            "chain": max_recovery[2],
            "branch": max_recovery[3],
            "ball": arb_record(max_recovery[1]),
        },
        "minimum_linear_modulo_margin": {
            "chain": min_modulo[2],
            "branch": min_modulo[3],
            "ball": arb_record(min_modulo[1]),
        },
        "maximum_pair_xx_exact_gap": {
            "sum_index": max_pair_x[1],
            "exact": fraction_record(max_pair_x[0]),
            "decimal": decimal_text(max_pair_x[0]),
        },
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    fixture = artifact["fixture"]
    return f"""# Hardy coefficient-to-physical-RAE roster gate

Date: 2026-08-06

Status: exact inversion and rigorous 32-pivot fixture recovery; not a proof and continuous q-cell coverage remains open

## What was proved

For either source branch, write the saved cubic phase as

```text
F(x) = a1*x + a2*x^2 + a3*x^3.
```

Before binary128 rounding the accepted source has

```text
branch +: a2 = xx/2 + 2*c,  a3 =  4*c/3,
branch -: a2 = xx/2 - 2*c,  a3 = -4*c/3.
```

Consequently both branches obey the same exact cancellation

```text
xx = 2*a2 - 3*a3.
```

The source also defines `xx=1/(pc-1)`.  If `r=rae/a>1`, then
`pc=2(r^2+r*sqrt(r^2-1))-1=(r+sqrt(r^2-1))^2`.  Since
`(r+s)(r-s)=1`, this gives the exact inverse

```text
pc = 1 + 1/xx,
r  = (sqrt(pc) + 1/sqrt(pc))/2.
```

Eight rational fixtures check these identities without floating-point arithmetic.

## Where pi comes from

The source does not insert an unexplained decimal.  Its ordinary real constant is
`p=4*atan(1)`.  For the more accurate central scale, PARI supplies `2*pi`, the
code divides that value by `2`, and then computes

```text
a = sqrt(8*t/pi).
```

Thus this is the usual circle constant characterized by `pi=4*atan(1)`, with an
independent high-precision PARI evaluation used by the physical scale.  It is not
an arbitrary circle or polygon chosen to fit the data.

## Independent physical roster

The checkpoint immediately before block `{fixture['block']}` records
`RN1={fixture['rn1_start']}`.  The completed block records `MT={fixture['mt']}`
and `{fixture['aenums']}` sums.  The accepted loop therefore gives

```text
rae_j = {fixture['rn1_start']} + {fixture['mt']} + 2
        + 2*({fixture['mt']}+1)*(j-1)
      = {fixture['first_rae']} + {fixture['spacing']}*(j-1).
```

The telemetry stores the first `{fixture['telemetry_sum_count']}` values.  From
each of their two coefficient branches independently, the inverse above recovers
the same unique integer pivot.  The recovered roster starts at
`{fixture['first_rae']}` and ends at `{fixture['last_telemetry_rae']}`.

At 320-bit Arb precision:

* maximum saved-coefficient versus ideal-formula gap:
  `{aggregate['maximum_coefficient_replay_gap']['ball']['upper_decimal']}`;
* maximum recovered-`rae` error:
  `{aggregate['maximum_recovered_rae_error']['ball']['upper_decimal']}`;
* minimum distance of the source linear modulo choice from a half-integer:
  `{aggregate['minimum_linear_modulo_margin']['ball']['lower_decimal']}`;
* maximum exact difference between the two branch reconstructions of `xx`:
  `{aggregate['maximum_pair_xx_exact_gap']['decimal']}`.

The calculation was repeated at 192 and 320 bits; all 64 modulo lifts and all
32 recovered integer pivots were identical.

## Boundary

This gate proves a discrete physical-coordinate adapter for the 32 saved block-20
pivots.  It does **not** yet prove that an interval of physical `rae` values maps
inside any saved q-selector cell, does not cover the remaining 180 pivots of the
block, and does not certify the numerical values returned by the source `PSI` or
`ERF` routines.  Those are separate obligations, so this is not a Riemann
Hypothesis or prize-level conclusion.

The next falsifiable test is to interval-evaluate the physical coefficient curve
between adjacent roster points and compare its four coordinate widths with the
already certified q-cell radii.  A failed containment is evidence that subdivision
or a different coordinate enclosure is required; it must not be relabelled as
coverage.
"""


def main() -> int:
    for path in (SOURCE, TELEMETRY, CHECKPOINT, INPUT, CHECKER):
        require(path.is_file(), f"missing RAE-map dependency: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")

    headers = load_chain_headers()
    checkpoint_rows = load_checkpoint_rows()
    before = checkpoint_rows[(0, 19)]
    after = checkpoint_rows[(0, 20)]
    t = Fraction(before["t"])
    require(t == 10_000_000_000, "fixture t drift")
    rn1_start_fraction = Fraction(before["rn1"])
    require(rn1_start_fraction.denominator == 1, "block RN1 is not integral")
    rn1_start = rn1_start_fraction.numerator
    mt = int(after["mt"])
    aenums = int(after["aenums"])
    rn1_end = Fraction(after["rn1"])
    require(mt == 209 and aenums == 212, "block-20 fixture controls drift")
    require(
        rn1_end == rn1_start + 2 * (mt + 1) * aenums,
        "checkpoint block increment does not match source roster",
    )

    precision_runs = [evaluate_precision(bits, t, rn1_start, mt, headers) for bits in PRECISION_LADDER]
    require(
        precision_runs[0]["selectors"] == precision_runs[1]["selectors"],
        "linear modulo selectors changed across precision ladder",
    )
    require(
        precision_runs[0]["recovered_roster"] == precision_runs[1]["recovered_roster"],
        "recovered physical roster changed across precision ladder",
    )
    final = precision_runs[-1]
    first_rae = rn1_start + mt + 2
    spacing = 2 * (mt + 1)
    full_last_rae = first_rae + spacing * (aenums - 1)
    telemetry_last_rae = first_rae + spacing * 31

    artifact: dict[str, Any] = {
        "kind": "jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate",
        "status": "rigorous_fixture_physical_roster_inversion_with_continuous_coverage_open",
        "source": {
            "path": relative(SOURCE),
            "sha256": SOURCE_SHA256,
            "locations": source_locations(),
            "pi_provenance": {
                "ordinary_real_constant": "p=4*atan(1)",
                "physical_scale_constant": "PARI Pi2n(1)=2*pi, divided by 2",
                "scale": "a=sqrt(8*t/pi)",
            },
        },
        "exact_inverse_theorem": {
            "domain": ["r=rae/a>1", "pc>1", "xx>0"],
            "branch_cancellation": "xx=2*a2-3*a3 on both source branches",
            "pc_from_xx": "pc=1+1/xx",
            "pc_factorization": "pc=(r+sqrt(r^2-1))^2",
            "r_from_pc": "r=(sqrt(pc)+1/sqrt(pc))/2",
            "reason": "(r+s)(r-s)=1 for s=sqrt(r^2-1)",
            "exact_rational_fixture_count": 8,
            "exact_rational_fixtures": exact_symbolic_fixtures(),
        },
        "fixture": {
            "t": fraction_record(t),
            "block": 20,
            "pre_block_checkpoint_stage_unit": [0, 19],
            "post_block_checkpoint_stage_unit": [0, 20],
            "rn1_start": rn1_start,
            "rn1_end": int(rn1_end),
            "mt": mt,
            "aenums": aenums,
            "first_rae": first_rae,
            "spacing": spacing,
            "full_last_rae": full_last_rae,
            "telemetry_sum_count": 32,
            "last_telemetry_rae": telemetry_last_rae,
            "telemetry_branch_count": 64,
        },
        "precision_ladder": {
            "bits": list(PRECISION_LADDER),
            "stable_linear_modulo_lifts": True,
            "stable_recovered_roster": True,
            "coefficient_replay_tolerance": fraction_record(COEFFICIENT_REPLAY_TOLERANCE),
            "rae_recovery_tolerance": fraction_record(RAE_RECOVERY_TOLERANCE),
        },
        "physical_scale": {"pi": final["pi"], "a": final["a"]},
        "aggregate": {
            "pair_count": 32,
            "branch_count": 64,
            "recovered_roster": final["recovered_roster"],
            "maximum_coefficient_replay_gap": final["maximum_coefficient_replay_gap"],
            "maximum_recovered_rae_error": final["maximum_recovered_rae_error"],
            "minimum_linear_modulo_margin": final["minimum_linear_modulo_margin"],
            "maximum_pair_xx_exact_gap": final["maximum_pair_xx_exact_gap"],
        },
        "rows": final["rows"],
        "next_handoff": {
            "first": "Interval-evaluate the physical rae-to-(a1,xr,a3,fracL) curve between adjacent pivots and test q-cell containment.",
            "second": "Extend telemetry or an equivalent source-faithful adapter from 32 to all 212 block-20 pivots only after the first containment test is understood.",
            "falsification_rule": "Any interval coordinate wider than its target q-cell radius is a failed coverage claim and requires subdivision or a new enclosure.",
        },
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
            "input": {"path": relative(INPUT), "sha256": file_hash(INPUT)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "This is a rigorous discrete fixture adapter for 32 saved physical pivots. It does not prove "
            "continuous q-cell coverage, the remaining 180 block pivots, special-function accuracy, RH, "
            "or a prize-level theorem."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built Hardy coefficient-to-RAE roster gate: "
        f"{artifact['aggregate']['pair_count']} pivots, {artifact['aggregate']['branch_count']} branches, "
        f"spacing {spacing}, continuous coverage open"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
