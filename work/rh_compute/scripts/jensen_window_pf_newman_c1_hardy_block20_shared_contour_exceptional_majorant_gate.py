#!/usr/bin/env python3
"""Certify gap-uniform shared-contour majorants for W2/W3 exceptional modes."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate as full
import jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate as pilot
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as rays
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate.py"
)
I = acb(0, 1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper_abs(value: acb) -> Fraction:
    return upper(abs(value))


def lower_abs(value: acb) -> Fraction:
    return lower(abs(value))


def b_majorant(phi2: Fraction, phi3: Fraction, length: int) -> tuple[arb, Fraction, str]:
    """Uniform bound after the W2 linear/quadratic reciprocal cancellation."""
    local_quadratic = phi2 + 3 * phi3 * length
    require(phi2 > 0 and local_quadratic > 0, "nonpositive b homotopy quadratic coefficient")
    if phi3 >= 0:
        angle = "7*pi/12"
        sine_two_lower = Fraction(1, 2)
        quadratic_floor = phi2
    else:
        angle = "3*pi/4"
        sine_two_lower = Fraction(1)
        quadratic_floor = local_quadratic
    pi = arb.pi()
    coefficient = cells.fraction_ball(abs(phi3))
    floor_ball = cells.fraction_ball(quadratic_floor)
    decay = 2 * pi * floor_ball * cells.fraction_ball(sine_two_lower)
    bound = 2 * pi * coefficient * (
        3 * length * pi.sqrt() / (4 * decay.sqrt() ** 3) + 1 / (2 * decay**2)
    )
    return bound, quadratic_floor, angle


def c_majorant(phi2: Fraction, phi3: Fraction) -> tuple[arb, str]:
    """Uniform bound after the W3 linear/quadratic reciprocal cancellation."""
    require(phi2 > 0, "nonpositive c quadratic coefficient")
    angle = "pi/6" if phi3 >= 0 else "5*pi/12"
    pi = arb.pi()
    coefficient = cells.fraction_ball(abs(phi3))
    decay = pi * cells.fraction_ball(phi2)
    return pi * coefficient / decay**2, angle


def main() -> int:
    ctx.dps = 110
    ctx.threads = 1
    recursive = native_q.load_recursive_chains()
    sector = json.loads(rays.SECTOR.read_text(encoding="utf-8"))
    sector_rows = {int(row["chain"]): row for row in sector["rows"]}
    source = json.loads(full.RESULT.read_text(encoding="utf-8"))
    source_rows = {int(row["chain"]): row for row in source["rows"]}
    require(set(recursive) == set(sector_rows) == set(source_rows), "shared-contour roster drift")

    rows: list[dict[str, Any]] = []
    max_b_ratio = (Fraction(0), 0)
    max_c_ratio = (Fraction(0), 0)
    min_b_slack: tuple[Fraction, int] | None = None
    min_c_slack: tuple[Fraction, int] | None = None
    max_b_bound = (Fraction(0), 0)
    max_c_bound = (Fraction(0), 0)
    w3_count = 0
    for chain in sorted(recursive):
        data = recursive[chain]
        parent = data["levels"][1]
        phi1, phi2, phi3 = (cells.binary128_fraction(value) for value in parent["coefficients_hex"])
        length = int(parent["length"])
        b_bound_ball, b_floor, b_angle = b_majorant(phi2, phi3, length)
        b_bound_lower = lower(b_bound_ball)
        b_bound_upper = upper(b_bound_ball)
        b_actual = native_q.complex_from_record(source_rows[chain]["b_exact_minus_quadratic"])
        b_actual_lower = lower_abs(b_actual)
        b_actual_upper = upper_abs(b_actual)
        require(b_bound_lower > b_actual_upper, f"chain {chain} W2 shared-contour majorant failure")
        require(b_actual_lower > 0, f"chain {chain} W2 exact-minus-quadratic contains zero")
        b_ratio = b_bound_upper / b_actual_lower
        b_slack = b_bound_lower - b_actual_upper
        max_b_ratio = max(max_b_ratio, (b_ratio, chain))
        max_b_bound = max(max_b_bound, (b_bound_upper, chain))
        min_b_slack = (b_slack, chain) if min_b_slack is None else min(min_b_slack, (b_slack, chain))

        row: dict[str, Any] = {
            "chain": chain,
            "sum_index": int(data["header"]["sum_index"]),
            "branch": int(data["header"]["branch"]),
            "selector": source_rows[chain]["selector"],
            "phi3_sign": "nonnegative" if phi3 >= 0 else "negative",
            "b_shared_angle": b_angle,
            "b_quadratic_floor": decimal(b_floor),
            "b_majorant_lower": decimal(b_bound_lower),
            "b_majorant_upper": decimal(b_bound_upper),
            "b_actual_upper": decimal(b_actual_upper),
            "b_majorant_to_actual_ratio_upper": decimal(b_ratio),
            "b_majorant_slack_lower": decimal(b_slack),
        }
        if phi1 < 0:
            w3_count += 1
            cutoff = Fraction(int(source_rows[chain]["ray_cutoff"]))
            c_low, _ = pilot.exceptional_mode_integral(
                data, sector_rows[chain], 0, "c", 1, pilot.PRECISIONS[0], cutoff
            )
            c_high, c_tail = pilot.exceptional_mode_integral(
                data, sector_rows[chain], 0, "c", 1, pilot.PRECISIONS[1], cutoff
            )
            require(c_low.overlaps(c_high), f"chain {chain} W3 exact-mode precision nonoverlap")
            eta = 1 + cells.fraction_ball(phi1)
            pi = arb.pi()
            c_quadratic = I / (2 * pi) - pilot.paper_fresnel(cells.fraction_ball(phi2), eta)
            c_actual = c_high - c_quadratic
            c_actual_lower = lower_abs(c_actual)
            c_actual_upper = upper_abs(c_actual)
            c_bound_ball, c_angle = c_majorant(phi2, phi3)
            c_bound_lower = lower(c_bound_ball)
            c_bound_upper = upper(c_bound_ball)
            require(c_bound_lower > c_actual_upper, f"chain {chain} W3 shared-contour majorant failure")
            require(c_actual_lower > 0, f"chain {chain} W3 exact-minus-quadratic contains zero")
            c_ratio = c_bound_upper / c_actual_lower
            c_slack = c_bound_lower - c_actual_upper
            max_c_ratio = max(max_c_ratio, (c_ratio, chain))
            max_c_bound = max(max_c_bound, (c_bound_upper, chain))
            min_c_slack = (c_slack, chain) if min_c_slack is None else min(min_c_slack, (c_slack, chain))
            row.update(
                {
                    "c_shared_angle": c_angle,
                    "c_majorant_lower": decimal(c_bound_lower),
                    "c_majorant_upper": decimal(c_bound_upper),
                    "c_actual_upper": decimal(c_actual_upper),
                    "c_majorant_to_actual_ratio_upper": decimal(c_ratio),
                    "c_majorant_slack_lower": decimal(c_slack),
                    "c_exact_mode_tail_upper": decimal(c_tail),
                }
            )
        rows.append(row)

    require(min_b_slack is not None and min_c_slack is not None, "shared-contour slack aggregate missing")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate",
        "status": "rigorous_gap_uniform_W2_W3_shared_contour_majorants_validated_all_374",
        "theorem": {
            "b": (
                "|M_b,exact-M_b,quad| <= 2*pi*|Phi3|*"
                "[3*N*sqrt(pi)/(4*a_b^(3/2))+1/(2*a_b^2)], "
                "a_b=2*pi*B_*sigma_b, B_*=min(Phi2,Phi2+3*Phi3*N), "
                "sigma_b=1/2 for Phi3>=0 and 1 for Phi3<0"
            ),
            "c_W3": "|M_c,exact-M_c,quad| <= pi*|Phi3|/(pi*Phi2)^2",
            "gap_independence": "Neither majorant contains delta=ceil(F'(N))-F'(N) nor eta=1+Phi1.",
        },
        "rows": rows,
        "aggregate": {
            "w2_call_count": len(rows),
            "w3_call_count": w3_count,
            "w2_majorant_count": sum(Fraction(row["b_majorant_slack_lower"]) > 0 for row in rows),
            "w3_majorant_count": sum(
                Fraction(row["c_majorant_slack_lower"]) > 0 for row in rows if row["selector"] == "W3"
            ),
            "minimum_w2_majorant_slack_lower": decimal(min_b_slack[0]),
            "minimum_w2_majorant_slack_witness": min_b_slack[1],
            "minimum_w3_majorant_slack_lower": decimal(min_c_slack[0]),
            "minimum_w3_majorant_slack_witness": min_c_slack[1],
            "maximum_w2_majorant_upper": decimal(max_b_bound[0]),
            "maximum_w2_majorant_witness": max_b_bound[1],
            "maximum_w3_majorant_upper": decimal(max_c_bound[0]),
            "maximum_w3_majorant_witness": max_c_bound[1],
            "maximum_w2_majorant_to_actual_ratio_upper": decimal(max_b_ratio[0]),
            "maximum_w2_majorant_to_actual_ratio_witness": max_b_ratio[1],
            "maximum_w3_majorant_to_actual_ratio_upper": decimal(max_c_ratio[0]),
            "maximum_w3_majorant_to_actual_ratio_witness": max_c_ratio[1],
        },
        "dependencies": {
            "exceptional_full": {"path": relative(full.RESULT), "sha256": file_hash(full.RESULT)},
            "sector": {"path": relative(rays.SECTOR), "sha256": file_hash(rays.SECTOR)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "The symbolic shared-contour inequalities are gap-uniform and their finite coefficient roster "
            "application is rigorous. They cover the W2 exceptional mode and the W3 lower exceptional mode, "
            "not the W4 zero-mode pair, generic endpoint tails, recurrence accumulation, outer Hardy control, "
            "Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built shared-contour exceptional majorant gate: 374/374 W2 and 165/165 W3 inequalities certified")
    return 0


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Shared-contour exceptional-mode majorants

Date: 2026-08-07

Status: rigorous symbolic W2/W3 inequalities with an all-374 finite-roster application; not a complete height-uniform correction theorem and not a proof of RH

## Common-contour lemma

After the reciprocal terms have canceled, integration by parts writes the W2
exceptional difference as a difference of phase integrals.  Interpolate
between the paper quadratic phase and the exact cubic phase.  For `Phi3>=0`,
deform both integrals to `arg(u)=7*pi/12`; for `Phi3<0`, use
`arg(u)=3*pi/4`.  Along the entire homotopy the quadratic decay is at least
`B_* sigma_b r^2`, where

```text
B_* = min(Phi2, Phi2+3*Phi3*N),
sigma_b = 1/2  if Phi3>=0,
sigma_b = 1    if Phi3<0.
```

The phase derivative with respect to the homotopy parameter has magnitude at
most `2*pi*|Phi3|*(3*N*r^2+r^3)`.  Dropping the nonnegative linear and cubic
decay terms and integrating the two Gaussian moments gives

```text
|M_b,exact-M_b,quad|
 <= 2*pi*|Phi3| [3*N*sqrt(pi)/(4*a_b^(3/2)) + 1/(2*a_b^2)],
a_b = 2*pi*B_*sigma_b.                                 (1)
```

For W3 use `arg(u)=pi/6` when `Phi3>=0` and `arg(u)=5*pi/12`
when `Phi3<0`.  Since `sin(2*arg(u))>=1/2`, the same homotopy gives

```text
|M_c,exact-M_c,quad| <= pi*|Phi3|/(pi*Phi2)^2.          (2)
```

Here `pi` is forced by the Fourier normalization `exp(2*pi*i*n*u)`.  It is
not an inserted geometric fitting constant.  Neither (1) nor (2) contains
`delta=ceil(F'(N))-F'(N)` or `eta=1+Phi1`; both remain finite as either gap
tends to zero.

## Certified application

Arb verifies (1) on all {aggregate['w2_call_count']} W2 modes and (2) on all
{aggregate['w3_call_count']} W3 modes:

```text
minimum W2 slack  >= {aggregate['minimum_w2_majorant_slack_lower']}  (chain {aggregate['minimum_w2_majorant_slack_witness']})
minimum W3 slack  >= {aggregate['minimum_w3_majorant_slack_lower']}  (chain {aggregate['minimum_w3_majorant_slack_witness']})
maximum W2 bound  <= {aggregate['maximum_w2_majorant_upper']}
maximum W3 bound  <= {aggregate['maximum_w3_majorant_upper']}
worst W2 bound/actual ratio <= {aggregate['maximum_w2_majorant_to_actual_ratio_upper']}
worst W3 bound/actual ratio <= {aggregate['maximum_w3_majorant_to_actual_ratio_upper']}
```

The large worst ratios show that these are conservative theorem bounds, not
fitted approximations.  They solve the small-gap singularity problem for W2
and W3, but do not yet bound the W4 zero-mode pair or the coupled generic
endpoint tails.  They imply no recurrence budget, outer Hardy estimate,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


if __name__ == "__main__":
    raise SystemExit(main())
