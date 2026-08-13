#!/usr/bin/env python3
"""Rigorous four-branch decomposition of the nonsaddle endpoint error."""

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

import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate as full
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as pilot
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate.py"
)
WITNESS_CHAINS = pilot.WITNESS_CHAINS
PRECISIONS = pilot.PRECISIONS
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


def upper_abs(value: acb) -> Fraction:
    return upper(abs(value))


def linear_ray(data: dict[str, Any], endpoint: int, family: str, dps: int) -> acb:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    phi1, phi2, phi3 = (pilot.exact_ball(value) for value in parent["coefficients_hex"])
    x = arb(endpoint)
    phase = phi1 * x + phi2 * x**2 + phi3 * x**3
    slope = phi1 + 2 * phi2 * x + 3 * phi3 * x**2
    start = int(data["levels"][2]["length"]) + 1 if family == "b" else 1
    start_ball = arb(start)
    pi = arb.pi()
    if family == "b":
        endpoint_factor = (-2 * pi * I * phase).exp()
        digamma_sum = start_ball.digamma() - (start_ball - slope).digamma()
    else:
        endpoint_factor = (2 * pi * I * phase).exp()
        digamma_sum = (start_ball + slope).digamma() - start_ball.digamma()
    return I * endpoint_factor * digamma_sum / (2 * pi)


def decompose(
    data: dict[str, Any],
    sector_row: dict[str, Any],
    eq69_row: dict[str, Any],
    w2_row: dict[str, Any],
    dps: int,
    cutoff: Fraction,
) -> dict[str, acb]:
    evaluated = pilot.evaluate(
        data,
        sector_row,
        eq69_row,
        w2_row,
        dps,
        ray_cutoff=cutoff,
        central_breaks=full.breaks_for_cutoff(cutoff),
    )
    length = int(data["levels"][1]["length"])
    b_zero, b_length, c_zero, c_length = evaluated["rays"]
    b_linear_zero = linear_ray(data, 0, "b", dps)
    b_linear_length = linear_ray(data, length, "b", dps)
    c_linear_zero = linear_ray(data, 0, "c", dps)
    c_linear_length = linear_ray(data, length, "c", dps)
    b_linear_source = b_linear_zero - b_linear_length
    c_linear_source = (c_linear_length - c_linear_zero).conjugate()
    b_remainder = (b_zero["total"] - b_linear_zero) - (b_length["total"] - b_linear_length)
    c_remainder = (
        (c_length["total"] - c_linear_length) - (c_zero["total"] - c_linear_zero)
    ).conjugate()
    paper_w2 = native_q.complex_from_record(w2_row["paper_w2"])
    paper_w34 = native_q.complex_from_record(w2_row["paper_w3_or_w4"])
    paper_w5 = native_q.complex_from_record(w2_row["paper_w5"])
    paper_half = native_q.complex_from_record(w2_row["paper_half_sum"])
    linear_aggregate = (
        evaluated["endpoint_half"]
        + evaluated["a_endpoint"]
        + b_linear_source
        + c_linear_source
    )
    generic_w5_gap = linear_aggregate - (paper_half + paper_w5)
    correction_prediction = (
        b_remainder
        + c_remainder
        + evaluated["zero_mode"]
        - paper_w2
        - paper_w34
    )
    correction_decomposition_gap = correction_prediction - evaluated["exact_minus_paper"]
    return {
        "b_linear_source": b_linear_source,
        "c_linear_source": c_linear_source,
        "b_remainder": b_remainder,
        "c_remainder": c_remainder,
        "zero_mode": evaluated["zero_mode"],
        "paper_w2": paper_w2,
        "paper_w34": paper_w34,
        "linear_aggregate": linear_aggregate,
        "generic_w5_gap": generic_w5_gap,
        "correction_prediction": correction_prediction,
        "exact_correction": evaluated["exact_minus_paper"],
        "correction_decomposition_gap": correction_decomposition_gap,
    }


def build_note(artifact: dict[str, Any]) -> str:
    rows = "\n".join(
        f"chain {row['chain']:>3}: |generic-W5 gap| <= {row['generic_w5_gap_upper']}, "
        f"|correction gap| <= {row['correction_decomposition_gap_upper']}"
        for row in artifact["rows"]
    )
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 endpoint-linearization decomposition pilot

Date: 2026-08-07

Status: rigorous four-branch finite-input decomposition; not an all-height proof and not a proof of RH

## Exact decomposition

At endpoint `x`, write `A_x=F'(x)` and

```text
F(x+u)=F(x)+A_x*u+B_x*u^2+Phi3*u^3,
B_x=Phi2+3*Phi3*x.
```

The exact ray is compared with its endpoint-linear model, not with a uniform
quadratic model.  Summing the linear rays analytically gives

```text
J_b,lin(x)=i*E_-(x)/(2*pi) [psi(m)-psi(m-A_x)],
J_c,lin(x)=i*E_+(x)/(2*pi) [psi(1+A_x)-psi(1)].       (1)
```

The four linear rays plus the exact harmonic endpoint term equal the paper's
endpoint half-sum plus W5.  Therefore

```text
Q_exact-Q_paper
 = R_b+R_c+zero_mode-W2-(W3 or W4),                  (2)
```

where `R_b,R_c` are the full-ray minus linear-ray endpoint differences with
the source orientations retained.  W2 and W3 are the paper's exceptional
quadratic-mode replacements; W4 is its positive-Phi1 zero-mode replacement.

## Four-branch certificate

The 70- and 110-digit Arb constructions overlap throughout.  Both (1) and
(2) enclose equality on chains 1, 2, 3, and 36:

```text
{rows}
```

Aggregate bounds:

```text
maximum |generic linear aggregate-(half+W5)| <= {aggregate['maximum_generic_w5_gap_upper']}
maximum |decomposition prediction-correction| <= {aggregate['maximum_correction_decomposition_gap_upper']}
```

This identifies the correct analytic remainder architecture and rejects the
naive replacement of every endpoint ray by one quadratic logarithmic ray.
The next step is an all-374 replay and explicit majorants for `R_b`, `R_c`,
and the exceptional W2--W4 differences.  This pilot proves no height-uniform
bound, recursive accumulation, outer Hardy control, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion.
"""


def main() -> int:
    recursive = native_q.load_recursive_chains()
    sector = json.loads(pilot.SECTOR.read_text(encoding="utf-8"))
    eq69 = json.loads(pilot.EQ69.read_text(encoding="utf-8"))
    w2 = json.loads(pilot.W2_W5.read_text(encoding="utf-8"))
    sector_rows = {int(row["chain"]): row for row in sector["rows"]}
    eq69_rows = {int(row["chain"]): row for row in eq69["rows"]}
    w2_rows = {int(row["chain"]): row for row in w2["rows"]}

    rows: list[dict[str, Any]] = []
    max_generic = Fraction(0)
    max_correction = Fraction(0)
    for chain in WITNESS_CHAINS:
        cutoff, _ = full.select_cutoff(recursive[chain], sector_rows[chain])
        low = decompose(
            recursive[chain], sector_rows[chain], eq69_rows[chain], w2_rows[chain], PRECISIONS[0], cutoff
        )
        high = decompose(
            recursive[chain], sector_rows[chain], eq69_rows[chain], w2_rows[chain], PRECISIONS[1], cutoff
        )
        for name in high:
            require(low[name].overlaps(high[name]), f"chain {chain} decomposition {name} precision nonoverlap")
        require(high["generic_w5_gap"].contains(0), f"chain {chain} generic/W5 identity excludes zero")
        require(
            high["correction_decomposition_gap"].contains(0),
            f"chain {chain} correction decomposition excludes zero",
        )
        generic_upper = upper_abs(high["generic_w5_gap"])
        correction_upper = upper_abs(high["correction_decomposition_gap"])
        max_generic = max(max_generic, generic_upper)
        max_correction = max(max_correction, correction_upper)
        parent = recursive[chain]["levels"][1]
        phi1 = cells.binary128_fraction(parent["coefficients_hex"][0])
        phi3 = cells.binary128_fraction(parent["coefficients_hex"][2])
        rows.append(
            {
                "chain": chain,
                "phi1_sign": "positive" if phi1 > 0 else "negative",
                "phi3_sign": "positive" if phi3 > 0 else "negative",
                "ray_cutoff": int(cutoff),
                "b_remainder": point_balls.acb_record(high["b_remainder"], 45),
                "c_remainder": point_balls.acb_record(high["c_remainder"], 45),
                "zero_mode": point_balls.acb_record(high["zero_mode"], 45),
                "paper_w2": point_balls.acb_record(high["paper_w2"], 45),
                "paper_w3_or_w4": point_balls.acb_record(high["paper_w34"], 45),
                "correction_prediction": point_balls.acb_record(high["correction_prediction"], 45),
                "exact_correction": point_balls.acb_record(high["exact_correction"], 45),
                "generic_w5_gap_upper": decimal(generic_upper),
                "correction_decomposition_gap_upper": decimal(correction_upper),
            }
        )

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate",
        "status": "rigorous_four_branch_endpoint_linearization_decomposition_validated",
        "precisions_decimal_digits": list(PRECISIONS),
        "witness_chains": list(WITNESS_CHAINS),
        "identities": {
            "generic": "endpoint_half+a_endpoint+B_linear+C_linear=paper_endpoint_half+W5",
            "correction": "Q_exact-Q_paper=R_b+R_c+zero_mode-W2-(W3_or_W4)",
        },
        "rows": rows,
        "aggregate": {
            "branch_count": len(rows),
            "maximum_generic_w5_gap_upper": decimal(max_generic),
            "maximum_correction_decomposition_gap_upper": decimal(max_correction),
        },
        "dependencies": {
            "sector": {"path": relative(pilot.SECTOR), "sha256": file_hash(pilot.SECTOR)},
            "logarithmic_reduction": {"path": relative(pilot.LOG_REDUCTION), "sha256": file_hash(pilot.LOG_REDUCTION)},
            "full_arb": {"path": relative(full.RESULT), "sha256": file_hash(full.RESULT)},
            "w2_w5": {"path": relative(pilot.W2_W5), "sha256": file_hash(pilot.W2_W5)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Rigorous four-branch finite-input decomposition only. No all-call or height-uniform remainder bound, "
            "recurrence accumulation, outer Hardy control, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built endpoint-linearization decomposition pilot: 4/4 branches, both exact identities enclose zero")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
