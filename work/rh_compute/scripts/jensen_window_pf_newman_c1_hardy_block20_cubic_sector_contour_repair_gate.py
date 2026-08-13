#!/usr/bin/env python3
"""Certify sign-aware decay sectors for the cubic nonsaddle contours."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import arb, ctx


PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
ADMISSIBILITY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate.json"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.py"
)
PRECISIONS = (90, 150)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def fraction_record(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": decimal(value),
    }


def evaluate(data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    child = data["levels"][2]
    require(int(parent["degree"]) == 3 and int(parent["length"]) == 104, "sector contour parent drift")
    phi1, phi2, phi3 = (cells.binary128_fraction(value) for value in parent["coefficients_hex"])
    require(phi3 != 0, "sector contour zero Phi3")
    length = int(parent["length"])
    xi = phi1 + 2 * phi2 * length + 3 * phi3 * length**2
    child_length = int(child["length"])
    require(xi.numerator // xi.denominator == child_length, "sector contour xi selector drift")
    delta = child_length + 1 - xi
    require(0 < delta < 1, "sector contour endpoint gap drift")

    second_zero = 2 * phi2
    second_length = 2 * phi2 + 6 * phi3 * length
    second_minimum = min(second_zero, second_length)
    require(second_minimum > 0, "sector contour lost strict convexity")
    c_gap = 1 + phi1
    require(c_gap > Fraction(1, 2), "sector contour 67c derivative gap drift")

    slanted_quadratic = arb(3).sqrt() * cells.fraction_ball(second_minimum) / 4
    if phi3 < 0:
        b_angle = "5*pi/6"
        b_linear = delta / 2
        b_quadratic = slanted_quadratic
        b_cubic = -phi3
        c_angle = "pi/2"
        c_linear = c_gap
        c_quadratic = arb(0)
        c_cubic = -phi3
    else:
        b_angle = "pi/2"
        b_linear = delta
        b_quadratic = arb(0)
        b_cubic = phi3
        c_angle = "pi/6"
        c_linear = c_gap / 2
        c_quadratic = slanted_quadratic
        c_cubic = phi3

    for name, value in (
        ("67b linear", b_linear),
        ("67b cubic", b_cubic),
        ("67c linear", c_linear),
        ("67c cubic", c_cubic),
    ):
        require(value > 0, f"sector contour {name} coefficient drift")
    require(lower(slanted_quadratic) > 0, "sector contour quadratic decay not positive")
    return {
        "phi1": phi1,
        "phi2": phi2,
        "phi3": phi3,
        "length": length,
        "xi": xi,
        "delta": delta,
        "second_zero": second_zero,
        "second_length": second_length,
        "second_minimum": second_minimum,
        "c_gap": c_gap,
        "slanted_quadratic": slanted_quadratic,
        "b_angle": b_angle,
        "b_linear": b_linear,
        "b_quadratic": b_quadratic,
        "b_cubic": b_cubic,
        "c_angle": c_angle,
        "c_linear": c_linear,
        "c_quadratic": c_quadratic,
        "c_cubic": c_cubic,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 cubic sector-contour repair gate

Date: 2026-08-07

Status: exact sign-aware contour lemma on the saved cubic roster; not a proof of the full recurrence or RH

## Sector choice

Let `F(z)=Phi1*z+Phi2*z^2+Phi3*z^3`, with `F''>0` on `[0,N]`.
For `(67b)` put `g_n(z)=n*z-F(z)` and for `(67c)` put
`h_n(z)=n*z+F(z)`.  The sign-aware rays are

```text
                         Phi3 < 0       Phi3 > 0
(67b), g_n                 5*pi/6          pi/2
(67c), h_n                   pi/2          pi/6
```

For the two slanted cases, exact cubic Taylor expansion at any
`x in [0,N]` gives

```text
Im[g_n(x+r*exp(5*pi*i/6))-g_n(x)]
 = (n-F'(x))*r/2 + sqrt(3)*F''(x)*r^2/4 + |Phi3|*r^3,  Phi3<0,

Im[h_n(x+r*exp(pi*i/6))-h_n(x)]
 = (n+F'(x))*r/2 + sqrt(3)*F''(x)*r^2/4 + Phi3*r^3,    Phi3>0.
```

On the retained vertical rays the corresponding identities are
`(n-F'(x))*r+Phi3*r^3` for positive `Phi3` in `(67b)` and
`(n+F'(x))*r+|Phi3|*r^3` for negative `Phi3` in `(67c)`.  Every displayed
coefficient is positive.  Thus the translated far connector is bounded by
a quadratic polynomial in `r` times `exp(-2*pi*|Phi3|*r^3)` and tends to
zero.

## Exact roster checks

All 374 recursive calls satisfy the hypotheses:

```text
strict F''>0 calls                         {aggregate['strict_convexity_count']}
negative / positive Phi3 calls            {aggregate['negative_phi3_count']} / {aggregate['positive_phi3_count']}
67b slanted / vertical contracts          {aggregate['b_slanted_count']} / {aggregate['b_vertical_count']}
67c slanted / vertical contracts          {aggregate['c_slanted_count']} / {aggregate['c_vertical_count']}
minimum F'' on [0,N]                      >= {aggregate['minimum_second_derivative']}
minimum 67b linear decay coefficient      >= {aggregate['minimum_b_linear_decay']}
minimum 67c linear decay coefficient      >= {aggregate['minimum_c_linear_decay']}
minimum cubic decay coefficient           >= {aggregate['minimum_cubic_decay']}
minimum slanted quadratic coefficient     >= {aggregate['minimum_slanted_quadratic_decay_lower']}
```

The 90- and 150-digit `sqrt(3)` interval evaluations overlap on every row.
Since the integrands are entire, the finite parallelogram deformation crosses
no singularity and creates no residue.  A saddle is not a singularity; Stokes
switching concerns a later asymptotic decomposition, not this exact contour
identity.

## Quadratic rotation and summation

The endpoint quadratic model is compatible with the paper's vertical
error-function contour.  For `(67b)`, rotate from `5*pi/6` to `pi/2`:
`sin(theta)>=1/2` and `-sin(2*theta)>=0`.  For `(67c)`, rotate from `pi/6`
to `pi/2`: `sin(theta)>=1/2` and `sin(2*theta)>=0`.  Hence the large
quadratic-model arc decays and the ray rotation is exact.

The outer `1/n` sums may also be interchanged with the repaired endpoint
rays on this finite roster.  After taking absolute values, their common
kernel is bounded near `r=0` by a constant times
`1+|log(r)|`, because `sum(exp(-pi*n*r)/n)=-log(1-exp(-pi*r))`; at infinity
the certified cubic factor dominates every polynomial amplitude.  The
logarithmic singularity is integrable.

## Pi provenance and remaining boundary

Here `pi` comes from the original Fourier exponential
`exp(2*pi*i*phase)`.  The angles are selected by the exact trigonometric
values at `pi/6`, `pi/2`, and `5*pi/6`; no circle-derived numerical fit is
inserted.

This repairs the existence of the exact nonsaddle contours and preserves the
paper's quadratic endpoint-ray value.  It does not bound the difference
between the full cubic endpoint integrals and their quadratic models, justify
the other amplitude truncations in W2--W4, establish a height-uniform
recurrence, control the outer Hardy remainder, or prove `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    for path in (PAPER, SOURCE, TELEMETRY, ADMISSIBILITY, CHECKER):
        require(path.is_file(), f"missing sector contour dependency: {path}")
    admissibility = json.loads(ADMISSIBILITY.read_text(encoding="utf-8"))
    require(
        admissibility["status"] == "exact_cubic_far_field_obstruction_to_stated_vertical_nonsaddle_contours",
        "sector contour obstruction gate not admitted",
    )

    recursive = native_q.load_recursive_chains()
    rows: list[dict[str, Any]] = []
    signs: Counter[str] = Counter()
    angle_counts: Counter[str] = Counter()
    minima: dict[str, tuple[Fraction, int] | None] = {
        "second": None,
        "b_linear": None,
        "c_linear": None,
        "cubic": None,
        "slanted_quadratic": None,
    }

    def retain_minimum(name: str, value: Fraction, chain: int) -> None:
        current = minima[name]
        if current is None or value < current[0]:
            minima[name] = (value, chain)

    for chain in sorted(recursive):
        low = evaluate(recursive[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], PRECISIONS[1])
        for name in ("slanted_quadratic", "b_quadratic", "c_quadratic"):
            require(low[name].overlaps(high[name]), f"chain {chain} sector contour {name} precision nonoverlap")
        for name in (
            "phi1",
            "phi2",
            "phi3",
            "length",
            "xi",
            "delta",
            "second_zero",
            "second_length",
            "second_minimum",
            "c_gap",
            "b_angle",
            "b_linear",
            "b_cubic",
            "c_angle",
            "c_linear",
            "c_cubic",
        ):
            require(low[name] == high[name], f"chain {chain} sector contour exact data drift: {name}")

        sign = "negative" if high["phi3"] < 0 else "positive"
        signs[sign] += 1
        angle_counts[f"b:{high['b_angle']}"] += 1
        angle_counts[f"c:{high['c_angle']}"] += 1
        retain_minimum("second", high["second_minimum"], chain)
        retain_minimum("b_linear", high["b_linear"], chain)
        retain_minimum("c_linear", high["c_linear"], chain)
        retain_minimum("cubic", abs(high["phi3"]), chain)
        retain_minimum("slanted_quadratic", lower(high["slanted_quadratic"]), chain)

        rows.append(
            {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "phi3_sign": sign,
                "second_derivative_at_zero": fraction_record(high["second_zero"]),
                "second_derivative_at_length": fraction_record(high["second_length"]),
                "second_derivative_minimum": fraction_record(high["second_minimum"]),
                "endpoint_delta": fraction_record(high["delta"]),
                "c_family_gap": fraction_record(high["c_gap"]),
                "b_contract": {
                    "angle": high["b_angle"],
                    "linear_decay": fraction_record(high["b_linear"]),
                    "quadratic_decay": point_balls.arb_record(high["b_quadratic"], 45),
                    "cubic_decay": fraction_record(high["b_cubic"]),
                },
                "c_contract": {
                    "angle": high["c_angle"],
                    "linear_decay": fraction_record(high["c_linear"]),
                    "quadratic_decay": point_balls.arb_record(high["c_quadratic"], 45),
                    "cubic_decay": fraction_record(high["c_cubic"]),
                },
                "precision_overlap": True,
            }
        )

    require(signs == Counter({"negative": 187, "positive": 187}), "sector contour sign roster drift")
    require(angle_counts["b:5*pi/6"] == angle_counts["b:pi/2"] == 187, "sector contour 67b angle drift")
    require(angle_counts["c:pi/6"] == angle_counts["c:pi/2"] == 187, "sector contour 67c angle drift")
    require(all(value is not None for value in minima.values()), "sector contour aggregate minimum missing")

    second = minima["second"]
    b_linear = minima["b_linear"]
    c_linear = minima["c_linear"]
    cubic = minima["cubic"]
    slanted_quadratic = minima["slanted_quadratic"]
    assert second and b_linear and c_linear and cubic and slanted_quadratic
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate",
        "status": "exact_sign_aware_cubic_nonsaddle_contour_existence_and_rotation_lemma",
        "scope": {
            "block": 20,
            "recursive_call_count": 374,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
        },
        "ray_contracts": {
            "67b_phi3_negative": {
                "angle": "5*pi/6",
                "identity": "Im(g_n(x+r*e^(5*pi*i/6))-g_n(x))=(n-F'(x))*r/2+sqrt(3)*F''(x)*r^2/4+abs(Phi3)*r^3",
            },
            "67b_phi3_positive": {
                "angle": "pi/2",
                "identity": "Im(g_n(x+i*r)-g_n(x))=(n-F'(x))*r+Phi3*r^3",
            },
            "67c_phi3_negative": {
                "angle": "pi/2",
                "identity": "Im(h_n(x+i*r)-h_n(x))=(n+F'(x))*r+abs(Phi3)*r^3",
            },
            "67c_phi3_positive": {
                "angle": "pi/6",
                "identity": "Im(h_n(x+r*e^(pi*i/6))-h_n(x))=(n+F'(x))*r/2+sqrt(3)*F''(x)*r^2/4+Phi3*r^3",
            },
        },
        "contour_lemma": {
            "analyticity": "The cubic phases and polynomial amplitudes are entire, so the finite parallelogram deformation has no residues.",
            "connector_bound": "polynomial(r)*exp(-2*pi*abs(Phi3)*r^3), uniformly in x and no worse for larger n",
            "quadratic_rotation_67b": "theta in [pi/2,5*pi/6]: sin(theta)>=1/2 and -sin(2*theta)>=0",
            "quadratic_rotation_67c": "theta in [pi/6,pi/2]: sin(theta)>=1/2 and sin(2*theta)>=0",
            "sum_kernel": "sum_{n>=1} exp(-pi*n*r)/n=-log(1-exp(-pi*r))=O(1+abs(log(r))) near zero",
            "sum_interchange": "absolute on each saved call because the logarithmic origin singularity is integrable and cubic decay controls infinity",
        },
        "aggregate": {
            "strict_convexity_count": 374,
            "negative_phi3_count": signs["negative"],
            "positive_phi3_count": signs["positive"],
            "b_slanted_count": angle_counts["b:5*pi/6"],
            "b_vertical_count": angle_counts["b:pi/2"],
            "c_slanted_count": angle_counts["c:pi/6"],
            "c_vertical_count": angle_counts["c:pi/2"],
            "minimum_second_derivative": decimal(second[0]),
            "minimum_second_derivative_witness": second[1],
            "minimum_b_linear_decay": decimal(b_linear[0]),
            "minimum_b_linear_decay_witness": b_linear[1],
            "minimum_c_linear_decay": decimal(c_linear[0]),
            "minimum_c_linear_decay_witness": c_linear[1],
            "minimum_cubic_decay": decimal(cubic[0]),
            "minimum_cubic_decay_witness": cubic[1],
            "minimum_slanted_quadratic_decay_lower": decimal(slanted_quadratic[0]),
            "minimum_slanted_quadratic_decay_witness": slanted_quadratic[1],
        },
        "rows": rows,
        "paper": {
            "path": relative(PAPER),
            "sha256": file_hash(PAPER),
            "pages": [23, 27, 28, 29, 30],
            "equations": [65, 67, 82, 83, 86, 89, 90, 91, 93, 94],
        },
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "admissibility_gate": {"path": relative(ADMISSIBILITY), "sha256": file_hash(ADMISSIBILITY)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Enclose the full-cubic endpoint ray integrals against the rotated quadratic models, including the exact "
            "F''(N)/2 versus Phi2 and polynomial-amplitude differences, before claiming W2--W4 error control."
        ),
        "proof_boundary": (
            "This proves exact sign-aware contour existence, connector decay, quadratic-ray rotation, and termwise "
            "summability on 374 saved cubic calls. It does not bound the full-cubic versus quadratic endpoint remainder, "
            "prove a height-uniform recurrence, control the outer Hardy representation, or prove Lambda<=0, PF-infinity, "
            "RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built cubic sector-contour repair gate: 374 strict-convexity calls, 187+187 sign-aware slanted rays")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
