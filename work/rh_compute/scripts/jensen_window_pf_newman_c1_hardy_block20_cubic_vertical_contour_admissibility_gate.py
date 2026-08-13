#!/usr/bin/env python3
"""Certify the cubic obstruction to the paper's vertical nonsaddle contours."""

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
CORRESPONDENCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate.py"
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


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def evaluate(data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    child = data["levels"][2]
    require(int(parent["degree"]) == 3 and int(parent["length"]) == 104, "cubic contour parent drift")
    a1, a2, a3 = (cells.binary128_fraction(value) for value in parent["coefficients_hex"])
    require(a3 != 0, "cubic contour zero Phi3")
    length = int(parent["length"])
    xi = a1 + 2 * a2 * length + 3 * a3 * length**2
    child_length = int(child["length"])
    require(xi.numerator // xi.denominator == child_length, "cubic contour xi selector drift")
    delta = child_length + 1 - xi
    require(0 < delta < 1, "cubic contour endpoint distance drift")

    if a3 < 0:
        family = "67b"
        phase = "g_n(z)=n*z-F(z)"
        endpoint = "z=N+iR"
        index = child_length + 1
        linear = delta
        threshold_squared = linear / (-a3)
        sign_formula = "Im g_n(N+iR)=R*(n-xi)+Phi3*R^3"
    else:
        family = "67c"
        phase = "h_n(z)=n*z+F(z)"
        endpoint = "z=iR"
        index = 1
        linear = 1 + a1
        require(linear > Fraction(1, 2), "cubic contour 67c linear drift")
        threshold_squared = linear / a3
        sign_formula = "Im h_1(iR)=R*(1+Phi1)-Phi3*R^3"

    threshold = cells.fraction_ball(threshold_squared).sqrt()
    witness_radius = 2 * threshold
    linear_ball = cells.fraction_ball(linear)
    a3_ball = cells.fraction_ball(a3)
    if family == "67b":
        imaginary_phase = witness_radius * linear_ball + a3_ball * witness_radius**3
    else:
        imaginary_phase = witness_radius * linear_ball - a3_ball * witness_radius**3
    identity_gap = imaginary_phase + 3 * linear_ball * witness_radius
    log_modulus = -2 * arb.pi() * imaginary_phase
    return {
        "family": family,
        "phase": phase,
        "endpoint": endpoint,
        "index": index,
        "phi3": a3,
        "linear": linear,
        "threshold_squared": threshold_squared,
        "threshold": threshold,
        "witness_radius": witness_radius,
        "imaginary_phase": imaginary_phase,
        "identity_gap": identity_gap,
        "log_modulus": log_modulus,
        "sign_formula": sign_formula,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 cubic vertical-contour admissibility gate

Date: 2026-08-07

Status: exact cubic far-field obstruction to the stated vertical nonsaddle contours; repair target, not a proof of RH

## Exact far-field calculation

Write

```text
F(z)=Phi1*z+Phi2*z^2+Phi3*z^3,   xi=F'(N).
```

For equation (82), the (67b) phase is `g_n(z)=n*z-F(z)`.  On the
paper's upper horizontal/vertical geometry,

```text
Im g_n(x+iR)=R*(n-F'(x))+Phi3*R^3.                    (1)
```

The bound asserted below equation (82) drops the cubic term.  If `Phi3<0`,
then at `x=N`, `n=ceil(xi)`, (1) becomes negative once
`R^2>(ceil(xi)-xi)/|Phi3|`; consequently
`|exp(2*pi*i*g_n)|=exp(-2*pi*Im g_n)` grows super-exponentially and the
vertical integrand in equation (83) does not even tend to zero.

For (67c), the phase in equation (90) is `h_n(z)=n*z+F(z)`, and

```text
Im h_n(x+iR)=R*(n+F'(x))-Phi3*R^3.                    (2)
```

If `Phi3>0`, then at `x=0`, `n=1`, (2) becomes negative once
`R^2>(1+Phi1)/Phi3`.  The stated vertical integral again cannot converge as
an ordinary improper integral.

## Exact roster

All 374 recursive block-20 cubics have nonzero `Phi3`:

```text
Phi3 < 0: {aggregate['negative_phi3_count']} calls; stated (67b) vertical contour fails
Phi3 > 0: {aggregate['positive_phi3_count']} calls; stated (67c) vertical contour fails
Phi3 = 0: {aggregate['zero_phi3_count']} calls
minimum crossover radius  >= {aggregate['minimum_crossover_radius_lower']}
maximum crossover radius  <= {aggregate['maximum_crossover_radius_upper']}
```

At the exact witness `R=2*R_cross`, both cases satisfy
`Im phase=-3*linear*R<0`.  The two independent 90/150-digit evaluations
overlap on every call, and the maximum identity gap is
`{aggregate['maximum_witness_identity_gap_upper']}`.

This does not show that the local W2--W4 formulas are unusable.  It shows
that their printed infinite vertical-contour derivation is incomplete for
every nonzero cubic on this roster.  In particular, the missing operation is
a contour deformation with a controlled connector and tail, not another
free fitted endpoint term.

## Repair target

The ray `z=z0+r*exp(i*pi/6)` is a valid cubic-only decay candidate for both
failing signs because `sin(3*pi/6)=1`.  Requiring the retained quadratic
endpoint model to decay as well sharpens the contract: use `5*pi/6` for the
negative-`Phi3` `(67b)` branch and `pi/6` for the positive-`Phi3` `(67c)`
branch, while the opposite-sign families retain `pi/2`.  The companion
sector-contour gate proves the resulting connector decay and summability;
the full-cubic versus quadratic endpoint remainder remains the next bound.

## Pi provenance and boundary

Pi appears only in the Fourier modulus identity
`|exp(2*pi*i*w)|=exp(-2*pi*Im(w))`; its positivity scales the growth or
decay and does not determine the cubic sign.  This gate identifies a missing
contour justification.  It does not yet supply the repaired contour, bound
the exact (67b)--(67c) aggregate, prove a height-uniform recurrence, control
the outer Hardy representation, or prove `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
"""


def main() -> int:
    for path in (PAPER, SOURCE, TELEMETRY, CORRESPONDENCE, CHECKER):
        require(path.is_file(), f"missing cubic contour dependency: {path}")
    correspondence = json.loads(CORRESPONDENCE.read_text(encoding="utf-8"))
    require(
        correspondence["status"] == "published_w2_w5_orientation_and_finite_component_residual_rigorously_enclosed",
        "cubic contour correspondence gate not admitted",
    )
    recursive = native_q.load_recursive_chains()
    rows: list[dict[str, Any]] = []
    signs: Counter[str] = Counter()
    minimum_threshold: tuple[Fraction, int] | None = None
    maximum_threshold: tuple[Fraction, int] | None = None
    maximum_identity: tuple[Fraction, int] | None = None
    maximum_imaginary_upper: tuple[Fraction, int] | None = None

    for chain in sorted(recursive):
        low = evaluate(recursive[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], PRECISIONS[1])
        require(low["family"] == high["family"], f"chain {chain} cubic contour family drift")
        for name in ("threshold", "witness_radius", "imaginary_phase", "identity_gap", "log_modulus"):
            require(low[name].overlaps(high[name]), f"chain {chain} cubic contour {name} precision nonoverlap")
        require(high["identity_gap"].contains(0), f"chain {chain} cubic contour witness identity drift")
        require(upper(high["imaginary_phase"]) < 0, f"chain {chain} cubic contour witness not growing")
        require(lower(high["log_modulus"]) > 0, f"chain {chain} cubic contour modulus witness drift")

        sign = "negative" if high["phi3"] < 0 else "positive"
        signs[sign] += 1
        threshold_lower = lower(high["threshold"])
        threshold_upper = upper(high["threshold"])
        identity_upper = cells.bound_fraction(abs(high["identity_gap"]).upper())
        imaginary_upper = upper(high["imaginary_phase"])
        if minimum_threshold is None or threshold_lower < minimum_threshold[0]:
            minimum_threshold = (threshold_lower, chain)
        if maximum_threshold is None or threshold_upper > maximum_threshold[0]:
            maximum_threshold = (threshold_upper, chain)
        if maximum_identity is None or identity_upper > maximum_identity[0]:
            maximum_identity = (identity_upper, chain)
        if maximum_imaginary_upper is None or imaginary_upper > maximum_imaginary_upper[0]:
            maximum_imaginary_upper = (imaginary_upper, chain)

        rows.append(
            {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "phi3_sign": sign,
                "failing_family": high["family"],
                "phase": high["phase"],
                "endpoint": high["endpoint"],
                "index": high["index"],
                "sign_formula": high["sign_formula"],
                "crossover_radius_squared_exact": {
                    "numerator": high["threshold_squared"].numerator,
                    "denominator": high["threshold_squared"].denominator,
                },
                "crossover_radius": point_balls.arb_record(high["threshold"], 45),
                "witness_radius": point_balls.arb_record(high["witness_radius"], 45),
                "witness_imaginary_phase": point_balls.arb_record(high["imaginary_phase"], 45),
                "witness_log_modulus": point_balls.arb_record(high["log_modulus"], 45),
                "witness_identity_gap_upper": decimal(identity_upper),
                "precision_overlap": True,
            }
        )

    require(signs == Counter({"negative": 187, "positive": 187}), "cubic contour sign roster drift")
    require(minimum_threshold and maximum_threshold and maximum_identity and maximum_imaginary_upper, "cubic contour aggregate missing")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate",
        "status": "exact_cubic_far_field_obstruction_to_stated_vertical_nonsaddle_contours",
        "scope": {
            "block": 20,
            "recursive_call_count": 374,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
        },
        "identities": {
            "67b": "Im(n*(x+iR)-F(x+iR))=R*(n-F'(x))+Phi3*R^3",
            "67c": "Im(n*(x+iR)+F(x+iR))=R*(n+F'(x))-Phi3*R^3",
            "modulus": "|exp(2*pi*i*w)|=exp(-2*pi*Im(w))",
            "witness": "At R=2*R_cross, Im(phase)=-3*linear*R<0",
        },
        "aggregate": {
            "negative_phi3_count": signs["negative"],
            "positive_phi3_count": signs["positive"],
            "zero_phi3_count": 0,
            "failing_67b_count": signs["negative"],
            "failing_67c_count": signs["positive"],
            "minimum_crossover_radius_lower": decimal(minimum_threshold[0]),
            "minimum_crossover_radius_witness": minimum_threshold[1],
            "maximum_crossover_radius_upper": decimal(maximum_threshold[0]),
            "maximum_crossover_radius_witness": maximum_threshold[1],
            "maximum_witness_identity_gap_upper": decimal(maximum_identity[0]),
            "maximum_witness_identity_gap_witness": maximum_identity[1],
            "maximum_witness_imaginary_phase_upper": decimal(maximum_imaginary_upper[0]),
            "maximum_witness_imaginary_phase_witness": maximum_imaginary_upper[1],
        },
        "rows": rows,
        "paper": {
            "path": relative(PAPER),
            "sha256": file_hash(PAPER),
            "pages": [27, 28, 29, 30],
            "equations": [82, 83, 84, 85, 87, 90, 91],
        },
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "correspondence_gate": {"path": relative(CORRESPONDENCE), "sha256": file_hash(CORRESPONDENCE)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "repair_target": {
            "candidate_ray_angle": "pi/6",
            "leading_decay": "sin(3*pi/6)=1",
            "obligations": [
                "construct a vertical-to-pi/6 homotopy for each failing family",
                "prove no saddle or Stokes contribution is crossed",
                "bound the connector and far tail uniformly in n",
                "sum the repaired 67b and 67c tails before comparing with W2--W4",
            ],
        },
        "proof_boundary": (
            "This proves a far-field obstruction to the stated vertical contours on 374 exact cubic inputs. It does not "
            "supply the repaired contour, bound the exact nonsaddle aggregate, prove a height-uniform recurrence, control "
            "the outer Hardy representation, or prove Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built cubic contour admissibility gate: 187 failing 67b, 187 failing 67c, zero admissible paired vertical cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
