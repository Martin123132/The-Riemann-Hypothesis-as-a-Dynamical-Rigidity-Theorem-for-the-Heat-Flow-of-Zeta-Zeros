#!/usr/bin/env python3
"""Certify the pole-free deformation induced by the diagnostic cutoff."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate as partition


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation62_pole_free_radius_saddle_collar_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation62_pole_free_radius_saddle_collar_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
LEGACY_PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_2015_1502.06903.pdf"

PRECISIONS = (70, 110)
EXPECTED_OUTPUTS = 15
TERMINAL_ALPHA = 159777
NEXT_ODD_ALPHA = 159779
SOURCE_TAIL_START = 37946
CELL_TAIL_START = 37941
COLLAR = tuple(range(CELL_TAIL_START, SOURCE_TAIL_START))
LOWER_NEIGHBOR = CELL_TAIL_START - 1
UPPER_NEIGHBOR = SOURCE_TAIL_START
FAR_ODD_BELOW = 6366037945
FAR_ODD_ABOVE = FAR_ODD_BELOW + 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def alpha_of_n(t: arb, n: int) -> arb:
    return 2 * arb(n) + t / (arb.pi() * n)


def lower_root(t: arb, alpha: arb) -> arb:
    saddle = (8 * t / arb.pi()).sqrt()
    require(alpha > saddle, "alpha must lie above the central saddle")
    return alpha * (1 - (1 - (saddle / alpha) ** 2).sqrt()) / 4


def u_saddle(t: arb, n: int) -> arb:
    """Appendix-B equation (B15a)."""
    return 1 / (1 - 1 / arb(n) - 2 * arb.pi() * n / t)


def radius_row(target_t: str, precision: int) -> dict[str, Any]:
    ctx.dps = precision
    t = arb(target_t)
    pi = arb.pi()
    alpha = {n: alpha_of_n(t, n) for n in range(LOWER_NEIGHBOR, UPPER_NEIGHBOR + 1)}
    require(alpha[UPPER_NEIGHBOR] < TERMINAL_ALPHA < alpha[UPPER_NEIGHBOR - 1], "source threshold ordering drift")
    require(TERMINAL_ALPHA < alpha[COLLAR[-1]], "first collar threshold not above terminal pole")
    for n in range(LOWER_NEIGHBOR, UPPER_NEIGHBOR):
        require(alpha[n] > alpha[n + 1], f"alpha(n) monotonicity drift at n={n}")
    require(alpha[LOWER_NEIGHBOR] < NEXT_ODD_ALPHA, "cell selector interval leaves the left pole gap")

    # These selectors are interior points of the two open threshold intervals.
    source_left = (arb(TERMINAL_ALPHA) + alpha[COLLAR[-1]]) / 2
    cell_left = (alpha[CELL_TAIL_START] + alpha[LOWER_NEIGHBOR]) / 2
    require(TERMINAL_ALPHA < source_left < alpha[COLLAR[-1]], "source selector interval failed")
    require(alpha[CELL_TAIL_START] < cell_left < alpha[LOWER_NEIGHBOR], "cell selector interval failed")
    require(source_left < cell_left < NEXT_ODD_ALPHA, "selectors leave the common left pole gap")

    source_root = lower_root(t, source_left)
    cell_root = lower_root(t, cell_left)
    require(SOURCE_TAIL_START - 1 < source_root < SOURCE_TAIL_START, "source selector tail start drift")
    require(CELL_TAIL_START - 1 < cell_root < CELL_TAIL_START, "cell selector tail start drift")

    source_radius = t / pi - source_left
    cell_radius = t / pi - cell_left
    require(0 < cell_radius < source_radius, "radius ordering drift")
    a = (8 * t / pi).sqrt()
    plus_branch_distance = t / pi - a
    minus_branch_distance = t / pi + a
    require(source_radius < plus_branch_distance < minus_branch_distance, "branch point enters radial shell")
    source_u0 = t / (pi * source_radius)
    cell_u0 = t / (pi * cell_radius)
    require(1 < source_u0 < cell_u0, "D lower-limit ordering drift")

    source_far = 2 * t / pi - source_left
    cell_far = 2 * t / pi - cell_left
    require(FAR_ODD_BELOW < cell_far < source_far < FAR_ODD_ABOVE, "far-side pole-free interval failed")

    saddle_rows: list[dict[str, Any]] = []
    for n in range(LOWER_NEIGHBOR, UPPER_NEIGHBOR + 1):
        us = u_saddle(t, n)
        inside_source = bool(us > source_u0)
        inside_cell = bool(us > cell_u0)
        if n == LOWER_NEIGHBOR:
            require(inside_source and inside_cell, "lower neighbor must remain inside both D contours")
            classification = "inside_both"
        elif n in COLLAR:
            require(inside_source and not inside_cell, f"collar saddle membership drift at n={n}")
            classification = "crossed_by_deformation"
        else:
            require(not inside_source and not inside_cell, "upper neighbor must remain outside both D contours")
            classification = "outside_both"
        # This is the exact B15a/124 threshold equivalence, not a numerical heuristic.
        threshold_identity = (1 / us) - (1 - pi * alpha[n] / t)
        require(threshold_identity.contains(0), f"saddle/alpha threshold identity failed at n={n}")
        saddle_rows.append(
            {
                "n": n,
                "alpha_of_n_ball": alpha[n].str(precision, more=True),
                "u_saddle_ball": us.str(precision, more=True),
                "inside_source_selector_D_range": inside_source,
                "inside_cell_selector_D_range": inside_cell,
                "classification": classification,
                "b15a_equation124_threshold_identity_ball": threshold_identity.str(precision, more=True),
            }
        )

    return {
        "target_t": target_t,
        "terminal_alpha": TERMINAL_ALPHA,
        "next_odd_alpha": NEXT_ODD_ALPHA,
        "source_selector": {
            "left_intercept_ball": source_left.str(precision, more=True),
            "radius_ball": source_radius.str(precision, more=True),
            "D_lower_limit_u0_ball": source_u0.str(precision, more=True),
            "continuous_lower_root_ball": source_root.str(precision, more=True),
            "classical_tail_start": SOURCE_TAIL_START,
            "far_intercept_ball": source_far.str(precision, more=True),
        },
        "cell_selector": {
            "left_intercept_ball": cell_left.str(precision, more=True),
            "radius_ball": cell_radius.str(precision, more=True),
            "D_lower_limit_u0_ball": cell_u0.str(precision, more=True),
            "continuous_lower_root_ball": cell_root.str(precision, more=True),
            "classical_tail_start": CELL_TAIL_START,
            "far_intercept_ball": cell_far.str(precision, more=True),
        },
        "far_pole_gap": {
            "odd_below": FAR_ODD_BELOW,
            "odd_above": FAR_ODD_ABOVE,
            "certificate": "odd_below < cell_far < source_far < odd_above",
        },
        "branch_shell": {
            "a_ball": a.str(precision, more=True),
            "plus_a_radial_distance_ball": plus_branch_distance.str(precision, more=True),
            "minus_a_radial_distance_ball": minus_branch_distance.str(precision, more=True),
            "certificate": "cell_radius < source_radius < distance(center,+a) < distance(center,-a)",
        },
        "saddle_rows": saddle_rows,
        "crossed_saddles": [row["n"] for row in saddle_rows if row["classification"] == "crossed_by_deformation"],
    }


def render_note(artifact: dict[str, Any]) -> str:
    center = artifact["output_rows"][7]
    source = center["source_selector"]
    cell = center["cell_selector"]
    return f"""# Pole-free radius and five-saddle collar gate

Date: 2026-08-09

Status: exact diagnostic contour-deformation target isolated; not a proof of RH

Lewis (2015), equation (8), uses a circle centred at `t/pi`.  Its left and right real
intercepts are `L=t/pi-R` and `U=2t/pi-L`.  The paper requires the radius to be
adjusted whenever an intercept would hit an odd-integer pole.

At the central output this gate chooses two explicit admissible left
intercepts inside the same gap `(159777,159779)`:

```text
source-tail selector L_s = {source['left_intercept_ball']}
diagnostic selector  L_c = {cell['left_intercept_ball']}
```

The equation-(124) root relation gives tail starts `37946` and `37941`,
respectively.  The first is the published/source endpoint convention; the
second is the introduced diagnostic midpoint convention.  On the far
side both corresponding intercepts lie strictly between the same consecutive
odd poles `{FAR_ODD_BELOW}` and `{FAR_ODD_ABOVE}`.  Both radii are also strictly
smaller than the radial distance to the branch point `z=+a`, with `z=-a`
farther away.  Hence the radial annulus contains neither an odd pole nor either
denominator branch point, and the enclosed alpha-pole roster is unchanged.

Appendix-B equation (B9) has lower limit

```text
u_0(R)=t/(pi R),
```

while its saddle is exactly, by (B15a),

```text
u_sad(N)=1/(1-1/N-2pi N/t).
```

The identity

```text
1/u_sad(N)=1-pi[2N+t/(pi N)]/t
```

shows that saddle crossing is governed by the same equation-(124) threshold.
At all fifteen outputs, `N=37940` remains inside both `D` ranges,
`N=37946` remains outside both, and exactly

```text
N=37941,37942,37943,37944,37945
```

cross the moving lower limit.  Therefore

```text
D(N,R_s)-D(N,R_c)
 = integral_[u_0(R_s)]^[u_0(R_c)] f_D(N,u) du
```

is an exact finite-interval identity for each collar mode, before saddle
approximation.

This does not permit adding the five classical main terms while keeping the
old error unchanged.  Because the two circles enclose the same poles, Cauchy
deformation makes the complete contour invariant.  If `A(R)` denotes the
first theta integral in equation (13), and

```text
S(R)=C(0,R)+D(0,R)+2 sum_(N>=1)(-1)^N[C(N,R)+D(N,R)],
```

then equations (13) and (24) give the exact compensation law

```text
R_s A(R_s)-i t S(R_s)/pi
 =R_c A(R_c)-i t S(R_c)/pi.
```

The diagnostic midpoint is therefore a mathematically defined finite
re-splitting candidate, not a paper-prescribed or proved sharper approximation.
Its transported `C`, `D`, first-integral, endpoint, and saddle remainders would
have to be derived as a new theorem.  This gate proves only the exact geometry;
it supplies no height-uniform remainder and no RH implication.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    low_rows: list[dict[str, Any]] = []
    high_rows: list[dict[str, Any]] = []
    for output_index in range(1, EXPECTED_OUTPUTS + 1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        low_rows.append(radius_row(target_t, PRECISIONS[0]))
        high_rows.append(radius_row(target_t, PRECISIONS[1]))

    output_rows: list[dict[str, Any]] = []
    for output_index, (low, high) in enumerate(zip(low_rows, high_rows), start=1):
        require(low["target_t"] == high["target_t"], "target precision roster drift")
        for selector in ("source_selector", "cell_selector"):
            require(low[selector]["classical_tail_start"] == high[selector]["classical_tail_start"], "selector tail precision drift")
            for key in ("left_intercept_ball", "radius_ball", "D_lower_limit_u0_ball", "continuous_lower_root_ball", "far_intercept_ball"):
                require(arb(low[selector][key]).overlaps(arb(high[selector][key])), f"selector precision drift: {selector}/{key}")
        require(low["crossed_saddles"] == high["crossed_saddles"] == list(COLLAR), "crossed-saddle precision drift")
        for low_saddle, high_saddle in zip(low["saddle_rows"], high["saddle_rows"]):
            for key in ("n", "inside_source_selector_D_range", "inside_cell_selector_D_range", "classification"):
                require(low_saddle[key] == high_saddle[key], f"saddle discrete precision drift: {key}")
            for key in ("alpha_of_n_ball", "u_saddle_ball", "b15a_equation124_threshold_identity_ball"):
                require(arb(low_saddle[key]).overlaps(arb(high_saddle[key])), f"saddle ball precision drift: {key}")
        output_rows.append({"output_index": output_index, "output_label": f"t{Decimal(output_index - 8) / Decimal(100):+.2f}", **high})

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation62_pole_free_radius_saddle_collar_gate",
        "status": "pole_free_radius_deformation_isolates_exact_five_saddle_collar",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "precisions_decimal_digits": list(PRECISIONS),
            "source_attribution": {
                "lewis_2015": [8, 13, 24, 26, "B9", "B15a"],
                "lewis_brereton_2026": [62, 124, 125, 126, 127],
            },
        },
        "exact_identities": {
            "D_finite_interval": "D(N,R_s)-D(N,R_c)=integral_(u0(R_s))^(u0(R_c)) f_D(N,u) du",
            "saddle_threshold": "1/u_sad(N)=1-pi*alpha(N)/t, alpha(N)=2N+t/(pi*N)",
            "full_contour_compensation": "R_s*A(R_s)-i*t*S(R_s)/pi=R_c*A(R_c)-i*t*S(R_c)/pi",
            "S_definition": "S(R)=C(0,R)+D(0,R)+2*sum_(N>=1)(-1)^N(C(N,R)+D(N,R))",
        },
        "output_rows": output_rows,
        "aggregate": {
            "output_count": EXPECTED_OUTPUTS,
            "crossed_saddles_per_output": len(COLLAR),
            "total_crossed_saddle_certificates": EXPECTED_OUTPUTS * len(COLLAR),
            "crossed_saddle_roster": list(COLLAR),
            "lower_neighbor_inside_both_count": EXPECTED_OUTPUTS,
            "upper_neighbor_outside_both_count": EXPECTED_OUTPUTS,
            "common_left_pole_gap": [TERMINAL_ALPHA, NEXT_ODD_ALPHA],
            "common_far_pole_gap": [FAR_ODD_BELOW, FAR_ODD_ABOVE],
        },
        "decision": {
            "pole_free_radius_deformation_exists_at_all_outputs": True,
            "denominator_branch_points_enter_deformation_annulus": False,
            "alpha_pole_roster_changes_under_deformation": False,
            "exactly_five_D_saddles_cross_lower_limit": True,
            "diagnostic_midpoint_defines_a_contour_resplitting_candidate": True,
            "diagnostic_midpoint_is_a_published_cutoff_rule": False,
            "five_classical_terms_can_be_added_without_remainder_transport": False,
            "height_uniform_transported_remainder_proved": False,
            "finite_gate_has_rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "cell_partition_gate": {"path": relative(cells.RESULT), "sha256": file_hash(cells.RESULT)},
            "hardy_2026_paper": {"path": relative(partition.PAPER), "sha256": file_hash(partition.PAPER)},
            "lewis_2015_paper": {"path": relative(LEGACY_PAPER), "sha256": file_hash(LEGACY_PAPER)},
            "source_file": {"path": relative(partition.SOURCE), "sha256": file_hash(partition.SOURCE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Prioritize a direct source-aligned bound for equations (126)--(127). Treat transport to the diagnostic midpoint only as a separate alternative-cutoff theorem requiring all contour channels jointly.",
        "proof_boundary": "Rigorous finite pole and saddle geometry at fifteen saved heights only. It provides an exact contour-resplitting identity for an introduced diagnostic radius, not a published cutoff repair, a bound for the transported remainder, a uniform hybrid theorem, Lambda<=0, PF-infinity, RH, or a prize-level result.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated pole-free radius saddle collar: "
        f"outputs={EXPECTED_OUTPUTS}, crossed={list(COLLAR)}, "
        f"certificates={artifact['aggregate']['total_crossed_saddle_certificates']}, pole-roster-unchanged=15"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
