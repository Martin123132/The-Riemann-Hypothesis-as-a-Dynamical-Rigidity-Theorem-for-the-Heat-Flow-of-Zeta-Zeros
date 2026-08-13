#!/usr/bin/env python3
"""Certify the moving contiguous mode roster induced by the y=64 charts."""

from __future__ import annotations

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
import sympy as sp


PARTITION_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_exact_normal_saddle_partition_gate.json"
PRIMITIVE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_three_chart_primitive_envelope_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_moving_roster_dirichlet_partition_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_moving_roster_dirichlet_partition_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 90
T = 10_000_000_000
C = 159577
MODE_LO = 39853
MODE_HI = 39936
Y0 = 64
RADIUS = 9
FRESNEL_WIDTHS = 4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def atanh_ball(value: arb) -> arb:
    return ((1 + value) / (1 - value)).log() / 2


def symbolic_roster() -> dict[str, str]:
    h, m, c_fourth, G, delta = sp.symbols("h m C4 G delta", positive=True, real=True)
    g = h * (m - c_fourth) + G
    require(sp.simplify(g.subs(m, m + 1) - g - h) == 0, "mode spacing identity failed")
    lower = sp.simplify(c_fourth + (-delta - G) / h)
    upper = sp.simplify(c_fourth + (delta - G) / h)
    require(sp.simplify(upper - lower - 2 * delta / h) == 0, "transition width identity failed")

    ratio, length = sp.symbols("r L", nonzero=True)
    geometric = (1 - ratio**length) / (1 - ratio)
    return {
        "affine_mode_function": "g_m(z)=h(m-C/4)+G(z), G(z)=beta*tanh(z/beta)-64x(z)/beta",
        "mode_spacing": "g_(m+1)(z)-g_m(z)=h>0",
        "real_cutoffs": "mu_-(z)=C/4+[-delta-G(z)]/h, mu_+(z)=C/4+[delta-G(z)]/h",
        "blocks": "m<=mu_-: nonstationary; mu_-<m<mu_+: Fresnel; m>=mu_+: Morse",
        "Fresnel_roster_width": "mu_+-mu_-=2delta/h",
        "geometric_kernel": str(geometric),
        "contiguous_Dirichlet_kernel": "D_[a,b](y)=exp[-i(d_a+d_b)y/2] sin((b-a+1)hy/2)/sin(hy/2)",
        "removable_values": "At sin(hy/2)=0 use the finite geometric sum by continuity.",
        "event_assignment": "g=-delta belongs to nonstationary; g=+delta belongs to Morse; open z cells carry fixed rosters.",
    }


def block_record(modes: list[int]) -> dict[str, Any]:
    if not modes:
        return {"count": 0, "first": None, "last": None}
    require(modes == list(range(modes[0], modes[-1] + 1)), "mode block is not contiguous")
    return {"count": len(modes), "first": modes[0], "last": modes[-1]}


def certified_moving_roster() -> dict[str, Any]:
    ctx.dps = PRECISION
    t = arb(T)
    endpoint_c = arb(C)
    pi = arb.pi()
    eta = pi * endpoint_c**2 / (8 * t)
    beta = (t * eta) ** (arb(1) / 3)
    h = 4 * beta / endpoint_c

    def tau(z: arb) -> arb:
        return (z / beta).tanh()

    def x(z: arb) -> arb:
        return (1 - tau(z)) / 2

    x_max = x(arb(-RADIUS))
    delta = arb(FRESNEL_WIDTHS) * (x_max / beta).sqrt()
    width = 2 * delta / h
    denominator = beta + arb(Y0) / (2 * beta)
    detunings = {
        mode: h * (arb(mode) - endpoint_c / 4)
        for mode in range(MODE_LO, MODE_HI + 1)
    }

    def root(level: arb, mode: int) -> arb:
        argument = (level - detunings[mode] + arb(Y0) / (2 * beta)) / denominator
        return beta * atanh_ball(argument)

    events: list[dict[str, Any]] = []
    event_balls: list[tuple[arb, int, str]] = []
    for mode in range(MODE_LO, MODE_HI + 1):
        left = root(-delta, mode)
        right = root(delta, mode)
        event_balls.append((left, mode, "enter_Fresnel"))
        event_balls.append((right, mode, "enter_Morse"))
    event_balls.sort(key=lambda row: float(row[0].mid()))

    gaps = []
    for index, (z, mode, event_type) in enumerate(event_balls):
        events.append(
            {
                "event_index": index,
                "z_ball": z.str(PRECISION, more=True),
                "mode": mode,
                "event": event_type,
            }
        )
        if index:
            gap = z - event_balls[index - 1][0]
            require(gap > 0, "event ordering is not rigorous")
            gaps.append(gap)

    boundaries = [arb(-RADIUS)] + [row[0].mid() for row in event_balls] + [arb(RADIUS)]
    cells = []
    histogram: dict[int, int] = {}
    for index, (left, right) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        sample = (left + right) / 2
        common = beta * tau(sample) - arb(Y0) * x(sample) / beta
        nonstationary = []
        fresnel = []
        morse = []
        for mode in range(MODE_LO, MODE_HI + 1):
            g = detunings[mode] + common
            if g < -delta:
                nonstationary.append(mode)
            elif g > delta:
                morse.append(mode)
            else:
                require(-delta < g < delta, "sample intersects an event")
                fresnel.append(mode)
        require(len(nonstationary) + len(fresnel) + len(morse) == MODE_HI - MODE_LO + 1, "cell roster lost a mode")
        n_block = block_record(nonstationary)
        f_block = block_record(fresnel)
        m_block = block_record(morse)
        occupied = nonstationary + fresnel + morse
        require(occupied == list(range(MODE_LO, MODE_HI + 1)), "cell blocks overlap or leave a gap")
        histogram[len(fresnel)] = histogram.get(len(fresnel), 0) + 1
        cells.append(
            {
                "cell_index": index,
                "z_left": left.str(PRECISION, more=True),
                "z_right": right.str(PRECISION, more=True),
                "nonstationary": n_block,
                "Fresnel": f_block,
                "Morse": m_block,
            }
        )

    distance_to_integer = min(abs(width - arb(2)), abs(arb(3) - width))
    require(arb(2) < width < arb(3), "Fresnel roster width is not between two and three")
    require(len(event_balls) == 168 and len(cells) == 169, "event/cell count drift")
    require(min(gaps) > arb("0.014"), "event separation below 0.014")
    require(max(histogram) == 3, "more than three simultaneous Fresnel modes")
    require(sum(histogram.values()) == 169, "cell histogram does not close")
    require(histogram.get(3) == 82 and histogram.get(2) == 83, "interior cell histogram drift")

    return {
        "height": T,
        "mode_range": [MODE_LO, MODE_HI],
        "mode_count": MODE_HI - MODE_LO + 1,
        "z_range": [-RADIUS, RADIUS],
        "h_ball": h.str(PRECISION, more=True),
        "delta_ball": delta.str(PRECISION, more=True),
        "real_Fresnel_roster_width_ball": width.str(PRECISION, more=True),
        "distance_of_width_to_nearest_integer_ball": distance_to_integer.str(PRECISION, more=True),
        "event_count": len(events),
        "open_cell_count": len(cells),
        "minimum_event_gap_ball": min(gaps).str(PRECISION, more=True),
        "maximum_simultaneous_Fresnel_modes": max(histogram),
        "Fresnel_count_cell_histogram": {str(key): value for key, value in sorted(histogram.items())},
        "events": events,
        "cells": cells,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_moving_roster"]
    return f"""# Moving y=64 roster and contiguous Dirichlet blocks

Date: 2026-08-11

Status: exact moving-roster partition validated; not a proof of the grouped
chart integral bound

The mode dependence of the boundary function is affine:

```text
g_m(z)=h(m-C/4)+G(z),
G(z)=beta*tanh(z/beta)-64x(z)/beta,                     (MR1)

g_(m+1)(z)-g_m(z)=h>0.                                 (MR2)
```

Consequently, at every fixed `z`, the 84 modes split into three contiguous
integer blocks:

```text
m<=mu_-(z):             nonstationary,
mu_-(z)<m<mu_+(z):      endpoint/Fresnel,
m>=mu_+(z):             ordinary Morse,                (MR3)

mu_+-mu_-=2delta/h
 ={c['real_Fresnel_roster_width_ball']}.                (MR4)
```

Since `2<2delta/h<3`, at most three integer modes can lie in the Fresnel
block at one `z`.  Thus the per-mode transition envelope from Section 11.320
is never paid 84 times simultaneously.

The 168 entry and exit events are all distinct.  Their minimum separation is

```text
{c['minimum_event_gap_ball']}>0.014,                    (MR5)
```

so they produce 169 open `z` cells with fixed rosters.  The cell histogram is

```text
Fresnel count 0: {c['Fresnel_count_cell_histogram'].get('0', 0)} cells,
Fresnel count 1: {c['Fresnel_count_cell_histogram'].get('1', 0)} cells,
Fresnel count 2: {c['Fresnel_count_cell_histogram'].get('2', 0)} cells,
Fresnel count 3: {c['Fresnel_count_cell_histogram'].get('3', 0)} cells. (MR6)
```

At `g=-delta` equality belongs to the nonstationary block; at `g=+delta`
it belongs to the Morse block.  The isolated event points have zero `z`
measure, and the half-open assignment makes the full partition exact.

Each nonempty block retains an exact finite Dirichlet kernel.  For consecutive
`a<=m<=b`,

```text
D_[a,b](y)
 =exp[-i(d_a+d_b)y/2]
  sin((b-a+1)hy/2)/sin(hy/2),                          (MR7)
```

with removable denominator zeros interpreted by the original finite sum.
This supplies a cancellation-preserving route for the nonstationary and
Morse blocks rather than summing their per-mode absolute envelopes.

Proof boundary: exact saved-height moving-roster algebra, event atlas, and
contiguous-kernel decomposition only.  No grouped chart integral estimate,
completed z integral, lower-interior join, `T_upper` theorem, `Lambda<=0`,
RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(PARTITION_GATE.is_file() and PRIMITIVE_GATE.is_file() and CHECKER.is_file(), "missing dependency or checker")
    symbolic = symbolic_roster()
    certified = certified_moving_roster()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_moving_roster_dirichlet_partition_gate",
        "status": "exact_y64_moving_roster_contiguous_Dirichlet_partition_complete",
        "passed": True,
        "symbolic_roster": symbolic,
        "certified_moving_roster": certified,
        "decision": {
            "chart_memberships_form_contiguous_mode_blocks": True,
            "Fresnel_block_has_at_most_three_modes": True,
            "all_168_events_distinct": True,
            "open_cell_rosters_exhaustive_and_disjoint": True,
            "contiguous_blocks_retain_exact_Dirichlet_kernels": True,
            "grouped_chart_integral_bound_proved": False,
            "lower_interior_join_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "On each fixed-roster z cell, apply summation by parts or the exact finite Dirichlet kernel to the nonstationary and Morse blocks, retain at most three transition modes explicitly, and rejoin all 169 cell boundaries without absolute-value duplication.",
        "proof_boundary": "Exact saved-height moving-roster algebra, event atlas, and contiguous-kernel decomposition only. No grouped chart integral estimate, completed z integral, lower-interior join, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "partition_gate": {"path": relative(PARTITION_GATE), "sha256": file_hash(PARTITION_GATE)},
            "primitive_gate": {"path": relative(PRIMITIVE_GATE), "sha256": file_hash(PRIMITIVE_GATE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"workers": 1, "process_priority": priority, "elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built moving y=64 roster: 168 events, 169 cells, at most 3 Fresnel modes")


if __name__ == "__main__":
    main()
