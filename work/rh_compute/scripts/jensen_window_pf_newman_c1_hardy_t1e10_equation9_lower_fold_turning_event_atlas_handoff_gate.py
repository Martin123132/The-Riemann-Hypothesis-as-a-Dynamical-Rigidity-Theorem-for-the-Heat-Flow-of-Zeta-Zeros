#!/usr/bin/env python3
"""Certify the fixed-selector turning-event atlas and exact mode handoff."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_displayed_height_window_gate as window_gate


WINDOW_GATE = window_gate.RESULT
PORTCULLIS_GATE = window_gate.PORTCULLIS_GATE
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 90
C = window_gate.cell_gate.ode_gate.C
T0 = window_gate.cell_gate.ode_gate.T
Q = (C - 1) // 4
EVENT_COUNT = 399
CURRENT_TRANSITION_COUNT = 84


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mode_for_event(index: int, q: int = Q, residue: int = 1) -> int:
    """Return the mode at the index-th descending event for C mod 4."""
    require(index >= 0, "event index must be nonnegative")
    if residue == 1:
        return q - index // 2 if index % 2 == 0 else q + (index + 1) // 2
    require(residue == 3, "odd selector residue must be 1 or 3")
    return q + 1 + index // 2 if index % 2 == 0 else q - (index - 1) // 2


def triangular(index: int) -> int:
    return index * (index + 1) // 2


def symbolic_identities() -> dict[str, str]:
    q, k = sp.symbols("q k", integer=True, nonnegative=True)
    pi = sp.pi
    endpoint = 4 * q + 1
    top = pi * endpoint**2 / 8

    n_even = 2 * k
    mode_even = q - k
    event_even = pi * mode_even * (endpoint - 2 * mode_even)
    target_even = top - pi * (n_even * (n_even + 1) / 2 + sp.Rational(1, 8))
    require(sp.simplify(event_even - target_even) == 0, "even event identity failed")

    n_odd = 2 * k + 1
    mode_odd = q + k + 1
    event_odd = pi * mode_odd * (endpoint - 2 * mode_odd)
    target_odd = top - pi * (n_odd * (n_odd + 1) / 2 + sp.Rational(1, 8))
    require(sp.simplify(event_odd - target_odd) == 0, "odd event identity failed")

    m, n, c = sp.symbols("m n C", integer=True)
    tau_m = m * (c - 2 * m)
    tau_n = n * (c - 2 * n)
    distinct_factor = sp.factor(tau_m - tau_n)
    require(
        sp.simplify(distinct_factor - (m - n) * (c - 2 * m - 2 * n)) == 0,
        "event distinctness factor failed",
    )
    return {
        "crossing_height": "tau_m=pi*m*(C-2m), obtained from alpha_m=2m+t/(pi*m)=C",
        "square_deficit": "pi*C^2/8-tau_m=pi*(C-4m)^2/8",
        "ordered_height": "tau_n=pi*C^2/8-pi*(n(n+1)/2+1/8)",
        "mode_order_for_C_4q_plus_1": "m_(2k)=q-k, m_(2k+1)=q+k+1",
        "successive_spacing": "tau_n-tau_(n+1)=pi*(n+1)",
        "distinctness_factor": str(distinct_factor),
        "odd_selector_distinctness": "For m!=n, equality would require 2(m+n)=C, impossible for odd C.",
        "exact_handoff": "sum_(m in O)J_m+sum_(m in T)J_m=sum_(m in O union {r})J_m+sum_(m in T minus {r})J_m",
    }


def event_height(index: int, endpoint: int = C) -> arb:
    return arb.pi() * arb(endpoint) ** 2 / 8 - arb.pi() * (arb(triangular(index)) + arb(1) / 8)


def contiguous_bounds(count: int, q: int = Q, residue: int = 1) -> tuple[int, int]:
    require(count > 0, "count must be positive")
    modes = [mode_for_event(index, q=q, residue=residue) for index in range(count)]
    require(set(modes) == set(range(min(modes), max(modes) + 1)), f"first {count} modes are not contiguous")
    return min(modes), max(modes)


def validate_all_handoffs() -> dict[str, Any]:
    sequence = [mode_for_event(index) for index in range(EVENT_COUNT)]
    require(len(sequence) == len(set(sequence)) == EVENT_COUNT, "event mode sequence repeats")
    checks = 0
    transition = set(sequence)
    ordinary: set[int] = set()
    master = set(sequence)
    for count in range(EVENT_COUNT, 0, -1):
        leaving = sequence[count - 1]
        require(leaving in transition and leaving not in ordinary, "handoff source ownership failed")
        before_union = ordinary | transition
        transition.remove(leaving)
        ordinary.add(leaving)
        require((ordinary | transition) == before_union == master, "handoff lost or duplicated a mode")
        require(ordinary.isdisjoint(transition), "handoff ownership overlaps")
        checks += 1
    require(transition == set() and ordinary == master, "handoff atlas did not exhaust")
    return {
        "checked_one_mode_handoffs": checks,
        "all_handoffs_preserve_disjoint_union": True,
        "identity_scope": "Exact finite mode sums J_m only; not equality of truncated Airy and ordinary-Morse approximations.",
    }


def certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    require(C == 4 * Q + 1, "saved endpoint is not 1 modulo 4")
    pi = arb.pi()
    top = pi * C**2 / 8
    lower = pi * (C - 2) ** 2 / 8
    selector_width = top - lower
    require(abs(selector_width - 2 * pi * Q).contains(0), "selector width identity failed")

    require(triangular(398) <= 2 * Q - 1, "event 398 is outside selector cell")
    require(triangular(399) >= 2 * Q, "event 399 unexpectedly enters selector cell")
    events = [event_height(index) for index in range(EVENT_COUNT)]
    require(all(lower < value < top for value in events), "saved event lies outside selector cell")
    require(event_height(EVENT_COUNT) < lower, "next event is not below selector cell")
    spacings = [events[index] - events[index + 1] for index in range(EVENT_COUNT - 1)]
    require(all(spacing > 0 for spacing in spacings), "event ordering failed")
    require(all((spacing - pi * (index + 1)).contains(0) for index, spacing in enumerate(spacings)), "event spacing law failed")

    current = arb(T0)
    tau_83 = event_height(83)
    tau_84 = event_height(84)
    require(tau_83 > current > tau_84, "saved height has the wrong event count")
    current_modes = [mode_for_event(index) for index in range(CURRENT_TRANSITION_COUNT)]
    require(set(current_modes) == set(range(39853, 39937)), "saved 84-mode roster mismatch")
    current_lo, current_hi = contiguous_bounds(CURRENT_TRANSITION_COUNT)
    window_radius = arb(window_gate.HEIGHT_RADIUS_TEXT)
    require(tau_83 > current + window_radius and current - window_radius > tau_84, "displayed window crosses an event")

    full_lo, full_hi = contiguous_bounds(EVENT_COUNT)
    next_endpoint = C + 2
    next_q = (next_endpoint - 3) // 4
    next_lo, next_hi = contiguous_bounds(EVENT_COUNT, q=next_q, residue=3)
    next_lower = pi * (next_endpoint - 2) ** 2 / 8
    require((next_lower - top).contains(0), "adjacent selector boundary mismatch")
    require((full_lo, full_hi) == (39695, 40093), "saved lower-edge roster mismatch")
    require((next_lo, next_hi) == (39696, 40094), "adjacent lower-edge roster mismatch")

    return {
        "fixed_selector_C": C,
        "selector_residue_mod_4": C % 4,
        "q": Q,
        "selector_cell_lower_ball": lower.str(PRECISION, more=True),
        "selector_cell_upper_ball": top.str(PRECISION, more=True),
        "selector_cell_width_ball": selector_width.str(PRECISION, more=True),
        "event_count_inside_selector_cell": EVENT_COUNT,
        "first_event_index": 0,
        "last_event_index": EVENT_COUNT - 1,
        "highest_event_mode": mode_for_event(0),
        "highest_event_distance_below_selector_top_ball": (top - events[0]).str(PRECISION, more=True),
        "minimum_event_spacing_ball": min(spacings, key=float).str(PRECISION, more=True),
        "lower_edge_transition_roster": [full_lo, full_hi],
        "lower_edge_transition_count": EVENT_COUNT,
        "saved_height": T0,
        "saved_transition_count": CURRENT_TRANSITION_COUNT,
        "saved_transition_roster": [current_lo, current_hi],
        "next_upward_event_index": 83,
        "next_upward_event_mode": mode_for_event(83),
        "next_upward_event_height_ball": tau_83.str(PRECISION, more=True),
        "next_upward_event_distance_ball": (tau_83 - current).str(PRECISION, more=True),
        "next_downward_event_index": 84,
        "next_downward_event_mode": mode_for_event(84),
        "next_downward_event_height_ball": tau_84.str(PRECISION, more=True),
        "next_downward_event_distance_ball": (current - tau_84).str(PRECISION, more=True),
        "displayed_window_upper_event_margin_ball": (tau_83 - current - window_radius).str(PRECISION, more=True),
        "displayed_window_lower_event_margin_ball": (current - window_radius - tau_84).str(PRECISION, more=True),
        "adjacent_selector_C": next_endpoint,
        "adjacent_selector_common_boundary_ball": top.str(PRECISION, more=True),
        "adjacent_selector_lower_edge_transition_roster": [next_lo, next_hi],
        "adjacent_selector_lower_edge_transition_count": EVENT_COUNT,
        "selector_jump_is_one_mode_handoff": False,
        "selector_jump_geometry": (
            "At C->C+2 the old fixed-C transition roster is empty near its upper edge, "
            "whereas the new selector starts with 399 modes at the removed lower alpha strip [C,C+2]."
        ),
        "handoff_validation": validate_all_handoffs(),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_event_atlas"]
    return f"""# Fixed-selector turning-event atlas and exact mode handoff

Date: 2026-08-12

Status: exact fixed-selector event ordering and finite-sum handoff validated;
not a proof of the Airy-to-Morse approximation join or selector jump

For a fixed odd lower endpoint `C`, a mode changes turning ownership when
its joint saddle reaches `alpha=C`:

```text
alpha_m=2m+t/(pi*m)=C,
tau_m=pi*m*(C-2m).                                    (TE1)
```

The distance from the selector top is the exact square

```text
pi*C^2/8-tau_m=pi*(C-4m)^2/8.                         (TE2)
```

Here `C=159577=4q+1`, `q=39894`.  Ordering the odd squares gives

```text
m_(2k)=q-k,       m_(2k+1)=q+k+1,
tau_n=pi*C^2/8-pi[n(n+1)/2+1/8],
tau_n-tau_(n+1)=pi(n+1).                              (TE3)
```

Thus all event heights are distinct and their minimum spacing is

```text
{c['minimum_event_spacing_ball']}.                     (TE4)
```

There are exactly {c['event_count_inside_selector_cell']} events in the fixed
selector cell.  Immediately above its lower edge their modes form the
contiguous roster `{c['lower_edge_transition_roster'][0]}..{c['lower_edge_transition_roster'][1]}`.
At `t=10^10`,

```text
tau_83 > t > tau_84,
transition roster =39853..39936,                       (TE5)

tau_83-t={c['next_upward_event_distance_ball']},
t-tau_84={c['next_downward_event_distance_ball']}.     (TE6)
```

The next upward event transfers mode `{c['next_upward_event_mode']}` from the
turning block to the ordinary block.  The next downward event transfers mode
`{c['next_downward_event_mode']}` in the reverse direction.  For every one of
the {c['event_count_inside_selector_cell']} events, direct set checks certify
the exact finite-sum identity

```text
sum_(m in O) J_m + sum_(m in T) J_m
 =sum_(m in O union {{r}}) J_m + sum_(m in T minus {{r}}) J_m. (TE7)
```

No mode is lost or counted twice.  This identity uses the same exact finite
Poisson mode integral `J_m` on both sides; it does not assert that separately
truncated ordinary-Morse and Airy approximations agree.

The odd-selector jump is a different event.  At `C -> C+2`, the old turning
roster is empty near its upper edge, while the new selector begins with the
399-mode roster `{c['adjacent_selector_lower_edge_transition_roster'][0]}..{c['adjacent_selector_lower_edge_transition_roster'][1]}`
at the removed lower-endpoint strip `[C,C+2]`.  It cannot be replaced by the
one-mode identity (TE7).

Proof boundary: exact event geometry and ownership conservation for exact
finite mode sums only.  No matched Airy/Morse remainder, 399-mode selector
strip reassembly, lower-interior join, complete `T_upper`, `Lambda<=0`, RH,
or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = window_gate.cell_gate.ode_gate.set_low_priority()
    for dependency in (WINDOW_GATE, PORTCULLIS_GATE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    ctx.dps = PRECISION
    certified = certificate()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate",
        "status": "exact_fixed_selector_turning_event_atlas_and_one_mode_handoff_complete",
        "passed": True,
        "symbolic_identities": symbolic_identities(),
        "certified_event_atlas": certified,
        "decision": {
            "all_fixed_selector_event_heights_exactly_ordered": True,
            "all_399_one_mode_handoffs_preserve_exact_finite_sum": True,
            "displayed_window_crosses_no_turning_event": True,
            "selector_jump_reducible_to_one_mode_handoff": False,
            "airy_morse_approximation_join_proved": False,
            "selector_strip_reassembly_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "displayed_window_gate": {"path": relative(WINDOW_GATE), "sha256": file_hash(WINDOW_GATE)},
            "portcullis_gate": {"path": relative(PORTCULLIS_GATE), "sha256": file_hash(PORTCULLIS_GATE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_obligation": (
            "Keep each crossing mode exact in a disjoint event buffer and derive explicit ordinary-Morse "
            "and Airy remainders outside that buffer; separately derive the 399-mode removed-strip identity at C->C+2."
        ),
        "proof_boundary": (
            "Exact fixed-selector event geometry and finite-sum ownership only. No Airy/Morse approximation join, "
            "selector-strip reassembly, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built turning-event atlas: 399 exact handoffs, current roster=84, selector jump is not one-mode")


if __name__ == "__main__":
    main()
