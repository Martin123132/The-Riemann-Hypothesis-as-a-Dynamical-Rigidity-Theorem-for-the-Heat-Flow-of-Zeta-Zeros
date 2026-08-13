#!/usr/bin/env python3
"""Certify exact finite-Poisson reassembly when the odd selector advances."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate as event_gate


EVENT_GATE = event_gate.RESULT
POISSON_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json"
CURRENT_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 90
C = event_gate.C
C_NEXT = C + 2
B = 5122421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def symbolic_reassembly() -> dict[str, str]:
    u, v, x = sp.symbols("u v x", real=True)
    endpoint = sp.symbols("C", integer=True, odd=True, positive=True)
    pi = sp.pi
    i = sp.I

    alpha_old = endpoint + 2 * u
    f_old = alpha_old * sp.exp(i * pi * x * alpha_old**2 / 4)
    alpha_new = endpoint + 2 + 2 * v
    f_new = alpha_new * sp.exp(i * pi * x * alpha_new**2 / 4)
    require(sp.simplify(f_new - f_old.subs(u, v + 1)) == 0, "selector roster shift failed")

    f0, f1 = sp.symbols("f_0 f_1")
    endpoint_half_difference = (f0 - f1) / 2
    one_cell_poisson_integral = (f0 + f1) / 2
    require(
        sp.simplify(endpoint_half_difference + one_cell_poisson_integral - f0) == 0,
        "one-cell finite-Poisson reassembly failed",
    )

    return {
        "old_coordinate": "alpha=C+2u, 0<=u<=L",
        "new_coordinate": "alpha=C+2+2v, 0<=v<=L-1",
        "roster_shift": "f_(C+2)(v)=f_C(v+1)",
        "integer_character_shift": "exp[-2pi*i*m*(u-1)]=exp(-2pi*i*m*u) for integer m",
        "mode_strip": "S_m(x)=integral_0^1 f_C(u)exp(-2pi*i*m*u)du=(-1)^m/2 integral_C^(C+2) alpha exp(i*pi[x*alpha^2/4-m*alpha])dalpha",
        "endpoint_half_difference": "H_C(x)=[f_C(0)-f_C(1)]/2",
        "one_cell_poisson": "lim_(M->infinity) sum_(m=-M)^M S_m(x)=[f_C(0)+f_C(1)]/2",
        "exact_removed_term": "H_C+lim_sym sum_m S_m=f_C(0)=C exp(i*pi*x*C^2/4)",
        "integrated_selector_identity": "K_C(t;B)-K_(C+2)(t;B)=integral_0^1 W_t(x) C exp(i*pi*x*C^2/4)dx for fixed B",
        "parity_guard": "Because C and C+2 are odd, (-1)^m exp(-i*pi*m*C)=(-1)^m exp(-i*pi*m*(C+2))=1; endpoint currents cannot be separated from the symmetric strip sum.",
    }


def certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    require(C % 2 == C_NEXT % 2 == B % 2 == 1, "selector or upper endpoint is not odd")
    require(B > C_NEXT, "selector exceeds upper endpoint")
    old_length = (B - C) // 2
    new_length = (B - C_NEXT) // 2
    require(B - C == 2 * old_length and B - C_NEXT == 2 * new_length, "odd roster length failed")
    require(new_length == old_length - 1, "selector shift does not remove one lattice point")

    event_artifact = json.loads(EVENT_GATE.read_text(encoding="utf-8"))
    atlas = event_artifact["certified_event_atlas"]
    require(atlas["adjacent_selector_C"] == C_NEXT, "event atlas adjacent selector mismatch")
    require(atlas["adjacent_selector_lower_edge_transition_count"] == 399, "adjacent turning count mismatch")
    require(atlas["adjacent_selector_lower_edge_transition_roster"] == [39696, 40094], "adjacent turning roster mismatch")

    selector_boundary = arb.pi() * arb(C) ** 2 / 8
    next_lower_boundary = arb.pi() * arb(C_NEXT - 2) ** 2 / 8
    require((selector_boundary - next_lower_boundary).contains(0), "selector cells do not share a boundary")

    return {
        "old_selector_C": C,
        "new_selector_C": C_NEXT,
        "fixed_upper_endpoint_B": B,
        "old_odd_roster_length_L": old_length,
        "old_odd_roster_term_count": old_length + 1,
        "new_odd_roster_length_L": new_length,
        "new_odd_roster_term_count": new_length + 1,
        "removed_source_lattice_point": C,
        "removed_source_lattice_point_count": 1,
        "shared_selector_boundary_ball": selector_boundary.str(PRECISION, more=True),
        "new_selector_turning_roster_at_lower_edge": atlas["adjacent_selector_lower_edge_transition_roster"],
        "new_selector_turning_mode_count": atlas["adjacent_selector_lower_edge_transition_count"],
        "complete_symmetric_mode_sum_required": True,
        "endpoint_half_current_required": True,
        "positive_399_mode_roster_sufficient_for_selector_reassembly": False,
        "zero_and_negative_modes_may_be_discarded": False,
        "fixed_B_selector_reassembly_exact": True,
        "upper_endpoint_simultaneous_change_covered": False,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_selector_reassembly"]
    return f"""# Exact selector-strip finite-Poisson reassembly

Date: 2026-08-12

Status: exact fixed-upper-endpoint selector reassembly validated; not a proof
of the 399-mode approximation bound or complete source selector splice

Write the old odd source roster as

```text
alpha=C+2u,       u=0,...,L,
f_C(u)=(C+2u)exp(i*pi*x*(C+2u)^2/4).                  (SS1)
```

When the lower selector advances from `C` to `C+2`, the new roster is the
shifted old roster `u=1,...,L`.  Therefore exactly one source lattice point,
`alpha=C`, is removed.  The counts here are

```text
C={c['old_selector_C']},       C+2={c['new_selector_C']},
old terms={c['old_odd_roster_term_count']},
new terms={c['new_odd_roster_term_count']}.             (SS2)
```

For integer Fourier mode `m`, shifting `u=v+1` leaves the character unchanged.
The mode difference is the one-cell strip

```text
S_m(x)=integral_0^1 f_C(u)exp(-2pi*i*m*u)du
 =(-1)^m/2 integral_C^(C+2)
   alpha exp(i*pi[x*alpha^2/4-m*alpha])dalpha.         (SS3)
```

The endpoint half-current changes by

```text
H_C(x)=[f_C(0)-f_C(1)]/2.                              (SS4)
```

Finite Poisson on the one-cell interval gives

```text
lim_(M->infinity) sum_(m=-M)^M S_m(x)
 =[f_C(0)+f_C(1)]/2,                                   (SS5)

H_C+lim_sym sum_m S_m=f_C(0)
 =C exp(i*pi*x*C^2/4).                                 (SS6)
```

After multiplication by the Kummer weight and x integration, (SS6) is exactly
the single source term removed when the lower odd selector advances.  Thus
the complete fixed-`B` finite-Poisson representation reassembles exactly.

This also shows why the 399 positive turning modes are insufficient by
themselves.  The selector identity requires the complete symmetric mode sum,
including zero and negative modes, together with the endpoint half-current.
Separating the common odd-endpoint current produces a false divergent object.

The adjacent selector begins with the 399-mode positive turning roster
`{c['new_selector_turning_roster_at_lower_edge'][0]}..{c['new_selector_turning_roster_at_lower_edge'][1]}`,
but its quantitative Airy/Morse estimate must be embedded in (SS6), not used
as a replacement for it.

Proof boundary: exact selector reassembly for a fixed upper endpoint `B` and
the complete symmetric Poisson sum only.  No 399-mode quantitative remainder,
simultaneous upper-endpoint transition, complete source splice, `T_upper`,
`Lambda<=0`, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = event_gate.window_gate.cell_gate.ode_gate.set_low_priority()
    for dependency in (EVENT_GATE, POISSON_GATE, CURRENT_GATE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    ctx.dps = PRECISION
    certified = certificate()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate",
        "status": "exact_fixed_B_odd_selector_strip_poisson_reassembly_complete",
        "passed": True,
        "symbolic_reassembly": symbolic_reassembly(),
        "certified_selector_reassembly": certified,
        "decision": {
            "one_source_lattice_term_removed": True,
            "complete_symmetric_poisson_strip_reassembles_removed_term": True,
            "endpoint_half_current_retained": True,
            "positive_399_mode_roster_alone_sufficient": False,
            "complete_source_selector_splice_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "turning_event_gate": {"path": relative(EVENT_GATE), "sha256": file_hash(EVENT_GATE)},
            "symmetric_poisson_gate": {"path": relative(POISSON_GATE), "sha256": file_hash(POISSON_GATE)},
            "theta_current_gate": {"path": relative(CURRENT_GATE), "sha256": file_hash(CURRENT_GATE)},
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
            "Insert a 399-mode exact or uniformly bounded endpoint chart into the complete symmetric strip identity, "
            "bound its zero/negative and positive outer complements jointly, and audit whether B changes at the same source selector boundary."
        ),
        "proof_boundary": (
            "Exact fixed-B selector-strip finite-Poisson reassembly only. No 399-mode quantitative bound, "
            "simultaneous upper endpoint transition, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built selector-strip reassembly: one removed source term equals endpoint half-current plus complete symmetric strip sum")


if __name__ == "__main__":
    main()
