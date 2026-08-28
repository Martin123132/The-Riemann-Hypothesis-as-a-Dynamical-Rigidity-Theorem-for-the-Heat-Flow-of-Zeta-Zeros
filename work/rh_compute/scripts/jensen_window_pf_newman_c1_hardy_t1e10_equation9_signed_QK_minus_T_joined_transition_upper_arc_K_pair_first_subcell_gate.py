#!/usr/bin/env python3
"""Certify the transition-plus-upper-arc K_T ownership pair on I_1."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from flint import arb

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_transition_height_derivative_scout as transition_scout
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_upper_arc_K_first_subcell_scout as arc_scout


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_transition_upper_arc_K_pair_first_subcell_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
TRANSITION_SCOUT = SCRIPT_DIR / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_transition_height_derivative_scout.py"
)
ARC_SCOUT = SCRIPT_DIR / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_upper_arc_K_first_subcell_scout.py"
)
DEPENDENCIES = {
    "transition_K_gate": ROOT / "work" / "rh_compute" / "results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_transition_K_first_subcell_gate.json"
    ),
    "first_height_subcell": ROOT / "work" / "rh_compute" / "results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_first_nonzero_height_subcell_gate.json"
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    transition = c["transition_Hardy_derivative_contribution"]
    arc = c["upper_arc_Hardy_derivative_contribution"]
    pair = c["paired_Hardy_derivative_contribution"]
    return f"""# Transition--upper-arc K pair on the first subcell

Date: 2026-08-28

Status: rigorous two-component Hardy-derivative interval; complete `K_T`
remains open

On `I_1=[10^10-10^-4,10^10+10^-4]`, let `K_tr` be the cancellation-stable
transition operator packet and `K_U` the cancellation-stable compact upper-arc
operator packet.  Their common Hardy projection gives

```text
Hardy_t[K_tr]/H={transition['ball']},               (TP1)

Hardy_t[K_U]/H ={arc['ball']}.                      (TP2)
```

The signs are opposite and the centres nearly cancel.  Joining before taking
an absolute value yields

```text
Hardy_t[K_tr+K_U]/H={pair['ball']},                 (TP3)

sup_I1 |Hardy_t[K_tr+K_U]/H| < 1.3e-6.             (TP4)
```

For `K_tr`, the operator is pushed through the alternating boundary-jet source
and through the 42-mode Gamma sum per mode.  For `K_U`, the rotated-arc weight
is

```text
delta+i(theta'-log(Y))-H'/H,
```

and the final projection transports the combined phase `theta+phi_U` with
derivative `theta'-log(Y)`.  Thus neither component is obtained by separately
widening its large derivative and phase terms.

The independent checker raises both calculations to 448 bits.  It changes the
transition jet order, derivative slabs, integral panels, and logarithm-series
order, while changing the arc Taylor degree, endpoint cutoff, initial panels,
and all main-panel breakpoints.  Both direct operator packets overlap their
dependency-heavy assembled forms.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
integer Fourier spacing, Gamma normalization, Riemann--Siegel phase, or
`p=t/(2*pi)`.  No circle fit, polygon fit, or visual pattern supplies `pi`.

Proof boundary: (TP1)--(TP4) cover only the transition and compact upper arc.
The positive-real tail derivative, lower finite cell, ordinary complement,
and tiny target correction remain outside the pair.  Therefore this does not
prove a complete nonzero-radius `K_T` bound, a wider `Q_K-T` sign interval,
the full event-cell sign theorem, an event-wall handoff, an all-height theorem,
`Lambda<=0`, RH, or a prize-level conclusion.
"""


def main() -> int:
    dependencies = {
        name: {"path": relative(path), "sha256": file_hash(path), "artifact": load_json(path)}
        for name, path in DEPENDENCIES.items()
    }
    require(
        dependencies["transition_K_gate"]["artifact"]["decision"]
        ["transition_K_nonzero_radius_interval_proved"]
        is True,
        "transition K dependency drift",
    )
    require(
        dependencies["first_height_subcell"]["artifact"]["decision"]
        ["first_nonzero_radius_Q_K_minus_T_sign_theorem_proved"]
        is True,
        "first height subcell drift",
    )

    transition = transition_scout.build("0.0001", 384, "production")
    arc = arc_scout.build("0.0001", 384, "production")
    transition_projection = arb(transition["transition_Q_derivative_contribution"]["ball"])
    arc_projection = arb(arc["upper_arc_Q_derivative_contribution"]["ball"])
    pair = transition_projection + arc_projection
    require(transition_projection.lower() > 0, "transition contribution lost positivity")
    require(arc_projection.upper() < 0, "upper-arc contribution lost negativity")
    require(abs(pair).upper() < arb("1.3e-6"), "paired derivative contribution misses 1.3e-6 cap")

    certificate = {
        "height_ball": transition["height_ball"],
        "radius": "0.0001",
        "precision_bits": 384,
        "transition_configuration": transition["configuration"],
        "upper_arc_configuration": arc["configuration"],
        "transition_Hardy_derivative_contribution": transition[
            "transition_Q_derivative_contribution"
        ],
        "upper_arc_Hardy_derivative_contribution": arc[
            "upper_arc_Q_derivative_contribution"
        ],
        "paired_Hardy_derivative_contribution": {
            "ball": pair.str(55, more=True),
            "absolute_upper": abs(pair).upper().str(40, more=True),
        },
        "transition_K_direct": transition["transition_K_direct"],
        "upper_arc_K_direct": arc["actual_K_direct"],
        "transition_direct_assembled_overlap": transition[
            "direct_and_assembled_transition_K_overlap"
        ],
        "upper_arc_direct_assembled_overlap": arc["direct_and_assembled_K_overlap"],
        "upper_arc_error_budgets": arc["error_budgets"],
        "upper_arc_phase_transport": arc["phase_transport"],
        "altered_transition_fixture": transition["altered_five_label_fixture"],
    }
    artifact = {
        "kind": "rh_c1_hardy_joined_transition_upper_arc_K_pair_first_subcell_gate",
        "date": "2026-08-28",
        "status": "joined_transition_upper_arc_K_pair_first_height_subcell_interval_certified",
        "passed": True,
        "certificate": certificate,
        "decision": {
            "transition_upper_arc_K_pair_interval_proved": True,
            "paired_Hardy_derivative_absolute_upper_below_1_3e_minus_6": True,
            "complete_K_T_interval_proved": False,
            "wider_Q_K_minus_T_sign_interval_proved": False,
            "maximal_event_cell_sign_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": value["path"], "sha256": value["sha256"]}
            for name, value in dependencies.items()
        },
        "sources": {
            relative(Path(__file__).resolve()): file_hash(Path(__file__).resolve()),
            relative(CHECKER): file_hash(CHECKER),
            relative(TRANSITION_SCOUT): file_hash(TRANSITION_SCOUT),
            relative(ARC_SCOUT): file_hash(ARC_SCOUT),
        },
        "next_obligation": (
            "Construct the same direct Hardy-operator packets for the lower finite cell, ordinary "
            "complement, positive-real tail, and tiny target correction, then assemble complete K_T."
        ),
        "proof_boundary": (
            "Transition-plus-compact-upper-arc pair only on |t-10^10|<=10^-4. No complete K_T, "
            "wider Q_K-T sign interval, full event-cell theorem, wall handoff, all-height theorem, "
            "Lambda<=0, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified transition-upper-arc K pair on first subcell: "
        f"pair={pair.str(20, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
