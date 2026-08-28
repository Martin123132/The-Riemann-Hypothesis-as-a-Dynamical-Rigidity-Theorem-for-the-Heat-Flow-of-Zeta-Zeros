#!/usr/bin/env python3
"""Certify the first nonzero height subcell with Q_K-T<0."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from flint import arb

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_first_subcell_interval_scout as scout


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_first_nonzero_height_subcell_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "joined_derivative": ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate.json",
    "event_atlas": ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_local_height_event_atlas_gate.json",
    "saved_closure": ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_one_sided_closure_gate.json",
}
DONOR_SOURCES = (
    Path(scout.__file__).resolve(),
    Path(scout.lower_cell.__file__).resolve(),
    Path(scout.upper_arc.__file__).resolve(),
    Path(scout.transition.__file__).resolve(),
    Path(scout.endpoint.__file__).resolve(),
    Path(scout.ibp.__file__).resolve(),
    Path(scout.upper_transition.__file__).resolve(),
    Path(scout.upper_far.__file__).resolve(),
    Path(scout.complement.__file__).resolve(),
    Path(scout.joined.__file__).resolve(),
)

HEIGHT = arb("10000000000")
RADIUS = arb("0.0001")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    q = c["Q_K_minus_T_height_box"]
    widths = c["component_widths"]
    return f"""# First nonzero joined Q_K-T height subcell

Date: 2026-08-28

Status: rigorous first nonzero-radius height subcell with `Q_K-T<0`
certified; extension to the maximal event cell remains open.

The exact A-free identity from Section 11.508 is

```text
Q_K(t)-T(t)=Hardy_t[Y_T(t)]/H(t),
Y_T=S_W-H D_T.                                      (HS1)
```

On the fixed-roster interval

```text
I_1=[10^10-10^(-4),10^10+10^(-4)],                 (HS2)
```

the same packet is rebuilt before projection as

```text
Y_T=U_unowned+O_join+I_transition
    -C_G D_(39853..39894)+(C_G-H)D_(622..39894).    (HS3)
```

Every scalar height occurrence in the lower cell, upper quarter arc,
ordinary complementary packet, transition packet, finite Dirichlet blocks,
phase, and prefactor is evaluated with Arb on the full box.  The Section
11.507 event atlas proves that every roster and endpoint used by (HS3) is
unchanged on `I_1`.

Two large-expression dependency losses are removed without approximation.
First,

```text
log H(t)-log H(t_0)=integral_(t_0)^t H'(u)/H(u) du. (HS4)
```

The direct and duplication formulas for `H'/H` overlap on the box, giving the
stable prefactor enclosure

```text
H(I_1)={c['H_ball']['ball']}.                        (HS5)
```

Second, the upper-arc action has one maximizer `delta_*(t)`.  The exact
envelope identity

```text
d A_max(t)/dt=delta_*(t)                            (HS6)
```

transports its point value across the box and gives

```text
A_max(I_1)={c['stable_arc_action_ball']}.            (HS7)
```

All original polynomial, omitted-numerator, compact-tail, and positive-tail
error guards remain active.  Representative component radii are

```text
U_unowned: real {widths['unowned_packet']['real_radius']},
           imag {widths['unowned_packet']['imag_radius']},
O_join:    real {widths['ordinary_packet']['real_radius']},
           imag {widths['ordinary_packet']['imag_radius']},
transition join:
           real {widths['joined_transition']['real_radius']},
           imag {widths['joined_transition']['imag_radius']}.
```

After the single final Hardy projection and division by the positive `H` box,

```text
Q_K-T={q['ball']},
lower={q['lower']},
upper={q['upper']}<0.                               (HS8)
```

The strict upper-endpoint negativity margin is
`{c['strict_negativity_margin_ball']}`.  The independent checker covers the
same interval by two adjacent half-boxes at 448 bits and changes the upper-arc
configuration, grouped panel counts, transition jet order, ordinary-tail
recurrence order, slab counts, and clustering power.

Pi provenance: every `pi` in (HS1)--(HS8) comes from the inherited
Riemann-Siegel phase, Gaussian/Fresnel kernel, Fourier characters, gamma
normalization, or the already-certified event equations.  No fitted geometric
constant is introduced.

Proof boundary: a rigorous fixed-roster sign theorem only on the compact
interval (HS2).  This does not extend the sign to the rest of the maximal
event cell, cross either event wall, establish an all-height theorem, prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.perf_counter()
    for path in (*DEPENDENCIES.values(), *DONOR_SOURCES, CHECKER):
        require(path.is_file(), f"missing source or dependency: {path}")
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    for name, payload in dependencies.items():
        require(payload.get("passed") is True, f"dependency did not pass: {name}")

    direct = scout.build("0.0001", 384, variant="production")
    qkt = arb(direct["projection"]["Q_K_minus_T_height_box"]["ball"])
    require(qkt.upper() < 0, "production height box lost Q_K-T negativity")
    require(qkt.upper() < arb("-0.00004"), "production sign margin fell below guard")

    t_box = arb(direct["height_ball"])
    event_cell = dependencies["event_atlas"]["certificate"]["primary_open_cell"]
    event_lower = arb(event_cell["lower_ball"]["ball"])
    event_upper = arb(event_cell["upper_ball"]["ball"])
    require(event_lower.upper() < t_box.lower(), "subcell crosses the lower event wall")
    require(t_box.upper() < event_upper.lower(), "subcell crosses the upper event wall")

    saved_qkt = arb(
        dependencies["saved_closure"]["certificate"]["Q_K_minus_T_ball"]["ball"]
    )
    require(qkt.overlaps(saved_qkt), "production subcell misses the saved-height closure")
    action_guard = direct["components"]
    stable_action = action_guard["upper_arc_action_guard"][
        "stable_envelope_action_ball"
    ]

    certificate = {
        "center_height": "10000000000",
        "radius": "0.0001",
        "height_interval": [
            (HEIGHT - RADIUS).str(45, more=True),
            (HEIGHT + RADIUS).str(45, more=True),
        ],
        "height_ball": direct["height_ball"],
        "event_cell_lower_ball": event_cell["lower_ball"],
        "event_cell_upper_ball": event_cell["upper_ball"],
        "fixed_roster_event_cell_containment": True,
        "exact_A_free_identity": "Q_K-T=Hardy_t[S_W-H*D_T]/H",
        "interval_packet_identity": (
            "Y_T=U_unowned+O_join+I_transition-C_G*D_(39853..39894)"
            "+(C_G-H)*D_(622..39894)"
        ),
        "H_transport": direct["projection"]["H_transport"],
        "H_ball": direct["projection"]["H_ball"],
        "stable_arc_action_ball": stable_action,
        "component_widths": action_guard,
        "complete_joined_packet_ball": direct["components"]["complete_joined_packet"],
        "theta_ball": direct["projection"]["theta_ball"],
        "Q_K_minus_T_height_box": direct["projection"]["Q_K_minus_T_height_box"],
        "strict_negativity_margin_ball": (-qkt.upper()).str(45, more=True),
        "saved_height_overlap": True,
        "all_donor_error_guards_retained": True,
        "production_configuration": direct["configuration"],
        "precision_bits": direct["precision_bits"],
    }
    artifact = {
        "kind": "rh_c1_hardy_joined_QK_minus_T_first_nonzero_height_subcell_gate",
        "date": "2026-08-28",
        "status": "first_nonzero_fixed_roster_joined_QK_minus_T_negative_height_subcell_certified",
        "passed": True,
        "certificate": certificate,
        "decision": {
            "first_nonzero_radius_Q_K_minus_T_sign_theorem_proved": True,
            "A_free_join_retained_before_projection": True,
            "stable_H_transport_used": True,
            "stable_upper_arc_envelope_transport_used": True,
            "maximal_event_cell_sign_proved": False,
            "event_wall_handoffs_proved": False,
            "all_height_theorem": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            relative(path): file_hash(path)
            for path in (*DONOR_SOURCES, BUILDER, CHECKER)
        },
        "runtime": {
            "elapsed_seconds": time.perf_counter() - started,
            "process_priority": direct["resource_mode"],
            "workers": 1,
        },
        "next_obligation": (
            "Sharpen the transition and upper-arc interval dependency, then enlarge the "
            "joined sign box within the maximal event cell before proving the two exact wall handoffs."
        ),
        "proof_boundary": (
            "Rigorous Q_K-T<0 only for |t-10^10|<=10^-4 in the fixed roster. "
            "No sign theorem on the rest of the maximal event cell, wall handoff, all-height "
            "theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified first nonzero joined Q_K-T height subcell: "
        f"upper={qkt.upper().str(20, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
