#!/usr/bin/env python3
"""Certify the joined-transition contribution to K_T on the first subcell."""

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

from flint import acb, arb

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_transition_height_derivative_scout as scout


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_transition_K_first_subcell_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
SCOUT = SCRIPT_DIR / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_transition_height_derivative_scout.py"
)
DEPENDENCIES = {
    "height_derivative_identity": ROOT / "work" / "rh_compute" / "results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_height_derivative_identity_gate.json"
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
    source_K = c["source_K"]
    gamma_K = c["Gamma_K_direct"]
    transition_K = c["transition_K_direct"]
    projected = c["transition_Q_derivative_contribution"]
    return f"""# Joined transition K on the first height subcell

Date: 2026-08-28

Status: rigorous transition-component interval certificate; the complete
`K_T` is not yet enclosed

On

```text
I_1=[10^10-10^-4,10^10+10^-4],
```

write

```text
J_tr=I_tr-C_G D_(39853..39894),
K_tr=J_tr'+(i theta'-H'/H)J_tr.                    (TK1)
```

The finite transition source is differentiated before its alternating
boundary-jet collapse.  Since

```text
d_t F(y,t)=-i log(y)F(y,t),
```

the complete Hardy transport operator has the cancellation-stable source
weight

```text
(d_t+i theta'-H'/H)F
  =[i(theta'-log(y))-H'/H]F.                       (TK2)
```

The 42-mode Gamma block is likewise summed with its operator weight per mode:

```text
C_G sum_m [i(theta'-log(m))+C_G'/C_G-H'/H]m^(-1/2-it).
                                                               (TK3)
```

No separately widened `D'` and `i theta'D` balls enter the direct result.
The source and Gamma operator packets are

```text
source K: real {source_K['real_ball']}
          imag {source_K['imag_ball']},

Gamma K:  real {gamma_K['real_ball']}
          imag {gamma_K['imag_ball']}.
```

Their direct joined enclosure is

```text
K_tr: real {transition_K['real_ball']}
      imag {transition_K['imag_ball']}.             (TK4)
```

It overlaps the deliberately dependency-heavy assembly from `J_tr'` and
`(i theta'-H'/H)J_tr`.  After the common Hardy projection and division by
positive `H`, this component contributes

```text
Hardy_t[K_tr]/H={projected['ball']} > 0.             (TK5)
```

The independent checker raises precision from 384 to 448 bits, changes jet
order `14 -> 16`, derivative slabs `256 -> 384`, integral panels `48 -> 64`,
and logarithm-series terms `34 -> 40`.  It also repeats the altered five-label
direct finite-source comparison.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
integer Fourier spacing, Gamma normalization, Riemann--Siegel phase, or
`p=t/(2*pi)`.  No fitted geometric constant supplies `pi`.

Proof boundary: (TK1)--(TK5) certify only the transition contribution to
`K_T` on `I_1`.  They do not enclose the lower-cell, upper-arc, ordinary-packet,
or tiny target-correction contributions, hence do not prove a complete
nonzero-radius `K_T` bound, a wider `Q_K-T` sign interval, the full event-cell
sign theorem, an event-wall handoff, an all-height theorem, `Lambda<=0`, RH,
or a prize-level conclusion.
"""


def main() -> int:
    dependencies = {
        name: {
            "path": relative(path),
            "sha256": file_hash(path),
            "artifact": load_json(path),
        }
        for name, path in DEPENDENCIES.items()
    }
    require(
        dependencies["height_derivative_identity"]["artifact"]["decision"]
        ["exact_joined_first_height_derivative_proved"]
        is True,
        "height derivative identity drift",
    )
    require(
        dependencies["first_height_subcell"]["artifact"]["decision"]
        ["first_nonzero_radius_Q_K_minus_T_sign_theorem_proved"]
        is True,
        "first height subcell drift",
    )

    run = scout.build("0.0001", 384, "production")
    direct_K = acb(
        arb(run["transition_K_direct"]["real_ball"]),
        arb(run["transition_K_direct"]["imag_ball"]),
    )
    assembled_K = acb(
        arb(run["transition_K_assembled"]["real_ball"]),
        arb(run["transition_K_assembled"]["imag_ball"]),
    )
    projected = arb(run["transition_Q_derivative_contribution"]["ball"])
    require(direct_K.overlaps(assembled_K), "two transition K assemblies miss")
    require(direct_K.real.rad() < arb("1e-6"), "direct transition K real radius widened")
    require(direct_K.imag.rad() < arb("1e-6"), "direct transition K imag radius widened")
    require(projected.lower() > 0, "transition Hardy derivative contribution lost positivity")
    require(run["altered_five_label_fixture"]["overlap"] is True, "altered fixture failed")

    certificate = {
        key: run[key]
        for key in (
            "height_ball",
            "radius",
            "precision_bits",
            "variant",
            "configuration",
            "exact_identities",
            "source_K",
            "Gamma_K_direct",
            "Gamma_K_assembled",
            "direct_and_assembled_Gamma_K_overlap",
            "transition_K_direct",
            "transition_K_assembled",
            "direct_and_assembled_transition_K_overlap",
            "transition_Q_derivative_contribution",
            "altered_five_label_fixture",
            "source_K_certificate",
        )
    }
    artifact = {
        "kind": "rh_c1_hardy_joined_transition_K_first_height_subcell_gate",
        "date": "2026-08-28",
        "status": "joined_transition_K_first_nonzero_height_subcell_interval_certified",
        "passed": True,
        "certificate": certificate,
        "decision": {
            "transition_K_nonzero_radius_interval_proved": True,
            "transition_Hardy_derivative_contribution_strictly_positive": True,
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
            relative(SCOUT): file_hash(SCOUT),
        },
        "next_obligation": (
            "Apply the same Hardy operator before interval widening to the upper arc, lower cell, "
            "ordinary complement, and tiny target correction; then assemble and independently "
            "check the complete K_T on I_1."
        ),
        "proof_boundary": (
            "Transition contribution to K_T only on |t-10^10|<=10^-4. No complete K_T bound, "
            "wider Q_K-T sign interval, full event-cell theorem, wall handoff, all-height theorem, "
            "Lambda<=0, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified joined-transition K on first height subcell: "
        f"Hardy_contribution={projected.str(20, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
