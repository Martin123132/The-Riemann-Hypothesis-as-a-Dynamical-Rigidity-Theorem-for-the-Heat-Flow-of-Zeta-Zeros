#!/usr/bin/env python3
"""Certify the lower finite-cell contribution to K_T on the first subcell."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_finite_cell_K_first_subcell_scout as scout


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_lower_finite_cell_K_first_subcell_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
SCOUT = SCRIPT_DIR / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_lower_finite_cell_K_first_subcell_scout.py"
)
DEPENDENCIES = {
    "lower_cell_value": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
        "lower_cell_reversed_roster_endpoint_gate.json"
    ),
    "height_derivative_identity": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_height_derivative_identity_gate.json"
    ),
    "first_height_subcell": ROOT / "work/rh_compute/results" / (
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


def complex_from(record: dict[str, Any]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    direct = c["lower_cell_K_direct"]
    assembled = c["lower_cell_K_assembled"]
    projected = c["lower_cell_Q_derivative_contribution"]
    remainder = c["direct_certificate"]["total_remainder_bound_ball"]
    return f"""# Joined lower finite-cell K on the first height subcell

Date: 2026-08-28

Status: rigorous lower finite-cell component certificate; the ordinary
complements and complete `K_T` are not yet enclosed

On

```text
I_1=[10^10-10^-4,10^10+10^-4],
```

write the lower finite cell as

```text
V_L(t)=integral_0^L F(y,t)dy,       L=621.5,
K_L=(d_t+i theta'-H'/H)V_L.                         (LK1)
```

Since `d_t F=-i log(y)F`, the complete Hardy operator is inserted before the
reversed 2481423-label endpoint collapse:

```text
K_L=integral_0^L [i(theta'-log(y))-H'/H]F(y,t)dy.   (LK2)
```

For `d_r=p/y+y-q_++r` and `T(f)=(f/d_r)'`, exact linearity gives

```text
T^n[(C-i log(y))y^(-1/2)]
 =C T^n[y^(-1/2)]-i T^n[log(y)y^(-1/2)],
C=i theta'-H'/H.                                   (LK3)
```

The second recurrence contains only monomials
`p^j y^(e/2)(log y)^ell d_r^(-k)` with `ell=0,1`.  Near zero its discarded
integral uses the exact elementary moments

```text
integral_0^S y^(b-1)dy=S^b/b,
integral_0^S y^(b-1)|log y|dy
 =S^b(-log(S)/b+1/b^2),             0<S<=1,
 =S^b(log(S)/b-1/b^2)+2/b^2,        S>=1.          (LK4)
```

On the upper slabs the complete signed complex coefficient of each
denominator power is interval-evaluated before its modulus is taken.  At
order 18, split `300`, and 4096 slabs the direct remainder is

```text
{remainder['ball']}.
```

The direct operator enclosure is

```text
K_L: real {direct['real_ball']}
     imag {direct['imag_ball']}.                   (LK5)
```

It overlaps the deliberately wider assembly

```text
V_L'+(i theta'-H'/H)V_L:
 real {assembled['real_ball']}
 imag {assembled['imag_ball']}.
```

After the common Hardy projection and stable positive `H` transport, the
lower finite-cell component contributes

```text
Hardy_t[K_L]/H={projected['ball']} < 0.             (LK6)
```

The independent checker raises precision from 384 to 448 bits, changes the
endpoint order `18 -> 20`, split `300 -> 280`, and slab count `4096 -> 6144`.
It repeats a changed three-label direct logarithmic-chart quadrature and
requires direct and assembled operator forms to overlap.

The next cancellation-preserving coordinate is

```text
K_(L+O)=(d_t+i theta'-H'/H)
        [V_L-T_lower-T_upper+Gamma_defect].         (LK7)
```

At `L`, the lower complement, finite source roster, and upper complement
partition the complete half-integer lattice.  Their endpoint currents must
therefore be joined before a final norm; (LK6) is not substituted for that
future joined calculation.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
half-integer Fourier spacing, Riemann--Siegel phase, Gamma normalization, or
`p=t/(2*pi)`.  No fitted geometric constant supplies `pi`.

Proof boundary: (LK1)--(LK6) certify only the lower finite-cell contribution
to `K_T` on `I_1`.  They do not enclose either ordinary complementary tail,
the complete lower-plus-ordinary join (LK7), the positive-real tail
derivative, tiny target correction, complete `K_T`, a wider `Q_K-T` sign
interval, the full event cell, wall handoffs, an all-height theorem,
`Lambda<=0`, RH, or a prize-level conclusion.
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
        dependencies["lower_cell_value"]["artifact"]["decision"]
        ["lower_cell_complex_ball_certified"]
        is True,
        "lower-cell value dependency drift",
    )
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
    direct = complex_from(run["lower_cell_K_direct"])
    assembled = complex_from(run["lower_cell_K_assembled"])
    projected = arb(run["lower_cell_Q_derivative_contribution"]["ball"])
    remainder = arb(run["direct_certificate"]["total_remainder_bound_ball"]["ball"])
    require(direct.overlaps(assembled), "direct and assembled lower-cell K miss")
    require(direct.real.rad() < arb("2e-8"), "lower-cell K real radius widened")
    require(direct.imag.rad() < arb("5e-8"), "lower-cell K imag radius widened")
    require(remainder.upper() < arb("1e-14"), "lower-cell weighted recurrence remainder widened")
    require(projected.upper() < 0, "lower finite-cell Hardy contribution lost negativity")
    fixture = run["altered_three_label_fixture"]
    require(fixture["overlap"] is True, "altered three-label fixture failed")
    require(float(fixture["direct_to_recurrence_discrepancy"]) < 1e-7, "fixture midpoint drift")

    certificate = {
        key: run[key]
        for key in (
            "height_ball",
            "radius",
            "precision_bits",
            "variant",
            "configuration",
            "exact_identities",
            "phase_and_prefactor",
            "lower_cell_value",
            "lower_cell_time_derivative",
            "lower_cell_K_direct",
            "lower_cell_K_assembled",
            "direct_and_assembled_overlap",
            "lower_cell_Q_derivative_contribution",
            "direct_certificate",
            "time_derivative_certificate",
            "altered_three_label_fixture",
        )
    }
    artifact = {
        "kind": "rh_c1_hardy_joined_lower_finite_cell_K_first_height_subcell_gate",
        "date": "2026-08-28",
        "status": "joined_lower_finite_cell_K_first_nonzero_height_subcell_interval_certified",
        "passed": True,
        "certificate": certificate,
        "decision": {
            "lower_finite_cell_K_nonzero_radius_interval_proved": True,
            "lower_finite_cell_Hardy_derivative_contribution_strictly_negative": True,
            "complete_lower_plus_ordinary_K_interval_proved": False,
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
            "Insert the same Hardy operator into the lower and upper ordinary complements, "
            "join all complete-half-lattice endpoint currents at L before a final norm, and "
            "independently certify K_(L+O) on I_1."
        ),
        "proof_boundary": (
            "Lower finite-cell contribution to K_T only on |t-10^10|<=10^-4. No complete "
            "lower-plus-ordinary K, complete K_T, wider Q_K-T sign interval, full event-cell "
            "theorem, wall handoff, all-height theorem, Lambda<=0, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified joined lower finite-cell K on first height subcell: "
        f"Hardy_contribution={projected.str(20, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
