#!/usr/bin/env python3
"""Certify the joined lower-cell plus ordinary K packet on the first subcell."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_ordinary_K_first_subcell_scout as scout


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_lower_ordinary_K_first_subcell_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
SCOUT = SCRIPT_DIR / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_lower_ordinary_K_first_subcell_scout.py"
)
DEPENDENCIES = {
    "lower_finite_K": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_lower_finite_cell_K_first_subcell_gate.json"
    ),
    "ordinary_value": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
        "lower_complement_signed_coefficient_interval_gate.json"
    ),
    "upper_group": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
        "upper_complement_closed_752_label_exterior_transition_layer_gate.json"
    ),
    "upper_far": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
        "upper_complement_closed_752_far_remainder_eight_round_gate.json"
    ),
    "endpoint_algebra": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
        "nonstationary_complement_root_unity_hurwitz_endpoint_gate.json"
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
    lower = c["lower_finite_cell_K"]
    ordinary_K = c["ordinary_K"]
    joined_K = c["joined_lower_ordinary_K"]
    ordinary_projection = c["ordinary_Q_derivative_contribution"]
    joined_projection = c["joined_lower_ordinary_Q_derivative_contribution"]
    return f"""# Joined lower-cell plus ordinary K on the first height subcell

Date: 2026-08-28

Status: rigorous lower-plus-ordinary component certificate; complete `K_T`
is not yet enclosed

On `I_1=[10^10-10^-4,10^10+10^-4]`, retain

```text
O_join=-T_lower-T_upper+Gamma_defect,
T_upper=G_752+F_C+I_L.                              (LOK1)
```

With `L_t=d_t+i theta'-H'/H`, the cancellation-preserving target is

```text
K_(L+O)=L_t[V_L-T_lower-G_752-F_C-I_L+Gamma_defect].
                                                               (LOK2)
```

Every endpoint expansion uses the exact weighted recurrence

```text
T^n[(C-i log y)y^(-1/2)]
 =C T^n[y^(-1/2)]-i T^n[(log y)y^(-1/2)],
C=i theta'-H'/H,       T(f)=(f/(q-y-p/y))'.        (LOK3)
```

For each denominator power, its complete signed complex coefficient is
interval-evaluated on a slab before taking a modulus.  The finite source,
lower complement, and remote upper complement endpoint partials are then
joined as complex balls before any recurrence remainder is added.  The
752-label endpoint-saddle collar is not differentiated after integration:
the weight `i(theta'-log y)-H'/H` is inserted inside its one grouped
geometric-amplitude integral.  Its exact zeros at `L` and `C` survive this
multiplication.

The retained component enclosures are

```text
K_L: real {lower['real_ball']}
     imag {lower['imag_ball']},

K_O: real {ordinary_K['real_ball']}
     imag {ordinary_K['imag_ball']},

K_(L+O): real {joined_K['real_ball']}
         imag {joined_K['imag_ball']}.              (LOK4)
```

The explicit sum of all endpoint-recurrence remainders is
`{c['joined_lower_ordinary_total_error']['ball']}`.  The Gamma defect is
retained through the exact rationalized expression for `1-C_G`; its Hardy
operator radius remains astronomically small.

After the common Hardy projection and stable positive `H` transport,

```text
Hardy_t[K_O]/H={ordinary_projection['ball']},
Hardy_t[K_(L+O)]/H={joined_projection['ball']} > 0. (LOK5)
```

The independent checker raises precision from 384 to 448 bits, changes the
lower finite-cell order, split, and slabs, changes the lower-complement order
`3 -> 4`, changes upper-tail order `8 -> 9`, changes all clustered slab
schemes, raises endpoint powers `16 -> 20`, and changes grouped panels
`192 -> 256`.  It also constructs changed four-label direct quadratures for
both signs of `q-y-p/y` and compares them with the weighted endpoint
recurrence.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
half-integer Fourier spacing, root-of-unity endpoint algebra,
Riemann--Siegel phase, Gamma normalization, or `p=t/(2*pi)`.  No fitted
geometric constant supplies `pi`.

Proof boundary: (LOK1)--(LOK5) certify only the joined lower finite-cell plus
ordinary contribution on `I_1`.  They omit the positive-real tail derivative
and tiny target correction and are not yet combined with the certified
transition--upper-arc pair.  Therefore they prove no complete `K_T` bound,
wider `Q_K-T` sign interval, full event-cell theorem, wall handoff,
all-height theorem, `Lambda<=0`, RH, or prize-level conclusion.
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
        dependencies["lower_finite_K"]["artifact"]["decision"]
        ["lower_finite_cell_K_nonzero_radius_interval_proved"]
        is True,
        "lower finite-cell K dependency drift",
    )
    require(
        dependencies["ordinary_value"]["artifact"]["decision"]
        ["complete_ordinary_packet_complex_ball_certified"]
        is True,
        "ordinary value dependency drift",
    )
    require(
        dependencies["upper_group"]["artifact"]["decision"]
        ["grouped_752_label_transition_enclosed"]
        is True,
        "upper grouped dependency drift",
    )
    require(
        dependencies["upper_far"]["artifact"]["decision"]
        ["complete_upper_complement_complex_ball_certified"]
        is True,
        "upper far dependency drift",
    )
    require(
        dependencies["endpoint_algebra"]["artifact"]["decision"]
        ["production_endpoint_denominator_sums_certified"]
        is True,
        "endpoint algebra dependency drift",
    )
    require(
        dependencies["first_height_subcell"]["artifact"]["decision"]
        ["first_nonzero_radius_Q_K_minus_T_sign_theorem_proved"]
        is True,
        "first height subcell drift",
    )

    run = scout.build("0.0001", 384, "production")
    joined_K = complex_from(run["joined_lower_ordinary_K"])
    component_sum = complex_from(run["lower_finite_cell_K"]) + complex_from(
        run["ordinary_K"]
    )
    projected = arb(run["joined_lower_ordinary_Q_derivative_contribution"]["ball"])
    explicit_error = arb(run["joined_lower_ordinary_total_error"]["ball"])
    require(joined_K.overlaps(component_sum), "joined/component lower-ordinary K miss")
    require(run["component_sum_overlap"] is True, "component overlap flag lost")
    require(explicit_error.upper() < arb("1e-5"), "joined explicit remainder widened")
    require(joined_K.real.rad() < arb("3e-5"), "joined real radius widened")
    require(joined_K.imag.rad() < arb("3e-5"), "joined imag radius widened")
    require(projected.lower() > 0, "joined lower-ordinary contribution lost positivity")
    require(len(run["complete_lattice_L_endpoint_currents"]) == 3, "L-current census drift")

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
            "lower_finite_cell_K",
            "lower_complement_K",
            "grouped_752_K",
            "finite_upper_K",
            "remote_upper_K",
            "Gamma_defect_K",
            "ordinary_K",
            "ordinary_Q_derivative_contribution",
            "joined_lower_ordinary_endpoint_partial",
            "joined_lower_ordinary_total_error",
            "joined_lower_ordinary_K",
            "component_sum_overlap",
            "joined_lower_ordinary_Q_derivative_contribution",
            "complete_lattice_L_endpoint_currents",
        )
    }
    artifact = {
        "kind": "rh_c1_hardy_joined_lower_ordinary_K_first_height_subcell_gate",
        "date": "2026-08-28",
        "status": "joined_lower_plus_ordinary_K_first_nonzero_height_subcell_interval_certified",
        "passed": True,
        "certificate": certificate,
        "decision": {
            "joined_lower_plus_ordinary_K_nonzero_radius_interval_proved": True,
            "complete_lower_plus_ordinary_K_interval_proved": True,
            "joined_lower_plus_ordinary_Hardy_contribution_strictly_positive": True,
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
            "Certify the positive-real tail Hardy derivative and exact tiny target correction, "
            "then assemble this lower-plus-ordinary packet with the transition--upper-arc pair "
            "and independently check complete K_T on I_1."
        ),
        "proof_boundary": (
            "Joined lower finite-cell plus ordinary contribution only on |t-10^10|<=10^-4. "
            "No positive-real tail derivative, tiny target correction, complete K_T, wider "
            "Q_K-T sign interval, full event-cell theorem, wall handoff, all-height theorem, "
            "Lambda<=0, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified joined lower-plus-ordinary K on first height subcell: "
        f"Hardy_contribution={projected.str(20, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
