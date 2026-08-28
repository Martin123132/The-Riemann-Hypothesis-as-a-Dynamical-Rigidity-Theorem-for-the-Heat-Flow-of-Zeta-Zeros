#!/usr/bin/env python3
"""Certify the complete K_T and its Hardy projection on I_1."""

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

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "complete_K_T_first_subcell_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
DEPENDENCIES = {
    "transition_upper_arc": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_transition_upper_arc_K_pair_first_subcell_gate.json"
    ),
    "lower_ordinary": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_lower_ordinary_K_first_subcell_gate.json"
    ),
    "positive_real_tail": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "positive_real_tail_K_first_subcell_gate.json"
    ),
    "tiny_target_correction": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "tiny_target_correction_K_first_subcell_gate.json"
    ),
    "height_derivative": ROOT / "work/rh_compute/results" / (
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
    return path.resolve().relative_to(ROOT).as_posix()


def complex_from(record: dict[str, Any]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def real_record(value: arb, digits: int = 70) -> dict[str, str]:
    return {
        "ball": value.str(digits, more=True),
        "lower": value.lower().str(55, more=True),
        "upper": value.upper().str(55, more=True),
    }


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(55, more=True),
    }


def assemble_certificate(dependencies: dict[str, dict[str, Any]]) -> dict[str, Any]:
    pair = dependencies["transition_upper_arc"]["certificate"]
    lower = dependencies["lower_ordinary"]["certificate"]
    tail = dependencies["positive_real_tail"]["certificate"]
    correction = dependencies["tiny_target_correction"]["certificate"]

    transition_K = complex_from(pair["transition_K_direct"])
    upper_arc_K = complex_from(pair["upper_arc_K_direct"])
    pair_K = transition_K + upper_arc_K
    lower_ordinary_K = complex_from(lower["joined_lower_ordinary_K"])
    correction_K = complex_from(correction["compact_weighted_mode_sum_K_corr"])
    tail_modulus = arb(tail["modulus_bound"]["upper"]).upper()
    tail_K = acb(arb(0, tail_modulus), arb(0, tail_modulus))

    # Join cancellation partners first. The tail modulus supplies a rigorous
    # centered rectangle for its otherwise unused phase.
    complete_K = lower_ordinary_K + pair_K + correction_K + tail_K

    theta = arb(lower["phase_and_prefactor"]["theta"]["ball"])
    H = arb(lower["phase_and_prefactor"]["H"]["ball"])
    require(H.lower() > 0, "H lost positivity")
    direct_projection = joined.hardy_projection(theta, complete_K) / H

    lower_projection = arb(
        lower["joined_lower_ordinary_Q_derivative_contribution"]["ball"]
    )
    pair_projection = arb(pair["paired_Hardy_derivative_contribution"]["ball"])
    correction_projection = arb(correction["Hardy_projection_over_H"]["ball"])
    tail_projection_bound = arb(
        tail["Hardy_projection_over_H_absolute_bound"]["upper"]
    ).upper()
    tail_projection = arb(0, tail_projection_bound)

    # Exact linearity of the common Hardy projection preserves the tighter
    # cancellation-aware component enclosures.
    projected_sum = (
        lower_projection
        + pair_projection
        + correction_projection
        + tail_projection
    )
    pair_absolute_upper = arb(
        pair["paired_Hardy_derivative_contribution"]["absolute_upper"]
    ).upper()
    conservative_projection = (
        lower_projection
        + correction_projection
        + arb(0, arb("1.3e-6") + tail_projection_bound)
    )

    require(pair_absolute_upper < arb("1.3e-6"), "pair envelope drift")
    require(direct_projection.lower() > 0, "direct complete K_T projection is not positive")
    require(projected_sum.lower() > 0, "componentwise complete projection is not positive")
    require(
        conservative_projection.lower() > 0,
        "envelope-only complete projection is not positive",
    )
    require(
        direct_projection.overlaps(projected_sum),
        "direct and componentwise complete projections do not overlap",
    )

    return {
        "height_ball": lower["height_ball"],
        "radius": lower["radius"],
        "precision_bits": 384,
        "transition_upper_arc_K": complex_record(pair_K),
        "lower_ordinary_K": complex_record(lower_ordinary_K),
        "positive_real_tail_K_rectangle": complex_record(tail_K),
        "positive_real_tail_K_modulus_upper": real_record(tail_modulus),
        "tiny_target_correction_K": complex_record(correction_K),
        "complete_K_T": complex_record(complete_K),
        "phase_and_prefactor": {
            "theta": real_record(theta),
            "H": real_record(H),
        },
        "component_Hardy_projections_over_H": {
            "lower_ordinary": real_record(lower_projection),
            "transition_upper_arc": real_record(pair_projection),
            "positive_real_tail_absolute_bound": real_record(tail_projection_bound),
            "tiny_target_correction": real_record(correction_projection),
        },
        "complete_Hardy_projection_over_H_direct": real_record(direct_projection),
        "complete_Hardy_projection_over_H_component_sum": real_record(projected_sum),
        "complete_Hardy_projection_over_H_pair_envelope": real_record(
            conservative_projection
        ),
        "direct_and_component_projection_overlap": True,
        "identities": {
            "complete_operator_split": (
                "K_T=K_(L+O)+(K_tr+K_U)+K_tail+K_corr"
            ),
            "height_derivative": "(Q_K-T)'=Hardy_t[K_T]/H",
            "projection_linearity": (
                "Hardy_t[K_T]/H=sum Hardy_t[K_component]/H"
            ),
            "tail_rectangle": (
                "|K_tail|<=B_tail implies Re(K_tail),Im(K_tail) in [-B_tail,B_tail]"
            ),
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    pieces = c["component_Hardy_projections_over_H"]
    return f"""# Complete K_T on the first height subcell

Date: 2026-08-28

Status: rigorous complete `K_T` and strictly positive joined height derivative
on the first nonzero-radius subcell; not a proof of any wider or global theorem

On `I_1=[10^10-10^-4,10^10+10^-4]`, exact split-contour ownership gives

```text
K_T=K_(L+O)+(K_tr+K_U)+K_tail+K_corr.              (KT1)
```

The first two parentheses are cancellation-preserving joins.  The positive-real
tail is inserted as a centered complex rectangle obtained from its rigorous
modulus bound, while the exact finite correction is retained with its phase.
Arb obtains the complete complex enclosure

```text
K_T: real {c['complete_K_T']['real_ball']}
     imag {c['complete_K_T']['imag_ball']}.         (KT2)
```

The exact A-free derivative identity is

```text
(Q_K-T)'=Hardy_t[K_T]/H.                            (KT3)
```

There are two production enclosures of its right side.  Direct projection of
(KT2) gives

```text
Hardy_t[K_T]/H={c['complete_Hardy_projection_over_H_direct']['ball']}.  (KT4)
```

Exact linearity also permits the four already stabilized projected packets to
be added without reintroducing the large shared phase dependency:

```text
lower plus ordinary : {pieces['lower_ordinary']['ball']}
transition plus arc : {pieces['transition_upper_arc']['ball']}
positive tail bound : {pieces['positive_real_tail_absolute_bound']['upper']}
tiny correction     : {pieces['tiny_target_correction']['ball']}

Hardy_t[K_T]/H={c['complete_Hardy_projection_over_H_component_sum']['ball']} > 0. (KT5)
```

As a deliberately coarser guard, the transition--arc midpoint is discarded
and only its certified absolute envelope `1.3e-6` is retained.  Even then,

```text
(Q_K-T)' in {c['complete_Hardy_projection_over_H_pair_envelope']['ball']} > 0. (KT6)
```

The independent checker changes precision and summation order, rebuilds the
phase and positive `H` transport, verifies direct/projected linearity on an
altered four-packet fixture, and repeats the envelope-only positivity test.

Pi provenance: every `pi` is inherited from the exact Mellin--Fresnel kernel,
Fourier spacing, Gamma normalization, Riemann--Siegel phase, or `p=t/(2*pi)`.
No circle, polygon, visual symmetry, or fitted geometric constant supplies
`pi`.

Proof boundary: (KT1)--(KT6) prove the complete nonzero-radius `K_T` enclosure
and `(Q_K-T)'>0` only on `I_1`.  The previously certified `Q_K-T<0` interval is
therefore increasing there, but no larger sign interval follows without a new
adaptive cover.  No full event-cell theorem, event-wall handoff, all-height
transport, equation-(4) global infinite-series identity, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    priority = joined.set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    ctx.prec = 384
    ctx.threads = 1
    require(CHECKER.is_file(), f"missing checker: {CHECKER}")

    dependencies: dict[str, dict[str, Any]] = {}
    dependency_records: dict[str, dict[str, str]] = {}
    for name, path in DEPENDENCIES.items():
        payload = load_json(path)
        require(payload.get("passed") is True, f"dependency did not pass: {name}")
        dependencies[name] = payload
        dependency_records[name] = {
            "path": relative(path),
            "sha256": file_hash(path),
        }

    require(
        dependencies["transition_upper_arc"]["decision"]
        ["transition_upper_arc_K_pair_interval_proved"]
        is True,
        "transition--upper-arc dependency drift",
    )
    require(
        dependencies["lower_ordinary"]["decision"]
        ["complete_lower_plus_ordinary_K_interval_proved"]
        is True,
        "lower-plus-ordinary dependency drift",
    )
    require(
        dependencies["positive_real_tail"]["decision"]
        ["positive_real_tail_K_interval_proved"]
        is True,
        "positive-real-tail dependency drift",
    )
    require(
        dependencies["tiny_target_correction"]["decision"]
        ["tiny_target_correction_K_interval_proved"]
        is True,
        "tiny-correction dependency drift",
    )
    require(
        dependencies["height_derivative"]["decision"]
        ["exact_joined_first_height_derivative_proved"]
        is True,
        "height-derivative identity drift",
    )
    require(
        dependencies["first_height_subcell"]["decision"]
        ["first_nonzero_radius_Q_K_minus_T_sign_theorem_proved"]
        is True,
        "first height-subcell sign theorem drift",
    )

    certificate = assemble_certificate(dependencies)
    artifact = {
        "kind": "rh_c1_hardy_complete_K_T_first_height_subcell_gate",
        "date": "2026-08-28",
        "status": "complete_K_T_and_positive_Q_K_minus_T_derivative_first_height_subcell_certified",
        "passed": True,
        "resource_mode": priority,
        "certificate": certificate,
        "decision": {
            "complete_K_T_interval_proved": True,
            "complete_K_T_Hardy_projection_interval_proved": True,
            "Q_K_minus_T_derivative_strictly_positive_on_I1": True,
            "Q_K_minus_T_strictly_increasing_on_I1": True,
            "first_nonzero_radius_Q_K_minus_T_sign_theorem_preserved": True,
            "wider_Q_K_minus_T_sign_interval_proved": False,
            "maximal_event_cell_sign_proved": False,
            "event_wall_handoff_proved": False,
            "all_height_transport_theorem_proved": False,
            "rh_implication": False,
        },
        "dependencies": dependency_records,
        "sources": {
            relative(Path(__file__).resolve()): file_hash(Path(__file__).resolve()),
            relative(CHECKER): file_hash(CHECKER),
        },
        "next_obligation": (
            "Use the complete positive K_T derivative enclosure to grow an adaptive, "
            "deterministic Q_K-T sign cover from I_1 across the maximal same-roster "
            "event cell, then prove exact one-mode handoffs at its two walls."
        ),
        "proof_boundary": (
            "Complete complex K_T and complete Hardy projection only on "
            "|t-10^10|<=10^-4. This proves Q_K-T is strictly increasing while remaining "
            "negative on that already-certified subcell. It does not prove a wider sign "
            "interval, maximal event-cell theorem, wall handoff, all-height transport, "
            "equation-(4) global infinite-series identity, Lambda<=0, PF-infinity, RH, "
            "or a prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    projection = certificate["complete_Hardy_projection_over_H_component_sum"]["ball"]
    print(
        "certified complete K_T and positive joined derivative on first subcell; "
        f"Hardy_projection={projection}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
