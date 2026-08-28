#!/usr/bin/env python3
"""Promote the offset-20 grouped-752 phase transport and wider derivative box."""

from __future__ import annotations

import hashlib
import json
import os
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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_grouped_752_panel_phase_transport_scout as transport
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_offset20_local_box_gate as value_gate


complete = transport.complete
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_offset20_grouped_752_phase_transport_monotonicity_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
CENTER = "10000000020"
DERIVATIVE_RADIUS = "0.07"
LOSS_RADIUS = "0.1"
COMPONENT_PRECISION_BITS = 384
PRODUCTION_PANELS = 192
PRODUCTION_PANEL_PRECISION = 384
INDEPENDENT_PANELS = 256
INDEPENDENT_PANEL_PRECISION = 448
EXTENDED_LEFT = "10000000019.93"
RIGHT_ANCHOR = "10000000020.00012"
FIRST_SUBCELL_RIGHT = "10000000000.0001"
REMAINING_GAP = "19.9299"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def real(record: dict[str, Any]) -> arb:
    return complete.real_from(record)


def complex_ball(record: dict[str, Any]) -> acb:
    return complete.complex_from(record)


def panel_manifest(rows: list[dict[str, Any]]) -> dict[str, Any]:
    ordered = sorted(rows, key=lambda row: int(row["index"]))
    row_hashes = [transport.canonical_hash(row) for row in ordered]
    digest = hashlib.sha256()
    for row, row_hash in zip(ordered, row_hashes):
        digest.update(f"{row['index']}:{row_hash}\n".encode("ascii"))
    return {
        "count": len(ordered),
        "first_index": int(ordered[0]["index"]),
        "last_index": int(ordered[-1]["index"]),
        "ordered_canonical_rows_sha256": digest.hexdigest(),
    }


def select_panel_rows(
    panel_cache: dict[str, dict[str, Any]],
    *,
    panels: int,
    precision_bits: int,
    variant: str,
) -> list[dict[str, Any]]:
    return transport.selected_panels(
        panel_cache, CENTER, panels, precision_bits, variant
    )


def selected_component_rows(
    component_cache: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    selected: dict[str, Any] = {}
    for radius, variant, role in (
        (DERIVATIVE_RADIUS, "production", "theorem_production"),
        (DERIVATIVE_RADIUS, "independent", "theorem_independent"),
        (LOSS_RADIUS, "production", "production_loss"),
    ):
        for component in complete.COMPONENTS:
            key = complete.cache_key(
                CENTER, radius, COMPONENT_PRECISION_BITS, variant, component
            )
            require(key in component_cache, f"missing selected component row: {key}")
            row = component_cache[key]
            require(row.get("passed") is True, f"selected component row failed: {key}")
            selected[f"{role}_{component}"] = {
                "cache_key": key,
                "canonical_row_sha256": complete.canonical_hash(row),
            }
    return selected


def evaluate(
    panel_rows: list[dict[str, Any]],
    component_cache: dict[str, dict[str, Any]],
    *,
    radius: str,
    variant: str,
) -> dict[str, Any]:
    return transport.evaluate_model(
        panel_rows,
        component_cache,
        CENTER,
        radius,
        COMPONENT_PRECISION_BITS,
        variant,
    )


def assemble_certificate(
    panel_cache: dict[str, dict[str, Any]],
    component_cache: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    production_panels = select_panel_rows(
        panel_cache,
        panels=PRODUCTION_PANELS,
        precision_bits=PRODUCTION_PANEL_PRECISION,
        variant="production",
    )
    independent_panels = select_panel_rows(
        panel_cache,
        panels=INDEPENDENT_PANELS,
        precision_bits=INDEPENDENT_PANEL_PRECISION,
        variant="independent",
    )
    production = evaluate(
        production_panels,
        component_cache,
        radius=DERIVATIVE_RADIUS,
        variant="production",
    )
    independent = evaluate(
        independent_panels,
        component_cache,
        radius=DERIVATIVE_RADIUS,
        variant="independent",
    )
    loss = evaluate(
        production_panels,
        component_cache,
        radius=LOSS_RADIUS,
        variant="production",
    )

    production_derivative = real(
        production["replacement_Q_K_minus_T_derivative_component_sum"]
    )
    independent_derivative = real(
        independent["replacement_Q_K_minus_T_derivative_component_sum"]
    )
    require(production_derivative.lower() > 0, "production derivative lost positivity")
    require(independent_derivative.lower() > 0, "independent derivative lost positivity")
    require(
        production_derivative.overlaps(independent_derivative),
        "production/independent derivative boxes miss",
    )
    require(
        real(production["replacement_lower_ordinary_Q_derivative"]).overlaps(
            real(independent["replacement_lower_ordinary_Q_derivative"])
        ),
        "production/independent lower-ordinary projections miss",
    )
    require(
        complex_ball(production["replacement_complete_K_T"]).overlaps(
            complex_ball(independent["replacement_complete_K_T"])
        ),
        "production/independent complete K_T boxes miss",
    )

    loss_derivative = real(loss["replacement_Q_K_minus_T_derivative_component_sum"])
    require(loss_derivative.lower() < 0 < loss_derivative.upper(), "loss row does not straddle zero")
    require(
        loss["strictly_positive_replacement_derivative"] is False,
        "loss row unexpectedly remains positive",
    )

    derivative_box = arb(arb(CENTER), arb(DERIVATIVE_RADIUS))
    event = load_json(complete.EVENT_ATLAS)["certificate"]["primary_open_cell"]
    event_lower = arb(event["lower_ball"]["ball"])
    event_upper = arb(event["upper_ball"]["ball"])
    require(event_lower.upper() < derivative_box.lower(), "derivative box crosses lower event wall")
    require(derivative_box.upper() < event_upper.lower(), "derivative box crosses upper event wall")

    value_artifact = load_json(value_gate.RESULT)
    require(value_artifact.get("passed") is True, "offset-20 value gate did not pass")
    value = value_artifact["certificate"]
    value_production = value_gate.real_from(value["production_Q_K_minus_T"])
    value_independent = value_gate.real_from(value["independent_Q_K_minus_T"])
    require(value_production.upper() < 0, "production anchor lost negativity")
    require(value_independent.upper() < 0, "independent anchor lost negativity")
    extended_left = arb(EXTENDED_LEFT)
    anchor = arb(RIGHT_ANCHOR)
    require(derivative_box.lower() < extended_left.lower(), "extension leaves derivative box")
    require(anchor.upper() < derivative_box.upper(), "anchor leaves derivative box")
    value_left = arb(value["closed_height_interval"][0])
    value_right = arb(value["closed_height_interval"][1])
    require(value_left.lower() < anchor.lower(), "anchor leaves value box on left")
    require(anchor.upper() < value_right.upper(), "anchor leaves value box on right")
    gap = extended_left - arb(FIRST_SUBCELL_RIGHT)
    require(gap.contains(arb(REMAINING_GAP)), "remaining-gap arithmetic drift")

    production_transport_error = real(production["grouped_K_transport_error"])
    independent_transport_error = real(independent["grouped_K_transport_error"])
    require(production_transport_error.upper() < arb("6e-6"), "production transport error too large")
    require(independent_transport_error.upper() < arb("6e-6"), "independent transport error too large")

    return {
        "center_height": CENTER,
        "derivative_radius": DERIVATIVE_RADIUS,
        "closed_derivative_interval": [
            derivative_box.lower().str(55, more=True),
            derivative_box.upper().str(55, more=True),
        ],
        "strictly_inside_fixed_roster_event_cell": True,
        "event_cell_open_walls": {
            "lower": event_lower.str(65, more=True),
            "upper": event_upper.str(65, more=True),
        },
        "panel_transport_identity": {
            "base_packet": "F_t(y)=F_c(y)*exp(-i*h*log(y))",
            "common_phase": "exp(-i*h*s0)",
            "panel_phase": "exp(-i*h*(s_j-s0))",
            "base_panel_error": "r*epsilon_j*B_j",
            "log_panel_error": "log(C)*r*epsilon_j*B_j",
            "l1_majorant": "B_j=752*(y_right-y_left)/sqrt(y_left)",
            "stable_phase": "theta(t_c)+h*(theta'(D)-s0)",
            "justification": (
                "|exp(-i*h*s)-exp(-i*h*s_j)|<=|h|*|s-s_j| and the "
                "real mean-value theorem for theta"
            ),
        },
        "production": {
            "panels": PRODUCTION_PANELS,
            "panel_precision_bits": PRODUCTION_PANEL_PRECISION,
            "panel_manifest": panel_manifest(production_panels),
            "grouped_K_transport_error": production["grouped_K_transport_error"],
            "replacement_lower_ordinary_Q_derivative": production[
                "replacement_lower_ordinary_Q_derivative"
            ],
            "replacement_complete_K_T": production["replacement_complete_K_T"],
            "Q_K_minus_T_derivative": production[
                "replacement_Q_K_minus_T_derivative_component_sum"
            ],
            "broad_direct_projection_overlap": production["replacement_old_overlaps"],
        },
        "independent": {
            "panels": INDEPENDENT_PANELS,
            "panel_precision_bits": INDEPENDENT_PANEL_PRECISION,
            "panel_order": "reverse",
            "panel_manifest": panel_manifest(independent_panels),
            "grouped_K_transport_error": independent["grouped_K_transport_error"],
            "replacement_lower_ordinary_Q_derivative": independent[
                "replacement_lower_ordinary_Q_derivative"
            ],
            "replacement_complete_K_T": independent["replacement_complete_K_T"],
            "Q_K_minus_T_derivative": independent[
                "replacement_Q_K_minus_T_derivative_component_sum"
            ],
            "broad_direct_projection_overlap": independent["replacement_old_overlaps"],
        },
        "production_independent_derivative_overlap": True,
        "production_independent_lower_ordinary_overlap": True,
        "production_independent_complete_K_T_overlap": True,
        "value_anchor": {
            "source_gate": relative(value_gate.RESULT),
            "source_gate_sha256": file_hash(value_gate.RESULT),
            "right_anchor": RIGHT_ANCHOR,
            "production_Q_K_minus_T": value["production_Q_K_minus_T"],
            "independent_Q_K_minus_T": value["independent_Q_K_minus_T"],
        },
        "monotonicity_sign_extension": {
            "closed_interval": [EXTENDED_LEFT, RIGHT_ANCHOR],
            "identity": "Q(t)=Q(a)-integral_t^a (Q_K-T)'(u) du",
            "strictly_increasing_on_derivative_box": True,
            "Q_K_minus_T_strictly_negative_on_closed_interval": True,
            "first_subcell_right_endpoint": FIRST_SUBCELL_RIGHT,
            "remaining_open_gap": REMAINING_GAP,
        },
        "radius_diagnostic": {
            "selected_two_configuration_radius": DERIVATIVE_RADIUS,
            "first_tested_production_zero_straddling_radius": LOSS_RADIUS,
            "production_loss_derivative": loss[
                "replacement_Q_K_minus_T_derivative_component_sum"
            ],
            "loss_kind": "interval_enclosure_width_only",
            "pointwise_derivative_zero_inferred_from_loss": False,
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    p = c["production"]["Q_K_minus_T_derivative"]
    i = c["independent"]["Q_K_minus_T_derivative"]
    extension = c["monotonicity_sign_extension"]
    loss = c["radius_diagnostic"]
    return f"""# Offset-20 grouped-752 phase-transport monotonicity gate

Date: 2026-08-28

Status: rigorous two-configuration local derivative and sign-extension
certificate; not a proof of the all-height theorem or RH

Let `t_c=10^10+20`, `h=t-t_c`, and

```text
D_20^phase=[t_c-0.07,t_c+0.07].                     (PT1)
```

The old direct grouped-collar integrator evaluated 192 height-dependent panel
phases separately.  This forgot their common rotation.  On panel `j`, choose
`s_j` as the midpoint of its log-height interval and put
`s_0=(log L+log C)/2`.  For the exact grouped integrand,

```text
F_t(y)=F_(t_c)(y) exp(-i h log y)
      =exp(-i h s_0) exp(-i h(s_j-s_0))F_(t_c)(y)+E_j. (PT2)
```

If `epsilon_j=max|log y-s_j|` and
`B_j=752*(y_right-y_left)/sqrt(y_left)`, then the finite 752-term geometric
sum gives

```text
|E_j|<=|h| epsilon_j B_j,                           (PT3)
```

because `|exp(-i x)-exp(-i y)|<=|x-y|` for real `x,y`.  The log-weighted
integral has the companion bound `log(C)|h|epsilon_j B_j`.  The common Hardy
phase is enclosed without losing `h` by

```text
theta(t_c)+h(theta'(D_20^phase)-s_0),               (PT4)
```

which follows from the real mean-value theorem.  Thus (PT2)--(PT4) are a
rigorous representation change, not sampled interpolation.

The 192-panel, 384-bit production configuration gives

```text
(Q_K-T)' in {p['ball']},
lower endpoint = {p['lower']} > 0.                  (PT5)
```

A 256-panel, 448-bit independent configuration integrates panels in reverse
order and changes the underlying lower-ordinary, transition, endpoint, and arc
settings.  It gives

```text
(Q_K-T)' in {i['ball']},
lower endpoint = {i['lower']} > 0.                  (PT6)
```

The derivative, lower-ordinary, and complete complex `K_T` enclosures overlap.
Each phase-transport result also overlaps its deliberately broader direct
interval evaluation.  Hence `(Q_K-T)'>0` throughout `D_20^phase`.

The prior value gate gives `Q_K-T<0` at
`a=10^10+20+0.00012`.  Therefore

```text
Q_K(t)-T(t)<0 for every
t in [{extension['closed_interval'][0]},
      {extension['closed_interval'][1]}].            (PT7)
```

The exact remaining gap to the first subcell is
`{extension['remaining_open_gap']}` height units.  At radius
`{loss['first_tested_production_zero_straddling_radius']}`, the production
phase-transport enclosure

```text
{loss['production_loss_derivative']['ball']}         (PT8)
```

straddles zero.  This is an enclosure-width diagnostic only and proves no
pointwise derivative zero or sign change.

## Where pi comes from

`pi` is the ordinary circumference-to-diameter constant.  It enters through
the exact Fourier kernel `exp(2*pi*i*q*y)`, quadratic kernel
`exp(-i*pi*y^2)`, Gamma normalization, and `p=t/(2*pi)`.  It is not fitted
from a plot, polygon, or numerical pattern.

## Proof boundary

This gate proves the complete derivative only on `D_20^phase` and negativity
only on (PT7).  It does not bridge the remaining gap, cover the event cell,
cross an event wall, prove an all-height theorem, `Lambda<=0`, PF-infinity,
RH, or a Clay-prize conclusion.
"""


def main() -> int:
    ctx.prec = 512
    ctx.threads = 1
    require(CHECKER.is_file(), f"missing checker: {CHECKER}")
    panel_cache = transport.load_jsonl(transport.CACHE)
    component_cache = transport.load_jsonl(complete.CACHE)
    certificate = assemble_certificate(panel_cache, component_cache)

    dependencies = {
        "event_atlas": {
            "path": relative(complete.EVENT_ATLAS),
            "sha256": file_hash(complete.EVENT_ATLAS),
        },
        "offset20_value_gate": {
            "path": relative(value_gate.RESULT),
            "sha256": file_hash(value_gate.RESULT),
        },
    }
    source_paths = (
        Path(__file__).resolve(),
        CHECKER,
        Path(transport.__file__).resolve(),
        Path(complete.__file__).resolve(),
        Path(transport.upper_group.__file__).resolve(),
        Path(transport.base.__file__).resolve(),
        Path(value_gate.__file__).resolve(),
    )
    artifact = {
        "kind": "rh_rigorous_offset20_grouped_752_phase_transport_monotonicity_gate",
        "date": "2026-08-28",
        "status": "offset20_grouped_752_phase_transport_derivative_and_sign_extension_certified",
        "passed": True,
        "precision_bits": 512,
        "certificate": certificate,
        "selected_component_rows": selected_component_rows(component_cache),
        "panel_cache_path": relative(transport.CACHE),
        "component_cache_path": relative(complete.CACHE),
        "dependencies": dependencies,
        "sources": {relative(path): file_hash(path) for path in source_paths},
        "reproduction": {
            "production": (
                "python work/rh_compute/scripts/"
                f"{transport.STEM}.py --radius 0.07 --variant production"
            ),
            "independent": (
                "python work/rh_compute/scripts/"
                f"{transport.STEM}.py --radius 0.07 --variant independent"
            ),
            "checker": f"python work/rh_compute/scripts/check_{STEM}.py",
        },
        "decision": {
            "grouped_752_common_phase_transport_proved": True,
            "offset20_radius_0p07_complete_derivative_positive": True,
            "production_independent_reproduction_passed": True,
            "offset20_leftward_sign_extension_to_19p93_proved": True,
            "connected_interval_from_first_subcell_proved": False,
            "maximal_event_cell_sign_proved": False,
            "event_wall_handoff_proved": False,
            "all_height_transport_theorem_proved": False,
            "rh_implication": False,
        },
        "proof_boundary": (
            "Only the grouped-752 phase transport, the radius-0.07 offset-20 complete "
            "derivative, and the resulting left-anchored negative interval are promoted. "
            "The exact 19.9299-unit gap and every global continuation step remain open."
        ),
    }

    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    result_tmp = RESULT.with_suffix(RESULT.suffix + ".tmp")
    note_tmp = NOTE.with_suffix(NOTE.suffix + ".tmp")
    result_tmp.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    note_tmp.write_text(render_note(artifact), encoding="utf-8")
    os.replace(result_tmp, RESULT)
    os.replace(note_tmp, NOTE)
    print(
        "certified offset-20 grouped-752 phase transport: "
        f"production_lower={certificate['production']['Q_K_minus_T_derivative']['lower']} "
        f"independent_lower={certificate['independent']['Q_K_minus_T_derivative']['lower']} "
        f"extended_left={EXTENDED_LEFT}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
