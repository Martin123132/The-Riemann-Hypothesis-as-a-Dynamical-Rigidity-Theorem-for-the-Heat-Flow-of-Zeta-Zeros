#!/usr/bin/env python3
"""Promote the offset-20 complete K_T derivative and its left sign extension."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_local_complete_K_T_scout as scout
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_offset20_local_box_gate as value_gate


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_offset20_complete_K_T_monotonicity_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
CACHE = scout.CACHE
CENTER = "10000000020"
OFFSET = "20"
DERIVATIVE_RADIUS = "0.006"
PRODUCTION_EDGE_RADIUS = "0.0061"
PRODUCTION_LOSS_RADIUS = "0.01"
VALUE_RADIUS = value_gate.PASS_RADIUS
PRECISION_BITS = 384
FIRST_SUBCELL_RIGHT = "10000000000.0001"
EXTENDED_LEFT = "10000000019.994"
RIGHT_ANCHOR = "10000000020.00012"
REMAINING_GAP = "19.9939"


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


def load_cache(path: Path = CACHE) -> dict[str, dict[str, Any]]:
    require(path.is_file(), f"missing component cache: {path}")
    rows: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            key = row.get("cache_key")
            require(key, f"component cache row {line_number} has no key")
            require(key not in rows, f"duplicate component cache key at row {line_number}")
            rows[key] = row
    return rows


def derivative_record(assembly: dict[str, Any], route: str) -> dict[str, Any]:
    return assembly[f"Q_K_minus_T_derivative_{route}"]


def real(record: dict[str, Any]) -> arb:
    return scout.real_from(record)


def complex_ball(record: dict[str, Any]) -> acb:
    return scout.complex_from(record)


def interval_record(value: arb) -> dict[str, str]:
    return scout.real_record(value)


def assembly(rows: dict[str, dict[str, Any]], radius: str, variant: str) -> dict[str, Any]:
    return scout.assemble(rows, CENTER, radius, PRECISION_BITS, variant)


def selected_rows(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    selections: dict[str, Any] = {}
    for radius, variant, role in (
        (DERIVATIVE_RADIUS, "production", "theorem_production"),
        (DERIVATIVE_RADIUS, "independent", "theorem_independent"),
        (PRODUCTION_EDGE_RADIUS, "production", "production_edge"),
        (PRODUCTION_LOSS_RADIUS, "production", "production_loss"),
    ):
        for component in scout.COMPONENTS:
            key = scout.cache_key(CENTER, radius, PRECISION_BITS, variant, component)
            require(key in rows, f"missing selected component row: {key}")
            row = rows[key]
            require(row.get("passed") is True, f"selected component row failed: {key}")
            selections[f"{role}_{component}"] = {
                "cache_key": key,
                "canonical_row_sha256": scout.canonical_hash(row),
            }
    return selections


def assemble_certificate(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    production = assembly(rows, DERIVATIVE_RADIUS, "production")
    independent = assembly(rows, DERIVATIVE_RADIUS, "independent")
    edge = assembly(rows, PRODUCTION_EDGE_RADIUS, "production")
    loss = assembly(rows, PRODUCTION_LOSS_RADIUS, "production")

    for name, row in (("production", production), ("independent", independent)):
        require(row["strictly_positive_derivative"] is True, f"{name} derivative lost positivity")
        direct = real(derivative_record(row, "direct"))
        components = real(derivative_record(row, "component_sum"))
        require(direct.lower() > 0, f"{name} direct projection is not positive")
        require(components.lower() > 0, f"{name} component projection is not positive")
        require(direct.overlaps(components), f"{name} direct/component projections miss")

    production_direct = real(derivative_record(production, "direct"))
    production_components = real(derivative_record(production, "component_sum"))
    independent_direct = real(derivative_record(independent, "direct"))
    independent_components = real(derivative_record(independent, "component_sum"))
    require(production_direct.overlaps(independent_direct), "production/independent direct projections miss")
    require(
        production_components.overlaps(independent_components),
        "production/independent component projections miss",
    )
    require(
        complex_ball(production["complete_K_T"]).overlaps(complex_ball(independent["complete_K_T"])),
        "production/independent complete K_T boxes miss",
    )

    edge_direct = real(derivative_record(edge, "direct"))
    edge_components = real(derivative_record(edge, "component_sum"))
    require(edge["strictly_positive_derivative"] is True, "production edge row lost positivity")
    require(edge_direct.lower() > 0 and edge_components.lower() > 0, "edge row is not positive")

    loss_direct = real(derivative_record(loss, "direct"))
    loss_components = real(derivative_record(loss, "component_sum"))
    require(loss["strictly_positive_derivative"] is False, "loss row unexpectedly stays positive")
    require(loss_direct.lower() < 0 < loss_direct.upper(), "direct loss row does not straddle zero")
    require(
        loss_components.lower() < 0 < loss_components.upper(),
        "component loss row does not straddle zero",
    )

    event = load_json(scout.EVENT_ATLAS)["certificate"]["primary_open_cell"]
    event_lower = arb(event["lower_ball"]["ball"])
    event_upper = arb(event["upper_ball"]["ball"])
    derivative_box = arb(arb(CENTER), arb(DERIVATIVE_RADIUS))
    require(event_lower.upper() < derivative_box.lower(), "derivative box crosses lower event wall")
    require(derivative_box.upper() < event_upper.lower(), "derivative box crosses upper event wall")

    value_artifact = load_json(value_gate.RESULT)
    require(value_artifact.get("passed") is True, "offset-20 value gate did not pass")
    require(
        value_artifact.get("status") == "offset20_closed_local_Q_K_minus_T_negative_interval_certified",
        "offset-20 value-gate status drift",
    )
    value = value_artifact["certificate"]
    require(value["center_height"] == CENTER, "value anchor center drift")
    require(value["radius"] == VALUE_RADIUS, "value anchor radius drift")
    value_production = value_gate.real_from(value["production_Q_K_minus_T"])
    value_independent = value_gate.real_from(value["independent_Q_K_minus_T"])
    require(value_production.upper() < 0, "production value anchor lost negativity")
    require(value_independent.upper() < 0, "independent value anchor lost negativity")

    extended_left = arb(EXTENDED_LEFT)
    anchor = arb(RIGHT_ANCHOR)
    require(derivative_box.lower() < extended_left.lower(), "left extension leaves derivative box")
    require(anchor.upper() < derivative_box.upper(), "right anchor leaves derivative box")
    value_left = arb(value["closed_height_interval"][0])
    value_right = arb(value["closed_height_interval"][1])
    require(value_left.lower() < anchor.lower(), "right anchor leaves value box on left")
    require(anchor.upper() < value_right.upper(), "right anchor leaves value box on right")
    gap = extended_left - arb(FIRST_SUBCELL_RIGHT)
    require(gap.contains(arb(REMAINING_GAP)), "remaining gap arithmetic drift")
    require(gap.lower() > 0, "offset-20 extension unexpectedly reaches the first subcell")

    return {
        "offset_from_t0": OFFSET,
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
        "production_derivative_direct": interval_record(production_direct),
        "production_derivative_component_sum": interval_record(production_components),
        "independent_derivative_direct": interval_record(independent_direct),
        "independent_derivative_component_sum": interval_record(independent_components),
        "production_independent_direct_overlap": True,
        "production_independent_component_overlap": True,
        "production_independent_complete_K_T_overlap": True,
        "value_anchor": {
            "source_gate": relative(value_gate.RESULT),
            "source_gate_sha256": file_hash(value_gate.RESULT),
            "closed_value_interval": value["closed_height_interval"],
            "right_anchor": RIGHT_ANCHOR,
            "production_Q_K_minus_T": value["production_Q_K_minus_T"],
            "independent_Q_K_minus_T": value["independent_Q_K_minus_T"],
        },
        "monotonicity_sign_extension": {
            "closed_interval": [
                EXTENDED_LEFT,
                RIGHT_ANCHOR,
            ],
            "identity": "Q(t)=Q(a)-integral_t^a (Q_K-T)'(u) du",
            "anchor": "a=10^10+20+0.00012",
            "strictly_increasing_on_derivative_box": True,
            "Q_K_minus_T_strictly_negative_on_closed_interval": True,
            "first_subcell_right_endpoint": FIRST_SUBCELL_RIGHT,
            "remaining_open_gap": REMAINING_GAP,
        },
        "radius_diagnostic": {
            "selected_two_configuration_radius": DERIVATIVE_RADIUS,
            "largest_tested_production_positive_radius": PRODUCTION_EDGE_RADIUS,
            "production_edge_direct": interval_record(edge_direct),
            "production_edge_component_sum": interval_record(edge_components),
            "first_tested_production_zero_straddling_radius": PRODUCTION_LOSS_RADIUS,
            "production_loss_direct": interval_record(loss_direct),
            "production_loss_component_sum": interval_record(loss_components),
            "loss_kind": "interval_enclosure_width_only",
            "pointwise_derivative_zero_inferred_from_loss": False,
        },
        "identities": {
            "complete_derivative": "(Q_K-T)'=Hardy_t[K_T]/H",
            "complete_ownership": "K_T=K_(L+O)+(K_tr+K_U)+K_tail+K_corr",
            "height_coordinate": "p=t/(2*pi)",
            "Fourier_kernel": "exp(2*pi*i*q*y)",
            "quadratic_kernel": "exp(-i*pi*y^2)",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    dp = c["production_derivative_component_sum"]
    di = c["independent_derivative_component_sum"]
    extension = c["monotonicity_sign_extension"]
    diagnostic = c["radius_diagnostic"]
    return f"""# Offset-20 complete `K_T` monotonicity bridge

Date: 2026-08-28

Status: rigorous two-configuration local derivative certificate and leftward
sign extension; not a proof of the all-height theorem or RH

Let

```text
D_20=[10^10+20-0.006,10^10+20+0.006].              (MB1)
```

The complete ownership identity is

```text
K_T=K_(L+O)+(K_tr+K_U)+K_tail+K_corr,
(Q_K-T)'=Hardy_t[K_T]/H.                            (MB2)
```

All five components are evaluated on the same Arb height box before the final
projection.  The production component sum gives

```text
(Q_K-T)' in {dp['ball']},
lower endpoint = {dp['lower']} > 0.                 (MB3)
```

An independent configuration reverses the ordinary expansion path and changes
the duplication, transition, endpoint, and upper-arc settings.  It gives

```text
(Q_K-T)' in {di['ball']},
lower endpoint = {di['lower']} > 0.                 (MB4)
```

The direct and component-sum projections overlap within both configurations;
the production and independent direct projections, component projections, and
complete complex `K_T` boxes also overlap.  Hence `Q_K-T` is strictly increasing
throughout `D_20`.

## Sign extension

The earlier value gate proves `Q_K-T<0` on
`I_20=[10^10+20-0.00012,10^10+20+0.00012]`.  Put
`a=10^10+20+0.00012`.  For every `t` between the left endpoint of `D_20` and
`a`, the fundamental theorem of calculus gives

```text
Q_K(t)-T(t)
 =Q_K(a)-T(a)-integral_t^a (Q_K-T)'(u) du
 <=Q_K(a)-T(a)<0.                                   (MB5)
```

Therefore

```text
Q_K(t)-T(t)<0 for every
t in [{extension['closed_interval'][0]},
      {extension['closed_interval'][1]}].            (MB6)
```

This extends the offset-20 sign island leftward by `0.00588` height units.  The
remaining open gap to the right endpoint `10000000000.0001` of the first
subcell is `{extension['remaining_open_gap']}` height units.

## Radius diagnostic

Production remains positive at radius
`{diagnostic['largest_tested_production_positive_radius']}`, but at radius
`{diagnostic['first_tested_production_zero_straddling_radius']}` its direct and
component enclosures both straddle zero.  This is an interval-enclosure width
obstruction only; it proves no pointwise derivative zero or sign change.

## Where pi comes from

Here `pi` is the ordinary circle constant, circumference divided by diameter.
It enters through the exact Fourier kernel `exp(2*pi*i*q*y)`, the quadratic
kernel `exp(-i*pi*y^2)`, Gamma normalization, and the height coordinate
`p=t/(2*pi)`.  It is not fitted from a plot, polygon, or numerical pattern.

## Proof boundary

This gate proves `(Q_K-T)'>0` only on `D_20` and `Q_K-T<0` only on the
left-anchored interval (MB6).  It does not bridge the remaining gap, certify
the whole event cell, cross an event wall, prove an all-height theorem,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize conclusion.
"""


def main() -> int:
    ctx.prec = 512
    ctx.threads = 1
    require(CHECKER.is_file(), f"missing checker: {CHECKER}")
    rows = load_cache()
    certificate = assemble_certificate(rows)

    dependencies = {
        "event_atlas": {
            "path": relative(scout.EVENT_ATLAS),
            "sha256": file_hash(scout.EVENT_ATLAS),
        },
        "offset20_value_gate": {
            "path": relative(value_gate.RESULT),
            "sha256": file_hash(value_gate.RESULT),
        },
    }
    source_paths = (
        Path(__file__).resolve(),
        CHECKER,
        Path(scout.__file__).resolve(),
        Path(scout.lower_ordinary.__file__).resolve(),
        Path(scout.transition.__file__).resolve(),
        Path(scout.positive_tail.__file__).resolve(),
        Path(scout.correction.__file__).resolve(),
        Path(scout.recentered_upper_arc.__file__).resolve(),
        Path(value_gate.__file__).resolve(),
    )
    artifact = {
        "kind": "rh_rigorous_upper_event_cell_offset20_complete_K_T_monotonicity_gate",
        "date": "2026-08-28",
        "status": "offset20_complete_K_T_positive_derivative_and_left_sign_extension_certified",
        "passed": True,
        "precision_bits": 512,
        "certificate": certificate,
        "selected_append_only_component_rows": selected_rows(rows),
        "cache_path": relative(CACHE),
        "dependencies": dependencies,
        "sources": {relative(path): file_hash(path) for path in source_paths},
        "reproduction": {
            "production": (
                "python work/rh_compute/scripts/"
                f"{scout.STEM}.py --radius 0.006 --variant production"
            ),
            "independent": (
                "python work/rh_compute/scripts/"
                f"{scout.STEM}.py --radius 0.006 --variant independent"
            ),
            "checker": f"python work/rh_compute/scripts/check_{STEM}.py",
        },
        "decision": {
            "offset20_complete_K_T_derivative_positive": True,
            "production_independent_reproduction_passed": True,
            "offset20_leftward_sign_extension_proved": True,
            "connected_interval_from_first_subcell_proved": False,
            "maximal_event_cell_sign_proved": False,
            "event_wall_handoff_proved": False,
            "all_height_transport_theorem_proved": False,
            "rh_implication": False,
        },
        "proof_boundary": (
            "Only the complete derivative on the radius-0.006 offset-20 box and the resulting "
            "left-anchored negative interval are promoted. The remaining gap to the first "
            "subcell is open, and the radius-0.01 loss proves only enclosure width."
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
    extension = certificate["monotonicity_sign_extension"]
    print(
        "certified offset-20 complete K_T monotonicity bridge: "
        f"production_lower={certificate['production_derivative_component_sum']['lower']} "
        f"independent_lower={certificate['independent_derivative_component_sum']['lower']} "
        f"extended_left={extension['closed_interval'][0]}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
