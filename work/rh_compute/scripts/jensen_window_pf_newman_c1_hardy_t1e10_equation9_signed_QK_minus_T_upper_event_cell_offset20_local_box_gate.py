#!/usr/bin/env python3
"""Promote the independently reproduced offset-20 local sign box."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_local_radius_ladder_scout as ladder


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_offset20_local_box_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
CACHE = ladder.CACHE
CENTER = "10000000020"
OFFSET = "20"
PASS_RADIUS = "0.00012"
LOSS_RADIUS = "0.00013"
PRECISION_BITS = 384
VERSION = ladder.VERSION
EVENT_ATLAS = ladder.point.EVENT_ATLAS
UPPER_ARC_SOURCE = Path(ladder.point.base.upper_arc.__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def row_hash(row: dict[str, Any]) -> str:
    encoded = json.dumps(row, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_cache(path: Path = CACHE) -> dict[str, dict[str, Any]]:
    require(path.is_file(), f"missing ladder cache: {path}")
    rows: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            key = row.get("cache_key")
            require(key, f"cache row {line_number} has no key")
            require(key not in rows, f"duplicate cache key at row {line_number}")
            rows[key] = row
    return rows


def key(radius: str, variant: str) -> str:
    return ladder.cache_key(CENTER, radius, PRECISION_BITS, variant)


def complex_from(record: dict[str, Any]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def real_from(record: dict[str, Any]) -> arb:
    """Reconstruct from directed endpoints without repeated display widening."""

    lower = arb(record["lower"]).lower()
    upper = arb(record["upper"]).upper()
    midpoint = (lower + upper) / 2
    radius = ((upper - lower) / 2).upper()
    return arb(midpoint, radius)


def real_record(value: arb, digits: int = 70) -> dict[str, str]:
    return {
        "ball": value.str(digits, more=True),
        "lower": value.lower().str(55, more=True),
        "upper": value.upper().str(55, more=True),
        "radius": value.rad().str(40, more=True),
    }


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(55, more=True),
    }


def semantic_arc_audit() -> dict[str, Any]:
    source = UPPER_ARC_SOURCE.read_text(encoding="utf-8")
    action_alarm = 'require(action_max < arb("0.01981"), "production action maximum exceeded its guard")'
    scale_alarm = 'require(rotated_arc_absolute < arb("0.057665"), "upper arc missed its certified scale guard")'
    require(source.count('arb("0.01981")') == 1, "0.01981 acquired another source use")
    require(source.count('arb("0.057665")') == 1, "0.057665 acquired another source use")
    require(action_alarm in source, "action calibration alarm changed")
    require(scale_alarm in source, "arc scale calibration alarm changed")
    require(source.count("action_max.exp()") == 2, "actual action is no longer used twice in errors")
    require(
        source.index("action_max = action(delta_star)")
        < source.index("polynomial_integral_error")
        < source.index('require(total_added_error < arb("1e-20")'),
        "action/error evaluation order changed",
    )
    require(
        source.index("rotated_arc_absolute = abs(rotated_arc)")
        < source.index(scale_alarm)
        < source.index("actual_arc ="),
        "arc scale alarm moved into the enclosure construction",
    )
    return {
        "action_alarm_literal_occurrences": 1,
        "arc_scale_alarm_literal_occurrences": 1,
        "actual_interval_action_used_in_error_bounds": True,
        "total_added_error_guard_retained": "total_added_error<1e-20",
        "arc_scale_alarm_not_substituted_into_actual_arc": True,
        "upper_arc_source": relative(UPPER_ARC_SOURCE),
        "upper_arc_source_sha256": file_hash(UPPER_ARC_SOURCE),
    }


def assemble_certificate(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    production = rows[key(PASS_RADIUS, "production")]
    independent = rows[key(PASS_RADIUS, "independent")]
    loss = rows[key(LOSS_RADIUS, "production")]

    for name, row in (("production", production), ("independent", independent)):
        require(row["analytic_guards_passed"] is True, f"{name} analytic guards failed")
        require(row["strictly_negative"] is True, f"{name} pass box is not negative")
        require(row["center_height"] == CENTER, f"{name} center drift")
        require(row["radius"] == PASS_RADIUS, f"{name} radius drift")
        require(row["precision_bits"] == PRECISION_BITS, f"{name} precision drift")
        require(row["version"] == VERSION, f"{name} ladder version drift")
        require(
            row["upper_arc_action_guard"]["all_substantive_donor_error_bounds_retained"] is True,
            f"{name} upper-arc error guard set was weakened",
        )

    q_production = real_from(production["Q_K_minus_T"])
    q_independent = real_from(independent["Q_K_minus_T"])
    packet_production = complex_from(production["complete_packet_width"])
    packet_independent = complex_from(independent["complete_packet_width"])
    transition_production = complex_from(production["joined_transition_width"])
    transition_independent = complex_from(independent["joined_transition_width"])
    require(q_production.upper() < 0, "production Q_K-T box lost negativity")
    require(q_independent.upper() < 0, "independent Q_K-T box lost negativity")
    require(q_production.overlaps(q_independent), "two Q_K-T boxes do not overlap")
    require(packet_production.overlaps(packet_independent), "complete packet boxes do not overlap")
    require(transition_production.overlaps(transition_independent), "transition boxes do not overlap")

    require(loss["analytic_guards_passed"] is True, "loss row is not a clean width diagnostic")
    require(loss["strictly_negative"] is False, "loss row unexpectedly remains negative")
    q_loss = real_from(loss["Q_K_minus_T"])
    require(q_loss.lower() < 0 < q_loss.upper(), "loss row does not straddle zero")

    event = load_json(EVENT_ATLAS)["certificate"]["primary_open_cell"]
    event_lower = arb(event["lower_ball"]["ball"])
    event_upper = arb(event["upper_ball"]["ball"])
    t_box = arb(arb(CENTER), arb(PASS_RADIUS))
    require(event_lower.upper() < t_box.lower(), "certified box crosses lower event wall")
    require(t_box.upper() < event_upper.lower(), "certified box crosses upper event wall")

    return {
        "offset_from_t0": OFFSET,
        "center_height": CENTER,
        "radius": PASS_RADIUS,
        "closed_height_interval": [
            t_box.lower().str(55, more=True),
            t_box.upper().str(55, more=True),
        ],
        "strictly_inside_fixed_roster_event_cell": True,
        "event_cell_open_walls": {
            "lower": event_lower.str(65, more=True),
            "upper": event_upper.str(65, more=True),
        },
        "production_Q_K_minus_T": real_record(q_production),
        "independent_Q_K_minus_T": real_record(q_independent),
        "production_independent_Q_K_minus_T_overlap": True,
        "production_complete_packet": complex_record(packet_production),
        "independent_complete_packet": complex_record(packet_independent),
        "production_independent_complete_packet_overlap": True,
        "production_joined_transition": complex_record(transition_production),
        "independent_joined_transition": complex_record(transition_independent),
        "production_independent_joined_transition_overlap": True,
        "direct_radius_bracket": {
            "largest_tested_certified_radius": PASS_RADIUS,
            "first_tested_zero_straddling_radius": LOSS_RADIUS,
            "loss_Q_K_minus_T": real_record(q_loss),
            "loss_kind": "interval_enclosure_width_only",
            "pointwise_nonnegativity_inferred_from_loss": False,
        },
        "upper_arc_calibration_audit": semantic_arc_audit(),
        "identities": {
            "height_coordinate": "p=t/(2*pi)",
            "Fourier_kernel": "exp(2*pi*i*q*y)",
            "quadratic_kernel": "exp(-i*pi*y^2)",
            "assembled_quantity": "Q_K(t)-T(t)=Hardy_t[P_W(t)]/H(t)",
        },
    }


def selected_rows(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    selected = {
        "production_pass": rows[key(PASS_RADIUS, "production")],
        "independent_pass": rows[key(PASS_RADIUS, "independent")],
        "production_first_loss": rows[key(LOSS_RADIUS, "production")],
    }
    return {
        name: {
            "cache_key": row["cache_key"],
            "canonical_row_sha256": row_hash(row),
        }
        for name, row in selected.items()
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    p = c["production_Q_K_minus_T"]
    i = c["independent_Q_K_minus_T"]
    loss = c["direct_radius_bracket"]
    interval = c["closed_height_interval"]
    return f"""# Offset-20 local `Q_K-T` sign box

Date: 2026-08-28

Status: rigorous two-configuration local interval certificate; disconnected
from the first height subcell and not a maximal-event-cell or RH theorem

Let

```text
I_20=[{interval[0]},
      {interval[1]}]
    =[10^10+20-0.00012, 10^10+20+0.00012].           (LB1)
```

The event-atlas wall enclosures place all of `I_20` strictly inside the same
fixed-roster open cell as `t=10^10`.  The complete joined contour packet is
therefore evaluated without changing labels or crossing an event wall.

The production packet gives

```text
Q_K-T in {p['ball']},
upper endpoint = {p['upper']} < 0.                    (LB2)
```

An independent configuration changes transition panels, endpoint-remainder
slabs, ordinary expansion order, and transition-jet order.  It gives

```text
Q_K-T in {i['ball']},
upper endpoint = {i['upper']} < 0.                    (LB3)
```

The two `Q_K-T`, complete-packet, and joined-transition boxes overlap.  Thus

```text
Q_K(t)-T(t)<0 for every t in I_20.                    (LB4)
```

## Radius barrier

At the same center, every analytic guard also passes at radius
`{loss['first_tested_zero_straddling_radius']}`, but the direct Arb enclosure is

```text
{loss['loss_Q_K_minus_T']['ball']}.                   (LB5)
```

It straddles zero.  This is an interval-dependency/width obstruction only; it
does not imply a pointwise sign change.

## Where pi comes from

Here `pi` is the ordinary circle constant, circumference divided by diameter.
It is not fitted or inserted from a plot.  It enters the exact Fourier kernel
`exp(2*pi*i*q*y)` and quadratic contour kernel `exp(-i*pi*y^2)`.  Differentiating
their phases produces the natural height coordinate `p=t/(2*pi)` used in the
stationary-point equations.

## Proof boundary

`I_20` is a second certified local island.  There remains an uncovered gap of
almost 20 height units between the first box around `10^10` and this box.  No
derivative sign on that gap, connected interval from `I_1`, maximal-event-cell
coverage, event-wall handoff, all-height theorem, RH, or Clay-prize conclusion
is asserted.
"""


def main() -> int:
    ctx.prec = 512
    ctx.threads = 1
    rows = load_cache()
    certificate = assemble_certificate(rows)
    require(CHECKER.is_file(), f"missing checker: {CHECKER}")

    dependencies = {
        "event_atlas": {
            "path": relative(EVENT_ATLAS),
            "sha256": file_hash(EVENT_ATLAS),
        }
    }
    source_paths = (
        Path(__file__).resolve(),
        CHECKER,
        Path(ladder.__file__).resolve(),
        Path(ladder.point.__file__).resolve(),
        Path(ladder.base.__file__).resolve(),
        UPPER_ARC_SOURCE,
    )
    artifact = {
        "kind": "rh_rigorous_upper_event_cell_offset20_local_box_gate",
        "date": "2026-08-28",
        "status": "offset20_closed_local_Q_K_minus_T_negative_interval_certified",
        "passed": True,
        "precision_bits": 512,
        "certificate": certificate,
        "selected_append_only_cache_rows": selected_rows(rows),
        "cache_path": relative(CACHE),
        "dependencies": dependencies,
        "sources": {relative(path): file_hash(path) for path in source_paths},
        "reproduction": {
            "production": (
                "python work/rh_compute/scripts/"
                f"{ladder.STEM}.py --offset 20 --radius 0.00012 --variant production --max-new 1"
            ),
            "independent": (
                "python work/rh_compute/scripts/"
                f"{ladder.STEM}.py --offset 20 --radius 0.00012 --variant independent --max-new 1"
            ),
            "checker": f"python work/rh_compute/scripts/check_{STEM}.py",
        },
        "decision": {
            "offset20_local_closed_box_sign_proved": True,
            "production_independent_reproduction_passed": True,
            "direct_radius_width_bracket_measured": True,
            "connected_interval_from_first_subcell_proved": False,
            "maximal_event_cell_sign_proved": False,
            "event_wall_handoff_proved": False,
            "all_height_transport_theorem_proved": False,
            "rh_implication": False,
        },
        "proof_boundary": (
            "Only Q_K-T<0 on the single closed offset-20 box is promoted. The radius-0.00013 "
            "row diagnoses direct-enclosure width and proves no pointwise sign failure."
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
        "certified offset-20 local box: "
        f"production_upper={certificate['production_Q_K_minus_T']['upper']} "
        f"independent_upper={certificate['independent_Q_K_minus_T']['upper']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
