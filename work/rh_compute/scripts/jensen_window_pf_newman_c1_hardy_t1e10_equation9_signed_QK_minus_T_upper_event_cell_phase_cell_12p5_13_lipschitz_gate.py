#!/usr/bin/env python3
"""Promote a two-configuration sign certificate on the 12.5--13 phase cell."""

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
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_anchored_grouped_752_phase_transport_scout as anchored
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_point_scout as point_scout


transport = anchored.transport
complete = anchored.complete
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_phase_cell_12p5_13_lipschitz_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
ANCHOR_OFFSET = "12.75"
ANCHOR = "10000000012.75"
ANCHOR_RADIUS = "0.0000001"
LEFT_OFFSET = "12.625"
LEFT_CENTER = "10000000012.625"
RIGHT_OFFSET = "12.875"
RIGHT_CENTER = "10000000012.875"
HALF_RADIUS = "0.125"
MAX_DISTANCE = "0.25"
FULL_LEFT = "10000000012.5"
FULL_RIGHT = "10000000013"
COMPONENT_PRECISION_BITS = 384
PRODUCTION_PANELS = 192
PRODUCTION_PANEL_PRECISION_BITS = 384
INDEPENDENT_PANELS = 256
INDEPENDENT_PANEL_PRECISION_BITS = 448


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def canonical_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


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


def select_panels(
    panel_cache: dict[str, dict[str, Any]],
    *,
    panels: int,
    precision_bits: int,
    variant: str,
) -> list[dict[str, Any]]:
    return transport.selected_panels(
        panel_cache, ANCHOR, panels, precision_bits, variant
    )


def select_point_row(
    point_cache: dict[str, dict[str, Any]], variant: str
) -> dict[str, Any]:
    key = point_scout.cache_key(
        ANCHOR, ANCHOR_RADIUS, COMPONENT_PRECISION_BITS, variant
    )
    require(key in point_cache, f"missing anchor point row: {key}")
    row = point_cache[key]
    require(row.get("passed") is True, f"anchor point row failed: {key}")
    require(
        row["record"]["projection"]["strictly_negative_on_box"] is True,
        f"anchor point row lost negativity: {key}",
    )
    return row


def component_manifest(
    component_cache: dict[str, dict[str, Any]], variant: str
) -> dict[str, Any]:
    selected: dict[str, Any] = {}
    for role, center in (("left", LEFT_CENTER), ("right", RIGHT_CENTER)):
        for component in complete.COMPONENTS:
            key = complete.cache_key(
                center,
                HALF_RADIUS,
                COMPONENT_PRECISION_BITS,
                variant,
                component,
            )
            require(key in component_cache, f"missing selected component row: {key}")
            row = component_cache[key]
            require(row.get("passed") is True, f"selected component row failed: {key}")
            selected[f"{role}_{component}"] = {
                "cache_key": key,
                "canonical_row_sha256": complete.canonical_hash(row),
            }
    return selected


def evaluate_half(
    panel_rows: list[dict[str, Any]],
    component_cache: dict[str, dict[str, Any]],
    *,
    target_offset: str,
    target_center: str,
    variant: str,
) -> dict[str, Any]:
    return anchored.evaluate_anchored_model(
        panel_rows,
        component_cache,
        anchor_text=ANCHOR,
        target_text=target_center,
        target_offset=target_offset,
        radius_text=HALF_RADIUS,
        precision_bits=COMPONENT_PRECISION_BITS,
        variant=variant,
    )


def model_summary(model: dict[str, Any]) -> dict[str, Any]:
    return {
        "target_offset": model["target_offset"],
        "target_center_height": model["target_center_height"],
        "radius": model["radius"],
        "closed_target_interval": model["closed_target_interval"],
        "grouped_K_transport_error": model["grouped_K_transport_error"],
        "replacement_lower_ordinary_Q_derivative": model[
            "replacement_lower_ordinary_Q_derivative"
        ],
        "replacement_complete_K_T": model["replacement_complete_K_T"],
        "Q_K_minus_T_derivative": model[
            "replacement_Q_K_minus_T_derivative_component_sum"
        ],
        "absolute_derivative_upper": model["absolute_derivative_upper"],
        "direct_component_overlap": model["replacement_direct_component_overlap"],
        "old_broad_enclosure_overlap": model["replacement_old_overlaps"],
    }


def configuration_certificate(
    *,
    variant: str,
    panel_rows: list[dict[str, Any]],
    panel_precision_bits: int,
    component_cache: dict[str, dict[str, Any]],
    point_row: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    left_model = evaluate_half(
        panel_rows,
        component_cache,
        target_offset=LEFT_OFFSET,
        target_center=LEFT_CENTER,
        variant=variant,
    )
    right_model = evaluate_half(
        panel_rows,
        component_cache,
        target_offset=RIGHT_OFFSET,
        target_center=RIGHT_CENTER,
        variant=variant,
    )
    q_anchor = real(point_row["record"]["projection"]["Q_K_minus_T_height_box"])
    left_derivative = real(
        left_model["replacement_Q_K_minus_T_derivative_component_sum"]
    )
    right_derivative = real(
        right_model["replacement_Q_K_minus_T_derivative_component_sum"]
    )
    left_modulus = abs(left_derivative).upper()
    right_modulus = abs(right_derivative).upper()
    distance = arb(MAX_DISTANCE)
    left_value_upper = (q_anchor.upper() + distance * left_modulus).upper()
    right_value_upper = (q_anchor.upper() + distance * right_modulus).upper()
    threshold = (-q_anchor.upper() / distance).lower()
    require(q_anchor.upper() < 0, f"{variant} anchor value lost negativity")
    require(left_modulus < threshold, f"{variant} left derivative modulus missed threshold")
    require(right_modulus < threshold, f"{variant} right derivative modulus missed threshold")
    require(left_value_upper < 0, f"{variant} left transported value lost negativity")
    require(right_value_upper < 0, f"{variant} right transported value lost negativity")

    certificate = {
        "variant": variant,
        "panels": len(panel_rows),
        "panel_precision_bits": panel_precision_bits,
        "panel_order": "forward" if variant == "production" else "reverse",
        "panel_manifest": panel_manifest(panel_rows),
        "component_manifest": component_manifest(component_cache, variant),
        "anchor_point_row": {
            "cache_key": point_row["cache_key"],
            "canonical_row_sha256": canonical_hash(point_row),
            "Q_K_minus_T": point_row["record"]["projection"][
                "Q_K_minus_T_height_box"
            ],
        },
        "mean_value_threshold": threshold.str(55, more=True),
        "left_half": {
            "model": model_summary(left_model),
            "derivative_modulus_upper": left_modulus.str(55, more=True),
            "transported_Q_K_minus_T_upper": left_value_upper.str(55, more=True),
            "strictly_negative_on_closed_half": True,
        },
        "right_half": {
            "model": model_summary(right_model),
            "derivative_modulus_upper": right_modulus.str(55, more=True),
            "transported_Q_K_minus_T_upper": right_value_upper.str(55, more=True),
            "strictly_negative_on_closed_half": True,
        },
        "strictly_negative_on_closed_phase_cell": True,
    }
    return certificate, {"left": left_model, "right": right_model}


def assemble_certificate() -> dict[str, Any]:
    ctx.prec = 512
    ctx.threads = 1
    panel_cache = transport.load_jsonl(transport.CACHE)
    component_cache = transport.load_jsonl(complete.CACHE)
    point_cache = point_scout.load_cache()

    production_panels = select_panels(
        panel_cache,
        panels=PRODUCTION_PANELS,
        precision_bits=PRODUCTION_PANEL_PRECISION_BITS,
        variant="production",
    )
    independent_panels = select_panels(
        panel_cache,
        panels=INDEPENDENT_PANELS,
        precision_bits=INDEPENDENT_PANEL_PRECISION_BITS,
        variant="independent",
    )
    production_point = select_point_row(point_cache, "production")
    independent_point = select_point_row(point_cache, "independent")
    production, production_models = configuration_certificate(
        variant="production",
        panel_rows=production_panels,
        panel_precision_bits=PRODUCTION_PANEL_PRECISION_BITS,
        component_cache=component_cache,
        point_row=production_point,
    )
    independent, independent_models = configuration_certificate(
        variant="independent",
        panel_rows=independent_panels,
        panel_precision_bits=INDEPENDENT_PANEL_PRECISION_BITS,
        component_cache=component_cache,
        point_row=independent_point,
    )

    production_anchor = real(production["anchor_point_row"]["Q_K_minus_T"])
    independent_anchor = real(independent["anchor_point_row"]["Q_K_minus_T"])
    require(production_anchor.overlaps(independent_anchor), "anchor value boxes miss")
    for side in ("left", "right"):
        production_derivative = real(
            production_models[side][
                "replacement_Q_K_minus_T_derivative_component_sum"
            ]
        )
        independent_derivative = real(
            independent_models[side][
                "replacement_Q_K_minus_T_derivative_component_sum"
            ]
        )
        require(
            production_derivative.overlaps(independent_derivative),
            f"{side} derivative boxes miss",
        )
        require(
            real(
                production_models[side][
                    "replacement_lower_ordinary_Q_derivative"
                ]
            ).overlaps(
                real(
                    independent_models[side][
                        "replacement_lower_ordinary_Q_derivative"
                    ]
                )
            ),
            f"{side} lower-ordinary projections miss",
        )
        require(
            complex_ball(
                production_models[side]["replacement_complete_K_T"]
            ).overlaps(
                complex_ball(
                    independent_models[side]["replacement_complete_K_T"]
                )
            ),
            f"{side} complete K_T boxes miss",
        )

    anchor = arb(ANCHOR)
    left_box = arb(arb(LEFT_CENTER), arb(HALF_RADIUS))
    right_box = arb(arb(RIGHT_CENTER), arb(HALF_RADIUS))
    require(left_box.contains(arb(FULL_LEFT)), "left phase-cell endpoint drift")
    require(left_box.contains(anchor), "left half misses anchor")
    require(right_box.contains(anchor), "right half misses anchor")
    require(right_box.contains(arb(FULL_RIGHT)), "right phase-cell endpoint drift")
    event = load_json(complete.EVENT_ATLAS)["certificate"]["primary_open_cell"]
    event_lower = arb(event["lower_ball"]["ball"])
    event_upper = arb(event["upper_ball"]["ball"])
    require(event_lower.upper() < left_box.lower(), "phase cell crosses lower event wall")
    require(right_box.upper() < event_upper.lower(), "phase cell crosses upper event wall")

    return {
        "anchor_height": ANCHOR,
        "anchor_probe_radius": ANCHOR_RADIUS,
        "closed_phase_cell": [FULL_LEFT, FULL_RIGHT],
        "left_closed_half": [FULL_LEFT, ANCHOR],
        "right_closed_half": [ANCHOR, FULL_RIGHT],
        "maximum_anchor_distance": MAX_DISTANCE,
        "strictly_inside_fixed_roster_event_cell": True,
        "event_cell_open_walls": {
            "lower": event_lower.str(65, more=True),
            "upper": event_upper.str(65, more=True),
        },
        "panel_transport_identity": {
            "height_transport": "F_t(y)=F_c(y)*exp(-i*(t-c)*log(y))",
            "common_phase": "exp(-i*(t-c)*s0)",
            "panel_error": "|t-c|*epsilon_j*B_j",
            "log_panel_error": "log(C)*|t-c|*epsilon_j*B_j",
            "mean_value_phase": "theta(t)=theta(c)+(t-c)*theta'(xi_t)",
        },
        "value_transport_identity": {
            "identity": "Q(t)=Q(a)+integral_a^t Q'(u) du",
            "bound": "Q(t)<=Q(a)+|t-a| sup_I |Q'|",
            "anchor": ANCHOR,
            "maximum_distance": MAX_DISTANCE,
            "no_stationary_point_uniqueness_required": True,
        },
        "production": production,
        "independent": independent,
        "production_independent_anchor_overlap": True,
        "production_independent_left_derivative_overlap": True,
        "production_independent_right_derivative_overlap": True,
        "production_independent_lower_ordinary_overlap": True,
        "production_independent_complete_K_T_overlap": True,
        "Q_K_minus_T_strictly_negative_on_closed_phase_cell": True,
        "scope_boundary": (
            "This is one closed phase cell inside one fixed-roster event cell. It neither "
            "covers the remaining height gap nor crosses an event wall, and it does not imply RH."
        ),
    }


def render_note(certificate: dict[str, Any]) -> str:
    production = certificate["production"]
    independent = certificate["independent"]
    return f"""# Phase-cell 12.5--13 anchored Lipschitz gate

Date: 2026-08-28

Status: rigorous two-configuration closed-cell sign certificate; not an
all-height theorem and not a proof of RH

Put `a=10^10+12.75` and split

```text
I_L=[10^10+12.5, 10^10+12.75],
I_R=[10^10+12.75, 10^10+13].                    (PC1)
```

Both intervals lie strictly inside the same fixed-roster event cell.  The
joined packet therefore has no roster jump on either half.

## Reused grouped atlas

For the grouped 752-label collar, a panel atlas is evaluated once at `a`.  If
`h=t-a`, the exact height transport is

```text
F_t(y)=F_a(y) exp(-i h log y).                  (PC2)
```

Factoring the common rotation at log-center `s_0` and one panel rotation at
`s_j` leaves the rigorous errors

```text
|E_j| <= |h| epsilon_j B_j,
|E_j^log| <= log(C)|h| epsilon_j B_j.           (PC3)
```

Here `epsilon_j` is the panel log half-width and
`B_j=752(y_right-y_left)/sqrt(y_left)`.  Production uses 192 forward panels
at 384 bits; the independent path uses 256 reverse panels at 448 bits.  No new
occurrence of pi is introduced in this phase-cell argument.

## Anchor and derivative norms

The production anchor upper endpoint is

```text
Q(a) <= {production['anchor_point_row']['Q_K_minus_T']['upper']}. (PC4)
```

The independent anchor upper endpoint is

```text
Q(a) <= {independent['anchor_point_row']['Q_K_minus_T']['upper']}. (PC5)
```

The transported complete derivative enclosures give

```text
production:  sup_I_L |Q'| <= {production['left_half']['derivative_modulus_upper']}
             sup_I_R |Q'| <= {production['right_half']['derivative_modulus_upper']}
independent: sup_I_L |Q'| <= {independent['left_half']['derivative_modulus_upper']}
             sup_I_R |Q'| <= {independent['right_half']['derivative_modulus_upper']}. (PC6)
```

For every `t` in either half, `|t-a|<=0.25`, so the real mean-value theorem
gives

```text
Q(t) = Q(a) + integral_a^t Q'(u) du
     <= Q(a) + 0.25 sup_I |Q'|.                 (PC7)
```

The resulting worst-case upper endpoints are

```text
production left:  {production['left_half']['transported_Q_K_minus_T_upper']}
production right: {production['right_half']['transported_Q_K_minus_T_upper']}
independent left: {independent['left_half']['transported_Q_K_minus_T_upper']}
independent right:{independent['right_half']['transported_Q_K_minus_T_upper']}. (PC8)
```

All four are strictly negative.  Since `I_L union I_R` is the full closed
interval in (PC1),

```text
Q_K(t)-T(t) < 0
for every t in [10^10+12.5, 10^10+13].          (PC9)
```

This proof does not assume or prove uniqueness of the stationary point seen
near offset 12.76.  Any additional hidden critical points are absorbed by the
rigorous derivative norm.  The certificate covers only this one half-unit
phase cell and makes no claim about the rest of the upper event cell.
"""


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def main() -> int:
    require(CHECKER.is_file(), f"missing checker: {CHECKER}")
    certificate = assemble_certificate()
    note = render_note(certificate)
    atomic_write(NOTE, note)
    source_paths = (
        Path(__file__).resolve(),
        CHECKER.resolve(),
        Path(anchored.__file__).resolve(),
        Path(transport.__file__).resolve(),
        Path(complete.__file__).resolve(),
        Path(point_scout.__file__).resolve(),
    )
    artifact = {
        "kind": "rh_c1_hardy_upper_event_cell_phase_cell_lipschitz_gate",
        "date": "2026-08-28",
        "status": "rigorous_two_configuration_closed_phase_cell_sign_certificate",
        "passed": True,
        "certificate": certificate,
        "proof_note": relative(NOTE),
        "proof_note_sha256": file_hash(NOTE),
        "checker": relative(CHECKER),
        "panel_cache": relative(transport.CACHE),
        "component_cache": relative(complete.CACHE),
        "point_cache": relative(point_scout.CACHE),
        "dependency_artifacts": {
            relative(complete.EVENT_ATLAS): file_hash(complete.EVENT_ATLAS),
        },
        "sources": {relative(path): file_hash(path) for path in source_paths},
    }
    atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        "promoted phase-cell sign gate: "
        f"interval=[{FULL_LEFT},{FULL_RIGHT}] "
        f"production_upper={certificate['production']['right_half']['transported_Q_K_minus_T_upper']} "
        f"independent_upper={certificate['independent']['right_half']['transported_Q_K_minus_T_upper']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
