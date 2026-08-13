#!/usr/bin/env python3
"""Certify a first nonzero fixed-selector height cell for the lower fold."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_airy_fresnel_detuning_ode_gate as ode_gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_height_transport_identity_gate as height_gate


HEIGHT_GATE = height_gate.RESULT
HEIGHT_CACHE = height_gate.CACHE
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_nonzero_height_cell_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_nonzero_height_cell_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = ode_gate.PRECISION
HEIGHT_RADIUS_TEXT = "0.01"
AI_ENVELOPE_TEXT = "0.54"
SECOND_DERIVATIVE_BOUND_TEXT = "0.146"
TAYLOR_REMAINDER_BOUND_TEXT = "0.0000073"
TOTAL_VARIATION_BOUND_TEXT = "0.001313"
BUILDER_PANELS = 2048


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def certify_airy_envelope(
    argument_left: arb,
    argument_right: arb,
    panels: int,
    envelope: arb,
) -> dict[str, Any]:
    """Prove the envelope on a covering interval partition."""
    require(argument_left < argument_right, "Airy argument interval is reversed")
    width = (argument_right - argument_left) / panels
    worst_upper = arb(0)
    worst_panel = -1
    for panel in range(panels):
        left = argument_left + panel * width
        right = argument_left + (panel + 1) * width
        cell = (left + right) / 2 + arb(0, (right - left) / 2)
        magnitude = abs(cell.airy_ai())
        require(magnitude < envelope, f"Airy envelope failed on panel {panel}")
        upper = magnitude.upper()
        if float(upper) > float(worst_upper):
            worst_upper = upper
            worst_panel = panel
    return {
        "panels": panels,
        "argument_left_ball": argument_left.str(PRECISION, more=True),
        "argument_right_ball": argument_right.str(PRECISION, more=True),
        "envelope_constant": envelope.str(PRECISION, more=True),
        "largest_observed_interval_upper": worst_upper.str(PRECISION, more=True),
        "largest_observed_panel": worst_panel,
        "all_panels_strictly_below_envelope": True,
    }


def grouped_values(rows: dict[int, dict[str, Any]], p: dict[str, arb]) -> tuple[acb, acb, acb]:
    i = acb(0, 1)
    beta = p["beta"]
    lam = p["lambda"]
    values: list[acb] = []
    first_derivatives: list[acb] = []
    second_derivatives: list[acb] = []
    for mode in range(ode_gate.MODE_LO, ode_gate.MODE_HI + 1):
        row = rows[mode]
        g = height_gate.parse_complex(row["G64"])
        g_prime = height_gate.parse_complex(row["G64_prime_detuning"])
        g_t = height_gate.parse_complex(row["G64_height_derivative_fixed_C"])
        values.append(g)
        first_derivatives.append(g_t)
        second_derivatives.append((-lam * g - i * g_prime) / beta**2)
    factor = 2 * arb.pi()
    return (
        factor * sum(values, acb(0)),
        factor * sum(first_derivatives, acb(0)),
        factor * sum(second_derivatives, acb(0)),
    )


def certificate(
    rows: dict[int, dict[str, Any]],
    p: dict[str, arb],
    panels: int,
    *,
    radius_text: str = HEIGHT_RADIUS_TEXT,
    second_derivative_bound_text: str = SECOND_DERIVATIVE_BOUND_TEXT,
    taylor_remainder_bound_text: str = TAYLOR_REMAINDER_BOUND_TEXT,
    total_variation_bound_text: str = TOTAL_VARIATION_BOUND_TEXT,
) -> dict[str, Any]:
    expected_modes = set(range(ode_gate.MODE_LO, ode_gate.MODE_HI + 1))
    require(set(rows) == expected_modes, "height-transport cache is incomplete")

    radius = arb(radius_text)
    envelope = arb(AI_ENVELOPE_TEXT)
    beta = p["beta"]
    lam = p["lambda"]
    lambda_min = lam - radius / beta
    lambda_max = lam + radius / beta
    require(lambda_min > 0, "lambda changed sign on the height cell")

    argument_left = -(lambda_max + ode_gate.Y)
    argument_right = -lambda_min
    airy = certify_airy_envelope(argument_left, argument_right, panels, envelope)

    mode_count = ode_gate.MODE_HI - ode_gate.MODE_LO + 1
    integral_weight_bound = ode_gate.Y * lambda_max + arb(ode_gate.Y) ** 2 / 2
    derived_second_bound = (
        2 * arb.pi() * mode_count * envelope * integral_weight_bound / beta**2
    )
    stated_second_bound = arb(second_derivative_bound_text)
    require(derived_second_bound < stated_second_bound, "stated second-derivative bound is too small")

    grouped, grouped_t, grouped_tt = grouped_values(rows, p)
    saved_transport = json.loads(HEIGHT_GATE.read_text(encoding="utf-8"))["certified_transport"]
    require(
        grouped.overlaps(height_gate.parse_complex(saved_transport["grouped_2pi_G64"])),
        "grouped value does not overlap the transport artifact",
    )
    require(
        grouped_t.overlaps(height_gate.parse_complex(saved_transport["grouped_2pi_height_derivative_fixed_C"])),
        "grouped first derivative does not overlap the transport artifact",
    )

    derived_remainder_bound = derived_second_bound * radius**2 / 2
    stated_remainder_bound = arb(taylor_remainder_bound_text)
    require(derived_remainder_bound <= stated_remainder_bound, "stated Taylor remainder is too small")
    linear_excursion = abs(grouped_t).upper() * radius
    derived_total_variation = linear_excursion + stated_remainder_bound
    stated_total_variation = arb(total_variation_bound_text)
    require(derived_total_variation < stated_total_variation, "stated total variation is too small")

    selector_upper = arb.pi() * arb(ode_gate.C) ** 2 / 8
    selector_lower = arb.pi() * arb(ode_gate.C - 2) ** 2 / 8
    height_lower = arb(ode_gate.T) - radius
    height_upper = arb(ode_gate.T) + radius
    require(height_lower > selector_lower and height_upper < selector_upper, "height cell crosses a selector fold")

    return {
        "height_center": ode_gate.T,
        "height_radius": radius_text,
        "height_lower_ball": height_lower.str(PRECISION, more=True),
        "height_upper_ball": height_upper.str(PRECISION, more=True),
        "fixed_odd_selector_C": ode_gate.C,
        "selector_cell_lower_ball": selector_lower.str(PRECISION, more=True),
        "selector_cell_upper_ball": selector_upper.str(PRECISION, more=True),
        "lambda_center_ball": lam.str(PRECISION, more=True),
        "lambda_min_ball": lambda_min.str(PRECISION, more=True),
        "lambda_max_ball": lambda_max.str(PRECISION, more=True),
        "airy_envelope": airy,
        "grouped_center_value": complex_record(grouped),
        "grouped_center_first_height_derivative": complex_record(grouped_t),
        "grouped_center_second_height_derivative": complex_record(grouped_tt),
        "derived_uniform_second_derivative_bound_ball": derived_second_bound.str(PRECISION, more=True),
        "stated_uniform_second_derivative_bound": second_derivative_bound_text,
        "derived_taylor_remainder_bound_ball": derived_remainder_bound.str(PRECISION, more=True),
        "stated_taylor_remainder_bound": taylor_remainder_bound_text,
        "linear_excursion_bound_ball": linear_excursion.str(PRECISION, more=True),
        "derived_total_variation_bound_ball": derived_total_variation.str(PRECISION, more=True),
        "stated_total_variation_bound": total_variation_bound_text,
        "uniform_cell_theorem": (
            f"For every real t with |t-10^10|<={radius_text} and fixed selector C=159577, "
            f"|S(t)-S(t0)-(t-t0)S'(t0)|<{taylor_remainder_bound_text} and "
            f"|S(t)-S(t0)|<{total_variation_bound_text}."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_height_cell"]
    return f"""# First nonzero fixed-selector lower-fold height cell

Date: 2026-08-11

Status: rigorous local height-cell enclosure validated; not a proof of a
selector-transition join or a complete source-height theorem

Put

```text
S(t)=2pi sum_m G_64(d_m;lambda(t)),
lambda(t)=[pi*C^2/8-t]/beta,
C=159577.
```

On the selector-stable cell

```text
|t-10^10| <= {c['height_radius']},
{c['height_lower_ball']} <= t <= {c['height_upper_ball']},
```

the Airy arguments lie in

```text
{c['airy_envelope']['argument_left_ball']}
    <= -lambda(t)-y <=
{c['airy_envelope']['argument_right_ball']},
0 <= y <= 64.
```

An Arb interval subdivision into {c['airy_envelope']['panels']} covering panels
proves `|Ai(x)|<0.54` throughout that interval.  Since

```text
G_tt=(-lambda*G-i*G')/beta^2,
```

the integral representation gives the uniform grouped bound

```text
|S''(t)| < {c['stated_uniform_second_derivative_bound']}.
```

Taylor's theorem therefore proves, for every real height in the cell,

```text
|S(t)-S(t0)-(t-t0)S'(t0)| < {c['stated_taylor_remainder_bound']},
|S(t)-S(t0)|                 < {c['stated_total_variation_bound']}.
```

This is the first certified nonzero-radius height theorem for the canonical
lower-fold lattice.  The radius is intentionally small: it establishes the
local mechanism and a reusable checker before any large subdivision campaign.

Proof boundary: this certificate covers only `t=10^10 +/- 0.01` with the odd
selector fixed.  It does not cross a selector fold, join adjacent rosters,
cover the source interval, prove `T_upper`, prove `Lambda<=0`, prove RH, or
establish a prize-level result.
"""


def main() -> None:
    started = time.perf_counter()
    priority = ode_gate.set_low_priority()
    for dependency in (HEIGHT_GATE, HEIGHT_CACHE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    ctx.dps = PRECISION
    p = ode_gate.parameters()
    rows = height_gate.load_cache()
    certified = certificate(rows, p, BUILDER_PANELS)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_nonzero_height_cell_gate",
        "status": "first_nonzero_fixed_selector_lower_fold_height_cell_complete",
        "passed": True,
        "exact_identity": {
            "second_height_derivative": "G_tt=(-lambda*G-i*G')/beta^2",
            "derivation": (
                "Because d lambda/dt=-1/beta is constant and "
                "partial_lambda^2 Ai(-lambda-y)=(-lambda-y)Ai(-lambda-y)."
            ),
        },
        "certified_height_cell": certified,
        "decision": {
            "selector_fixed_throughout_cell": True,
            "airy_envelope_interval_certified": True,
            "uniform_second_height_derivative_bound_proved": True,
            "nonzero_radius_height_cell_proved": True,
            "selector_transition_join_proved": False,
            "source_height_interval_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "height_transport_gate": {"path": relative(HEIGHT_GATE), "sha256": file_hash(HEIGHT_GATE)},
            "height_transport_cache": {"path": relative(HEIGHT_CACHE), "sha256": file_hash(HEIGHT_CACHE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_obligation": (
            "Replace the deliberately local radius by a resumable selector-cell subdivision, "
            "then prove adjacent-roster endpoint reassembly at every crossed odd selector fold."
        ),
        "proof_boundary": (
            "A first nonzero fixed-selector height cell only. No selector transition join, "
            "source-height theorem, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built first nonzero lower-fold height cell: |S''|<0.146 and radius=0.01 certified")


if __name__ == "__main__":
    main()
