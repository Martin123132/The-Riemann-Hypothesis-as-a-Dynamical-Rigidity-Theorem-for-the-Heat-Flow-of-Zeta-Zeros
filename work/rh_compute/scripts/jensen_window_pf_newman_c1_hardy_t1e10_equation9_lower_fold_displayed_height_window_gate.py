#!/usr/bin/env python3
"""Cover the source evaluator's complete displayed lower-fold height window."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import re
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

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_nonzero_height_cell_gate as cell_gate


SOURCE_OUTPUT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/t1e10/zeta14cubicmult_et005_output.txt"
PORTCULLIS_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json"
CELL_GATE = cell_gate.RESULT
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_displayed_height_window_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_displayed_height_window_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = cell_gate.PRECISION
HEIGHT_RADIUS_TEXT = "0.07"
SECOND_DERIVATIVE_BOUND_TEXT = "0.146"
TAYLOR_REMAINDER_BOUND_TEXT = "0.000358"
TOTAL_VARIATION_BOUND_TEXT = "0.009494"
BUILDER_PANELS = 2560


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_height_roster() -> dict[str, Any]:
    text = SOURCE_OUTPUT.read_text(encoding="utf-8")
    center_match = re.search(r"^t=\s*([0-9.]+e[+-]?\d+)\s*$", text, re.MULTILINE | re.IGNORECASE)
    require(center_match is not None, "source central height was not found")
    center = Decimal(center_match.group(1))
    offsets = [
        Decimal(match.group(1))
        for match in re.finditer(
            r"^Grand total of Hardy function Z\(t([+-]\d+\.\d+)\)=",
            text,
            re.MULTILINE,
        )
    ]
    expected = [Decimal("-0.07") + Decimal(index) * Decimal("0.01") for index in range(15)]
    require(offsets == expected, "displayed source-height offsets are not the expected 15-point roster")
    require(center == Decimal("1e10"), "source central height is not 10^10")
    return {
        "central_height": str(center),
        "offset_count": len(offsets),
        "offsets": [f"{offset:+.2f}" for offset in offsets],
        "minimum_offset": f"{offsets[0]:+.2f}",
        "maximum_offset": f"{offsets[-1]:+.2f}",
        "step": "0.01",
        "printed_hardy_values_used_as_proof_data": False,
    }


def turning_roster_stability() -> dict[str, Any]:
    radius = arb(HEIGHT_RADIUS_TEXT)
    center = arb(cell_gate.ode_gate.T)
    endpoint = arb(cell_gate.ode_gate.C)
    pi = arb.pi()

    def roots(height: arb) -> tuple[arb, arb]:
        discriminant = (1 - 8 * height / (pi * endpoint**2)).sqrt()
        return endpoint * (1 - discriminant) / 4, endpoint * (1 + discriminant) / 4

    height_low = center - radius
    height_high = center + radius
    lower_at_low, upper_at_low = roots(height_low)
    lower_at_high, upper_at_high = roots(height_high)
    nt_low = (height_low / (2 * pi)).sqrt()
    nt_high = (height_high / (2 * pi)).sqrt()

    require(arb(39852) < lower_at_low and lower_at_high < arb(39853), "lower turning root changes integer cell")
    require(arb(39936) < upper_at_high and upper_at_low < arb(39937), "upper turning root changes integer cell")
    require(arb(39894) < nt_low and nt_high < arb(39895), "Riemann-Siegel split changes integer cell")
    margins = [
        lower_at_low - 39852,
        39853 - lower_at_high,
        upper_at_high - 39936,
        39937 - upper_at_low,
        nt_low - 39894,
        39895 - nt_high,
    ]
    require(all(margin > 0 for margin in margins), "turning-roster margin is not positive")
    minimum_margin = min(margins, key=float)
    return {
        "lower_root_at_height_low_ball": lower_at_low.str(PRECISION, more=True),
        "lower_root_at_height_high_ball": lower_at_high.str(PRECISION, more=True),
        "upper_root_at_height_high_ball": upper_at_high.str(PRECISION, more=True),
        "upper_root_at_height_low_ball": upper_at_low.str(PRECISION, more=True),
        "riemann_siegel_split_at_height_low_ball": nt_low.str(PRECISION, more=True),
        "riemann_siegel_split_at_height_high_ball": nt_high.str(PRECISION, more=True),
        "minimum_integer_boundary_margin_ball": minimum_margin.str(PRECISION, more=True),
        "lower_turning_modes": "39853..39894",
        "upper_turning_modes": "39895..39936",
        "complete_turning_roster": "39853..39936",
        "mode_count": 84,
        "roster_constant_throughout_window": True,
    }


def certificate(panels: int) -> dict[str, Any]:
    roster = source_height_roster()
    turning_roster = turning_roster_stability()
    p = cell_gate.ode_gate.parameters()
    rows = cell_gate.height_gate.load_cache()
    cell = cell_gate.certificate(
        rows,
        p,
        panels,
        radius_text=HEIGHT_RADIUS_TEXT,
        second_derivative_bound_text=SECOND_DERIVATIVE_BOUND_TEXT,
        taylor_remainder_bound_text=TAYLOR_REMAINDER_BOUND_TEXT,
        total_variation_bound_text=TOTAL_VARIATION_BOUND_TEXT,
    )
    require(arb(cell["height_radius"]) == arb(HEIGHT_RADIUS_TEXT) or cell["height_radius"] == HEIGHT_RADIUS_TEXT, "height radius mismatch")
    return {
        "source_height_roster": roster,
        "turning_roster_stability": turning_roster,
        "canonical_lower_fold_window": cell,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_displayed_window"]["canonical_lower_fold_window"]
    r = artifact["certified_displayed_window"]["source_height_roster"]
    tr = artifact["certified_displayed_window"]["turning_roster_stability"]
    return f"""# Complete displayed lower-fold height window

Date: 2026-08-11

Status: rigorous canonical lower-fold enclosure on all displayed source
heights; not a proof of the complete source formula or a selector-fold join

The saved evaluator output requests the exact 15-point roster

```text
t-0.07, t-0.06, ..., t, ..., t+0.06, t+0.07,
t=10^10.                                                (DW1)
```

The roster is parsed from the pinned source output.  Its printed Hardy values
are not used as proof data.

Across the same height cell, the exact inverse roots remain in

```text
39852 < N_-(C;t) < 39853,
39936 < N_+(C;t) < 39937,
39894 < sqrt(t/(2pi)) < 39895,                         (DW2)
```

with minimum integer-boundary margin
`{tr['minimum_integer_boundary_margin_ball']}`.  Hence the source turning
roster stays exactly `39853..39936` with 84 modes throughout the window.

Throughout the one selector-stable interval `|t-10^10|<=0.07`, the exact
fixed-selector identity remains

```text
G_tt=(-lambda G-i G')/beta^2.                          (DW3)
```

An Arb cover with {c['airy_envelope']['panels']} intervals proves
`|Ai(x)|<0.54` on the complete Airy argument range.  A differently partitioned
independent checker repeats the cover.  Therefore

```text
|S''(t)|<0.146,                                         (DW4)
|S(t)-S(t0)-(t-t0)S'(t0)|<0.000358,                   (DW5)
|S(t)-S(t0)|<0.009494.                                 (DW6)
```

All {r['offset_count']} displayed heights lie in this single certified cell.
This proves height variation only for the canonical grouped lower-fold
transform `S`; it does not certify the other source blocks or their final
assembly.

Proof boundary: complete displayed-height window for the canonical lower-fold
transform only.  No adjacent-selector crossing, complete source evaluator
error theorem, lower-interior join, `T_upper`, `Lambda<=0`, RH, or prize-level
conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = cell_gate.ode_gate.set_low_priority()
    for dependency in (SOURCE_OUTPUT, PORTCULLIS_GATE, CELL_GATE, cell_gate.HEIGHT_CACHE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    ctx.dps = PRECISION
    certified = certificate(BUILDER_PANELS)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_displayed_height_window_gate",
        "status": "canonical_lower_fold_complete_15_point_displayed_height_window_certified",
        "passed": True,
        "certified_displayed_window": certified,
        "decision": {
            "source_displayed_height_roster_parsed": True,
            "all_15_displayed_heights_inside_one_selector_stable_cell": True,
            "all_15_displayed_heights_preserve_84_mode_turning_roster": True,
            "canonical_lower_fold_window_enclosed": True,
            "printed_hardy_values_used_as_proof_data": False,
            "complete_source_formula_enclosed": False,
            "selector_transition_join_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "source_output": {"path": relative(SOURCE_OUTPUT), "sha256": file_hash(SOURCE_OUTPUT)},
            "portcullis_gate": {"path": relative(PORTCULLIS_GATE), "sha256": file_hash(PORTCULLIS_GATE)},
            "first_nonzero_cell_gate": {"path": relative(CELL_GATE), "sha256": file_hash(CELL_GATE)},
            "height_transport_cache": {"path": relative(cell_gate.HEIGHT_CACHE), "sha256": file_hash(cell_gate.HEIGHT_CACHE)},
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
            "Derive the adjacent-roster identity at an odd-selector crossing, then decide whether "
            "the full selector cell needs coverage or only source-requested windows around each central height."
        ),
        "proof_boundary": (
            "Canonical lower-fold control over all 15 displayed heights only. No complete source formula, "
            "selector-fold join, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built complete displayed lower-fold height window: 15 offsets inside radius=0.07")


if __name__ == "__main__":
    main()
