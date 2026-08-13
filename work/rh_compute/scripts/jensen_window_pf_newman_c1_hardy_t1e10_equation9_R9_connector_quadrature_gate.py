#!/usr/bin/env python3
"""Certify the R=4 to R=9 connector and the true outgoing contour tail."""

from __future__ import annotations

import argparse
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

from jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate import (
    CharacteristicCore,
    parse_complex,
)


CORE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.json"
GEOMETRY_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_saddle_migration_contour_geometry_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_connector_quadrature_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_connector_quadrature_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_connector_quadrature_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 35
PANELS_PER_SIDE = 20
PANEL_WIDTH = arb(1) / 4
FORMULA_VERSION = "R9_connector_degree17_tanh_v1"
TOLERANCE = "1e-14"
Y_BOUND = 64
MODE_COUNT = 84


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def panel_bounds(index: int) -> tuple[arb, arb, str]:
    if index < PANELS_PER_SIDE:
        left = arb(-9) + PANEL_WIDTH * index
        side = "negative"
    else:
        left = arb(4) + PANEL_WIDTH * (index - PANELS_PER_SIDE)
        side = "positive"
    return left, left + PANEL_WIDTH, side


def load_cache() -> dict[tuple[str, int], dict[str, Any]]:
    rows: dict[tuple[str, int], dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    for line in CACHE.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("formula_version") == FORMULA_VERSION:
            rows[(row["integral_kind"], int(row["panel_index"]))] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def integrate_panel(core: CharacteristicCore, kind: str, index: int) -> acb:
    left, right, _ = panel_bounds(index)
    function = core.canonical_integrand if kind == "canonical" else core.finite_t_surrogate_integrand
    return acb.integral(
        function,
        left,
        right,
        abs_tol=arb(TOLERANCE),
        rel_tol=arb(TOLERANCE),
        eval_limit=500_000,
        depth_limit=50,
    )


def fill_cache(core: CharacteristicCore, rows: dict[tuple[str, int], dict[str, Any]]) -> None:
    total_panels = 2 * PANELS_PER_SIDE
    for kind in ("canonical", "finite_t_degree17"):
        for index in range(total_panels):
            key = (kind, index)
            if key in rows:
                continue
            started = time.perf_counter()
            value = integrate_panel(core, kind, index)
            left, right, side = panel_bounds(index)
            row = {
                "formula_version": FORMULA_VERSION,
                "integral_kind": kind,
                "panel_index": index,
                "side": side,
                "z_left": str(left),
                "z_right": str(right),
                "value": complex_record(value),
                "precision_decimal_digits": PRECISION,
                "tolerance": TOLERANCE,
                "elapsed_seconds": round(time.perf_counter() - started, 3),
            }
            append_cache(row)
            rows[key] = row
            print(f"cached {kind} connector {index + 1}/{total_panels}: {row['elapsed_seconds']:.3f}s", flush=True)


def sum_rows(rows: dict[tuple[str, int], dict[str, Any]], kind: str) -> acb:
    total_panels = 2 * PANELS_PER_SIDE
    require(all((kind, index) in rows for index in range(total_panels)), f"incomplete {kind} connector")
    return sum((parse_complex(rows[(kind, index)]["value"]) for index in range(total_panels)), acb(0))


def certified_values(core: CharacteristicCore, rows: dict[tuple[str, int], dict[str, Any]]) -> dict[str, Any]:
    prior = json.loads(CORE_GATE.read_text(encoding="utf-8"))["certified_values"]
    canonical_R4 = parse_complex(prior["canonical_rectangle_ball"])
    finite_R4 = parse_complex(prior["finite_t_exact_rectangle_ball"])
    full_line = parse_complex(prior["full_z_line_compact_y_ball"])
    canonical_connector = sum_rows(rows, "canonical")
    finite_connector_surrogate = sum_rows(rows, "finite_t_degree17")

    radius = arb(9) / core.beta
    tanh_tail = arb(3) * radius**18 / (1 - radius)
    phase_tail = (
        core.beta**3 * tanh_tail
        + arb(Y_BOUND) * core.beta * tanh_tail
        + arb(Y_BOUND) ** 2 * tanh_tail / (4 * core.beta)
    )
    connector_integral_tail = (
        arb(2).sqrt() / core.pi
        * (1 + core.cy * arb(Y_BOUND))
        * arb(10) * arb(Y_BOUND) * arb(MODE_COUNT) * phase_tail
    )
    finite_connector = acb(
        arb(finite_connector_surrogate.real, connector_integral_tail),
        arb(finite_connector_surrogate.imag, connector_integral_tail),
    )

    canonical_R9 = canonical_R4 + canonical_connector
    finite_R9 = finite_R4 + finite_connector
    finite_correction = finite_R9 - canonical_R9
    true_tail = full_line - canonical_R9
    connector_fraction = abs(canonical_connector) / abs(canonical_R9)
    correction_fraction = abs(finite_correction) / abs(finite_R9)
    tail_fraction = abs(true_tail) / abs(canonical_R9)

    require(connector_integral_tail < arb("1e-27"), "connector tanh tail exceeds 1e-27")
    require(tail_fraction < arb("0.01"), "true R9 tail exceeds one percent")
    require(correction_fraction < arb("0.01"), "expanded finite-t correction exceeds one percent")

    return {
        "canonical_connector_ball": complex_record(canonical_connector),
        "finite_t_connector_surrogate_ball": complex_record(finite_connector_surrogate),
        "finite_t_connector_exact_ball": complex_record(finite_connector),
        "canonical_R9_rectangle_ball": complex_record(canonical_R9),
        "finite_t_R9_rectangle_ball": complex_record(finite_R9),
        "finite_t_R9_minus_canonical_ball": complex_record(finite_correction),
        "finite_t_R9_relative_correction_ball": correction_fraction.str(PRECISION, more=True),
        "connector_relative_to_R9_ball": connector_fraction.str(PRECISION, more=True),
        "true_R9_outgoing_tail_ball": complex_record(true_tail),
        "true_R9_outgoing_tail_absolute_ball": abs(true_tail).str(PRECISION, more=True),
        "true_R9_outgoing_tail_relative_ball": tail_fraction.str(PRECISION, more=True),
        "degree17_connector_tanh_tail_ball": tanh_tail.str(PRECISION, more=True),
        "degree17_connector_phase_tail_ball": phase_tail.str(PRECISION, more=True),
        "degree17_connector_integral_tail_ball": connector_integral_tail.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_values"]
    return f"""# R=9 characteristic connector quadrature

Date: 2026-08-11

Status: finite R=9 connector quadrature validated; not a proof of a
height-uniform contour-tail theorem or the `y>64` saddle join

Section 11.316 shows that the prior `|z|>4` contribution contains migrating
Airy saddles once `y>16-lambda`.  The stationary core must therefore include

```text
-9<=z<=9,       0<=y<=64.                               (R9Q1)
```

Forty quarter-unit connector panels on `[-9,-4] union [4,9]` give

```text
I_connector={c['canonical_connector_ball']['real_ball']}
            +i*{c['canonical_connector_ball']['imag_ball']}.           (R9Q2)
```

Combining (R9Q2) with the certified `R=4` core gives

```text
I_canonical(9,64)={c['canonical_R9_rectangle_ball']['real_ball']}
                  +i*{c['canonical_R9_rectangle_ball']['imag_ball']}.  (R9Q3)
```

The exact finite-`t` connector again uses the degree-17 stable tanh form and
closed Fresnel first moment.  Its omitted Taylor contribution to the complete
connector is below

```text
{c['degree17_connector_integral_tail_ball']} <1e-27.     (R9Q4)
```

The exact finite-`t` expanded rectangle divided by `sqrt(2)/pi` is

```text
I_exact(9,64)={c['finite_t_R9_rectangle_ball']['real_ball']}
              +i*{c['finite_t_R9_rectangle_ball']['imag_ball']},       (R9Q5)
```

and its relative correction to (R9Q3) is

```text
{c['finite_t_R9_relative_correction_ball']}.             (R9Q6)
```

Subtracting (R9Q3) from the independent full-line Airy representation leaves
the true outgoing `R=9` contour tail

```text
I_tail={c['true_R9_outgoing_tail_ball']['real_ball']}
       +i*{c['true_R9_outgoing_tail_ball']['imag_ball']},               (R9Q7)

|I_tail|/|I_canonical(9,64)|
 ={c['true_R9_outgoing_tail_relative_ball']}.                           (R9Q8)
```

Unlike the old `|z|>4` quantity, (R9Q7) contains no canonical stationary
point for `0<=y<=64`.  It is the term that must now be certified directly on
the outgoing `pi/6` and `5pi/6` rays.

Pi provenance is unchanged: the rays are forced by the cubic Airy phase and
all scale factors retain the Kummer/Fourier--Poisson normalization.

Proof boundary: rigorous finite connector quadrature and a saved-height true
tail measurement only.  No direct ray-integral bound, height-uniform theorem,
`y>64` partition, complete `T_upper`, `Lambda<=0`, RH, or prize-level
conclusion is proved.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-cache", action="store_true")
    args = parser.parse_args()
    started = time.perf_counter()
    priority = set_low_priority()
    require(CORE_GATE.is_file() and GEOMETRY_GATE.is_file() and CHECKER.is_file(), "missing dependency or checker")
    if args.rebuild_cache and CACHE.exists():
        CACHE.unlink()
    ctx.dps = PRECISION
    core = CharacteristicCore()
    rows = load_cache()
    fill_cache(core, rows)
    values = certified_values(core, rows)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_connector_quadrature_gate",
        "status": "R9_stationary_connector_and_true_outgoing_tail_interval_quadrature_complete",
        "passed": True,
        "scope": {"connector": "[-9,-4] union [4,9]", "y_range": "0<=y<=64", "panels_per_integral": 40},
        "certified_values": values,
        "decision": {
            "migrating_saddle_connector_rigorously_integrated": True,
            "finite_t_R9_rectangle_rigorously_integrated": True,
            "true_R9_tail_measured": True,
            "true_R9_tail_direct_ray_bound_proved": False,
            "y_gt_64_join_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Evaluate and bound the true R=9 tail directly on the outgoing pi/6 and 5pi/6 rays, then solve and partition the y>64 normal stationary locus for the endpoint/interior handoff.",
        "proof_boundary": "Rigorous finite connector quadrature and a saved-height true-tail measurement only. No direct ray-integral bound, height-uniform theorem, y>64 partition, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "core_gate": {"path": relative(CORE_GATE), "sha256": file_hash(CORE_GATE)},
            "geometry_gate": {"path": relative(GEOMETRY_GATE), "sha256": file_hash(GEOMETRY_GATE)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"workers": 1, "process_priority": priority, "cache_rows": len(rows), "elapsed_seconds": round(time.perf_counter()-started,3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print("built R9 connector quadrature and true outgoing tail measurement")


if __name__ == "__main__":
    main()
