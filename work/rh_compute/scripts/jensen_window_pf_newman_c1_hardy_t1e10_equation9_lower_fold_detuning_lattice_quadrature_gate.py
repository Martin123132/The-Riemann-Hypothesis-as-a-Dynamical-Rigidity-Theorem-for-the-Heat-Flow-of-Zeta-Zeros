#!/usr/bin/env python3
"""Certify all 84 finite Airy-Fresnel detuning values and their grouped sum."""

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


ODE_GATE = ode_gate.RESULT
CORE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_lattice_quadrature_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_lattice_quadrature_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_lattice_quadrature_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = ode_gate.PRECISION
MODE_LO = ode_gate.MODE_LO
MODE_HI = ode_gate.MODE_HI
FORMULA_VERSION = "finite_Airy_Fresnel_G64_detuning_lattice_v1"


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


def load_cache() -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    for line in CACHE.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("formula_version") != FORMULA_VERSION:
            continue
        mode = int(row["mode"])
        require(mode not in rows, f"duplicate cached mode {mode}")
        rows[mode] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def compute_row(mode: int, p: dict[str, arb]) -> dict[str, Any]:
    started = time.perf_counter()
    d = p["h"] * (arb(mode) - arb(ode_gate.C) / 4)
    value = ode_gate.moment_integral(d, 0, p)
    return {
        "formula_version": FORMULA_VERSION,
        "mode": mode,
        "detuning_ball": d.str(PRECISION, more=True),
        "G64": complex_record(value),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }


def fill_cache(rows: dict[int, dict[str, Any]], p: dict[str, arb]) -> None:
    for mode in range(MODE_LO, MODE_HI + 1):
        if mode in rows:
            continue
        row = compute_row(mode, p)
        append_cache(row)
        rows[mode] = row


def target_from_core() -> acb:
    core = json.loads(CORE_GATE.read_text(encoding="utf-8"))
    record = core["certified_values"]["full_z_line_compact_y_ball"]
    return parse_complex(record)


def build_certificate(rows: dict[int, dict[str, Any]], p: dict[str, arb]) -> dict[str, Any]:
    require(set(rows) == set(range(MODE_LO, MODE_HI + 1)), "detuning cache is incomplete")
    values = [parse_complex(rows[mode]["G64"]) for mode in range(MODE_LO, MODE_HI + 1)]
    grouped = 2 * arb.pi() * sum(values, acb(0))
    target = target_from_core()
    require(grouped.real.overlaps(target.real), "grouped real component misses Airy-first target")
    require(grouped.imag.overlaps(target.imag), "grouped imaginary component misses Airy-first target")
    midpoint_differences = {
        "real": (grouped.real.mid() - target.real.mid()).str(PRECISION, more=True),
        "imag": (grouped.imag.mid() - target.imag.mid()).str(PRECISION, more=True),
    }
    magnitudes = [abs(value) for value in values]
    maximum_index = max(range(len(values)), key=lambda index: float(magnitudes[index].mid()))
    adjacent_variation = sum((abs(values[index + 1] - values[index]) for index in range(len(values) - 1)), arb(0))
    return {
        "height": ode_gate.T,
        "Y": ode_gate.Y,
        "mode_range": [MODE_LO, MODE_HI],
        "row_count": len(rows),
        "grouped_2pi_sum": complex_record(grouped),
        "independent_Airy_first_target": complex_record(target),
        "midpoint_differences": midpoint_differences,
        "maximum_single_G64_mode": MODE_LO + maximum_index,
        "maximum_single_G64_magnitude_ball": magnitudes[maximum_index].str(PRECISION, more=True),
        "discrete_adjacent_variation_ball": adjacent_variation.str(PRECISION, more=True),
        "cache_total_elapsed_seconds": round(sum(float(row["elapsed_seconds"]) for row in rows.values()), 3),
        "detuning_spacing_ball": p["h"].str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_lattice"]
    return f"""# Complete lower-fold detuning-lattice quadrature

Date: 2026-08-11

Status: all 84 finite Airy-Fresnel detuning values and their grouped sum
validated; not a proof of the height-uniform fold theorem

For every mode `m=39853..39936`, the cache contains an independent complex
ball for

```text
G_64(d_m)=integral_0^64 Ai(-lambda-y)
 exp(i*y^2/(4beta)-i*d_m*y)dy.                         (DL1)
```

All 84 rows are deterministic, append-only, and individually resumable.
Their exact grouped canonical value is

```text
2pi sum_m G_64(d_m)
 ={c['grouped_2pi_sum']['real_ball']}
  +i*{c['grouped_2pi_sum']['imag_ball']}.              (DL2)
```

Both components overlap the independently computed Airy-first value

```text
{c['independent_Airy_first_target']['real_ball']}
 +i*{c['independent_Airy_first_target']['imag_ball']}. (DL3)
```

The midpoint differences are

```text
real: {c['midpoint_differences']['real']},
imag: {c['midpoint_differences']['imag']}.              (DL4)
```

This cross-check uses two different organizations of the same canonical
fold: 84 separate detuning transforms in (DL2), and one Airy integral with
the complete finite Dirichlet kernel in (DL3).  It supplies rigorous lattice
state data for the detuning ODE from Section 11.324.

The discrete adjacent variation is recorded as

```text
sum_m |G_64(d_(m+1))-G_64(d_m)|
 ={c['discrete_adjacent_variation_ball']}.              (DL5)
```

It is diagnostic state information, not a termwise final-error bound.  The
next step is to propagate the ODE over detuning and height cells while
preserving the grouped finite sum, then match the outer-Morse limit.

Proof boundary: complete rigorous canonical detuning lattice at one saved
height and `Y=64` only.  No finite-t outer-chart join, complete `T_upper`,
height-uniform source theorem, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = ode_gate.set_low_priority()
    for dependency in (ODE_GATE, CORE_GATE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    ctx.dps = PRECISION
    p = ode_gate.parameters()
    rows = load_cache()
    fill_cache(rows, p)
    rows = load_cache()
    certificate = build_certificate(rows, p)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_lattice_quadrature_gate",
        "status": "complete_saved_height_lower_fold_detuning_lattice_certified",
        "passed": True,
        "certified_lattice": certificate,
        "decision": {
            "all_84_detuning_values_rigorously_integrated": True,
            "grouped_sum_overlaps_independent_Airy_first_value": True,
            "detuning_ODE_state_lattice_complete": True,
            "height_uniform_propagation_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "ode_gate": {"path": relative(ODE_GATE), "sha256": file_hash(ODE_GATE)},
            "core_gate": {"path": relative(CORE_GATE), "sha256": file_hash(CORE_GATE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_obligation": "Use the certified 84-point state as anchors for an interval detuning/height propagation of the finite ODE, preserving the grouped sum and matching the ordinary Morse carrier at the chart boundary.",
        "proof_boundary": "Complete rigorous canonical detuning lattice at t=10^10 and Y=64 only. No height-uniform ODE propagation, outer join, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built complete detuning lattice: 84 rows, grouped sum overlaps Airy-first target")


if __name__ == "__main__":
    main()
