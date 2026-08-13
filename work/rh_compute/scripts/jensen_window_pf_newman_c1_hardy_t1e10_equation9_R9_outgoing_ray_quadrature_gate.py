#!/usr/bin/env python3
"""Certify the canonical R=9 outgoing Airy-ray integral directly."""

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


GEOMETRY_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_saddle_migration_contour_geometry_gate.json"
CONNECTOR_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_connector_quadrature_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_outgoing_ray_quadrature_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_outgoing_ray_quadrature_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_outgoing_ray_quadrature_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 30
RADIUS = 9
S_CUTOFF = 2
S_PANELS = 16
Y_BOUND = 64
Y_PANELS = 8
MODE_COUNT = 84
S_PANEL_WIDTH = arb(S_CUTOFF) / S_PANELS
Y_PANEL_WIDTH = arb(Y_BOUND) / Y_PANELS
FORMULA_VERSION = "R9_direct_Airy_rays_combined_phase_v1"
OUTER_TOLERANCE = "1e-11"
INNER_TOLERANCE = "1e-12"


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


def ray_data(core: CharacteristicCore, side: str) -> tuple[acb, acb, acb]:
    if side == "positive":
        direction = (core.i * core.pi / 6).exp()
        return acb(RADIUS), direction, direction
    require(side == "negative", f"unknown ray side {side}")
    direction = (core.i * 5 * core.pi / 6).exp()
    return acb(-RADIUS), direction, -direction


def combined_ray_integrand(core: CharacteristicCore, side: str, s: acb) -> acb:
    base, direction, orientation = ray_data(core, side)
    z = base + direction * s

    def y_integrand(y: acb, _: bool) -> acb:
        phase = z**3 / 3 - core.lam * z - z * y + y**2 / (4 * core.beta)
        kernel = sum(
            ((-core.i * detuning * y).exp() for detuning in core.detunings),
            acb(0),
        )
        return orientation * (core.i * phase).exp() * kernel

    total = acb(0)
    for index in range(Y_PANELS):
        left = Y_PANEL_WIDTH * index
        right = Y_PANEL_WIDTH * (index + 1)
        total += acb.integral(
            y_integrand,
            left,
            right,
            abs_tol=arb(INNER_TOLERANCE),
            rel_tol=arb(INNER_TOLERANCE),
            eval_limit=100_000,
            depth_limit=30,
        )
    return total


def panel_bounds(index: int) -> tuple[arb, arb]:
    return S_PANEL_WIDTH * index, S_PANEL_WIDTH * (index + 1)


def integrate_ray_panel(core: CharacteristicCore, side: str, index: int) -> acb:
    left, right = panel_bounds(index)
    return acb.integral(
        lambda s, _: combined_ray_integrand(core, side, s),
        left,
        right,
        abs_tol=arb(OUTER_TOLERANCE),
        rel_tol=arb(OUTER_TOLERANCE),
        eval_limit=100_000,
        depth_limit=30,
    )


def load_cache() -> dict[tuple[str, int], dict[str, Any]]:
    rows: dict[tuple[str, int], dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    for line in CACHE.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("formula_version") == FORMULA_VERSION:
            rows[(row["side"], int(row["panel_index"]))] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def fill_cache(core: CharacteristicCore, rows: dict[tuple[str, int], dict[str, Any]]) -> None:
    for side in ("positive", "negative"):
        for index in range(S_PANELS):
            key = (side, index)
            if key in rows:
                continue
            started = time.perf_counter()
            value = integrate_ray_panel(core, side, index)
            left, right = panel_bounds(index)
            row = {
                "formula_version": FORMULA_VERSION,
                "side": side,
                "panel_index": index,
                "s_left": str(left),
                "s_right": str(right),
                "value": complex_record(value),
                "precision_decimal_digits": PRECISION,
                "outer_tolerance": OUTER_TOLERANCE,
                "inner_tolerance": INNER_TOLERANCE,
                "elapsed_seconds": round(time.perf_counter() - started, 3),
            }
            append_cache(row)
            rows[key] = row
            print(
                f"cached {side} ray {index + 1}/{S_PANELS}: "
                f"{row['elapsed_seconds']:.3f}s",
                flush=True,
            )


def sum_rows(rows: dict[tuple[str, int], dict[str, Any]], side: str) -> acb:
    require(all((side, index) in rows for index in range(S_PANELS)), f"incomplete {side} ray")
    return sum(
        (parse_complex(rows[(side, index)]["value"]) for index in range(S_PANELS)),
        acb(0),
    )


def tail_envelope(core: CharacteristicCore) -> dict[str, arb]:
    s = arb(S_CUTOFF)
    linear = (arb(RADIUS) ** 2 - core.lam - arb(Y_BOUND)) / 2
    quadratic = arb(3).sqrt() * arb(RADIUS) / 2
    q = linear * s + quadratic * s**2 + s**3 / 3
    derivative = linear + 2 * quadratic * s + s**2
    bound = 2 * arb(MODE_COUNT) * arb(Y_BOUND) * (-q).exp() / derivative
    return {
        "linear": linear,
        "quadratic": quadratic,
        "q_cutoff": q,
        "q_prime_cutoff": derivative,
        "both_rays_bound": bound,
    }


def certified_values(core: CharacteristicCore, rows: dict[tuple[str, int], dict[str, Any]]) -> dict[str, Any]:
    positive = sum_rows(rows, "positive")
    negative = sum_rows(rows, "negative")
    truncated = positive + negative
    envelope = tail_envelope(core)
    bound = envelope["both_rays_bound"]
    direct = acb(arb(truncated.real, bound), arb(truncated.imag, bound))

    prior = json.loads(CONNECTOR_GATE.read_text(encoding="utf-8"))["certified_values"]
    measured = parse_complex(prior["true_R9_outgoing_tail_ball"])
    real_overlap = bool(direct.real.overlaps(measured.real))
    imag_overlap = bool(direct.imag.overlaps(measured.imag))
    require(envelope["linear"] > arb("5.9"), "ray linear decay margin lost")
    require(bound < arb("1e-14"), "analytic ray tail exceeds 1e-14")
    require(real_overlap and imag_overlap, "direct ray enclosure misses independent tail measurement")

    return {
        "positive_ray_truncated_ball": complex_record(positive),
        "negative_ray_truncated_ball": complex_record(negative),
        "combined_ray_truncated_ball": complex_record(truncated),
        "combined_ray_complete_enclosure_ball": complex_record(direct),
        "independent_real_axis_tail_ball": complex_record(measured),
        "direct_minus_independent_ball": complex_record(direct - measured),
        "componentwise_overlap_with_independent": {
            "real": real_overlap,
            "imag": imag_overlap,
        },
        "linear_decay_margin_ball": envelope["linear"].str(PRECISION, more=True),
        "quadratic_decay_coefficient_ball": envelope["quadratic"].str(PRECISION, more=True),
        "q_at_s_cutoff_ball": envelope["q_cutoff"].str(PRECISION, more=True),
        "q_prime_at_s_cutoff_ball": envelope["q_prime_cutoff"].str(PRECISION, more=True),
        "both_rays_s_gt_2_absolute_bound_ball": bound.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_values"]
    return f"""# R=9 outgoing Airy-ray quadrature

Date: 2026-08-11

Status: direct saved-height outgoing-ray enclosure validated; not a proof of
a height-uniform contour theorem or the `y>64` saddle join

For `a=lambda+y`, deform the positive and negative real tails to

```text
z_+(s)=9+exp(i*pi/6)s,
z_-(s)=-9+exp(5i*pi/6)s,       s>=0.                    (R9R1)
```

The negative real tail has the opposite outward orientation, so the complete
ray integral is

```text
I_ray= exp(i*pi/6) integral_0^infinity F(z_+(s))ds
      -exp(5i*pi/6) integral_0^infinity F(z_-(s))ds.    (R9R2)
```

Here `F` retains the complete 84-mode kernel inside the `y` integral.  The
phase is combined before exponentiation; no exponentially large closed-form
Fresnel terms are subtracted numerically.

On either ray,

```text
|exp(i phi_(lambda+y))|
 =exp(-[(81-lambda-y)s/2+(9sqrt(3))s^2/2+s^3/3]).       (R9R3)
```

After taking absolute values only for `s>2`, convexity gives the two-ray
envelope

```text
2*84*64*exp(-q(2))/q'(2)
 ={c['both_rays_s_gt_2_absolute_bound_ball']} <1e-14,   (R9R4)

q(s)=(81-lambda-64)s/2+(9sqrt(3))s^2/2+s^3/3.          (R9R5)
```

Thirty-two resumable interval panels on `0<=s<=2` give

```text
I_ray,trunc={c['combined_ray_truncated_ball']['real_ball']}
            +i*{c['combined_ray_truncated_ball']['imag_ball']}.        (R9R6)
```

Adding (R9R4) componentwise gives the complete enclosure

```text
I_ray={c['combined_ray_complete_enclosure_ball']['real_ball']}
      +i*{c['combined_ray_complete_enclosure_ball']['imag_ball']}.     (R9R7)
```

Both components of (R9R7) overlap the independently obtained full-line minus
`R=9` real-axis ball.  This closes the direct contour check without using that
independent value as numerical input to the ray integration.

Pi provenance: the angles in (R9R1) are forced by the cubic Airy phase.  The
remaining pi factors retain the Kummer/Fourier--Poisson normalization; no
circle, polygon, or fitted geometric constant is introduced.

Proof boundary: rigorous direct contour quadrature and analytic `s>2` tail at
one saved height only.  No height-uniform ray theorem, `y>64` endpoint/interior
partition, complete `T_upper`, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-cache", action="store_true")
    args = parser.parse_args()
    started = time.perf_counter()
    priority = set_low_priority()
    require(GEOMETRY_GATE.is_file() and CONNECTOR_GATE.is_file() and CHECKER.is_file(), "missing dependency or checker")
    if args.rebuild_cache and CACHE.exists():
        CACHE.unlink()
    ctx.dps = PRECISION
    core = CharacteristicCore()
    rows = load_cache()
    fill_cache(core, rows)
    values = certified_values(core, rows)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_outgoing_ray_quadrature_gate",
        "status": "R9_direct_outgoing_Airy_ray_interval_enclosure_complete",
        "passed": True,
        "scope": {
            "rays": ["9+exp(i*pi/6)s", "-9+exp(5i*pi/6)s"],
            "s_numerical_range": "0<=s<=2",
            "y_range": "0<=y<=64",
            "ray_panels": 32,
            "inner_y_panels_per_evaluation": Y_PANELS,
        },
        "exact_geometry": {
            "positive_orientation": "+exp(i*pi/6)",
            "negative_orientation": "-exp(5i*pi/6)",
            "ray_decay": "Im phi_(lambda+y)=(81-lambda-y)s/2+(9sqrt(3))s^2/2+s^3/3",
            "arc_guard": "On the connecting sectors based at +/-9, every linear, quadratic, and cubic imaginary-phase term is nonnegative and the cubic term is positive away from the real boundary.",
        },
        "certified_values": values,
        "decision": {
            "combined_phase_ray_quadrature_rigorous": True,
            "analytic_s_gt_2_tail_bound_rigorous": True,
            "direct_enclosure_matches_independent_real_axis_tail": True,
            "saved_height_R9_tail_directly_certified": True,
            "height_uniform_ray_theorem_proved": False,
            "y_gt_64_join_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Promote the ray enclosure, then solve the y>64 moving normal stationary locus and define the no-overlap endpoint/interior handoff before attempting a height-uniform assembly.",
        "proof_boundary": "Rigorous direct contour quadrature and analytic s>2 tail at one saved height only. No height-uniform ray theorem, y>64 endpoint/interior partition, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "geometry_gate": {"path": relative(GEOMETRY_GATE), "sha256": file_hash(GEOMETRY_GATE)},
            "connector_gate": {"path": relative(CONNECTOR_GATE), "sha256": file_hash(CONNECTOR_GATE)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "cache_rows": len(rows),
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built direct R9 outgoing-ray enclosure")


if __name__ == "__main__":
    main()
