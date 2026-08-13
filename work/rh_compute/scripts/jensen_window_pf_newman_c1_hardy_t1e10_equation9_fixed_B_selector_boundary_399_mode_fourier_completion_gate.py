#!/usr/bin/env python3
"""Certify the fixed-B 399-mode selector-boundary canonical Fourier completion."""

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
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate as scheduling_gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate as selector_gate


SCHEDULING_GATE = scheduling_gate.RESULT
SELECTOR_GATE = selector_gate.RESULT
EVENT_GATE = selector_gate.EVENT_GATE
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 70
ABS_TOL = "1e-20"
FORMULA_VERSION = "fixed_B_selector_boundary_399_mode_fourier_completion_v1"
LOWER_ALPHA = 159_577
UPPER_ALPHA = LOWER_ALPHA + 2
FIXED_B = 5_122_423
CENTER_MODE = (LOWER_ALPHA + 3) // 4
HALF_WIDTH = 199
MODE_LO = CENTER_MODE - HALF_WIDTH
MODE_HI = CENTER_MODE + HALF_WIDTH
MODE_COUNT = MODE_HI - MODE_LO + 1


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


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def parameters() -> dict[str, arb]:
    ctx.dps = PRECISION
    pi = arb.pi()
    lower = arb(LOWER_ALPHA)
    height = pi * lower**2 / 8
    beta = height ** (arb(1) / 3)
    sigma = 4 * beta / (pi * lower)
    h = 4 * beta / lower
    y_width = 2 / sigma
    return {
        "pi": pi,
        "lower": lower,
        "upper": lower + 2,
        "height": height,
        "beta": beta,
        "sigma": sigma,
        "h": h,
        "Y": y_width,
    }


def symbolic_identities() -> dict[str, str]:
    a, beta, k, s = sp.symbols("A beta k s", positive=True, real=True)
    pi = sp.pi
    sigma = 4 * beta / (pi * a)
    h = 4 * beta / a
    y_width = 2 / sigma
    require(sp.simplify(h * y_width - 2 * pi) == 0, "one-period identity failed")
    beta_cube = pi * a**2 / 8
    endpoint_phase = sp.simplify((y_width**2 / (4 * beta) - 3 * pi / 2).subs(beta**3, beta_cube))
    require(sp.simplify(endpoint_phase + pi) == 0, "endpoint phase identity failed")
    mode = (a + 3) / 4 + k
    detuning = sp.simplify(h * (mode - a / 4))
    require(sp.simplify(detuning - h * (k + sp.Rational(3, 4))) == 0, "centered detuning failed")
    return {
        "selector_height": "t*=pi*A^2/8",
        "strip_coordinate": "alpha=A+sigma*y, 0<=y<=Y=2/sigma",
        "scales": "beta=t*^(1/3), sigma=4beta/(pi*A), h=4beta/A",
        "one_period": "hY=2pi",
        "symmetric_mode_reindex": "m=(A+3)/4+k, -199<=k<=199",
        "detuning_reindex": "d_m=h(k+3/4)",
        "coefficient": "G_m=int_0^Y Ai(-y) exp(i*y^2/(4beta)-i*d_m*y)dy=hat(g)(k)",
        "period_function": "g(s)=Y Ai(-Ys) exp(i*Y^2*s^2/(4beta)-3pi*i*s/2)",
        "endpoint_phase": "Y^2/(4beta)-3pi/2=-pi",
        "endpoint_values": "g(0)=Y Ai(0), g(1)=-Y Ai(-Y)",
        "fourier_midpoint": "lim_(M->infinity) sum_(k=-M)^M hat(g)(k)=[g(0)+g(1)]/2",
        "endpoint_half_current": "H=[g(0)-g(1)]/2",
        "canonical_completion": "H+sum_core+(midpoint-sum_core)=g(0)",
        "fixed_shift_guard": "A fixed shift of a symmetric Fourier cutoff has the same limit because hat(g)(k)->0.",
    }


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


def canonical_coefficient(mode: int, p: dict[str, arb], tolerance: str = ABS_TOL) -> acb:
    i = acb(0, 1)
    d = p["h"] * (arb(mode) - p["lower"] / 4)

    def integrand(y: acb, _: bool) -> acb:
        return (-y).airy_ai() * (i * (y**2 / (4 * p["beta"]) - d * y)).exp()

    return acb.integral(
        integrand,
        arb(0),
        p["Y"],
        abs_tol=arb(tolerance),
        rel_tol=arb(tolerance),
        eval_limit=600_000,
        depth_limit=60,
    )


def fill_cache(rows: dict[int, dict[str, Any]], p: dict[str, arb]) -> None:
    missing = [mode for mode in range(MODE_LO, MODE_HI + 1) if mode not in rows]
    for index, mode in enumerate(missing, 1):
        started = time.perf_counter()
        d = p["h"] * (arb(mode) - p["lower"] / 4)
        value = canonical_coefficient(mode, p)
        row = {
            "formula_version": FORMULA_VERSION,
            "mode": mode,
            "centered_index": mode - CENTER_MODE,
            "detuning_ball": d.str(PRECISION, more=True),
            "coefficient": complex_record(value),
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        }
        append_cache(row)
        rows[mode] = row
        if index % 20 == 0 or index == len(missing):
            print(f"cached selector coefficients {index}/{len(missing)}", flush=True)


def direct_grouped_integral(p: dict[str, arb], panels: int = MODE_COUNT, tolerance: str = ABS_TOL) -> acb:
    i = acb(0, 1)
    center_offset = arb(CENTER_MODE) - p["lower"] / 4

    def direct_kernel(y: acb) -> acb:
        return sum(
            ((-i * p["h"] * (arb(mode) - p["lower"] / 4) * y).exp() for mode in range(MODE_LO, MODE_HI + 1)),
            acb(0),
        )

    def kernel(y: acb) -> acb:
        if y.real.contains(0) or (y.real - p["Y"]).contains(0):
            return direct_kernel(y)
        angle = p["h"] * y / 2
        return (-i * p["h"] * center_offset * y).exp() * (MODE_COUNT * angle).sin() / angle.sin()

    def integrand(y: acb, _: bool) -> acb:
        return (-y).airy_ai() * (i * y**2 / (4 * p["beta"])).exp() * kernel(y)

    total = acb(0)
    for index in range(panels):
        left = p["Y"] * arb(index) / panels
        right = p["Y"] * arb(index + 1) / panels
        total += acb.integral(
            integrand,
            left,
            right,
            abs_tol=arb(tolerance),
            rel_tol=arb(tolerance),
            eval_limit=200_000,
            depth_limit=40,
        )
    return 2 * p["pi"] * total


def geometry_certificate(p: dict[str, arb]) -> dict[str, Any]:
    require(LOWER_ALPHA % 4 == 1, "lower endpoint is not 1 modulo 4")
    require(MODE_COUNT == 399 and (MODE_LO, MODE_HI) == (39696, 40094), "399-mode roster drift")
    require(CENTER_MODE == 39895, "center mode drift")
    require((p["h"] * p["Y"] - 2 * p["pi"]).contains(0), "hY=2pi failed")

    lower_root = (p["upper"] - (p["upper"] ** 2 - p["lower"] ** 2).sqrt()) / 4
    upper_root = (p["upper"] + (p["upper"] ** 2 - p["lower"] ** 2).sqrt()) / 4
    require(arb(MODE_LO - 1) < lower_root < arb(MODE_LO), "lower stationary root bracket failed")
    require(arb(MODE_HI) < upper_root < arb(MODE_HI + 1), "upper stationary root bracket failed")

    rows: list[tuple[int, arb, arb]] = []
    for mode in range(MODE_LO, MODE_HI + 1):
        alpha = 2 * arb(mode) + p["lower"] ** 2 / (8 * mode)
        lower_clearance = alpha - p["lower"]
        upper_clearance = p["upper"] - alpha
        require(lower_clearance.lower() > 0 and upper_clearance.lower() > 0, f"mode {mode} leaves strip")
        rows.append((mode, lower_clearance, upper_clearance))
    min_lower = min(rows, key=lambda row: row[1].lower())
    min_upper = min(rows, key=lambda row: row[2].lower())

    outside_low = 2 * arb(MODE_LO - 1) + p["lower"] ** 2 / (8 * (MODE_LO - 1)) - p["upper"]
    outside_high = 2 * arb(MODE_HI + 1) + p["lower"] ** 2 / (8 * (MODE_HI + 1)) - p["upper"]
    require(outside_low.lower() > 0 and outside_high.lower() > 0, "first outer mode enters strip")

    threshold_low = -p["Y"].sqrt() + p["Y"] / (2 * p["beta"])
    threshold_high = p["Y"].sqrt() + p["Y"] / (2 * p["beta"])
    d_lo = p["h"] * (arb(MODE_LO) - p["lower"] / 4)
    d_hi = p["h"] * (arb(MODE_HI) - p["lower"] / 4)
    d_outer_low = p["h"] * (arb(MODE_LO - 1) - p["lower"] / 4)
    d_outer_high = p["h"] * (arb(MODE_HI + 1) - p["lower"] / 4)
    core_low_margin = d_lo - threshold_low
    core_high_margin = threshold_high - d_hi
    first_outer_low_margin = threshold_low - d_outer_low
    first_outer_high_margin = d_outer_high - threshold_high
    for margin in (core_low_margin, core_high_margin, first_outer_low_margin, first_outer_high_margin):
        require(margin.lower() > 0, "canonical branch margin lost positivity")
    require((1 / first_outer_low_margin).lower() > 1200, "near-tangent outer-mode guard weakened")

    return {
        "lower_strip_alpha": LOWER_ALPHA,
        "upper_strip_alpha": UPPER_ALPHA,
        "fixed_upper_endpoint_B": FIXED_B,
        "selector_boundary_height_ball": p["height"].str(PRECISION, more=True),
        "beta_ball": p["beta"].str(PRECISION, more=True),
        "sigma_ball": p["sigma"].str(PRECISION, more=True),
        "normal_width_Y_ball": p["Y"].str(PRECISION, more=True),
        "detuning_spacing_h_ball": p["h"].str(PRECISION, more=True),
        "hY_minus_2pi_ball": (p["h"] * p["Y"] - 2 * p["pi"]).str(30, more=True),
        "mode_center": CENTER_MODE,
        "centered_index_range": [-HALF_WIDTH, HALF_WIDTH],
        "stationary_mode_range": [MODE_LO, MODE_HI],
        "stationary_mode_count": MODE_COUNT,
        "continuous_stationary_root_lower_ball": lower_root.str(PRECISION, more=True),
        "continuous_stationary_root_upper_ball": upper_root.str(PRECISION, more=True),
        "minimum_lower_alpha_clearance_mode": min_lower[0],
        "minimum_lower_alpha_clearance_ball": min_lower[1].str(PRECISION, more=True),
        "minimum_upper_alpha_clearance_mode": min_upper[0],
        "minimum_upper_alpha_clearance_ball": min_upper[2].str(PRECISION, more=True),
        "first_lower_outer_alpha_excess_ball": outside_low.str(PRECISION, more=True),
        "first_upper_outer_alpha_excess_ball": outside_high.str(PRECISION, more=True),
        "core_lower_branch_margin_ball": core_low_margin.str(PRECISION, more=True),
        "core_upper_branch_margin_ball": core_high_margin.str(PRECISION, more=True),
        "first_outer_lower_branch_margin_ball": first_outer_low_margin.str(PRECISION, more=True),
        "first_outer_upper_branch_margin_ball": first_outer_high_margin.str(PRECISION, more=True),
        "reciprocal_first_outer_lower_margin_ball": (1 / first_outer_low_margin).str(PRECISION, more=True),
        "near_tangent_outer_mode_requires_grouped_completion": True,
    }


def numerical_certificate(rows: dict[int, dict[str, Any]], p: dict[str, arb]) -> dict[str, Any]:
    require(set(rows) == set(range(MODE_LO, MODE_HI + 1)), "coefficient cache is incomplete")
    coefficients = [parse_complex(rows[mode]["coefficient"]) for mode in range(MODE_LO, MODE_HI + 1)]
    coefficient_sum = sum(coefficients, acb(0))
    core_integral = 2 * p["pi"] * coefficient_sum
    direct_integral = direct_grouped_integral(p)
    require(core_integral.real.overlaps(direct_integral.real), "direct grouped real component misses cache sum")
    require(core_integral.imag.overlaps(direct_integral.imag), "direct grouped imaginary component misses cache sum")

    ai0 = arb(0).airy_ai()
    ai_y = (-p["Y"]).airy_ai()
    endpoint_phase = p["Y"] ** 2 / (4 * p["beta"]) - 3 * p["pi"] / 2
    require((endpoint_phase + p["pi"]).contains(0), "upper period endpoint phase is not -pi")
    g0 = acb(p["Y"] * ai0)
    g1 = acb(-p["Y"] * ai_y)
    midpoint = (g0 + g1) / 2
    endpoint_half_current = (g0 - g1) / 2
    complement = midpoint - coefficient_sum
    closure = endpoint_half_current + coefficient_sum + complement - g0
    require(closure.real.contains(0) and closure.imag.contains(0), "canonical Fourier completion does not close")

    magnitudes = [abs(value) for value in coefficients]
    max_index = max(range(len(coefficients)), key=lambda index: magnitudes[index].upper())
    adjacent_variation = sum((abs(coefficients[index + 1] - coefficients[index]) for index in range(len(coefficients) - 1)), arb(0))
    physical_factor = arb(2).sqrt() / p["pi"]
    return {
        "cached_coefficient_count": len(coefficients),
        "canonical_399_coefficient_sum": complex_record(coefficient_sum),
        "canonical_399_integral_2pi_sum": complex_record(core_integral),
        "independent_direct_kernel_integral": complex_record(direct_integral),
        "endpoint_g0": complex_record(g0),
        "endpoint_g1": complex_record(g1),
        "complete_symmetric_midpoint": complex_record(midpoint),
        "endpoint_half_current": complex_record(endpoint_half_current),
        "zero_negative_outer_positive_complement": complex_record(complement),
        "scaled_complete_symmetric_midpoint_2pi": complex_record(2 * p["pi"] * midpoint),
        "scaled_endpoint_half_current_2pi": complex_record(2 * p["pi"] * endpoint_half_current),
        "scaled_complement_2pi": complex_record(2 * p["pi"] * complement),
        "scaled_closure_residual_2pi": complex_record(2 * p["pi"] * closure),
        "physical_amplitude_factor_sqrt2_over_pi": physical_factor.str(PRECISION, more=True),
        "physical_399_core": complex_record(physical_factor * core_integral),
        "physical_joint_complement": complex_record(physical_factor * 2 * p["pi"] * complement),
        "physical_endpoint_half_current": complex_record(physical_factor * 2 * p["pi"] * endpoint_half_current),
        "physical_completed_removed_term": complex_record(physical_factor * 2 * p["pi"] * g0),
        "maximum_single_coefficient_mode": MODE_LO + max_index,
        "maximum_single_coefficient_magnitude_ball": magnitudes[max_index].str(PRECISION, more=True),
        "adjacent_coefficient_variation_ball": adjacent_variation.str(PRECISION, more=True),
        "cache_total_elapsed_seconds": round(sum(float(row["elapsed_seconds"]) for row in rows.values()), 3),
    }


def render_note(artifact: dict[str, Any]) -> str:
    g = artifact["geometry"]
    n = artifact["numerical_completion"]
    return f"""# Fixed-B selector-boundary 399-mode Fourier completion

Date: 2026-08-12

Status: exact canonical one-period completion and rigorous 399-mode interval
quadrature validated; finite-t Airy-model transport remains open

Hold the theorem endpoint fixed at `B={FIXED_B}` and consider the selector
strip

```text
A={LOWER_ALPHA} <= alpha <= A+2={UPPER_ALPHA},
t*=pi*A^2/8={g['selector_boundary_height_ball']}.       (FC1)
```

Use the exact lower-fold scales

```text
beta=t*^(1/3), sigma=4beta/(pi A), h=4beta/A,
alpha=A+sigma*y, 0<=y<=Y=2/sigma.                      (FC2)
```

They satisfy the exact one-period identity

```text
hY=2pi.                                                 (FC3)
```

The joint saddle is `alpha_m=2m+t*/(pi m)`.  Its two roots at
`alpha=A+2` lie in

```text
{g['continuous_stationary_root_lower_ball']},
{g['continuous_stationary_root_upper_ball']}.           (FC4)
```

Consequently the strip contains exactly the 399 integer saddle modes
`{MODE_LO}..{MODE_HI}`.  They are the symmetric block

```text
m={CENTER_MODE}+k, -199<=k<=199,
d_m=h(k+3/4).                                           (FC5)
```

The closest lower-fold saddle clearance is only
`{g['minimum_lower_alpha_clearance_ball']}` and the first excluded lower
outer mode misses the upper endpoint by only
`{g['first_lower_outer_alpha_excess_ball']}`.  Its reciprocal canonical
branch margin exceeds `{g['reciprocal_first_outer_lower_margin_ball']}`.
Thus a per-mode nonstationary triangle bound is structurally unsuitable.

For each mode define the canonical coefficient

```text
G_m=int_0^Y Ai(-y) exp(i*y^2/(4beta)-i*d_m*y)dy.        (FC6)
```

All 399 coefficients were integrated as independent complex balls in an
append-only cache.  Their grouped value is

```text
2pi sum_(m={MODE_LO})^{MODE_HI} G_m
 ={n['canonical_399_integral_2pi_sum']['real_ball']}
  +i*{n['canonical_399_integral_2pi_sum']['imag_ball']}. (FC7)
```

A second interval calculation using the exact 399-term Dirichlet kernel gives

```text
{n['independent_direct_kernel_integral']['real_ball']}
 +i*{n['independent_direct_kernel_integral']['imag_ball']},       (FC8)
```

and both components overlap.

The key completion is exact.  Put `s=y/Y`.  By (FC3), `G_(39895+k)` is the
`k`th Fourier coefficient of

```text
g(s)=Y Ai(-Ys) exp(i*Y^2*s^2/(4beta)-3pi*i*s/2).       (FC9)
```

The endpoint phase simplifies exactly to `-pi`, so

```text
g(0)=Y Ai(0), g(1)=-Y Ai(-Y).                          (FC10)
```

Dirichlet-Jordan therefore gives the complete symmetric mode sum as
`[g(0)+g(1)]/2`.  The zero, negative, and outer-positive modes, retained as
one cancellation object, are exactly

```text
midpoint-sum_399
 ={n['zero_negative_outer_positive_complement']['real_ball']}
  +i*{n['zero_negative_outer_positive_complement']['imag_ball']}. (FC11)
```

After the common `2pi` factor this complement is

```text
{n['scaled_complement_2pi']['real_ball']}
 +i*{n['scaled_complement_2pi']['imag_ball']}.          (FC12)
```

Adding the endpoint half-current, the 399 block, and (FC11) reconstructs
`g(0)` exactly.  No outer mode is bounded separately and no divergent
endpoint current is detached.

This closes the complete canonical Fourier bookkeeping at the selector
boundary.  It does not yet prove that the exact finite-`t` Kummer strip is
within a useful explicit error of this canonical completion.  That
cancellation-preserving finite-t comparison is the next theorem target.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(CACHE)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No complete finite-t strip error, `T_upper`, `Lambda<=0`, RH, or prize-level
conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for dependency in (SCHEDULING_GATE, SELECTOR_GATE, EVENT_GATE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    p = parameters()
    symbolic = symbolic_identities()
    geometry = geometry_certificate(p)
    rows = load_cache()
    fill_cache(rows, p)
    rows = load_cache()
    numerical = numerical_certificate(rows, p)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate",
        "status": "exact_canonical_399_mode_fourier_completion_and_joint_complement_certified",
        "passed": True,
        "symbolic_identities": symbolic,
        "geometry": geometry,
        "numerical_completion": numerical,
        "decision": {
            "fixed_theorem_endpoint_B": FIXED_B,
            "all_399_joint_saddles_inside_one_cell_strip": True,
            "stationary_roster_is_symmetric_fourier_block": True,
            "canonical_399_mode_integral_rigorously_computed": True,
            "zero_negative_outer_positive_complement_completed_jointly": True,
            "endpoint_half_current_retained": True,
            "per_mode_outer_triangle_bound_rejected": True,
            "finite_t_canonical_error_bounded": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "source_transition_scheduling_gate": {"path": relative(SCHEDULING_GATE), "sha256": file_hash(SCHEDULING_GATE)},
            "selector_strip_gate": {"path": relative(SELECTOR_GATE), "sha256": file_hash(SELECTOR_GATE)},
            "turning_event_gate": {"path": relative(EVENT_GATE), "sha256": file_hash(EVENT_GATE)},
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
        "next_obligation": (
            "Compare the exact finite-t one-cell Kummer chart with the completed canonical Fourier object before absolute values, "
            "using a compact logit core plus contour tails and preserving the endpoint-midpoint identity."
        ),
        "proof_boundary": (
            "Exact one-cell scale/roster algebra, rigorous canonical 399-mode quadrature, and exact joint canonical Fourier complement only. "
            "No complete finite-t canonical error, T_upper theorem, Lambda<=0, RH, or prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built fixed-B selector chart: 399 canonical modes plus exact joint Fourier complement")


if __name__ == "__main__":
    main()
