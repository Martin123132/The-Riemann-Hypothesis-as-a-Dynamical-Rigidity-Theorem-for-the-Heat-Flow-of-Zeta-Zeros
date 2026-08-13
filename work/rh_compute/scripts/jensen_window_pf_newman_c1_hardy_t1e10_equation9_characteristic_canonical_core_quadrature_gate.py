#!/usr/bin/env python3
"""Certify the characteristic canonical core and its finite-t comparison."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx


NORMAL_FORM_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 35
T = 10_000_000_000
C = 159577
MODE_LO = 39853
MODE_HI = 39936
Z_BOUND = 4
Y_BOUND = 64
PANELS = 32
FORMULA_VERSION = "characteristic_core_degree17_tanh_v1"
TOLERANCE = "1e-14"

TANH_COEFFICIENTS = {
    1: Fraction(1),
    3: Fraction(-1, 3),
    5: Fraction(2, 15),
    7: Fraction(-17, 315),
    9: Fraction(62, 2835),
    11: Fraction(-1382, 155925),
    13: Fraction(21844, 6081075),
    15: Fraction(-929569, 638512875),
    17: Fraction(6404582, 10854718875),
}


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


def load_cache() -> dict[tuple[str, int], dict[str, Any]]:
    rows: dict[tuple[str, int], dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    for line in CACHE.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("formula_version") != FORMULA_VERSION:
            continue
        rows[(row["integral_kind"], int(row["panel_index"]))] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


class CharacteristicCore:
    def __init__(self) -> None:
        self.i = acb(0, 1)
        self.pi = arb.pi()
        self.t = arb(T)
        self.c = arb(C)
        self.y = arb(Y_BOUND)
        self.eta = self.pi * self.c**2 / (8 * self.t)
        self.beta = (self.t * self.eta) ** (arb(1) / 3)
        self.lam = (self.eta - 1) * self.t / self.beta
        self.sigma = 4 * self.beta / (self.pi * self.c)
        self.h = 4 * self.beta / self.c
        self.cy = self.sigma / self.c
        self.detunings = [
            self.h * (arb(mode) - self.c / 4)
            for mode in range(MODE_LO, MODE_HI + 1)
        ]
        self.f0_constant = (self.i * self.pi / 4).exp() * (self.pi / 2).sqrt()
        self.f0_rotation = (-self.i * self.pi / 4).exp() / arb(2).sqrt()

    def f0(self, q: acb) -> acb:
        return self.f0_constant * (self.f0_rotation * q).erf()

    def normal_integral(self, a: acb, p: acb, analytic: bool) -> acb:
        root = (2 * a).sqrt()
        q0 = -p * root / (2 * a)
        q1 = (self.y - p / (2 * a)) * root
        return (
            (-self.i * p**2 / (4 * a)).exp()
            / root
            * (self.f0(q1) - self.f0(q0))
        )

    def first_normal_moment(self, a: acb, p: acb, h0: acb) -> acb:
        endpoint = (self.i * (a * self.y**2 - p * self.y)).exp()
        return p * h0 / (2 * a) - self.i * (endpoint - 1) / (2 * a)

    def tanh_polynomial(self, value: acb) -> acb:
        total = acb(0)
        for degree, coefficient in TANH_COEFFICIENTS.items():
            total += arb(coefficient.numerator) / coefficient.denominator * value**degree
        return total

    def argument_minus_tanh_polynomial(self, value: acb) -> acb:
        total = acb(0)
        for degree, coefficient in TANH_COEFFICIENTS.items():
            if degree >= 3:
                total -= arb(coefficient.numerator) / coefficient.denominator * value**degree
        return total

    def canonical_integrand(self, z: acb, analytic: bool) -> acb:
        a = 1 / (4 * self.beta)
        inner = acb(0)
        for detuning in self.detunings:
            inner += self.normal_integral(a, z + detuning, analytic)
        carrier = (self.i * (z**3 / 3 - self.lam * z)).exp()
        return carrier * inner

    def finite_t_surrogate_integrand(self, z: acb, analytic: bool) -> acb:
        scaled = z / self.beta
        tau = self.tanh_polynomial(scaled)
        x = (1 - tau) / 2
        a = x / (2 * self.beta)
        common = self.beta * tau
        inner = acb(0)
        for detuning in self.detunings:
            p = common + detuning
            h0 = self.normal_integral(a, p, analytic)
            inner += h0 + self.cy * self.first_normal_moment(a, p, h0)
        amplitude = scaled.cosh() ** (-arb(3) / 2)
        endpoint_phase = (
            -self.lam * z
            + self.beta**3 * self.argument_minus_tanh_polynomial(scaled)
        )
        return amplitude * (self.i * endpoint_phase).exp() * inner

    def airy_integrand(self, y: acb, _: bool) -> acb:
        airy = (-(self.lam + y)).airy_ai()
        carrier = (self.i * y**2 / (4 * self.beta)).exp()
        kernel = sum(
            ((-self.i * detuning * y).exp() for detuning in self.detunings),
            acb(0),
        )
        return airy * carrier * kernel


def panel_bounds(index: int) -> tuple[arb, arb]:
    left = arb(-Z_BOUND) + arb(2 * Z_BOUND * index) / PANELS
    right = arb(-Z_BOUND) + arb(2 * Z_BOUND * (index + 1)) / PANELS
    return left, right


def integrate_panel(core: CharacteristicCore, kind: str, index: int) -> acb:
    left, right = panel_bounds(index)
    function = (
        core.canonical_integrand
        if kind == "canonical"
        else core.finite_t_surrogate_integrand
    )
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
    for kind in ("canonical", "finite_t_degree17"):
        for index in range(PANELS):
            key = (kind, index)
            if key in rows:
                continue
            started = time.perf_counter()
            value = integrate_panel(core, kind, index)
            left, right = panel_bounds(index)
            row = {
                "formula_version": FORMULA_VERSION,
                "integral_kind": kind,
                "panel_index": index,
                "z_left": str(left),
                "z_right": str(right),
                "value": complex_record(value),
                "precision_decimal_digits": PRECISION,
                "tolerance": TOLERANCE,
                "elapsed_seconds": round(time.perf_counter() - started, 3),
            }
            append_cache(row)
            rows[key] = row
            print(
                f"cached {kind} panel {index + 1}/{PANELS}: "
                f"{row['elapsed_seconds']:.3f}s",
                flush=True,
            )


def sum_panels(rows: dict[tuple[str, int], dict[str, Any]], kind: str) -> acb:
    require(all((kind, index) in rows for index in range(PANELS)), f"incomplete {kind} cache")
    return sum((parse_complex(rows[(kind, index)]["value"]) for index in range(PANELS)), acb(0))


def full_line_compact_y(core: CharacteristicCore) -> acb:
    total = acb(0)
    for left in range(0, Y_BOUND, 4):
        total += acb.integral(
            core.airy_integrand,
            arb(left),
            arb(left + 4),
            abs_tol=arb("1e-24"),
            rel_tol=arb("1e-24"),
            eval_limit=200_000,
            depth_limit=40,
        )
    return 2 * core.pi * total


def certified_values(core: CharacteristicCore, rows: dict[tuple[str, int], dict[str, Any]]) -> dict[str, Any]:
    canonical = sum_panels(rows, "canonical")
    surrogate = sum_panels(rows, "finite_t_degree17")
    full_line = full_line_compact_y(core)

    radius = arb(Z_BOUND) / core.beta
    tanh_tail = arb(3) * radius**18 / (1 - radius)
    phase_tail = (
        core.beta**3 * tanh_tail
        + arb(Y_BOUND) * core.beta * tanh_tail
        + arb(Y_BOUND) ** 2 * tanh_tail / (4 * core.beta)
    )
    integral_tail = (
        arb(2).sqrt()
        / core.pi
        * (1 + core.cy * arb(Y_BOUND))
        * arb(2 * Z_BOUND)
        * arb(Y_BOUND)
        * arb(MODE_HI - MODE_LO + 1)
        * phase_tail
    )
    finite_t = acb(
        arb(surrogate.real, integral_tail),
        arb(surrogate.imag, integral_tail),
    )

    correction = finite_t - canonical
    z_tail = full_line - canonical
    relative_correction = abs(correction) / abs(finite_t)
    relative_z_tail = abs(z_tail) / abs(canonical)
    sqrt_lambda = core.lam.sqrt()
    plus_gap = sqrt_lambda - core.detunings[-1]
    minus_gap = sqrt_lambda + core.detunings[0]

    require(plus_gap > arb("0.0058"), "positive Airy branch margin failed")
    require(minus_gap > arb("0.032"), "negative Airy branch margin failed")
    require(abs(correction) < arb("0.00031"), "finite-t correction exceeds 0.00031")
    require(relative_correction < arb("0.0000066"), "relative correction exceeds 6.6e-6")
    require(relative_z_tail < arb("0.064"), "canonical z-tail ratio exceeds 0.064")
    require(integral_tail < arb("5e-35"), "degree-17 tanh tail exceeds 5e-35")

    return {
        "canonical_rectangle_ball": complex_record(canonical),
        "finite_t_degree17_rectangle_ball": complex_record(surrogate),
        "finite_t_exact_rectangle_ball": complex_record(finite_t),
        "finite_t_minus_canonical_ball": complex_record(correction),
        "finite_t_minus_canonical_absolute_ball": abs(correction).str(PRECISION, more=True),
        "finite_t_relative_correction_ball": relative_correction.str(PRECISION, more=True),
        "full_z_line_compact_y_ball": complex_record(full_line),
        "canonical_abs_z_gt_4_contribution_ball": complex_record(z_tail),
        "canonical_abs_z_gt_4_absolute_ball": abs(z_tail).str(PRECISION, more=True),
        "canonical_abs_z_gt_4_relative_to_rectangle_ball": relative_z_tail.str(PRECISION, more=True),
        "sqrt_lambda_ball": sqrt_lambda.str(PRECISION, more=True),
        "positive_Airy_branch_stationary_gap_ball": plus_gap.str(PRECISION, more=True),
        "negative_Airy_branch_stationary_gap_ball": minus_gap.str(PRECISION, more=True),
        "degree17_tanh_tail_ball": tanh_tail.str(PRECISION, more=True),
        "exact_phase_tail_ball": phase_tail.str(PRECISION, more=True),
        "exact_integral_tail_ball": integral_tail.str(PRECISION, more=True),
        "physical_amplitude_factor_ball": (arb(2).sqrt() / core.pi).str(PRECISION, more=True),
        "physical_finite_t_rectangle_ball": complex_record(arb(2).sqrt() / core.pi * finite_t),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_values"]
    return f"""# Characteristic canonical-core interval quadrature

Date: 2026-08-10

Status: finite characteristic-core quadrature validated; not a proof of the
outer normal-saddle join or the complete Hardy transformation

For the grouped canonical phase from Section 11.314, define

```text
I_M(Z,Y)=sum_m integral_(-Z)^Z integral_0^Y
 exp(i[z^3/3-lambda*z-(z+d_m)y+y^2/(4beta)]) dy dz.     (CQ1)
```

Two exact one-dimensional representations are used.  Integrating `z` on the
full real line first gives

```text
I_M(infinity,Y)=2pi integral_0^Y Ai(-lambda-y)
 exp(i*y^2/(4beta))D_M(y)dy.                            (CQ2)
```

Integrating `y` first on `[0,Y]`, complete the square:

```text
H_beta(p;Y)=exp(-i*beta*p^2)sqrt(2beta)
 integral_(-sqrt(2beta)p)^((Y-2beta*p)/sqrt(2beta))
 exp(i*q^2/2)dq.                                        (CQ3)
```

Then `I_M(Z,Y)` is the `z` integral of
`exp(i[z^3/3-lambda*z]) sum_m H_beta(z+d_m;Y)`.

The detuning lattice gives strict branch margins

```text
sqrt(lambda)-d_max={c['positive_Airy_branch_stationary_gap_ball']} >0.0058,
sqrt(lambda)+d_min={c['negative_Airy_branch_stationary_gap_ball']} >0.032.
```

Thus neither leading Airy branch has a normal stationary point on
`0<=y<=64`.  This is an arithmetic lattice effect; the dangerous continuous
threshold lies just outside the final integer mode.

At `Z=4`, `Y=64`, 32 deterministic interval panels give

```text
I_M(4,64)={c['canonical_rectangle_ball']['real_ball']}
          +i*{c['canonical_rectangle_ball']['imag_ball']}.             (CQ4)
```

The exact finite-`t` logit phase remains quadratic in `y`.  A degree-17 odd
Taylor form for `tanh(z/beta)` removes interval dependency cancellation; a
Cauchy estimate on `|w|=1` bounds its contribution to the full rectangle by

```text
{c['exact_integral_tail_ball']} <5e-35.                 (CQ5)
```

After restoring the exact amplitude, the finite-`t` rectangle divided by the
constant `sqrt(2)/pi` is

```text
{c['finite_t_exact_rectangle_ball']['real_ball']}
+i*{c['finite_t_exact_rectangle_ball']['imag_ball']}.                  (CQ6)
```

Its direct correction to (CQ4) has magnitude and relative size

```text
{c['finite_t_minus_canonical_absolute_ball']} <0.00031,
{c['finite_t_relative_correction_ball']} <6.6e-6.       (CQ7)
```

This cancellation-preserving comparison is much sharper than integrating the
pointwise phase-error majorant absolutely.

The full-`z` value from (CQ2) is

```text
{c['full_z_line_compact_y_ball']['real_ball']}
+i*{c['full_z_line_compact_y_ball']['imag_ball']}.                     (CQ8)
```

Consequently the canonical `|z|>4`, `0<=y<=64` contribution has relative
size

```text
{c['canonical_abs_z_gt_4_relative_to_rectangle_ball']} <0.064.         (CQ9)
```

It is controlled and explicitly measured, but it is not negligible.  The
next theorem must place (CQ9) on a legal steepest-descent contour and join the
`y>64` component to every displaced alpha saddle.

Pi provenance: `pi` in (CQ2) is the standard Airy Fourier inversion factor;
`beta`, `d_m`, and the physical `sqrt(2)/pi` factor retain the Kummer and
Fourier--Poisson normalization traced in Sections 11.308 and 11.314.  No
geometric fit is introduced.

Proof boundary: rigorous finite quadrature and exact-transform algebra at
`t=10^10` only.  No height-uniform canonical bound, `y>64` saddle partition,
complete `T_upper` assembly, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-cache", action="store_true")
    args = parser.parse_args()
    started = time.perf_counter()
    priority = set_low_priority()
    require(NORMAL_FORM_GATE.is_file() and CHECKER.is_file(), "missing dependency or checker")
    if args.rebuild_cache and CACHE.exists():
        CACHE.unlink()

    ctx.dps = PRECISION
    core = CharacteristicCore()
    rows = load_cache()
    fill_cache(core, rows)
    values = certified_values(core, rows)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate",
        "status": "characteristic_canonical_core_and_finite_t_interval_comparison_complete",
        "passed": True,
        "scope": {
            "height": T,
            "modes": f"{MODE_LO}..{MODE_HI}",
            "z_rectangle": f"|z|<={Z_BOUND}",
            "y_rectangle": f"0<=y<={Y_BOUND}",
            "panels_per_integral": PANELS,
        },
        "representations": {
            "Airy_first": "2pi*int_0^Y Ai(-lambda-y) exp(i*y^2/(4beta)) D_M(y)dy",
            "Fresnel_first": "int_-Z^Z exp(i[z^3/3-lambda*z]) sum_m H_beta(z+d_m;Y)dz",
            "finite_t_normal_integral": "The exact logit phase is quadratic in y and the exact amplitude is linear in y, so both are evaluated by H(a,p;Y) and its first p-moment.",
        },
        "certified_values": values,
        "decision": {
            "two_one_dimensional_representations_proved": True,
            "discrete_near_endpoint_stationary_point_excluded": True,
            "canonical_rectangle_rigorously_integrated": True,
            "finite_t_rectangle_rigorously_integrated": True,
            "relative_finite_t_correction_below_6_6e_6": True,
            "canonical_z_tail_negligible": False,
            "outer_normal_saddle_join_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Put the non-negligible canonical |z|>4 term on an explicit steepest-descent contour, then partition y>64 so that every displaced alpha saddle is assigned to either the endpoint chart or the ordinary interior Morse chart without overlap or loss.",
        "proof_boundary": "Exact transform algebra and rigorous finite interval quadrature at t=10^10 only. No height-uniform canonical theorem, y>64 saddle join, complete T_upper assembly, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "normal_form_gate": {"path": relative(NORMAL_FORM_GATE), "sha256": file_hash(NORMAL_FORM_GATE)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
            "cache_rows": len(rows),
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built characteristic canonical-core quadrature: correction<6.6e-6, z-tail<0.064")


if __name__ == "__main__":
    main()
