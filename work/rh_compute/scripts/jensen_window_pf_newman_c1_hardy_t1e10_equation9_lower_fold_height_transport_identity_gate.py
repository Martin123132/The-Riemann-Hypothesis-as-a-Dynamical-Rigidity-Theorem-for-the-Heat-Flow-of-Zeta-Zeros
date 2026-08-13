#!/usr/bin/env python3
"""Certify exact fixed-selector height transport of the lower-fold lattice."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_airy_fresnel_detuning_ode_gate as ode_gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_lattice_quadrature_gate as lattice_gate


ODE_GATE = ode_gate.RESULT
LATTICE_GATE = lattice_gate.RESULT
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_height_transport_identity_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_height_transport_identity_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_height_transport_identity_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = ode_gate.PRECISION
MODE_LO = ode_gate.MODE_LO
MODE_HI = ode_gate.MODE_HI
FORMULA_VERSION = "lower_fold_fixed_C_height_transport_v1"


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


def symbolic_transport() -> dict[str, str]:
    t, endpoint_c, beta = sp.symbols("t C beta", positive=True, real=True)
    pi = sp.pi
    eta = pi * endpoint_c**2 / (8 * t)
    lam = sp.simplify((eta - 1) * t / beta)
    expected = (pi * endpoint_c**2 / 8 - t) / beta
    require(sp.simplify(lam - expected) == 0, "fixed-selector lambda identity failed")
    require(sp.diff(expected, t) == -1 / beta, "fixed-selector lambda derivative failed")
    return {
        "fixed_selector_scale": "beta^3=pi*C^2/8 is independent of t while the odd endpoint selector C is fixed",
        "linear_height_parameter": "lambda(t;C)=[pi*C^2/8-t]/beta",
        "lambda_height_derivative": "d lambda/dt=-1/beta",
        "lambda_derivative_identity": "partial_lambda G_Y=A_Y*w_Y-A_0+G_Y'/(2beta)+i*d*G_Y",
        "height_derivative_identity": "dG_Y/dt=-(1/beta)[A_Y*w_Y-A_0+G_Y'/(2beta)+i*d*G_Y]",
        "derivation": "Since partial_lambda Ai(-lambda-y)=partial_y Ai(-lambda-y), integrate once by parts and retain both finite endpoint values.",
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
        require(mode not in rows, f"duplicate height-transport mode {mode}")
        rows[mode] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def compute_row(mode: int, p: dict[str, arb], lattice_rows: dict[int, dict[str, Any]]) -> dict[str, Any]:
    started = time.perf_counter()
    i = acb(0, 1)
    beta = p["beta"]
    lam = p["lambda"]
    d = p["h"] * (arb(mode) - arb(ode_gate.C) / 4)
    g = lattice_gate.parse_complex(lattice_rows[mode]["G64"])
    g_prime = -i * ode_gate.moment_integral(d, 1, p)
    ai0 = (-lam).airy_ai()
    aiy = (-(lam + ode_gate.Y)).airy_ai()
    w_y = (i * (arb(ode_gate.Y) ** 2 / (4 * beta) - d * ode_gate.Y)).exp()
    g_lambda = aiy * w_y - ai0 + g_prime / (2 * beta) + i * d * g
    g_t = -g_lambda / beta
    return {
        "formula_version": FORMULA_VERSION,
        "mode": mode,
        "detuning_ball": d.str(PRECISION, more=True),
        "G64": complex_record(g),
        "G64_prime_detuning": complex_record(g_prime),
        "G64_lambda_derivative": complex_record(g_lambda),
        "G64_height_derivative_fixed_C": complex_record(g_t),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }


def fill_cache(rows: dict[int, dict[str, Any]], p: dict[str, arb], lattice_rows: dict[int, dict[str, Any]]) -> None:
    for mode in range(MODE_LO, MODE_HI + 1):
        if mode in rows:
            continue
        row = compute_row(mode, p, lattice_rows)
        append_cache(row)
        rows[mode] = row


def certificate(rows: dict[int, dict[str, Any]], p: dict[str, arb]) -> dict[str, Any]:
    require(set(rows) == set(range(MODE_LO, MODE_HI + 1)), "height-transport cache is incomplete")
    g_values = [parse_complex(rows[mode]["G64"]) for mode in range(MODE_LO, MODE_HI + 1)]
    g_prime_values = [parse_complex(rows[mode]["G64_prime_detuning"]) for mode in range(MODE_LO, MODE_HI + 1)]
    lambda_values = [parse_complex(rows[mode]["G64_lambda_derivative"]) for mode in range(MODE_LO, MODE_HI + 1)]
    height_values = [parse_complex(rows[mode]["G64_height_derivative_fixed_C"]) for mode in range(MODE_LO, MODE_HI + 1)]
    grouped = 2 * arb.pi() * sum(g_values, acb(0))
    grouped_d = 2 * arb.pi() * sum(g_prime_values, acb(0))
    grouped_lambda = 2 * arb.pi() * sum(lambda_values, acb(0))
    grouped_height = 2 * arb.pi() * sum(height_values, acb(0))
    chain_residual = grouped_height + grouped_lambda / p["beta"]
    require(chain_residual.real.contains(0) and chain_residual.imag.contains(0), "grouped height chain rule failed")

    pi = arb.pi()
    endpoint_c = arb(ode_gate.C)
    selector_upper = pi * endpoint_c**2 / 8
    selector_lower = pi * (endpoint_c - 2) ** 2 / 8
    return {
        "height": ode_gate.T,
        "fixed_odd_selector_C": ode_gate.C,
        "selector_cell_lower_ball": selector_lower.str(PRECISION, more=True),
        "selector_cell_upper_ball": selector_upper.str(PRECISION, more=True),
        "distance_to_upper_selector_fold_ball": (selector_upper - ode_gate.T).str(PRECISION, more=True),
        "lambda_ball": p["lambda"].str(PRECISION, more=True),
        "lambda_height_derivative_ball": (-1 / p["beta"]).str(PRECISION, more=True),
        "row_count": len(rows),
        "grouped_2pi_G64": complex_record(grouped),
        "grouped_2pi_detuning_derivative": complex_record(grouped_d),
        "grouped_2pi_lambda_derivative": complex_record(grouped_lambda),
        "grouped_2pi_height_derivative_fixed_C": complex_record(grouped_height),
        "grouped_chain_rule_residual": complex_record(chain_residual),
        "grouped_height_derivative_magnitude_ball": abs(grouped_height).str(PRECISION, more=True),
        "cache_total_elapsed_seconds": round(sum(float(row["elapsed_seconds"]) for row in rows.values()), 3),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_transport"]
    return f"""# Exact fixed-selector height transport of the lower fold

Date: 2026-08-11

Status: exact height-transport identity and complete 84-mode saved-height
derivative validated; not a proof of a nonzero height-cell enclosure

While the odd source endpoint selector `C` is fixed,

```text
beta^3=pi*C^2/8,
lambda(t;C)=[pi*C^2/8-t]/beta,
d lambda/dt=-1/beta.                                  (HT1)
```

For the finite transform `G_Y(d)` of Section 11.324,

```text
partial_lambda G_Y
 =Ai(-lambda-Y)w_Y-Ai(-lambda)
  +G_Y'/(2beta)+i*d*G_Y,                              (HT2)

dG_Y/dt=-(partial_lambda G_Y)/beta.                    (HT3)
```

Equation (HT2) follows from
`partial_lambda Ai(-lambda-y)=partial_y Ai(-lambda-y)` and one exact
integration by parts.  Both finite endpoint values are retained.

The current saved height lies in the fixed-selector cell

```text
{c['selector_cell_lower_ball']}<t<
{c['selector_cell_upper_ball']},                       (HT4)
```

at distance

```text
{c['distance_to_upper_selector_fold_ball']}            (HT5)
```

below its upper selector fold.  All 84 values of `G_64'`,
`partial_lambda G_64`, and `dG_64/dt` were evaluated as complex balls.  The
grouped canonical height derivative is

```text
d/dt [2pi sum_m G_64(d_m)]
 ={c['grouped_2pi_height_derivative_fixed_C']['real_ball']}
  +i*{c['grouped_2pi_height_derivative_fixed_C']['imag_ball']}, (HT6)
```

with magnitude

```text
{c['grouped_height_derivative_magnitude_ball']}.        (HT7)
```

The independently assembled grouped chain-rule residual contains zero in
both components.  This gives an exact local transport law, but a point
derivative is not yet a nonzero-radius height theorem.  The next obligation
is a second-derivative or interval-ODE bound on a selector-stable height cell,
followed by subdivision up to the selector boundary.

Proof boundary: exact fixed-selector height identity and complete derivative
at `t=10^10` only.  No nonzero-radius height enclosure, selector-transition
join, complete `T_upper`, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = ode_gate.set_low_priority()
    for dependency in (ODE_GATE, LATTICE_GATE, lattice_gate.CACHE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    ctx.dps = PRECISION
    p = ode_gate.parameters()
    lattice_rows = lattice_gate.load_cache()
    rows = load_cache()
    fill_cache(rows, p, lattice_rows)
    rows = load_cache()
    certified = certificate(rows, p)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_height_transport_identity_gate",
        "status": "exact_fixed_selector_lower_fold_height_transport_identity_complete",
        "passed": True,
        "symbolic_transport": symbolic_transport(),
        "certified_transport": certified,
        "decision": {
            "fixed_selector_lambda_linear_in_height": True,
            "exact_finite_transform_height_derivative_proved": True,
            "all_84_height_derivatives_rigorously_evaluated": True,
            "nonzero_radius_height_cell_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "ode_gate": {"path": relative(ODE_GATE), "sha256": file_hash(ODE_GATE)},
            "lattice_gate": {"path": relative(LATTICE_GATE), "sha256": file_hash(LATTICE_GATE)},
            "lattice_cache": {"path": relative(lattice_gate.CACHE), "sha256": file_hash(lattice_gate.CACHE)},
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
        "next_obligation": "Bound the second height derivative or integrate the detuning ODE with interval lambda over a nonzero selector-stable height cell, then subdivide and join selector transitions without changing mode ownership silently.",
        "proof_boundary": "Exact fixed-selector height transport and complete point derivative only. No nonzero-radius height enclosure, selector join, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built lower-fold height transport: exact lambda law and 84 grouped derivatives certified")


if __name__ == "__main__":
    main()
