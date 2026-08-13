#!/usr/bin/env python3
"""Certify exact endpoint-current reassembly across the moving y=64 rosters."""

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


PRIMITIVE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_three_chart_primitive_envelope_gate.json"
ROSTER_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_moving_roster_dirichlet_partition_gate.json"
POISSON_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_endpoint_roster_reassembly_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_endpoint_roster_reassembly_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 90
T = 10_000_000_000
C = 159577
B = 5122421
MODE_LO = 39853
MODE_HI = 39936
Y0 = 64


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


def block_modes(block: dict[str, Any]) -> range:
    if int(block["count"]) == 0:
        return range(0)
    return range(int(block["first"]), int(block["last"]) + 1)


def certified_reassembly() -> dict[str, Any]:
    ctx.dps = PRECISION
    pi = arb.pi()
    endpoint_c = arb(C)
    beta = (pi * endpoint_c**2 / 8) ** (arb(1) / 3)
    sigma = 4 * beta / (pi * endpoint_c)
    h = 4 * beta / endpoint_c
    y_b = arb(B - C) / sigma
    detunings = {
        mode: h * (arb(mode) - endpoint_c / 4)
        for mode in range(MODE_LO, MODE_HI + 1)
    }
    i = acb(0, 1)

    h_over_sigma = h / sigma
    hy_b = h * y_b
    require(h_over_sigma.overlaps(pi), "h/sigma does not overlap pi")
    require(hy_b.overlaps(pi * (B - C)), "upper endpoint lattice phase failed")

    difference = B - C
    require(difference % 4 == 0, "B-C is not divisible by four")
    quarter_gap = difference // 4
    require(C % 2 == 1 and B % 2 == 1, "source endpoints are not odd")
    require(quarter_gap % 2 == 1, "upper common phase is not odd parity")
    upper_character = -1

    y0 = arb(Y0)
    full_y64 = sum(
        ((-i * detunings[mode] * y0).exp() for mode in range(MODE_LO, MODE_HI + 1)),
        acb(0),
    )
    d_lo = detunings[MODE_LO]
    d_hi = detunings[MODE_HI]
    closed_y64 = (
        (-i * (d_lo + d_hi) * y0 / 2).exp()
        * (arb((MODE_HI - MODE_LO + 1)) * h * y0 / 2).sin()
        / (h * y0 / 2).sin()
    )
    require(full_y64.overlaps(closed_y64), "full y=64 Dirichlet form failed")
    y64_magnitude = abs(full_y64)
    require(y64_magnitude < arb("0.612"), "y=64 full-kernel magnitude exceeds 0.612")

    roster = json.loads(ROSTER_GATE.read_text(encoding="utf-8"))["certified_moving_roster"]
    cells = roster["cells"]
    require(len(cells) == 169, "unexpected moving-roster cell count")
    max_reassembly_error = arb(0)
    for cell in cells:
        blocks = (cell["nonstationary"], cell["Fresnel"], cell["Morse"])
        modes = [mode for block in blocks for mode in block_modes(block)]
        require(modes == list(range(MODE_LO, MODE_HI + 1)), "cell modes do not reassemble in order")
        pieces = []
        for block in blocks:
            pieces.append(
                sum(
                    ((-i * detunings[mode] * y0).exp() for mode in block_modes(block)),
                    acb(0),
                )
            )
        error = abs(sum(pieces, acb(0)) - full_y64)
        if error > max_reassembly_error:
            max_reassembly_error = error
    require(max_reassembly_error.contains(0), "cell endpoint kernels do not reassemble exactly")

    upper_full_kernel = upper_character * (MODE_HI - MODE_LO + 1)
    require(upper_full_kernel == -84, "upper full kernel is not -84")

    return {
        "height": T,
        "source_endpoints": [C, B],
        "mode_range": [MODE_LO, MODE_HI],
        "mode_count": MODE_HI - MODE_LO + 1,
        "h_over_sigma_ball": h_over_sigma.str(PRECISION, more=True),
        "pi_ball": pi.str(PRECISION, more=True),
        "y_B_ball": y_b.str(PRECISION, more=True),
        "h_y_B_ball": hy_b.str(PRECISION, more=True),
        "B_minus_C": difference,
        "B_minus_C_over_4": quarter_gap,
        "upper_common_character": upper_character,
        "full_D84_y_B_exact": upper_full_kernel,
        "full_D84_y64": complex_record(full_y64),
        "full_D84_y64_magnitude_ball": y64_magnitude.str(PRECISION, more=True),
        "fixed_roster_cells_checked": len(cells),
        "maximum_y64_cell_reassembly_error_ball": max_reassembly_error.str(PRECISION, more=True),
        "endpoint_current_identity": "sum_m 2E_m(U)/(i*pi*x)=2 exp(i common(U)) D_84(U)/(i*pi*x)",
        "artificial_split_identity": "sum_m K_m(0,y_B)=sum_m K_m(0,64)+sum_m K_m(64,y_B)",
    }


def symbolic_checks() -> dict[str, str]:
    beta, endpoint_c = sp.symbols("beta C", positive=True, real=True)
    pi = sp.pi
    sigma = 4 * beta / (pi * endpoint_c)
    h = 4 * beta / endpoint_c
    require(sp.simplify(h / sigma - pi) == 0, "symbolic h/sigma identity failed")
    mode, endpoint_b = sp.symbols("m B", integer=True)
    y_b = (endpoint_b - endpoint_c) / sigma
    detuning = h * (mode - endpoint_c / 4)
    require(
        sp.simplify(detuning * y_b - pi * (endpoint_b - endpoint_c) * (mode - endpoint_c / 4)) == 0,
        "symbolic endpoint phase failed",
    )
    return {
        "scale_identity": "h/sigma=pi",
        "upper_lattice_identity": "h*y_B=pi*(B-C)",
        "mode_character": "exp(-i*d_m*y_B)=exp(i*pi*C*(B-C)/4) when B-C is even",
        "saved_endpoint_value": "C=159577, B=5122421 imply exp(-i*d_m*y_B)=-1 for every integer m",
        "full_kernel_value": "D_84(y_B)=-84",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_reassembly"]
    return f"""# Exact y=64 endpoint-current roster reassembly

Date: 2026-08-11

Status: exact endpoint-current reassembly and upper lattice character
validated; not a proof of the grouped bulk integrals or `T_upper`

The exact scales satisfy

```text
h/sigma=pi,
y_B=(B-C)/sigma,
h*y_B=pi(B-C).                                         (ER1)
```

Both source endpoints are odd.  At the saved endpoints,

```text
B-C=4962844=4*1240711,                                 (ER2)
```

and `1240711` is odd.  Therefore, for every integer mode,

```text
exp(-i d_m y_B)
 =exp(i*pi*C(B-C)/4)=-1,                               (ER3)

D_84(y_B)=-84.                                         (ER4)
```

This is an exact lattice character, not a numerical large-argument phase
evaluation.  At the artificial split endpoint the complete kernel instead
satisfies

```text
|D_84(64)|={c['full_D84_y64_magnitude_ball']}<0.612.   (ER5)
```

On every one of the 169 open moving-roster cells, the nonstationary,
Fresnel, and Morse mode blocks reassemble in order to `39853..39936`.
Consequently their endpoint currents sum pointwise to the complete kernel:

```text
sum_chart D_chart(U)=D_84(U),       U=64 or y_B.        (ER6)
```

The `U=64` current cancels exactly against the opposite current from
`K(0,64)` by the primitive telescope of Section 11.320.  The `U=y_B`
current is the genuine source `B` endpoint current and must remain paired
with the grouped Fresnel/bulk contribution; (ER4) evaluates its mode
character but does not bound or discard it.

Proof boundary: exact saved-height scale algebra, endpoint characters, and
169-cell endpoint reassembly only.  No grouped bulk cell integral,
lower-interior join, source `T_upper` identification, complete error theorem,
`Lambda<=0`, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for dependency in (PRIMITIVE_GATE, ROSTER_GATE, POISSON_GATE):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    certified = certified_reassembly()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_endpoint_roster_reassembly_gate",
        "status": "exact_y64_endpoint_roster_reassembly_and_upper_lattice_character_complete",
        "passed": True,
        "symbolic": symbolic_checks(),
        "certified_reassembly": certified,
        "decision": {
            "all_169_cell_endpoint_rosters_reassemble": True,
            "artificial_y64_current_cancels_by_telescope": True,
            "upper_D84_character_equals_minus_84": True,
            "upper_source_current_bounded": False,
            "grouped_bulk_cell_integrals_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "primitive_gate": {"path": relative(PRIMITIVE_GATE), "sha256": file_hash(PRIMITIVE_GATE)},
            "roster_gate": {"path": relative(ROSTER_GATE), "sha256": file_hash(ROSTER_GATE)},
            "poisson_gate": {"path": relative(POISSON_GATE), "sha256": file_hash(POISSON_GATE)},
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
        "next_obligation": "Estimate only the grouped bulk/Fresnel pieces on fixed-roster cells, carry the exact B-endpoint current into the source T_upper comparison, and use the y=64 telescope before taking absolute values.",
        "proof_boundary": "Exact saved-height endpoint characters and roster reassembly only. No grouped bulk cell integral, lower-interior join, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved.",
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built endpoint roster reassembly: 169 cells, D84(yB)=-84, |D84(64)|<0.612")


if __name__ == "__main__":
    main()
