#!/usr/bin/env python3
"""Select the exact extracted Hankel carrier at every lower-fold event cell."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_hankel_carrier_event_selection_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "exact_Hankel_split": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_exact_hankel_branch_factorization_gate.json",
    "all_event_height_cells": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_continuous_height_gate.json",
    "turning_event_geometry": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_fresnel_retention_gate.json",
}

C = 159_577
Y = 64
MAX_ODD = 797
PRECISION = 100


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


def symbolic_certificate() -> dict[str, str]:
    beta, d, delta, y, sigma = sp.symbols("beta d delta y sigma", positive=True, real=True)
    lam = d**2 - delta / beta
    selected = sigma * (sp.sqrt(lam) - d)
    rationalized = -sigma * delta / (beta * (sp.sqrt(lam) + d))
    require(sp.simplify(selected - rationalized) == 0, "selected endpoint slope rationalization failed")

    x = sp.symbols("X", positive=True)
    xi = sp.Rational(2, 3) * x ** sp.Rational(3, 2)
    require(sp.diff(xi, x) == sp.sqrt(x), "cubic carrier derivative drift")
    return {
        "event_height": "t=tau_m+delta; lambda=d_m^2-delta/beta",
        "branch_label": "sigma=sign(d_m); C_sigma is H^(1) for sigma=+1 and H^(2) for sigma=-1",
        "selected_endpoint_slope": "F_sigma(0)=-d_m+sigma*sqrt(lambda)=-sigma*delta/[beta*(sqrt(lambda)+abs(d_m))]",
        "opposite_endpoint_slope": "F_-sigma(0)=-sigma*(abs(d_m)+sqrt(lambda))",
        "carrier_derivatives": "F_+ = y/(2beta)-d_m+sqrt(lambda+y); F_- = y/(2beta)-d_m-sqrt(lambda+y)",
        "extracted_Hankel_amplitude": "R_nu^+(xi)=sqrt(pi*xi/2)*exp(-i*(xi-pi*nu/2-pi/4))*H_nu^(1)(xi); R_nu^-=conjugate(R_nu^+)",
        "exact_branch_form": "C_+=exp(i*(xi-pi/4))/(2sqrt(pi))*[U*X^(-1/4)R_(1/3)^+ + i*V*X^(1/4)R_(2/3)^+]; C_-=conjugate(C_+)",
    }


def interval_certificate() -> dict[str, str]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    beta = (pi * arb(C) ** 2 / 8) ** (arb(1) / 3)
    half_cell_lambda_shift = pi / (16 * beta)
    d_min = beta / C
    d_max = arb(MAX_ODD) * beta / C
    lambda_min = d_min**2 - half_cell_lambda_shift
    lambda_min_exact = pi / (16 * beta)
    require(lambda_min.overlaps(lambda_min_exact), "minimum lambda identity drift")
    require(lambda_min > arb("9e-5"), "event-cell positive-X floor failed")

    x_max = d_max**2 + half_cell_lambda_shift + Y
    require(x_max < arb(180), "fold-strip X ceiling failed")
    require(beta**2 > arb("4.6e6"), "beta-square floor failed")
    require(x_max < beta**2, "carrier monotonicity inequality failed")

    opposite_endpoint_gap = d_min + lambda_min.sqrt()
    selected_face_slope = half_cell_lambda_shift / (d_min + lambda_min.sqrt())
    require(opposite_endpoint_gap > arb("0.023"), "opposite-branch endpoint gap failed")
    require(selected_face_slope < arb("0.004"), "selected branch face slope drift")

    # d=(4m-C)beta/C and 4m-C ranges over the 399 nonzero odd labels.
    labels = [(-1) ** (n + 1) * (2 * n + 1) for n in range(399)]
    require(len(labels) == 399 and len(set(labels)) == 399, "event label roster drift")
    require(min(abs(value) for value in labels) == 1, "minimum odd event label drift")
    require(max(abs(value) for value in labels) == MAX_ODD, "maximum odd event label drift")

    return {
        "beta_ball": beta.str(PRECISION, more=True),
        "d_min_ball": d_min.str(PRECISION, more=True),
        "d_max_ball": d_max.str(PRECISION, more=True),
        "lambda_min_ball": lambda_min.str(PRECISION, more=True),
        "lambda_min_exact_formula": "pi/(16*beta)",
        "fold_strip_X_max_ball": x_max.str(PRECISION, more=True),
        "beta_squared_ball": (beta**2).str(PRECISION, more=True),
        "opposite_branch_endpoint_gap_lower_ball": opposite_endpoint_gap.str(PRECISION, more=True),
        "selected_branch_face_slope_upper_ball": selected_face_slope.str(PRECISION, more=True),
        "event_count": "399",
        "event_abs_odd_label_range": "1..797",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Beta^-4 Hankel carrier: exact event selection

Date: 2026-08-13
Status: exact carrier-selection gate; not a proof of the amplitude splice

Factor the positive-`X` Hankel branches from Section 11.378 as

```text
C_+(X)=e^(i[xi-pi/4]) B_+(X),
C_-(X)=e^(-i[xi-pi/4])B_-(X),       xi=2X^(3/2)/3,    (ES1)
```

where `B_+` contains the exact scaled Hankel amplitudes of orders `1/3` and
`2/3`, and `B_-=conj(B_+)` for real data.  This factorization is exact: no
large-`X` replacement is made.  Restoring the normal carrier
`exp(i[y^2/(4beta)-d_m*y])` gives extracted phase derivatives

```text
F_+(y)=y/(2beta)-d_m+sqrt(lambda+y),
F_-(y)=y/(2beta)-d_m-sqrt(lambda+y).                  (ES2)
```

At event `tau_m`, write `t=tau_m+delta`.  The exact event identity is

```text
lambda=d_m^2-delta/beta.                              (ES3)
```

Set `sigma=sign(d_m)`, choosing `C_+` for positive `d_m` and `C_-` for
negative `d_m`.  Rationalization gives the exact endpoint orientation

```text
F_sigma(0)
 =-d_m+sigma*sqrt(lambda)
 =-sigma*delta/[beta(sqrt(lambda)+|d_m|)],             (ES4)

F_-sigma(0)=-sigma(|d_m|+sqrt(lambda)).                (ES5)
```

Thus the selected carrier crosses the lower endpoint exactly once, at the
event center, while the opposite carrier does not approach that endpoint
saddle.

For the 399 cells, `|4m-C|` runs through the nonzero odd labels of absolute
size `1,...,797`.  Since `beta^3=pi*C^2/8` and `|delta|<=pi/16`,

```text
lambda >= (beta/C)^2-pi/(16beta)
        = pi/(16beta)
        = {c['lambda_min_ball']} >9e-5.                (ES6)
```

Hence every point of every certified event cell and every `0<=y<=64` lies
strictly in `X=lambda+y>0`.  The entire fold strip also obeys

```text
X<{c['fold_strip_X_max_ball']}<180<beta^2.             (ES7)
```

Consequently the opposite extracted carrier has fixed derivative sign on
the fold strip, with endpoint magnitude greater than
`{c['opposite_branch_endpoint_gap_lower_ball']}`.  The selected carrier is
strictly monotone in `y` and has at most one carrier saddle.  Its endpoint
slope at either cell face has magnitude below
`{c['selected_branch_face_slope_upper_ball']}`.

Equations (ES1)--(ES7) select the correct cubic carrier and orient every
event cell.  They do not bound derivatives of the exact scaled Hankel
amplitudes `B_+/-`; therefore they do not yet license integration by parts on
the opposite branch or identify the selected branch with the exact
logistic/Gamma carrier.  No corridor remainder, `Q_K-T`, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level result is
proved.
"""


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    dependencies = {}
    for name, path in DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {name}")
        dependencies[name] = json.loads(path.read_text(encoding="utf-8"))
    require(dependencies["exact_Hankel_split"]["decision"]["exact_Hankel_branch_factorization_proved"] is True, "Hankel split drift")
    require(dependencies["all_event_height_cells"]["claims"]["all_399_exact_event_cells_uniformly_certified"] is True, "event cell atlas drift")
    require(dependencies["turning_event_geometry"]["decision"]["exact_Airy_event_Fresnel_coordinate_join_proved"] is True, "turning-event identity drift")

    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_Hankel_carrier_selection_and_all_event_positive_X_floor_proved",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "all_399_event_cells_strictly_positive_X": True,
            "selected_extracted_carrier_crosses_at_event_center": True,
            "opposite_extracted_carrier_nonstationary_on_fold_strip": True,
            "selected_extracted_carrier_has_at_most_one_fold_strip_saddle": True,
            "exact_scaled_Hankel_amplitude_derivative_bound_proved": False,
            "ordinary_logistic_Gamma_carrier_identification_proved": False,
            "corridor_remainder_bound_proved": False,
            "complete_Q_K_minus_T_bound_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "next_action": "Bound the exact scaled Hankel amplitude and its first derivatives on the positive-X cell range, then compare the selected branch with the endpoint-retaining logistic/Gamma carrier and integrate the opposite branch nonstationarily.",
        "proof_boundary": "Exact extracted-carrier selection and positive-X event-cell geometry only. No scaled-Hankel amplitude derivative bound, quantitative ordinary-carrier identification, corridor splice, Q_K-T bound, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
        "dependencies": {name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "resource_policy": {"workers": 1, "process_priority": priority},
        "runtime": {"elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    artifact["sources"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("proved exact Hankel carrier event selection: 399 positive-X cells, opposite carrier gap >0.023", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
