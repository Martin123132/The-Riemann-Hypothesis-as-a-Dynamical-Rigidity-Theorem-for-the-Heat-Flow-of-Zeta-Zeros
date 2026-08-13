#!/usr/bin/env python3
"""Certify the analytic pole-subtracted cotangent current on the B window."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_current_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "step_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate.json",
    "B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
WINDOW_XI = 70
LOCAL_MODES = (621, 622)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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
    c, mode = sp.symbols("c m", positive=True, real=True)
    summand = 1 / (4 * (c**2 - mode**2))
    derivative = sp.factor(sp.diff(summand, c))
    require(
        sp.simplify(derivative + c / (2 * (c**2 - mode**2) ** 2)) == 0,
        "summand derivative failed",
    )

    k, j = sp.symbols("k j", positive=True, integer=True)
    # Laurent-expand the classical cotangent sum with its kth pole removed.
    delta = sp.symbols("delta", real=True)
    full = sp.pi * sp.cot(sp.pi * (k + delta)) / (8 * (k + delta)) - 1 / (8 * (k + delta) ** 2)
    pole = 1 / (4 * ((k + delta) ** 2 - k**2))
    removed_limit = sp.simplify(sp.limit(full - pole, delta, 0))
    require(removed_limit == -sp.Rational(3, 16) / k**2, "one-pole removable limit failed")
    two_removed = sp.simplify(removed_limit - 1 / (4 * (k**2 - j**2)))

    return {
        "all_pair_sum": "S_all(c)=sum_(m>=1)1/[4(c^2-m^2)]=pi*cot(pi*c)/(8c)-1/(8c^2)",
        "local_subtraction": "S_hat(c)=S_all(c)-1/[4(c^2-621^2)]-1/[4(c^2-622^2)]",
        "summand_derivative": "d/dc {1/[4(c^2-m^2)]}=-c/[2(c^2-m^2)^2]<0",
        "strict_monotonicity": "S_hat is strictly decreasing on every interval containing no unremoved integer pole",
        "one_pole_limit": "lim_(c->k){S_all(c)-1/[4(c^2-k^2)]}=-3/(16k^2)",
        "two_pole_limit": f"lim_(c->k)S_hat(c)={sp.sstr(two_removed)} with (k,j)=(621,622) or (622,621)",
        "paired_boundary_current": "C_hat(x)=-2(2+i*pi*B^2*x)S_hat(Bx/2)/pi^2",
    }


def cotangent_sum(c: arb, pi: arb) -> arb:
    value = pi * (pi * c).cot() / (8 * c) - 1 / (8 * c**2)
    for mode in LOCAL_MODES:
        value -= 1 / (4 * (c**2 - arb(mode) ** 2))
    return value


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t = arb(T)
    endpoint = arb(B)
    pi = arb.pi()
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_h = hessian.sqrt()
    spacing = 2 * root_h / endpoint
    c0 = endpoint * x0 / 2
    c_low = endpoint * (x0 - arb(WINDOW_XI) / root_h) / 2
    c_high = endpoint * (x0 + arb(WINDOW_XI) / root_h) / 2

    s_low = cotangent_sum(c_low, pi)
    s_center = cotangent_sum(c0, pi)
    s_high = cotangent_sum(c_high, pi)
    require(s_low > 0, "left cotangent background is not positive")
    require(s_center < 0, "central cotangent background is not negative")
    require(s_high < 0, "right cotangent background is not negative")
    require(max(abs(s_low), abs(s_high)) < arb("0.000288"), "cotangent background exceeds 0.000288")

    root_low = arb("621.49827876663")
    root_high = arb("621.49827876665")
    root_low_value = cotangent_sum(root_low, pi)
    root_high_value = cotangent_sum(root_high, pi)
    require(root_low_value > 0 and root_high_value < 0, "cotangent root bracket failed")
    xi_root_low = spacing * (root_low - c0)
    xi_root_high = spacing * (root_high - c0)
    require(-arb("6.568") < xi_root_low < xi_root_high < -arb("6.566"), "standardized root bracket drift")

    def removable_value(mode: int, other: int) -> arb:
        k = arb(mode)
        j = arb(other)
        return (-3 / (16 * k**2)) - 1 / (4 * (k**2 - j**2))

    removable_621 = removable_value(621, 622)
    removable_622 = removable_value(622, 621)
    require(removable_621 > 0 and removable_622 < 0, "removable-value signs failed")

    center_current_abs = (
        2
        * (4 + (pi * endpoint**2 * x0) ** 2).sqrt()
        * abs(s_center)
        / pi**2
    )
    one_xi_l1_lower = center_current_abs / root_h
    require(center_current_abs > arb("88000"), "central analytic current below 88000")
    require(one_xi_l1_lower > arb("0.000302"), "one-xi triangle lower bound below 0.000302")

    return {
        "height": T,
        "endpoint": B,
        "window_xi": WINDOW_XI,
        "c_window_low_ball": c_low.str(PRECISION, more=True),
        "c_trace_ball": c0.str(PRECISION, more=True),
        "c_window_high_ball": c_high.str(PRECISION, more=True),
        "S_hat_left_ball": s_low.str(PRECISION, more=True),
        "S_hat_trace_ball": s_center.str(PRECISION, more=True),
        "S_hat_right_ball": s_high.str(PRECISION, more=True),
        "S_hat_uniform_absolute_upper_bound": "0.000288",
        "root_c_bracket": ["621.49827876663", "621.49827876665"],
        "root_endpoint_sign_balls": [root_low_value.str(PRECISION, more=True), root_high_value.str(PRECISION, more=True)],
        "root_xi_bracket_balls": [xi_root_low.str(PRECISION, more=True), xi_root_high.str(PRECISION, more=True)],
        "removable_value_at_621_ball": removable_621.str(PRECISION, more=True),
        "removable_value_at_622_ball": removable_622.str(PRECISION, more=True),
        "central_analytic_current_absolute_ball": center_current_abs.str(PRECISION, more=True),
        "absolute_integral_lower_bound_on_xi_0_to_1_ball": one_xi_l1_lower.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Pole-subtracted cotangent current on the B window

Date: 2026-08-13

Status: exact meromorphic reduction plus interval certificate; not a proof of
the complete B-face estimate

For `c=Bx/2`, the paired first boundary current from Section 11.362 contains

```text
S_all(c)=sum_(m>=1)1/[4(c^2-m^2)]
        =pi*cot(pi*c)/(8c)-1/(8c^2).                 (PC1)
```

Remove the two exact local modes before estimating:

```text
S_hat(c)=S_all(c)-1/[4(c^2-621^2)]
                  -1/[4(c^2-622^2)].                (PC2)
```

The apparent singularities at 621 and 622 are removable.  In general,

```text
lim_(c->k){{S_all(c)-1/[4(c^2-k^2)]}}=-3/(16k^2),    (PC3)
```

so the second local term can simply be subtracted at each limit.  The
certified values have opposite signs:

```text
S_hat(621)={c['removable_value_at_621_ball']},
S_hat(622)={c['removable_value_at_622_ball']}.         (PC4)
```

Every retained summand satisfies

```text
d/dc {{1/[4(c^2-m^2)]}}=-c/[2(c^2-m^2)^2]<0.         (PC5)
```

Thus `S_hat` is strictly decreasing on the entire `|xi|<=70` B window.  Its
endpoint and trace values are

```text
S_hat(c_low)={c['S_hat_left_ball']},
S_hat(c_B)  ={c['S_hat_trace_ball']},
S_hat(c_high)={c['S_hat_right_ball']},                (PC6)
```

and `|S_hat|<0.000288` throughout.  There is exactly one zero, with

```text
c_* in [{c['root_c_bracket'][0]}, {c['root_c_bracket'][1]}],
xi_* in {c['root_xi_bracket_balls']}.                 (PC7)
```

The corresponding analytic first boundary current is

```text
C_hat(x)=-2(2+i*pi*B^2*x)S_hat(Bx/2)/pi^2.           (PC8)
```

Its absolute value at the trace saddle exceeds `88000`.  Since `S_hat` is
negative and decreasing for `0<=xi<=1`, even that one-unit strip has raw
absolute integral above

```text
{c['absolute_integral_lower_bound_on_xi_0_to_1_ball']} > 0.000302. (PC9)
```

Equation (PC9) is a cancellation guard, not a lower bound for the signed
residual.  It proves that the analytic cotangent background cannot be
triangled away: its oscillatory integral must be combined with the exact
621/622 currents and the endpoint/negative/outer completion.

Pi provenance: `pi` comes from the equation-(9) face current and the integer
cotangent partial fraction.  No fitted constant is used.

Proof boundary: exact pole removal, monotonicity, unique-zero enclosure, and
saved-height size bounds for the first analytic boundary current only.  The
exact local currents, the other rational terms, outside-window estimate, and
A-fold splice remain open.  This result does not establish the complete
paired residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or a prize-level
conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["step_tail"]["decision"]["nonlocal_leading_pair_sum_has_cotangent_form"] is True,
        "cotangent dependency drift",
    )
    require(
        dependencies["B_window"]["decision"]["xi_70_window_contains_exactly_B_crossings_621_622"] is True,
        "B-window dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "B_face_first_paired_boundary_current_poles_removed_and_analytic_background_certified",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "local_621_622_poles_removed_exactly": True,
            "pole_subtracted_cotangent_background_analytic_on_xi_70_window": True,
            "pole_subtracted_background_strictly_decreasing": True,
            "pole_subtracted_background_has_one_certified_window_zero": True,
            "raw_absolute_value_bound_is_cancellation_compatible": False,
            "exact_local_modes_621_622_enclosed": False,
            "complete_B_face_estimate_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Construct the exact continuous 621/622 two-mode replacement in the same face normalization, combine it with C_hat before integration, then extend the pole subtraction to the next two rational Fresnel currents. Do not triangle C_hat or separate endpoint continuations.",
        "proof_boundary": "Exact pole-subtracted first cotangent current only. No exact local two-mode enclosure, higher rational-current sum, outside-window estimate, A-fold splice, complete paired residual, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified pole-subtracted B cotangent current and unique window zero", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
