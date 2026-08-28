#!/usr/bin/env python3
"""Certify the affine-weighted Mordell recursion and two full second steps."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb
from mpmath import mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_affine_weighted_second_step_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PRIOR_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_interval_one_step_current_gate.py"
)
PRIOR_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_interval_one_step_current_gate.json"
)
BOX_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate.json"
)

A = 159_577
L = 2_481_422
K = L - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def arb_upper_text(value: arb, digits: int = 40) -> str:
    midpoint, radius, exponent = value.mid_rad_10exp()
    upper = mp.mpf(abs(int(midpoint)) + int(radius)) * mp.power(10, int(exponent))
    return mp.nstr(upper, digits, min_fixed=-10, max_fixed=10)


def acb_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": str(value.real),
        "imag_ball": str(value.imag),
        "radius_upper": arb_upper_text(value.rad()),
        "absolute_upper": arb_upper_text(value.abs_upper()),
    }


def minimal_period(a: Fraction, tau: Fraction, limit: int = 131_072) -> int:
    for period in range(1, limit + 1):
        if (2 * tau * period).denominator == 1 and (
            a * period + tau * period * period
        ).denominator == 1:
            return period
    raise RuntimeError("affine-weighted theta period not found")


def weighted_period_ball(
    p: Fraction,
    a: Fraction,
    tau: Fraction,
    n: int,
    prior: Any,
) -> tuple[acb, int]:
    """Exact rational-period compression of sum_(k=0)^n (p+k)e(...)."""

    require(n >= 0, "weighted theta length must be nonnegative")
    period = minimal_period(a, tau)
    cycle = [prior.acb_cis_pi(2 * (a * k + tau * k * k)) for k in range(period)]
    z0 = sum(cycle, acb(0))
    z1 = sum((k * cycle[k] for k in range(period)), acb(0))
    cycles, remainder = divmod(n + 1, period)
    r0 = sum(cycle[:remainder], acb(0))
    r1 = sum((k * cycle[k] for k in range(remainder)), acb(0))
    s0 = cycles * z0 + r0
    s1 = period * cycles * (cycles - 1) * z0 / 2 + cycles * z1
    s1 += cycles * period * r0 + r1
    return acb(prior.arb_rational(p)) * s0 + s1, period


def normalize_weighted_parameters(a: Fraction, tau: Fraction, prior: Any) -> dict[str, Any]:
    """Normalize a real quadratic phase to 0<=tau<=1/4 exactly."""

    original_a = a
    original_tau = tau
    operations: list[dict[str, str]] = []

    integer_tau_shift = prior.floor_fraction(tau + Fraction(1, 2))
    if integer_tau_shift:
        tau -= integer_tau_shift
        operations.append({"operation": "integer_tau_period", "shift": str(integer_tau_shift)})

    if tau > Fraction(1, 4):
        tau -= Fraction(1, 2)
        a += Fraction(1, 2)
        operations.append({"operation": "parity_half_shift", "tau_shift": "-1/2", "a_shift": "1/2"})
    elif tau < Fraction(-1, 4):
        tau += Fraction(1, 2)
        a -= Fraction(1, 2)
        operations.append({"operation": "parity_half_shift", "tau_shift": "1/2", "a_shift": "-1/2"})

    conjugated = tau < 0
    if conjugated:
        tau = -tau
        a = -a
        operations.append({"operation": "conjugation", "a_map": "a->-a", "tau_map": "tau->-tau"})

    a_shift = prior.floor_fraction(a + Fraction(1, 2))
    if a_shift:
        a -= a_shift
        operations.append({"operation": "integer_a_period", "shift": str(a_shift)})
    require(Fraction(-1, 2) <= a < Fraction(1, 2), "normalized a escaped")
    require(Fraction(0) <= tau <= Fraction(1, 4), "normalized tau escaped")
    return {
        "original_a": original_a,
        "original_tau": original_tau,
        "a": a,
        "tau": tau,
        "conjugated": conjugated,
        "operations": operations,
    }


def affine_weighted_mordell_step(
    p: Fraction,
    a: Fraction,
    tau: Fraction,
    n: int,
    settings: Any,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    """Apply the exact endpoint-complete Mordell identity to W_n(p;a,tau)."""

    require(tau > 0 and n >= 0, "invalid affine Mordell inputs")
    flint.ctx.prec = settings.precision_bits
    shift, z = prior.nearest_integer_shift(a)
    m = prior.floor_fraction(2 * n * tau)
    w = z / (2 * tau)
    sigma = -Fraction(1, 1) / (4 * tau)
    child_p = 2 * tau * p - z
    child, child_period = weighted_period_ball(child_p, w, sigma, m, prior)
    c = prior.acb_cis_pi(Fraction(1, 4) - z * z / (2 * tau)) / acb(
        2 * prior.arb_rational(tau)
    ).sqrt()
    main = c * child / (2 * prior.arb_rational(tau))

    lower_argument = z - tau + Fraction(1, 2)
    upper_argument = z + (2 * n + 1) * tau - m - Fraction(1, 2)
    lower_h, lower_hz, lower_record = prior.mordell_h_and_derivative_balls(
        lower_argument, -2 * tau, settings
    )
    upper_h, upper_hz, upper_record = prior.mordell_h_and_derivative_balls(
        upper_argument, -2 * tau, settings
    )
    lower_phase = -(z - tau / 2)
    upper_phase = 2 * (Fraction(n) + Fraction(1, 2)) * (
        z + tau * (Fraction(n) + Fraction(1, 2))
    )
    pi_i = acb.pi() * acb(0, 1)
    endpoint = -acb(0, 1) * (
        prior.acb_cis_pi(lower_phase)
        * (acb(prior.arb_rational(2 * p - 1)) * lower_h + lower_hz / pi_i)
        + ((-1) ** m)
        * prior.acb_cis_pi(upper_phase)
        * (acb(prior.arb_rational(2 * p + 2 * n + 1)) * upper_h + upper_hz / pi_i)
    ) / 4
    complete = main + endpoint
    direct, direct_period = weighted_period_ball(p, a, tau, n, prior)
    difference = complete - direct
    require(complete.overlaps(direct) and difference.contains(0), "affine Mordell step mismatch")
    return complete, {
        "p": str(p),
        "a": str(a),
        "tau": str(tau),
        "n": n,
        "nearest_integer_shift": shift,
        "z": str(z),
        "child_m": m,
        "child_a": str(w),
        "child_tau": str(sigma),
        "child_p": str(child_p),
        "child_period": child_period,
        "direct_period": direct_period,
        "active_index_contraction_ratio": str(Fraction(m, max(1, n))),
        "lower_argument": str(lower_argument),
        "upper_argument": str(upper_argument),
        "lower_Mordell": lower_record,
        "upper_Mordell": upper_record,
        "child_weighted_current": acb_record(child),
        "main_current": acb_record(main),
        "joined_endpoint_current": acb_record(endpoint),
        "complete_current": acb_record(complete),
        "direct_weighted_current": acb_record(direct),
        "complete_minus_direct": acb_record(difference),
        "difference_contains_zero": True,
        "passed": True,
    }


def normalized_full_second_step(
    x: Fraction,
    s: Fraction,
    settings: Any,
    prior: Any,
) -> dict[str, Any]:
    """Replace the first transformed current by a normalized second step."""

    first = prior.compute_full_row(x, s, settings)
    r = first["integer_shift_r"]
    n = first["transformed_m"]
    original_a = Fraction(first["transformed_w"])
    original_tau = Fraction(first["transformed_sigma"])
    normalization = normalize_weighted_parameters(original_a, original_tau, prior)
    normalized_a = normalization["a"]
    normalized_tau = normalization["tau"]
    p = Fraction(r)

    direct_original, original_period = weighted_period_ball(p, original_a, original_tau, n, prior)
    direct_normalized, normalized_period = weighted_period_ball(p, normalized_a, normalized_tau, n, prior)
    reconstructed = direct_normalized.conjugate() if normalization["conjugated"] else direct_normalized
    normalization_difference = reconstructed - direct_original
    require(normalization_difference.contains(0), "weighted normalization mismatch")

    if normalized_tau == 0:
        recursive_value = reconstructed
        step_kind = "exact_tau_zero_terminal"
        second_step = None
        child_m = None
    else:
        normalized_complete, second_step = affine_weighted_mordell_step(
            p, normalized_a, normalized_tau, n, settings, prior
        )
        recursive_value = normalized_complete.conjugate() if normalization["conjugated"] else normalized_complete
        step_kind = "contracting_affine_Mordell_step"
        child_m = second_step["child_m"]

    recursive_difference = recursive_value - direct_original
    require(recursive_difference.contains(0), "recursive transformed current mismatch")

    tau1 = x / 2
    c = Fraction(A) + 2 * s
    a1 = x * c / 2
    _, z1 = prior.nearest_integer_shift(a1)
    main_phase = Fraction(1, 4) + (a1 * a1 - z1 * z1) / (2 * tau1)
    multiplier = prior.acb_cis_pi(main_phase) / (
        acb(prior.arb_rational(tau1)) * acb(2 * prior.arb_rational(tau1)).sqrt()
    )
    recursive_main = multiplier * recursive_value
    first_endpoint = acb(
        arb(first["joined_endpoint_current"]["real_ball"]),
        arb(first["joined_endpoint_current"]["imag_ball"]),
    )
    recursive_complete = recursive_main + first_endpoint
    source = acb(
        arb(first["exact_period_source_current"]["real_ball"]),
        arb(first["exact_period_source_current"]["imag_ball"]),
    )
    source_difference = recursive_complete - source
    require(recursive_complete.overlaps(source), "second-step source overlap failed")
    require(source_difference.contains(0), "second-step source difference excludes zero")
    return {
        "x": str(x),
        "s": str(s),
        "first_transformed_n": n,
        "first_transformed_p": str(p),
        "first_transformed_a": str(original_a),
        "first_transformed_tau": str(original_tau),
        "normalization": {
            **normalization,
            "original_a": str(normalization["original_a"]),
            "original_tau": str(normalization["original_tau"]),
            "a": str(normalization["a"]),
            "tau": str(normalization["tau"]),
        },
        "original_period": original_period,
        "normalized_period": normalized_period,
        "normalization_difference": acb_record(normalization_difference),
        "normalization_difference_contains_zero": True,
        "step_kind": step_kind,
        "second_step": second_step,
        "second_child_m": child_m,
        "recursive_transformed_current": acb_record(recursive_value),
        "recursive_minus_direct_transformed": acb_record(recursive_difference),
        "recursive_transformed_difference_contains_zero": True,
        "recursive_first_main": acb_record(recursive_main),
        "inherited_first_endpoint": acb_record(first_endpoint),
        "recursive_complete_source_current": acb_record(recursive_complete),
        "exact_period_source_current": acb_record(source),
        "recursive_complete_minus_source": acb_record(source_difference),
        "recursive_complete_difference_contains_zero": True,
        "passed": True,
    }


SHORT_CASES = (
    {"p": Fraction(3, 2), "a": Fraction(1, 6), "tau": Fraction(1, 4), "n": 7},
    {"p": Fraction(7, 3), "a": Fraction(-1, 5), "tau": Fraction(1, 6), "n": 12},
    {"p": Fraction(-2), "a": Fraction(2, 7), "tau": Fraction(1, 8), "n": 9},
    {"p": Fraction(5), "a": Fraction(7, 6), "tau": Fraction(1, 4), "n": 15},
)


def render_note(artifact: dict[str, Any]) -> str:
    summary = artifact["summary"]
    interior = artifact["full_second_step_rows"][0]
    corner = artifact["full_second_step_rows"][1]
    return f"""# Affine-weighted Mordell recursion and full second steps

Date: 2026-08-24

Status: exact-recursion and pointwise interval-certificate note; the general
affine-weighted Mordell step and two full-roster second-step rows are
certified, but recursive parameter boxes and the small-tau neighbourhood
branch remain open

## Affine-weighted current

Define

```text
W_n(p;a,tau)=sum_(k=0)^n (p+k)exp(2pi i[ak+tau k^2]). (AW1)
```

For `a=z+r`, `-1/2<=z<1/2`, and `tau>0`, put

```text
m=floor(2n tau),  w=z/(2tau),  sigma=-1/(4tau),
C=exp(i*pi/4-i*pi*z^2/(2tau))/sqrt(2tau).
```

Differentiating Kuznetsov's exact theta identity with
`D_a=(2pi i)^(-1)partial_a` gives

```text
W_n(p;a,tau)
 =C/(2tau) W_m(2tau p-z;w,sigma)+E_n(p;a,tau),       (AW2)
```

where the endpoint-complete affine current is

```text
E_n=-i/4{{
 E_-[(2p-1)h(u_-,-2tau)+h_z(u_-,-2tau)/(pi i)]
 +(-1)^m E_+[(2p+2n+1)h(u_+,-2tau)
                    +h_z(u_+,-2tau)/(pi i)]}}.       (AW3)
```

Thus the child affine weight is exactly `p'=2tau p-z`.  The special physical
identity in Section 11.463 is the case `p=a/(2tau)`, for which `p'=r`.

## Exact normalization

For integer `j,l`, the exact symmetries used before each recursive step are

```text
W(p;a+j,tau+l)=W(p;a,tau),
W(p;a+1/2,tau-1/2)=W(p;a,tau),
W(p;a,-tau)=conjugate(W(p;-a,tau)).                  (AW4)
```

The second line follows from `(k-k^2)/2` being integral.  These operations
normalize every point row to `-1/2<=a<1/2`, `0<=tau<=1/4` without changing
`p`.

## Certified rows

Four short rational rows enclose their independent exact-period currents.
For the full `x=2/5,s=1/3` row, the first child

```text
W_992568(31916;-7/6,-5/4)
```

normalizes by conjugation to `W_992568(31916;1/6,1/4)`.  Equation (AW2)
then contracts the active index from `992568` to
`{interior['second_child_m']}` and produces

```text
p'=95747/6,  w'=1/3,  sigma'=-1.                     (AW5)
```

The child is period three.  Its interval second step overlaps the original
period-six transformed current, and after restoring the first-step multiplier
and joined endpoints the complete source difference contains zero.

At `x=1/2,s=0`, the first child normalizes exactly to

```text
W_1240710(39894;-1/2,0),                             (AW6)
```

an alternating affine sum with period two.  The recursion therefore
terminates exactly at the centre; it does not divide by `tau=0`.  A nonzero
box around this centre contains very small positive `tau`, so a rigorous
small-tau neighbourhood branch is still required.

The largest full-row complete/source difference-ball radius is
`{summary['maximum_full_source_difference_radius_upper']}`.  Both differences
contain zero.  An altered 384-bit, `Y=11` replay independently reproduces all
short and full rows.

## Decision

The affine weight is now closed under the exact Mordell step, and the first
interior recursion contracts by one half before reaching a period-three
child.  The corner centre has an exact terminal branch.  What remains is to
lift (AW2)--(AW4) from points to branch-partitioned parameter boxes and add a
uniform small-tau enclosure around the terminal face.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Primary-source boundary: Kuznetsov's exact theta/Mordell identity supplies
the unweighted transform and endpoint functions.  Equations (AW2)--(AW4),
the affine endpoint coefficients, and their recursion are derived here.
Practical Gauss--Laguerre quadrature is not used.

Pi provenance: every `pi` comes from the inherited quadratic Fourier phase,
the differential normalization `D_a`, or the exact Mordell representation.
No fitted circle or polygon constant is introduced.

This gate proves the affine-weighted exact recurrence, four short interval
rows, one contracting full second step, and one exact tau-zero terminal full
row.  It does not prove recursive parameter boxes, a small-tau neighbourhood
branch, physical quadrature, the non-A bound, joined `R_after_A`, `R_Dir`,
`Q_K-T`, an all-height theorem, `Lambda<=0`, PF-infinity, RH, or a prize-level
conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(PRIOR_RESULT.is_file() and BOX_RESULT.is_file(), "Mordell dependencies missing")
    prior_artifact = json.loads(PRIOR_RESULT.read_text(encoding="utf-8"))
    box_artifact = json.loads(BOX_RESULT.read_text(encoding="utf-8"))
    require(prior_artifact.get("passed") is True and box_artifact.get("passed") is True, "dependency not passed")
    prior = load_module("mordell_affine_second_step_dependency", PRIOR_BUILDER)
    settings = prior.PRODUCTION

    short_rows = []
    for case in SHORT_CASES:
        complete, record = affine_weighted_mordell_step(settings=settings, prior=prior, **case)
        require(complete.is_finite(), "nonfinite short affine step")
        short_rows.append(record)

    full_rows = [
        normalized_full_second_step(Fraction(2, 5), Fraction(1, 3), settings, prior),
        normalized_full_second_step(Fraction(1, 2), Fraction(0), settings, prior),
    ]
    require(full_rows[0]["step_kind"] == "contracting_affine_Mordell_step", "interior row did not recurse")
    require(full_rows[0]["second_child_m"] == 496_284, "interior contraction drift")
    require(full_rows[1]["step_kind"] == "exact_tau_zero_terminal", "corner row did not terminate")
    source_difference_radii = [
        arb(row["recursive_complete_minus_source"]["radius_upper"]) for row in full_rows
    ]

    artifact = {
        "kind": STEM,
        "status": "affine_weighted_Mordell_recurrence_four_short_rows_one_contracting_full_second_step_and_one_tau_zero_terminal_certified_recursive_boxes_small_tau_open",
        "passed": True,
        "scope": {"height": 10_000_000_000, "A": A, "L": L, "K": K, "workers": 1},
        "production_settings": settings.__dict__,
        "exact_theorem": {
            "weighted_current": "W_n(p;a,tau)=sum_(k=0)^n (p+k) exp(2*pi*i*(a*k+tau*k^2))",
            "child_weight": "p_prime=2*tau*p-z",
            "main_transform": "C/(2*tau)*W_m(p_prime;z/(2*tau),-1/(4*tau))",
            "lower_endpoint_coefficient": "2*p-1",
            "upper_endpoint_coefficient": "2*p+2*n+1",
        },
        "short_interval_rows": short_rows,
        "full_second_step_rows": full_rows,
        "summary": {
            "short_row_count": len(short_rows),
            "full_row_count": len(full_rows),
            "interior_first_child_n": full_rows[0]["first_transformed_n"],
            "interior_second_child_m": full_rows[0]["second_child_m"],
            "interior_child_period": full_rows[0]["second_step"]["child_period"],
            "corner_terminal_period": full_rows[1]["normalized_period"],
            "all_full_source_differences_contain_zero": all(
                row["recursive_complete_difference_contains_zero"] for row in full_rows
            ),
            "maximum_full_source_difference_radius_upper": arb_upper_text(max(source_difference_radii)),
        },
        "decision": {
            "general_affine_weight_closed_under_Mordell_step": True,
            "exact_phase_normalization_proved": True,
            "interior_full_second_step_contracts": True,
            "corner_center_has_exact_tau_zero_terminal": True,
            "two_full_recursive_source_reassemblies_certified": True,
            "recursive_parameter_boxes_built": False,
            "small_tau_neighbourhood_branch_built": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Lift the affine recursion and normalization to branch-partitioned parameter boxes. "
            "At the x=2/5 normalization wall, hull the parity/conjugation and floor branches; near "
            "x=1/2, derive a uniform small-tau weighted-sum branch that joins the exact tau-zero terminal."
        ),
        "proof_boundary": (
            "Exact affine-weighted Mordell identity, exact point normalization, four short interval rows, "
            "one contracting full pointwise second step, and one exact tau-zero terminal row only. No "
            "recursive parameter boxes, small-tau neighbourhood branch, physical quadrature, non-A bound, "
            "joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level "
            "conclusion is proved."
        ),
        "primary_source": prior_artifact["primary_source"],
        "dependencies": {
            "prior_result": {"path": relative(PRIOR_RESULT), "sha256": file_hash(PRIOR_RESULT)},
            "prior_builder": {"path": relative(PRIOR_BUILDER), "sha256": file_hash(PRIOR_BUILDER)},
            "parameter_box_result": {"path": relative(BOX_RESULT), "sha256": file_hash(BOX_RESULT)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "certified affine-weighted Mordell recurrence, one contracting full second step, and one tau-zero terminal",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
