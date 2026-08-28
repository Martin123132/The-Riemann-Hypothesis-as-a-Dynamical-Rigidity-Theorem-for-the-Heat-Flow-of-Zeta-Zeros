#!/usr/bin/env python3
"""Certify the finite quarter-disk vertical/arc/tail representation."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from typing import Callable


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_quarter_disk_vertical_arc_tail_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
ARB_HELPER = BUILDER.with_name(STEM.removesuffix("_gate") + "_arb.py")
DEPENDENCIES = {
    "mellin_hurwitz": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_gamma_balanced_mellin_hurwitz_finite_difference_gate.json",
    "pole_safe_join": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_PW_unowned_cell_pole_safe_two_boundary_contour_join_gate.json",
}

HEIGHT = 10_000_000_000
Y_UPPER = "39936.5"
Q_MINUS = "79788.5"
LABEL_COUNT = 2_481_423


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
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def fmt(value: mp.mpf | mp.mpc, digits: int = 45) -> str:
    return mp.nstr(value, digits)


def delta_q(u: mp.mpf | mp.mpc, q_minus: mp.mpf, count: int) -> mp.mpc:
    if u == 0:
        return mp.mpc(count)
    return mp.exp(-q_minus * u) * mp.expm1(-count * u) / mp.expm1(-u)


def roster(y: mp.mpf, q_minus: mp.mpf, count: int) -> mp.mpc:
    return mp.fsum(mp.exp(2j * mp.pi * (q_minus + index) * y) for index in range(count))


def scaled_prefactor(s: mp.mpc) -> mp.mpc:
    return mp.exp(mp.pi * mp.im(s) / 2 + 0.25j * mp.pi) * mp.power(2 * mp.pi, s - 1)


def integrate_zero_to_bound(
    s: mp.mpc,
    rest: Callable[[mp.mpf], mp.mpc],
    bound: mp.mpf,
    log_cuts: tuple[mp.mpf, ...],
    direct_cuts: tuple[mp.mpf, ...],
) -> mp.mpc:
    require(bound > 0, "integration bound must be positive")
    split = min(bound, mp.mpf(1))
    v0 = -mp.log(split)

    def logarithmic(v: mp.mpf) -> mp.mpc:
        if not mp.isfinite(v):
            return mp.mpc(0)
        y = mp.exp(-v)
        return mp.exp(-(1 - s) * v) * rest(y)

    shifted = tuple(v0 + cut for cut in log_cuts)
    low = mp.quad(logarithmic, shifted[:-1])
    low += mp.quadosc(logarithmic, [shifted[-2], mp.inf], omega=abs(mp.im(s)))
    if bound <= 1:
        return low
    require(direct_cuts[0] == 1 and direct_cuts[-1] == bound, "direct cuts do not span [1,bound]")
    return low + mp.quad(lambda y: mp.power(y, -s) * rest(y), direct_cuts)


def altered_contour_replay() -> dict[str, str]:
    mp.mp.dps = 68
    s = mp.mpc("0.5", "3.25")
    q_minus = mp.mpf("2.5")
    count = 5
    y_upper = mp.mpf("3.5")
    radius = 2 * mp.pi * y_upper
    prefactor = scaled_prefactor(s)
    log_cuts = tuple(mp.mpf(text) for text in ("0", "1", "3", "7", "15", "30", "+inf"))
    direct_cuts = (
        mp.mpf(1),
        mp.mpf(2),
        mp.mpf(4),
        mp.mpf(8),
        mp.mpf(12),
        mp.mpf(16),
        radius,
    )

    finite_mellin = integrate_zero_to_bound(
        s,
        lambda u: mp.exp(1j * u * u / (4 * mp.pi)) * delta_q(u, q_minus, count),
        radius,
        log_cuts,
        direct_cuts,
    )
    positive_tail = mp.quad(
        lambda u: mp.power(u, -s)
        * mp.exp(1j * u * u / (4 * mp.pi))
        * delta_q(u, q_minus, count),
        (radius, mp.mpf(28), mp.mpf(36), mp.mpf(48), mp.mpf(64), mp.mpf(90)),
    )
    tail_cut = mp.mpf(90)
    omitted_tail_bound = (
        abs(prefactor)
        * mp.power(tail_cut, mp.mpf("-0.5"))
        * mp.exp(-q_minus * tail_cut)
        / (q_minus * (1 - mp.exp(-tail_cut)))
    )

    vertical = integrate_zero_to_bound(
        s,
        lambda y: mp.exp(-1j * mp.pi * y * y) * roster(y, q_minus, count),
        y_upper,
        log_cuts,
        (mp.mpf(1), mp.mpf(2), y_upper),
    )

    u0 = -1j * radius
    delta0 = delta_q(u0, q_minus, count)
    endpoint = (
        prefactor
        * mp.power(u0, -s)
        * mp.exp(1j * u0 * u0 / (4 * mp.pi))
        * delta0
        * radius
    )

    def normalized_arc(delta: mp.mpf) -> mp.mpc:
        w_exponent = 1j * radius * mp.expm1(1j * delta)
        w = mp.exp(w_exponent)
        exponent = (
            mp.im(s) * delta
            + 0.5j * delta
            - 1j * mp.pi * y_upper * y_upper * mp.expm1(2j * delta)
            + 1j * q_minus * radius * mp.expm1(1j * delta)
        )
        quotient = (1 + mp.power(w, count)) / (1 + w)
        return mp.exp(exponent) * quotient

    arc = endpoint * mp.quad(
        normalized_arc,
        (
            mp.mpf(0),
            mp.mpf("1e-8"),
            mp.mpf("1e-6"),
            mp.mpf("1e-4"),
            mp.mpf("0.001"),
            mp.mpf("0.01"),
            mp.mpf("0.1"),
            mp.mpf("0.5"),
            mp.pi / 2,
        ),
    )

    finite_discrepancy = abs(prefactor * finite_mellin - vertical - arc)
    full_discrepancy = abs(
        prefactor * (finite_mellin + positive_tail)
        - vertical
        - arc
        - prefactor * positive_tail
    )
    endpoint_scale_discrepancy = abs(abs(endpoint) - mp.sqrt(y_upper))
    require(finite_discrepancy < mp.mpf("1e-34"), "altered finite contour identity failed")
    require(full_discrepancy < mp.mpf("1e-34"), "altered full contour identity failed")
    require(endpoint_scale_discrepancy < mp.mpf("1e-54"), "endpoint scale identity failed")
    require(omitted_tail_bound < mp.mpf("1e-70"), "altered omitted positive tail is too large")
    return {
        "s": fmt(s),
        "q_minus": fmt(q_minus),
        "count": str(count),
        "Y": fmt(y_upper),
        "finite_mellin": fmt(finite_mellin),
        "vertical_current": fmt(vertical),
        "compact_arc": fmt(arc),
        "scaled_positive_tail": fmt(prefactor * positive_tail),
        "scaled_positive_tail_omission_bound": fmt(omitted_tail_bound),
        "finite_identity_discrepancy": fmt(finite_discrepancy),
        "full_identity_discrepancy": fmt(full_discrepancy),
        "endpoint_scale_discrepancy": fmt(endpoint_scale_discrepancy),
    }


def run_arb_helper(variant: str) -> dict[str, object]:
    command = ["py", "-3.13", str(ARB_HELPER), "--variant", variant]
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=1800,
        check=False,
    )
    require(completed.returncode == 0, f"Arb helper failed: {completed.stderr[-2000:]}")
    payload = json.loads(completed.stdout)
    require(payload.get("passed") is True, "Arb helper did not pass")
    return payload


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    for path in (CHECKER, ARB_HELPER):
        require(path.is_file(), f"missing gate source: {path.name}")

    dependencies: dict[str, dict[str, str]] = {}
    for name, path in DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {name}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        require(payload.get("passed") is True, f"dependency did not pass: {name}")
        dependencies[name] = {"path": relative(path), "sha256": file_hash(path)}

    altered = altered_contour_replay()
    arb_certificate = run_arb_helper("production")
    require(arb_certificate["rotated_arc_absolute"]["upper_float"] < 0.057665, "arc guard drift")
    require(
        arb_certificate["error_budgets"]["positive_real_tail_log10_upper"]["upper_float"]
        < -100_000_000,
        "positive tail guard drift",
    )

    artifact = {
        "kind": "rh_c1_hardy_equation4_joined_packet_quarter_disk_vertical_arc_tail_gate",
        "date": "2026-08-27",
        "status": "exact_quarter_disk_identity_and_production_upper_arc_arb_enclosure_certified_lower_cell_and_complete_DK_open",
        "passed": True,
        "production": {
            "t": str(HEIGHT),
            "s": "1/2+i*t",
            "Y": Y_UPPER,
            "R": "2*pi*Y",
            "q_minus": Q_MINUS,
            "q_plus": "2561211.5",
            "label_count": LABEL_COUNT,
        },
        "definitions": {
            "F_s": "F_s(u)=u^(-s)*exp(i*u^2/(4*pi))*Delta_q(u)",
            "M_s": "M_s=integral_0^infinity F_s(u)du",
            "vertical_current": "V_Y=integral_0^Y y^(-s)*exp(-i*pi*y^2)*D_W(y)dy",
            "clockwise_arc": "A_Y^cw=integral_[R to -iR clockwise] F_s(u)du",
            "positive_tail": "T_Y=integral_R^infinity F_s(u)du",
        },
        "exact_identities": {
            "orientation": "M_s=integral_[0 to -iR]F_s(u)du-A_Y^cw+T_Y",
            "pointwise_vertical_map": (
                "E_s*F_s(-2*pi*i*y)*(-2*pi*i*dy)="
                "y^(-s)*exp(-i*pi*y^2)*D_W(y)dy"
            ),
            "scaled_quarter_disk": "S_W=V_Y+C_Y+E_s*T_Y, C_Y=-E_s*A_Y^cw",
            "unowned_join": "U_unowned=integral_0^L g_W(y)dy+C_U+E_s*T_U",
            "endpoint_carrier": (
                "C_U=exp(i*(pi-t*log(U)-pi*U^2))*C_U^rot for U=39936+1/2"
            ),
            "rotated_arc_integrand": (
                "C_U^rot=sqrt(U)*integral_0^(pi/2) exp(X(delta))"
                "*(1+w(delta)^N)/(1+w(delta)) ddelta"
            ),
            "stable_X": arb_certificate["stable_endpoint_exponent"],
            "stable_w": arb_certificate["stable_geometric_coordinate"],
        },
        "altered_replay": altered,
        "arb_certificate": arb_certificate,
        "route_decision": {
            "absolute_arc_scale": "the signed upper arc has certified modulus below 0.057665",
            "positive_tail": "the positive-real tail is exponentially negligible at production height",
            "comparison_only": (
                "the arc scale lies inside the numerical width of the saved D_K corridor, but the arc is not D_K"
            ),
            "next_action": (
                "derive a cancellation-preserving lower-cell evaluator and couple its signed value with C_U, "
                "the ordinary Gamma-subtracted packet, and the grouped natural A lift before taking a norm"
            ),
        },
        "dependencies": dependencies,
        "source_hashes": {
            "builder": file_hash(BUILDER),
            "checker": file_hash(CHECKER),
            "arb_helper": file_hash(ARB_HELPER),
        },
        "validation": {
            "resource_mode": priority,
            "elapsed_seconds": time.perf_counter() - started,
            "independent_checker_required": True,
        },
        "scope": {
            "certified": (
                "Exact finite quarter-disk orientation and vertical-current map, the cancellation-preserving "
                "unowned-cell representation, and a rigorous Arb enclosure of the complete production upper arc "
                "with explicit polynomial, geometric-numerator, compact-tail, and positive-real-tail errors."
            ),
            "not_certified": (
                "No lower-cell enclosure, complete joined packet, J_Z or D_K enclosure, non-A theorem, all-height "
                "result, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
            ),
        },
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    arc_abs = arb_certificate["rotated_arc_absolute"]
    action_max = arb_certificate["action_maximum"]
    action_delta = arb_certificate["action_stationary_delta"]
    tail_log = arb_certificate["error_budgets"]["positive_real_tail_log10_upper"]
    NOTE.write_text(
        f"""# Quarter-disk vertical/arc/tail representation for the joined packet

Date: 2026-08-27

Status: exact contour identity and rigorous production upper-arc enclosure certified; lower-cell and complete `D_K` assembly open.

Let

```text
F_s(u)=u^(-s)exp(i*u^2/(4*pi))Delta_q(u),
M_s=integral_0^infinity F_s(u)du,
E_s=exp(pi*t/2+i*pi/4)(2*pi)^(s-1),
R=2*pi*Y.
```

Close the fourth-quadrant quarter disk with the principal branch.  The small
circle at zero vanishes because `Re(1-s)=1/2`.  If `A_Y^cw` is the clockwise
arc from `R` to `-iR` and `T_Y` is the positive-real tail from `R`, contour
orientation gives

```text
M_s=integral_[0 to -iR]F_s(u)du-A_Y^cw+T_Y.         (QD1)
```

On `u=-2*pi*i*y`, with `0<y<Y`, the principal logarithm gives exactly

```text
E_s u^(-s)(-2*pi*i*dy)=y^(-s)dy,
exp(i*u^2/(4*pi))=exp(-i*pi*y^2),
Delta_q(-2*pi*i*y)=D_W(y).                          (QD2)
```

Thus, writing `C_Y=-E_s A_Y^cw`,

```text
S_W=V_Y+C_Y+E_s T_Y,
V_Y=integral_0^Y y^(-s)exp(-i*pi*y^2)D_W(y)dy.      (QD3)
```

Combining (QD3) with the already-certified two-boundary join at
`L=621.5`, `U=39936.5` cancels the whole interior current before any norm:

```text
U_unowned=integral_0^L g_W(y)dy+C_U+E_s T_U.        (QD4)
```

For the production upper boundary put `u=-2*pi*i*U*exp(i*delta)`.  Since `U`
and `q_-` are half integers and the label count `N=2481423` is odd, the
endpoint carrier separates exactly:

```text
C_U=exp(i*[pi-t*log(U)-pi*U^2]) C_U^rot,

C_U^rot=sqrt(U) integral_0^(pi/2)
  exp(X(delta)) (1+w(delta)^N)/(1+w(delta)) ddelta. (QD5)
```

The stable coordinates used by Arb are

```text
X(delta)=(t+2*pi*U^2-q_-R+i/2)delta
 -i*pi*U^2[expm1(2i*delta)-2i*delta]
 +i*q_-R[expm1(i*delta)-i*delta],

w(delta)=exp(i*R*expm1(i*delta)).                   (QD6)
```

They expose the small real endpoint slope
`206.8358840647606...` instead of subtracting three unresolved
`10^10`-scale interval terms.

Arb proves that the real action has its unique maximum at

```text
delta_* in {action_delta['ball']},
A(delta_*) in {action_max['ball']}.
```

The complete signed arc, including the tiny initial geometric-numerator
layer, the finite Taylor-coordinate error, and the omitted compact tail,
satisfies

```text
C_U^rot in
  {arb_certificate['rotated_arc']['real_ball']}
 +i*{arb_certificate['rotated_arc']['imag_ball']},

|C_U|=|C_U^rot| in {arc_abs['ball']} < 0.057665.    (QD7)
```

The separate positive-real tail has

```text
log10 |E_s T_U| <= {tail_log['ball']},              (QD8)
```

so it is negligible on every displayed corridor scale.  An altered
`t=3.25`, five-label contour replay checks (QD1)--(QD3) directly.  The
independent checker changes the altered height, roster, endpoint, precision,
Taylor degree, initial split, and Arb paneling, and requires rigorous overlap
with (QD7).

The value `0.057665` is usefully comparable with the width of the stored
`D_K` corridor, but `C_U` is not `D_K`.  The lower finite cell, ordinary
Gamma-subtracted packet, grouped natural A lift, and all remaining equation-
(4) ownership terms must still be assembled with their signs intact.

Pi provenance: every `pi` in (QD1)--(QD8) comes from the inherited
Riemann-Siegel Gaussian, the fixed Mellin scaling `R=2*pi*Y`, the principal
quarter turn, or the odd Fourier roster.  No fitted geometric constant is
inserted.

## Proof boundary

Exact finite quarter-disk identity, exact vertical-current normalization,
cancellation-preserving unowned-cell reduction, and rigorous production
upper-arc/tail enclosure only.  No lower-cell enclosure, complete joined
packet, `J_Z`, `D_K`, non-A, all-height, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion.
""",
        encoding="utf-8",
    )
    print("certified quarter-disk vertical arc tail representation and production upper arc")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
