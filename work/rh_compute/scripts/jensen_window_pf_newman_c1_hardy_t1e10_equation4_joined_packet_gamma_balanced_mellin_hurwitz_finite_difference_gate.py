#!/usr/bin/env python3
"""Certify a gamma-balanced Mellin-Hurwitz coordinate for the joined packet."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any, Callable


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_gamma_balanced_mellin_hurwitz_finite_difference_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "entire_kernel": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_riemann_auxiliary_entire_kernel_finite_difference_gate.json",
    "geometric_half_line": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_QK_geometric_half_line_and_saddle_route_gate.json",
}

HEIGHT = 10_000_000_000
N_MINUS = 79_788
N_PLUS = 2_561_211
LABEL_COUNT = N_PLUS - N_MINUS
Q_MINUS = Fraction(2 * N_MINUS + 1, 2)
Q_PLUS = Fraction(2 * N_PLUS + 1, 2)
PI_DIGITS = "3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067"
PI_DENOMINATOR = 10 ** (len(PI_DIGITS) - 1)
PI_LOWER = Fraction(int(PI_DIGITS), PI_DENOMINATOR)
PI_UPPER = Fraction(int(PI_DIGITS) + 1, PI_DENOMINATOR)


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


def fmt(value: mp.mpf | mp.mpc, digits: int = 40) -> str:
    return mp.nstr(value, digits)


def geometric_difference(u: mp.mpf | mp.mpc, q_minus: mp.mpf, count: int) -> mp.mpc:
    if u == 0:
        return mp.mpc(count)
    return mp.exp(-q_minus * u) * mp.expm1(-count * u) / mp.expm1(-u)


def geometric_roster(x: mp.mpf, first_label: int, count: int) -> mp.mpc:
    omega = mp.exp(0.25j * mp.pi)
    return mp.fsum(
        mp.exp(-mp.pi * omega * (first_label + 2 * index) * x)
        for index in range(count)
    )


def integrate_zero_to_infinity(
    s: mp.mpc,
    remainder: Callable[[mp.mpf], mp.mpc],
    low_cuts: tuple[mp.mpf, ...],
    high_cuts: tuple[mp.mpf, ...],
) -> mp.mpc:
    def logarithmic(v: mp.mpf) -> mp.mpc:
        if not mp.isfinite(v):
            return mp.mpc(0)
        u = mp.exp(-v)
        return mp.exp(-(1 - s) * v) * remainder(u)

    require(low_cuts[-1] == mp.inf, "logarithmic cuts must end at infinity")
    require(high_cuts[-1] != mp.inf, "upper cutoff must be finite and explicit")
    low = mp.quad(logarithmic, low_cuts[:-1])
    low += mp.quadosc(logarithmic, [low_cuts[-2], mp.inf], omega=abs(mp.im(s)))
    high = mp.quad(lambda u: mp.power(u, -s) * remainder(u), high_cuts)
    return low + high


def scaled_prefactor(s: mp.mpc) -> mp.mpc:
    t = mp.im(s)
    return mp.exp(mp.pi * t / 2 + 0.25j * mp.pi) * mp.power(2 * mp.pi, s - 1)


def direct_source(
    s: mp.mpc,
    first_label: int,
    count: int,
    low_cuts: tuple[mp.mpf, ...],
    high_cuts: tuple[mp.mpf, ...],
) -> mp.mpc:
    t = mp.im(s)
    k0 = mp.exp(3 * mp.pi * t / 4 + 3j * mp.pi / 8)
    integral = integrate_zero_to_infinity(
        s,
        lambda x: mp.exp(-mp.pi * x * x) * geometric_roster(x, first_label, count),
        low_cuts,
        high_cuts,
    )
    return k0 * integral


def mellin_source(
    s: mp.mpc,
    first_label: int,
    count: int,
    low_cuts: tuple[mp.mpf, ...],
    high_cuts: tuple[mp.mpf, ...],
) -> tuple[mp.mpc, mp.mpc]:
    q_minus = mp.mpf(first_label) / 2
    integral = integrate_zero_to_infinity(
        s,
        lambda u: mp.exp(1j * u * u / (4 * mp.pi))
        * geometric_difference(u, q_minus, count),
        low_cuts,
        high_cuts,
    )
    return scaled_prefactor(s) * integral, integral


def hurwitz_difference(w: mp.mpc, q_minus: mp.mpf, q_plus: mp.mpf) -> mp.mpc:
    return mp.zeta(w, q_minus) - mp.zeta(w, q_plus)


def mellin_taylor_main(s: mp.mpc, q_minus: mp.mpf, q_plus: mp.mpf, order: int) -> mp.mpc:
    return mp.fsum(
        mp.power(1j / (4 * mp.pi), index)
        / mp.factorial(index)
        * mp.gamma(1 - s + 2 * index)
        * hurwitz_difference(1 - s + 2 * index, q_minus, q_plus)
        for index in range(order)
    )


def absolute_taylor_bound_log10(order: int) -> mp.mpf:
    q_minus = mp.mpf(Q_MINUS.numerator) / Q_MINUS.denominator
    q_plus = mp.mpf(Q_PLUS.numerator) / Q_PLUS.denominator
    w = mp.mpf(2 * order) + mp.mpf("0.5")
    delta = hurwitz_difference(w, q_minus, q_plus)
    require(mp.im(delta) == 0 and mp.re(delta) > 0, "real Hurwitz difference lost positivity")
    return (
        mp.pi * HEIGHT / (2 * mp.log(10))
        - mp.log(2 * mp.pi) / (2 * mp.log(10))
        + mp.loggamma(w) / mp.log(10)
        - order * mp.log(4 * mp.pi) / mp.log(10)
        - mp.loggamma(order + 1) / mp.log(10)
        + mp.log(mp.re(delta)) / mp.log(10)
    )


def upper_ratio_coordinate(order: int) -> Fraction:
    return Fraction(order, 1) + Fraction(3, 16 * (order + 1))


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    require(CHECKER.is_file(), "missing independent checker")
    dependency_rows: dict[str, dict[str, str]] = {}
    for name, path in DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {name}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        require(payload.get("passed") is True, f"dependency did not pass: {name}")
        dependency_rows[name] = {"path": relative(path), "sha256": file_hash(path)}

    mp.mp.dps = 90
    altered_s = mp.mpc("0.5", "3.25")
    altered_first = 5
    altered_count = 5
    low_cuts = (mp.mpf(0), mp.mpf(1), mp.mpf(3), mp.mpf(7), mp.mpf(15), mp.inf)
    high_cuts = (
        mp.mpf(1),
        mp.mpf(2),
        mp.mpf(4),
        mp.mpf(8),
        mp.mpf(16),
        mp.mpf(32),
        mp.mpf(64),
    )

    direct = direct_source(altered_s, altered_first, altered_count, low_cuts, high_cuts)
    transformed, transformed_integral = mellin_source(
        altered_s, altered_first, altered_count, low_cuts, high_cuts
    )
    transform_discrepancy = abs(direct - transformed)
    require(transform_discrepancy < mp.mpf("1e-55"), "scaled Mellin transform mismatch")

    q_minus = mp.mpf(altered_first) / 2
    q_plus = q_minus + altered_count
    anchor_integral = integrate_zero_to_infinity(
        altered_s,
        lambda u: geometric_difference(u, q_minus, altered_count),
        low_cuts,
        high_cuts,
    )
    anchor_hurwitz = mp.gamma(1 - altered_s) * hurwitz_difference(
        1 - altered_s, q_minus, q_plus
    )
    anchor_discrepancy = abs(anchor_integral - anchor_hurwitz)
    require(anchor_discrepancy < mp.mpf("1e-55"), "Hurwitz Mellin anchor mismatch")

    order = 3
    main_terms = mellin_taylor_main(altered_s, q_minus, q_plus, order)

    def taylor_remainder(u: mp.mpf) -> mp.mpc:
        x = u * u / (4 * mp.pi)
        polynomial = mp.fsum(mp.power(1j * x, index) / mp.factorial(index) for index in range(order))
        return geometric_difference(u, q_minus, altered_count) * (mp.exp(1j * x) - polynomial)

    remainder_integral = integrate_zero_to_infinity(
        altered_s, taylor_remainder, low_cuts, high_cuts
    )
    expansion_discrepancy = abs(transformed_integral - main_terms - remainder_integral)
    require(expansion_discrepancy < mp.mpf("1e-54"), "finite Mellin-Hurwitz expansion mismatch")

    explicit_bound = (
        abs(scaled_prefactor(altered_s))
        * mp.gamma(2 * order + mp.mpf("0.5"))
        * mp.re(hurwitz_difference(2 * order + mp.mpf("0.5"), q_minus, q_plus))
        / (mp.power(4 * mp.pi, order) * mp.factorial(order))
    )
    scaled_remainder = abs(scaled_prefactor(altered_s) * remainder_integral)
    require(scaled_remainder <= explicit_bound, "explicit Taylor remainder bound failed")

    balance_left = abs(scaled_prefactor(altered_s) * mp.gamma(1 - altered_s)) ** 2
    balance_right = 1 / (1 + mp.exp(-2 * mp.pi * mp.im(altered_s)))
    balance_discrepancy = abs(balance_left - balance_right)
    require(balance_discrepancy < mp.mpf("1e-75"), "gamma balance identity failed")

    threshold = 20_000_022_018
    q_squared = Q_MINUS * Q_MINUS
    coordinate_at_threshold = upper_ratio_coordinate(threshold)
    coordinate_after_threshold = upper_ratio_coordinate(threshold + 1)
    require(
        coordinate_at_threshold < PI_LOWER * q_squared,
        "pi lower bound does not certify the decreasing threshold",
    )
    require(
        coordinate_after_threshold > PI_UPPER * q_squared,
        "pi upper bound does not certify the first failed upper-ratio test",
    )
    require(
        upper_ratio_coordinate(threshold) - upper_ratio_coordinate(threshold - 1) > 0,
        "upper-ratio coordinate is not increasing",
    )

    scale_rows = [
        {"order": item, "absolute_bound_log10": fmt(absolute_taylor_bound_log10(item), 45)}
        for item in (1, 4, 8, 16, 32)
    ]

    artifact: dict[str, Any] = {
        "kind": STEM,
        "status": "gamma_balanced_mellin_hurwitz_finite_difference_and_absolute_taylor_barrier_certified",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "s": "1/2+i*t",
            "q_minus": f"{Q_MINUS.numerator}/{Q_MINUS.denominator}",
            "q_plus": f"{Q_PLUS.numerator}/{Q_PLUS.denominator}",
            "label_count": LABEL_COUNT,
        },
        "dependencies": dependency_rows,
        "scaled_mellin_coordinate": {
            "substitution": "u=2*pi*exp(i*pi/4)*x, followed within 0<=arg(u)<=pi/4 to the positive real axis",
            "finite_difference": "Delta_q(u)=(exp(-q_minus*u)-exp(-q_plus*u))/(1-exp(-u))",
            "removable_value": f"Delta_q(0)={LABEL_COUNT}",
            "integral": "M_s=integral_0^infinity u^(-s)*exp(i*u^2/(4*pi))*Delta_q(u)du",
            "source_identity": "S_W=exp(pi*t/2+i*pi/4)*(2*pi)^(s-1)*M_s",
            "contour_sector_has_gaussian_nongrowth": True,
            "altered_validation": {
                "t": "3.25",
                "first_label": altered_first,
                "label_count": altered_count,
                "direct_source": fmt(direct),
                "mellin_source": fmt(transformed),
                "discrepancy_absolute": fmt(transform_discrepancy),
            },
        },
        "gamma_balance": {
            "zero_gaussian_anchor": "integral u^(-s)*Delta_q(u)du=Gamma(1-s)*(zeta(1-s,q_minus)-zeta(1-s,q_plus))",
            "prefactor_identity": "|exp(pi*t/2+i*pi/4)*(2*pi)^(s-1)*Gamma(1-s)|^2=1/(1+exp(-2*pi*t))",
            "anchor_discrepancy_absolute": fmt(anchor_discrepancy),
            "balance_discrepancy_absolute": fmt(balance_discrepancy),
        },
        "finite_mellin_hurwitz_expansion": {
            "main": "sum_(k=0)^(K-1) (i/(4*pi))^k/k!*Gamma(1-s+2k)*Delta_zeta(1-s+2k)",
            "remainder": "integral u^(-s)*Delta_q(u)*(exp(i*x)-sum_(k<K)(i*x)^k/k!)du, x=u^2/(4*pi)",
            "absolute_bound": "|E_s*R_K|<=|E_s|*Gamma(2K+1/2)*Delta_zeta(2K+1/2)/((4*pi)^K*K!)",
            "altered_order": order,
            "altered_expansion_discrepancy_absolute": fmt(expansion_discrepancy),
            "altered_scaled_remainder_absolute": fmt(scaled_remainder),
            "altered_explicit_bound": fmt(explicit_bound),
        },
        "actual_height_scale_audit": {
            "bound_rows": scale_rows,
            "successive_bound_ratio": "B_(K+1)/B_K=[(2K+1/2)(2K+3/2)/(4*pi*(K+1))]*Delta_zeta(2K+5/2)/Delta_zeta(2K+1/2)",
            "ratio_upper_bound": "[K+3/(16(K+1))]/(pi*q_minus^2)",
            "last_order_with_certified_ratio_upper_bound_below_one": threshold,
            "first_order_where_that_upper_bound_exceeds_one": threshold + 1,
            "consequence": "this absolute Taylor majorant is strictly decreasing through more than twenty billion orders, so minimizing it is not a scalable evaluator",
            "pi_interval_digits": len(PI_DIGITS) - 1,
        },
        "decision": {
            "gamma_balanced_mellin_coordinate_selected": True,
            "absolute_taylor_remainder_route_rejected": True,
            "next_action": "retain the oscillatory Mellin remainder and derive a saddle-aligned contour or endpoint-complete Mordell representation; do not norm its Taylor remainder on the positive axis",
            "joined_packet_enclosed": False,
            "J_Z_enclosed": False,
            "D_K_enclosed": False,
            "rh_implication": False,
        },
        "runtime": {
            "seconds": time.perf_counter() - started,
            "active_compute_workers": 1,
            "priority": priority,
            "thread_caps": 1,
        },
        "source_hashes": {"builder": file_hash(BUILDER), "checker": file_hash(CHECKER)},
        "proof_boundary": (
            "Exact gamma-balanced Mellin-Hurwitz finite-difference coordinate, finite expansion with explicit "
            "remainder, and a rigorous non-scalability result for its positive-axis absolute Taylor majorant only. "
            "No saddle-aligned remainder bound, quantitative joined-packet, J_Z, or D_K enclosure, and no non-A, "
            "all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }

    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(
        f"""# Gamma-balanced Mellin-Hurwitz coordinate for the joined packet

Date: 2026-08-27

Status: exact scaled coordinate and explicit absolute-majorant barrier certified; saddle-aligned enclosure open.

Put

```text
omega=exp(i*pi/4),
q_-=79788+1/2,       q_+=2561211+1/2,
Delta_q(u)=[exp(-q_-u)-exp(-q_+u)]/[1-exp(-u)].
```

The denominator is removable at zero and `Delta_q(0)=2481423`.  In the
decaying half-line source use `u=2*pi*omega*x`.  The Gaussian becomes
`exp(i*u^2/(4*pi))`.  Rotating the resulting ray from `arg(u)=pi/4` to the
positive real axis is legal: the Gaussian does not grow in that sector, the
finite difference is regular at zero, and its real part decays at infinity.
Therefore

```text
M_s(q_-,q_+)
 =integral_0^infinity u^(-s)exp(i*u^2/(4*pi))Delta_q(u)du,

S_W=E_s M_s(q_-,q_+),
E_s=exp(pi*t/2+i*pi/4)(2*pi)^(s-1).                 (MH1)
```

Removing only the Gaussian chirp gives the exact finite Hurwitz anchor

```text
integral_0^infinity u^(-s)Delta_q(u)du
 =Gamma(1-s)[zeta(1-s,q_-)-zeta(1-s,q_+)].          (MH2)
```

The apparently large factor in (MH1) is exactly gamma-balanced:

```text
|E_s Gamma(1-s)|^2=1/[1+exp(-2*pi*t)].              (MH3)
```

Thus no exponential loss is intrinsic to this coordinate.

For every integer `K>=1`, Taylor's formula gives the exact finite expansion

```text
M_s=sum_(k=0)^(K-1) (i/(4*pi))^k/k!
       Gamma(1-s+2k)Delta_zeta(1-s+2k)+R_K,         (MH4)

Delta_zeta(w)=zeta(w,q_-)-zeta(w,q_+),

R_K=integral_0^infinity u^(-s)Delta_q(u)
    [exp(ix)-sum_(k<K)(ix)^k/k!]du,
x=u^2/(4*pi).                                       (MH5)
```

The elementary real-axis remainder estimate is explicit:

```text
|E_s R_K| <= |E_s| Gamma(2K+1/2)Delta_zeta(2K+1/2)
              /[(4*pi)^K K!].                       (MH6)
```

It is also unusable at the production height.  If `B_K` denotes the right
side of (MH6), positivity of the finite Hurwitz sum gives

```text
B_(K+1)/B_K
 <= [K+3/(16(K+1))]/(pi*q_-^2).                    (MH7)
```

Exact rational arithmetic with a 99-decimal enclosing interval for `pi`
proves that the right side of (MH7) is below one through

```text
K=20,000,022,018
```

and exceeds one at the next integer.  Consequently this particular absolute
Taylor majorant keeps decreasing for more than twenty billion orders before
it can even reach a turning point.  It is not a viable evaluator.  The useful
part is (MH1)--(MH4): the remainder must remain oscillatory and be moved onto
a saddle-aligned contour or an endpoint-complete Mordell representation.

An altered `t=3.25`, five-label row compares the original decaying half-line,
the scaled Mellin integral, the Hurwitz anchor, and the order-three finite
expansion.  The independent checker changes the height, roster, order,
quadrature partition, and precision.

Primary-source provenance: the entire auxiliary-kernel contour move follows
the no-residue mechanism in Arias de Reyna, https://arxiv.org/abs/2407.02016.
The unequal-truncation/Mordell route context is O'Sullivan,
https://arxiv.org/abs/1811.01130.  No implicit asymptotic constant is used.

Pi provenance: every `pi` comes from the inherited Riemann-Siegel Gaussian,
the fixed quarter-turn, the Mellin scaling, or the standard gamma identity.

## Proof boundary

{artifact['proof_boundary']}
""",
        encoding="utf-8",
    )
    print("certified gamma-balanced Mellin-Hurwitz finite difference and absolute Taylor barrier")
    print(f"result: {relative(RESULT)}")
    print(f"note: {relative(NOTE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
