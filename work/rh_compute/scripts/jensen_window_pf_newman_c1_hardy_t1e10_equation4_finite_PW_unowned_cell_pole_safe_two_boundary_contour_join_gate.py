#!/usr/bin/env python3
"""Join the finite RSI prefix to both unowned Fresnel tails pole-safely."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp
from sympy import Rational, euler


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_PW_unowned_cell_pole_safe_two_boundary_contour_join_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "branch_alignment": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_label_branch_correction_QK_alignment_gate.json",
    "fresnel_cells": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_common_carrier_subtraction_reduction_gate.json",
    "transition_packet": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_alternating_boundary_jet_packet_gate.json",
}

HEIGHT = 10_000_000_000
ALPHA_MIN = 159_577
ALPHA_MAX = 5_122_421
N_MINUS = 79_788
N_PLUS = 2_561_211
M = 2_481_423
LOWER_BOUNDARY = "621.5"
UPPER_BOUNDARY = "39936.5"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {relative(path)}")
    return json.loads(path.read_text(encoding="utf-8"))


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


def finite_source(z: mp.mpc, alpha_min: int, label_count: int) -> mp.mpc:
    return mp.fsum(mp.exp(1j * mp.pi * (alpha_min + 2 * j) * z) for j in range(label_count))


def paired_csc_source(z: mp.mpc, n_minus: int, n_plus: int) -> mp.mpc:
    return (
        0.5j
        * (mp.exp(2j * mp.pi * n_minus * z) - mp.exp(2j * mp.pi * n_plus * z))
        / mp.sin(mp.pi * z)
    )


def unregularized_csc_packet(z: mp.mpc, s: mp.mpc, terms: int) -> mp.mpc:
    return (
        mp.power(z, -s)
        * mp.exp(-1j * mp.pi * z * z + 2j * mp.pi * terms * z)
        / (2j * mp.sin(mp.pi * z))
    )


def regularized_prefix_packet(z: mp.mpc, s: mp.mpc, terms: int) -> mp.mpc:
    return unregularized_csc_packet(z, s, terms) - unregularized_csc_packet(z, s, 0)


def direct_prefix_packet(z: mp.mpc, s: mp.mpc, terms: int) -> mp.mpc:
    source = mp.fsum(mp.exp(1j * mp.pi * (2 * j + 1) * z) for j in range(terms))
    return mp.power(z, -s) * mp.exp(-1j * mp.pi * z * z) * source


def half_integer_source(alpha_min: int, label_count: int, k: int) -> mp.mpc:
    alternating_sum = 1 if label_count % 2 else 0
    phase = mp.exp(1j * mp.pi * alpha_min * (mp.mpf(k) + mp.mpf("0.5")))
    return phase * alternating_sum


def alternating_source_moment(alpha_min: int, label_count: int, order: int) -> int:
    x = Rational(alpha_min, 2)
    value = Rational(2) ** (order - 1) * (
        euler(order, x) - ((-1) ** label_count) * euler(order, x + label_count)
    )
    require(value.q == 1, "alternating source moment is not integral")
    return int(value)


def source_and_correction(t: mp.mpf, alpha: int) -> tuple[mp.mpc, mp.mpc]:
    pi = mp.pi
    a = mp.mpf("0.25") - 0.5j * t
    b = mp.mpf("0.75") - 0.5j * t
    xi = 1j * pi * alpha * alpha / 4
    c = pi * alpha * (1 + 1j) / mp.sqrt(2)
    d = mp.exp(-pi * t / 4 - 1j * pi / 8)
    even = mp.gamma(a) * pi ** (-a) * mp.hyp1f1(a, mp.mpf("0.5"), xi) / 2
    odd = c * mp.gamma(b) * pi ** (-b) * mp.hyp1f1(b, mp.mpf("1.5"), xi) / 2
    source = 1j * mp.exp(pi * t) * d * (even - odd)
    correction = d * (even + odd)
    return source, correction


def direct_contour_prefix(t: mp.mpf, alpha_min: int, label_count: int) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * t
    rotation = mp.exp(-1j * mp.pi / 4)

    def integrand(q: mp.mpf) -> mp.mpc:
        z = mp.mpf("0.5") + q * rotation
        return (
            -mp.power(z, -s)
            * mp.exp(-1j * mp.pi * z * z)
            * finite_source(z, alpha_min, label_count)
            * rotation
        )

    return mp.quad(integrand, [-mp.inf, -6, -1, 0, 1, 3, 6, 9, 13, 20, mp.inf])


def integrated_source_lower_series(
    t: mp.mpf,
    alpha_min: int,
    label_count: int,
    u_cut: mp.mpf,
    term_count: int,
) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * t
    x_cut = mp.exp(u_cut)
    lam = mp.pi * (1 + 1j) / mp.sqrt(2)
    total = mp.mpc(0)
    for j in range(label_count):
        linear = -lam * (alpha_min + 2 * j)
        quadratic = -mp.pi
        previous_two = mp.mpc(0)
        previous_one = mp.mpc(1)
        for order in range(term_count):
            if order == 0:
                coefficient = mp.mpc(1)
            elif order == 1:
                coefficient = linear
            else:
                coefficient = (
                    linear * previous_one + 2 * quadratic * previous_two
                ) / order
            total += coefficient * mp.power(x_cut, order + 1 - s) / (order + 1 - s)
            previous_two, previous_one = previous_one, coefficient
    return total


def direct_source_half_line(
    t: mp.mpf,
    alpha_min: int,
    label_count: int,
) -> tuple[mp.mpc, mp.mpf, mp.mpf]:
    s = mp.mpf("0.5") + 1j * t
    lam = mp.pi * (1 + 1j) / mp.sqrt(2)
    k_zero = mp.exp(3 * mp.pi * t / 4 + 3j * mp.pi / 8)
    u_lower = mp.mpf("-10")
    u_upper = mp.mpf("2")

    def integrand(u: mp.mpf) -> mp.mpc:
        x = mp.exp(u)
        source = mp.fsum(mp.exp(-lam * (alpha_min + 2 * j) * x) for j in range(label_count))
        return mp.exp((1 - s) * u) * mp.exp(-mp.pi * x * x) * source

    lower_36 = integrated_source_lower_series(t, alpha_min, label_count, u_lower, 36)
    lower_44 = integrated_source_lower_series(t, alpha_min, label_count, u_lower, 44)
    lower_series_delta = abs(k_zero * (lower_44 - lower_36))
    middle = mp.quad(integrand, [u_lower, -7, -4, -2, 0, u_upper])

    x_upper = mp.exp(u_upper)
    decay_rate = mp.pi / mp.sqrt(2)
    exponent = mp.pi * x_upper * x_upper + decay_rate * alpha_min * x_upper
    derivative_floor = 2 * mp.pi * x_upper + decay_rate * alpha_min
    upper_tail_bound = (
        abs(k_zero)
        * label_count
        * mp.power(x_upper, mp.mpf("-0.5"))
        * mp.exp(-exponent)
        / derivative_floor
    )
    return k_zero * (lower_44 + middle), lower_series_delta, upper_tail_bound


def direct_finite_fresnel(
    t: mp.mpf,
    alpha_min: int,
    label_count: int,
    lower: mp.mpf,
    upper: mp.mpf,
) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * t

    def integrand(y: mp.mpf) -> mp.mpc:
        return (
            mp.power(y, -s)
            * mp.exp(-1j * mp.pi * y * y)
            * finite_source(y, alpha_min, label_count)
        )

    points = [lower]
    cursor = lower + mp.mpf("0.25")
    while cursor < upper:
        points.append(cursor)
        cursor += mp.mpf("0.25")
    points.append(upper)
    return mp.quad(integrand, points)


def altered_roster_row(
    t_text: str,
    alpha_min: int,
    alpha_max: int,
    lower_text: str,
    upper_text: str,
) -> dict[str, Any]:
    mp.mp.dps = 75
    t = mp.mpf(t_text)
    lower = mp.mpf(lower_text)
    upper = mp.mpf(upper_text)
    label_count = (alpha_max - alpha_min) // 2 + 1

    closed_source = mp.mpc(0)
    closed_correction = mp.mpc(0)
    for alpha in range(alpha_min, alpha_max + 1, 2):
        source, correction = source_and_correction(t, alpha)
        closed_source += source
        closed_correction += correction

    direct_prefix = direct_contour_prefix(t, alpha_min, label_count)
    direct_source, lower_series_delta, upper_tail_bound = direct_source_half_line(
        t, alpha_min, label_count
    )
    finite_interior = direct_finite_fresnel(t, alpha_min, label_count, lower, upper)
    joined_closed = closed_source - finite_interior
    joined_contour = direct_prefix + closed_correction - finite_interior

    n_minus = (alpha_min - 1) // 2
    n_plus = (alpha_max + 1) // 2
    regularized_j = {}
    for terms in (n_minus, n_plus):
        prefix_source = mp.mpc(0)
        for alpha in range(1, 2 * terms, 2):
            source, _ = source_and_correction(t, alpha)
            prefix_source += source
        prefix_interior = (
            direct_finite_fresnel(t, 1, terms, lower, upper) if terms else mp.mpc(0)
        )
        regularized_j[terms] = prefix_interior - prefix_source
    regularized_join = regularized_j[n_minus] - regularized_j[n_plus]

    source_discrepancy = abs(closed_source - direct_source)
    branch_discrepancy = abs(closed_source - direct_prefix - closed_correction)
    join_discrepancy = abs(joined_closed - joined_contour)
    regularized_join_discrepancy = abs(joined_closed - regularized_join)
    require(
        source_discrepancy < mp.mpf("1e-52"),
        f"altered source half-line replay failed: {mp.nstr(source_discrepancy, 18)}",
    )
    require(
        lower_series_delta < mp.mpf("1e-64"),
        f"altered lower-series stabilization failed: {mp.nstr(lower_series_delta, 18)}",
    )
    require(
        upper_tail_bound < mp.mpf("1e-80"),
        f"altered upper-tail bound is too large: {mp.nstr(upper_tail_bound, 18)}",
    )
    require(
        branch_discrepancy < mp.mpf("1e-52"),
        f"altered branch/contour replay failed: {mp.nstr(branch_discrepancy, 18)}",
    )
    require(
        join_discrepancy < mp.mpf("1e-52"),
        f"altered two-boundary join replay failed: {mp.nstr(join_discrepancy, 18)}",
    )
    require(
        regularized_join_discrepancy < mp.mpf("1e-52"),
        f"altered regularized J_N replay failed: {mp.nstr(regularized_join_discrepancy, 18)}",
    )

    return {
        "t": t_text,
        "alpha_min": alpha_min,
        "alpha_max": alpha_max,
        "label_count": label_count,
        "n_minus": n_minus,
        "n_plus": n_plus,
        "lower_boundary": lower_text,
        "upper_boundary": upper_text,
        "closed_source": mp.nstr(closed_source, 58),
        "direct_rotated_negative_boundary": mp.nstr(direct_source, 58),
        "direct_prefix_contour": mp.nstr(direct_prefix, 58),
        "closed_full_branch_correction": mp.nstr(closed_correction, 58),
        "finite_interior_fresnel": mp.nstr(finite_interior, 58),
        "joined_unowned_packet": mp.nstr(joined_closed, 58),
        "regularized_J_n_minus": mp.nstr(regularized_j[n_minus], 58),
        "regularized_J_n_plus": mp.nstr(regularized_j[n_plus], 58),
        "source_half_line_discrepancy_absolute": mp.nstr(source_discrepancy, 18),
        "source_lower_series_36_to_44_delta_absolute": mp.nstr(lower_series_delta, 18),
        "source_upper_tail_omission_bound_absolute": mp.nstr(upper_tail_bound, 18),
        "branch_identity_discrepancy_absolute": mp.nstr(branch_discrepancy, 18),
        "two_boundary_join_discrepancy_absolute": mp.nstr(join_discrepancy, 18),
        "regularized_J_difference_discrepancy_absolute": mp.nstr(regularized_join_discrepancy, 18),
    }


def render_note(artifact: dict[str, Any]) -> str:
    replay = artifact["altered_roster_replay"]
    return f"""# Pole-safe two-boundary contour join for the unowned cells

Date: 2026-08-27

Status: exact paired-contour join certified; quantitative Mordell compression open

Put `s=1/2+it`, use the principal logarithm, and define

```text
D_W(z)=sum_(j=0)^(M-1) exp[i*pi*(A+2j)z],
g_W(z)=z^(-s)exp(-i*pi*z^2)D_W(z),
A=159577, M=2481423.                                (UJ1)
```

The finite RSI prefix and the complete branch correction are

```text
P_W=-integral_C g_W(z)dz,
B_W=integral_0^infinity[Fresnel] g_W(y)dy,
S_W=P_W+B_W.                                        (UJ2)
```

For `L=621.5` and `U=39936.5`, the lower and upper unowned cells therefore
join before any norm as

```text
U_unowned
 =P_W+integral_0^L g_W+integral_U^infinity[Fresnel]g_W
 =S_W-integral_L^U g_W.                             (UJ3)
```

Deforming `C` to the upper boundary of the real axis is legal for the finite
source.  The indentation at zero vanishes because `Re(s)=1/2`.  Hence

```text
U_unowned
 =-integral_(-infinity)^0[upper] g_W(z)dz
  -integral_L^U g_W(y)dy.                           (UJ4)
```

Writing `z=-y+i0` in the first term and rotating
`y=exp(-i*pi/4)x` through the decaying quadrant gives exactly

```text
-integral_(-infinity)^0[upper]g_W(z)dz
 =K_0 integral_0^infinity x^(-s)exp(-pi*x^2)G_W(x)dx
 =S_W,

K_0=exp(3*pi*t/4+3*pi*i/8),
G_W(x)=sum_(j=0)^(M-1)exp[-pi(1+i)(A+2j)x/sqrt(2)]. (UJ5)
```

Thus (UJ3)--(UJ5) are the same identity in cell, real-boundary, and decaying
half-line normalizations.

The finite source also has the exact csc form

```text
D_W(z)=(i/2)[exp(2*pi*i*n_-*z)-exp(2*pi*i*n_+*z)]csc(pi*z),
n_-=79788, n_+=2561211.                             (UJ6)
```

If

```text
G_N(z)=(1/(2i))z^(-s)exp(-i*pi*z^2+2*pi*i*N*z)csc(pi*z),
```

then `g_W=G_(n_+)-G_(n_-)`.  This difference must remain inside every real-
boundary integral.  At each nonzero integer `k`, both `G_N` terms have the
same residue `k^(-s)/(2*pi*i)`, so the poles cancel.  At zero each separate
term behaves as `z^(-s-1)/(2*pi*i)`, whereas the paired difference behaves
as `M z^(-s)` and is locally integrable.  Consequently separate objects
formed directly from the unregularized `G_N` are not legitimate ordinary
real-tail integrals.

There is, however, a canonical common subtraction.  Put

```text
R_N(z)=G_N(z)-G_0(z)
      =z^(-s)exp(-i*pi*z^2)
       sum_(r=0)^(N-1)exp[i*pi*(2r+1)z].             (UJ7)
```

Every `R_N` is a finite-source packet: its integer poles are removable and it
behaves as `N z^(-s)` at zero.  Therefore the individual Fresnel objects

```text
mathcal J_N(L,U)
 =integral_C R_N(z)dz
  -integral_0^L R_N(y)dy
  -integral_U^infinity[Fresnel]R_N(y)dy
```

are legitimate, and the common subtraction cancels in their difference:

```text
U_unowned=mathcal J_(n_-)(L,U)-mathcal J_(n_+)(L,U). (UJ8)
```

If `S_N` is the finite physical Kummer prefix and `B_N[L,U]` is the
corresponding finite-source interior Fresnel integral, then the same
definitions give the endpoint-complete identity

```text
mathcal J_N(L,U)=B_N[L,U]-S_N.                       (UJ8a)
```

This is the pole-safe finite/infinite contour normalization required for a
Mordell transform.  It is exact, but it is not yet a compressed or bounded
evaluation of either `mathcal J_N`.

Because `A=1 mod 4`, `M` is odd, and both production boundaries are half
integers,

```text
D_W(k+1/2)=i(-1)^k,

D_W^(r)(k+1/2)
 =i(-1)^k(i*pi)^r 2^(r-1)
   [E_r(A/2)+E_r(A/2+M)].                           (UJ9)
```

Here `E_r` is the Euler polynomial; the plus sign uses odd `M`.  This gives
all source endpoint jets in constant time for the next endpoint-complete
transformation.  If `f(z)=z^(-s)exp(-i*pi*z^2)`, then

```text
g_W^(r)(h)=sum_(q=0)^r binom(r,q)f^(r-q)(h)D_W^(q)(h),
                                      h in {{L,U}},  (UJ10)

(log f)' =-s/z-2i*pi*z,
(log f)''= s/z^2-2i*pi,
(log f)^(q)=(-1)^q s(q-1)!/z^q,       q>=3.
```

Together with (UJ9), this is a constant-work recurrence for every full
integrand endpoint jet.  The altered `{replay['label_count']}`-label replay at
`t={replay['t']}` compares the hypergeometric source, direct RSI contour,
rotated negative boundary, branch correction, and finite interior integral.
Its largest join discrepancy is
`{replay['two_boundary_join_discrepancy_absolute']}`.

There is also an exact cancellation guard for the final assembly.  With
`a=39852.5`, write

```text
O_join= integral_L^a g_W
        -C_G sum_(m=622)^39852 m^(-s),

T_join= integral_a^U g_W
        -C_G sum_(m=39853)^39936 m^(-s)-mathcal A_A^nat.
```

Then (UJ3) gives, before projection,

```text
U_unowned+O_join+T_join
 =S_W-C_G sum_(m=622)^39936 m^(-s)-mathcal A_A^nat. (UJ11)
```

The two positive Fresnel subintervals cancel the `-integral_L^U g_W` in
`U_unowned` exactly.  Therefore the three packets must be assembled before a
final norm; (UJ11) is the cancellation-preserving complete target, not three
unrelated absolute-value estimates.

Pi provenance: every `pi` in (UJ1)--(UJ11) descends from the original RSI
Gaussian/sine kernel, the odd Fourier roster, the fixed quarter-turn, or the
Riemann-Siegel branch phase.  No fitted circle constant is inserted.

Proof boundary: (UJ3)--(UJ11) certify the exact pole-safe contour join and
remove the separate infinite upper-tail obligation.  They do not compress or
bound the two-boundary packet, bound the ordinary Gamma-subtracted packet,
enclose `J_Z` or `D_K`, prove a non-A bound or all-height theorem, or prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")

    require(ALPHA_MIN == 2 * N_MINUS + 1, "lower source/prefix map drift")
    require(ALPHA_MAX == 2 * N_PLUS - 1, "upper source/prefix map drift")
    require(N_PLUS - N_MINUS == M, "source count drift")
    require(ALPHA_MIN % 4 == 1 and M % 2 == 1, "half-integer phase hypotheses drift")

    mp.mp.dps = 80
    source_rows = []
    regularized_rows = []
    surrogate_s = mp.mpf("0.5") + mp.mpf("3.125") * 1j
    for z in (mp.mpc("0.37", "0.21"), mp.mpc("1.25", "-0.17"), mp.mpc("2.5", "0.13")):
        direct = finite_source(z, 5, 9)
        paired = paired_csc_source(z, 2, 11)
        discrepancy = abs(direct - paired)
        require(discrepancy < mp.mpf("1e-70"), "paired csc source identity failed")
        source_rows.append({"z": mp.nstr(z, 20), "discrepancy_absolute": mp.nstr(discrepancy, 18)})
        for terms in (1, 2, 5, 11):
            regularized = regularized_prefix_packet(z, surrogate_s, terms)
            finite_prefix = direct_prefix_packet(z, surrogate_s, terms)
            regularized_discrepancy = abs(regularized - finite_prefix)
            require(regularized_discrepancy < mp.mpf("1e-69"), "regularized finite-prefix identity failed")
            regularized_rows.append(
                {
                    "z": mp.nstr(z, 20),
                    "N": terms,
                    "discrepancy_absolute": mp.nstr(regularized_discrepancy, 18),
                }
            )

    endpoint_rows = []
    for k in (-2, -1, 0, 1, 4, 7):
        value = half_integer_source(ALPHA_MIN, M, k)
        target = 1j * (1 if k % 2 == 0 else -1)
        discrepancy = abs(value - target)
        require(discrepancy < mp.mpf("1e-60"), "production half-integer endpoint phase failed")
        endpoint_rows.append({"k": k, "value": mp.nstr(value, 22), "discrepancy_absolute": mp.nstr(discrepancy, 18)})

    source_moments = []
    for order in range(9):
        source_moments.append(
            {
                "order": order,
                "alternating_power_sum": str(alternating_source_moment(ALPHA_MIN, M, order)),
                "jet_formula": "i*(-1)^k*(i*pi)^r*alternating_power_sum",
            }
        )

    replay = altered_roster_row("3.25", 5, 13, "1.5", "3.5")

    artifact: dict[str, Any] = {
        "kind": STEM,
        "status": "exact_PW_unowned_cells_pole_safe_two_boundary_contour_join_certified",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "s": "1/2+i*t",
            "contour": "C: z=1/2+q*exp(-i*pi/4), q from -infinity to infinity",
            "log_branch": "principal; negative real boundary approached from above",
        },
        "actual_roster": {
            "alpha_min": ALPHA_MIN,
            "alpha_max": ALPHA_MAX,
            "label_count": M,
            "n_minus": N_MINUS,
            "n_plus": N_PLUS,
            "lower_boundary": LOWER_BOUNDARY,
            "upper_boundary": UPPER_BOUNDARY,
        },
        "exact_join": {
            "finite_source": "D_W(z)=sum_(j=0)^(M-1)exp(i*pi*(A+2j)*z)",
            "integrand": "g_W(z)=z^(-s)exp(-i*pi*z^2)D_W(z)",
            "unowned_cell_definition": "U_unowned=P_W+integral_0^L g_W+integral_U^infinity[Fresnel]g_W",
            "cell_join": "U_unowned=S_W-integral_L^U g_W",
            "two_boundary_join": "U_unowned=-integral_(-infinity)^0[upper]g_W-integral_L^U g_W",
            "rotated_negative_boundary": "-integral_(-infinity)^0[upper]g_W=K_0*integral_0^infinity x^(-s)exp(-pi*x^2)G_W(x)dx=S_W",
            "separate_upper_infinite_tail_remaining": False,
        },
        "complete_packet_reassembly": {
            "ordinary_packet": "O_join=integral_L^a g_W-C_G*sum_(m=622)^39852m^(-s)",
            "transition_packet": "T_join=integral_a^U g_W-C_G*sum_(m=39853)^39936m^(-s)-mathcal_A_A^nat",
            "identity": "U_unowned+O_join+T_join=S_W-C_G*sum_(m=622)^39936m^(-s)-mathcal_A_A^nat",
            "positive_fresnel_strip_cancels_exactly": True,
            "final_norm_requires_complete_reassembly": True,
        },
        "paired_csc_packet": {
            "source_identity": "D_W=(i/2)(exp(2*pi*i*n_minus*z)-exp(2*pi*i*n_plus*z))*csc(pi*z)",
            "G_N": "(1/(2*i))*z^(-s)*exp(-i*pi*z^2+2*pi*i*N*z)*csc(pi*z)",
            "paired_identity": "g_W=G_(n_plus)-G_(n_minus)",
            "common_nonzero_integer_residue": "k^(-s)/(2*pi*i)",
            "zero_behavior_separate": "G_N(z)~z^(-s-1)/(2*pi*i)",
            "zero_behavior_paired": "g_W(z)~M*z^(-s)",
            "separate_real_tail_J_N_ordinary_integrals_permitted": False,
            "pair_before_real_boundary_integration_required": True,
            "numerical_source_rows": source_rows,
        },
        "regularized_prefix_packets": {
            "definition": "R_N=G_N-G_0",
            "finite_prefix_identity": "R_N=z^(-s)exp(-i*pi*z^2)*sum_(r=0)^(N-1)exp(i*pi*(2r+1)z)",
            "integer_poles_removable": True,
            "zero_behavior": "R_N(z)~N*z^(-s)",
            "joined_packet": "mathcal_J_N(L,U)=integral_C R_N-integral_0^L R_N-integral_U^infinity[Fresnel]R_N",
            "unowned_identity": "U_unowned=mathcal_J_(n_minus)(L,U)-mathcal_J_(n_plus)(L,U)",
            "endpoint_complete_identity": "mathcal_J_N(L,U)=B_N[L,U]-S_N",
            "individual_regularized_J_N_integrals_legitimate": True,
            "numerical_rows": regularized_rows,
        },
        "half_integer_endpoints": {
            "hypotheses": "A=1 mod 4 and M odd",
            "formula": "D_W(k+1/2)=i*(-1)^k",
            "all_order_formula": "D_W^(r)(k+1/2)=i*(-1)^k*(i*pi)^r*2^(r-1)*(E_r(A/2)+E_r(A/2+M))",
            "full_integrand_jet": "g_W^(r)(h)=sum_(q=0)^r binom(r,q)f^(r-q)(h)D_W^(q)(h), f=z^(-s)exp(-i*pi*z^2)",
            "log_f_derivatives": [
                "(log f)'=-s/z-2*i*pi*z",
                "(log f)''=s/z^2-2*i*pi",
                "(log f)^(q)=(-1)^q*s*(q-1)!/z^q for q>=3",
            ],
            "constant_work_full_integrand_endpoint_jets": True,
            "production_rows": endpoint_rows,
            "production_alternating_moments_orders_0_through_8": source_moments,
        },
        "altered_roster_replay": replay,
        "decision": {
            "P_W_joined_exactly_to_lower_and_upper_unowned_cells": True,
            "paired_integer_pole_cancellation_proved": True,
            "paired_zero_integrability_proved": True,
            "separate_csc_tail_packets_rejected": True,
            "common_G0_regularization_proved": True,
            "regularized_finite_prefix_J_N_join_proved": True,
            "complete_three_packet_fresnel_cancellation_proved": True,
            "quantitative_endpoint_complete_Mordell_compression_proved": False,
            "joined_unowned_packet_enclosed_at_actual_height": False,
            "ordinary_owned_packet_enclosed": False,
            "D_K_enclosed": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "references": {
            "gabcke": "https://arxiv.org/abs/1512.01186",
            "kuznetsov": "https://arxiv.org/abs/1306.4081",
        },
        "runtime": {
            "seconds": time.time() - started,
            "priority": priority,
            "active_compute_workers": 1,
            "thread_caps": 1,
        },
        "next_obligation": "Derive an endpoint-complete Mordell transformation or rigorous grouped quadrature for the regularized difference mathcal_J_(n_minus)-mathcal_J_(n_plus), or work in the completely reassembled form (UJ11). Use the constant-work endpoint jets (UJ9)--(UJ10), retain the common G_0 subtraction through every contour move, and do not norm the three physical packets separately.",
        "proof_boundary": "Exact pole-safe P_W/unowned-cell two-boundary contour join, common-G_0 regularized finite-prefix packets, all-order half-integer endpoint jets, and complete cancellation-preserving packet reassembly only. No quantitative Mordell compression, actual-height joined-packet bound, ordinary packet bound, J_Z or D_K enclosure, non-A or all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }

    note = render_note(artifact)
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(note, encoding="utf-8")
    print("certified pole-safe P_W/unowned-cell two-boundary contour join", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
