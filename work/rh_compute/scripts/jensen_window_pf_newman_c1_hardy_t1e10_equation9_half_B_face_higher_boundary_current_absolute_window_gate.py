#!/usr/bin/env python3
"""Bound the next two pole-subtracted B face currents on the xi window."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_higher_boundary_current_absolute_window_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "signed_first_current": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_signed_window_gate.json",
    "triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
}

PRECISION = 90
T = 10_000_000_000
B = 5_122_421
WINDOW_XI = 70
SUM_CUTOFF = 100_000
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
    z, x, endpoint, mode, c = sp.symbols("z x B m c", positive=True, real=True)
    phase = sp.pi * mode**2 / z + sp.pi * endpoint**2 * z / 4
    driver = z ** (-sp.Rational(1, 2)) * (2 + sp.I * sp.pi * endpoint**2 * z)
    i_phase = sp.I * sp.diff(phase, z)
    g0 = driver / i_phase
    g1 = sp.diff(g0, z) / i_phase
    g2 = sp.diff(g1, z) / i_phase
    prefactor = x ** (-sp.Rational(3, 2)) / (2 * sp.I * sp.pi)
    terms = [
        prefactor * g0.subs(z, x),
        -prefactor * g1.subs(z, x),
        prefactor * g2.subs(z, x),
    ]
    reduced = [sp.factor(term.subs(x, 2 * c / endpoint)) for term in terms]
    expected0 = -(1 + sp.I * sp.pi * endpoint * c) / (sp.pi**2 * (c**2 - mode**2))
    require(sp.simplify(reduced[0] - expected0) == 0, "first boundary term drift")

    r2, r3, r4, r5 = sp.symbols("R_2 R_3 R_4 R_5")
    current1 = c / (sp.pi**3 * endpoint) * (
        (-4 * sp.pi * endpoint * c**3 + 4 * sp.I * c**2) * r3
        + (5 * sp.pi * endpoint * c - 3 * sp.I) * r2
    )
    current2 = -sp.I * c**2 / (sp.pi**4 * endpoint**2) * (
        (-48 * sp.pi * endpoint * c**5 + 48 * sp.I * c**4) * r5
        + (84 * sp.pi * endpoint * c**3 - 60 * sp.I * c**2) * r4
        + (-35 * sp.pi * endpoint * c + 15 * sp.I) * r3
    )
    return {
        "IBP_operator": "g_0=h/(i*phi'), g_(j+1)=g_j'/(i*phi'); boundary terms alternate +g_0,-g_1,+g_2",
        "first_current": str(reduced[0]),
        "second_current_sum": str(current1),
        "third_current_sum": str(current2),
        "lattice_sums": "R_k(c)=sum_(m>=1,m notin {621,622})(c^2-m^2)^(-k), k=2,...,5",
        "derivative_recurrence": "R_(k+1)=-R_k'/(2*k*c)",
        "tail_rule": "For m>M, |m^2-c^2|^(-k)<=rho^(-k)m^(-2k), rho=1-(c_high/(M+1))^2",
    }


def absolute_lattice_sums(c_low: arb, c_high: arb) -> dict[int, arb]:
    sums = {order: arb(0) for order in range(2, 6)}
    for mode in range(1, 621):
        denominator = c_low**2 - arb(mode) ** 2
        require(denominator > 0, "left denominator lost positivity")
        for order in sums:
            sums[order] += denominator ** (-order)
    for mode in range(623, SUM_CUTOFF + 1):
        denominator = arb(mode) ** 2 - c_high**2
        require(denominator > 0, "right denominator lost positivity")
        for order in sums:
            sums[order] += denominator ** (-order)

    next_mode = arb(SUM_CUTOFF + 1)
    rho = 1 - (c_high / next_mode) ** 2
    require(rho > arb("0.9999"), "tail rho unexpectedly small")
    cutoff = arb(SUM_CUTOFF)
    for order in sums:
        sums[order] += rho ** (-order) * cutoff ** (1 - 2 * order) / (2 * order - 1)
    return sums


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t, endpoint, pi = arb(T), arb(B), arb.pi()
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = hessian.sqrt()
    x_low = x0 - arb(WINDOW_XI) / root_hessian
    x_high = x0 + arb(WINDOW_XI) / root_hessian
    c_low, c_high = endpoint * x_low / 2, endpoint * x_high / 2
    sums = absolute_lattice_sums(c_low, c_high)

    c = c_high
    current1_sup = c / (pi**3 * endpoint) * (
        (4 * pi * endpoint * c**3 + 4 * c**2) * sums[3]
        + (5 * pi * endpoint * c + 3) * sums[2]
    )
    current2_sup = c**2 / (pi**4 * endpoint**2) * (
        (48 * pi * endpoint * c**5 + 48 * c**4) * sums[5]
        + (84 * pi * endpoint * c**3 + 60 * c**2) * sums[4]
        + (35 * pi * endpoint * c + 15) * sums[3]
    )

    # x(1-x) is increasing on this window, which lies below 1/2.
    weight_max = (x_low * (1 - x_low)) ** (-arb(1) / 4)
    width = 2 * arb(WINDOW_XI) / root_hessian
    physical_factor = 2 * (pi / (32 * t)) ** (arb(1) / 4)
    current1_physical = physical_factor * width * weight_max * current1_sup
    current2_physical = physical_factor * width * weight_max * current2_sup
    combined = current1_physical + current2_physical
    require(current1_physical < arb("1.49e-6"), "second face current exceeds 1.49e-6")
    require(current2_physical < arb("2.1e-10"), "third face current exceeds 2.1e-10")
    require(combined < arb("1.491e-6"), "combined higher face currents exceed 1.491e-6")

    return {
        "height": T,
        "endpoint": B,
        "window_xi": [-WINDOW_XI, WINDOW_XI],
        "sum_cutoff": SUM_CUTOFF,
        "x_window": [x_low.str(PRECISION, more=True), x_high.str(PRECISION, more=True)],
        "c_window": [c_low.str(PRECISION, more=True), c_high.str(PRECISION, more=True)],
        "kummer_weight_max_ball": weight_max.str(PRECISION, more=True),
        "kummer_weight_max_location": "x_window_low",
        "absolute_lattice_sum_balls": {str(order): value.str(PRECISION, more=True) for order, value in sums.items()},
        "second_current_uniform_absolute_ball": current1_sup.str(PRECISION, more=True),
        "third_current_uniform_absolute_ball": current2_sup.str(PRECISION, more=True),
        "second_current_physical_absolute_integral_ball": current1_physical.str(PRECISION, more=True),
        "third_current_physical_absolute_integral_ball": current2_physical.str(PRECISION, more=True),
        "combined_higher_current_physical_absolute_integral_ball": combined.str(PRECISION, more=True),
        "combined_higher_current_physical_absolute_upper_bound": "1.491e-6",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Higher pole-subtracted B face currents on the Gaussian window

Date: 2026-08-13

Status: exact boundary-current algebra plus saved-height absolute enclosure;
not a proof of the complete B-face expansion

For the exact B endpoint triangle put

```text
phi(z)=pi*m^2/z+pi*B^2*z/4,
h(z)=z^(-1/2)(2+i*pi*B^2*z).
```

Repeated integration by parts with
`g_0=h/(i*phi')`, `g_(j+1)=g_j'/(i*phi')` produces alternating face
currents `+g_0,-g_1,+g_2`.  After `c=Bx/2` and summation over all positive
modes except 621 and 622, define

```text
R_k(c)=sum_(m>=1,m notin {{621,622}})(c^2-m^2)^(-k).  (HC1)
```

The next two exact currents after Section 11.373 are

```text
C_1=c/(pi^3 B){{(-4pi Bc^3+4i c^2)R_3
                 +(5pi Bc-3i)R_2}},                  (HC2)

C_2=-i c^2/(pi^4 B^2){{(-48pi Bc^5+48i c^4)R_5
                 +(84pi Bc^3-60i c^2)R_4
                 +(-35pi Bc+15i)R_3}}.               (HC3)
```

These formulas retain the signed pair before lattice summation.  Since the
dangerous first current has already been integrated with its trace phase,
the two higher currents may now be bounded absolutely.  On `|xi|<=70`, a
finite Arb sum through `{c['sum_cutoff']}` plus the analytic tail

```text
sum_(m>M)|m^2-c^2|^(-k)
 <=rho^(-k) M^(1-2k)/(2k-1),
rho=1-[c_high/(M+1)]^2,                               (HC4)
```

gives

```text
sup |C_1| <= {c['second_current_uniform_absolute_ball']},
sup |C_2| <= {c['third_current_uniform_absolute_ball']}. (HC5)
```

Restoring the complete Kummer weight, window width, real projection factor,
and `(pi/(32t))^(1/4)` normalization yields

```text
E_1 <= {c['second_current_physical_absolute_integral_ball']},
E_2 <= {c['third_current_physical_absolute_integral_ball']},
E_1+E_2 <= {c['combined_higher_current_physical_absolute_integral_ball']}
          <1.491e-6.                                  (HC6)
```

The `C_2` cost is negligible; the deliberately coarse absolute enclosure of
`C_1` is below `1.49e-6`.  No oscillatory gain is claimed for either.

This does not certify the remainder after the third face current, the exact
local replacements, or either outside-window tail.  It also must not be
added to the first-current complex modulus: the certified first-current
quantity is its signed physical projection.

Pi provenance: all pi factors are differentiated from the equation-(9)
triangle phase or inherited from its paper normalization.  No fitted
constant is used.

Proof boundary: the second and third pole-subtracted B face currents on the
saved-height `|xi|<=70` window only.  No complete boundary expansion,
complete B estimate, A-fold splice, complete paired residual, complete
`T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["signed_first_current"]["decision"]["first_pole_subtracted_current_integrated_before_absolute_value"] is True, "signed-current dependency drift")
    require(dependencies["triangle"]["decision"]["pair_Volterra_integral_reduced_to_separable_triangle"] is True, "triangle dependency drift")
    require(dependencies["B_window"]["decision"]["xi_70_window_contains_exactly_B_crossings_621_622"] is True, "window dependency drift")

    artifact = {
        "kind": STEM,
        "status": "second_and_third_pole_subtracted_B_face_window_currents_below_1p491e_minus_6_physical",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "second_and_third_face_current_formulas_derived_exactly": True,
            "local_621_622_poles_excluded_before_summation": True,
            "infinite_lattice_tails_bounded_analytically": True,
            "combined_higher_window_current_below_1p491e_minus_6": True,
            "remainder_after_third_face_current_bounded": False,
            "exact_local_621_622_replacements_added": False,
            "outside_window_tail_bounded": False,
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
        "next_obligation": "Derive a uniform remainder after the third face integration by parts for all nonlocal modes and join it to the exact local 621/622 curved-domain replacements. Then combine the signed first current, higher-current bounds, and Fresnel remainders in one B-window budget before treating |xi|>70.",
        "proof_boundary": "Only the second and third nonlocal B face currents on the saved-height |xi|<=70 window. No expansion remainder, local replacement, outside-window, complete B, A-fold, complete paired residual, complete T_upper, height-uniform, Lambda<=0, PF-infinity, RH, or prize-level theorem is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified higher B face window currents below 1.491e-6", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
