#!/usr/bin/env python3
"""Certify the all-order signed-pair endpoint expansion and its route guard."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_all_order_endpoint_jet_route_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__)
CHECKER = Path(__file__).with_name("check_" + STEM + ".py")
DEPENDENCIES = {
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "two_jet_reassembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate.json",
    "A_transition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate.json",
}

A = 159_577
B = 5_122_421
T_LO = 622
T_HI = 39_894


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


def symbolic_certificate() -> dict[str, Any]:
    d = sp.symbols("d", nonzero=True)
    jumps = sp.symbols("j0:8")
    r_plus, r_minus = sp.symbols("r_plus r_minus")
    checked_orders: list[int] = []
    for order in (2, 4, 6, 8):
        plus = -sum(jumps[j] / d ** (j + 1) for j in range(order)) + r_plus / d**order
        minus = -sum(jumps[j] / (-d) ** (j + 1) for j in range(order)) + r_minus / (-d) ** order
        expected = -2 * sum(jumps[2 * k + 1] / d ** (2 * k + 2) for k in range(order // 2))
        expected += (r_plus + r_minus) / d**order
        require(sp.simplify(plus + minus - expected) == 0, f"pair expansion drift at order {order}")
        checked_orders.append(order)

    y, q = sp.symbols("y q")
    polynomial = y
    leading_coefficients: list[str] = []
    for n in range(9):
        poly = sp.Poly(sp.expand(polynomial), y)
        require(poly.degree() == n + 1, f"derivative degree drift at order {n}")
        require(sp.simplify(poly.LC() - q**n) == 0, f"leading coefficient drift at order {n}")
        leading_coefficients.append(str(poly.LC()))
        polynomial = sp.expand(2 * sp.diff(polynomial, y) + q * y * polynomial)

    return {
        "checked_even_orders": checked_orders,
        "all_order_pair_identity": "I_m+I_-m=-2*sum_(k=0)^(r-1) Delta f_x^(2k+1)/(2*pi*i*m)^(2k+2)+(R_(m,2r)+R_(-m,2r))/(2*pi*i*m)^(2r)",
        "remainder_definition": "R_(+/-m,2r)=integral_0^L f_x^(2r)(u)*exp(-/+2*pi*i*m*u)du",
        "remainder_bound": "absolute remainder <=2*integral_0^L abs(f_x^(2r)(u))du/(2*pi*m)^(2r)",
        "endpoint_derivative_recurrence": "p_0(y)=y; p_(n+1)(y)=2*p_n'(y)+q*y*p_n(y); q=i*pi*x",
        "endpoint_derivative_form": "f_x^(n)(u)=p_n(A+2u)*exp(i*pi*x*(A+2u)^2/4)",
        "highest_monomial": "p_n(y)=q^n*y^(n+1)+lower powers",
        "checked_leading_coefficients_n0_to_n8": leading_coefficients,
        "successive_top_odd_jet_ratio": "rho_D(x,m)^2=[x*D/(2*m)]^2",
    }


def exact_turning_geometry() -> dict[str, Any]:
    rho_b_target = Fraction(B, 4 * T_HI)
    rho_a_target = Fraction(A, 4 * T_HI)
    rho_a_outer = Fraction(A, 4 * (T_HI + 1))
    x_b_622 = Fraction(2 * T_LO, B)
    require(rho_b_target > 32, "B target ratio is not above 32")
    require(rho_a_target > 1, "A target edge is not on the non-descending side")
    require(rho_a_outer < 1, "first outer A mode is not beyond the half-boundary")
    require(Fraction(0) < x_b_622 < Fraction(1, 2), "B crossing is outside the half-domain")
    require(x_b_622 * B == 2 * T_LO, "B crossing identity drift")
    return {
        "B_half_cell_ratio_at_target_top": {
            "exact": f"{rho_b_target.numerator}/{rho_b_target.denominator}",
            "decimal": mp.nstr(mp.mpf(rho_b_target.numerator) / rho_b_target.denominator, 24),
            "strict_lower_bound": "32",
        },
        "A_half_cell_ratio_at_target_top": {
            "exact": f"{rho_a_target.numerator}/{rho_a_target.denominator}",
            "decimal": mp.nstr(mp.mpf(rho_a_target.numerator) / rho_a_target.denominator, 24),
            "relation": ">1",
        },
        "A_half_cell_ratio_at_first_outer_mode": {
            "exact": f"{rho_a_outer.numerator}/{rho_a_outer.denominator}",
            "decimal": mp.nstr(mp.mpf(rho_a_outer.numerator) / rho_a_outer.denominator, 24),
            "relation": "<1",
        },
        "B_mode_622_normal_crossing": {
            "x_exact": f"{x_b_622.numerator}/{x_b_622.denominator}",
            "x_decimal": mp.nstr(mp.mpf(x_b_622.numerator) / x_b_622.denominator, 24),
            "rho_B": "1 exactly",
        },
        "target_band_consequence": "At x=1/2 every m in 622..39894 has rho_A>1 and rho_B>32; at x=2m/D the corresponding endpoint ratio is exactly one.",
    }


def endpoint_polynomials(q: mp.mpc, maximum_order: int) -> list[list[mp.mpc]]:
    coefficients = [mp.mpc(0), mp.mpc(1)]
    rows = [coefficients]
    for _ in range(maximum_order):
        derivative = [(j + 1) * coefficients[j + 1] for j in range(len(coefficients) - 1)]
        following = [mp.mpc(0)] * max(len(derivative), len(coefficients) + 1)
        for j, value in enumerate(derivative):
            following[j] += 2 * value
        for j, value in enumerate(coefficients):
            following[j + 1] += q * value
        coefficients = following
        rows.append(coefficients)
    return rows


def evaluate_polynomial(coefficients: list[mp.mpc], y: mp.mpf) -> mp.mpc:
    value = mp.mpc(0)
    for coefficient in reversed(coefficients):
        value = value * y + coefficient
    return value


def diagnostic_scout() -> list[dict[str, Any]]:
    mp.mp.dps = 80
    cases = (
        ("B_normal_crossing_m622", mp.mpf(T_LO), 2 * mp.mpf(T_LO) / B, "B"),
        ("A_half_boundary_m39894", mp.mpf(T_HI), mp.mpf("0.5"), "A"),
        ("A_first_outer_m39895", mp.mpf(T_HI + 1), mp.mpf("0.5"), "A"),
    )
    output: list[dict[str, Any]] = []
    for name, mode, x, endpoint_name in cases:
        q = mp.j * mp.pi * x
        polynomials = endpoint_polynomials(q, 7)
        endpoint = mp.mpf(B if endpoint_name == "B" else A)
        fourier_d = 2 * mp.pi * mp.j * mode
        previous: mp.mpf | None = None
        terms: list[dict[str, str | int | None]] = []
        for n in (1, 3, 5, 7):
            derivative = evaluate_polynomial(polynomials[n], endpoint)
            derivative *= mp.exp(q * endpoint**2 / 4)
            magnitude = abs(2 * derivative / fourier_d ** (n + 1))
            terms.append(
                {
                    "odd_derivative_order": n,
                    "endpoint_term_modulus": mp.nstr(magnitude, 22),
                    "ratio_to_previous": None if previous is None else mp.nstr(magnitude / previous, 22),
                }
            )
            previous = magnitude
        rho = x * endpoint / (2 * mode)
        output.append(
            {
                "case": name,
                "endpoint": endpoint_name,
                "mode": int(mode),
                "x": mp.nstr(x, 24),
                "rho_squared": mp.nstr(rho**2, 22),
                "terms": terms,
            }
        )
    return output


def render_note(geometry: dict[str, Any], scout: list[dict[str, Any]]) -> str:
    b_ratio = geometry["B_half_cell_ratio_at_target_top"]["decimal"]
    a_ratio = geometry["A_half_cell_ratio_at_target_top"]["decimal"]
    a_outer = geometry["A_half_cell_ratio_at_first_outer_mode"]["decimal"]
    b_cross = geometry["B_mode_622_normal_crossing"]["x_decimal"]
    scout_lines: list[str] = []
    for row in scout:
        ratios = [term["ratio_to_previous"] for term in row["terms"] if term["ratio_to_previous"] is not None]
        scout_lines.append(f"{row['case']}: rho^2={row['rho_squared']}; observed ratios={', '.join(str(value) for value in ratios)}")
    telemetry = "\n".join(scout_lines)
    return f"""# All-order signed-pair endpoint jets and the turning-ratio route guard

Date: 2026-08-23

Status: exact all-order pair identity and exact turning-ratio obstruction
certified; raw higher-jet route rejected, phase-adapted joined bound open

Let

```text
f_x(u)=(A+2u)exp[i*pi*x*(A+2u)^2/4],
I_m(x)=integral_0^L f_x(u)exp(-2*pi*i*m*u)du,
D_m=2*pi*i*m.
```

Repeated integration by parts at any even order `2r` and pairing `+m` with
`-m` before taking absolute values gives the exact identity

```text
I_m+I_-m
 =-2 sum_(k=0)^(r-1) Delta f_x^(2k+1)/D_m^(2k+2)
  +[R_(m,2r)+R_(-m,2r)]/D_m^(2r),                    (AJ1)

R_(+/-m,2r)
 =integral_0^L f_x^(2r)(u)exp(-/+2*pi*i*m*u)du.       (AJ2)
```

Thus every even endpoint derivative cancels from the signed pair; only odd
endpoint jumps remain.  The elementary absolute remainder is

```text
|remainder|<=2 integral_0^L |f_x^(2r)(u)|du
                    /(2*pi*m)^(2r).                   (AJ3)
```

Put `y=A+2u` and `q=i*pi*x`.  The endpoint derivatives have the exact
polynomial recurrence

```text
p_0(y)=y,
p_(n+1)(y)=2p_n'(y)+q*y*p_n(y),
f_x^(n)(u)=p_n(y)exp(q*y^2/4),                        (AJ4)

p_n(y)=q^n y^(n+1)+lower powers.                     (AJ5)
```

For endpoint `D` in `{{A,B}}`, the ratio of successive highest odd-jet
monomials is therefore exactly

```text
rho_D(x,m)^2=[x*D/(2m)]^2.                            (AJ6)
```

This is the same normal turning coordinate already resolved by the exact
Fresnel/Morse transition charts.  At a normal endpoint crossing `x=2m/D`,
`rho_D=1`; repeated endpoint extraction is not a descending geometric
hierarchy there.

The saved integers make the obstruction exact.  At `x=1/2` and the top
target mode,

```text
rho_B={b_ratio}>32,
rho_A={a_ratio}>1.                                    (AJ7)
```

For the first outer mode `m=39895`,

```text
rho_A={a_outer}<1,                                    (AJ8)
```

so the A half-boundary lies precisely between modes 39894 and 39895.  The B
normal crossing for mode 622 is at

```text
x=1244/5122421={b_cross},    rho_B=1.                 (AJ9)
```

An 80-digit floating scout, used only to choose the route, evaluates the
separate endpoint terms of (AJ1) through derivative order seven:

```text
{telemetry}
```

The exact conclusion is limited but decisive.  Raising the raw endpoint-jet
order cannot give a uniform descending absolute hierarchy across the B
crossing or the A half-boundary.  This does not rule out a signed resummation,
a contour argument, or the already-certified transition functions.  It says
that the next `R_after_A` estimate must retain those phase-adapted A/B
extractions and act on the remaining joined projector current; it must not
replace them by a global high-order Bernoulli or endpoint-derivative norm.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (AJ1)--(AJ9) is inherited from the Kummer
quadratic phase and integer Fourier character.  No fitted or geometric
occurrence is introduced.

Proof boundary: exact all-order paired integration by parts, derivative
recurrence, highest-monomial ratio, and rational turning geometry only.  The
floating order scout is diagnostic.  No lower bound on the true signed
remainder, no exclusion of cancellation-aware resummation, no quantitative
`R_after_A` or `R_Dir` bound, no complete `Q_K-T`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["symmetric_poisson"]["decision"]["paired_one_over_m_endpoint_current_cancels"] is True, "pair cancellation dependency drift")
    require(dependencies["two_jet_reassembly"]["decision"]["common_kernel_reassembled_as_symmetric_pair_series"] is True, "two-jet dependency drift")
    require(dependencies["A_transition"]["decision"]["projector_completed_exact_A_transition_certified"] is True, "A-transition dependency drift")

    symbolic = symbolic_certificate()
    geometry = exact_turning_geometry()
    scout = diagnostic_scout()
    artifact = {
        "kind": STEM,
        "status": "exact_all_order_signed_pair_endpoint_expansion_and_turning_ratio_route_guard_certified_joined_bound_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "target_positive_modes": [T_LO, T_HI],
            "kummer_half_domain": [0, "1/2"],
        },
        "symbolic_certificate": symbolic,
        "exact_turning_geometry": geometry,
        "diagnostic_order_scout": {
            "rigorous_input_to_theorem": False,
            "precision_decimal_digits": 80,
            "rows": scout,
        },
        "decision": {
            "all_order_signed_pair_identity_proved": True,
            "only_odd_endpoint_jumps_survive_pairing": True,
            "endpoint_derivative_recurrence_proved": True,
            "highest_jet_ratio_proved": True,
            "B_target_half_cell_ratio_above_32": True,
            "A_half_boundary_bracketed_by_39894_39895": True,
            "raw_higher_jet_absolute_hierarchy_uniformly_descending": False,
            "phase_adapted_transition_charts_still_required": True,
            "cancellation_aware_resummation_excluded": False,
            "quantitative_R_after_A_bound_proved": False,
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
            "process_priority": priority,
        },
        "next_obligation": "Keep the certified B and A transition extractions. Derive a phase-adapted evaluator for the remaining joined post-A projector current, preferably by applying the existing Volterra/common-phase transport after algebraic deletion of the A and B atoms, rather than by increasing a global endpoint-jet order.",
        "proof_boundary": "Exact all-order paired integration by parts, derivative recurrence, highest-monomial ratio, and rational turning geometry only. The floating order scout is diagnostic. No lower bound on the true signed remainder, exclusion of cancellation-aware resummation, quantitative R_after_A or R_Dir bound, complete Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(geometry, scout), encoding="utf-8")
    print("certified all-order endpoint-jet route guard", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
