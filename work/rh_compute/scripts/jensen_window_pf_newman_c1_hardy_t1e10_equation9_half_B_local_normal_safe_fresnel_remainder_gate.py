#!/usr/bin/env python3
"""Certify the normalized Fresnel remainders on the two B-local normal-safe arcs."""

from __future__ import annotations

from fractions import Fraction
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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_normal_safe_fresnel_remainder_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "local_gap": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_exact_quadratic_domain_gap_gate.json",
    "Morse_normalization": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
    "half_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
}

PRECISION = 110
T = 10_000_000_000
B = 5_122_421
SPLITS = {
    621: Fraction("0.00024238523659"),
    622: Fraction("0.000242916437095"),
}


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


def arb_fraction(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def symbolic_certificate() -> dict[str, str]:
    x, endpoint, mode, d = sp.symbols("x B m d", positive=True, real=True)
    weight = x ** sp.Rational(7, 4) * (1 - x) ** (-sp.Rational(1, 4))
    derivative = sp.factor(sp.diff(weight, x))
    require(
        sp.simplify(derivative - weight * (7 - 6 * x) / (4 * x * (1 - x))) == 0,
        "transformed-weight derivative failed",
    )
    remainder = 15 * mode * x**2 / (4 * sp.pi**4 * d**7)
    weighted = sp.simplify((x * (1 - x)) ** (-sp.Rational(1, 4)) * remainder)
    expected = 15 * mode * weight / (4 * sp.pi**4 * d**7)
    require(sp.simplify(weighted - expected) == 0, "weighted remainder failed")
    require(sp.integrate(d**-7, d) == -1 / (6 * d**6), "d^-7 primitive failed")

    normalizer = (sp.pi / (32 * sp.symbols("t", positive=True))) ** sp.Rational(1, 4)
    return {
        "completed_remainder": "|R_ext,m(x)|<=15*m*x^2/[4*pi^4*|Bx/2-m|^7]",
        "Kummer_weight": "|W_t(x)|=[x(1-x)]^(-1/4)",
        "weighted_remainder": "15*m*x^(7/4)*(1-x)^(-1/4)/[4*pi^4*d^7], d=|Bx/2-m|",
        "weight_monotonicity": "h(x)=x^(7/4)(1-x)^(-1/4) is strictly increasing on 0<x<=1/2 because h'/h=(7-6x)/[4x(1-x)]>0",
        "distance_change": "dx=2*dd/B on either monotone crossing side",
        "panel_primitive": "integral_d0^d1 d^(-7) dd=(d0^(-6)-d1^(-6))/6",
        "paper_normalizer": str(normalizer),
        "physical_projection_factor": "2*(pi/(32t))^(1/4)",
        "pi_provenance": "pi comes from the equation-(9) Fresnel phase and its paper normalization; Arb pi is used directly",
    }


def panel_bound(mode: int, split: Fraction) -> tuple[arb, list[dict[str, str]]]:
    endpoint = Fraction(B)
    if mode == 621:
        distance = Fraction(mode) - endpoint * split / 2
        distance_max = Fraction(mode)
        side = "left"
    else:
        distance = endpoint * split / 2 - Fraction(mode)
        distance_max = endpoint / 4 - Fraction(mode)
        side = "right"
    require(0 < distance < distance_max, f"invalid distance interval for mode {mode}")

    panels: list[dict[str, str]] = []
    total = arb(0)
    current = distance
    while current < distance_max:
        following = min(2 * current, distance_max)
        if side == "left":
            x_for_weight = 2 * (Fraction(mode) - current) / endpoint
        else:
            x_for_weight = 2 * (Fraction(mode) + following) / endpoint
        x_ball = arb_fraction(x_for_weight)
        d0 = arb_fraction(current)
        d1 = arb_fraction(following)
        h_upper = x_ball ** (arb(7) / 4) * (1 - x_ball) ** (-arb(1) / 4)
        inverse_power_integral = (d0**-6 - d1**-6) / 6
        contribution = h_upper * inverse_power_integral
        total += contribution
        panels.append(
            {
                "d_left": str(current),
                "d_right": str(following),
                "x_weight_endpoint": str(x_for_weight),
                "h_upper_ball": h_upper.str(PRECISION, more=True),
                "dimensionless_panel_ball": contribution.str(PRECISION, more=True),
            }
        )
        current = following

    pi = arb.pi()
    raw_bound = arb(15 * mode) * total / (2 * pi**4 * B)
    return raw_bound, panels


def interval_certificate(local_gap: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    normalizer = (pi / (32 * arb(T))) ** (arb(1) / 4)
    physical_factor = 2 * normalizer
    modes: dict[str, Any] = {}
    combined = arb(0)

    for mode in (621, 622):
        split = SPLITS[mode]
        bracket = local_gap["interval_certificate"][f"mode_{mode}"]["x_bracket"]
        require(Fraction(bracket[0]) < split < Fraction(bracket[1]), f"mode {mode} split escaped balance bracket")
        x = arb_fraction(split)
        q = (arb(B) * x - 2 * arb(mode)) / (2 * x).sqrt()
        if mode == 621:
            require(q < -arb("18.20"), "mode-621 split q margin failed")
            arc = "0<x<=x_split"
            monotone_q_guard = "q_B is increasing and negative at x_split, hence |q_B(x)|>=|q_B(x_split)| on this arc"
        else:
            require(q > arb("14.52"), "mode-622 split q margin failed")
            arc = "x_split<=x<=1/2"
            monotone_q_guard = "q_B is increasing and positive at x_split, hence q_B(x)>=q_B(x_split) on this arc"

        raw, panels = panel_bound(mode, split)
        physical = physical_factor * raw
        combined += physical
        limit = arb("4.02e-11") if mode == 621 else arb("1.56e-10")
        require(physical < limit, f"mode-{mode} normalized remainder exceeds limit")
        modes[str(mode)] = {
            "mode": mode,
            "arc": arc,
            "split_x_exact": str(split),
            "certified_balance_bracket": bracket,
            "q_at_split_ball": q.str(PRECISION, more=True),
            "monotone_q_guard": monotone_q_guard,
            "dyadic_panel_count": len(panels),
            "raw_half_integral_remainder_ball": raw.str(PRECISION, more=True),
            "physically_normalized_remainder_ball": physical.str(PRECISION, more=True),
            "panels": panels,
        }

    require(combined < arb("1.97e-10"), "combined normalized local normal-safe remainder exceeds limit")
    return {
        "height": T,
        "endpoint": B,
        "paper_normalizer_ball": normalizer.str(PRECISION, more=True),
        "physical_projection_factor_ball": physical_factor.str(PRECISION, more=True),
        "modes": modes,
        "combined_physically_normalized_remainder_ball": combined.str(PRECISION, more=True),
        "combined_physically_normalized_remainder_upper_bound": "1.97e-10",
        "reference_physical_target": "8.6e-6",
        "target_fraction_upper_bound": "2.30e-5",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    m621 = c["modes"]["621"]
    m622 = c["modes"]["622"]
    return f"""# Normal-safe Fresnel remainder for local B modes

Date: 2026-08-13

Status: saved-height interval certificate for two specified arcs; not a
proof of the complete local B current

The exact local-domain gate locates the balance points between the outer
Gaussian and normal Fresnel derivatives.  Choose the exact rational split
points

```text
x_621={m621['split_x_exact']} inside {m621['certified_balance_bracket']},
x_622={m622['split_x_exact']} inside {m622['certified_balance_bracket']}. (NF1)
```

On `0<x<=x_621`, `q_B(621,x)` is negative and increasing, while on
`x_622<=x<=1/2`, `q_B(622,x)` is positive and increasing.  The certified
split values are

```text
q_621={m621['q_at_split_ball']},
q_622={m622['q_at_split_ball']}.                       (NF2)
```

Thus the third-order Fresnel expansion is uniformly valid on both arcs.
Its completed-current remainder from Section 11.364 is

```text
|R_ext,m(x)|<=15*m*x^2/[4*pi^4*|Bx/2-m|^7].          (NF3)
```

Restoring the exact Kummer weight and the real half-domain projection gives

```text
E_m^N=2*(pi/(32t))^(1/4)
 integral_arc [x(1-x)]^(-1/4)|R_ext,m(x)|dx.          (NF4)
```

This is the normalization used by equation (9), not the raw completed
current.  Put `d=|Bx/2-m|` and
`h(x)=x^(7/4)(1-x)^(-1/4)`.  Since

```text
h'(x)/h(x)=(7-6x)/[4x(1-x)]>0 on 0<x<=1/2,           (NF5)
```

each dyadic `d` panel is bounded at the correct monotone endpoint, while
`integral d^(-7)dd` is evaluated exactly.  Arb outward rounding gives

```text
E_621^N <= {m621['physically_normalized_remainder_ball']},
E_622^N <= {m622['physically_normalized_remainder_ball']},

E_621^N+E_622^N
 <= {c['combined_physically_normalized_remainder_ball']} < 1.97e-10. (NF6)
```

This uses only {m621['dyadic_panel_count']} and {m622['dyadic_panel_count']}
analytic panels respectively; there is no oscillatory numerical quadrature.
The combined bound is below `2.30e-5` of the reference physical target
`8.6e-6`.

The result closes only the positive local B-tail remainders on the two
normal-safe arcs.  On the complementary outer-safe arcs, the exact 621/622
currents must remain combined with the pole-subtracted cotangent background;
its raw absolute value is known to be too large to triangle away.  Negative,
zero, outer-positive, and A-fold sectors are not reassigned here.

Pi provenance: every `pi` in (NF3)--(NF6) comes from the equation-(9)
Fresnel phase or the paper normalization `(pi/(32t))^(1/4`; Arb evaluates
that same constant directly.  No geometric or fitted value is inserted.

Proof boundary: the two stated saved-height positive-B normal-safe remainder
arcs only.  No complementary outer-safe local-current estimate, complete B
face estimate, A-fold splice, complete paired residual, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["B_window"]["decision"]["three_term_Fresnel_tail_with_explicit_remainder_proved"] is True,
        "Fresnel-remainder dependency drift",
    )
    require(
        dependencies["local_gap"]["decision"]["local_exact_max_gradient_component_above_45p63"] is True,
        "local-gap dependency drift",
    )
    require(
        dependencies["Morse_normalization"]["decision"]["classical_inverse_sqrt_mode_carrier_recovered"] is True,
        "paper-normalization dependency drift",
    )
    require(
        dependencies["half_reflection"]["decision"]["full_source_main_equals_twice_real_half_integral"] is True,
        "half-reflection dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "two_local_B_normal_safe_positive_tail_remainders_below_1p97e_minus_10_in_physical_normalization",
        "passed": True,
        "scope": {
            "height": T,
            "endpoint": B,
            "modes": [621, 622],
            "mode_621_arc": "0<x<=0.00024238523659",
            "mode_622_arc": "0.000242916437095<=x<=1/2",
        },
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(dependencies["local_gap"]),
        "decision": {
            "exact_equation9_normalization_restored": True,
            "mode_621_normal_safe_positive_B_tail_remainder_bounded": True,
            "mode_622_normal_safe_positive_B_tail_remainder_bounded": True,
            "combined_physical_remainder_below_1p97e_minus_10": True,
            "oscillatory_quadrature_used": False,
            "complementary_outer_safe_local_current_bounded": False,
            "negative_zero_outer_positive_completion_bounded": False,
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
        "next_obligation": "On the complementary outer-safe arcs, form the exact smooth sum of the 621/622 crossing currents and the pole-subtracted first, third, and fifth rational cotangent currents in one outer-Morse amplitude. Bound that signed amplitude and its derivative before outer Gaussian integration by parts.",
        "proof_boundary": "Only the positive local B-tail Fresnel remainders on two explicitly stated saved-height normal-safe arcs. No complementary outer-safe current estimate, negative/zero/outer-positive completion, complete B estimate, A-fold splice, complete paired residual, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified normalized local B normal-safe Fresnel remainders below 1.97e-10", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
