#!/usr/bin/env python3
"""Certify the exact signed-pair logistic-Morse transform on x<=1/2."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_signed_pair_logistic_morse_transform_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "paired_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "half_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
    "logistic_transform": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate.json",
    "universal_morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "midpoint_guard": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_boundary_linear_roster_pair_anchor_gate.json",
}

T = 10_000_000_000
A = 159_577
B = 5_122_421


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
    alpha, x, m, t, endpoint = sp.symbols("alpha x m t D", positive=True, real=True)
    q_plus = sp.sqrt(x / 2) * (alpha - 2 * m / x)
    original_plus = sp.pi * x * alpha**2 / 4 - sp.pi * m * alpha
    completed_plus = -sp.pi * m**2 / x + sp.pi * q_plus**2 / 2
    require(sp.simplify(original_plus - completed_plus) == 0, "positive completion failed")

    q_minus = sp.sqrt(x / 2) * (alpha + 2 * m / x)
    original_minus = sp.pi * x * alpha**2 / 4 + sp.pi * m * alpha
    completed_minus = -sp.pi * m**2 / x + sp.pi * q_minus**2 / 2
    require(sp.simplify(original_minus - completed_minus) == 0, "negative completion failed")
    require(
        sp.simplify(completed_plus.subs(q_plus, q_minus) - completed_minus) != 0,
        "signed Fresnel coordinates were incorrectly identified",
    )

    r, v = sp.symbols("r v", positive=True, real=True)
    logistic_x = 1 / (1 + r * v)
    r_value = t / (2 * sp.pi * m**2)
    psi = -sp.pi * m**2 / logistic_x + t * sp.log((1 - logistic_x) / logistic_x) / 2
    psi_at_one = sp.simplify(psi.subs(v, 1))
    universal = sp.simplify(psi - psi_at_one)
    expected = -t * (v - 1 - sp.log(v)) / 2
    require(
        sp.simplify(sp.expand_log(universal.subs(r, r_value) - expected, force=True)) == 0,
        "universal logistic phase failed",
    )

    v_half = sp.simplify(((1 - x) / (r_value * x)).subs(x, sp.Rational(1, 2)))
    require(sp.simplify(v_half - 2 * sp.pi * m**2 / t) == 0, "half-boundary coordinate failed")
    threshold = sp.sqrt(t / (2 * sp.pi))
    require(sp.simplify(v_half.subs(m, threshold) - 1) == 0, "half-boundary threshold failed")

    # The common outer phase and parity are the mechanism that permits the
    # signed pair to be formed inside one Morse amplitude.
    require(sp.simplify((-sp.pi * (-m) ** 2 / x) - (-sp.pi * m**2 / x)) == 0, "m-squared phase drift")

    delta_f, delta_f1, rem_plus, rem_minus = sp.symbols("Delta_f Delta_f1 R_plus R_minus")
    denominator = 2 * sp.pi * sp.I * m
    i_plus = -delta_f / denominator - delta_f1 / denominator**2 + rem_plus / denominator**2
    i_minus = i_plus.subs({m: -m, rem_plus: rem_minus})
    pair = sp.simplify(i_plus + i_minus)
    require(sp.simplify(sp.diff(pair, delta_f)) == 0, "one-over-m current survived")

    return {
        "signed_Fresnel_coordinate": "q_D^(sigma)(m,x)=sqrt(x/2)[D-2 sigma m/x], sigma in {+1,-1}",
        "signed_current": "P_(sigma m)=[E_B^(sigma)-E_A^(sigma)]/(i*pi)+sigma*m*sqrt(2/x)[F(q_B^(sigma))-F(q_A^(sigma))]",
        "signed_mode": "J_(sigma m)=(-1)^m exp(-i*pi*m^2/x) P_(sigma m)/x",
        "pair_current": "P_m^pair=P_m+P_(-m)",
        "pair_mode": "J_m+J_(-m)=(-1)^m exp(-i*pi*m^2/x)P_m^pair/x",
        "logistic_variables": "r=t/(2*pi*m^2), v=(1-x)/(r*x), x=1/(1+r*v)",
        "Morse_coordinate": "s=sgn(v-1)sqrt(2[v-1-log v])",
        "universal_phase": "psi_m(x)-psi_m(x_m)=-t*s^2/4",
        "half_boundary": "v_H=2*pi*m^2/t, s_H=sgn(v_H-1)sqrt(2[v_H-1-log v_H])",
        "paired_half_transform": "K_m+K_(-m)=(-1)^m exp(i psi_m(x_m)) integral_(s_H)^infinity exp(-i*t*s^2/4) A_m^pair(s)ds",
        "paired_amplitude": "A_m^pair=r*x(s)*(dv/ds)*[x(s)(1-x(s))]^(-1/4)[P_m(x(s))+P_(-m)(x(s))]",
        "orientation": "x:0->1/2 maps s:+infinity->s_H; -dx/ds=r*x^2*dv/ds gives the displayed positive orientation.",
        "pairing_guard": "P_m and P_(-m) have distinct signed Fresnel coordinates but exactly the same outer Morse phase; pair amplitudes, not Fresnel endpoints, are combined.",
    }


def finite_height_ledger() -> dict[str, Any]:
    mp.mp.dps = 80
    t = mp.mpf(T)
    pi = mp.pi

    def values(mode: int) -> dict[str, str | int]:
        v_half = 2 * pi * mode**2 / t
        s_half = mp.sign(v_half - 1) * mp.sqrt(2 * (v_half - 1 - mp.log(v_half)))
        x_mode = 1 / (1 + t / (2 * pi * mode**2))
        return {
            "mode": mode,
            "v_H": mp.nstr(v_half, 70),
            "s_H": mp.nstr(s_half, 70),
            "x_m": mp.nstr(x_mode, 70),
        }

    nt = mp.sqrt(t / (2 * pi))
    rows = [values(mode) for mode in (1, 621, 622, 39_694, 39_894, 39_895)]
    require(mp.mpf(39_894) < nt < mp.mpf(39_895), "threshold interval drift")
    require(all(mp.mpf(row["s_H"]) < 0 for row in rows[:-1]), "interior half-boundary sign drift")
    require(mp.mpf(rows[-1]["s_H"]) > 0, "outer half-boundary sign drift")
    return {
        "height": T,
        "N_t": mp.nstr(nt, 70),
        "N_t_interval": [39_894, 39_895],
        "rows": rows,
        "classification": {
            "m_1_through_39894": "s_H<0, so the half-domain Morse interval contains s=0",
            "m_at_least_39895": "s_H>0, so the half-domain interval starts beyond the Morse saddle",
            "negative_partner": "uses the same s_H and outer Gaussian phase because the transform depends on m^2",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = artifact["finite_height_ledger"]["rows"]
    near = {row["mode"]: row for row in rows}
    return f"""# Exact signed-pair half-domain logistic-Morse transform

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of a quantitative amplitude enclosure

For `m>0`, introduce two genuinely different Fresnel coordinates

```text
q_D^(sigma)(m,x)=sqrt(x/2)(D-2 sigma m/x),
sigma in {{+1,-1}}, D in {{A,B}}.                       (PM1)
```

Completing the square for the `+m` and `-m` Fourier characters gives

```text
P_(sigma m)=[E_B^(sigma)-E_A^(sigma)]/(i*pi)
 +sigma m sqrt(2/x)[F(q_B^(sigma))-F(q_A^(sigma))],

J_(sigma m)=(-1)^m exp(-i*pi*m^2/x)P_(sigma m)/x.     (PM2)
```

Although the Fresnel coordinates in (PM1) are different, the parity and
outer phase in (PM2) are identical.  Hence the symmetric pair is exactly

```text
J_m+J_-m=(-1)^m exp(-i*pi*m^2/x)(P_m+P_-m)/x.         (PM3)
```

Use the global logistic coordinate

```text
r=t/(2*pi*m^2), v=(1-x)/(r*x),
s=sgn(v-1)sqrt(2[v-1-log v]).                          (PM4)
```

The phase remains globally exact,

```text
psi_m(x)-psi_m(x_m)=-t*s^2/4.                         (PM5)
```

At the half-boundary `x=1/2`,

```text
v_H=2*pi*m^2/t,
s_H=sgn(v_H-1)sqrt(2[v_H-1-log v_H]).                 (PM6)
```

Since `x:0->1/2` maps monotonically to `s:+infinity->s_H`, the pair in the
normal form of Section 11.354 has the exact one-phase representation

```text
K_m+K_-m=(-1)^m exp(i psi_m(x_m))
 integral_(s_H)^infinity exp(-i*t*s^2/4)A_m^pair(s)ds,

A_m^pair=r*x*(dv/ds)[x(1-x)]^(-1/4)(P_m+P_-m).        (PM7)
```

This is the phase-adapted cancellation object that was missing from the
tangent B calculation.  It keeps the negative completion in the same
Gaussian integral as its positive partner; the `1/m` endpoint current has
already cancelled before an absolute value is taken.

At `t=10^10`, `sqrt(t/(2*pi))` lies strictly between `39894` and `39895`.
The boundary coordinates are

```text
s_H(39894)={near[39894]['s_H']},
s_H(39895)={near[39895]['s_H']}.                       (PM8)
```

Thus modes `1..39894` integrate across the Gaussian saddle, while every
outer pair `m>=39895` starts on its nonstationary side.  The negative partner
does not create another saddle; it changes only the exact pair amplitude.

Equation (PM7) also explains why midpoint roster parity could not cancel the
Fourier pairs coefficientwise.  The cancellation is distributed through
the signed Fresnel currents along the complete `s` contour.

Pi provenance: every `pi` comes from completing the equation-(9)
Kummer/Fourier phase, the logistic saddle scale, or the exact reflection
phase.  No fitted constant is used.

Proof boundary: exact signed-current algebra, half-domain orientation, and
one-phase pair transform only.  No bound for `A_m^pair`, nonlinear B
crossing, outer-pair sum, or A-fold splice is proved.  No theorem establishing
`T_upper`, `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion follows here.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["paired_residual"]["decision"]["finite_paired_half_domain_target_residual_identity_proved"] is True,
        "paired residual dependency drift",
    )
    require(
        dependencies["logistic_transform"]["decision"]["global_logistic_Morse_phase_is_exact"] is True,
        "logistic phase dependency drift",
    )
    require(
        dependencies["universal_morse"]["decision"]["grouped_boundary_fresnel_current_preserved"] is True,
        "grouped-current dependency drift",
    )
    require(
        dependencies["symmetric_poisson"]["decision"]["paired_one_over_m_endpoint_current_cancels"] is True,
        "pair cancellation dependency drift",
    )
    require(
        dependencies["midpoint_guard"]["decision"]["all_nonzero_symmetric_Poisson_pairs_cancel_at_midpoint"] is False,
        "midpoint guard drift",
    )

    artifact = {
        "kind": STEM,
        "status": "exact_signed_pair_half_domain_logistic_morse_transform_proved_amplitude_bound_open",
        "passed": True,
        "scope": {"height": T, "source_endpoints": [A, B], "mode_condition": "integer m>0"},
        "symbolic_certificate": symbolic_certificate(),
        "finite_height_ledger": finite_height_ledger(),
        "decision": {
            "positive_and_negative_modes_share_outer_Morse_phase": True,
            "signed_Fresnel_coordinates_are_identical": False,
            "signed_pair_formed_inside_one_half_domain_Gaussian_integral": True,
            "half_domain_lower_Morse_limit_derived_exactly": True,
            "outer_positive_pairs_begin_beyond_saddle_from_39895": True,
            "negative_partner_is_second_stationary_branch": False,
            "pair_amplitude_quantitatively_bounded": False,
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
            "sympy_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Analyze the exact pair amplitude A_m^pair in the low/B-crossing and outer ranges. Seek a uniform integration-by-parts or contour majorant on s>=s_H that retains P_m+P_-m, and compare its signed aggregate with the full q-density tangent barrier before touching the fold-owned block.",
        "proof_boundary": "Exact signed-current algebra and half-domain one-phase logistic-Morse pair transform only. No quantitative pair-amplitude bound, nonlinear B crossing, outer-pair sum, A-fold splice, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact signed-pair half-domain logistic-Morse transform", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
