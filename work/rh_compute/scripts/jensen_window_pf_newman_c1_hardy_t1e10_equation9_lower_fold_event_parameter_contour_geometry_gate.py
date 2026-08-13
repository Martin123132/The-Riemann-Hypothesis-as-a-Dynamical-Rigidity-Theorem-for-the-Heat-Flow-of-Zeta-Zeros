#!/usr/bin/env python3
"""Certify the event-parameter range and a common enlarged contour."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_contour_geometry_gate"
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.json"
CONTINUOUS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_continuous_height_cell_gate.json"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 100
C = 159_577
Q = (C - 1) // 4
EVENT_COUNT = 399
Y_MAX = 64
OLD_RADIUS = 9
CONTOUR_RADIUS = 14
CONTOUR_HEIGHT = 1
HORIZONTAL_END = 21
CONNECTOR_CAUCHY_RADIUS = arb("0.225")
HORIZONTAL_CAUCHY_RADIUS = arb("0.3")


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


def mode_for_event(index: int) -> int:
    return Q - index // 2 if index % 2 == 0 else Q + (index + 1) // 2


def symbolic_identities() -> dict[str, str]:
    c, m, beta, h, z, y = sp.symbols("C m beta h z y", real=True)
    epsilon = (4 * m - c) / c
    d = beta * epsilon
    tstar = beta**3
    tau = sp.pi * m * (c - 2 * m)
    tstar_c = sp.pi * c**2 / 8
    require(sp.simplify((tstar_c - tau) - sp.pi * (c - 4 * m) ** 2 / 8) == 0, "event square failed")
    require(sp.simplify((tstar_c - tau) - tstar_c * epsilon**2) == 0, "normalized event square failed")
    require(sp.simplify(beta * d**2 - beta**3 * epsilon**2) == 0, "detuning scale failed")
    lam = (tstar - (tstar - beta * d**2 + h)) / beta
    require(sp.simplify(lam - (d**2 - h / beta)) == 0, "lambda cell law failed")
    canonical = z**3 / 3 - (d**2 - h / beta) * z - (z + d) * y + y**2 / (4 * beta)
    require(sp.simplify(sp.diff(canonical, z) - (z**2 - d**2 + h / beta - y)) == 0, "canonical saddle equation failed")
    return {
        "event_parameter": "epsilon_m=(4m-C)/C",
        "detuning": "d_m=beta*epsilon_m",
        "event_square": "t*-tau_m=pi(C-4m)^2/8=beta*d_m^2",
        "cell_lambda": "lambda(tau_m+h)=d_m^2-h/beta",
        "height_transport": "D_(m,tau_m+h)=exp(i*h*z/beta)D_(m,tau_m)",
        "canonical_saddles": "z=+/-sqrt(d_m^2-h/beta+y)",
        "exact_saddles": "z=+/-beta*acosh(sqrt((beta^2+y+y^2/(4beta^2))/(beta^2-d_m^2+h/beta)))",
    }


def certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    pi = arb.pi()
    c = arb(C)
    tstar = pi * c**2 / 8
    beta = tstar ** (arb(1) / 3)
    kappa = pi / (16 * beta)
    roster = []
    old_failures = []
    for index in range(EVENT_COUNT):
        mode = mode_for_event(index)
        signed_odd = 4 * mode - C
        expected = (-1 if index % 2 == 0 else 1) * (2 * index + 1)
        require(signed_odd == expected, "signed event law failed")
        epsilon = arb(signed_odd) / c
        d = beta * epsilon
        canonical_outer = (d**2 + kappa + Y_MAX).sqrt()
        if canonical_outer.lower() > OLD_RADIUS:
            old_failures.append(index)
        roster.append({
            "event_index": index,
            "mode": mode,
            "signed_odd": signed_odd,
            "epsilon_ball": epsilon.str(PRECISION, more=True),
            "detuning_ball": d.str(PRECISION, more=True),
            "canonical_outer_saddle_ball": canonical_outer.str(PRECISION, more=True),
        })

    d_values = [arb(row["detuning_ball"]) for row in roster]
    d_abs_max = max((abs(value).upper() for value in d_values), key=float)
    lambda_abs_max = d_abs_max**2 + kappa
    canonical_saddle_max = (lambda_abs_max + Y_MAX).sqrt()
    exact_ratio = (
        beta**2 + Y_MAX + arb(Y_MAX) ** 2 / (4 * beta**2)
    ) / (beta**2 - d_abs_max**2 - kappa)
    exact_saddle_max = beta * exact_ratio.sqrt().acosh()

    r = arb(CONTOUR_RADIUS)
    v = arb(CONTOUR_HEIGHT)
    canonical_path_bracket = r**2 - v**2 / 3 - lambda_abs_max - Y_MAX
    scaled_gap = (r**2 - v**2) / beta**2
    exact_path_bracket = (
        beta**2 * scaled_gap / (1 + scaled_gap)
        - lambda_abs_max
        - Y_MAX
        - arb(Y_MAX) ** 2 / (4 * beta**2)
    )

    connector_real = r - CONNECTOR_CAUCHY_RADIUS
    connector_imag = arb("0.875") + CONNECTOR_CAUCHY_RADIUS
    horizontal_real = r + arb("0.25") - HORIZONTAL_CAUCHY_RADIUS
    horizontal_imag = v + HORIZONTAL_CAUCHY_RADIUS

    def cauchy_brackets(real_lower: arb, imag_upper: arb) -> tuple[arb, arb]:
        canonical = real_lower**2 - imag_upper**2 / 3 - lambda_abs_max - Y_MAX
        gap = (real_lower**2 - imag_upper**2) / beta**2
        exact = beta**2 * gap / (1 + gap) - lambda_abs_max - Y_MAX - arb(Y_MAX) ** 2 / (4 * beta**2)
        return canonical, exact

    connector_canonical, connector_exact = cauchy_brackets(connector_real, connector_imag)
    horizontal_canonical, horizontal_exact = cauchy_brackets(horizontal_real, horizontal_imag)
    minimum_cauchy_bracket = min(connector_canonical, connector_exact, horizontal_canonical, horizontal_exact)

    end = arb(HORIZONTAL_END)
    exact_far_exponent_margin = end**2 / 2 - (
        Y_MAX + arb(1) / 2 + arb(Y_MAX) ** 2 / (4 * beta**2) + lambda_abs_max
    )
    canonical_far_exponent_margin = end**2 - (Y_MAX + arb(1) / 3 + lambda_abs_max)
    maximum_local_q_radius = ((end + HORIZONTAL_CAUCHY_RADIUS) ** 2 + horizontal_imag**2).sqrt() / beta
    cosh_modulus_lower = (horizontal_imag / beta).cos()

    require(canonical_saddle_max < r, "canonical saddle escaped R=14")
    require(exact_saddle_max < r, "exact saddle escaped R=14")
    require(canonical_path_bracket > 0 and exact_path_bracket > 0, "lifted path lost decay")
    require(minimum_cauchy_bracket > 0, "lifted Cauchy disk lost decay")
    require(exact_far_exponent_margin > 0 and canonical_far_exponent_margin > 0, "far-tail exponent lost decay")
    require(maximum_local_q_radius < arb("0.02"), "local tanh chart radius failed")
    require(cosh_modulus_lower > arb("0.999"), "cosh zero-free margin failed")
    require(len(old_failures) > 0, "old contour unexpectedly covers every event")

    return {
        "C": C,
        "event_count": EVENT_COUNT,
        "mode_range": [min(row["mode"] for row in roster), max(row["mode"] for row in roster)],
        "signed_odd_range": [min(row["signed_odd"] for row in roster), max(row["signed_odd"] for row in roster)],
        "beta_ball": beta.str(PRECISION, more=True),
        "cell_kappa_pi_over_16beta_ball": kappa.str(PRECISION, more=True),
        "maximum_abs_detuning_ball": d_abs_max.str(PRECISION, more=True),
        "maximum_abs_lambda_on_cells_ball": lambda_abs_max.str(PRECISION, more=True),
        "maximum_canonical_saddle_abs_ball": canonical_saddle_max.str(PRECISION, more=True),
        "maximum_exact_saddle_abs_ball": exact_saddle_max.str(PRECISION, more=True),
        "old_R9_failure_count": len(old_failures),
        "old_R9_first_failure_event": old_failures[0],
        "old_R9_last_failure_event": old_failures[-1],
        "common_contour": {
            "radius": CONTOUR_RADIUS,
            "height": CONTOUR_HEIGHT,
            "horizontal_end": HORIZONTAL_END,
            "canonical_saddle_clearance_ball": (r - canonical_saddle_max).str(PRECISION, more=True),
            "exact_saddle_clearance_ball": (r - exact_saddle_max).str(PRECISION, more=True),
            "canonical_lifted_path_bracket_ball": canonical_path_bracket.str(PRECISION, more=True),
            "exact_lifted_path_bracket_ball": exact_path_bracket.str(PRECISION, more=True),
            "minimum_lifted_Cauchy_bracket_ball": minimum_cauchy_bracket.str(PRECISION, more=True),
            "exact_far_exponent_margin_ball": exact_far_exponent_margin.str(PRECISION, more=True),
            "canonical_far_exponent_margin_ball": canonical_far_exponent_margin.str(PRECISION, more=True),
            "maximum_local_tanh_q_radius_ball": maximum_local_q_radius.str(PRECISION, more=True),
            "cosh_modulus_lower_ball": cosh_modulus_lower.str(PRECISION, more=True),
        },
        "roster": roster,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    contour = c["common_contour"]
    return f"""# Event-parameter contour geometry

Date: 2026-08-12

Status: exact 399-event parameter range and one common enlarged contour
certified; not a proof of a uniform full-line remainder estimate

For event index `0<=n<=398`,

```text
4m_n-C=(-1)^(n+1)(2n+1),
epsilon_n=(4m_n-C)/C,       d_n=beta epsilon_n.       (EPC1)
```

The event square and every `pi/16` height cell give

```text
t*-tau_n=beta d_n^2,
lambda(tau_n+h)=d_n^2-h/beta,
D_(n,tau_n+h)=exp(i h z/beta)D_(n,tau_n).             (EPC2)
```

The full roster has modes `{c['mode_range'][0]}..{c['mode_range'][1]}` and
`{c['signed_odd_range'][0]}<=4m-C<={c['signed_odd_range'][1]}`.  Its largest
absolute detuning is

```text
{c['maximum_abs_detuning_ball']}.                     (EPC3)
```

The canonical and exact event-cell saddle envelopes are, respectively,

```text
|z_can|<={c['maximum_canonical_saddle_abs_ball']},
|z_exact|<={c['maximum_exact_saddle_abs_ball']}.       (EPC4)
```

Therefore the prototype `R=9` contour is not atlas-uniform: it misses the
outer saddle envelope for {c['old_R9_failure_count']} events, beginning at
event {c['old_R9_first_failure_event']}.  A fixed contour with vertical
connectors at `Re z=+/-14`, lift `Im z=1`, and signed horizontal integration
through `|Re z|=21` instead has certified saddle clearances

```text
canonical: {contour['canonical_saddle_clearance_ball']},
exact:     {contour['exact_saddle_clearance_ball']}.   (EPC5)
```

The minimum factored phase bracket on every lifted Cauchy disk is
`{contour['minimum_lifted_Cauchy_bracket_ball']}>0`.  At `|Re z|=21`, the
exact logistic-Gaussian and canonical Gaussian exponent margins are
`{contour['exact_far_exponent_margin_ball']}` and
`{contour['canonical_far_exponent_margin_ball']}`.  The entire contour stays
inside the local hyperbolic chart with `|q|<0.02` and a zero-free cosh
margin.

This establishes common contour geometry only.  It does not show that the
signed exact-minus-canonical integral stays within the prototype numerical
budget at extreme detuning.  It does not propagate all 399 event remainders.
It does not establish complete `T_upper` or prove `Lambda<=0` or RH.  It does
not establish a prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for dependency in (ATLAS, CONTINUOUS, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    artifact = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "event_parameter_range_and_common_contour_geometry_certified",
        "passed": True,
        "symbolic_identities": symbolic_identities(),
        "certificate": certificate(),
        "claims": {
            "exact_399_event_parameter_roster_proved": True,
            "prototype_R9_atlas_uniform": False,
            "common_R14_lifted_contour_geometry_proved": True,
            "uniform_full_line_remainder_proved": False,
            "all_399_event_propagation_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": (
            "Exact event-parameter, saddle-envelope, and common-contour geometry only. No uniform full-line remainder, "
            "399-event propagation theorem, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
        "next_action": "Run the signed full-line enclosure at both extreme detunings on the certified R=14 contour before choosing a global interval or paired-event route.",
        "dependencies": {
            "event_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "continuous_prototype": {"path": relative(CONTINUOUS), "sha256": file_hash(CONTINUOUS)},
        },
        "artifacts": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"workers": 1, "process_priority": priority, "elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    artifact["artifacts"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(
        f"certified event contour geometry: R9 failures={artifact['certificate']['old_R9_failure_count']}, "
        f"R14 saddle max={artifact['certificate']['maximum_canonical_saddle_abs_ball']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
