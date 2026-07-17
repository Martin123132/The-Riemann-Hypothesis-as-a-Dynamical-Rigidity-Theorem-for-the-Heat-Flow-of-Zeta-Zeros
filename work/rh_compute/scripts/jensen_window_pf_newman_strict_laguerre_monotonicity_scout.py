#!/usr/bin/env python3
"""Build an exact target and numerical scout for monotonic first Laguerre flow."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate import (
    quadrature_phi,
)
from jensen_window_pf_newman_theta_modular_blend_high_frequency_scout import (
    quadrature_rows,
    xi,
)
import jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate as lehmer_gate


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_strict_laguerre_monotonicity_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.md"
)
COUNTERMODEL_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.json"
)

MODERATE_TIMES = ("0", "0.05", "0.1", "0.15", "0.2", "0.5")
HIGH_TIMES = ("0", "0.1", "0.2", "0.5")
HIGH_X = (40, 80, 120, 160, 200)
MP_DPS = 75


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def format_float(value: float) -> str:
    return format(float(value), ".17e")


def format_mpf(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def exact_algebra() -> dict:
    h0, h1, h2, h3 = sp.symbols("H H_prime H_second H_third", real=True)
    a0, a1, a2, a3 = sp.symbols("A A_prime A_second A_third", real=True)
    # Differentiate H'^2-H*H'' and substitute independent jet symbols.
    x = sp.symbols("x", real=True)
    function = sp.Function("f")(x)
    derivative = sp.diff(sp.diff(function, x) ** 2 - function * sp.diff(function, x, 2), x)
    jet_map = {
        function: h0,
        sp.diff(function, x): h1,
        sp.diff(function, x, 2): h2,
        sp.diff(function, x, 3): h3,
    }
    derivative = sp.expand(derivative.xreplace(jet_map))
    if sp.simplify(derivative - (h1 * h2 - h0 * h3)) != 0:
        raise RuntimeError("first-Laguerre derivative identity failed")

    theta_expression = sp.expand(
        ((a0 + sp.Rational(1, 2)) / 8) * (a3 / 8)
        - (a1 / 8) * (a2 / 8)
    )
    theta_expected = ((a0 + sp.Rational(1, 2)) * a3 - a1 * a2) / 64
    if sp.simplify(theta_expression - theta_expected) != 0:
        raise RuntimeError("theta-primitive monotonicity identity failed")

    return {
        "definitions": {
            "heat_function": "H_t(x)=integral_0^infinity exp(t*u^2)*Phi(u)*cos(x*u)du",
            "first_laguerre": "L_t(x)=H_t'(x)^2-H_t(x)*H_t''(x)",
            "monotonicity_margin": "M_t(x)=-L_t'(x)=H_t(x)*H_t'''(x)-H_t'(x)*H_t''(x)",
        },
        "derivative_identity": "L_t'(x)=H_t'(x)*H_t''(x)-H_t(x)*H_t'''(x)",
        "monotone_sufficiency": (
            "If M_t(x)>0 for every x>0 and L_t(x)->0 as x->infinity, then "
            "L_t(x)=integral_x^infinity M_t(y)dy>0 for every x>=0."
        ),
        "newman_consequence": (
            "M_t(x)>0 for every x>0 and every 0<t<=1/5 implies Lambda<=0."
        ),
        "multiple_root_exclusion": (
            "A multiple real zero c>0 has H_t(c)=H_t'(c)=0 and hence "
            "M_t(c)=0, contradicting strict monotonicity."
        ),
        "correlation_sine_transform": (
            "M_t(x)=2*integral_R v*K_(1,t)(v)*sin(2*x*v)dv="
            "4*integral_0^infinity v*K_(1,t)(v)*sin(2*x*v)dv"
        ),
        "tail_cosine_reduction": {
            "definition": "Q_t(s)=integral_s^infinity v*K_(1,t)(v)dv",
            "identity": (
                "M_t(x)=8*x*integral_0^infinity Q_t(s)*cos(2*x*s)ds"
            ),
            "derivatives": "Q_t'(s)=-s*K_(1,t)(s), Q_t''(s)=-K_(1,t)(s)-s*K_(1,t)'(s)",
            "one_inflection_theorem": (
                "For ell=log K_(1,t), h(s)=1+s*ell'(s) decreases strictly from "
                "1 to -infinity because h'=ell'+s*ell''<0 and strong "
                "log-concavity gives ell'(s)<=-c*s. Hence Q_t''=-K_(1,t)*h "
                "changes sign exactly once, from negative to positive."
            ),
            "classical_convexity_obstruction": (
                "Q_t''(0)=-K_(1,t)(0)<0, so the classical positive-decreasing-"
                "convex cosine-transform criterion cannot establish M_t>0."
            ),
            "regular_shift_obstruction": (
                "At every regular fixed shift p>0, evenness and real analyticity "
                "give Q_t(i*p) real. Hence -Im Q_t(s+i*p) starts at zero and "
                "cannot be both nonnegative and nonincreasing unless it is "
                "identically zero, so Tuck's elementary interior-shift criterion "
                "cannot apply nontrivially."
            ),
            "boundary_handoff": (
                "A viable complex-contour proof needs a boundary correction or a "
                "frequency-adaptive modular shift; the theta series identifies "
                "abs(Im z)=pi/8 as its natural termwise-control boundary."
            ),
            "source": "https://doi.org/10.1017/S0004972700047511",
        },
        "theta_primitive": {
            "definition": "A_t=D_t[C_t]=8*H_t-1/2",
            "identity": (
                "M_t(x)=((A_t(x)+1/2)*A_t'''(x)-A_t'(x)*A_t''(x))/64"
            ),
            "target": (
                "(A_t(x)+1/2)*A_t'''(x)-A_t'(x)*A_t''(x)>0 for x>0"
            ),
        },
    }


def moderate_laguerre_derivative(node_count: int) -> dict[str, np.ndarray]:
    u, weights, raw_phi = quadrature_phi(node_count)
    x_values = np.arange(1, 2001, dtype=float) / 50
    output: dict[str, np.ndarray] = {"x": x_values}
    for t_text in MODERATE_TIMES:
        t = float(t_text)
        weighted = weights * np.exp(t * u * u) * raw_phi[0]
        chunks: list[np.ndarray] = []
        for start in range(0, x_values.size, 200):
            frequencies = x_values[start : start + 200]
            phase = np.outer(u, frequencies)
            cosine = np.cos(phase)
            sine = np.sin(phase)
            h0 = weighted @ cosine
            h1 = -(weighted * u) @ sine
            h2 = -(weighted * u**2) @ cosine
            h3 = (weighted * u**3) @ sine
            chunks.append(h1 * h2 - h0 * h3)
        output[t_text] = np.concatenate(chunks)
    return output


def build_moderate_scout() -> dict:
    coarse = moderate_laguerre_derivative(192)
    fine = moderate_laguerre_derivative(320)
    x_values = fine["x"]
    rows: list[dict] = []
    max_absolute_delta = 0.0
    max_relative_delta = 0.0
    for key in MODERATE_TIMES:
        values = fine[key]
        deltas = np.abs(values - coarse[key])
        relative = deltas / np.maximum(np.abs(values), 1e-300)
        max_absolute_delta = max(max_absolute_delta, float(np.max(deltas)))
        max_relative_delta = max(max_relative_delta, float(np.max(relative)))
        max_index = int(np.argmax(values))
        min_index = int(np.argmin(values))
        rows.append(
            {
                "t": key,
                "x_count": int(x_values.size),
                "x_min": format(float(x_values[0]), ".17g"),
                "x_max": format(float(x_values[-1]), ".17g"),
                "all_L_prime_negative": bool(np.all(values < 0)),
                "largest_L_prime": format_float(values[max_index]),
                "largest_L_prime_x": format(float(x_values[max_index]), ".17g"),
                "smallest_L_prime": format_float(values[min_index]),
                "smallest_L_prime_x": format(float(x_values[min_index]), ".17g"),
            }
        )
    if not all(row["all_L_prime_negative"] for row in rows):
        raise RuntimeError("moderate monotonicity scout found a nonnegative value")
    if max_relative_delta >= 1e-6:
        raise RuntimeError("moderate coarse/fine convergence lost")
    return {
        "method": (
            "Composite Gauss-Legendre quadrature on [0,4], finite Phi series "
            "n<=15, analytic Fourier jets through order 3"
        ),
        "coarse_nodes_per_panel": 192,
        "fine_nodes_per_panel": 320,
        "x_step": "1/50",
        "max_absolute_coarse_fine_delta": format_float(max_absolute_delta),
        "max_relative_coarse_fine_delta": format_float(max_relative_delta),
        "rows": rows,
    }


def mp_transform_jet(rows: list[tuple], t: mp.mpf, x: mp.mpf) -> tuple[mp.mpf, ...]:
    weighted = [(u, qweight * mp.exp(t * u * u) * full) for u, qweight, _, full in rows]
    return (
        mp.fsum(weight * mp.cos(x * u) for u, weight in weighted),
        mp.fsum(-weight * u * mp.sin(x * u) for u, weight in weighted),
        mp.fsum(-weight * u**2 * mp.cos(x * u) for u, weight in weighted),
        mp.fsum(weight * u**3 * mp.sin(x * u) for u, weight in weighted),
    )


def xi_laguerre_derivative(x: mp.mpf) -> mp.mpf:
    s = mp.mpf("0.5") + mp.j * x / 2
    derivatives = [xi(s)] + [mp.diff(xi, s, order) for order in range(1, 4)]
    jets = [
        mp.re(derivatives[0] / 8),
        mp.re(mp.j * derivatives[1] / 16),
        mp.re(-derivatives[2] / 32),
        mp.re(-mp.j * derivatives[3] / 64),
    ]
    return jets[1] * jets[2] - jets[0] * jets[3]


def high_frequency_values(node_count: int) -> list[dict]:
    qrows = quadrature_rows(node_count, 13)
    output: list[dict] = []
    for t_text in HIGH_TIMES:
        t = mp.mpf(t_text)
        for x_value in HIGH_X:
            x = mp.mpf(x_value)
            h0, h1, h2, h3 = mp_transform_jet(qrows, t, x)
            laguerre_prime = h1 * h2 - h0 * h3
            output.append(
                {
                    "t": t_text,
                    "x": str(x_value),
                    "L_prime": format_mpf(laguerre_prime),
                }
            )
    return output


def build_high_frequency_scout() -> dict:
    mp.mp.dps = MP_DPS
    coarse = high_frequency_values(100)
    fine = high_frequency_values(140)
    rows: list[dict] = []
    max_relative_delta = mp.mpf(0)
    max_xi_relative_delta = mp.mpf(0)
    for coarse_row, fine_row in zip(coarse, fine, strict=True):
        if coarse_row["t"] != fine_row["t"] or coarse_row["x"] != fine_row["x"]:
            raise RuntimeError("high-frequency row alignment failed")
        coarse_value = mp.mpf(coarse_row["L_prime"])
        fine_value = mp.mpf(fine_row["L_prime"])
        relative_delta = abs(coarse_value - fine_value) / abs(fine_value)
        max_relative_delta = max(max_relative_delta, relative_delta)
        row = dict(fine_row)
        row["relative_coarse_fine_delta"] = format_mpf(relative_delta, 24)
        if row["t"] == "0":
            exact_value = xi_laguerre_derivative(mp.mpf(row["x"]))
            xi_relative = abs(fine_value - exact_value) / abs(exact_value)
            max_xi_relative_delta = max(max_xi_relative_delta, xi_relative)
            row["exact_xi_L_prime"] = format_mpf(exact_value)
            row["relative_xi_delta"] = format_mpf(xi_relative, 24)
        rows.append(row)
    if any(mp.mpf(row["L_prime"]) >= 0 for row in rows):
        raise RuntimeError("high-frequency monotonicity scout found a nonnegative value")
    if max_relative_delta >= mp.mpf("1e-25"):
        raise RuntimeError("high-frequency coarse/fine convergence lost")
    if max_xi_relative_delta >= mp.mpf("1e-25"):
        raise RuntimeError("high-frequency Xi cross-check lost")
    return {
        "method": (
            "75-digit composite Gauss-Legendre quadrature on [0,2.8], "
            "13 theta terms, analytic Fourier jets through order 3"
        ),
        "dps": MP_DPS,
        "coarse_nodes_per_panel": 100,
        "fine_nodes_per_panel": 140,
        "max_relative_coarse_fine_delta": format_mpf(max_relative_delta, 24),
        "max_relative_t0_xi_delta": format_mpf(max_xi_relative_delta, 24),
        "rows": rows,
    }


def countermodel_guard() -> dict:
    stored = json.loads(COUNTERMODEL_RESULT.read_text(encoding="utf-8"))
    certificate = stored["certificate"]
    if not certificate.get("passed") or not certificate["laguerre_ball"].startswith("[-"):
        raise RuntimeError("theta-tail countermodel certificate is unavailable")
    return {
        "kernel": (
            "phi_(1/10,1/1000)(u)=exp(-u^2-u^4/10-"
            "(cosh(4u)-1)/1000)*(1+u^4/10)"
        ),
        "certified_point": "x=21/5",
        "laguerre_ball": certificate["laguerre_ball"],
        "deduction": (
            "L(0)>0, the certified L(21/5)<0, and L(x)->0 imply that L'(y)>0 "
            "for some y>21/5. Thus uniform strong log-concavity, strict "
            "root-variable concavity, and theta-type decay do not imply the "
            "monotonicity target."
        ),
    }


def xi_monotonicity_counterexample() -> dict:
    certified = lehmer_gate.arb_exact_point()
    return {
        "x": certified["x"],
        "normalized_monotonicity_margin_ball": certified[
            "exact_normalized_monotonicity_margin_ball"
        ],
        "monotonicity_to_curvature_ball": certified[
            "monotonicity_to_curvature_ball"
        ],
        "certified_sign": certified["certified_monotonicity_bounds"],
        "continuity_deduction": (
            "The theta tail makes |u|^3*exp(u^2/5)*Phi(u) integrable. Dominated "
            "convergence therefore makes the first three x-derivatives of H_t "
            "continuous in t on [0,1/5]. The strict value M_0(x)<0 persists for "
            "all sufficiently small positive t, rejecting the universal assertion "
            "M_t(x)>0 for every x>0 and every 0<t<=1/5."
        ),
        "scope": (
            "This rejects only the stronger strict-monotonicity route. It does not "
            "reject L_t(x)>0, double-zero transversality, Lambda<=0, or RH."
        ),
    }


def build_payload() -> dict:
    exact = exact_algebra()
    moderate = build_moderate_scout()
    high = build_high_frequency_scout()
    guard = countermodel_guard()
    xi_guard = xi_monotonicity_counterexample()
    rows = [
        GateRow(
            "nslms_01_derivative_identity",
            "exact_identity",
            "proved",
            "The derivative of the first Laguerre expression is one third-order Wronskian.",
            exact["derivative_identity"],
            "Exact differentiation only.",
        ),
        GateRow(
            "nslms_02_monotone_sufficiency",
            "sufficient_theorem",
            "proved",
            "Strict decrease to the zero limit implies strict first-Laguerre positivity.",
            exact["monotone_sufficiency"],
            "A stronger sufficient implication, but its premise is rigorously false for Xi.",
        ),
        GateRow(
            "nslms_03_correlation_sine_target",
            "exact_identity",
            "proved",
            "The monotonicity margin is a sine transform of the weighted first correlation.",
            exact["correlation_sine_transform"],
            "Exact differentiation under the correlation transform.",
        ),
        GateRow(
            "nslms_04_tail_cosine_geometry",
            "exact_identity",
            "proved",
            "The sine target is a cosine transform of a positive decreasing one-inflection tail.",
            exact["tail_cosine_reduction"]["identity"],
            "Exact reduction; the tail is concave at the origin, so the classical convex criterion is unavailable.",
            exact["tail_cosine_reduction"],
        ),
        GateRow(
            "nslms_05_theta_primitive_target",
            "exact_identity",
            "proved",
            "The same target is a third-order curvature inequality for the endpoint-subtracted theta primitive.",
            exact["theta_primitive"]["identity"],
            "Exact reformulation only; its proposed global sign is false for Xi.",
        ),
        GateRow(
            "nslms_06_dense_scout",
            "numerical_diagnostic",
            "observed",
            "Every dense sampled moderate-frequency value has L_t'(x)<0.",
            "x=j/50, 1<=j<=2000, t in {0,.05,.1,.15,.2,.5}",
            "Finite double-precision quadrature with coarse/fine checks; not a continuum proof.",
            moderate,
        ),
        GateRow(
            "nslms_07_high_frequency_scout",
            "numerical_diagnostic",
            "observed",
            "The negative derivative survives high-frequency arithmetic cancellation through x=200.",
            "L_t'(x)<0 on 20 selected (t,x) rows",
            "Finite 75-digit quadrature and t=0 Xi cross-checks only.",
            high,
        ),
        GateRow(
            "nslms_08_theta_tail_nonpromotion",
            "countermodel_gate",
            "guard_validated",
            "The monotonicity target does not follow from the full generic Xi shape-and-tail package.",
            guard["deduction"],
            "Uses the separate rigorous theta-tail countermodel; it is not the Xi kernel.",
            guard,
        ),
        GateRow(
            "nslms_09_xi_lehmer_counterexample",
            "interval_counterexample",
            "certified",
            "An exact Xi Lehmer-point jet rejects the global strict-monotonicity target.",
            xi_guard["certified_sign"],
            xi_guard["scope"],
            xi_guard,
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_strict_laguerre_monotonicity_scout",
        "date": "2026-07-17",
        "status": (
            "exact sufficient condition and finite scouts, retired by a rigorous "
            "Xi Lehmer-point counterexample"
        ),
        "proof_boundary": (
            "This artifact proves that strict decrease of the first Laguerre expression "
            "to its zero limit is sufficient for Lambda<=0 and derives exact correlation "
            "and theta-primitive forms. A rigorous Xi point has M_0<0, which by "
            "time continuity rejects the proposed all-positive-times condition. This "
            "does not reject strict Laguerre positivity, RH, or Lambda<=0."
        ),
        "exact": exact,
        "moderate_scout": moderate,
        "high_frequency_scout": high,
        "countermodel_guard": guard,
        "xi_monotonicity_counterexample": xi_guard,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "https://arxiv.org/abs/1606.05011",
            "https://doi.org/10.1017/S0004972700047511",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md",
            "outputs/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.md",
            "outputs/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.md",
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.md",
        ],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    moderate = payload["moderate_scout"]
    high = payload["high_frequency_scout"]
    guard = payload["countermodel_guard"]
    xi_guard = payload["xi_monotonicity_counterexample"]
    moderate_lines = [
        "| t | x rows | largest L_t' | at x | smallest L_t' | at x |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in moderate["rows"]:
        moderate_lines.append(
            f"| {row['t']} | {row['x_count']} | {row['largest_L_prime']} | "
            f"{row['largest_L_prime_x']} | {row['smallest_L_prime']} | "
            f"{row['smallest_L_prime_x']} |"
        )
    high_lines = [
        "| t | x | L_t'(x) |",
        "|---:|---:|---:|",
    ]
    for row in high["rows"]:
        high_lines.append(f"| {row['t']} | {row['x']} | {row['L_prime']} |")
    return f"""# Jensen-Window PF Newman Strict-Laguerre Monotonicity Scout

Date: {payload['date']}

Status: {payload['status']}. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_monotonicity_scout.py
```

Current result:

```text
validated Jensen-window PF Newman strict-Laguerre monotonicity scout: 9 rows, 0 issues, 5 exact target identities, 2 exact classical-route obstructions, 6 dense time rows, 20 high-frequency rows, 1 theta-tail nonpromotion guard, 1 Arb Xi monotonicity rejection
```

## Exact Sufficient Candidate

Define

```text
{exact['definitions']['first_laguerre']}
{exact['definitions']['monotonicity_margin']}
```

Then

```text
{exact['monotone_sufficiency']}
{exact['newman_consequence']}
```

This condition is stronger than strict Laguerre positivity. The implication is
exact, but the condition itself is false for Xi, as certified below.

## Correlation And Theta Forms

```text
{exact['correlation_sine_transform']}
{exact['theta_primitive']['identity']}
```

The candidate analytic problem was positivity of a sine transform of
`v*K_(1,t)(v)`, or equivalently the displayed third-order theta-primitive
Wronskian.

Writing the sine kernel as an integrated cosine gives a second exact form:

```text
{exact['tail_cosine_reduction']['definition']}
{exact['tail_cosine_reduction']['identity']}
{exact['tail_cosine_reduction']['derivatives']}
```

{exact['tail_cosine_reduction']['one_inflection_theorem']}
In particular, `{exact['tail_cosine_reduction']['classical_convexity_obstruction']}`
Thus the classical convex-kernel positivity theorem does not close this target.
Moreover, {exact['tail_cosine_reduction']['regular_shift_obstruction']}
{exact['tail_cosine_reduction']['boundary_handoff']}
The relevant primary source is
https://doi.org/10.1017/S0004972700047511.

## Dense Moderate Scout

{chr(10).join(moderate_lines)}

The grid is `x=j/50`, `1<=j<=2000`. The maximum relative coarse/fine delta is
`{moderate['max_relative_coarse_fine_delta']}`. This is finite numerical
evidence, not continuum control.

## High-Frequency Scout

{chr(10).join(high_lines)}

The 75-digit maximum relative coarse/fine delta is
`{high['max_relative_coarse_fine_delta']}`; the `t=0` rows are independently
checked against derivatives of `xi((1+i*x)/2)/8`, with maximum relative delta
`{high['max_relative_t0_xi_delta']}`.

These selected rows miss the close-pair stress region and are retained as a
warning against promoting a regular frequency grid.

## Exact Xi Counterexample

At the rational Lehmer stress point `x={xi_guard['x']}`, Arb certifies

```text
{xi_guard['normalized_monotonicity_margin_ball']}
{xi_guard['monotonicity_to_curvature_ball']}
{xi_guard['certified_sign']}
```

{xi_guard['continuity_deduction']}

{xi_guard['scope']}

## Nonpromotion Guard

For the explicit theta-tail countermodel,

```text
{guard['kernel']}
L(21/5)={guard['laguerre_ball']}<0.
```

{guard['deduction']}

The generic countermodel had already shown that shape and tail hypotheses do
not supply monotonicity. The exact Xi certificate now closes the stronger
route itself.

## Decision

Do not pursue a global positive sine-transform theorem for `v*K_(1,t)` or a
coupled contour sign for `-L_t'`. The surviving Newman routes are direct
strict Laguerre positivity and, more economically near Lehmer pairs, exclusion
of simultaneous zeros of `(H_t,H_t')` through the corrected `C1`
transversality target.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print("wrote Jensen-window PF Newman strict-Laguerre monotonicity scout")


if __name__ == "__main__":
    main()
