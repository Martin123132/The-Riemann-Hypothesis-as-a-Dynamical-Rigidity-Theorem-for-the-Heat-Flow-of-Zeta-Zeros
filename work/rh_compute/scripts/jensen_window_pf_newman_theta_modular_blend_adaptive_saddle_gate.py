#!/usr/bin/env python3
"""Build the frequency-adaptive saddle gate for the modular theta blend."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.md"
)
PANELS = (
    ("0", ".04"),
    (".04", ".08"),
    (".08", ".14"),
    (".14", ".22"),
    (".22", ".34"),
    (".34", ".5"),
    (".5", ".72"),
    (".72", "1"),
    ("1", "1.4"),
    ("1.4", "2"),
    ("2", "2.8"),
)
X_VALUES = (80, 100, 120, 150, 200)
TIMES = ("0", "0.5")


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def theta_summand(u: mp.mpf, n: int) -> mp.mpf:
    return (
        2 * mp.pi**2 * n**4 * mp.exp(9 * u)
        - 3 * mp.pi * n**2 * mp.exp(5 * u)
    ) * mp.exp(-mp.pi * n**2 * mp.exp(4 * u))


def blend_weight(u: mp.mpf) -> mp.mpf:
    return (1 + mp.erf(3 * mp.sinh(4 * u))) / 2


def quadrature_rows(node_count: int, theta_terms: int) -> list[tuple]:
    nodes, weights = mp.gauss_quadrature(node_count, "legendre")
    rows: list[tuple] = []
    for left_text, right_text in PANELS:
        left = mp.mpf(left_text)
        right = mp.mpf(right_text)
        for node, weight in zip(nodes, weights):
            u = (right - left) * node / 2 + (right + left) / 2
            qweight = (right - left) * weight / 2
            omega = blend_weight(u)
            forward = [theta_summand(u, n) for n in range(1, theta_terms + 1)]
            blocks = [
                omega * forward[n - 1]
                + (1 - omega) * theta_summand(-u, n)
                for n in range(1, theta_terms + 1)
            ]
            rows.append((u, qweight, tuple(blocks), mp.fsum(forward)))
    return rows


def jet_from_values(
    rows: list[tuple], values: list[mp.mpf], t: mp.mpf, x: mp.mpf
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    weighted = [
        (row[0], row[1] * mp.exp(t * row[0] ** 2) * value)
        for row, value in zip(rows, values, strict=True)
    ]
    return (
        mp.fsum(weight * mp.cos(x * u) for u, weight in weighted),
        mp.fsum(-weight * u * mp.sin(x * u) for u, weight in weighted),
        mp.fsum(-weight * u * u * mp.cos(x * u) for u, weight in weighted),
        mp.fsum(weight * u**3 * mp.sin(x * u) for u, weight in weighted),
    )


def laguerre(jet: tuple[mp.mpf, ...]) -> mp.mpf:
    return jet[1] ** 2 - jet[0] * jet[2]


def monotonicity_margin(jet: tuple[mp.mpf, ...]) -> mp.mpf:
    """Return -L'(x)=H*H'''-H'*H''."""
    return jet[0] * jet[3] - jet[1] * jet[2]


def transition_index(a: int, t: mp.mpf, x: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    y = mp.atan(x / a) / 4
    for _ in range(30):
        y = mp.atan((x + 2 * t * y) / a) / 4
    n_star = mp.sqrt(mp.sqrt(a * a + (x + 2 * t * y) ** 2) / (4 * mp.pi))
    return y, n_star


def compute_adaptive_rows(
    node_count: int, dps: int, theta_terms: int
) -> list[dict]:
    mp.mp.dps = dps
    rows = quadrature_rows(node_count, theta_terms)
    block_columns = [
        [row[2][n] for row in rows] for n in range(theta_terms)
    ]
    full_column = [row[3] for row in rows]
    output: list[dict] = []
    for t_text in TIMES:
        t = mp.mpf(t_text)
        for x_int in X_VALUES:
            x = mp.mpf(x_int)
            full_jet = jet_from_values(rows, full_column, t, x)
            full_laguerre = laguerre(full_jet)
            full_margin = monotonicity_margin(full_jet)
            cumulative = [mp.mpf(0) for _ in rows]
            partial_rows: list[dict] = []
            first_half_relative: int | None = None
            first_margin_half_relative: int | None = None
            for n, block_column in enumerate(block_columns[:10], start=1):
                cumulative = [
                    left + right
                    for left, right in zip(cumulative, block_column, strict=True)
                ]
                partial_jet = jet_from_values(rows, cumulative, t, x)
                partial_laguerre = laguerre(partial_jet)
                partial_margin = monotonicity_margin(partial_jet)
                ratio = partial_laguerre / full_laguerre
                relative_error = abs(ratio - 1)
                margin_ratio = partial_margin / full_margin
                margin_relative_error = abs(margin_ratio - 1)
                if first_half_relative is None and relative_error < mp.mpf("0.5"):
                    first_half_relative = n
                if (
                    first_margin_half_relative is None
                    and margin_relative_error < mp.mpf("0.5")
                ):
                    first_margin_half_relative = n
                partial_rows.append(
                    {
                        "n": n,
                        "partial_to_full_ratio": mp.nstr(
                            ratio, dps - 8, strip_zeros=False
                        ),
                        "relative_laguerre_error": mp.nstr(
                            relative_error, dps - 8, strip_zeros=False
                        ),
                        "partial_margin_to_full_ratio": mp.nstr(
                            margin_ratio, dps - 8, strip_zeros=False
                        ),
                        "relative_margin_error": mp.nstr(
                            margin_relative_error, dps - 8, strip_zeros=False
                        ),
                    }
                )
            y5, n5 = transition_index(5, t, x)
            y9, n9 = transition_index(9, t, x)
            recommended = int(mp.ceil(max(n5, n9))) + 2
            output.append(
                {
                    "t": t_text,
                    "x": x_int,
                    "full_laguerre": mp.nstr(
                        full_laguerre, dps - 8, strip_zeros=False
                    ),
                    "full_monotonicity_margin": mp.nstr(
                        full_margin, dps - 8, strip_zeros=False
                    ),
                    "first_n_relative_error_below_half": first_half_relative,
                    "first_n_margin_relative_error_below_half": first_margin_half_relative,
                    "principal_crossing_y_a5": mp.nstr(y5, 25),
                    "principal_crossing_y_a9": mp.nstr(y9, 25),
                    "transition_index_a5": mp.nstr(n5, 25),
                    "transition_index_a9": mp.nstr(n9, 25),
                    "ceil_max_transition_plus_two": recommended,
                    "adaptive_bound_covers_observed_n": (
                        first_half_relative is not None
                        and first_half_relative <= recommended
                    ),
                    "adaptive_bound_covers_observed_margin_n": (
                        first_margin_half_relative is not None
                        and first_margin_half_relative <= recommended
                    ),
                    "partials": partial_rows,
                }
            )
    return output


def compare_rows(coarse: list[dict], fine: list[dict]) -> dict:
    max_transition_delta = mp.mpf(0)
    stable_first_n = True
    stable_margin_first_n = True
    max_near_ratio_delta = mp.mpf(0)
    max_near_margin_ratio_delta = mp.mpf(0)
    for left, right in zip(coarse, fine, strict=True):
        if (left["t"], left["x"]) != (right["t"], right["x"]):
            raise RuntimeError("adaptive coarse/fine key mismatch")
        stable_first_n &= (
            left["first_n_relative_error_below_half"]
            == right["first_n_relative_error_below_half"]
        )
        stable_margin_first_n &= (
            left["first_n_margin_relative_error_below_half"]
            == right["first_n_margin_relative_error_below_half"]
        )
        n = right["first_n_relative_error_below_half"]
        if n is None:
            continue
        left_ratio = mp.mpf(left["partials"][n - 1]["partial_to_full_ratio"])
        right_ratio = mp.mpf(right["partials"][n - 1]["partial_to_full_ratio"])
        max_near_ratio_delta = max(
            max_near_ratio_delta,
            abs(left_ratio - right_ratio) / max(abs(right_ratio), mp.mpf("1e-100")),
        )
        margin_n = right["first_n_margin_relative_error_below_half"]
        if margin_n is not None:
            left_margin_ratio = mp.mpf(
                left["partials"][margin_n - 1]["partial_margin_to_full_ratio"]
            )
            right_margin_ratio = mp.mpf(
                right["partials"][margin_n - 1]["partial_margin_to_full_ratio"]
            )
            max_near_margin_ratio_delta = max(
                max_near_margin_ratio_delta,
                abs(left_margin_ratio - right_margin_ratio)
                / max(abs(right_margin_ratio), mp.mpf("1e-100")),
            )
        for field in ("transition_index_a5", "transition_index_a9"):
            max_transition_delta = max(
                max_transition_delta,
                abs(mp.mpf(left[field]) - mp.mpf(right[field])),
            )
    return {
        "stable_first_n": bool(stable_first_n),
        "stable_margin_first_n": bool(stable_margin_first_n),
        "max_relative_near_transition_ratio_delta": mp.nstr(
            max_near_ratio_delta, 20
        ),
        "max_relative_near_transition_margin_ratio_delta": mp.nstr(
            max_near_margin_ratio_delta, 20
        ),
        "max_absolute_transition_index_delta": mp.nstr(max_transition_delta, 20),
    }


def build_exact() -> dict:
    u, t, a, x, n = sp.symbols("u t a x n", positive=True, real=True)
    phase = t * u**2 + (a + sp.I * x) * u - sp.pi * n**2 * sp.exp(4 * u)
    derivative = sp.diff(phase, u)
    expected = 2 * t * u + a + sp.I * x - 4 * sp.pi * n**2 * sp.exp(4 * u)
    if sp.simplify(derivative - expected) != 0:
        raise RuntimeError("theta saddle derivative identity failed")
    saddle_zero = sp.log((a + sp.I * x) / (4 * sp.pi * n**2)) / 4
    if sp.simplify(expected.subs({t: 0, u: saddle_zero})) != 0:
        raise RuntimeError("zero-time saddle identity failed")
    y = sp.symbols("y", real=True)
    crossing = sp.expand_complex(expected.subs(u, sp.I * y))
    real_part = sp.simplify(sp.re(crossing))
    imag_part = sp.simplify(sp.im(crossing))
    if real_part != a - 4 * sp.pi * n**2 * sp.cos(4 * y):
        raise RuntimeError("positive-time crossing real part failed")
    if imag_part != 2 * t * y + x - 4 * sp.pi * n**2 * sp.sin(4 * y):
        raise RuntimeError("positive-time crossing imaginary part failed")
    return {
        "component_phase": (
            "F_(a,n,t,x)(u)=t*u^2+(a+i*x)*u-pi*n^2*exp(4u), a in {5,9}"
        ),
        "saddle_equation": (
            "2*t*u+a+i*x-4*pi*n^2*exp(4u)=0"
        ),
        "zero_time_saddle": (
            "u_*=(1/4)*Log((a+i*x)/(4*pi*n^2))"
        ),
        "zero_time_real_part": (
            "Re(u_*)=(1/4)*log(sqrt(a^2+x^2)/(4*pi*n^2))"
        ),
        "zero_time_transition": (
            "n_*(a,0,x)=(a^2+x^2)^(1/4)/(2*sqrt(pi))="
            "sqrt(x/(4*pi))*(1+O_a(x^-2))"
        ),
        "positive_time_crossing": (
            "At u=i*y: 4*pi*n^2*cos(4y)=a and "
            "4*pi*n^2*sin(4y)=x+2*t*y."
        ),
        "positive_time_fixed_point": (
            "4y=atan((x+2*t*y)/a), 0<y<pi/8; "
            "n_*^2=sqrt(a^2+(x+2*t*y)^2)/(4*pi)"
        ),
        "uniform_asymptotic": (
            "For 0<=t<=1/2 and a in {5,9}, n_*(a,t,x)="
            "sqrt(x/(4*pi))*(1+O(x^-1)) as x->infinity."
        ),
        "interpretation": (
            "The modular switch changes side near Re(u_*)=0. A fixed arithmetic "
            "truncation eventually misses the saddle transition; a viable "
            "truncation must grow on the Riemann-Siegel scale sqrt(x/(4*pi))."
        ),
    }


def build_artifact(
    dps: int, coarse_nodes: int, fine_nodes: int, theta_terms: int
) -> dict:
    exact = build_exact()
    coarse = compute_adaptive_rows(coarse_nodes, dps, theta_terms)
    fine = compute_adaptive_rows(fine_nodes, dps, theta_terms)
    convergence = compare_rows(coarse, fine)
    if not convergence["stable_first_n"]:
        raise RuntimeError("adaptive first-N diagnostic is not quadrature-stable")
    if not convergence["stable_margin_first_n"]:
        raise RuntimeError("adaptive margin first-N diagnostic is not quadrature-stable")
    if mp.mpf(convergence["max_relative_near_transition_ratio_delta"]) >= mp.mpf("1e-20"):
        raise RuntimeError("adaptive near-transition ratio convergence failed")
    if mp.mpf(
        convergence["max_relative_near_transition_margin_ratio_delta"]
    ) >= mp.mpf("1e-20"):
        raise RuntimeError("adaptive near-transition margin ratio convergence failed")
    if any(not row["adaptive_bound_covers_observed_n"] for row in fine):
        raise RuntimeError("saddle-index collar failed to cover an observed row")
    if any(not row["adaptive_bound_covers_observed_margin_n"] for row in fine):
        raise RuntimeError("saddle-index collar failed to cover a monotonicity row")
    rows = [
        GateRow(
            id="ntmbasg_01_component_phase",
            role="exact_identity",
            readiness="available_exact",
            claim="Each theta component has an explicit deformed complex phase.",
            formula=exact["component_phase"],
            proof_boundary="Exact component algebra only.",
        ),
        GateRow(
            id="ntmbasg_02_zero_time_saddle",
            role="exact_identity",
            readiness="available_exact",
            claim="At zero Newman time the component saddle is explicit.",
            formula=exact["zero_time_saddle"],
            proof_boundary="Principal saddle geometry, not a contour theorem.",
        ),
        GateRow(
            id="ntmbasg_03_transition_scale",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The modular side transition occurs on the Riemann-Siegel square-root scale.",
            formula=exact["zero_time_transition"],
            proof_boundary="Exact saddle threshold; no Laguerre sign follows.",
        ),
        GateRow(
            id="ntmbasg_04_positive_time_crossing",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The positive-time transition obeys an exact scalar fixed-point law.",
            formula=exact["positive_time_fixed_point"],
            proof_boundary="Saddle crossing only; no steepest-descent remainder is bounded.",
        ),
        GateRow(
            id="ntmbasg_05_uniform_scale",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Bounded positive Newman time preserves the leading square-root transition scale.",
            formula=exact["uniform_asymptotic"],
            proof_boundary="Leading scale only; constants and a collar theorem remain open.",
        ),
        GateRow(
            id="ntmbasg_06_adaptive_scout",
            role="finite_diagnostic",
            readiness="diagnostic_only",
            claim="High-precision rows show that the block counts needed for scale-level Laguerre and monotonicity-margin agreement grow together with the saddle index.",
            formula="N_L=N_M=3,4,4,5,6 at x=80,100,120,150,200 for both sampled times",
            proof_boundary=(
                "Ten point diagnostics only; the exact Lehmer-point counterexample "
                "shows that the sampled monotonicity signs are not globally representative."
            ),
            diagnostics=fine,
        ),
        GateRow(
            id="ntmbasg_07_collar_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="A two-index collar covers all observed transition rows but is not promoted to either sign theorem.",
            formula="max(N_L,N_M)<=ceil(max(n_*5,n_*9))+2 on all ten rows",
            proof_boundary="Finite empirical collar; no global tail sign is asserted.",
            diagnostics=convergence,
        ),
        GateRow(
            id="ntmbasg_08_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The surviving architecture applies the frequency-adaptive modular remainder only to direct Laguerre positivity.",
            formula=(
                "Choose N(x,t) on the saddle scale, retain both a=5 and a=9 "
                "transitions, and prove the Laguerre remainder sign; the global "
                "monotonicity-margin sign is rejected by the Xi Lehmer counterexample."
            ),
            proof_boundary=(
                "Direct Laguerre control remains open; the rejected stronger sign is "
                "not needed for the corrected C1 transversality route."
            ),
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate",
        "date": "2026-07-16",
        "status": (
            "exact modular-blend saddle geometry with finite frequency-adaptive "
            "diagnostics and a retired monotonicity branch; not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact proves the exact component saddle and transition-index "
            "formulas and records ten high-precision adaptive-truncation diagnostics "
            "for both the Laguerre value and its sampled monotonicity margin. The "
            "latter global sign is false by the separate Xi Lehmer certificate. It "
            "does not prove a contour deformation, a uniform two-index collar, the "
            "Laguerre remainder sign, strict Laguerre positivity, Lambda<=0, or RH."
        ),
        "parameters": {
            "dps": dps,
            "coarse_nodes_per_panel": coarse_nodes,
            "fine_nodes_per_panel": fine_nodes,
            "theta_terms": theta_terms,
            "x_values": list(X_VALUES),
            "times": list(TIMES),
            "panels": [list(panel) for panel in PANELS],
        },
        "exact": exact,
        "convergence": convergence,
        "diagnostics": fine,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_theta_modular_blend_gate.md",
            "outputs/jensen_window_pf_newman_theta_modular_blend_high_frequency_scout.md",
            "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    table = [
        "| t | x | n* (a=5) | n* (a=9) | first N for L | first N for -L' | ceil(max n*)+2 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in artifact["diagnostics"]:
        table.append(
            f"| {row['t']} | {row['x']} | "
            f"{mp.nstr(mp.mpf(row['transition_index_a5']), 7)} | "
            f"{mp.nstr(mp.mpf(row['transition_index_a9']), 7)} | "
            f"{row['first_n_relative_error_below_half']} | "
            f"{row['first_n_margin_relative_error_below_half']} | "
            f"{row['ceil_max_transition_plus_two']} |"
        )
    return "\n".join(
        [
            "# Jensen-Window PF Newman Modular-Blend Adaptive-Saddle Gate",
            "",
            "Date: 2026-07-16",
            "",
            "Status: exact saddle geometry with finite frequency-adaptive diagnostics.",
            "This is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.py",
            "```",
            "",
            "## Saddle Law",
            "",
            "The two exponential components of each theta summand have `a=9` and",
            "`a=5`. After the Newman weight and spectral phase,",
            "",
            "```text",
            exact["component_phase"],
            exact["saddle_equation"],
            "```",
            "",
            "At `t=0`,",
            "",
            "```text",
            exact["zero_time_saddle"],
            exact["zero_time_real_part"],
            exact["zero_time_transition"],
            "```",
            "",
            "Thus the saddle crosses the modular blend's `u=0` transition when",
            "`n` is on the Riemann-Siegel scale `sqrt(x/(4*pi))`.",
            "",
            "At positive Newman time the crossing saddle is `u=i*y` and obeys",
            "",
            "```text",
            exact["positive_time_crossing"],
            exact["positive_time_fixed_point"],
            exact["uniform_asymptotic"],
            "```",
            "",
            "## Adaptive Scout",
            "",
            "For each sampled point, the table compares the exact transition index",
            "with the first blended block counts whose Laguerre expression and",
            "strict monotonicity margin `-L_t'` are within 50% of their full",
            "high-precision values:",
            "",
            *table,
            "",
            "Both observed counts agree row by row and grow from three to six, while",
            "a two-index collar beyond the larger `a=5,9` transition covers all ten",
            "rows. This is a stable truncation diagnostic only. The exact Xi Lehmer",
            "counterexample shows that the sampled `-L_t'` signs are not globally",
            "representative.",
            "",
            "## Decision",
            "",
            "The fixed finite truncation is closed, but the blend identifies its",
            "correct replacement: choose `N=N(x,t)` on the transition scale and",
            "prove a sign-aware remainder after retaining both theta components.",
            "That task now concerns `L_t` only. Arb certifies `-L_0'<0` at the",
            "Lehmer stress point, and time continuity rejects the proposed global",
            "monotonicity condition. Near close pairs, corrected `C1` double-zero",
            "transversality remains the lower-derivative Newman target.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--dps", type=int, default=75)
    parser.add_argument("--coarse-nodes", type=int, default=100)
    parser.add_argument("--fine-nodes", type=int, default=140)
    parser.add_argument("--theta-terms", type=int, default=13)
    args = parser.parse_args()
    artifact = build_artifact(
        args.dps, args.coarse_nodes, args.fine_nodes, args.theta_terms
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman modular-blend adaptive-saddle gate: "
        f"{len(artifact['rows'])} rows, {len(artifact['diagnostics'])} diagnostics"
    )


if __name__ == "__main__":
    main()
