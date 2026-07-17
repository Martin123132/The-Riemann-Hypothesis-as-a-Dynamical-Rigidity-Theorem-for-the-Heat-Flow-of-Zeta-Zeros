#!/usr/bin/env python3
"""Build the Polymath-15 normalized Laguerre bridge for Newman heat flow."""

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
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.md"
)
SOURCE_URL = "https://arxiv.org/abs/1904.12438"
X_VALUES = (200, 1000, 10_000, 1_000_000)
TIME_VALUES = ("0", "0.01", "0.1", "0.5")


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def alpha(s: mp.mpc) -> mp.mpc:
    return 1 / (2 * s) + 1 / (s - 1) + mp.log(s / (2 * mp.pi)) / 2


def alpha_prime(s: mp.mpc) -> mp.mpc:
    return -1 / (2 * s**2) - 1 / (s - 1) ** 2 + 1 / (2 * s)


def alpha_second(s: mp.mpc) -> mp.mpc:
    return 1 / s**3 + 2 / (s - 1) ** 3 - 1 / (2 * s**2)


def phase_geometry(x: mp.mpf, t: mp.mpf) -> dict:
    s = (1 - mp.j * x) / 2
    a = alpha(s)
    a1 = alpha_prime(s)
    a2 = alpha_second(s)
    log_m_first = a + (t / 2) * a * a1
    log_m_second = a1 + (t / 2) * (a1**2 + a * a2)
    beta_first = -mp.re(log_m_first) / 2
    beta_second = -mp.im(log_m_second) / 4
    potential = mp.re(log_m_second) / 4
    lower_without_potential = 4 * beta_first**2 - 2 * abs(beta_second)
    return {
        "beta_first": beta_first,
        "beta_second": beta_second,
        "potential": potential,
        "single_saddle_lower_without_potential": lower_without_potential,
    }


def phase_diagnostics(dps: int = 60) -> list[dict]:
    mp.mp.dps = dps
    rows: list[dict] = []
    for t_text in TIME_VALUES:
        t = mp.mpf(t_text)
        for x_int in X_VALUES:
            values = phase_geometry(mp.mpf(x_int), t)
            rows.append(
                {
                    "t": t_text,
                    "x": x_int,
                    **{
                        key: mp.nstr(value, 35, strip_zeros=False)
                        for key, value in values.items()
                    },
                }
            )
    return rows


def build_exact() -> dict:
    a0, a1, a2, z0, z1, z2, v = sp.symbols(
        "a0 a1 a2 z0 z1 z2 v", real=True
    )
    h0 = a0 * z0
    h1 = a1 * z0 + a0 * z1
    h2 = a2 * z0 + 2 * a1 * z1 + a0 * z2
    normalized = a0**2 * (
        z1**2 - z0 * z2 + ((a1 / a0) ** 2 - a2 / a0) * z0**2
    )
    if sp.simplify(h1**2 - h0 * h2 - normalized) != 0:
        raise RuntimeError("normalized Laguerre identity failed")

    p0, p1, p2, e0, e1, e2 = sp.symbols(
        "p0 p1 p2 e0 e1 e2", real=True
    )
    curvature = lambda q0, q1, q2: q1**2 - q0 * q2 + v * q0**2
    perturbation = sp.expand(
        curvature(p0 + e0, p1 + e1, p2 + e2)
        - curvature(p0, p1, p2)
    )
    expected_perturbation = (
        2 * p1 * e1
        + e1**2
        - p0 * e2
        - e0 * p2
        - e0 * e2
        + 2 * v * p0 * e0
        + v * e0**2
    )
    if sp.simplify(perturbation - expected_perturbation) != 0:
        raise RuntimeError("curvature perturbation identity failed")

    beta, beta1, beta2 = sp.symbols("beta beta1 beta2", real=True)
    single0 = 2 * sp.cos(beta)
    single1 = -2 * beta1 * sp.sin(beta)
    single2 = -2 * beta2 * sp.sin(beta) - 2 * beta1**2 * sp.cos(beta)
    single_curvature = sp.expand_trig(
        sp.trigsimp(curvature(single0, single1, single2))
    )
    expected_single = (
        4 * beta1**2
        + 2 * beta2 * sp.sin(2 * beta)
        + 4 * v * sp.cos(beta) ** 2
    )
    if sp.trigsimp(single_curvature - expected_single) != 0:
        raise RuntimeError("single-saddle curvature identity failed")

    return {
        "published_input": {
            "region": "0<t<=1/2, 0<=y<=1, x>=200",
            "normalizer": "B_t(x+i*y)=M_t((1+y-i*x)/2)",
            "M_t": "M_t(s)=exp((t/4)*alpha(s)^2)*M_0(s)",
            "alpha": "alpha(s)=1/(2s)+1/(s-1)+(1/2)Log(s/(2*pi))",
            "cutoff": "N=floor(sqrt(x/(4*pi)+t/16))",
            "coefficient": "b_n^t=exp((t/4)*log(n)^2)",
            "estimate": "H_t(x+i*y)/B_t(x+i*y)=f_t(x+i*y)+r_t(x+i*y), |r_t|<=e_A+e_B+e_C0",
        },
        "real_axis": {
            "s_minus": "s=(1-i*x)/2",
            "shifted_exponent": "s_*=s+(t/2)*alpha(s)",
            "dirichlet_sum": "D_(N,t)(x)=sum_(n<=N) exp((t/4)*log(n)^2)*n^(-s_*)",
            "collapse": "f_t(x)=D_(N,t)(x)+(conj(B_t(x))/B_t(x))*conj(D_(N,t)(x))",
            "real_main": "P_(N,t)(x)=2*Re(exp(i*beta_t(x))*D_(N,t)(x)), B_t(x)=A_t(x)*exp(i*beta_t(x))",
        },
        "normalized_curvature": (
            "L[H_t]=A_t^2*C_t[Z_t], Z_t=H_t/A_t, "
            "C_t[Z]=Z'^2-Z*Z''+V_t*Z^2, V_t=-(log A_t)''"
        ),
        "phase_geometry": {
            "log_derivatives": (
                "g_t'=alpha+(t/2)*alpha*alpha'; "
                "g_t''=alpha'+(t/2)*(alpha'^2+alpha*alpha'')"
            ),
            "beta_first": "beta_t'=-(1/2)*Re(g_t'(s))",
            "beta_second": "beta_t''=-(1/4)*Im(g_t''(s))",
            "potential": "V_t=(1/4)*Re(g_t''(s))",
        },
        "cauchy_transfer": {
            "hypothesis": (
                "r=H_t/B_t-f_(N,t) is holomorphic on |z-x|<=rho and "
                "sup|r|<=E_rho"
            ),
            "epsilon_0": "eps_0=E_rho",
            "epsilon_1": "eps_1=(1/rho+|beta'|)*E_rho",
            "epsilon_2": (
                "eps_2=(2/rho^2+2*|beta'|/rho+|beta''|+beta'^2)*E_rho"
            ),
            "conclusion": "|(Z_t-P_(N,t))^(j)|<=eps_j for j=0,1,2",
        },
        "curvature_error": (
            "Err=2*|P'|*eps_1+eps_1^2+|P|*eps_2+eps_0*|P''|+"
            "eps_0*eps_2+|V|*(2*|P|*eps_0+eps_0^2); "
            "C[P]>Err implies L[H_t]>0"
        ),
        "single_saddle": (
            "P_0=2*cos(beta); C[P_0]=4*beta'^2+2*beta''*sin(2*beta)+"
            "4*V*cos(beta)^2"
        ),
        "cell": (
            "For fixed N, 4*pi*(N^2-t/16)<=x<"
            "4*pi*((N+1)^2-t/16); a Cauchy disk must remain in one cell."
        ),
        "boundary_layer": (
            "|b_n^t*n^(-s_*)|<=n^(-1/2)*exp(-(t/4)*log(n)*"
            "log((x/(4*pi))/n)+(t/(2*x^2))*log(n)); hence n>=2 "
            "dominance is controlled by t*log(x), and the nonuniform regime is "
            "t*log(x)=O(1)."
        ),
    }


def build_artifact(dps: int) -> dict:
    exact = build_exact()
    diagnostics = phase_diagnostics(dps)
    if any(mp.mpf(row["potential"]) <= 0 for row in diagnostics):
        raise RuntimeError("sampled normalizer potential lost positivity")
    if any(
        mp.mpf(row["single_saddle_lower_without_potential"]) <= 0
        for row in diagnostics
    ):
        raise RuntimeError("sampled single-saddle lower bound lost positivity")
    rows = [
        GateRow(
            id="np15nlb_01_published_input",
            role="published_theorem_input",
            readiness="available_published",
            claim="Polymath 15 supplies an effective positive-time Riemann-Siegel approximation with the saddle-scale cutoff.",
            formula=exact["published_input"]["estimate"],
            proof_boundary="Imported theorem; this artifact does not reprove its analytic estimates.",
        ),
        GateRow(
            id="np15nlb_02_real_axis_collapse",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="On the real axis the two complex saddle sums collapse to one real Hardy-type finite sum.",
            formula=exact["real_axis"]["real_main"],
            proof_boundary="Exact conjugation algebra for a fixed cutoff N.",
        ),
        GateRow(
            id="np15nlb_03_normalized_curvature",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Positive Laguerre curvature is invariant under the nonvanishing amplitude after adding an explicit potential.",
            formula=exact["normalized_curvature"],
            proof_boundary="Exact product differentiation; no sign is assumed.",
        ),
        GateRow(
            id="np15nlb_04_phase_potential",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The normalizer phase derivatives and curvature potential are explicit elementary functions.",
            formula=exact["phase_geometry"]["log_derivatives"],
            proof_boundary="Exact differentiation of the published normalizer.",
        ),
        GateRow(
            id="np15nlb_05_cauchy_jet_transfer",
            role="exact_lemma",
            readiness="conditional_ready",
            claim="A uniform holomorphic remainder bound on one cutoff cell gives explicit zeroth-through-second real jet errors.",
            formula=exact["cauchy_transfer"]["conclusion"],
            proof_boundary="Requires a disk staying inside the theorem region with one fixed N.",
        ),
        GateRow(
            id="np15nlb_06_curvature_error",
            role="exact_lemma",
            readiness="conditional_ready",
            claim="The finite main curvature has an explicit sufficient margin over the transferred remainder.",
            formula=exact["curvature_error"],
            proof_boundary="A sufficient criterion; its uniform main-sum margin remains to be proved.",
        ),
        GateRow(
            id="np15nlb_07_single_saddle",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The n=1 saddle has a phase-speed curvature floor before all arithmetic perturbations.",
            formula=exact["single_saddle"],
            proof_boundary="Exact formula; global positivity of its lower bound is not inferred from samples.",
        ),
        GateRow(
            id="np15nlb_08_phase_scout",
            role="finite_diagnostic",
            readiness="diagnostic_only",
            claim="Sampled normalizer potential and phase-only curvature lower bounds are positive throughout the tested range.",
            formula="V_t>0 and 4*beta_t'^2-2*|beta_t''|>0 on all 16 sampled points",
            proof_boundary="Point diagnostics only, not an interval or ray proof.",
            diagnostics=diagnostics,
        ),
        GateRow(
            id="np15nlb_09_scale_match",
            role="exact_handoff",
            readiness="ready_to_apply",
            claim="The discovered modular-blend saddle scale matches the cutoff in the published theorem.",
            formula=exact["published_input"]["cutoff"],
            proof_boundary="Leading-scale match; the two decompositions have different remainders.",
        ),
        GateRow(
            id="np15nlb_10_transition_collar",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Cauchy differentiation leaves discrete collars at cutoff transitions.",
            formula=exact["cell"],
            proof_boundary="Needs overlapping adjacent-N estimates or a direct transition-collar theorem.",
        ),
        GateRow(
            id="np15nlb_11_boundary_layer",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="Single-saddle dominance is nonuniform only in the shrinking boundary layer t*log(x)=O(1).",
            formula=exact["boundary_layer"],
            proof_boundary="Coefficient reduction only; no uniform tail constant or Laguerre sign is proved.",
        ),
        GateRow(
            id="np15nlb_12_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The live tail programme splits into a dominant-saddle ray and a near-zero-time arithmetic boundary layer.",
            formula=(
                "Prove C[P_0+q]>Err for t*log(x)>=C, certify cutoff collars, "
                "then retain the full Dirichlet coupling for 0<t*log(x)<C."
            ),
            proof_boundary="Open; not strict Laguerre positivity, Lambda<=0, or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_normalized_laguerre_bridge",
        "date": "2026-07-16",
        "status": (
            "exact normalized-curvature transfer from the published Polymath-15 "
            "effective Riemann-Siegel theorem; not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact proves the real-axis collapse, amplitude-normalized "
            "Laguerre identity, Cauchy jet-transfer lemma, perturbation criterion, "
            "and the t*log(x) boundary-layer reduction. It does not prove a "
            "uniform finite-main curvature margin, cutoff-transition collars, "
            "strict Laguerre positivity, Lambda<=0, or RH."
        ),
        "source": {
            "title": "Effective approximation of heat flow evolution of the Riemann xi function, and a new upper bound for the de Bruijn-Newman constant",
            "authors": "D. H. J. Polymath",
            "url": SOURCE_URL,
            "theorem": "Theorem 1.3 (effective Riemann-Siegel approximation)",
        },
        "parameters": {
            "diagnostic_dps": dps,
            "diagnostic_x_values": list(X_VALUES),
            "diagnostic_t_values": list(TIME_VALUES),
        },
        "exact": exact,
        "diagnostics": diagnostics,
        "rows": [asdict(row) for row in rows],
        "sources": [
            SOURCE_URL,
            "outputs/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    table = [
        "| t | x | beta' | beta'' | V | phase-only lower |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in artifact["diagnostics"]:
        table.append(
            "| {t} | {x} | {b1} | {b2} | {v} | {lower} |".format(
                t=row["t"],
                x=row["x"],
                b1=mp.nstr(mp.mpf(row["beta_first"]), 9),
                b2=mp.nstr(mp.mpf(row["beta_second"]), 8),
                v=mp.nstr(mp.mpf(row["potential"]), 8),
                lower=mp.nstr(
                    mp.mpf(row["single_saddle_lower_without_potential"]), 9
                ),
            )
        )
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Normalized Laguerre Bridge",
            "",
            "Date: 2026-07-16",
            "",
            "Status: exact normalized-curvature transfer from a published",
            "positive-time Riemann-Siegel theorem. This is not a proof of",
            "`Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.py",
            "```",
            "",
            "Published input: D. H. J. Polymath, [Effective approximation of heat",
            f"flow evolution of the Riemann xi function]({SOURCE_URL}), Theorem 1.3.",
            "",
            "## Real Normalization",
            "",
            "For `s=(1-i*x)/2`, put",
            "",
            "```text",
            exact["published_input"]["normalizer"],
            exact["published_input"]["M_t"],
            exact["published_input"]["alpha"],
            exact["published_input"]["cutoff"],
            exact["real_axis"]["dirichlet_sum"],
            "```",
            "",
            "On the real axis the published two-sum approximation has `kappa=0`",
            "and its second sum is the phase-adjusted conjugate of the first:",
            "",
            "```text",
            exact["real_axis"]["collapse"],
            exact["real_axis"]["real_main"],
            "```",
            "",
            "Writing `Z_t=H_t/A_t`, exact product differentiation gives",
            "",
            "```text",
            exact["normalized_curvature"],
            "```",
            "",
            "Thus the complex approximation becomes a real finite curvature",
            "problem with an explicit normalizer potential.",
            "",
            "## Jet Transfer",
            "",
            "If one fixed-`N` approximation is holomorphic on a disk of radius",
            "`rho` and its remainder is bounded there by `E_rho`, Cauchy's",
            "estimate and the phase product give",
            "",
            "```text",
            exact["cauchy_transfer"]["epsilon_0"],
            exact["cauchy_transfer"]["epsilon_1"],
            exact["cauchy_transfer"]["epsilon_2"],
            exact["curvature_error"],
            "```",
            "",
            "This is the missing derivative bridge: a scalar Riemann-Siegel",
            "remainder can certify the Laguerre expression once the finite main",
            "sum clears the explicit error.",
            "",
            "## Phase Floor",
            "",
            "The `n=1` term is exactly `P_0=2*cos(beta_t)`, with",
            "",
            "```text",
            exact["single_saddle"],
            "```",
            "",
            *table,
            "",
            "These 16 rows are diagnostics, not a ray proof. They show that the",
            "normalizer itself supplies a large positive phase-speed floor; the",
            "remaining issue is controlling the arithmetic terms and the theorem",
            "remainder in `C^2`.",
            "",
            "## Remaining Geometry",
            "",
            "The published cutoff exactly matches the square-root scale found by",
            "the modular-blend saddle audit. Cauchy disks cannot cross its discrete",
            "transition points without an overlapping adjacent-`N` estimate:",
            "",
            "```text",
            exact["cell"],
            "```",
            "",
            "For the nonconstant arithmetic terms, the published lower bound on",
            "`Re(s_*)` yields",
            "",
            "```text",
            exact["boundary_layer"],
            "```",
            "",
            "So the global tail now has two sharply separated tasks: prove",
            "single-saddle dominance when `t*log(x)` is large, and handle the",
            "shrinking but RH-critical boundary layer `t*log(x)=O(1)` without",
            "discarding the coupled Dirichlet cancellation.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--dps", type=int, default=60)
    args = parser.parse_args()
    artifact = build_artifact(args.dps)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 normalized Laguerre bridge: "
        f"{len(artifact['rows'])} rows, {len(artifact['diagnostics'])} diagnostics"
    )


if __name__ == "__main__":
    main()
