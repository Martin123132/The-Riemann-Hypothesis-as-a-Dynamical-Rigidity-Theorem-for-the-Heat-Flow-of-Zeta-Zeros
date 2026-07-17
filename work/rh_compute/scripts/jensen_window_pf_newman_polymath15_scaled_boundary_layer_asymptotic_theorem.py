#!/usr/bin/env python3
"""Build the c=t*log(x/4pi) asymptotic Newman boundary-layer theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.md"
)
SOURCE_URL = "https://arxiv.org/abs/1904.12438"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    amplitude, u, u1, v, v1, theta, potential = sp.symbols(
        "amplitude u u1 v v1 theta potential", real=True
    )
    value = 2 * amplitude * sp.cos(theta)
    first = 2 * amplitude * (
        u * sp.cos(theta) - v * sp.sin(theta)
    )
    second = 2 * amplitude * (
        (u**2 + u1 - v**2) * sp.cos(theta)
        - (2 * u * v + v1) * sp.sin(theta)
    )
    curvature = sp.trigsimp(first**2 - value * second + potential * value**2)
    expected = 4 * amplitude**2 * (
        v**2
        + (potential - u1) * sp.cos(theta) ** 2
        + v1 * sp.sin(2 * theta) / 2
    )
    if sp.trigsimp(sp.expand_trig(curvature - expected)) != 0:
        raise RuntimeError("phase-amplitude curvature identity failed")

    return {
        "scale": (
            "L=log(x/(4*pi)), c=t*L, 4+epsilon<=c<=25, "
            "t=c/L and L->infinity"
        ),
        "main_sum": (
            "D_(N,t)=sum_(n<=N)w_n*n^(-s_*), "
            "w_n=exp((t/4)log(n)^2), N=floor(sqrt(exp(L)+t/16))"
        ),
        "sigma": (
            "Re(s_*)=1/2+c/4+o_epsilon(1)>="
            "3/2+epsilon/8 for sufficiently large L"
        ),
        "coefficient_envelope": (
            "w_n*n^(-Re(s_*))<=n^(-p_epsilon), "
            "p_epsilon=1+epsilon/16, 1<=n<=N"
        ),
        "moments": (
            "D_k=sum_(n<=N)w_n*log(n)^k*n^(-s_*), "
            "Z_k=sum_(n>=1)log(n)^k*n^(-s_*); "
            "|D_k-Z_k|<=25/(4L)*sum_(n>=1)log(n)^(k+2)n^(-p_epsilon)+"
            "sum_(n>N)log(n)^k*n^(-Re(s_*))=O_epsilon(1/L), k=0,1,2"
        ),
        "zeta_floor": (
            "|zeta(s_*)|>=1/zeta(3/2+epsilon/8)>0"
        ),
        "jet_convergence": (
            "D_(N,t)^(j)=(zeta(s_*(x)))^(j)+O_epsilon(1/L), j=0,1,2"
        ),
        "phase_amplitude": (
            "For D=a*exp(i*psi), theta=beta+psi, u=(log a)', v=theta': "
            "C[P]=4a^2*(v^2+(V-u')cos(theta)^2+(v'/2)sin(2theta))"
        ),
        "phase_bounds": (
            "a>=m_epsilon>0, u'=O_epsilon(1), "
            "v=-L/4+O_epsilon(1), v'=O_epsilon(1), V=o(1)"
        ),
        "main_curvature": (
            "C_t[P_(N,t)]>=c_epsilon*L^2 for all sufficiently large L"
        ),
        "raw_remainder": (
            "On fixed radius-1/4 collars, the Polymath-15 raw remainder plus "
            "one adjacent-cutoff block is O_epsilon(exp(-kappa_epsilon*L)); "
            "its normalized C2 curvature cost is o_epsilon(L^2)"
        ),
        "theorem": (
            "For every epsilon>0 there exists L_epsilon such that "
            "0<t<=1/2, L>=L_epsilon, and t*L>=4+epsilon imply L_t(x)>0"
        ),
        "critical_layer": "0<t*log(x/(4*pi))<=4+o(1)",
        "current_global_boundary": (
            "The stronger oscillatory-zeta handoff lowers the corpus-wide "
            "asymptotic boundary to c_*=4911678521/1933561194; the threshold "
            "four remains only the boundary of this absolute-value method"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15sblat_01_scaled_regime",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The nonuniform positive-time tail is governed by the single scaled parameter c=tL.",
            formula=exact["scale"],
            proof_boundary="The branch c>=25 is already covered by the explicit global-ray certificate.",
        ),
        GateRow(
            id="np15sblat_02_power_envelope",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="Above c=4+epsilon all weighted cutoff coefficients share a summable power envelope.",
            formula=exact["coefficient_envelope"],
            proof_boundary="Uniform for 4+epsilon<=c<=25 after increasing L_epsilon.",
        ),
        GateRow(
            id="np15sblat_03_zeta_moments",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="The finite weighted Dirichlet moments converge through second order to ordinary zeta moments.",
            formula=exact["moments"],
            proof_boundary="Uses e^q-1<=q*e^q, t<=25/L, and the summable envelope.",
        ),
        GateRow(
            id="np15sblat_04_euler_floor",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="The limiting zeta factor is uniformly separated from zero by its Euler product.",
            formula=exact["zeta_floor"],
            proof_boundary="The line Re(s_*)>=3/2+epsilon/8 lies in absolute convergence.",
        ),
        GateRow(
            id="np15sblat_05_c2_convergence",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="The finite main sum and its first two x derivatives inherit uniform zeta control.",
            formula=exact["jet_convergence"],
            proof_boundary="Uses s_*'=-i/2+o(1) and s_*''=o(1).",
        ),
        GateRow(
            id="np15sblat_06_phase_amplitude_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="A nonvanishing complex Dirichlet amplitude converts the real Laguerre expression into phase speed plus bounded logarithmic corrections.",
            formula=exact["phase_amplitude"],
            proof_boundary="Exact symbolic identity.",
        ),
        GateRow(
            id="np15sblat_07_main_asymptotic_curvature",
            role="proved_asymptotic_theorem",
            readiness="ready_to_apply",
            claim="The finite Riemann-Siegel main curvature is uniformly positive above the absolute-convergence boundary.",
            formula=exact["main_curvature"],
            proof_boundary="The threshold L_epsilon is existential, not numerically optimized.",
        ),
        GateRow(
            id="np15sblat_08_remainder_transfer",
            role="proved_asymptotic_theorem",
            readiness="ready_to_apply",
            claim="The published analytic remainder and cutoff jumps are exponentially smaller than the main curvature.",
            formula=exact["raw_remainder"],
            proof_boundary="Uses the same raw Schwarz-Cauchy architecture as the explicit c>=25 theorem.",
        ),
        GateRow(
            id="np15sblat_09_exact_h_asymptotic",
            role="proved_theorem",
            readiness="ready_to_apply",
            claim="Exact Newman first-Laguerre positivity holds asymptotically for every scaled time strictly above four.",
            formula=exact["theorem"],
            proof_boundary="Quantified separately for each epsilon>0.",
        ),
        GateRow(
            id="np15sblat_10_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="This absolute-value Dirichlet method alone leaves the critical scaled strip at and below four.",
            formula=exact["critical_layer"],
            proof_boundary=(
                "Historical method handoff; the stronger oscillatory-zeta "
                "theorem supersedes it as the corpus-wide boundary. Not "
                "Lambda<=0 or RH."
            ),
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem",
        "date": "2026-07-17",
        "status": (
            "historical quantified asymptotic strict Laguerre theorem for exact "
            "H_t whenever t*log(x/(4*pi))>=4+epsilon; the theorem remains "
            "valid but its global boundary is superseded; not a proof of "
            "Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "For every fixed epsilon>0 this artifact proves the existence of an "
            "L_epsilon beyond which L_t(x)>0 whenever 0<t<=1/2 and "
            "t*log(x/(4*pi))>=4+epsilon. It does not supply a practical numerical "
            "L_epsilon. Its t*log(x/(4*pi))<=4+o(1) handoff is the historical "
            "boundary of the absolute-value method, not the current global "
            "boundary: the oscillatory-zeta theorem lowers that boundary to "
            "c_*=4911678521/1933561194. This artifact does not prove positivity "
            "for all x and t, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.md",
            "outputs/jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.md",
            "outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Scaled Boundary-Layer Asymptotic Theorem",
            "",
            "Date: 2026-07-17",
            "",
            "Status: historical quantified asymptotic strict first-Laguerre",
            "positivity for the exact Newman heat flow above the scaled threshold",
            "four. The theorem remains valid, but its corpus-wide boundary has",
            "been superseded. This is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.py",
            "```",
            "",
            "The effective finite approximation and error bounds come from D. H. J.",
            f"Polymath, [Effective approximation of heat flow evolution of the Riemann xi function]({SOURCE_URL}),",
            "Theorem 1.3. The c>=25 branch is supplied by the explicit global-ray",
            "certificate; this note's absolute-value argument treats the historical",
            "bounded-c branch 4+epsilon<=c<=25.",
            "",
            "## Scaled Regime",
            "",
            "Fix `epsilon>0` and put",
            "",
            "```text",
            exact["scale"],
            exact["main_sum"],
            "```",
            "",
            "The exact critical-axis formula for the published exponent gives",
            "",
            "```text",
            exact["sigma"],
            "```",
            "",
            "At the largest retained index, `log(N)=L/2+o(1)`. Hence",
            "",
            "```text",
            exact["coefficient_envelope"],
            "```",
            "",
            "The number four appears here structurally: the worst exponent tends",
            "to `1/2+c/8`, which is summable exactly when `c>4`.",
            "",
            "## Zeta Limit",
            "",
            "For `k=0,1,2`, compare the finite weighted moments with the ordinary",
            "Dirichlet moments of zeta. Since `exp(q)-1<=q*exp(q)` for `q>=0`,",
            "",
            "```text",
            exact["moments"],
            "```",
            "",
            "The first infinite sum is finite for every fixed epsilon; the second",
            "is exponentially small because `Re(s_*)>1`. Therefore",
            "",
            "```text",
            exact["jet_convergence"],
            "```",
            "",
            "The Euler product supplies a noncircular lower bound independent of",
            "the spectral phase:",
            "",
            "```text",
            exact["zeta_floor"],
            "```",
            "",
            "Indeed `|1-p^(-s)|<=1+p^(-Re(s))`, and the resulting product is at",
            "least `1/zeta(3/2+epsilon/8)`.",
            "",
            "## Curvature",
            "",
            "Write the nonzero finite main sum as `D=a*exp(i*psi)` and set",
            "`theta=beta+psi`. Exact symbolic differentiation gives",
            "",
            "```text",
            exact["phase_amplitude"],
            "```",
            "",
            "The zeta lower bound and moment convergence imply",
            "",
            "```text",
            exact["phase_bounds"],
            "```",
            "",
            "so the `v^2` term is of order `L^2`, while every possible negative",
            "term is bounded. Thus, for a positive epsilon-dependent constant,",
            "",
            "```text",
            exact["main_curvature"],
            "```",
            "",
            "## Exact Heat Flow",
            "",
            "On a radius-`1/4` complex collar, the published `eA+eB` term is",
            "exponentially small because its bracket is `O(L^2*exp(-L))`; the",
            "`eC0` term remains `O(exp(-L/4))`; and an adjacent cutoff contributes",
            "one `n^(-p_epsilon)` block at `n` of order `exp(L/2)`. Restoring the",
            "raw normalizer costs only `exp(L/16+O(1))`. Consequently",
            "",
            "```text",
            exact["raw_remainder"],
            "```",
            "",
            "Cauchy differentiation and the exact normalized-curvature perturbation",
            "identity preserve the main sign. Combining this bounded-c branch with",
            "the already certified c>=25 ray proves",
            "",
            "```text",
            exact["theorem"],
            "```",
            "",
            "## Historical Method Boundary",
            "",
            "This absolute-value argument stops at",
            "",
            "```text",
            exact["critical_layer"],
            "```",
            "",
            "At this threshold the uniform coefficient exponent reaches one, so",
            "absolute Dirichlet control genuinely stops. Any continuation through",
            "that layer must use cancellation rather than another absolute-tail",
            "majorant.",
            "",
            "The later oscillatory-zeta handoff does exactly that and lowers the",
            "current corpus-wide boundary to",
            "",
            "```text",
            exact["current_global_boundary"],
            "```",
            "",
            "Thus the strip at and below four is not the live global gap; it is the",
            "historical boundary of this particular method.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 scaled boundary-layer asymptotic theorem: "
        f"{len(artifact['rows'])} rows, critical scaled threshold 4"
    )


if __name__ == "__main__":
    main()
