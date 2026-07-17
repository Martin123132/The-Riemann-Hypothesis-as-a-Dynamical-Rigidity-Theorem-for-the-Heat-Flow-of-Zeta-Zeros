#!/usr/bin/env python3
"""Promote uniform strong log-concavity through the Newman target window."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402
import sympy as sp  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_time_strong_logconcavity_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.md"
)

CERTIFICATE_PRECISION_BITS = 256
TAIL_START = 9
TARGET_TIME = sp.Rational(1, 5)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def arb_text(value: flint.arb, digits: int = 55) -> str:
    return value.str(digits).replace("e", "E")


def curvature_certificate(
    *,
    precision_bits: int = CERTIFICATE_PRECISION_BITS,
    tail_start: int = TAIL_START,
) -> dict:
    """Certify Phi(0), Phi''(0), and kappa=-Phi''(0)/Phi(0)."""
    if tail_start < 3:
        raise ValueError("tail_start must be at least three")

    old_precision = flint.ctx.prec
    flint.ctx.prec = precision_bits
    try:
        arb = flint.arb
        pi = arb.pi()
        finite_phi0 = arb(0)
        finite_phi2 = arb(0)
        for n in range(1, tail_start):
            a = pi * n * n
            exponential = (-a).exp()
            finite_phi0 += a * (2 * a - 3) * exponential
            finite_phi2 += (
                a
                * (32 * a**3 - 224 * a**2 + 330 * a - 75)
                * exponential
            )

        # For n >= tail_start, the Phi terms are positive and bounded by
        # 2*pi^2*n^4*exp(-pi*n^2).  The Phi'' terms are positive and bounded
        # by 362*pi^4*n^8*exp(-pi*n^2).  Consecutive ratios decrease in n.
        n0 = arb(tail_start)
        ratio_phi0 = (
            ((n0 + 1) / n0) ** 4 * (-pi * (2 * n0 + 1)).exp()
        )
        ratio_phi2 = (
            ((n0 + 1) / n0) ** 8 * (-pi * (2 * n0 + 1)).exp()
        )
        if not (ratio_phi0 < 1 and ratio_phi2 < 1):
            raise RuntimeError("theta-tail geometric ratio did not certify below one")

        phi0_tail_upper = (
            2
            * pi**2
            * n0**4
            * (-pi * n0**2).exp()
            / (1 - ratio_phi0)
        )
        phi2_tail_upper = (
            362
            * pi**4
            * n0**8
            * (-pi * n0**2).exp()
            / (1 - ratio_phi2)
        )

        # Symmetric error balls are slightly weaker than the proved one-sided
        # tails, but keep the interval construction simple and fully rigorous.
        phi0 = finite_phi0 + arb(0, phi0_tail_upper.upper())
        phi2 = finite_phi2 + arb(0, phi2_tail_upper.upper())
        if not (phi0 > 0 and phi2 < 0):
            raise RuntimeError("origin Phi signs did not certify")

        kappa = -phi2 / phi0
        threshold = kappa / 2
        target = arb(1) / 5
        margin = threshold - target
        if not margin > 0:
            raise RuntimeError("strong-log-concavity threshold missed 1/5")

        return {
            "arithmetic": "python-flint Arb balls",
            "precision_bits": precision_bits,
            "tail_start": tail_start,
            "phi0_finite_ball": arb_text(finite_phi0),
            "phi2_finite_ball": arb_text(finite_phi2),
            "phi0_tail_upper": arb_text(phi0_tail_upper, 30),
            "phi2_tail_upper": arb_text(phi2_tail_upper, 30),
            "phi0_ball": arb_text(phi0),
            "phi2_ball": arb_text(phi2),
            "kappa_ball": arb_text(kappa),
            "threshold_kappa_over_2_ball": arb_text(threshold),
            "target_time": "1/5",
            "threshold_margin_ball": arb_text(margin),
            "target_strong_curvature_lower_ball": arb_text(kappa - 2 * target),
            "correlation_strong_curvature_lower_ball": arb_text(
                2 * (kappa - 2 * target)
            ),
            "passed": True,
        }
    finally:
        flint.ctx.prec = old_precision


def exact_algebra() -> dict:
    u, a = sp.symbols("u a", positive=True)
    summand = (
        2 * a**2 * sp.exp(9 * u) - 3 * a * sp.exp(5 * u)
    ) * sp.exp(-a * sp.exp(4 * u))
    value_at_zero = sp.factor(summand.subs(u, 0))
    second_at_zero = sp.factor(sp.diff(summand, u, 2).subs(u, 0))
    if value_at_zero != a * (2 * a - 3) * sp.exp(-a):
        raise RuntimeError("Phi summand value formula failed")
    if second_at_zero != (
        a * (32 * a**3 - 224 * a**2 + 330 * a - 75) * sp.exp(-a)
    ):
        raise RuntimeError("Phi summand second-derivative formula failed")

    ell1, ell2 = sp.symbols("ell1 ell2", real=True)
    q_second = (u * ell2 - ell1) / (4 * u**3)
    if sp.simplify(4 * u**3 * q_second - (u * ell2 - ell1)) != 0:
        raise RuntimeError("root-variable curvature identity failed")

    p, q, alpha, beta, n, s = sp.symbols(
        "p q alpha beta n s", real=True, positive=True
    )
    hessian_form = (
        (alpha + beta - 2 * n / s**2) * p**2
        + 2 * (alpha - beta) * p * q
        + (alpha + beta) * q**2
    )
    pair_form = alpha * (p + q) ** 2 + beta * (p - q) ** 2 - 2 * n * p**2 / s**2
    if sp.expand(hessian_form - pair_form) != 0:
        raise RuntimeError("correlation Hessian identity failed")

    return {
        "published_input": {
            "statement": (
                "q(r)=log(Phi(sqrt(r))) is strictly concave for r>0"
            ),
            "primary_source": "https://doi.org/10.1007/BF02075457",
            "secondary_source": "https://arxiv.org/abs/1309.0055",
        },
        "root_variable_identity": (
            "q''(u^2)=(u*(log Phi)''(u)-(log Phi)'(u))/(4*u^3)<0"
        ),
        "uniform_phi_curvature": (
            "With kappa=-(log Phi)''(0)=-Phi''(0)/Phi(0), "
            "(log Phi)''(u)<=-kappa for every real u"
        ),
        "deformed_kernel": {
            "definition": "phi_t(u)=exp(t*u^2)*Phi(u)",
            "curvature": (
                "(log phi_t)''(u)<=(log Phi)''(0)+2*t=-(kappa-2*t)"
            ),
            "admissibility": (
                "For 0<=t<kappa/2, phi_t is positive, even, strictly decreasing "
                "on (0,infinity), strongly log-concave, and super-exponentially decaying"
            ),
        },
        "correlation": {
            "definition": (
                "K_(n,t)(v)=integral_R phi_t(s+v)*phi_t(s-v)*s^(2n)ds"
            ),
            "half_line": (
                "K_(n,t)(v)=2*integral_(s>0) "
                "phi_t(s+v)*phi_t(s-v)*s^(2n)ds"
            ),
            "joint_hessian": (
                "D^2 log integrand on (s,v) has quadratic form "
                "a*(p+q)^2+b*(p-q)^2-2*n*p^2/s^2, "
                "a=(log phi_t)''(s+v), b=(log phi_t)''(s-v)"
            ),
            "joint_bound": (
                "D^2 log integrand <=-2*(kappa-2*t)*I for 0<=t<kappa/2"
            ),
            "marginal_theorem": (
                "Prekopa marginalization gives (log K_(n,t))''(v)"
                "<=-2*(kappa-2*t) for every n>=0"
            ),
            "marginal_source": "https://doi.org/10.1016/0022-1236(76)90004-5",
        },
        "zeroth_correlation_guard": {
            "identity": "Fourier[K_(0,t)](xi)=2*H_t(xi/2)^2",
            "properties_at_t0": (
                "K_(0,0) is positive, even, strongly log-concave, positive definite, "
                "and has the same Xi super-Gaussian correlation tail"
            ),
            "zero_source": (
                "Hardy's theorem gives infinitely many real zeros of H_0, so the "
                "Fourier transform of K_(0,0) has real double zeros"
            ),
            "hardy_source": "https://arxiv.org/abs/1712.08435",
        },
        "nonpromotion_decision": (
            "Uniform strong log-concavity, admissibility, positive definiteness, and "
            "the exact Xi correlation tail do not force Fourier zero-freeness: K_(0,0) "
            "has all of them and still has double Fourier zeros. Any proof for K_(1,t) "
            "must use the s^2 weighting or a coupling in the correlation hierarchy."
        ),
        "open_handoff": (
            "Exploit structure that distinguishes K_(1,t) from K_(0,t), such as a "
            "relative F_2/F_1 curvature inequality, an s^2-weighted modular square, "
            "or a hierarchy coercivity estimate. Shape-only Fourier criteria are closed."
        ),
        "symbolic_checks": {
            "phi_summand_at_zero": str(value_at_zero),
            "phi_summand_second_at_zero": str(second_at_zero),
            "root_curvature_numerator": str(sp.factor(4 * u**3 * q_second)),
            "correlation_hessian_pair_form": str(pair_form),
        },
    }


def build_payload() -> dict:
    exact = exact_algebra()
    certificate = curvature_certificate()
    rows = [
        GateRow(
            id="npslc_01_published_root_concavity",
            role="published_theorem",
            readiness="ready_to_apply",
            claim="A published theorem gives strict concavity of log(Phi(sqrt(r))).",
            formula=exact["published_input"]["statement"],
            proof_boundary="Published input theorem; not a new proof of that source result.",
            diagnostics=exact["published_input"],
        ),
        GateRow(
            id="npslc_02_uniform_phi_curvature",
            role="exact_lemma",
            readiness="available_exact",
            claim="The root-variable theorem yields a uniform strong-log-concavity constant for Phi.",
            formula=exact["uniform_phi_curvature"],
            proof_boundary="Exact calculus consequence of the published theorem only.",
        ),
        GateRow(
            id="npslc_03_origin_curvature_certificate",
            role="interval_certificate",
            readiness="interval_validated",
            claim="Arb certifies the origin curvature and its positive margin beyond the target heat time.",
            formula=(
                "kappa=-Phi''(0)/Phi(0), kappa/2>1/5"
            ),
            proof_boundary="Finite theta sum plus explicit geometric tail enclosure only.",
            diagnostics=certificate,
        ),
        GateRow(
            id="npslc_04_positive_time_admissibility",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every deformed Xi kernel in the target window is uniformly strongly log-concave and admissible.",
            formula=exact["deformed_kernel"]["curvature"],
            proof_boundary="Kernel-shape theorem only; it does not control Fourier zeros.",
            diagnostics=exact["deformed_kernel"],
        ),
        GateRow(
            id="npslc_05_correlation_marginal",
            role="published_theorem_composition",
            readiness="ready_to_apply",
            claim="Prekopa marginalization propagates a uniform strong-log-concavity bound to every Xi correlation.",
            formula=exact["correlation"]["marginal_theorem"],
            proof_boundary="Exact composition with Prekopa; no Fourier zero-freeness follows.",
            diagnostics=exact["correlation"],
        ),
        GateRow(
            id="npslc_06_zeroth_square_transform",
            role="exact_identity",
            readiness="available_exact",
            claim="The zeroth correlation transform is exactly the square of the Newman heat function.",
            formula=exact["zeroth_correlation_guard"]["identity"],
            proof_boundary="Exact Fourier normalization identity only.",
            diagnostics=exact["zeroth_correlation_guard"],
        ),
        GateRow(
            id="npslc_07_xi_shape_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Even Xi-specific strong shape, positive definiteness, and super-Gaussian tail do not force transform zero-freeness.",
            formula=exact["nonpromotion_decision"],
            proof_boundary="K_(0,0) counterexample; it does not decide the s^2-weighted K_(1,t) target.",
            diagnostics=exact["zeroth_correlation_guard"],
        ),
        GateRow(
            id="npslc_08_target_window_margin",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The full published target interval 0<=t<=1/5 lies deep inside the strong-log-concavity regime.",
            formula=(
                "kappa-2/5>0 and (log K_(n,t))''<=-2*(kappa-2*t)"
            ),
            proof_boundary="Uniform shape margin only; not strict Fourier positivity.",
            diagnostics={
                "target_time": certificate["target_time"],
                "kernel_margin": certificate["target_strong_curvature_lower_ball"],
                "correlation_margin": certificate[
                    "correlation_strong_curvature_lower_ball"
                ],
            },
        ),
        GateRow(
            id="npslc_09_weighted_hierarchy_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The surviving route must exploit the s^2 weight or correlation-hierarchy coupling.",
            formula=exact["open_handoff"],
            proof_boundary="Open coercivity target; not RH or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_positive_time_strong_logconcavity_gate",
        "date": "2026-07-17",
        "status": (
            "exact positive-time strong-log-concavity theorem with Xi-specific nonpromotion gate"
        ),
        "proof_boundary": (
            "This artifact composes the published strict concavity of log(Phi(sqrt(r))) "
            "with an Arb origin-curvature certificate and Prekopa marginalization. It "
            "proves uniform strong log-concavity and admissibility of phi_t and every "
            "K_(n,t) throughout 0<=t<=1/5. The exact K_(0,0) square-transform guard "
            "shows that these shape and tail properties do not force Fourier "
            "zero-freeness. It does not prove strict positivity of Fourier[K_(1,t)], "
            "Wiener density, Lambda<=0, or RH."
        ),
        "sources": [
            "https://doi.org/10.1007/BF02075457",
            "https://arxiv.org/abs/1309.0055",
            "https://doi.org/10.1016/0022-1236(76)90004-5",
            "https://arxiv.org/abs/1712.08435",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.md",
        ],
        "exact": exact,
        "certificate": certificate,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = payload["certificate"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Positive-Time Strong Log-Concavity Gate",
            "",
            "Date: 2026-07-17",
            "",
            "Status: exact positive-time strong-log-concavity theorem with an",
            "Xi-specific nonpromotion gate. This is not a proof of RH or",
            "`Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_positive_time_strong_logconcavity_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Newman positive-time strong-log-concavity gate: 9 rows, 0 issues, 1 published input, 2 exact curvature/admissibility theorems, 1 Arb threshold certificate, 1 Prekopa correlation theorem, 1 square-transform identity, 1 Xi-specific nonpromotion gate, 1 target-window margin, 1 weighted-hierarchy handoff",
            "```",
            "",
            "## Published Input",
            "",
            "Csordas and Varga proved",
            "",
            "```text",
            exact["published_input"]["statement"],
            "```",
            "",
            "Primary source: https://doi.org/10.1007/BF02075457.",
            "",
            "## Uniform Curvature",
            "",
            "Writing `u=sqrt(r)` gives",
            "",
            "```text",
            exact["root_variable_identity"],
            exact["uniform_phi_curvature"],
            exact["deformed_kernel"]["curvature"],
            "```",
            "",
            "Thus `phi_t` is positive, even, strictly decreasing,",
            "super-exponentially decaying, and uniformly strongly log-concave",
            "whenever `0<=t<kappa/2`.",
            "",
            "## Arb Threshold",
            "",
            "The origin values use the exact theta-summand formulas",
            "",
            "```text",
            "phi_n(0)=a*(2*a-3)*exp(-a)",
            "phi_n''(0)=a*(32*a^3-224*a^2+330*a-75)*exp(-a)",
            "a=pi*n^2",
            f"Phi(0)={certificate['phi0_ball']}",
            f"Phi''(0)={certificate['phi2_ball']}",
            f"kappa={certificate['kappa_ball']}",
            f"kappa/2={certificate['threshold_kappa_over_2_ball']}",
            f"kappa/2-1/5={certificate['threshold_margin_ball']}",
            "```",
            "",
            "The finite sum ends at `n=8`; explicit decreasing-ratio geometric",
            "bounds enclose every omitted `n>=9` term.",
            "",
            "## Correlation Propagation",
            "",
            "On the positive half-line in `s`, the logarithmic Hessian has",
            "quadratic form",
            "",
            "```text",
            exact["correlation"]["joint_hessian"],
            exact["correlation"]["joint_bound"],
            exact["correlation"]["marginal_theorem"],
            "```",
            "",
            "This follows by Prekopa marginalization. In particular, throughout",
            "`0<=t<=1/5`, every correlation is much more strongly log-concave",
            "than a merely qualitative admissibility argument records.",
            "",
            "## Nonpromotion Gate",
            "",
            "The same shape package does not imply Fourier zero-freeness:",
            "",
            "```text",
            exact["zeroth_correlation_guard"]["identity"],
            exact["zeroth_correlation_guard"]["properties_at_t0"],
            exact["zeroth_correlation_guard"]["zero_source"],
            "```",
            "",
            exact["nonpromotion_decision"],
            "",
            "## Live Handoff",
            "",
            exact["open_handoff"],
            "",
        ]
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Newman positive-time strong-log-concavity gate: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
