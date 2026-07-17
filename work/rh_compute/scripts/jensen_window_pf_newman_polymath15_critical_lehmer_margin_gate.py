#!/usr/bin/env python3
"""Stress the critical Newman curvature margin at an explicit Lehmer point."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout import (  # noqa: E402
    compute_point,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.md"
)
X0 = Fraction(1_401_016_343, 100_000)
C_VALUES = ("0", "0.01", "0.1", "0.5")
PRECISION_BITS = 256


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def m_zero(s: flint.acb) -> flint.acb:
    pi = flint.arb.pi()
    return (
        (s * (s - 1) / 2)
        * (-s * pi.log() / 2).exp()
        * (2 * pi).sqrt()
        * ((s / 2 - flint.acb(0.5)) * (s / 2).log() - s / 2).exp()
        / 8
    )


def arb_exact_point() -> dict:
    """Enclose exact H_0 curvature at the explicit rational point X0."""

    flint.ctx.prec = PRECISION_BITS
    x = flint.arb(flint.fmpq(X0.numerator, X0.denominator))
    imag = flint.acb(0, 1)
    s0 = (flint.acb(1) + imag * flint.acb(x)) / 2
    s = flint.acb_series([s0, imag / 2, 0, 0], 4)
    xi = (
        s
        * (s - 1)
        * (-s * flint.arb.pi().log() / 2).exp()
        * (s / 2).gamma()
        * s.zeta()
        / 2
    )
    h0 = xi[0].real / 8
    h1 = xi[1].real / 8
    h2 = 2 * xi[2].real / 8
    h3 = 6 * xi[3].real / 8
    normalizer = m_zero((flint.acb(1) - imag * flint.acb(x)) / 2)
    amplitude_squared = (normalizer * normalizer.conjugate()).real
    if not amplitude_squared > 0:
        raise RuntimeError("exact Lehmer-point normalizer amplitude is not positive")
    curvature = (h1 * h1 - h0 * h2) / amplitude_squared
    monotonicity_margin = (h0 * h3 - h1 * h2) / amplitude_squared
    monotonicity_to_curvature = monotonicity_margin / curvature
    ell = (x / (4 * flint.arb.pi())).log()
    scaled = curvature / (ell * ell)
    if not curvature > 0:
        raise RuntimeError("exact Lehmer-point curvature positivity failed")
    if not scaled < flint.arb(1) / 1000:
        raise RuntimeError("exact Lehmer-point small-margin bound failed")
    if not scaled < flint.arb(3) / 40:
        raise RuntimeError("dominant-ray margin separation failed")
    if not monotonicity_margin < 0:
        raise RuntimeError("exact Lehmer-point monotonicity counterexample failed")
    return {
        "precision_bits": PRECISION_BITS,
        "x": f"{X0.numerator}/{X0.denominator}",
        "L_ball": arb_text(ell),
        "H0_ball": arb_text(h0),
        "H1_ball": arb_text(h1),
        "H2_ball": arb_text(h2),
        "H3_ball": arb_text(h3),
        "amplitude_squared_ball": arb_text(amplitude_squared),
        "exact_normalized_curvature_ball": arb_text(curvature),
        "exact_normalized_curvature_mid": arb_text(curvature.mid()),
        "scaled_curvature_ball": arb_text(scaled),
        "exact_normalized_monotonicity_margin_ball": arb_text(monotonicity_margin),
        "monotonicity_to_curvature_ball": arb_text(monotonicity_to_curvature),
        "certified_bounds": "0<C[H_0/A_0]/L^2<1/1000<3/40",
        "certified_monotonicity_bounds": (
            "M_0(x)=H_0(x)H_0'''(x)-H_0'(x)H_0''(x)=-L_0'(x)<0"
        ),
    }


def build_exact() -> dict:
    a, g0, g1, g2, g3 = sp.symbols("a g0 g1 g2 g3", real=True)
    f0 = -a**2 * g0
    f1 = -a**2 * g1
    f2 = 2 * g0 - a**2 * g2
    laguerre = sp.expand(f1**2 - f0 * f2)
    pair_expected = 2 * a**2 * g0**2 + a**4 * (g1**2 - g0 * g2)
    if sp.simplify(laguerre - pair_expected) != 0:
        raise RuntimeError("close-pair Laguerre identity failed")
    f3 = 6 * g1 - a**2 * g3
    monotonicity = sp.expand(f0 * f3 - f1 * f2)
    monotonicity_expected = -4 * a**2 * g0 * g1 + a**4 * (
        g0 * g3 - g1 * g2
    )
    if sp.simplify(monotonicity - monotonicity_expected) != 0:
        raise RuntimeError("close-pair monotonicity identity failed")

    h0, h1, h2, h3, h4 = sp.symbols("h0 h1 h2 h3 h4", real=True)
    heat_derivative = h2**2 - 2 * h1 * h3 + h0 * h4
    if sp.simplify(heat_derivative.subs({h0: 0, h1: 0}) - h2**2) != 0:
        raise RuntimeError("double-zero heat transversality failed")

    y, tau = sp.symbols("y tau", real=True)
    model = y**2 - 2 * tau
    model_laguerre = sp.expand(sp.diff(model, y) ** 2 - model * sp.diff(model, y, 2))
    if sp.simplify(model_laguerre - (2 * y**2 + 4 * tau)) != 0:
        raise RuntimeError("quadratic heat model failed")
    return {
        "close_pair": (
            "If F(m+y)=(y^2-a^2)G(m+y), then "
            "L[F](m)=2*a^2*G(m)^2+a^4*(G'(m)^2-G(m)G''(m))"
        ),
        "gap_scale": "a=(gamma_(j+1)-gamma_j)/2, so the leading margin is quadratic in the gap",
        "close_pair_monotonicity": (
            "If F(m+y)=(y^2-a^2)G(m+y), then M[F](m)="
            "-4*a^2*G(m)*G'(m)+a^4*(G(m)*G'''(m)-G'(m)*G''(m))"
        ),
        "heat_derivative": (
            "For partial_t H=-partial_x^2 H, partial_t L="
            "H_xx^2-2H_x H_xxx+H H_xxxx; at H=H_x=0 this equals H_xx^2"
        ),
        "quadratic_model": (
            "H_tau(y)=y^2-2*tau has L[H_tau](y)=2*y^2+4*tau"
        ),
        "forbidden_promotion": (
            "No c-independent positive multiple of L^2 obtained on the dominant "
            "ray may be extended unchanged through the critical layer"
        ),
        "monotonicity_continuity": (
            "The theta tail makes |u|^3*exp(u^2/5)*Phi(u) integrable, so "
            "dominated convergence makes H_t and its first three x-derivatives "
            "continuous in t on [0,1/5]. Hence M_0(x)<0 implies M_t(x)<0 "
            "for every sufficiently small positive t."
        ),
    }


def numerical_scout(coarse_dps: int, fine_dps: int) -> dict:
    x_text = f"{X0.numerator / X0.denominator:.5f}"
    # Let compute_point parse the decimal after setting its working precision.
    coarse = [compute_point(x_text, c, coarse_dps) for c in C_VALUES]
    fine = [compute_point(x_text, c, fine_dps) for c in C_VALUES]
    rows = []
    max_delta = mp.mpf(0)
    for left, right in zip(coarse, fine, strict=True):
        lc = mp.mpf(left["corrected_curvature"])
        rc = mp.mpf(right["corrected_curvature"])
        max_delta = max(max_delta, abs(lc - rc) / abs(rc))
        ell = mp.mpf(right["L"])
        jet = [mp.mpf(value) for value in right["corrected_main_jet"]]
        transversality = jet[0] ** 2 + (jet[1] / ell) ** 2
        row = {
            "c": right["c"],
            "t": right["t"],
            "N": right["N"],
            "p": right["p"],
            "corrected_curvature": right["corrected_curvature"],
            "corrected_curvature_over_L2": mp.nstr(rc / ell**2, 35),
            "corrected_value": right["corrected_main_jet"][0],
            "corrected_first_over_L": mp.nstr(jet[1] / ell, 35),
            "transversality_norm": mp.nstr(transversality, 35),
        }
        if right["c"] == "0":
            row["exact_xi_curvature"] = right["exact_xi_curvature"]
            row["relative_corrected_to_exact_delta"] = right[
                "relative_corrected_to_exact_delta"
            ]
        rows.append(row)
    if max_delta >= mp.mpf("1e-25"):
        raise RuntimeError("Lehmer stress scout precision failed")
    if any(mp.mpf(row["corrected_curvature"]) <= 0 for row in rows):
        raise RuntimeError("Lehmer stress scout sign failed")
    return {
        "coarse_dps": coarse_dps,
        "fine_dps": fine_dps,
        "max_relative_corrected_curvature_delta": mp.nstr(max_delta, 25),
        "rows": rows,
    }


def build_artifact(coarse_dps: int, fine_dps: int) -> dict:
    exact = build_exact()
    arb_data = arb_exact_point()
    scout = numerical_scout(coarse_dps, fine_dps)
    rows = [
        GateRow(
            id="np15clmg_01_exact_pair_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The first Laguerre margin near a close real pair is quadratically gap-sensitive.",
            formula=exact["close_pair"],
            proof_boundary="Exact local factorization identity; the regular factor must be retained.",
        ),
        GateRow(
            id="np15clmg_02_heat_transversality",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="A nondegenerate double zero is crossed transversally by Newman heat time.",
            formula=exact["heat_derivative"],
            proof_boundary="Uses the Newman sign partial_t H=-partial_x^2 H.",
        ),
        GateRow(
            id="np15clmg_03_quadratic_model",
            role="exact_model",
            readiness="ready_to_apply",
            claim="The universal quadratic collision model has a time-dependent, not uniform, curvature floor.",
            formula=exact["quadratic_model"],
            proof_boundary="Local model only; it does not exclude a positive collision time.",
        ),
        GateRow(
            id="np15clmg_04_arb_exact_point",
            role="interval_certificate",
            readiness="certified",
            claim="Arb certifies a strict but very small exact xi curvature margin at one explicit rational point.",
            formula=arb_data["certified_bounds"],
            proof_boundary="A single exact t=0 point; no interval or RH promotion.",
            diagnostics=arb_data,
        ),
        GateRow(
            id="np15clmg_05_fixed_floor_rejected",
            role="counter_gate",
            readiness="closed",
            claim="The dominant-ray 3/40 times L-squared margin cannot be continued unchanged to c=0.",
            formula="C[H_0/A_0]/L^2<1/1000<3/40 at x=14010.16343",
            proof_boundary="Rejects only fixed-floor continuation, not a gap- or time-adaptive estimate.",
        ),
        GateRow(
            id="np15clmg_06_corrected_scout",
            role="high_precision_diagnostic",
            readiness="diagnostic_only",
            claim="The corrected Polymath-15 finite main tracks the small-margin regime stably.",
            formula="c in {0,0.01,0.1,0.5}; corrected curvature and first-jet norm recorded",
            proof_boundary="Point diagnostics; the finite main is not the exact heat flow for c>0.",
            diagnostics=scout,
        ),
        GateRow(
            id="np15clmg_07_target_refinement",
            role="strategy_gate",
            readiness="ready_to_apply",
            claim="Critical-layer transfer must be adaptive to near-collisions or target double-zero transversality directly.",
            formula="control (J,J'/L) against a C1 remainder instead of demanding C[J]>=kappa*L^2",
            proof_boundary="A strategy consequence, not a proof of the required lower bound.",
        ),
        GateRow(
            id="np15clmg_08_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The positive point values do not certify a uniform sign or exclude a positive Newman boundary.",
            formula=exact["forbidden_promotion"],
            proof_boundary="Not Lambda<=0 or RH.",
        ),
        GateRow(
            id="np15clmg_09_monotonicity_counterexample",
            role="interval_counterexample",
            readiness="certified",
            claim="Arb rejects global strict decrease of the Xi first-Laguerre expression at the Lehmer stress point.",
            formula=arb_data["certified_monotonicity_bounds"],
            proof_boundary=(
                "Continuity rejects the universal assertion M_t>0 for all x and "
                "all 0<t<=1/5, but does not reject strict Laguerre positivity or "
                "the transversality route."
            ),
            diagnostics={
                "x": arb_data["x"],
                "normalized_margin_ball": arb_data[
                    "exact_normalized_monotonicity_margin_ball"
                ],
                "margin_to_curvature_ball": arb_data[
                    "monotonicity_to_curvature_ball"
                ],
                "continuity_deduction": exact["monotonicity_continuity"],
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate",
        "date": "2026-07-17",
        "status": (
            "rigorous explicit small-margin and strict-monotonicity counter-gates "
            "plus corrected finite diagnostics; not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact rigorously proves a small positive normalized t=0 "
            "curvature at one explicit rational x and exact local close-pair "
            "identities. It rejects an unchanged fixed L-squared margin in the "
            "critical layer. It also certifies M_0<0 there, which by continuity "
            "rejects global strict decrease for all positive target times. It does "
            "not reject strict Laguerre positivity, certify a frequency interval, "
            "exclude double zeros, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "arb": arb_data,
        "scout": scout,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.md",
        ],
    }


def render_note(artifact: dict) -> str:
    arb_data = artifact["arb"]
    table = [
        "| c=tL | N | C[J]/L^2 | J | J'/L | J^2+(J'/L)^2 |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in artifact["scout"]["rows"]:
        table.append(
            "| {c} | {N} | {curv} | {j0} | {j1} | {trans} |".format(
                c=row["c"],
                N=row["N"],
                curv=mp.nstr(mp.mpf(row["corrected_curvature_over_L2"]), 9),
                j0=mp.nstr(mp.mpf(row["corrected_value"]), 9),
                j1=mp.nstr(mp.mpf(row["corrected_first_over_L"]), 9),
                trans=mp.nstr(mp.mpf(row["transversality_norm"]), 9),
            )
        )
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical Lehmer Margin Gate",
            "",
            "Date: 2026-07-17",
            "",
            "Status: rigorous small-margin counter-gate and corrected finite",
            "diagnostics. This is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.py",
            "```",
            "",
            "## Exact Stress Point",
            "",
            "At the explicit rational point `x=14010.16343`, Arb evaluates the",
            "exact `H_0=xi((1+i*x)/2)/8` through third order and divides its",
            "Laguerre expression by the exact normalizer amplitude. It certifies",
            "",
            "```text",
            arb_data["exact_normalized_curvature_ball"],
            arb_data["scaled_curvature_ball"],
            arb_data["certified_bounds"],
            "```",
            "",
            "Thus the `3/40*L^2` margin proved on `tL>=25` cannot simply be",
            "continued into the critical layer.",
            "",
            "## Strict Monotonicity Counter-Gate",
            "",
            "The same Arb jet certifies",
            "",
            "```text",
            arb_data["exact_normalized_monotonicity_margin_ball"],
            arb_data["monotonicity_to_curvature_ball"],
            arb_data["certified_monotonicity_bounds"],
            "```",
            "",
            "The theta tail makes `|u|^3*exp(u^2/5)*Phi(u)` integrable, so",
            "dominated convergence makes the heat-flow jet through third order",
            "continuous in `t` on `[0,1/5]`. The strict negative value at `t=0` therefore",
            "persists for all sufficiently small positive `t`. This rigorously",
            "retires the proposed global condition `M_t(x)>0` for every `x>0` and",
            "every `0<t<=1/5`. It does not challenge `L_t(x)>0` itself.",
            "",
            "## Why The Margin Shrinks",
            "",
            "For a close pair with local factorization",
            "`F(m+y)=(y^2-a^2)G(m+y)`, exact differentiation gives",
            "",
            "```text",
            artifact["exact"]["close_pair"],
            artifact["exact"]["close_pair_monotonicity"],
            "```",
            "",
            "The leading curvature is quadratic in the half-gap, while the stronger",
            "monotonicity sign also depends on the slope of the regular factor and has",
            "no fixed close-pair sign. Under the Newman",
            "heat sign, a double zero is nevertheless transversal:",
            "",
            "```text",
            artifact["exact"]["heat_derivative"],
            artifact["exact"]["quadratic_model"],
            "```",
            "",
            "## Corrected Main",
            "",
            "Independent 60- and 80-digit runs of the corrected Polymath-15 main",
            "give:",
            "",
            *table,
            "",
            "The numerical rows are diagnostics. Their role is to show the shape",
            "of the viable target: a fixed curvature floor is too blunt near close",
            "pairs, while the two-component first-jet norm directly measures",
            "distance from a double zero and requires only a `C1` error transfer.",
            "",
            "No finite or point computation here is promoted to a uniform sign,",
            "absence of positive-time collisions, `Lambda <= 0`, or RH.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--coarse-dps", type=int, default=60)
    parser.add_argument("--fine-dps", type=int, default=80)
    args = parser.parse_args()
    artifact = build_artifact(args.coarse_dps, args.fine_dps)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 critical Lehmer margin gate: "
        f"{len(artifact['rows'])} rows, 1 Arb small-margin point, "
        "3 exact local identities, 1 fixed-floor rejection, "
        "1 Arb strict-monotonicity rejection"
    )


if __name__ == "__main__":
    main()
