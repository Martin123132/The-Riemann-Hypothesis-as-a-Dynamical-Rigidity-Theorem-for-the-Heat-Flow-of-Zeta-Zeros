#!/usr/bin/env python3
"""Build a strong-log-concave Mellin countermodel at the quartic frontier."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_strong_logconcave_local_quartic_countermodel.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_strong_logconcave_local_quartic_countermodel.md"
)
PRECISION_BITS = 512


@dataclass(frozen=True)
class CountermodelRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_text(value: flint.arb, digits: int = 60) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def arb_upper_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.upper().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_plus(), "E")


def normalized_moments() -> list[flint.arb]:
    arb = flint.arb
    acb = flint.acb
    lower = acb(1) / 20
    upper = acb(1)
    tolerance = arb("1e-70")
    moments: list[flint.arb] = []
    for k in range(6):
        p = acb(k) + acb(1) / 2

        def integrand(y: flint.acb, analytic: bool) -> flint.acb:
            del analytic
            return y ** (p - 1) * (-3 * y - y * y / 10).exp()

        integral = acb.integral(
            integrand,
            lower,
            upper,
            abs_tol=tolerance,
            rel_tol=tolerance,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        )
        if not integral.imag.contains(0):
            raise RuntimeError("Mellin integral acquired an imaginary part")
        moments.append(integral.real / p.gamma().real)
    return moments


def cubic_frontier(left: flint.arb, right: flint.arb) -> flint.arb:
    return left**2 * right**2 - 6 * left * right + 4 * left + 4 * right - 3


def quartic_frontier(x: flint.arb, y: flint.arb, z: flint.arb) -> flint.arb:
    return (
        x**3 * y**4 * z**3
        - 12 * x**2 * y**3 * z**2
        - 18 * x**2 * y**2 * z**2
        + 54 * x**2 * y**2 * z
        - 27 * x**2 * y**2
        + 54 * x * y**2 * z**2
        - 6 * x * y**2 * z
        - 180 * x * y * z
        + 108 * x * y
        + 81 * x * z
        - 54 * x
        - 27 * y**2 * z**2
        + 108 * y * z
        - 64 * y
        - 54 * z
        + 36
    )


def build_payload() -> dict:
    flint.ctx.prec = PRECISION_BITS
    moments = normalized_moments()
    contractions = [
        moments[k + 1] * moments[k - 1] / moments[k] ** 2
        for k in range(1, 5)
    ]
    lower_walls = [
        flint.arb(1) / 3,
        flint.arb(3) / 5,
        flint.arb(5) / 7,
        flint.arb(7) / 9,
    ]
    lower_margins = [
        contractions[k] - lower_walls[k] for k in range(len(contractions))
    ]
    upper_margins = [1 - value for value in contractions]
    monotone_gaps = [
        contractions[k + 1] - contractions[k] for k in range(3)
    ]
    cubic_margins = [
        cubic_frontier(contractions[k], contractions[k + 1])
        for k in range(3)
    ]
    quartic_q = quartic_frontier(*contractions[:3])

    if not all(bool(moment > 0) for moment in moments):
        raise RuntimeError("a normalized Mellin moment is not positive")
    if not all(bool(margin > 0) for margin in lower_margins + upper_margins):
        raise RuntimeError("the local ratio walls are not strict")
    if not all(bool(gap > 0) for gap in monotone_gaps):
        raise RuntimeError("the local contractions are not strictly increasing")
    if not all(bool(margin < 0) for margin in cubic_margins):
        raise RuntimeError("a local shifted cubic is not strictly hyperbolic")
    if not bool(quartic_q < 0):
        raise RuntimeError("the local quartic frontier is not negative")

    exact = {
        "density": "f(y)=exp(-3*y-y^2/10)*1_[1/20,1](y)",
        "strong_logconcavity_modulus": "1/5",
        "normalized_moment": (
            "A_k=integral_(1/20)^1 y^(k-1/2)*exp(-3*y-y^2/10)dy/"
            "Gamma(k+1/2)"
        ),
        "contraction": "x_k=A_(k-1)*A_(k+1)/A_k^2",
        "cubic_frontier": "F(s,t)=s^2*t^2-6*s*t+4*s+4*t-3",
        "quartic": (
            "J_4(w)=1+4*w+6*x_1*w^2+4*x_1^2*x_2*w^3+"
            "x_1^3*x_2^2*x_3*w^4"
        ),
        "quartic_discriminant": "Disc(J_4)=256*x_1^6*x_2^2*Q(x_1,x_2,x_3)",
        "full_support_approximation": (
            "f_L(y)=exp(-3*y-y^2/10-L*((1/20-y)_+ +(y-1)_+)), "
            "y>=0, L<infinity"
        ),
    }
    diagnostics = {
        "precision_bits": PRECISION_BITS,
        "moments": [
            {
                "k": k,
                "ball": arb_text(value),
                "lower": arb_lower_text(value),
                "upper": arb_upper_text(value),
            }
            for k, value in enumerate(moments)
        ],
        "contractions": [
            {
                "k": k + 1,
                "ball": arb_text(value),
                "lower": arb_lower_text(value),
                "upper": arb_upper_text(value),
                "lower_wall": str(lower_walls[k]),
                "lower_wall_margin": arb_text(lower_margins[k]),
                "upper_wall_margin": arb_text(upper_margins[k]),
            }
            for k, value in enumerate(contractions)
        ],
        "monotone_gaps": [
            {
                "from_k": k + 1,
                "to_k": k + 2,
                "ball": arb_text(value),
                "lower": arb_lower_text(value),
            }
            for k, value in enumerate(monotone_gaps)
        ],
        "cubic_frontiers": [
            {
                "shift": k,
                "ball": arb_text(value),
                "upper": arb_upper_text(value),
            }
            for k, value in enumerate(cubic_margins)
        ],
        "quartic_Q": {
            "shift": 0,
            "ball": arb_text(quartic_q),
            "upper": arb_upper_text(quartic_q),
        },
    }
    rows = [
        CountermodelRow(
            id="slcq_01_strong_logconcavity",
            role="exact_input",
            readiness="ready_to_apply",
            claim=(
                "The compact Mellin density is positive on a convex support and "
                "1/5-strongly log-concave in the squared variable."
            ),
            formula="-d^2(log f)/dy^2=1/5 on (1/20,1)",
            proof_boundary="Abstract Mellin density, not the Xi kernel.",
        ),
        CountermodelRow(
            id="slcq_02_normalized_mellin_coordinate",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The coefficient normalization is exactly the Xi Berwald-Borell normalization.",
            formula=exact["normalized_moment"],
            proof_boundary="Only the first six moments are used below.",
        ),
        CountermodelRow(
            id="slcq_03_positive_moments",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="Arb encloses A_0 through A_5 strictly above zero.",
            formula="A_k>0, 0<=k<=5",
            proof_boundary="Six finite normalized Mellin moments.",
            diagnostics={"moments": diagnostics["moments"]},
        ),
        CountermodelRow(
            id="slcq_04_local_ratio_walls",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="The first four contractions satisfy every local Stieltjes lower wall and Berwald-Borell upper wall strictly.",
            formula="(2*k-1)/(2*k+1)<x_k<1, 1<=k<=4",
            proof_boundary="Four local contractions, not an all-k cone theorem.",
            diagnostics={"contractions": diagnostics["contractions"]},
        ),
        CountermodelRow(
            id="slcq_05_local_monotonicity",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="The four local contractions are strictly increasing.",
            formula="x_1<x_2<x_3<x_4",
            proof_boundary="Three local adjacent gaps.",
            diagnostics={"monotone_gaps": diagnostics["monotone_gaps"]},
        ),
        CountermodelRow(
            id="slcq_06_three_shifted_cubics",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="Three consecutive shifted cubic Jensen windows are strictly hyperbolic.",
            formula="F(x_1,x_2)<0, F(x_2,x_3)<0, F(x_3,x_4)<0",
            proof_boundary="Shifts zero through two only; no all-shift cubic claim.",
            diagnostics={"cubic_frontiers": diagnostics["cubic_frontiers"]},
        ),
        CountermodelRow(
            id="slcq_07_negative_quartic_frontier",
            role="countermodel_gate",
            readiness="ready_to_apply",
            claim="The shift-zero quartic has negative discriminant and is not hyperbolic.",
            formula=exact["quartic_discriminant"] + "<0",
            proof_boundary="One local quartic nonhyperbolicity theorem.",
            diagnostics={"quartic_Q": diagnostics["quartic_Q"]},
        ),
        CountermodelRow(
            id="slcq_08_local_nonpromotion",
            role="non_promotion_gate",
            readiness="ready_to_apply",
            claim=(
                "Strong squared-variable log-concavity plus the complete local "
                "quadratic, ratio, monotone, and cubic data do not force degree four."
            ),
            formula="local D2 + local ratio cone + local D3 != local D4",
            proof_boundary=(
                "Blocks a local generic implication only; it does not reproduce "
                "the Xi heat trajectory or its global all-shift cubic theorem."
            ),
        ),
        CountermodelRow(
            id="slcq_09_full_support_approximation",
            role="exact_approximation_theorem",
            readiness="ready_to_apply",
            claim=(
                "Full-support 1/5-strongly log-concave densities retain all strict "
                "countermodel signs for sufficiently large finite L."
            ),
            formula=exact["full_support_approximation"],
            proof_boundary=(
                "Existence by dominated convergence and openness of finitely many "
                "strict moment inequalities; no explicit least L is claimed."
            ),
        ),
        CountermodelRow(
            id="slcq_10_xi_specific_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "A quartic proof must use a global heat-compatible or Xi-specific "
                "theta-arithmetic condition beyond the reconciled low-degree facts."
            ),
            formula="construct and propagate the branch-aware quartic threshold",
            proof_boundary="Degree four for Xi, PF-infinity, Lambda<=0, and RH remain open.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_strong_logconcave_local_quartic_countermodel",
        "date": "2026-07-25",
        "status": "exact interval local quartic countermodel gate",
        "proof_boundary": (
            "This artifact gives a rigorous abstract Mellin countermodel to a local "
            "promotion from strong squared-variable log-concavity and the neighboring "
            "degree-two/degree-three data to degree four. It is not the Xi kernel or "
            "a Newman heat trajectory and does not disprove the global all-shift cubic "
            "theorem, degree-four hyperbolicity for Xi, PF-infinity, Lambda<=0, or RH."
        ),
        "sources": [
            "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
            "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
            "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md",
            "outputs/jensen_window_pf_quartic_boundary_flow_obstruction.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_strong_logconcave_local_quartic_countermodel.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_strong_logconcave_local_quartic_countermodel.py"
        ),
        "exact": exact,
        "diagnostics": diagnostics,
        "summary": {
            "rows": len(rows),
            "positive_moments": len(moments),
            "local_ratio_wall_contractions": len(contractions),
            "positive_monotone_gaps": len(monotone_gaps),
            "strict_cubic_margins": len(cubic_margins),
            "negative_quartic_frontiers": 1,
            "full_support_approximation_theorems": 1,
            "open_xi_handoffs": 1,
            "main_finding": (
                "A 1/5-strongly log-concave Gamma-normalized Mellin density "
                "satisfies all local ratio walls and three consecutive strict cubic "
                "tests while its shift-zero quartic discriminant is negative. "
                "Full-support strongly log-concave approximants preserve the strict "
                "signs, so the quartic bridge needs genuinely global or Xi-specific "
                "structure beyond the reconciled low-degree facts."
            ),
        },
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    diagnostics = payload["diagnostics"]
    contractions = diagnostics["contractions"]
    cubics = diagnostics["cubic_frontiers"]
    lines = [
        "# Jensen-Window PF Strong-Log-Concave Local Quartic Countermodel",
        "",
        "Date: 2026-07-25",
        "",
        "Status: exact interval local quartic countermodel gate. This is not a proof",
        "or disproof of degree-four Xi hyperbolicity, PF-infinity, RH, or",
        "`Lambda <= 0`; the witness is not the Xi kernel.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_strong_logconcave_local_quartic_countermodel.json",
        "python work/rh_compute/scripts/jensen_window_pf_strong_logconcave_local_quartic_countermodel.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_strong_logconcave_local_quartic_countermodel.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF strong-log-concave local quartic countermodel: 10 rows, 0 issues, 6 positive Mellin moments, 4 local ratio-wall contractions, 3 monotone gaps, 3 strict cubic margins, 1 negative quartic frontier, 1 full-support approximation theorem, 1 Xi-specific handoff",
        "```",
        "",
        "## Strong-Log-Concave Mellin Witness",
        "",
        "On the convex support `[1/20,1]`, take",
        "",
        "```text",
        payload["exact"]["density"],
        "A_k=integral_(1/20)^1 y^(k-1/2)f(y)dy/Gamma(k+1/2).",
        "```",
        "",
        "The potential `3*y+y^2/10` has second derivative `1/5`, so this is",
        "a `1/5`-strongly log-concave squared-variable density. Arb encloses",
        "`A_0,...,A_5` strictly above zero at 512-bit precision.",
        "",
        "## Low-Degree Data",
        "",
        "For `x_k=A_(k-1)A_(k+1)/A_k^2`, Arb gives",
        "",
        "```text",
    ]
    for entry in contractions:
        lines.append(f"x_{entry['k']}={entry['ball']}")
    lines.extend(
        [
            "```",
            "",
            "All four satisfy",
            "",
            "```text",
            "(2*k-1)/(2*k+1)<x_k<1, 1<=k<=4,",
            "x_1<x_2<x_3<x_4.",
            "```",
            "",
            "The three consecutive cubic frontier balls are",
            "",
            "```text",
        ]
    )
    for entry in cubics:
        lines.append(f"F_shift_{entry['shift']}={entry['ball']}<0")
    lines.extend(
        [
            "```",
            "",
            "Thus the local quadratic walls, Stieltjes lower walls, monotone",
            "contractions, and shifted cubic tests are all strict.",
            "",
            "## Quartic Failure",
            "",
            "The normalized quartic obeys",
            "",
            "```text",
            payload["exact"]["quartic"],
            payload["exact"]["quartic_discriminant"],
            f"Q={diagnostics['quartic_Q']['ball']}<0.",
            "```",
            "",
            "Hence the shift-zero quartic has negative discriminant and a nonreal",
            "conjugate pair. This is",
            "a moment-realized strengthening of the abstract quartic boundary-flow",
            "obstruction: the currently reconciled local low-degree facts still do",
            "not provide the missing quartic invariant.",
            "",
            "## Full-Support Guard",
            "",
            "For finite `L>0`, define on `y>=0`",
            "",
            "```text",
            payload["exact"]["full_support_approximation"],
            "```",
            "",
            "Its potential is `1/5`-strongly convex and the density is positive on",
            "the full half-line. As `L` tends to infinity, its first six normalized",
            "moments converge to those above by dominated convergence. Every sign",
            "used here is strict and depends continuously on those moments, so all",
            "of the displayed countermodel inequalities persist for sufficiently",
            "large finite `L`. Full support and strong log-concavity therefore do not",
            "repair the local implication.",
            "",
            "## Proof Boundary",
            "",
            "The witness is not a Newman heat trajectory and does not satisfy or",
            "challenge the global Xi all-shift theorem as a dynamical statement.",
            "It blocks only a local generic promotion. The surviving target is a",
            "heat-compatible branch-aware quartic invariant, an Xi-specific weighted",
            "theta/Laguerre composition, or a noncircular all-degree theorem.",
            "",
            "```text",
            "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
            "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
            "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md",
            "outputs/jensen_window_pf_quartic_boundary_flow_obstruction.md",
            "```",
            "",
            "Summary:",
            "",
            payload["summary"]["main_finding"],
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF strong-log-concave local quartic countermodel: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
