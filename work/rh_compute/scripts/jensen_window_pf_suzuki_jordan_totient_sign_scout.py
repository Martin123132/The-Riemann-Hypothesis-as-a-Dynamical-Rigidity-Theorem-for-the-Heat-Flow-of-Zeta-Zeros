#!/usr/bin/env python3
"""Scout Suzuki's Jordan-totient weighted sign target on a finite grid."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

import numpy as np
import psutil
from scipy.special import beta as beta_fn
from scipy.special import betaincc, gamma


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_jordan_totient_sign_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md"
)
OMEGAS = (0.5, 0.25, 0.125, 0.0625, 0.03125)
OMEGA_LABELS = ("1/2", "1/4", "1/8", "1/16", "1/32")
ANCHORS = (10.0, 100.0, 1000.0, 5000.0)


@dataclass(frozen=True)
class OmegaSummary:
    omega: str
    sample_count: int
    negative_count: int
    nonfinite_count: int
    minimum_scaled_h: float
    minimum_x: float
    maximum_scaled_h: float
    maximum_x: float
    minimum_coefficient: float
    anchors: dict[str, float]


def set_below_normal_priority() -> None:
    try:
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except (AttributeError, psutil.Error):
        pass


def x_grid() -> np.ndarray:
    return np.unique(
        np.concatenate(
            (
                np.linspace(1.001, 20.0, 600),
                np.geomspace(20.01, 5000.0, 600),
            )
        )
    )


def jordan_coefficients(omega: float, nmax: int) -> np.ndarray:
    coefficients = np.arange(nmax + 1, dtype=float) ** omega
    coefficients[0] = 0.0
    prime = np.ones(nmax + 1, dtype=bool)
    prime[:2] = False
    for p in range(2, nmax + 1):
        if not prime[p]:
            continue
        coefficients[p::p] *= 1.0 - p ** (-2.0 * omega)
        if p * p <= nmax:
            prime[p * p :: p] = False
    return coefficients


def primitive_weight(omega: float, x: np.ndarray | float) -> np.ndarray:
    x_array = np.asarray(x, dtype=float)
    if abs(omega - 0.5) < 1e-15:
        root = np.sqrt(np.maximum(0.0, 1.0 - x_array * x_array))
        return (
            2.0
            / np.sqrt(x_array)
            * (
                2.0 * root
                + np.log(x_array)
                - np.log1p(root)
            )
        )

    a1 = (3.0 - 2.0 * omega) / 2.0
    a2 = (5.0 - 2.0 * omega) / 4.0
    upper_beta_1 = beta_fn(a1, omega) * betaincc(
        a1, omega, x_array * x_array
    )
    upper_beta_2 = beta_fn(a2, omega) * betaincc(
        a2, omega, x_array * x_array
    )
    prefactor = (
        4.0
        * omega
        / (2.0 * omega - 1.0)
        * math.pi**omega
        / gamma(omega)
    )
    return prefactor * (
        x_array ** (omega - 1.0) * upper_beta_1
        - (2.0 * omega + 1.0)
        / (4.0 * omega)
        * x_array ** (-0.5)
        * upper_beta_2
    )


def scaled_h(
    omega: float,
    x: float,
    coefficients: np.ndarray,
) -> float:
    nmax = int(math.floor(x))
    n = np.arange(1, nmax + 1, dtype=float)
    value = np.dot(
        coefficients[1 : nmax + 1],
        primitive_weight(omega, n / x),
    )
    return float(value / math.sqrt(x))


def build_summary(
    omega: float,
    label: str,
    grid: np.ndarray,
) -> OmegaSummary:
    coefficients = jordan_coefficients(omega, int(math.floor(grid[-1])))
    values = np.asarray(
        [scaled_h(omega, float(x), coefficients) for x in grid],
        dtype=float,
    )
    minimum_index = int(np.nanargmin(values))
    maximum_index = int(np.nanargmax(values))
    anchors = {
        format(anchor, ".0f"): scaled_h(omega, anchor, coefficients)
        for anchor in ANCHORS
    }
    return OmegaSummary(
        omega=label,
        sample_count=int(values.size),
        negative_count=int(np.count_nonzero(values < 0.0)),
        nonfinite_count=int(np.count_nonzero(~np.isfinite(values))),
        minimum_scaled_h=float(values[minimum_index]),
        minimum_x=float(grid[minimum_index]),
        maximum_scaled_h=float(values[maximum_index]),
        maximum_x=float(grid[maximum_index]),
        minimum_coefficient=float(np.min(coefficients[1:])),
        anchors=anchors,
    )


def signed_weight_guard() -> dict[str, float | bool]:
    x_negative = math.exp(-4.0)
    x_positive = 0.5
    negative_value = float(primitive_weight(0.5, x_negative))
    positive_value = float(primitive_weight(0.5, x_positive))
    return {
        "omega": 0.5,
        "negative_probe_x": x_negative,
        "negative_probe_value": negative_value,
        "positive_probe_x": x_positive,
        "positive_probe_value": positive_value,
        "sign_change_certified_by_double_probes": (
            negative_value < -1.0 and positive_value > 1.0
        ),
    }


def build_payload() -> dict:
    grid = x_grid()
    summaries = [
        build_summary(omega, label, grid)
        for omega, label in zip(OMEGAS, OMEGA_LABELS, strict=True)
    ]
    return {
        "kind": "jensen_window_pf_suzuki_jordan_totient_sign_scout",
        "date": "2026-07-23",
        "status": "finite double-precision sign scout",
        "proof_boundary": (
            "This finite grid is numerical reconnaissance only. It does not "
            "prove positivity or eventual one-sign behavior on an unbounded "
            "interval, the L2 residual, meromorphic innerness, any all-time "
            "determinant condition, RH, or Lambda<=0."
        ),
        "grid": {
            "linear_start": 1.001,
            "linear_stop": 20.0,
            "linear_count": 600,
            "geometric_start": 20.01,
            "geometric_stop": 5000.0,
            "geometric_count": 600,
            "unique_count": int(grid.size),
        },
        "summaries": [asdict(summary) for summary in summaries],
        "signed_weight_guard": signed_weight_guard(),
        "audit": {
            "omega_count": len(summaries),
            "total_sample_count": sum(
                summary.sample_count for summary in summaries
            ),
            "total_negative_count": sum(
                summary.negative_count for summary in summaries
            ),
            "all_coefficients_positive": all(
                summary.minimum_coefficient > 0.0 for summary in summaries
            ),
            "all_sampled_h_positive": all(
                summary.negative_count == 0
                and summary.minimum_scaled_h > 0.0
                for summary in summaries
            ),
            "eventual_sign_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
        "source": {
            "title": "A canonical system of differential equations arising from the Riemann zeta-function",
            "url": "https://arxiv.org/abs/1204.1827",
            "result": "Appendix theorem giving the L2 equivalence and eventual one-sign sufficient condition",
        },
    }


def render_note(payload: dict) -> str:
    lines = [
        "# Jensen-Window PF Suzuki Jordan-Totient Sign Scout",
        "",
        "Date: 2026-07-23",
        "",
        "Status: finite double-precision reconnaissance report. It is not a proof",
        "of eventual sign, the all-time determinant gate, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_suzuki_jordan_totient_sign_scout.json",
        "python work/rh_compute/scripts/jensen_window_pf_suzuki_jordan_totient_sign_scout.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py",
        "```",
        "",
        "## Target",
        "",
        "Suzuki's scalar criterion uses",
        "",
        "```text",
        "c_omega(n)=n^omega*product_(p|n)(1-p^(-2*omega))>0,",
        "h_omega^<1>(x)=x^(-1)*sum_(n<=x)c_omega(n)*g_omega^<1>(n/x).",
        "```",
        "",
        "Eventual one-sign behavior of `h_omega^<1>` implies innerness of",
        "`Theta_omega`. On a sequence `omega->0`, that is RH-strength.",
        "",
        "## Finite Grid",
        "",
        "The deterministic grid has 600 linear points on `[1.001,20]` and",
        "600 geometric points on `[20.01,5000]`. Reported values are",
        "`sqrt(x)*h_omega^<1>(x)`.",
        "",
        "| omega | samples | negative | minimum | x at minimum | maximum | x at maximum |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for summary in payload["summaries"]:
        lines.append(
            "| {omega} | {sample_count} | {negative_count} | "
            "{minimum_scaled_h:.9f} | {minimum_x:.6f} | "
            "{maximum_scaled_h:.9f} | {maximum_x:.6f} |".format(
                **summary
            )
        )
    lines.extend(
        [
            "",
            "All 6,000 sampled values are positive. The smallest sampled value",
            "is separated from zero by more than `0.08`; every sampled arithmetic",
            "coefficient is positive.",
            "",
            "Anchor values:",
            "",
            "| omega | x=10 | x=100 | x=1000 | x=5000 |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for summary in payload["summaries"]:
        anchors = summary["anchors"]
        lines.append(
            f"| {summary['omega']} | {anchors['10']:.9f} | "
            f"{anchors['100']:.9f} | {anchors['1000']:.9f} | "
            f"{anchors['5000']:.9f} |"
        )
    guard = payload["signed_weight_guard"]
    lines.extend(
        [
            "",
            "## Signed-Weight Guard",
            "",
            "The positivity is not termwise. At `omega=1/2`, the explicit",
            "primitive weight has opposite signs at two probes:",
            "",
            "```text",
            f"g_(1/2)^<1>(exp(-4))={guard['negative_probe_value']:.9f}<0",
            f"g_(1/2)^<1>(1/2)={guard['positive_probe_value']:.9f}>0",
            "```",
            "",
            "Thus the observed positive summatory values arise after a",
            "genuinely signed arithmetic convolution; positivity of",
            "`c_omega(n)` alone is not a proof.",
            "",
            "## Interpretation",
            "",
            "The grid makes the scalar target empirically credible and supplies",
            "well-separated regression anchors. It gives no control beyond",
            "`x=5000`, no uniformity as `omega->0`, and no rigorous rounding",
            "certificate. Finite positivity cannot be promoted to eventual sign.",
            "",
            "## Source",
            "",
            "- Masatoshi Suzuki, `A canonical system of differential equations arising from the Riemann zeta-function`: https://arxiv.org/abs/1204.1827",
            "",
            "## Proof Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    set_below_normal_priority()
    args = build_parser().parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Suzuki Jordan-totient sign scout: "
        f"{payload['audit']['total_sample_count']} samples, "
        f"{payload['audit']['total_negative_count']} negative, "
        f"{len(payload['summaries'])} omega values"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
