#!/usr/bin/env python3
"""Absorb every theta summand n>=6 into the local degree-zero theorem."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import jensen_window_pf_newman_theta_fifth_summand_positive_time_boundary_scout as scout  # noqa: E402


continuation = scout.continuation
STEM = "jensen_window_pf_newman_theta_complete_tail_positive_time_degree_certificate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
FIFTH_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_fifth_summand_"
    "positive_time_degree_certificate.json"
)
INDEX_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_first_jet_winding_gate.json"
)

DATE = "2026-08-04"
PRECISION_BITS = 288
ABSOLUTE_TOLERANCE = "1e-44"
TIME_HIGH = Fraction(9, 20)
TAIL_RATIO_BOUND = "1e-9"
MOMENT_SAFETY_FACTOR = "2"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_component_ratio_audit() -> dict:
    require(7**4 < 2 * 6**4, "successive polynomial ratio bound failed")
    require(6 < 2**38, "geometric denominator bound failed")
    require(6 * 10**9 < 2**33, "decimal tail-ratio conversion failed")
    return {
        "theta_pi_provenance": (
            "The pi is the canonical Jacobi-theta normalization in phi_n, inherited "
            "from sum_n exp(-pi*n^2*y); no arbitrary circle constant is inserted."
        ),
        "pointwise_reduction": (
            "For z_n=pi*n^2*exp(4u), phi_n=exp(u)*z_n*(2z_n-3)*exp(-z_n). "
            "Since z_5>25*3>3, 2z_5-3>=z_5 and 2z_n-3<=2z_n. Thus for n>=6, "
            "phi_n/phi_5 <= 2*(n/5)^4*exp(-pi*(n^2-25)*exp(4u))."
        ),
        "uniform_exponential_reduction": (
            "Using pi>3 and exp(4u)>=1 on u>=0 gives "
            "phi_n/phi_5 < b_n=2*(n/5)^4*exp(-3*(n^2-25))."
        ),
        "geometric_tail": (
            "For n>=6, b_(n+1)/b_n <= (7/6)^4*exp(-3*(2n+1)) "
            "<2*2^-39=2^-38, using e>2. Also b_6<5*2^-33. Hence "
            "sum_(n>=6)b_n < (5*2^-33)/(1-2^-38) < 6*2^-33 < 10^-9."
        ),
        "conclusion": (
            "For every u>=0, sum_(n>=6) phi_n(u) < 10^-9*phi_5(u)."
        ),
        "tail_ratio_upper": TAIL_RATIO_BOUND,
    }


def render_note(payload: dict) -> str:
    lines = [
        "# Complete-Theta Local Positive-Time Degree Certificate",
        "",
        f"Date: {DATE}",
        "",
        "Status: rigorous complete-theta local contact exclusion. This is not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Theorem",
        "",
        "The complete theta-kernel heat transform and its first x derivative have no",
        "common zero in",
        "",
        "```text",
        "0<=t<=0.45,  135.5<=x<=136.",
        "```",
        "",
        f"The certified finite boundary separation is `{payload['finite_margin_lower']}`.",
        f"The complete n>=6 tail budget is `{payload['maximum_tail_derivative_bound']}`.",
        f"Their ratio is at most `{payload['tail_to_margin_ratio_upper']}`.",
        "",
        "The boundary homotopy from the five-term endpoint to the complete theta sum",
        "therefore stays nonzero and preserves winding zero. The all-multiplicity",
        "same-sign index theorem then excludes every interior contact.",
        "",
        "## Proof Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build(args: argparse.Namespace, baseline: list[float], priority_lowered: bool) -> dict:
    started = perf_counter()
    ratio_audit = exact_component_ratio_audit()
    n5_tail_audit = scout.extended_tail_audit()
    fifth_payload = json.loads(FIFTH_PARENT.read_text(encoding="utf-8"))
    index_payload = json.loads(INDEX_PARENT.read_text(encoding="utf-8"))
    require(
        fifth_payload["status"]
        == "rigorous fifth-summand positive-time exclusion certificate complete",
        "fifth-summand parent is incomplete",
    )
    require(
        fifth_payload["mu_zero_winding"]["polygon"]["winding_t_x"] == 0,
        "fifth-summand parent winding is not zero",
    )
    require(
        "+floor(m/2)" in json.dumps(index_payload, sort_keys=True),
        "all-multiplicity index parent is incomplete",
    )

    continuation.flint.ctx.prec = args.precision_bits
    continuation.ABS_TOL = args.absolute_tolerance
    fifth = scout.FifthTransformFamily(
        Fraction(0),
        TIME_HIGH,
        TIME_HIGH,
        Fraction(0),
        weighted=False,
    )
    raw_moment_uppers = [fifth.moment_upper(order, TIME_HIGH) for order in range(2)]
    moment_safety_factor = continuation.arb(MOMENT_SAFETY_FACTOR)
    moment_uppers = [moment_safety_factor * value for value in raw_moment_uppers]
    tail_ratio = continuation.arb(TAIL_RATIO_BOUND)
    tail_bounds = [tail_ratio * moment for moment in moment_uppers]
    maximum_tail = max(tail_bounds)
    finite_margin = continuation.arb(
        fifth_payload["boundary_certificate"]["minimum_separation_lower"]
    ).lower()
    require(maximum_tail.upper() < finite_margin.lower(), "complete tail exceeds finite margin")
    tail_to_margin = maximum_tail / finite_margin
    require(tail_to_margin.upper() < 1, "tail-to-margin ratio is not below one")

    return {
        "kind": STEM,
        "schema_version": 1,
        "date": DATE,
        "status": "rigorous complete-theta local positive-time exclusion complete",
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "parents": {
            "fifth_summand_degree": {
                "path": str(FIFTH_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(FIFTH_PARENT),
            },
            "all_multiplicity_index": {
                "path": str(INDEX_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(INDEX_PARENT),
            },
        },
        "resource_policy": {
            "mode": "daytime",
            "active_compute_workers": 1,
            "thread_caps": 1,
            "below_normal_priority_applied": priority_lowered,
            "baseline_cpu_percent": baseline,
            "baseline_mean_percent": sum(baseline) / len(baseline),
            "elapsed_seconds": perf_counter() - started,
        },
        "rectangle": {
            "time": ["0", "0.45"],
            "x": ["135.5", "136"],
            "coordinate_orientation": "(t,x)",
        },
        "component_ratio_audit": ratio_audit,
        "n5_positive_moment_tail_audit": n5_tail_audit,
        "positive_moment_budget": {
            "reference_component": "phi_5",
            "time_upper": "0.45",
            "derivative_orders": [0, 1],
            "precision_bits": args.precision_bits,
            "absolute_integration_tolerance": args.absolute_tolerance,
            "safety_factor": MOMENT_SAFETY_FACTOR,
            "raw_moment_uppers": [
                continuation.arb_text(value) for value in raw_moment_uppers
            ],
            "moment_uppers": [continuation.arb_text(value) for value in moment_uppers],
            "tail_derivative_bounds": [
                continuation.arb_text(value) for value in tail_bounds
            ],
            "uniformity": (
                "For 0<=t<=0.45 and every real x, |partial_x^j H_n(t,x)| is "
                "bounded by the positive j-th moment at t=0.45."
            ),
        },
        "finite_margin_lower": continuation.arb_text(finite_margin),
        "maximum_tail_derivative_bound": continuation.arb_text(maximum_tail),
        "tail_to_margin_ratio_upper": continuation.arb_text(tail_to_margin.upper()),
        "boundary_homotopy": (
            "For F_s=sum_(n=1)^5 H_n+s*sum_(n>=6)H_n, 0<=s<=1, each "
            "selected boundary component retains at least the five-term certified "
            "margin minus the complete-tail bound. This is strictly positive."
        ),
        "degree_composition": (
            "The five-term endpoint has boundary winding zero. The complete-tail "
            "homotopy preserves it. Every multiplicity-m contact has local index "
            "-floor(m/2) in (t,x), so total degree zero excludes every contact."
        ),
        "theorem": (
            "The complete theta-kernel heat transform and its first x derivative "
            "have no common zero in [0,0.45]x[135.5,136]."
        ),
        "proof_boundary": (
            "This is a rigorous complete-theta contact exclusion only in the declared "
            "local positive-time rectangle. It does not control other frequency "
            "rectangles, negative heat time, all real x, the global de Bruijn-Newman "
            "constant, a degree-uniform Jensen remainder, Lambda<=0, RH, or a "
            "prize-level conclusion."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--precision-bits", type=int, default=PRECISION_BITS)
    parser.add_argument("--absolute-tolerance", default=ABSOLUTE_TOLERANCE)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    priority_lowered = continuation.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "complete-theta tail baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("complete-theta tail run deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "could not apply below-normal process priority")
    payload = build(args, baseline, priority_lowered)
    continuation.write_json_atomic(args.out, payload)
    continuation.write_text_atomic(args.note, render_note(payload))
    print(
        "built complete-theta local degree certificate: "
        f"tail_to_margin={payload['tail_to_margin_ratio_upper']}, "
        "winding=0, no_contacts=True"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
