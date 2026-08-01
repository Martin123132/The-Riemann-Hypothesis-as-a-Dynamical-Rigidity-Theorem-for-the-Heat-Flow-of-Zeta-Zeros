#!/usr/bin/env python3
"""Scout the explicit Suzuki k=2 central-polynomial residual."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
import json
import math
from pathlib import Path
import time

import mpmath as mp
import numpy as np
import psutil
from scipy.integrate import cumulative_trapezoid
from scipy.special import beta as beta_fn
from scipy.special import betaincc, gamma

try:
    from threadpoolctl import threadpool_limits
except ImportError:
    threadpool_limits = None


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_k2_residual_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_suzuki_k2_residual_scout.md"
)
OMEGAS = (
    ("1/4", 0.25),
    ("1/8", 0.125),
    ("1/16", 0.0625),
    ("1/32", 0.03125),
)
ANCHORS = (10.0, 100.0, 1000.0, 5000.0)


def set_below_normal_priority() -> None:
    try:
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except (AttributeError, psutil.Error):
        pass


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


def primitive_weight(omega: float, x: np.ndarray) -> np.ndarray:
    a1 = (3.0 - 2.0 * omega) / 2.0
    a2 = (5.0 - 2.0 * omega) / 4.0
    upper_beta_1 = beta_fn(a1, omega) * betaincc(
        a1, omega, x * x
    )
    upper_beta_2 = beta_fn(a2, omega) * betaincc(
        a2, omega, x * x
    )
    prefactor = (
        4.0
        * omega
        / (2.0 * omega - 1.0)
        * math.pi**omega
        / gamma(omega)
    )
    return prefactor * (
        x ** (omega - 1.0) * upper_beta_1
        - (2.0 * omega + 1.0)
        / (4.0 * omega)
        * x ** (-0.5)
        * upper_beta_2
    )


def h1_value(
    omega: float,
    x: float,
    coefficients: np.ndarray,
) -> float:
    nmax = int(math.floor(x))
    if nmax < 1:
        return 0.0
    n = np.arange(1, nmax + 1, dtype=float)
    value = np.dot(
        coefficients[1 : nmax + 1],
        primitive_weight(omega, n / x),
    )
    return float(value / math.sqrt(x))


def xi(value: mp.mpf | mp.mpc) -> mp.mpf | mp.mpc:
    return (
        value
        * (value - 1)
        * mp.power(mp.pi, -value / 2)
        * mp.gamma(value / 2)
        * mp.zeta(value)
    )


def q1_value(omega: float) -> float:
    mp.mp.dps = 60
    point = mp.mpf("0.5") + mp.mpf(str(omega))
    return float(-2 * mp.diff(lambda s: mp.log(xi(s)), point))


def build_g2_kernel(
    omega: float,
    t_max: float,
    fine_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # g_(omega,1)(exp(-v))*exp(-v/2) behaves like v^omega at
    # v=0. This power map resolves that endpoint without a huge uniform grid.
    power = 4.0 / (1.0 + omega)
    unit_grid = np.linspace(0.0, 1.0, fine_count)
    v_fine = t_max * unit_grid**power
    y_fine = np.exp(-v_fine)
    integrand_fine = np.exp(-v_fine / 2.0) * primitive_weight(
        omega, y_fine
    )
    integral_fine = cumulative_trapezoid(
        integrand_fine, v_fine, initial=0.0
    )
    v_coarse = v_fine[::2]
    integral_coarse = cumulative_trapezoid(
        integrand_fine[::2], v_coarse, initial=0.0
    )
    return v_fine, integral_fine, v_coarse, integral_coarse


def h2_value(
    x: float,
    weighted_coefficients: np.ndarray,
    v_grid: np.ndarray,
    g2_integral: np.ndarray,
) -> float:
    nmax = int(math.floor(x))
    if nmax < 1:
        return 0.0
    n = np.arange(1, nmax + 1, dtype=float)
    v = np.log(x / n)
    kernel = np.interp(v, v_grid, g2_integral)
    return float(np.dot(weighted_coefficients[1 : nmax + 1], kernel))


def write_checkpoint(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def build_summary(
    label: str,
    omega: float,
    t_grid: np.ndarray,
    kernel_fine_count: int,
    cpu_threshold: float,
    cpu_samples: list[float],
) -> tuple[dict, bool]:
    coefficients = jordan_coefficients(
        omega, int(math.floor(math.exp(float(t_grid[-1]))))
    )
    denominator = np.sqrt(np.arange(len(coefficients), dtype=float))
    denominator[0] = 1.0
    weighted_coefficients = coefficients / denominator
    weighted_coefficients[0] = 0.0
    (
        v_fine,
        g2_fine,
        v_coarse,
        g2_coarse,
    ) = build_g2_kernel(
        omega, float(t_grid[-1]), kernel_fine_count
    )
    h2_values = np.empty_like(t_grid)
    high_cpu_streak = 0
    park_requested = False
    for index, t_value in enumerate(t_grid):
        h2_values[index] = h2_value(
            math.exp(float(t_value)),
            weighted_coefficients,
            v_fine,
            g2_fine,
        )
        if index and index % 200 == 0:
            sample = float(psutil.cpu_percent(interval=0.2))
            cpu_samples.append(sample)
            if sample > cpu_threshold:
                high_cpu_streak += 1
            else:
                high_cpu_streak = 0
            if high_cpu_streak >= 2:
                park_requested = True

    q1 = q1_value(omega)
    residual = h2_values - t_grid - q1
    coarse_t = t_grid[::2]
    coarse_residual = residual[::2]
    fine_energy = float(np.trapezoid(residual * residual, t_grid))
    coarse_energy = float(
        np.trapezoid(coarse_residual * coarse_residual, coarse_t)
    )
    tail_start = 3 * (len(t_grid) - 1) // 4
    tail = residual[tail_start:]
    anchor_residuals = {
        format(anchor, ".0f"): float(
            np.interp(math.log(anchor), t_grid, residual)
        )
        for anchor in ANCHORS
    }
    anchor_h1 = {
        format(anchor, ".0f"): h1_value(
            omega, anchor, coefficients
        )
        for anchor in ANCHORS
    }
    h2_anchor_kernel_differences = []
    for anchor in ANCHORS:
        fine_anchor = h2_value(
            anchor, weighted_coefficients, v_fine, g2_fine
        )
        coarse_anchor = h2_value(
            anchor, weighted_coefficients, v_coarse, g2_coarse
        )
        h2_anchor_kernel_differences.append(
            abs(fine_anchor - coarse_anchor)
        )
    cumulative_discrepancy = h2_values - t_grid
    return (
        {
            "omega": label,
            "omega_float": omega,
            "q1": q1,
            "initial_residual": float(residual[0]),
            "final_residual": float(residual[-1]),
            "minimum_residual": float(np.min(residual)),
            "minimum_t": float(t_grid[int(np.argmin(residual))]),
            "maximum_residual": float(np.max(residual)),
            "maximum_t": float(t_grid[int(np.argmax(residual))]),
            "fine_finite_energy": fine_energy,
            "coarse_finite_energy": coarse_energy,
            "coarse_fine_energy_relative_difference": abs(
                fine_energy - coarse_energy
            )
            / max(np.finfo(float).tiny, abs(fine_energy)),
            "kernel_fine_count": int(len(v_fine)),
            "kernel_coarse_count": int(len(v_coarse)),
            "kernel_resolution_anchor_max_abs_difference": float(
                max(h2_anchor_kernel_differences)
            ),
            "tail_quarter_rms": float(
                math.sqrt(np.mean(tail * tail))
            ),
            "tail_quarter_max_abs": float(np.max(np.abs(tail))),
            "tail_quarter_mean": float(np.mean(tail)),
            "final_cumulative_h1_minus_one": float(
                cumulative_discrepancy[-1]
            ),
            "recurrence_identity_error": abs(
                float(residual[-1])
                - (float(cumulative_discrepancy[-1]) - q1)
            ),
            "anchor_residuals": anchor_residuals,
            "anchor_h1": anchor_h1,
        },
        park_requested,
    )


def render_note(payload: dict) -> str:
    lines = [
        "# Jensen-Window PF Suzuki K2 Residual Scout",
        "",
        "Date: 2026-07-23",
        "",
        "Status: finite double-precision reconnaissance. This is not a proof",
        "of an L2 tail, fixed-shift innerness, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_suzuki_k2_residual_scout.json",
        "python work/rh_compute/scripts/jensen_window_pf_suzuki_k2_residual_scout.py",
        "python work/rh_compute/scripts/jensen_window_pf_suzuki_k2_residual_scout.py --resume",
        "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_k2_residual_scout.py",
        "```",
        "",
        "## Target",
        "",
        "For `t=log x`, the exact level-two residual is",
        "",
        "```text",
        "r_(omega,2)(t)",
        " =-q_(omega,1)",
        "  +integral_0^t[H_(omega,1)(exp(u))-1]du,",
        "q_(omega,1)=-2*xi'(1/2+omega)/xi(1/2+omega).",
        "```",
        "",
        "A cofinal family of genuine `L2(0,infinity)` estimates would be",
        "RH-equivalent. This scout integrates only through `x=5000`.",
        "",
        "## Finite Results",
        "",
        "| omega | q1 | r2(10) | r2(100) | r2(1000) | r2(5000) | finite energy | coarse/fine rel. |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in payload["summaries"]:
        anchors = item["anchor_residuals"]
        lines.append(
            f"| {item['omega']} | {item['q1']:.8g} | "
            f"{anchors['10']:.8g} | {anchors['100']:.8g} | "
            f"{anchors['1000']:.8g} | {anchors['5000']:.8g} | "
            f"{item['fine_finite_energy']:.8g} | "
            f"{item['coarse_fine_energy_relative_difference']:.3g} |"
        )
    lines.extend(
        [
            "",
            "The fine grid is uniform in `t` with "
            f"`{payload['grid']['fine_count']}` points; its even-index",
            "subgrid provides the independent coarse quadrature comparison.",
            "",
            "## Interpretation",
            "",
            "These values test normalization and reveal finite-range shape.",
            "They cannot establish convergence, square integrability,",
            "positive-time Hardy support, uniformity as `omega->0`, or any",
            "unbounded arithmetic estimate. Boundary spectral energy is",
            "automatic and is not what this time-domain scout certifies.",
            "",
            "## Runtime",
            "",
            f"Mode: `{payload['runtime']['resource_mode']}`; one below-normal",
            "worker, numerical thread pools limited to one.",
            f"Completed shifts: `{len(payload['summaries'])}/4`; parked:",
            f"`{str(payload['runtime']['parked']).lower()}`.",
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
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--x-max", type=float, default=5000.0)
    parser.add_argument("--fine-count", type=int, default=2401)
    parser.add_argument("--kernel-fine-count", type=int, default=68001)
    parser.add_argument("--cpu-threshold", type=float, default=75.0)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="validate and continue the checkpoint at --out",
    )
    return parser


def expected_grid(args: argparse.Namespace) -> dict:
    return {
        "x_min": 1.0,
        "x_max": args.x_max,
        "t_max": math.log(args.x_max),
        "fine_count": args.fine_count,
        "coarse_count": (args.fine_count + 1) // 2,
        "kernel_fine_count": args.kernel_fine_count,
        "kernel_coarse_count": (args.kernel_fine_count + 1) // 2,
        "uniform_in": "t=log(x)",
    }


def new_payload(args: argparse.Namespace, cpu_samples: list[float]) -> dict:
    return {
        "kind": "jensen_window_pf_suzuki_k2_residual_scout",
        "date": "2026-07-23",
        "status": "running finite double-precision reconnaissance",
        "proof_boundary": (
            "Finite double-precision time-domain quadrature only. It does "
            "not prove convergence as x tends to infinity, an L2 residual, "
            "Hardy causality, fixed-shift innerness, RH, or Lambda<=0."
        ),
        "grid": expected_grid(args),
        "summaries": [],
        "runtime": {
            "resource_mode": "daytime_one_worker",
            "priority": "below_normal",
            "thread_limit": 1,
            "cpu_threshold": args.cpu_threshold,
            "cpu_samples": cpu_samples,
            "parked": False,
            "completed": False,
        },
        "audit": {
            "l2_tail_proved": False,
            "hardy_causality_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def load_resume_payload(
    args: argparse.Namespace,
) -> tuple[dict, list[float], float]:
    if not args.out.is_file():
        raise ValueError(f"resume checkpoint does not exist: {args.out}")
    payload = json.loads(args.out.read_text(encoding="utf-8"))
    if payload.get("kind") != "jensen_window_pf_suzuki_k2_residual_scout":
        raise ValueError("resume checkpoint has the wrong kind")

    actual_grid = payload.get("grid", {})
    for key, expected in expected_grid(args).items():
        actual = actual_grid.get(key)
        if isinstance(expected, float):
            matches = isinstance(actual, (int, float)) and math.isclose(
                float(actual), expected, rel_tol=0.0, abs_tol=1e-14
            )
        else:
            matches = actual == expected
        if not matches:
            raise ValueError(
                f"resume grid mismatch for {key}: "
                f"checkpoint={actual!r}, requested={expected!r}"
            )

    summaries = payload.get("summaries")
    if not isinstance(summaries, list) or len(summaries) > len(OMEGAS):
        raise ValueError("resume checkpoint has an invalid summary list")
    expected_prefix = tuple(label for label, _ in OMEGAS[: len(summaries)])
    actual_prefix = tuple(item.get("omega") for item in summaries)
    if actual_prefix != expected_prefix:
        raise ValueError(
            "resume checkpoint omega order mismatch: "
            f"{actual_prefix!r} != {expected_prefix!r}"
        )

    runtime = payload.setdefault("runtime", {})
    cpu_samples = runtime.setdefault("cpu_samples", [])
    if not isinstance(cpu_samples, list):
        raise ValueError("resume checkpoint CPU samples are not a list")
    previous_elapsed = float(
        runtime.get(
            "elapsed_seconds_total",
            runtime.get("elapsed_seconds", 0.0),
        )
    )
    runtime.update(
        {
            "resource_mode": "daytime_one_worker",
            "priority": "below_normal",
            "thread_limit": 1,
            "cpu_threshold": args.cpu_threshold,
            "parked": False,
            "completed": False,
            "resume_count": int(runtime.get("resume_count", 0)) + 1,
        }
    )
    payload["status"] = "running finite double-precision reconnaissance"
    return payload, cpu_samples, previous_elapsed


def main() -> int:
    args = build_parser().parse_args()
    set_below_normal_priority()
    started = time.monotonic()
    if args.resume:
        try:
            payload, cpu_samples, previous_elapsed = load_resume_payload(
                args
            )
        except (OSError, ValueError, json.JSONDecodeError) as error:
            raise SystemExit(f"cannot resume scout: {error}") from error
    else:
        cpu_samples = []
        previous_elapsed = 0.0
        payload = new_payload(args, cpu_samples)

    t_grid = np.linspace(
        0.0, math.log(args.x_max), args.fine_count
    )
    limiter = (
        threadpool_limits(limits=1)
        if threadpool_limits is not None
        else nullcontext()
    )
    with limiter:
        for label, omega in OMEGAS[len(payload["summaries"]) :]:
            summary, park_requested = build_summary(
                label,
                omega,
                t_grid,
                args.kernel_fine_count,
                args.cpu_threshold,
                cpu_samples,
            )
            payload["summaries"].append(summary)
            cycle_elapsed = time.monotonic() - started
            payload["runtime"]["cycle_elapsed_seconds"] = cycle_elapsed
            payload["runtime"]["elapsed_seconds"] = (
                previous_elapsed + cycle_elapsed
            )
            payload["runtime"]["elapsed_seconds_total"] = (
                previous_elapsed + cycle_elapsed
            )
            write_checkpoint(args.out, payload)
            if park_requested:
                payload["runtime"]["parked"] = True
                break

    payload["runtime"]["completed"] = (
        len(payload["summaries"]) == len(OMEGAS)
    )
    if payload["runtime"]["completed"]:
        payload["runtime"]["parked"] = False
    payload["status"] = (
        "finite double-precision reconnaissance complete"
        if payload["runtime"]["completed"]
        else "finite double-precision reconnaissance parked"
    )
    cycle_elapsed = time.monotonic() - started
    payload["runtime"]["cycle_elapsed_seconds"] = cycle_elapsed
    payload["runtime"]["elapsed_seconds"] = previous_elapsed + cycle_elapsed
    payload["runtime"]["elapsed_seconds_total"] = (
        previous_elapsed + cycle_elapsed
    )
    write_checkpoint(args.out, payload)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Suzuki k=2 residual scout: "
        f"{len(payload['summaries'])}/{len(OMEGAS)} shifts, "
        f"parked={payload['runtime']['parked']}, "
        f"elapsed={payload['runtime']['elapsed_seconds']:.1f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
