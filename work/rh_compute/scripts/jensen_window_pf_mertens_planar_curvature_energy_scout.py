#!/usr/bin/env python3
"""Build a finite scaling scout for the planar curvature-energy gate."""

from __future__ import annotations

import argparse
import ctypes
import json
import math
import os
from pathlib import Path
from time import perf_counter


for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_curvature_energy_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_curvature_energy_scout.md"
)
DEFAULT_ALPHAS = (0.125, 0.25, 0.5, 0.75)


def set_below_normal_priority() -> None:
    if os.name != "nt":
        return
    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = (
            ctypes.c_void_p,
            ctypes.c_uint32,
        )
        kernel32.SetPriorityClass.restype = ctypes.c_int
        kernel32.SetPriorityClass(
            kernel32.GetCurrentProcess(),
            0x00004000,
        )
    except (AttributeError, OSError):
        pass


def mobius_values(limit: int) -> np.ndarray:
    values = np.ones(limit + 1, dtype=np.int8)
    values[0] = 0
    prime = np.ones(limit + 1, dtype=bool)
    prime[:2] = False
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        prime[p::p] = False
        values[p::p] *= -1
        values[p * p :: p * p] = 0
    return values


def dyadic_sizes(max_size: int) -> list[int]:
    if max_size < 16:
        raise ValueError("max_size must be at least 16")
    sizes: list[int] = []
    size = 16
    while size <= max_size:
        sizes.append(size)
        size *= 2
    return sizes


def arithmetic_arrays(
    size: int,
    mobius: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    side = 3 * size - 2
    arithmetic = np.zeros((side, side), dtype=np.int8)
    shifts = np.arange(1, side + 1)
    for row, n in enumerate(range(size + 1, 4 * size - 1)):
        m = n + shifts
        valid = m <= 4 * size - 1
        arithmetic[row, valid] = mobius[n] * mobius[m[valid]]
    prefixes = arithmetic.cumsum(
        axis=0,
        dtype=np.int64,
    ).cumsum(
        axis=1,
        dtype=np.int64,
    )
    return arithmetic, prefixes


def kernel_surface(
    size: int,
    alpha: float,
) -> tuple[int, np.ndarray, np.ndarray]:
    ambient = 2 * size + 1
    rank = min(
        size,
        math.ceil(size ** (1.0 - alpha / 2.0)),
    )
    odd = np.arange(1, 2 * rank, 2, dtype=np.float64)
    grid = (
        np.arange(0, 2 * size + 1, dtype=np.float64)
        * math.pi
        / ambient
    )
    cosine = np.cos(np.outer(grid, odd)) @ (1.0 / (odd * odd))
    indices = np.arange(1, size + 1)
    current_kernel = 2.0 / ambient * (
        cosine[np.abs(indices[:, None] - indices[None, :])]
        - cosine[indices[:, None] + indices[None, :]]
    )
    endpoint = current_kernel[:, -1]
    corner = current_kernel[-1, -1]
    future_indices = np.arange(
        2 * size,
        4 * size,
        dtype=np.float64,
    )
    future = (
        (2.0 * size / future_indices) ** (1.0 + alpha)
        * (4.0 * size - future_indices)
        / (2.0 * size)
    )
    side = 3 * size - 2
    extended = np.zeros((side + 1, side + 1), dtype=np.float64)
    shifts = np.arange(1, side + 1)
    for row, n in enumerate(range(size + 1, 4 * size - 1)):
        m = n + shifts
        valid = m <= 4 * size - 1
        if n < 2 * size:
            current = valid & (m < 2 * size)
            future_mask = valid & (m >= 2 * size)
            extended[row, :side][current] = current_kernel[
                n - size - 1,
                m[current] - size - 1,
            ]
            extended[row, :side][future_mask] = (
                endpoint[n - size - 1]
                * future[m[future_mask] - 2 * size]
            )
        else:
            extended[row, :side][valid] = (
                corner
                * future[n - 2 * size]
                * future[m[valid] - 2 * size]
            )
    mixed = (
        extended[:-1, :-1]
        - extended[1:, :-1]
        - extended[:-1, 1:]
        + extended[1:, 1:]
    )
    return rank, extended[:-1, :-1], mixed


def compute_row(
    size: int,
    alpha: float,
    arithmetic: np.ndarray,
    prefixes: np.ndarray,
) -> dict:
    started = perf_counter()
    rank, weights, mixed = kernel_surface(size, alpha)
    prefix_float = prefixes.astype(np.float64)
    absolute_mixed = np.abs(mixed)
    variation = float(absolute_mixed.sum())
    energy = float(
        np.sum(absolute_mixed * prefix_float * prefix_float)
    )
    direct = float(np.sum(arithmetic * weights))
    abel = float(np.sum(prefix_float * mixed))
    analytic_bound = (
        3.0
        * math.pi**2
        * (1.0 + math.log(2.0 * rank))
        / size
    )
    return {
        "K": size,
        "alpha": alpha,
        "R": rank,
        "mixed_variation": variation,
        "mixed_variation_times_K": variation * size,
        "analytic_variation_bound": analytic_bound,
        "variation_to_bound_ratio": variation / analytic_bound,
        "curvature_energy": energy,
        "energy_per_K": energy / size,
        "energy_per_K_log_2K": (
            energy / (size * math.log(2.0 * size))
        ),
        "max_planar_prefix": int(np.max(np.abs(prefixes))),
        "max_prefix_per_K": (
            int(np.max(np.abs(prefixes))) / size
        ),
        "signed_offdiagonal": direct,
        "double_abel_value": abel,
        "direct_abel_abs_error": abs(direct - abel),
        "elapsed_seconds": perf_counter() - started,
    }


def add_doubling_exponents(rows: list[dict]) -> None:
    previous: dict[float, dict] = {}
    for row in rows:
        alpha = float(row["alpha"])
        earlier = previous.get(alpha)
        if earlier is None:
            row["energy_doubling_exponent"] = None
            row["variation_doubling_exponent"] = None
        else:
            row["energy_doubling_exponent"] = math.log(
                row["curvature_energy"] / earlier["curvature_energy"],
                2.0,
            )
            row["variation_doubling_exponent"] = math.log(
                row["mixed_variation"] / earlier["mixed_variation"],
                2.0,
            )
        previous[alpha] = row


def build_payload(
    max_size: int = 1024,
    alphas: tuple[float, ...] = DEFAULT_ALPHAS,
) -> dict:
    set_below_normal_priority()
    sizes = dyadic_sizes(max_size)
    mobius = mobius_values(4 * max(sizes) + 4)
    rows: list[dict] = []
    for size in sizes:
        arithmetic, prefixes = arithmetic_arrays(size, mobius)
        for alpha in alphas:
            rows.append(
                compute_row(
                    size,
                    alpha,
                    arithmetic,
                    prefixes,
                )
            )
    add_doubling_exponents(rows)
    return {
        "kind": "jensen_window_pf_mertens_planar_curvature_energy_scout",
        "date": "2026-07-24",
        "status": (
            "finite double-precision scaling diagnostic for the open "
            "planar curvature-energy gate"
        ),
        "parameters": {
            "sizes": sizes,
            "alphas": list(alphas),
            "floating_point": "numpy float64",
            "worker_policy": (
                "one process, BLAS thread caps set to one, "
                "below-normal Windows priority"
            ),
        },
        "rows": rows,
        "audit": {
            "row_count": len(rows),
            "max_K": max(sizes),
            "alpha_count": len(alphas),
            "max_direct_abel_abs_error": max(
                row["direct_abel_abs_error"] for row in rows
            ),
            "max_variation_times_K": max(
                row["mixed_variation_times_K"] for row in rows
            ),
            "max_energy_per_K": max(
                row["energy_per_K"] for row in rows
            ),
            "max_prefix_per_K": max(
                row["max_prefix_per_K"] for row in rows
            ),
            "all_analytic_variation_bounds_pass": all(
                row["mixed_variation"]
                < row["analytic_variation_bound"]
                for row in rows
            ),
            "finite_energy_target_proved": False,
            "asymptotic_scaling_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
        "proof_boundary": (
            "Finite non-rigorous float64 diagnostics only. The grid does "
            "not prove an all-K estimate, a curvature-energy theorem, RH, "
            "PF-infinity, or Lambda <= 0."
        ),
    }


def render_note(payload: dict) -> str:
    lines = [
        "# Jensen-Window PF Mertens Planar Curvature-Energy Scout",
        "",
        "Date: 2026-07-24",
        "",
        "Status: finite double-precision scaling diagnostic for the open",
        "planar curvature-energy gate; not a proof artifact.",
        "This diagnostic does not prove an asymptotic theorem, RH,",
        "PF-infinity, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/"
        "jensen_window_pf_mertens_planar_curvature_energy_scout.json",
        "python work/rh_compute/scripts/"
        "jensen_window_pf_mertens_planar_curvature_energy_scout.py",
        "python work/rh_compute/scripts/"
        "check_jensen_window_pf_mertens_planar_curvature_energy_scout.py",
        "```",
        "",
        "## Purpose",
        "",
        "The exact handoff",
        "`outputs/jensen_window_pf_mertens_planar_abel_handoff.md` proves",
        "",
        "```text",
        "V_(alpha,K)<3*pi^2*(1+log(2R))/K,",
        "|O_(alpha,K)|^2<=V_(alpha,K)E_(alpha,K),",
        "```",
        "",
        "and leaves `E_(alpha,K)=O_epsilon(K^(1+epsilon))` open. This",
        "scout asks only whether the finite actual-Mobius values immediately",
        "falsify that scale.",
        "",
        "## Grid",
        "",
        "| K | alpha | R | K V_K | E_K/K | M_K/K | O_K | Abel error |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["rows"]:
        lines.append(
            "| {K} | {alpha:g} | {R} | {kv:.9g} | {ek:.9g} | "
            "{mk:.9g} | {off:.9g} | {error:.3g} |".format(
                K=row["K"],
                alpha=row["alpha"],
                R=row["R"],
                kv=row["mixed_variation_times_K"],
                ek=row["energy_per_K"],
                mk=row["max_prefix_per_K"],
                off=row["signed_offdiagonal"],
                error=row["direct_abel_abs_error"],
            )
        )
    audit = payload["audit"]
    lines.extend(
        [
            "",
            "## Finite Observation",
            "",
            "Across this recorded grid:",
            "",
            f"- maximum `K*V_K`: `{audit['max_variation_times_K']:.12g}`",
            f"- maximum `E_K/K`: `{audit['max_energy_per_K']:.12g}`",
            f"- maximum `M_K/K`: `{audit['max_prefix_per_K']:.12g}`",
            "- maximum direct/double-Abel discrepancy: "
            f"`{audit['max_direct_abel_abs_error']:.12g}`",
            "",
            "Thus the tested values through `K=1024` are compatible with",
            "`V_K=O(1/K)` and `E_K=O(K)` on this grid. This is a finite",
            "observation, not evidence for a uniform constant and not an",
            "extrapolation theorem. The proved analytic statement remains",
            "`V_K=O((1+log R)/K)`, and the all-scale Mobius energy estimate",
            "remains open.",
            "",
            "## Boundary",
            "",
            "The calculations use ordinary float64 arithmetic and finitely",
            "many dyadic blocks. Passing the checker confirms reproduction",
            "on selected small rows, the exact finite Abel identity to",
            "floating tolerance, and honest status language. It does not",
            "prove asymptotic scaling, the curvature-energy target, the full",
            "Burnol bound, RH, PF-infinity, or `Lambda <= 0`.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_alphas(value: str) -> tuple[float, ...]:
    alphas = tuple(float(item) for item in value.split(",") if item)
    if not alphas or any(not 0.0 < alpha < 1.0 for alpha in alphas):
        raise argparse.ArgumentTypeError(
            "alphas must be comma-separated values in (0,1)"
        )
    return alphas


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-k", type=int, default=1024)
    parser.add_argument(
        "--alphas",
        type=parse_alphas,
        default=DEFAULT_ALPHAS,
    )
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload(args.max_k, args.alphas)
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.result.write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built Mertens planar curvature-energy scout: "
        f"{len(payload['rows'])} rows, K<={payload['audit']['max_K']}, "
        f"{payload['audit']['alpha_count']} alpha values"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
