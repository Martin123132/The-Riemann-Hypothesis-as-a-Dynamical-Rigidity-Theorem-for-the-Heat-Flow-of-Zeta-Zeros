#!/usr/bin/env python3
"""Track a three-theta-summand tangency as the fourth coefficient turns on."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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

import mpmath as mp
import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_theta_fourth_summand_tangency_homotopy_scout"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
PARENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout.json"
)
SCHEMA_VERSION = 1
DATE = "2026-08-03"
DPS = 50
CUTOFF = mp.mpf(2)
LAMBDA_STEP_TEXT = "0.02"
LAMBDA_LAST_TEXT = "0.22"
T_MAX_TEXT = "0.5"


class ResourcePark(RuntimeError):
    """Raised only between completed continuation rows."""


class CpuMonitor:
    def __init__(self) -> None:
        self.samples: list[float] = []
        self.consecutive_high = 0

    def sample(self) -> None:
        value = float(psutil.cpu_percent(interval=0.5))
        self.samples.append(value)
        self.consecutive_high = self.consecutive_high + 1 if value > 75.0 else 0
        if self.consecutive_high >= 2:
            raise ResourcePark("two consecutive continuation CPU samples exceeded 75%")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def request_below_normal_priority() -> bool:
    try:
        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        else:
            process.nice(5)
        return True
    except (psutil.Error, PermissionError, OSError):
        return False


def decimal(value: mp.mpf, digits: int = 44) -> str:
    return mp.nstr(value, digits)


def integration_breaks(index: int) -> list[mp.mpf]:
    scale = 1 / (4 * mp.pi * index * index)
    candidates = (
        mp.mpf(0),
        scale,
        4 * scale,
        16 * scale,
        mp.mpf("0.02"),
        mp.mpf("0.05"),
        mp.mpf("0.1"),
        mp.mpf("0.2"),
        mp.mpf("0.5"),
        mp.mpf(1),
        CUTOFF,
    )
    result: list[mp.mpf] = []
    for value in sorted(set(candidates)):
        if 0 <= value <= CUTOFF and (not result or value > result[-1]):
            result.append(value)
    return result


def component_jets(index: int, time: mp.mpf, x: mp.mpf) -> list[mp.mpf]:
    n = mp.mpf(index)
    pi_n2 = mp.pi * n * n

    def kernel(u: mp.mpf) -> mp.mpf:
        exp4u = mp.exp(4 * u)
        return (
            pi_n2
            * mp.exp(5 * u + time * u * u)
            * (2 * pi_n2 * exp4u - 3)
            * mp.exp(-pi_n2 * exp4u)
        )

    derivatives = (
        lambda u: mp.cos(x * u),
        lambda u: -u * mp.sin(x * u),
        lambda u: -(u**2) * mp.cos(x * u),
        lambda u: (u**3) * mp.sin(x * u),
    )
    breaks = integration_breaks(index)
    return [
        mp.quad(lambda u, derivative=derivative: kernel(u) * derivative(u), breaks)
        for derivative in derivatives
    ]


def homotopy_field(
    coefficient: mp.mpf, time: mp.mpf, x: mp.mpf
) -> tuple[list[mp.mpf], list[mp.mpf]]:
    rows = [component_jets(index, time, x) for index in range(1, 5)]
    field = [
        sum(rows[index][order] for index in range(3))
        + coefficient * rows[3][order]
        for order in range(4)
    ]
    return field, rows[3]


def implicit_tangent(field: list[mp.mpf], fourth: list[mp.mpf]) -> mp.matrix:
    jacobian = mp.matrix(
        [[-field[2], field[1]], [-field[3], field[2]]]
    )
    return mp.lu_solve(jacobian, mp.matrix([-fourth[0], -fourth[1]]))


def correct_contact(
    coefficient: mp.mpf,
    time: mp.mpf,
    x: mp.mpf,
    maximum_iterations: int = 12,
) -> tuple[mp.mpf, mp.mpf, list[mp.mpf], list[mp.mpf], int]:
    tolerance = mp.mpf("1e-36")
    for iteration in range(maximum_iterations):
        field, fourth = homotopy_field(coefficient, time, x)
        jacobian = mp.matrix(
            [[-field[2], field[1]], [-field[3], field[2]]]
        )
        step = mp.lu_solve(jacobian, mp.matrix([-field[0], -field[1]]))
        size = max(abs(step[0]), abs(step[1]))
        if size > mp.mpf("0.25"):
            step *= mp.mpf("0.25") / size
        time += step[0]
        x += step[1]
        if size < mp.mpf("1e-30"):
            break
    field, fourth = homotopy_field(coefficient, time, x)
    require(
        max(abs(field[0]), abs(field[1])) < tolerance,
        f"contact Newton solve did not converge at lambda={coefficient}",
    )
    return time, x, field, fourth, iteration + 1


def solve_boundary_crossing(
    coefficient: mp.mpf, x: mp.mpf
) -> tuple[mp.mpf, mp.mpf, list[mp.mpf], list[mp.mpf], int]:
    tolerance = mp.mpf("1e-40")
    time = mp.mpf(0)
    for iteration in range(10):
        field, fourth = homotopy_field(coefficient, time, x)
        jacobian = mp.matrix(
            [[fourth[0], field[1]], [fourth[1], field[2]]]
        )
        step = mp.lu_solve(jacobian, mp.matrix([-field[0], -field[1]]))
        coefficient += step[0]
        x += step[1]
        if max(abs(step[0]), abs(step[1])) < mp.mpf("1e-42"):
            break
    field, fourth = homotopy_field(coefficient, time, x)
    require(
        max(abs(field[0]), abs(field[1])) < tolerance,
        "t=0 crossing solve did not converge",
    )
    return coefficient, x, field, fourth, iteration + 1


def continuation_row(
    coefficient: mp.mpf,
    time: mp.mpf,
    x: mp.mpf,
    field: list[mp.mpf],
    fourth: list[mp.mpf],
    iterations: int,
) -> dict:
    tangent = implicit_tangent(field, fourth)
    return {
        "lambda": decimal(coefficient),
        "time": decimal(time),
        "x": decimal(x),
        "field_jets_0_to_3": [decimal(value) for value in field],
        "fourth_component_jets_0_to_3": [decimal(value) for value in fourth],
        "jacobian_determinant": decimal(-(field[2] ** 2)),
        "dt_dlambda": decimal(tangent[0]),
        "dx_dlambda": decimal(tangent[1]),
        "newton_iterations": iterations,
    }


def cutoff_tail_log10_bound(index: int, derivative_order: int) -> float:
    epsilon = float(T_MAX_TEXT) / (4.0 * math.e * math.e)
    alpha = math.pi * index * index - epsilon
    y0 = math.exp(4.0 * float(CUTOFF))
    log_bound = (
        2.0 * math.log(math.pi)
        + 4.0 * math.log(index)
        + math.lgamma(derivative_order + 1.0)
        - math.log(2.0)
        + 1.5 * math.log(y0)
        - alpha * y0
        - math.log(alpha - 1.5 / y0)
    )
    return log_bound / math.log(10.0)


def build_payload(baseline: list[float], priority_lowered: bool) -> dict:
    started = perf_counter()
    monitor = CpuMonitor()
    lambda_step = mp.mpf(LAMBDA_STEP_TEXT)
    lambda_last = mp.mpf(LAMBDA_LAST_TEXT)
    seed_time = mp.mpf("0.432004524260958805434747496521")
    seed_x = mp.mpf("135.561632042676054682681116537")
    time, x, field, fourth, iterations = correct_contact(
        mp.mpf(0), seed_time, seed_x
    )
    rows = [continuation_row(mp.mpf(0), time, x, field, fourth, iterations)]
    tangent = implicit_tangent(field, fourth)
    parked = False

    coefficient = mp.mpf(0)
    try:
        while coefficient < lambda_last:
            next_coefficient = coefficient + lambda_step
            guess_time = time + lambda_step * tangent[0]
            guess_x = x + lambda_step * tangent[1]
            time, x, field, fourth, iterations = correct_contact(
                next_coefficient, guess_time, guess_x
            )
            coefficient = next_coefficient
            tangent = implicit_tangent(field, fourth)
            rows.append(
                continuation_row(
                    coefficient, time, x, field, fourth, iterations
                )
            )
            monitor.sample()
    except ResourcePark:
        parked = True

    crossing: dict | None = None
    if not parked:
        cross_lambda, cross_x, cross_field, cross_fourth, cross_iterations = (
            solve_boundary_crossing(
                mp.mpf("0.2017"), mp.mpf("135.862")
            )
        )
        cross_tangent = implicit_tangent(cross_field, cross_fourth)
        crossing = {
            "lambda": decimal(cross_lambda),
            "time": "0",
            "x": decimal(cross_x),
            "field_jets_0_to_3": [decimal(value) for value in cross_field],
            "fourth_component_jets_0_to_3": [
                decimal(value) for value in cross_fourth
            ],
            "jacobian_determinant": decimal(-(cross_field[2] ** 2)),
            "dt_dlambda": decimal(cross_tangent[0]),
            "dx_dlambda": decimal(cross_tangent[1]),
            "newton_iterations": cross_iterations,
        }

    completed_times = [mp.mpf(row["time"]) for row in rows]
    completed_xs = [mp.mpf(row["x"]) for row in rows]
    monotone_time = all(
        completed_times[index + 1] < completed_times[index]
        for index in range(len(completed_times) - 1)
    )
    monotone_x = all(
        completed_xs[index + 1] > completed_xs[index]
        for index in range(len(completed_xs) - 1)
    )
    tails = {
        str(index): {
            str(order): cutoff_tail_log10_bound(index, order)
            for order in range(4)
        }
        for index in range(1, 5)
    }
    elapsed = perf_counter() - started
    return {
        "kind": STEM,
        "schema_version": SCHEMA_VERSION,
        "date": DATE,
        "status": (
            "parked after resource threshold"
            if parked
            else "high-precision fourth-summand tangency homotopy scout complete"
        ),
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "parent": {
            "result": str(PARENT_RESULT.relative_to(REPO_ROOT)).replace("\\", "/"),
            "result_sha256": sha256_path(PARENT_RESULT),
        },
        "resource_policy": {
            "mode": "daytime",
            "active_compute_workers": 1,
            "thread_caps": 1,
            "below_normal_priority_applied": priority_lowered,
            "baseline_cpu_percent": baseline,
            "baseline_mean_percent": sum(baseline) / len(baseline),
            "runtime_cpu_percent": monitor.samples,
            "resource_parked": parked,
            "elapsed_seconds": elapsed,
        },
        "precision": {
            "decimal_digits": DPS,
            "integration_interval": [0, 2],
            "cutoff_tail_log10_bounds": tails,
            "interval_certified": False,
        },
        "exact_setup": {
            "component": (
                "phi_n(u)=pi*n^2*exp(5u)*(2*pi*n^2*exp(4u)-3)"
                "*exp(-pi*n^2*exp(4u))"
            ),
            "component_transform": (
                "H_(n,t)(x)=integral_0^infinity exp(tu^2)phi_n(u)cos(xu)du"
            ),
            "homotopy": (
                "F_lambda=H_1+H_2+H_3+lambda*H_4, 0<=lambda<=1"
            ),
            "contact": "F_lambda(t,x)=0 and partial_x F_lambda(t,x)=0",
            "heat_identities": (
                "partial_t F=-partial_x^2 F and "
                "partial_t partial_x F=-partial_x^3 F"
            ),
            "contact_jacobian": "[[-F_xx,0],[-F_xxx,F_xx]]",
            "contact_jacobian_determinant": "-F_xx^2",
            "implicit_response": (
                "[dt/dlambda,dx/dlambda]^T=-J^(-1)[H_4,H_4']^T"
            ),
            "pi_provenance": (
                "pi is the canonical Jacobi-theta constant in "
                "sum_n exp(-pi*n^2*y), here y=exp(4u). It is fixed by the "
                "theta modular normalization inherited by Xi; it is not inserted "
                "from an arbitrary circle or polygon. Replacing it changes the "
                "theta kernel and no longer represents the same zeta problem."
            ),
        },
        "continuation": {
            "lambda_step": decimal(lambda_step),
            "rows": rows,
            "time_strictly_decreasing_on_rows": monotone_time,
            "x_strictly_increasing_on_rows": monotone_x,
            "tracked_last_lambda": rows[-1]["lambda"],
        },
        "boundary_crossing": crossing,
        "diagnostic_decision": {
            "tracked_branch_exit": (
                None
                if crossing is None
                else (
                    "The regular branch issuing from the lambda=0 three-summand "
                    "contact crosses t=0 transversely before lambda=1."
                )
            ),
            "fixed_root_comparison": (
                "The fixed-root tangency-chain theorem does not apply: x drifts "
                "strictly along the sampled branch, and the initial response is order one."
            ),
            "next_theorem_target": (
                "Interval-certify branch existence, uniqueness, F_xx nonvanishing, "
                "and dt/dlambda<0 through the t=0 crossing; then count all contact "
                "branches in the compact heat rectangle before adding later summands."
            ),
            "proof_boundary": (
                "This is a high-precision, analytically tailed point continuation, "
                "not interval continuation. It tracks one branch only. It does not "
                "exclude another four-summand branch, prove no re-entry for lambda>"
                "the crossing, control the complete theta sum, exclude an Xi contact, "
                "prove Lambda<=0, RH, or a prize-level conclusion."
            ),
        },
    }


def render_note(payload: dict) -> str:
    rows = payload["continuation"]["rows"]
    crossing = payload["boundary_crossing"]
    lines = [
        "# Fourth Theta-Summand Tangency Homotopy Scout",
        "",
        f"Date: {DATE}",
        "",
        "Status: high-precision branch diagnostic. This is not a proof of RH or contact exclusion.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Exact Homotopy",
        "",
        "```text",
        "F_lambda(t,x)=H_(1,t)(x)+H_(2,t)(x)+H_(3,t)(x)+lambda H_(4,t)(x),",
        "F_lambda=0,  partial_x F_lambda=0.",
        "```",
        "",
        "At a contact, heat flow gives the exact Jacobian and response",
        "",
        "```text",
        "J=[[-F_xx,0],[-F_xxx,F_xx]],  det J=-F_xx^2,",
        "(t',x')=-J^(-1)(H_4,H_4').",
        "```",
        "",
        "## Pi Provenance",
        "",
        payload["exact_setup"]["pi_provenance"],
        "",
        "## Continuation",
        "",
        "| lambda | t | x | dt/dlambda | dx/dlambda | F_xx |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {mp.nstr(mp.mpf(row['lambda']), 5)} | "
            f"{mp.nstr(mp.mpf(row['time']), 12)} | "
            f"{mp.nstr(mp.mpf(row['x']), 14)} | "
            f"{mp.nstr(mp.mpf(row['dt_dlambda']), 10)} | "
            f"{mp.nstr(mp.mpf(row['dx_dlambda']), 10)} | "
            f"{mp.nstr(mp.mpf(row['field_jets_0_to_3'][2]), 8)} |"
        )
    if crossing is not None:
        lines.extend(
            [
                "",
                "## Boundary Crossing",
                "",
                f"The tracked branch reaches `t=0` at `lambda="
                f"{mp.nstr(mp.mpf(crossing['lambda']), 18)}` and "
                f"`x={mp.nstr(mp.mpf(crossing['x']), 18)}`.",
                f"The transverse derivative is `dt/dlambda="
                f"{mp.nstr(mp.mpf(crossing['dt_dlambda']), 14)}`.",
                "Thus this branch exits nonnegative heat time after only about one",
                "fifth of the true fourth coefficient has been switched on.",
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation",
            "",
            "This supplies a concrete drifting-root mechanism. The fixed-root",
            "tangency-chain theorem does not cover it, because the contact location",
            "moves and the fourth-summand response is order one despite its small mass.",
            "The next serious theorem is an interval branch-and-count argument, not a",
            "claim that this one tracked branch exhausts the four-summand problem.",
            "",
            "## Proof Boundary",
            "",
            payload["diagnostic_decision"]["proof_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    parser.add_argument("--dps", type=int, default=DPS)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    priority_lowered = request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "theta homotopy baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, mean={baseline_mean:.1f}%"
    )
    if baseline_mean > args.baseline_max_percent:
        print("theta homotopy scout deferred: daytime baseline is already busy")
        return 2
    mp.mp.dps = args.dps
    payload = build_payload(baseline, priority_lowered)
    write_json_atomic(args.out, payload)
    write_text_atomic(args.note, render_note(payload))
    crossing = payload["boundary_crossing"]
    crossing_summary = "unfinished" if crossing is None else crossing["lambda"]
    print(
        "built fourth-summand tangency homotopy: "
        f"rows={len(payload['continuation']['rows'])}, "
        f"t=0 crossing lambda={crossing_summary}, "
        f"parked={payload['resource_policy']['resource_parked']}, "
        f"elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 3 if payload["resource_policy"]["resource_parked"] else 0


if __name__ == "__main__":
    sys.exit(main())
