#!/usr/bin/env python3
"""Project physical shifted-Hardy kernels onto exact low-mode nullspaces."""

from __future__ import annotations

import argparse
import ctypes
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import mpmath as mp
import numpy as np

from jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate import (
    CHEBYSHEV_POINTS,
    DEFAULT_SHIFT_STEP,
    MIDPOINT_POINTS,
    evaluate_kernel,
    hardy_phases,
    parse_complex,
    parse_real_part,
    polynomial,
    shifts,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_hardy_exact_low_mode_projected_kernel_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{KIND}.md"
BRIDGE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate.json"
)
JOINT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_hardy_joint_correlated_error_ledger_gate.json"
)
HANDOFF_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_hardy_modewise_derivative_handoff_gate.json"
)
ROOTS_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_calibrated_roster_interior_phase_root_scout.json"
)
ADAPTER_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_physical_scalar_edge_adapter_gate.json"
)
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts"
    / "check_jensen_window_pf_newman_c1_hardy_exact_low_mode_projected_kernel_gate.py"
)

ACTIVE_LOW_MODES = {
    ("base_rows", "P_V"): (0, 2, 4),
    ("base_rows", "P_A"): (1, 3, 5),
    ("base_rows", "P_Q"): (1, 3, 5),
    ("base_rows", "P_N"): (2, 4),
    ("tangent_rows", "P_V"): (1, 3, 5),
    ("tangent_rows", "P_A"): (2, 4),
    ("tangent_rows", "P_Q"): (2, 4),
    ("tangent_rows", "P_N"): (3, 5),
}
HEIGHT_KEYS = ("t1e10_et005", "t1e12_et005")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def set_below_normal_priority() -> None:
    if os.name != "nt":
        return
    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint32)
        kernel32.SetPriorityClass.restype = ctypes.c_int
        kernel32.SetPriorityClass(kernel32.GetCurrentProcess(), 0x00004000)
    except (AttributeError, OSError):
        pass


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative_path(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT)).replace("\\", "/")


def load_decimal_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"), parse_float=Decimal)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_project(
    coefficients: list[Fraction],
    inactive_modes: tuple[int, ...],
    modes: list[list[int]],
) -> list[Fraction]:
    projected = coefficients[:]
    for mode_index in inactive_modes:
        mode = modes[mode_index]
        norm_squared = sum(value * value for value in mode)
        moment = sum(
            coefficient * value
            for coefficient, value in zip(projected, mode, strict=True)
        )
        projected = [
            coefficient - moment * Fraction(value, norm_squared)
            for coefficient, value in zip(projected, mode, strict=True)
        ]
    return projected


def exact_factor(record: dict) -> float:
    return float(Fraction(record["numerator"])) / math.sqrt(
        int(record["norm_squared"])
    )


def derivative_budget(
    coefficients: list[Fraction], handoff: dict
) -> tuple[list[float], list[Fraction]]:
    mode_rows = handoff["orthogonal_modes"]
    moments: list[Fraction] = []
    normalized_moments: list[float] = []
    for row in mode_rows:
        mode = row["integer_vector"]
        moment = sum(
            coefficient * value
            for coefficient, value in zip(coefficients, mode, strict=True)
        )
        moments.append(moment)
        normalized_moments.append(
            float(moment) / math.sqrt(int(row["norm_squared"]))
        )

    budgets = []
    for derivative in range(6):
        budgets.append(
            math.fsum(
                abs(normalized_moments[index])
                * (
                    exact_factor(row["derivative_weights"][str(derivative)])
                    if str(derivative) in row["derivative_weights"]
                    else 0.0
                )
                for index, row in enumerate(mode_rows)
            )
        )
    rough = math.fsum(
        abs(normalized_moments[index])
        * exact_factor(row["sixth_derivative_remainder_weight"])
        for index, row in enumerate(mode_rows)
    )
    coefficient_l2 = math.sqrt(
        math.fsum(float(value) ** 2 for value in coefficients)
    )
    tail = math.sqrt(Fraction(handoff["rough_tail_weight"]["exact_squared"]))
    budgets.append(rough + coefficient_l2 * tail)
    return budgets, moments


def held_out_grid(upper: float) -> np.ndarray:
    midpoint = (np.arange(MIDPOINT_POINTS) + 0.5) * upper / MIDPOINT_POINTS
    chebyshev = 0.5 * upper * (
        1.0 - np.cos(np.linspace(0.0, np.pi, CHEBYSHEV_POINTS))
    )
    return np.unique(np.concatenate((midpoint, chebyshev, [0.0, upper])))


def source_audit() -> dict:
    paths = {
        "shifted_hardy_bridge": BRIDGE_RESULT,
        "joint_error_ledger": JOINT_RESULT,
        "derivative_handoff": HANDOFF_RESULT,
        "calibrated_roots": ROOTS_RESULT,
        "physical_adapter": ADAPTER_RESULT,
        "builder": BUILDER,
        "checker": CHECKER,
    }
    return {
        key: {"path": relative_path(path), "sha256": file_hash(path)}
        for key, path in paths.items()
    }


def build_rows_and_summary() -> tuple[list[dict], dict]:
    bridge = load_decimal_json(BRIDGE_RESULT)
    joint = load_decimal_json(JOINT_RESULT)
    handoff = load_decimal_json(HANDOFF_RESULT)["exact_mode_handoff"]
    roots_payload = json.loads(ROOTS_RESULT.read_text(encoding="utf-8"))
    adapter = json.loads(ADAPTER_RESULT.read_text(encoding="utf-8"))
    roots = [
        roots_payload["precision_ladder"][-1],
        *roots_payload["additional_roots"],
    ]
    physical_rows = adapter["physical_rows"]
    modes = [row["integer_vector"] for row in handoff["orthogonal_modes"]]
    errors = {
        key: [float(value) for value in joint["residual_modes"][key]["values"]]
        for key in HEIGHT_KEYS
    }
    require(len(roots) == len(physical_rows) == len(bridge["physical_kernels"]) == 3, "root alignment failed")

    mp.mp.dps = 180
    shift_values = shifts(DEFAULT_SHIFT_STEP)
    rows = []
    for root_index, (root_row, physical, bridge_root) in enumerate(
        zip(roots, physical_rows, bridge["physical_kernels"], strict=True)
    ):
        root = root_row["root"]
        t_value = mp.mpf(physical["chart"]["T_0"])
        n_value = int(physical["N"])
        heat_time = float(physical["chart"]["t"])
        sigma = parse_real_part(root["s_star"])
        omega = parse_complex(root["eta"])
        if n_value % 2:
            omega = -omega
        upper = math.log(n_value)
        test = held_out_grid(upper)
        test_weight = np.exp(
            heat_time * test**2 / 4.0 - (sigma - 0.5) * test
        )
        phases = hardy_phases(t_value, shift_values)
        stored = {
            (row["family"], row["name"]): row for row in bridge_root["rows"]
        }
        for family in ("base_rows", "tangent_rows"):
            for name, values in physical["polynomials"][family].items():
                source_row = stored[(family, name)]
                coefficients = [
                    Fraction(value) for value in source_row["real_hardy_coefficients"]
                ]
                active = ACTIVE_LOW_MODES[(family, name)]
                inactive = tuple(index for index in range(6) if index not in active)
                projected = exact_project(coefficients, inactive, modes)
                for mode_index in inactive:
                    require(
                        sum(
                            coefficient * value
                            for coefficient, value in zip(
                                projected, modes[mode_index], strict=True
                            )
                        )
                        == 0,
                        "exact inactive-mode projection failed",
                    )

                projected_float = np.asarray(
                    [float(value) for value in projected], dtype=np.float64
                )
                target_polynomial = polynomial(values)
                truth = (
                    omega
                    * test_weight
                    * np.polynomial.polynomial.polyval(test, target_polynomial)
                )
                prediction = evaluate_kernel(
                    phases, shift_values, test, projected_float
                )
                target_scale = max(1.0, float(np.max(np.abs(truth))))
                held_out_error = float(
                    np.max(np.abs(prediction - truth)) / target_scale
                )
                budgets, moments = derivative_budget(projected, handoff)
                original_float = np.asarray(
                    [float(value) for value in coefficients], dtype=np.float64
                )
                projections = {
                    key: abs(
                        math.fsum(
                            coefficient * error
                            for coefficient, error in zip(
                                projected_float, errors[key], strict=True
                            )
                        )
                    )
                    for key in HEIGHT_KEYS
                }
                ratio = projections[HEIGHT_KEYS[1]] / max(
                    projections[HEIGHT_KEYS[0]], 1.0e-300
                )
                rows.append(
                    {
                        "root_index": root_index,
                        "N": n_value,
                        "family": family,
                        "name": name,
                        "active_low_modes": list(active),
                        "exact_annihilated_low_modes": list(inactive),
                        "exact_projected_coefficients": [
                            fraction_text(value) for value in projected
                        ],
                        "projected_coefficients_decimal": [
                            format(float(value), ".17g") for value in projected
                        ],
                        "exact_low_mode_moments": [
                            fraction_text(value) for value in moments
                        ],
                        "coefficient_l1": math.fsum(
                            abs(float(value)) for value in projected
                        ),
                        "coefficient_l2": math.sqrt(
                            math.fsum(float(value) ** 2 for value in projected)
                        ),
                        "coefficient_change_l2": float(
                            np.linalg.norm(projected_float - original_float)
                        ),
                        "original_held_out_relative_error": float(
                            source_row["held_out_relative_error"]
                        ),
                        "projected_held_out_relative_error": held_out_error,
                        "target_scale": target_scale,
                        "derivative_budget_multipliers": budgets,
                        "normalized_derivative_budget_multipliers": [
                            value / target_scale for value in budgets
                        ],
                        "residual_projections": projections,
                        "absolute_projection_height_ratio": ratio,
                        "projection_increased_at_t1e12": ratio > 1.0,
                    }
                )

    require(len(rows) == 24, "physical projected-row count failed")
    exact_annihilations = sum(
        len(row["exact_annihilated_low_modes"]) for row in rows
    )
    require(exact_annihilations == 84, "exact low-mode constraint count failed")
    derivative_maxima = []
    for derivative in range(7):
        maximum_row = max(
            rows,
            key=lambda row: row["normalized_derivative_budget_multipliers"][
                derivative
            ],
        )
        derivative_maxima.append(
            {
                "derivative_order": derivative,
                "maximum_absolute_multiplier": max(
                    row["derivative_budget_multipliers"][derivative]
                    for row in rows
                ),
                "maximum_normalized_multiplier": maximum_row[
                    "normalized_derivative_budget_multipliers"
                ][derivative],
                "attained_at": {
                    "root_index": maximum_row["root_index"],
                    "family": maximum_row["family"],
                    "name": maximum_row["name"],
                },
            }
        )
    summary = {
        "physical_rows": len(rows),
        "exact_low_mode_annihilations": exact_annihilations,
        "maximum_coefficient_change_l2": max(
            row["coefficient_change_l2"] for row in rows
        ),
        "maximum_projected_coefficient_l1": max(
            row["coefficient_l1"] for row in rows
        ),
        "maximum_projected_held_out_relative_error": max(
            row["projected_held_out_relative_error"] for row in rows
        ),
        "increased_projection_rows": sum(
            row["projection_increased_at_t1e12"] for row in rows
        ),
        "maximum_absolute_projection_height_ratio": max(
            row["absolute_projection_height_ratio"] for row in rows
        ),
        "maximum_residual_projection": {
            key: max(row["residual_projections"][key] for row in rows)
            for key in HEIGHT_KEYS
        },
        "normalized_derivative_budget_maxima": derivative_maxima,
    }
    return rows, summary


def conceptual_rows() -> list[dict]:
    return [
        {
            "id": "helpk_01_saved_kernels",
            "role": "source_audit",
            "readiness": "available_diagnostic",
            "claim": "All twenty-four accepted h=0.01 physical kernel vectors are retained as the starting point.",
            "proof_boundary": "The starting vectors remain floating-point fits, not interval objects.",
        },
        {
            "id": "helpk_02_integer_modes",
            "role": "exact_input",
            "readiness": "available_exact",
            "claim": "The first six shift modes are primitive integer Gram vectors.",
            "proof_boundary": "Finite-grid algebra only.",
        },
        {
            "id": "helpk_03_rational_projection",
            "role": "exact_lemma",
            "readiness": "available_exact",
            "claim": "Orthogonal subtraction projects each finite-decimal coefficient vector to a rational nullspace vector.",
            "proof_boundary": "This changes the diagnostic kernel representation; it does not alter the evaluator error.",
        },
        {
            "id": "helpk_04_exact_annihilation",
            "role": "finite_certificate",
            "readiness": "available_exact",
            "claim": "Eighty-four selected low-mode moments vanish as exact rational identities.",
            "proof_boundary": "Only the listed modes and rows are covered.",
        },
        {
            "id": "helpk_05_fit_replay",
            "role": "finite_diagnostic",
            "readiness": "diagnostic_validated",
            "claim": "Endpoint-inclusive held-out replay remains below 1e-10 after projection.",
            "proof_boundary": "Floating-point held-out diagnostics are not interval approximation theorems.",
        },
        {
            "id": "helpk_06_derivative_budget",
            "role": "conditional_handoff",
            "readiness": "available_conditional",
            "claim": "The exact moment zeros remove the selected derivative channels from the Section 11.237 bound.",
            "proof_boundary": "The external derivative bounds remain open.",
        },
        {
            "id": "helpk_07_scalar_guard",
            "role": "nonpromotion_guard",
            "readiness": "guard_validated",
            "claim": "The projected kernels preserve the 21/24 two-height scalar-extrapolation obstruction.",
            "proof_boundary": "Exact mode cleanup does not validate extrapolation.",
        },
        {
            "id": "helpk_08_physical_guard",
            "role": "nonpromotion_guard",
            "readiness": "mandatory",
            "claim": "No physical-height value is evaluated or inferred.",
            "proof_boundary": "The physical run remains behind the external error theorem.",
        },
    ]


def build_artifact() -> dict:
    rows, summary = build_rows_and_summary()
    success = (
        "built Newman C1 Hardy exact low-mode projected-kernel gate: "
        f"{summary['physical_rows']} rows, "
        f"{summary['exact_low_mode_annihilations']} exact annihilations, "
        f"max held-out error {summary['maximum_projected_held_out_relative_error']:.12g}, "
        f"{summary['increased_projection_rows']}/24 projected errors still increase, "
        "0 physical values and 1 open derivative theorem"
    )
    return {
        "kind": KIND,
        "status": "exact_low_mode_projection_with_diagnostic_kernel_replay",
        "active_low_mode_registry": {
            f"{family}:{name}": list(modes)
            for (family, name), modes in ACTIVE_LOW_MODES.items()
        },
        "rows": rows,
        "summary": summary,
        "conceptual_rows": conceptual_rows(),
        "source_audit": source_audit(),
        "route_decision": (
            "Use the exact projected coefficient vectors for subsequent error budgeting. "
            "Their selected low-mode cancellations are identities, while the projected "
            "kernel approximation remains a diagnostic until interval-certified."
        ),
        "proof_boundary": (
            "This gate proves eighty-four exact rational low-mode annihilations and "
            "replays finite kernel and two-height diagnostics. It does not prove the "
            "external evaluator derivative theorem, an interval kernel approximation, "
            "a physical carrier value, retained observation, determinant sign or bound, "
            "complete-current inequality, all-q transport, contact exclusion, Lambda<=0, "
            "PF-infinity, RH, or a prize-level conclusion."
        ),
        "success": success,
    }


def render_note(artifact: dict) -> str:
    summary = artifact["summary"]
    derivative_rows = summary["normalized_derivative_budget_maxima"]
    lines = [
        "# Newman C1 Hardy Exact Low-Mode Projected-Kernel Gate",
        "",
        "Date: 2026-08-05",
        "",
        f"Status: `{artifact['status']}`; not a proof of the external error theorem, RH, or a prize-level conclusion.",
        "",
        "## Exact projection",
        "",
        "The twenty-four saved physical coefficient vectors already showed a sparse "
        "low-mode pattern, but only to floating precision.  For each row, this gate "
        "defines the replacement vector",
        "",
        "```text",
        "c_tilde=c-sum_(k in Z) <c,p_k> p_k/||p_k||_2^2,",
        "```",
        "",
        "where the `p_k` are the primitive integer Gram vectors and `Z` is the row's "
        "listed inactive-mode set.  The saved decimals are finite rationals, so every "
        "entry of `c_tilde` and every identity `<c_tilde,p_k>=0` is exact.",
        "",
        f"Across `{summary['physical_rows']}` rows this enforces "
        f"`{summary['exact_low_mode_annihilations']}` exact annihilations.  The largest "
        f"coefficient-vector L2 change is `{summary['maximum_coefficient_change_l2']:.12g}`.",
        "",
        "## Kernel replay",
        "",
        "The complete midpoint, Chebyshev, and endpoint held-out grids were replayed at "
        "180-digit phase precision.  The maximum projected-kernel relative error is "
        f"`{summary['maximum_projected_held_out_relative_error']:.12g}`, while the largest "
        f"projected coefficient L1 norm remains `{summary['maximum_projected_coefficient_l1']:.12g}`. "
        "This is a floating diagnostic, not an interval certificate.",
        "",
        "## Derivative sensitivity",
        "",
        "Substituting the exact moment zeros into Section 11.237 gives the following "
        "largest coefficient multipliers after normalization by each target scale:",
        "",
        "| derivative order | maximum normalized multiplier | row |",
        "|---:|---:|---|",
    ]
    for row in derivative_rows:
        where = row["attained_at"]
        lines.append(
            f"| {row['derivative_order']} | {row['maximum_normalized_multiplier']:.12g} | "
            f"root {where['root_index']} `{where['family']}:{where['name']}` |"
        )
    lines.extend(
        [
            "",
            "Thus the million-scale independent L1 amplification is not the relevant "
            "structured derivative budget.  The exact low-mode design reduces the worst "
            "normalized multipliers for orders zero through six to the values above.  It "
            "still does not supply the missing external `D_r(T)` bounds.",
            "",
            "## Falsification retained",
            "",
            f"At the two calibration heights, `{summary['increased_projection_rows']}/24` "
            "absolute projected residuals still increase.  The largest increase factor is "
            f"`{summary['maximum_absolute_projection_height_ratio']:.12g}`.  Exact moment "
            "cleanup therefore does not revive the rejected scalar-extrapolation route.",
            "",
            "No new `pi` is introduced.  The grid is the evaluator's literal `0.01` shift "
            "grid, and the projection uses rational finite-grid inner products only.",
            "",
            artifact["proof_boundary"],
            "",
            "```text",
            artifact["success"],
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    set_below_normal_priority()
    artifact = build_artifact()
    atomic_write(args.out, json.dumps(artifact, indent=2) + "\n")
    atomic_write(args.note, render_note(artifact))
    print(artifact["success"])


if __name__ == "__main__":
    main()
