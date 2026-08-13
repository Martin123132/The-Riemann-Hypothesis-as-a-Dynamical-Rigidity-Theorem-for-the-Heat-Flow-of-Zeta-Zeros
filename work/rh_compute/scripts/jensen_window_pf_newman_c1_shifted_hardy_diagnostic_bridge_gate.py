#!/usr/bin/env python3
"""Build the shifted-Hardy diagnostic bridge and stability gate."""

from __future__ import annotations

import ctypes
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re

import mpmath as mp
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
EXTERNAL_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode"
FIXTURE_ROOT = EXTERNAL_ROOT / "fixtures"

SOURCE_PATHS = {
    "calibrated_roots": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_calibrated_roster_interior_phase_root_scout.json",
    "physical_adapter": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_physical_scalar_edge_adapter_gate.json",
    "low_n_oracle": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_fixed_cell_low_n_determinant_evaluator_gate.json",
    "ridge_inversion": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_ideal_cubic_"
        "double_morse_ridge_carrier_inversion_gate.json"
    ),
}

PATCH_PATHS = {
    "single_mpi_sum": EXTERNAL_ROOT / "patches/0001-initialize-mpi-sum.patch",
    "single_finalize": EXTERNAL_ROOT / "patches/0002-finalize-mpi-before-stop.patch",
    "multi_finalize": EXTERNAL_ROOT / "patches/0003-finalize-multi-mpi-before-stop.patch",
    "multi_openmpi": EXTERNAL_ROOT / "patches/0004-fix-multi-openmpi-interface.patch",
    "spacing_004": EXTERNAL_ROOT / "patches/0005-multi-shift-step-004.patch",
    "spacing_002": EXTERNAL_ROOT / "patches/0006-multi-shift-step-002.patch",
}

FIXTURE_PATHS = {
    "single_t1e10_et05": FIXTURE_ROOT
    / "t1e10/zeta14cubic_single_et05_output.txt",
    "multi_t1e10_et05": FIXTURE_ROOT
    / "t1e10/zeta14cubicmult_et05_output.txt",
    "multi_t1e10_et005": FIXTURE_ROOT
    / "t1e10/zeta14cubicmult_et005_output.txt",
    "multi_t1e10_et005_h002": FIXTURE_ROOT
    / "t1e10/zeta14cubicmult_et005_h002_output.txt",
    "multi_t1e10_et005_h004": FIXTURE_ROOT
    / "t1e10/zeta14cubicmult_et005_h004_output.txt",
    "multi_t1e12_et005": FIXTURE_ROOT
    / "t1e12/zeta14cubicmult_et005_output.txt",
}

UPSTREAM_REPOSITORY = "https://github.com/dml2391/Hardy-function-fastcodes"
UPSTREAM_COMMIT = "2e16dac3206b707052c3ac4cacdf3d1a2325e636"
UPSTREAM_LICENSE = "GPL-3.0"
SHIFT_COUNT = 15
DEFAULT_SHIFT_STEP = 0.01
TERMINAL_PEEL = 0
TRAIN_POINTS = 2401
MIDPOINT_POINTS = 4800
CHEBYSHEV_POINTS = 2001


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


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


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def relative_path(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT)).replace("\\", "/")


def parse_complex(text: str) -> complex:
    match = re.fullmatch(
        r"\(([^ ]+)\s*([+-])\s*([^j]+)j\)", text.strip()
    )
    require(match is not None, f"cannot parse complex value: {text}")
    sign = 1.0 if match.group(2) == "+" else -1.0
    return complex(float(match.group(1)), sign * float(match.group(3)))


def parse_real_part(text: str) -> float:
    match = re.match(r"\(([^ ]+)", text.strip())
    require(match is not None, f"cannot parse real part: {text}")
    return float(match.group(1))


def parse_grand_totals(path: Path) -> np.ndarray:
    values = re.findall(
        r"Grand total of Hardy function Z\(t[^=]*=\s*([-+0-9.Ee]+)",
        path.read_text(encoding="utf-8"),
    )
    require(len(values) == SHIFT_COUNT, f"expected 15 Hardy values in {path}")
    return np.asarray([float(value) for value in values], dtype=np.float64)


def parse_single_total(path: Path) -> float:
    values = re.findall(
        r"Grand total of Hardy function Z\(t\)=\s*([-+0-9.Ee]+)",
        path.read_text(encoding="utf-8"),
    )
    require(len(values) == 1, "single fixture Hardy value count")
    return float(values[0])


def load_and_audit_sources() -> tuple[dict, dict, dict[str, dict[str, str]]]:
    payloads = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }
    require(
        payloads["calibrated_roots"]["summary"]["calibrated_roots"] == 3,
        "physical-root count drifted",
    )
    require(
        payloads["physical_adapter"]["summary"]["compiled_base_polynomials"]
        == 12,
        "physical base-row count drifted",
    )
    require(
        payloads["physical_adapter"]["summary"]["compiled_tangent_polynomials"]
        == 12,
        "physical tangent-row count drifted",
    )
    require(
        payloads["low_n_oracle"]["summary"]["fixtures"] == 3,
        "low-N oracle drifted",
    )
    ridge_text = json.dumps(payloads["ridge_inversion"], sort_keys=True)
    require(
        "2mathscr M_N-mathscr F_N=mathscr M_N-tau_band+E_band" in ridge_text,
        "carrier inversion identity drifted",
    )

    audit: dict[str, dict[str, str]] = {}
    for group, paths in (
        ("source", SOURCE_PATHS),
        ("patch", PATCH_PATHS),
        ("fixture", FIXTURE_PATHS),
    ):
        for key, path in paths.items():
            require(path.is_file(), f"missing {group} artifact: {path}")
            audit[f"{group}:{key}"] = {
                "path": relative_path(path),
                "sha256": file_hash(path),
            }
    return payloads["calibrated_roots"], payloads["physical_adapter"], audit


def shifts(step: float) -> np.ndarray:
    return np.arange(-7, 8, dtype=np.float64) * step


def hardy_phases(t_value: mp.mpf, shift_values: np.ndarray) -> np.ndarray:
    phases = []
    for delta in shift_values:
        shifted = t_value + mp.mpf(str(float(delta)))
        phases.append(complex(mp.e ** (-1j * mp.siegeltheta(shifted))))
    return np.asarray(phases, dtype=np.complex128)


def leading_rs_correction(t_value: mp.mpf) -> float:
    a_value = mp.sqrt(t_value / (2 * mp.pi))
    n_value = int(mp.floor(a_value))
    fraction = a_value - n_value
    numerator = mp.cos(
        2 * mp.pi * (fraction * (fraction - 1) - mp.mpf(1) / 16)
    )
    denominator = mp.cos(2 * mp.pi * fraction)
    value = (-1) ** (n_value - 1) * a_value ** (-mp.mpf("0.5"))
    return float(value * numerator / denominator)


def direct_carrier(t_text: str) -> dict:
    t_value = mp.mpf(t_text)
    n_value = int(mp.floor(mp.sqrt(t_value / (2 * mp.pi))))
    integers = np.arange(1, n_value + 1, dtype=np.float64)
    logs = np.log(integers)
    base = np.empty(n_value, dtype=np.complex128)
    for index in range(1, n_value + 1):
        phase = mp.e ** (1j * t_value * mp.log(index))
        base[index - 1] = complex(phase) / math.sqrt(index)
    return {
        "T": t_value,
        "N": n_value,
        "logs": logs,
        "base": base,
    }


def direct_main_samples(carrier: dict, step: float) -> dict:
    shift_values = shifts(step)
    phases = hardy_phases(carrier["T"], shift_values)
    values = []
    corrections = []
    for delta, phase in zip(shift_values, phases, strict=True):
        twist = np.exp(1j * delta * carrier["logs"])
        values.append(2.0 * float(np.real(phase * np.sum(carrier["base"] * twist))))
        corrections.append(
            leading_rs_correction(carrier["T"] + mp.mpf(str(float(delta))))
        )
    return {
        "shifts": shift_values,
        "phases": phases,
        "main": np.asarray(values),
        "correction": np.asarray(corrections),
    }


def error_summary(fast_main: np.ndarray, direct_main: np.ndarray) -> dict:
    error = fast_main - direct_main
    return {
        "errors": [float(value) for value in error],
        "mean": float(np.mean(error)),
        "standard_deviation": float(np.std(error)),
        "maximum_absolute": float(np.max(np.abs(error))),
    }


def fit_real_hardy_kernel(
    t_value: mp.mpf,
    shift_values: np.ndarray,
    train: np.ndarray,
    target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    phases = hardy_phases(t_value, shift_values)
    matrix = 2.0 * phases[None, :] * np.exp(
        1j * train[:, None] * shift_values[None, :]
    )
    real_matrix = np.vstack((matrix.real, matrix.imag))
    real_target = np.concatenate((target.real, target.imag))
    coefficients = np.linalg.lstsq(real_matrix, real_target, rcond=1.0e-14)[0]
    singular_values = np.linalg.svd(real_matrix, compute_uv=False)
    return coefficients, phases, singular_values


def evaluate_kernel(
    phases: np.ndarray,
    shift_values: np.ndarray,
    points: np.ndarray,
    coefficients: np.ndarray,
) -> np.ndarray:
    matrix = 2.0 * phases[None, :] * np.exp(
        1j * points[:, None] * shift_values[None, :]
    )
    return matrix @ coefficients


def polynomial(array: list[dict[str, str]]) -> np.ndarray:
    return np.asarray(
        [complex(float(value["re"]), float(value["im"])) for value in array],
        dtype=np.complex128,
    )


def physical_kernel_certificate(
    roots_payload: dict, adapter_payload: dict
) -> tuple[list[dict], dict]:
    roots = [
        roots_payload["precision_ladder"][-1],
        *roots_payload["additional_roots"],
    ]
    physical_rows = adapter_payload["physical_rows"]
    require(len(roots) == len(physical_rows) == 3, "physical row alignment")
    shift_values = shifts(DEFAULT_SHIFT_STEP)
    rows = []
    endpoint_obstructions = []
    for index, (root_row, physical) in enumerate(
        zip(roots, physical_rows, strict=True)
    ):
        root = root_row["root"]
        t_value = mp.mpf(physical["chart"]["T_0"])
        n_value = int(physical["N"])
        require(n_value == int(root["N"]), "physical N mismatch")
        heat_time = float(physical["chart"]["t"])
        sigma = parse_real_part(root["s_star"])
        omega = parse_complex(root["eta"])
        if n_value % 2:
            omega = -omega
        endpoint_line = complex(mp.e ** (-1j * mp.siegeltheta(t_value)))
        ratio = omega / endpoint_line
        endpoint_obstructions.append(abs(ratio.imag))

        upper = math.log(n_value - TERMINAL_PEEL)
        train = np.linspace(0.0, upper, TRAIN_POINTS)
        midpoint = (np.arange(MIDPOINT_POINTS) + 0.5) * upper / MIDPOINT_POINTS
        chebyshev = 0.5 * upper * (
            1.0 - np.cos(np.linspace(0.0, np.pi, CHEBYSHEV_POINTS))
        )
        test = np.unique(np.concatenate((midpoint, chebyshev, [0.0, upper])))
        train_weight = np.exp(
            heat_time * train**2 / 4.0 - (sigma - 0.5) * train
        )
        test_weight = np.exp(
            heat_time * test**2 / 4.0 - (sigma - 0.5) * test
        )

        root_rows = []
        for family in ("base_rows", "tangent_rows"):
            for name, values in physical["polynomials"][family].items():
                coefficients = polynomial(values)
                target = (
                    omega
                    * train_weight
                    * np.polynomial.polynomial.polyval(train, coefficients)
                )
                fit, phases, singular_values = fit_real_hardy_kernel(
                    t_value, shift_values, train, target
                )
                truth = (
                    omega
                    * test_weight
                    * np.polynomial.polynomial.polyval(test, coefficients)
                )
                prediction = evaluate_kernel(phases, shift_values, test, fit)
                scale = max(1.0, float(np.max(np.abs(truth))))
                error = float(np.max(np.abs(prediction - truth)) / scale)
                root_rows.append(
                    {
                        "family": family,
                        "name": name,
                        "polynomial_degree": len(coefficients) - 1,
                        "real_hardy_coefficients": [float(value) for value in fit],
                        "held_out_relative_error": error,
                        "coefficient_l1": float(np.sum(np.abs(fit))),
                        "coefficient_l2": float(np.linalg.norm(fit)),
                        "target_scale": scale,
                        "design_condition_number": float(
                            singular_values[0] / singular_values[-1]
                        ),
                    }
                )
        rows.append(
            {
                "root_index": index,
                "N": n_value,
                "T_0": str(t_value),
                "shift_step": DEFAULT_SHIFT_STEP,
                "shifts": [float(value) for value in shift_values],
                "terminal_peel": TERMINAL_PEEL,
                "bulk_log_upper": upper,
                "endpoint_phase_ratio": {
                    "re": ratio.real,
                    "im": ratio.imag,
                },
                "unpeeled_endpoint_quadrature_obstruction": abs(ratio.imag),
                "rows": root_rows,
            }
        )

    all_rows = [row for root in rows for row in root["rows"]]
    summary = {
        "roots": len(rows),
        "physical_rows": len(all_rows),
        "terminal_peel": TERMINAL_PEEL,
        "maximum_held_out_relative_error": max(
            row["held_out_relative_error"] for row in all_rows
        ),
        "maximum_coefficient_l1": max(row["coefficient_l1"] for row in all_rows),
        "maximum_design_condition_number": max(
            row["design_condition_number"] for row in all_rows
        ),
        "maximum_unpeeled_endpoint_obstruction": max(endpoint_obstructions),
    }
    return rows, summary


def low_height_reconstruction(
    carrier: dict,
    samples: dict,
    fast_main: np.ndarray,
    physical_phase_ratio: complex,
) -> dict:
    shift_values = samples["shifts"]
    t_value = carrier["T"]
    n_value = carrier["N"]
    log_a = float(mp.log(mp.sqrt(t_value / (2 * mp.pi))))
    upper = math.log(n_value)
    train = np.linspace(0.0, upper, 1601)
    test = (np.arange(3200) + 0.5) * upper / 3200
    heat_time = 0.0002
    sigma = 0.5025
    unit = complex(mp.e ** (-1j * mp.siegeltheta(t_value))) * physical_phase_ratio
    train_weight = np.exp(
        heat_time * train**2 / 4.0 - (sigma - 0.5) * train
    )
    test_weight = np.exp(
        heat_time * test**2 / 4.0 - (sigma - 0.5) * test
    )
    source_weight = np.exp(
        heat_time * carrier["logs"] ** 2 / 4.0
        - (sigma - 0.5) * carrier["logs"]
    )
    rows = []
    for degree in range(6):
        target = unit * train_weight * (1j * (train - log_a)) ** degree
        coefficients, phases, _ = fit_real_hardy_kernel(
            t_value, shift_values, train, target
        )
        truth = unit * test_weight * (1j * (test - log_a)) ** degree
        prediction = evaluate_kernel(phases, shift_values, test, coefficients)
        grid_scale = max(1.0, float(np.max(np.abs(truth))))
        grid_error = float(np.max(np.abs(prediction - truth)) / grid_scale)

        terms = (
            carrier["base"]
            * unit
            * source_weight
            * (1j * (carrier["logs"] - log_a)) ** degree
        )
        endpoint = 0.5 * float(terms[0].real + terms[-1].real)
        direct = float(np.sum(terms.real) - endpoint)
        exact_reconstruction = float(coefficients @ samples["main"] - endpoint)
        fast_reconstruction = float(coefficients @ fast_main - endpoint)
        observation_scale = max(1.0, abs(direct))
        rows.append(
            {
                "degree": degree,
                "grid_relative_error": grid_error,
                "coefficient_l1": float(np.sum(np.abs(coefficients))),
                "direct_half_endpoint_observation": direct,
                "exact_main_relative_error": abs(exact_reconstruction - direct)
                / observation_scale,
                "fast_main_relative_error": abs(fast_reconstruction - direct)
                / observation_scale,
                "fast_minus_exact_reconstruction": fast_reconstruction
                - exact_reconstruction,
            }
        )
    return {
        "T": str(t_value),
        "N": n_value,
        "shift_step": DEFAULT_SHIFT_STEP,
        "rows": rows,
        "maximum_grid_relative_error": max(row["grid_relative_error"] for row in rows),
        "maximum_exact_main_relative_error": max(
            row["exact_main_relative_error"] for row in rows
        ),
        "maximum_fast_main_relative_error": max(
            row["fast_main_relative_error"] for row in rows
        ),
    }


def build_rows(exact: dict) -> list[GateRow]:
    return [
        GateRow("shd_01_pin", "upstream pin", "guard_validated", "The external Hardy implementation is pinned by repository, commit, and license.", f"{UPSTREAM_REPOSITORY} at {UPSTREAM_COMMIT}; {UPSTREAM_LICENSE}.", "The third-party source is not copied into this corpus."),
        GateRow("shd_02_port", "portability", "diagnostic_validated", "Four minimal patches make the single and multi programs compile and exit cleanly under OpenMPI.", "MPI_SUM initialization, reachable MPI_Finalize calls, removal of one dead undefined callback, and default-integer reduction counts.", "These patches certify interoperability, not numerical correctness."),
        GateRow("shd_03_single", "single fixture", "diagnostic_validated", "The patched single evaluator is independently checked at T=10^10.", str(exact["single_fixture"]), "The observed error is not an interval bound."),
        GateRow("shd_04_multi", "multi fixture", "diagnostic_validated", "The upstream 0.01 multi grid is checked against direct Riemann-Siegel main sums.", str(exact["accepted_grid"]), "Only the tested low heights are covered."),
        GateRow("shd_05_scale", "height scaling", "diagnostic_validated", "The error decreases between T=10^10 and T=10^12 in the direction predicted by the documented asymptotic floor.", str(exact["scaling"]), "Two heights do not establish an error theorem or its constant."),
        GateRow("shd_06_phase", "phase lock", "diagnostic_validated", "The physical normalizer nearly aligns with the stationary Hardy endpoint direction.", str(exact["phase_lock"]), "The residual endpoint quadrature is retained, not rounded to zero."),
        GateRow("shd_07_kernel", "kernel identity", "proved", "A real combination of corrected Hardy main sums implements the fitted complex carrier kernel.", str(exact["kernel_identity"]), "The identity is conditional on the supplied main-sum values."),
        GateRow("shd_08_endpoint", "endpoint ownership", "diagnostic_validated", "The physical phase lock removes the stationary endpoint quadrature obstruction.", str(exact["terminal_peel"]), "The two starred endpoint half-weights are still adjusted explicitly."),
        GateRow("shd_09_physical", "physical kernel", "diagnostic_validated", "All 24 physical base and tangent kernels pass independent held-out bulk grids.", str(exact["physical_kernel"]), "This checks interpolation only, not Hardy-sample accuracy."),
        GateRow("shd_10_condition", "conditioning", "guard_validated", "The coefficient norms make independent per-sample error bounds unusable at proof scale.", str(exact["conditioning"]), "Empirical correlated cancellation cannot replace a certified joint error bound."),
        GateRow("shd_11_h002", "spacing control", "falsified", "The 0.02 spacing shortcut is rejected.", str(exact["rejected_002"]), "It must not be used for a physical run."),
        GateRow("shd_12_h004", "spacing control", "falsified", "The 0.04 spacing shortcut is rejected.", str(exact["rejected_004"]), "It must not be used for a physical run."),
        GateRow("shd_13_runtime", "runtime guard", "guard_validated", "The physical multi-value job is not launched until it is resumable within the session contract.", str(exact["runtime_guard"]), "A non-resumable run may not be allowed to cross the four-hour turn ceiling."),
        GateRow("shd_14_boundary", "proof boundary", "guard_validated", "This bridge is a numerical diagnostic and not proof evidence for the determinant sign or RH.", str(exact["proof_boundary"]), "No interval enclosure or asymptotic error constant is supplied."),
        GateRow("shd_15_handoff", "next implementation", "open", "The next gate is a resumable 0.01 multi evaluator with an explicit joint-error ledger.", str(exact["next_action"]), "No physical value is claimed before that gate closes."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    summary = payload["physical_kernel_summary"]
    calibration = payload["calibration"]
    lines = [
        "# Newman C1 Shifted-Hardy Diagnostic Bridge Gate",
        "",
        "Date: 2026-08-05",
        "",
        "Status: pinned external evaluator, low-height calibration, and physical kernel interpolation validated; physical carrier values and proof error bounds remain open; not a proof of RH.",
        "",
        "## Structural decision",
        "",
        str(exact["structural_obstruction"]),
        "",
        "The alternate diagnostic bridge uses",
        "",
        "```text",
        str(exact["kernel_identity"]),
        "```",
        "",
        "The common phase is essential. Without the physical normalizer, the stationary endpoint loses one quadrature. The retained mismatch is",
        "",
        "```text",
        str(exact["phase_lock"]),
        "```",
        "",
        "## External calibration",
        "",
        f"- Upstream: `{UPSTREAM_REPOSITORY}`.",
        f"- Commit: `{UPSTREAM_COMMIT}`.",
        f"- License: `{UPSTREAM_LICENSE}`.",
        f"- Single T=10^10 absolute error: `{calibration['single_t1e10']['absolute_error']:.12g}`.",
        f"- Accepted 0.01-grid maximum main error at T=10^10: `{calibration['t1e10_et005']['maximum_absolute']:.12g}`.",
        f"- Accepted 0.01-grid maximum main error at T=10^12: `{calibration['t1e12_et005']['maximum_absolute']:.12g}`.",
        f"- Rejected 0.02-grid maximum main error: `{calibration['t1e10_h002']['maximum_absolute']:.12g}`.",
        f"- Rejected 0.04-grid maximum main error: `{calibration['t1e10_h004']['maximum_absolute']:.12g}`.",
        "",
        "The wider grids are falsification controls. They are not admissible physical samplers.",
        "",
        "## Physical kernel",
        "",
        f"No terminal band peel is required. Across `{summary['roots']}` roots and `{summary['physical_rows']}` base/tangent rows, including both interval endpoints, the maximum held-out relative error is `{summary['maximum_held_out_relative_error']:.12g}`. The maximum coefficient l1 norm is `{summary['maximum_coefficient_l1']:.12g}`.",
        "",
        str(exact["conditioning"]),
        "",
        "## Runtime boundary",
        "",
        str(exact["runtime_guard"]),
        "",
        "## Proof boundary",
        "",
        str(exact["proof_boundary"]),
        "",
        str(exact["next_action"]),
        "",
        payload["success"],
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    set_below_normal_priority()
    mp.mp.dps = 80
    roots_payload, adapter_payload, source_audit = load_and_audit_sources()

    carrier_10 = direct_carrier("1e10")
    samples_10_001 = direct_main_samples(carrier_10, 0.01)
    samples_10_002 = direct_main_samples(carrier_10, 0.02)
    samples_10_004 = direct_main_samples(carrier_10, 0.04)
    carrier_12 = direct_carrier("1e12")
    samples_12_001 = direct_main_samples(carrier_12, 0.01)

    grand_10_et05 = parse_grand_totals(FIXTURE_PATHS["multi_t1e10_et05"])
    grand_10_et005 = parse_grand_totals(FIXTURE_PATHS["multi_t1e10_et005"])
    grand_10_h002 = parse_grand_totals(FIXTURE_PATHS["multi_t1e10_et005_h002"])
    grand_10_h004 = parse_grand_totals(FIXTURE_PATHS["multi_t1e10_et005_h004"])
    grand_12_et005 = parse_grand_totals(FIXTURE_PATHS["multi_t1e12_et005"])

    fast_10_et05 = grand_10_et05 - samples_10_001["correction"]
    fast_10_et005 = grand_10_et005 - samples_10_001["correction"]
    fast_10_h002 = grand_10_h002 - samples_10_002["correction"]
    fast_10_h004 = grand_10_h004 - samples_10_004["correction"]
    fast_12_et005 = grand_12_et005 - samples_12_001["correction"]

    calibration = {
        "t1e10_et05": error_summary(fast_10_et05, samples_10_001["main"]),
        "t1e10_et005": error_summary(fast_10_et005, samples_10_001["main"]),
        "t1e10_h002": error_summary(fast_10_h002, samples_10_002["main"]),
        "t1e10_h004": error_summary(fast_10_h004, samples_10_004["main"]),
        "t1e12_et005": error_summary(fast_12_et005, samples_12_001["main"]),
    }
    single_value = parse_single_total(FIXTURE_PATHS["single_t1e10_et05"])
    exact_single = mp.siegelz(mp.mpf("1e10"))
    calibration["single_t1e10"] = {
        "fast_value": single_value,
        "mpmath_value": str(exact_single),
        "absolute_error": float(abs(exact_single - single_value)),
        "relative_error": float(abs(exact_single - single_value) / abs(exact_single)),
    }

    first_root = roots_payload["precision_ladder"][-1]["root"]
    first_t = mp.mpf(adapter_payload["physical_rows"][0]["chart"]["T_0"])
    first_omega = parse_complex(first_root["eta"])
    if int(first_root["N"]) % 2:
        first_omega = -first_omega
    first_line = complex(mp.e ** (-1j * mp.siegeltheta(first_t)))
    physical_phase_ratio = first_omega / first_line

    low_reconstruction = {
        "t1e10_et05": low_height_reconstruction(
            carrier_10, samples_10_001, fast_10_et05, physical_phase_ratio
        ),
        "t1e10_et005": low_height_reconstruction(
            carrier_10, samples_10_001, fast_10_et005, physical_phase_ratio
        ),
        "t1e12_et005": low_height_reconstruction(
            carrier_12, samples_12_001, fast_12_et005, physical_phase_ratio
        ),
    }

    physical_kernels, physical_summary = physical_kernel_certificate(
        roots_payload, adapter_payload
    )

    require(calibration["single_t1e10"]["absolute_error"] < 0.01, "single fixture")
    require(calibration["t1e10_et005"]["maximum_absolute"] < 0.01, "0.01 low fixture")
    require(calibration["t1e12_et005"]["maximum_absolute"] < 0.005, "0.01 scaling fixture")
    require(calibration["t1e10_h002"]["maximum_absolute"] > 1.0, "0.02 control did not fail")
    require(calibration["t1e10_h004"]["maximum_absolute"] > 3.0, "0.04 control did not fail")
    require(physical_summary["physical_rows"] == 24, "physical kernel row count")
    require(
        physical_summary["maximum_held_out_relative_error"] < 1.0e-8,
        "physical held-out kernel error",
    )
    require(
        physical_summary["maximum_coefficient_l1"] > 1.0e5,
        "conditioning guard unexpectedly disappeared",
    )
    require(
        low_reconstruction["t1e12_et005"]["maximum_fast_main_relative_error"]
        < 0.01,
        "t1e12 reconstructed fixture",
    )

    exact = {
        "structural_obstruction": "The exact ridge identity 2 M_N-F_N=M_N-tau_band+E_band leaves one raw carrier. The diagonal-plus-adjacent collar is locally exact but does not by itself give an O(polylog N) physical evaluator.",
        "single_fixture": f"At T=10^10 the patched single code gives {single_value:.10f}, while mpmath gives {mp.nstr(exact_single, 30)}; absolute error {calibration['single_t1e10']['absolute_error']:.12g}.",
        "accepted_grid": f"For the upstream shift_step=0.01 grid, the maximum direct-main errors are {calibration['t1e10_et005']['maximum_absolute']:.12g} at T=10^10 and {calibration['t1e12_et005']['maximum_absolute']:.12g} at T=10^12.",
        "scaling": f"The accepted-grid mean error magnitude decreases from {abs(calibration['t1e10_et005']['mean']):.12g} to {abs(calibration['t1e12_et005']['mean']):.12g}. This is compatible with, but does not prove, the documented O(T^(-1/4)) floor.",
        "phase_lock": f"At the first physical root omega_c/exp(-i theta_RS(T_0))={physical_phase_ratio.real:.16g}{physical_phase_ratio.imag:+.16g}i, so the retained perpendicular endpoint fraction is {abs(physical_phase_ratio.imag):.12g}.",
        "kernel_identity": "If K_j(lambda)=2 exp[-i theta_RS(T+delta_j)]exp(i delta_j lambda) and g(lambda) is fitted as sum_j c_jK_j(lambda) with real c_j, then sum_j c_j H_main(T+delta_j)=Re sum_(n<=N)n^(-1/2+iT)g(log n). The half-endpoint carrier convention is restored by explicit endpoint subtraction.",
        "terminal_peel": "At the actual frequency T_0=Omega, the parity-adjusted source phase agrees with the stationary Hardy endpoint direction to the recorded precision. No terminal band is peeled: the fit includes lambda=log N, while the starred carrier's half-weights at n=1 and n=N are restored by explicit subtraction.",
        "physical_kernel": f"Three roots times eight base/tangent rows give 24 fitted kernels. The maximum held-out relative bulk error is {physical_summary['maximum_held_out_relative_error']:.12g}.",
        "conditioning": f"The largest physical coefficient l1 norm is {physical_summary['maximum_coefficient_l1']:.12g}. A worst-case independent sample-error estimate is therefore amplified by that factor; only empirical low-height correlation currently reduces it.",
        "rejected_002": f"At T=10^10, shift_step=0.02 has maximum direct-main error {calibration['t1e10_h002']['maximum_absolute']:.12g} and is rejected.",
        "rejected_004": f"At T=10^10, shift_step=0.04 has maximum direct-main error {calibration['t1e10_h004']['maximum_absolute']:.12g} and is rejected.",
        "runtime_guard": "The upstream multi evaluator has no restart facility. Measured one-rank runtime rises from about 9 seconds at T=10^10 to 153 seconds at T=10^12, so an uncheckpointed physical T about 3.26e22 run cannot be assumed to finish inside one four-hour cycle even with four night workers.",
        "proof_boundary": "The external algorithm documents only an asymptotic relative-error order, with no validated constant or interval enclosure here. Kernel interpolation, raw ten-decimal outputs, and correlated low-height errors prove no physical observation, determinant sign, flow inequality, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "next_action": "Add append-only block checkpoints and deterministic resume to the accepted shift_step=0.01 multi evaluator, preserve each of the 15 partial sums and the RS-tail state, and validate stop/resume equivalence at T=10^10 and T=10^12. Then formulate a joint correlated-error diagnostic; do not launch the physical job or call it proof evidence before both gates close.",
    }
    rows = build_rows(exact)
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-05",
        "status": "shifted-Hardy diagnostic bridge calibrated; physical kernel fits retained; wider grids rejected; resumable physical evaluator and proof error bound open",
        "upstream": {
            "repository": UPSTREAM_REPOSITORY,
            "commit": UPSTREAM_COMMIT,
            "license": UPSTREAM_LICENSE,
        },
        "source_audit": source_audit,
        "source_sha256": {
            key: value["sha256"] for key, value in source_audit.items()
        },
        "calibration": calibration,
        "low_height_reconstruction": low_reconstruction,
        "physical_kernels": physical_kernels,
        "physical_kernel_summary": physical_summary,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "source_artifacts": len(source_audit),
            "accepted_shift_steps": [DEFAULT_SHIFT_STEP],
            "rejected_shift_steps": [0.02, 0.04],
            "low_height_main_calibrations": 5,
            "low_height_reconstructions": 3,
            "physical_roots": physical_summary["roots"],
            "physical_kernel_rows": physical_summary["physical_rows"],
            "terminal_peel": TERMINAL_PEEL,
            "physical_values_evaluated": 0,
            "interval_enclosures": 0,
            "open_rows": sum(row.readiness == "open" for row in rows),
        },
        "next_action": exact["next_action"],
        "pi_provenance": "Every pi in the bridge comes from T_0=2pi a^2, the Riemann-Siegel theta normalization, or the Fourier character exp(i delta log n). No circle or polygon is inserted into the arithmetic.",
        "proof_boundary": exact["proof_boundary"],
        "success": "built Newman C1 shifted-Hardy diagnostic bridge gate: 15 rows, 0 issues, 4 portability patches, 6 raw fixtures, 5 direct-main calibrations, 3 low-height reconstructions, 24 physical kernel fits, shift_step 0.01 retained, 0.02 and 0.04 rejected, 0 physical values and 1 open resumable-evaluator row",
    }
    atomic_write(RESULT_PATH, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(payload["success"])


if __name__ == "__main__":
    main()
