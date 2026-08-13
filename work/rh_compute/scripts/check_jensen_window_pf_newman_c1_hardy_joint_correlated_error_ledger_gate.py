#!/usr/bin/env python3
"""Independently check the shifted-Hardy joint correlated-error ledger."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import mpmath as mp
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_hardy_joint_correlated_error_ledger_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
HEIGHT_KEYS = ("t1e10_et005", "t1e12_et005")
RCOND_GRID = (1.0e-8, 1.0e-10, 1.0e-12, 1.0e-14, 1.0e-15, 1.0e-16)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_dot(left: np.ndarray, right: np.ndarray) -> float:
    return math.fsum(float(a) * float(b) for a, b in zip(left, right, strict=True))


def close(left: float, right: float, *, atol: float = 1.0e-12, rtol: float = 1.0e-10) -> bool:
    return abs(left - right) <= atol + rtol * max(abs(left), abs(right))


def parse_complex(text: str) -> complex:
    match = re.fullmatch(r"\(([^ ]+)\s*([+-])\s*([^j]+)j\)", text.strip())
    require(match is not None, f"cannot parse complex value: {text}")
    sign = 1.0 if match.group(2) == "+" else -1.0
    return complex(float(match.group(1)), sign * float(match.group(3)))


def parse_real_part(text: str) -> float:
    match = re.match(r"\(([^ ]+)", text.strip())
    require(match is not None, f"cannot parse real part: {text}")
    return float(match.group(1))


def polynomial(values: list[dict[str, str]]) -> np.ndarray:
    return np.asarray(
        [complex(float(value["re"]), float(value["im"])) for value in values],
        dtype=np.complex128,
    )


def check_sources(payload: dict) -> dict[str, Path]:
    require(len(payload["source_audit"]) == 9, "source audit count")
    paths = {}
    for key, audit in payload["source_audit"].items():
        path = REPO_ROOT / audit["path"]
        require(path.is_file(), f"missing source {key}")
        require(file_hash(path) == audit["sha256"], f"source hash drift {key}")
        paths[key] = path
    return paths


def independent_polynomial_tails(error: np.ndarray) -> dict[str, float]:
    x = np.linspace(-1.0, 1.0, 15)
    norm = float(np.linalg.norm(error))
    tails = {}
    for degree in range(7):
        design = np.vander(x, degree + 1, increasing=True)
        fit = np.linalg.lstsq(design, error, rcond=None)[0]
        residual = error - design @ fit
        tails[str(degree)] = float(np.linalg.norm(residual)) / norm
    return tails


def check_residual_modes(payload: dict, bridge: dict) -> dict[str, np.ndarray]:
    errors = {}
    for key in HEIGHT_KEYS:
        source = np.asarray(bridge["calibration"][key]["errors"], dtype=float)
        stored = payload["residual_modes"][key]
        require(np.array_equal(source, np.asarray(stored["values"])), f"residual values {key}")
        require(close(float(np.mean(source)), stored["mean"]), f"residual mean {key}")
        require(close(float(np.std(source)), stored["standard_deviation"]), f"residual std {key}")
        require(close(float(np.linalg.norm(source)), stored["l2_norm"]), f"residual norm {key}")
        tails = independent_polynomial_tails(source)
        for degree, value in tails.items():
            require(
                close(value, stored["tail_energy_fraction_after_degree"][degree], atol=1.0e-11),
                f"polynomial tail {key}:{degree}",
            )
        errors[key] = source

    low = errors[HEIGHT_KEYS[0]]
    high = errors[HEIGHT_KEYS[1]]
    low_norm = float(np.linalg.norm(low))
    high_norm = float(np.linalg.norm(high))
    cosine = stable_dot(low, high) / (low_norm * high_norm)
    scale = stable_dot(low, high) / stable_dot(low, low)
    relative_residual = float(np.linalg.norm(high - scale * low)) / high_norm
    cross = payload["cross_height"]
    require(close(cosine, cross["cosine_similarity"]), "cross-height cosine")
    require(close(scale, cross["best_scalar_t1e10_to_t1e12"]), "cross-height scale")
    require(close(relative_residual, cross["relative_residual_after_best_scalar"]), "cross-height residual")
    return errors


def check_saved_projections(
    payload: dict, bridge: dict, errors: dict[str, np.ndarray]
) -> None:
    source_rows = {
        (int(root["root_index"]), row["family"], row["name"]): row
        for root in bridge["physical_kernels"]
        for row in root["rows"]
    }
    stored_rows = payload["physical_projections"]
    require(len(source_rows) == len(stored_rows) == 24, "physical row count")
    increased = 0
    constant_annihilating = 0
    ratios = []
    for stored in stored_rows:
        key = (stored["root_index"], stored["family"], stored["name"])
        require(key in source_rows, f"unknown physical row {key}")
        source = source_rows[key]
        coefficients = np.asarray(source["real_hardy_coefficients"], dtype=float)
        l1 = math.fsum(abs(float(value)) for value in coefficients)
        coefficient_sum = math.fsum(float(value) for value in coefficients)
        require(close(l1, stored["coefficient_l1"], atol=1.0e-9), f"l1 {key}")
        require(close(coefficient_sum, stored["coefficient_sum"], atol=1.0e-12), f"sum {key}")
        absolute = {}
        for height in HEIGHT_KEYS:
            total = stable_dot(coefficients, errors[height])
            projection = stored["projections"][height]
            require(close(total, projection["projection"], atol=1.0e-12), f"projection {key}:{height}")
            independent = l1 * float(np.max(np.abs(errors[height])))
            require(close(independent, projection["independent_l1_bound"], atol=1.0e-9), f"bound {key}:{height}")
            common = float(np.mean(errors[height])) * coefficient_sum
            require(close(common, projection["constant_mode_projection"], atol=1.0e-12), f"common {key}:{height}")
            require(close(total - common, projection["centered_projection"], atol=1.0e-12), f"centered {key}:{height}")
            absolute[height] = abs(total)
        ratio = absolute[HEIGHT_KEYS[1]] / max(absolute[HEIGHT_KEYS[0]], 1.0e-300)
        require(close(ratio, stored["absolute_projection_height_ratio"], atol=1.0e-10), f"ratio {key}")
        ratios.append(ratio)
        increased += ratio > 1.0
        constant_annihilating += abs(coefficient_sum) / l1 < 1.0e-12

    summary = payload["projection_summary"]
    require(increased == summary["increased_projection_rows"] == 21, "increased projections")
    require(
        constant_annihilating == summary["constant_annihilating_rows"] == 21,
        "constant annihilation",
    )
    require(close(max(ratios), summary["maximum_absolute_projection_height_ratio"]), "maximum ratio")


def replay_rcond_ensemble(
    payload: dict,
    roots_payload: dict,
    adapter_payload: dict,
    errors: dict[str, np.ndarray],
) -> None:
    roots = [roots_payload["precision_ladder"][-1], *roots_payload["additional_roots"]]
    physical_rows = adapter_payload["physical_rows"]
    require(len(roots) == len(physical_rows) == 3, "root/adapter alignment")
    shifts = np.arange(-7, 8, dtype=float) * 0.01
    stored_ensemble = payload["rcond_ensemble"]
    require(len(stored_ensemble) == len(RCOND_GRID), "rcond row count")

    for rcond, stored in zip(RCOND_GRID, stored_ensemble, strict=True):
        maximum_error = 0.0
        maximum_l1 = 0.0
        projections = {key: [] for key in HEIGHT_KEYS}
        increased = 0
        for root_row, physical in zip(roots, physical_rows, strict=True):
            root = root_row["root"]
            t_value = mp.mpf(physical["chart"]["T_0"])
            n_value = int(physical["N"])
            heat_time = float(physical["chart"]["t"])
            sigma = parse_real_part(root["s_star"])
            omega = parse_complex(root["eta"])
            if n_value % 2:
                omega = -omega
            phases = np.asarray(
                [
                    complex(
                        mp.e
                        ** (
                            -1j
                            * mp.siegeltheta(
                                t_value + mp.mpf(str(float(delta)))
                            )
                        )
                    )
                    for delta in shifts
                ]
            )
            upper = math.log(n_value)
            train = np.linspace(0.0, upper, 2401)
            midpoint = (np.arange(4800) + 0.5) * upper / 4800
            chebyshev = 0.5 * upper * (
                1.0 - np.cos(np.linspace(0.0, np.pi, 2001))
            )
            test = np.unique(np.concatenate((midpoint, chebyshev, [0.0, upper])))
            train_weight = np.exp(
                heat_time * train**2 / 4.0 - (sigma - 0.5) * train
            )
            test_weight = np.exp(
                heat_time * test**2 / 4.0 - (sigma - 0.5) * test
            )
            train_matrix = 2.0 * phases[None, :] * np.exp(
                1j * train[:, None] * shifts[None, :]
            )
            real_matrix = np.vstack((train_matrix.real, train_matrix.imag))
            test_matrix = 2.0 * phases[None, :] * np.exp(
                1j * test[:, None] * shifts[None, :]
            )

            for family in ("base_rows", "tangent_rows"):
                for values in physical["polynomials"][family].values():
                    source_polynomial = polynomial(values)
                    target = (
                        omega
                        * train_weight
                        * np.polynomial.polynomial.polyval(train, source_polynomial)
                    )
                    real_target = np.concatenate((target.real, target.imag))
                    coefficients = np.linalg.lstsq(
                        real_matrix, real_target, rcond=rcond
                    )[0]
                    truth = (
                        omega
                        * test_weight
                        * np.polynomial.polynomial.polyval(test, source_polynomial)
                    )
                    prediction = test_matrix @ coefficients
                    scale = max(1.0, float(np.max(np.abs(truth))))
                    maximum_error = max(
                        maximum_error,
                        float(np.max(np.abs(prediction - truth)) / scale),
                    )
                    maximum_l1 = max(
                        maximum_l1,
                        math.fsum(abs(float(value)) for value in coefficients),
                    )
                    values_by_height = {
                        key: abs(stable_dot(coefficients, errors[key]))
                        for key in HEIGHT_KEYS
                    }
                    for key in HEIGHT_KEYS:
                        projections[key].append(values_by_height[key])
                    increased += (
                        values_by_height[HEIGHT_KEYS[1]]
                        > values_by_height[HEIGHT_KEYS[0]]
                    )

        require(close(rcond, stored["rcond"], atol=0.0), f"rcond value {rcond}")
        require(close(maximum_error, stored["maximum_held_out_relative_error"], atol=1.0e-12), f"held error {rcond}")
        require(close(maximum_l1, stored["maximum_coefficient_l1"], atol=1.0e-5), f"l1 replay {rcond}")
        require(increased == stored["increased_projection_rows"] == 21, f"increase replay {rcond}")
        for key in HEIGHT_KEYS:
            require(
                close(
                    max(projections[key]),
                    stored["maximum_absolute_projection"][key],
                    atol=1.0e-10,
                ),
                f"projection replay {rcond}:{key}",
            )


def main() -> None:
    mp.mp.dps = 80
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    paths = check_sources(payload)
    bridge = json.loads(paths["shifted_hardy_result"].read_text(encoding="utf-8"))
    roots = json.loads(paths["calibrated_roots"].read_text(encoding="utf-8"))
    adapter = json.loads(paths["physical_adapter"].read_text(encoding="utf-8"))
    errors = check_residual_modes(payload, bridge)
    check_saved_projections(payload, bridge, errors)
    replay_rcond_ensemble(payload, roots, adapter, errors)

    require(payload["physical_values_evaluated"] == 0, "physical values")
    require(payload["interval_enclosures"] == 0, "intervals")
    require("Retire scalar" in payload["route_decision"], "route decision")
    require("mode" in payload["missing_theorem"].lower(), "missing theorem")
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "success marker")
    require("21/24" in note, "projection obstruction marker")
    require("not a proof" in note.lower(), "proof boundary marker")

    print(
        "validated Newman C1 Hardy joint correlated-error ledger gate: "
        "2 residual vectors, 15 shifts, 24 physical rows, "
        "21 constant-annihilating rows, 21/24 projected errors increase, "
        "6 cutoff replays, scalar extrapolation retired, "
        "0 physical values and 1 open modewise-error theorem"
    )


if __name__ == "__main__":
    main()
