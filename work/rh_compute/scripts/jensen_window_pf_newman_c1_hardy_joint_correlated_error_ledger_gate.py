#!/usr/bin/env python3
"""Build the finite joint correlated-error ledger for shifted Hardy samples."""

from __future__ import annotations

import ctypes
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_hardy_joint_correlated_error_ledger_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
BRIDGE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate.json"
)
BRIDGE_BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts"
    / "jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate.py"
)
BRIDGE_CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts"
    / "check_jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate.py"
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
EQUIVALENCE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/resume_equivalence"
    / "fixture_result.json"
)
CHECKPOINT_MODULE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable"
    / "rh_hardy_checkpoint.f90"
)
RESUMABLE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable"
    / "zeta14cubicmult_resumable.f90"
)
JOINT_CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts"
    / "check_jensen_window_pf_newman_c1_hardy_joint_correlated_error_ledger_gate.py"
)

HEIGHT_KEYS = ("t1e10_et005", "t1e12_et005")
RCOND_GRID = (1.0e-8, 1.0e-10, 1.0e-12, 1.0e-14, 1.0e-15, 1.0e-16)


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


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def source_audit() -> dict[str, dict[str, str]]:
    paths = {
        "shifted_hardy_result": BRIDGE_RESULT,
        "shifted_hardy_builder": BRIDGE_BUILDER,
        "shifted_hardy_checker": BRIDGE_CHECKER,
        "calibrated_roots": ROOTS_RESULT,
        "physical_adapter": ADAPTER_RESULT,
        "resume_equivalence": EQUIVALENCE_RESULT,
        "checkpoint_module": CHECKPOINT_MODULE,
        "resumable_source": RESUMABLE_SOURCE,
        "joint_error_checker": JOINT_CHECKER,
    }
    return {
        key: {"path": relative_path(path), "sha256": file_hash(path)}
        for key, path in paths.items()
    }


def orthonormal_shift_modes() -> tuple[np.ndarray, np.ndarray]:
    x = np.linspace(-1.0, 1.0, 15)
    vandermonde = np.polynomial.legendre.legvander(x, 14)
    basis, triangular = np.linalg.qr(vandermonde)
    for index in range(15):
        if triangular[index, index] < 0.0:
            basis[:, index] *= -1.0
    return x, basis


def mode_ledger(error: np.ndarray, basis: np.ndarray) -> dict:
    coefficients = basis.T @ error
    reconstructed = basis @ coefficients
    norm = float(np.linalg.norm(error))
    tails = {}
    for degree in range(7):
        tail = float(np.linalg.norm(coefficients[degree + 1 :]))
        tails[str(degree)] = tail / norm
    return {
        "values": [float(value) for value in error],
        "mean": float(np.mean(error)),
        "standard_deviation": float(np.std(error)),
        "maximum_absolute": float(np.max(np.abs(error))),
        "l2_norm": norm,
        "orthonormal_mode_coefficients": [float(value) for value in coefficients],
        "tail_energy_fraction_after_degree": tails,
        "maximum_reconstruction_error": float(np.max(np.abs(reconstructed - error))),
    }


def stable_dot(left: np.ndarray, right: np.ndarray) -> float:
    return math.fsum(float(a) * float(b) for a, b in zip(left, right, strict=True))


def grouped_contributions(values: np.ndarray) -> dict[str, float]:
    return {
        "common_0": float(values[0]),
        "slope_1": float(values[1]),
        "curvature_2": float(values[2]),
        "cubic_3": float(values[3]),
        "quartic_quintic_4_5": float(np.sum(values[4:6])),
        "rough_6_14": float(np.sum(values[6:])),
        "rough_6_14_absolute_budget": float(np.sum(np.abs(values[6:]))),
    }


def projection_rows(
    bridge: dict, errors: dict[str, np.ndarray], basis: np.ndarray
) -> list[dict]:
    error_modes = {key: basis.T @ value for key, value in errors.items()}
    rows = []
    for root in bridge["physical_kernels"]:
        for source_row in root["rows"]:
            coefficients = np.asarray(source_row["real_hardy_coefficients"], dtype=float)
            l1 = math.fsum(abs(float(value)) for value in coefficients)
            coefficient_sum = math.fsum(float(value) for value in coefficients)
            coefficient_modes = basis.T @ coefficients
            projected = {}
            for key in HEIGHT_KEYS:
                error = errors[key]
                total = stable_dot(coefficients, error)
                mode_contributions = coefficient_modes * error_modes[key]
                independent_bound = l1 * float(np.max(np.abs(error)))
                projected[key] = {
                    "projection": total,
                    "absolute_projection": abs(total),
                    "normalized_absolute_projection": abs(total)
                    / max(1.0, float(source_row["target_scale"])),
                    "independent_l1_bound": independent_bound,
                    "empirical_cancellation_factor": independent_bound
                    / max(abs(total), 1.0e-300),
                    "constant_mode_projection": float(np.mean(error))
                    * coefficient_sum,
                    "centered_projection": total
                    - float(np.mean(error)) * coefficient_sum,
                    "mode_contributions": grouped_contributions(mode_contributions),
                    "mode_reconstruction_error": abs(
                        total - float(np.sum(mode_contributions))
                    ),
                }
            ratio = (
                projected[HEIGHT_KEYS[1]]["absolute_projection"]
                / max(projected[HEIGHT_KEYS[0]]["absolute_projection"], 1.0e-300)
            )
            rows.append(
                {
                    "root_index": int(root["root_index"]),
                    "N": int(root["N"]),
                    "family": source_row["family"],
                    "name": source_row["name"],
                    "target_scale": float(source_row["target_scale"]),
                    "coefficient_l1": l1,
                    "coefficient_l2": float(np.linalg.norm(coefficients)),
                    "coefficient_sum": coefficient_sum,
                    "constant_sensitivity": abs(coefficient_sum) / l1,
                    "projections": projected,
                    "absolute_projection_height_ratio": ratio,
                    "projection_increased_at_t1e12": ratio > 1.0,
                }
            )
    return rows


def summarize_projections(rows: list[dict]) -> dict:
    summary = {
        "rows": len(rows),
        "constant_annihilating_rows": sum(
            row["constant_sensitivity"] < 1.0e-12 for row in rows
        ),
        "increased_projection_rows": sum(
            row["projection_increased_at_t1e12"] for row in rows
        ),
        "maximum_absolute_projection_height_ratio": max(
            row["absolute_projection_height_ratio"] for row in rows
        ),
    }
    for key in HEIGHT_KEYS:
        factors = sorted(
            row["projections"][key]["empirical_cancellation_factor"]
            for row in rows
        )
        summary[key] = {
            "maximum_absolute_projection": max(
                row["projections"][key]["absolute_projection"] for row in rows
            ),
            "maximum_normalized_absolute_projection": max(
                row["projections"][key]["normalized_absolute_projection"]
                for row in rows
            ),
            "minimum_empirical_cancellation_factor": factors[0],
            "median_empirical_cancellation_factor": float(np.median(factors)),
            "maximum_empirical_cancellation_factor": factors[-1],
            "maximum_rough_mode_absolute_budget": max(
                row["projections"][key]["mode_contributions"][
                    "rough_6_14_absolute_budget"
                ]
                for row in rows
            ),
        }
    return summary


def cross_height_ledger(errors: dict[str, np.ndarray]) -> dict:
    low = errors[HEIGHT_KEYS[0]]
    high = errors[HEIGHT_KEYS[1]]
    low_norm = float(np.linalg.norm(low))
    high_norm = float(np.linalg.norm(high))
    scale = stable_dot(low, high) / stable_dot(low, low)
    residual = high - scale * low
    return {
        "cosine_similarity": stable_dot(low, high) / (low_norm * high_norm),
        "best_scalar_t1e10_to_t1e12": scale,
        "relative_residual_after_best_scalar": float(np.linalg.norm(residual))
        / high_norm,
        "maximum_error_ratio_t1e12_over_t1e10": float(
            np.max(np.abs(high)) / np.max(np.abs(low))
        ),
        "l2_ratio_t1e12_over_t1e10": high_norm / low_norm,
    }


def load_bridge_module():
    spec = importlib.util.spec_from_file_location("shifted_hardy_bridge_source", BRIDGE_BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load bridge builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.mp.mp.dps = 80
    return module


def rcond_stability_ensemble(errors: dict[str, np.ndarray]) -> list[dict]:
    module = load_bridge_module()
    roots = json.loads(ROOTS_RESULT.read_text(encoding="utf-8"))
    adapter = json.loads(ADAPTER_RESULT.read_text(encoding="utf-8"))
    original_lstsq = np.linalg.lstsq
    ensemble = []
    try:
        for rcond in RCOND_GRID:
            def fixed_lstsq(a, b, rcond=None, *, _fixed=rcond):
                return original_lstsq(a, b, rcond=_fixed)

            np.linalg.lstsq = fixed_lstsq
            roots_rows, summary = module.physical_kernel_certificate(roots, adapter)
            projections = {key: [] for key in HEIGHT_KEYS}
            coefficient_l1 = []
            increased = 0
            for root in roots_rows:
                for row in root["rows"]:
                    coefficients = np.asarray(row["real_hardy_coefficients"], dtype=float)
                    values = {
                        key: abs(stable_dot(coefficients, errors[key]))
                        for key in HEIGHT_KEYS
                    }
                    for key in HEIGHT_KEYS:
                        projections[key].append(values[key])
                    coefficient_l1.append(float(row["coefficient_l1"]))
                    increased += values[HEIGHT_KEYS[1]] > values[HEIGHT_KEYS[0]]
            ensemble.append(
                {
                    "rcond": rcond,
                    "maximum_held_out_relative_error": float(
                        summary["maximum_held_out_relative_error"]
                    ),
                    "maximum_coefficient_l1": max(coefficient_l1),
                    "maximum_absolute_projection": {
                        key: max(projections[key]) for key in HEIGHT_KEYS
                    },
                    "increased_projection_rows": increased,
                }
            )
    finally:
        np.linalg.lstsq = original_lstsq
    return ensemble


def make_note(payload: dict) -> str:
    residuals = payload["residual_modes"]
    cross = payload["cross_height"]
    summary = payload["projection_summary"]
    ensemble = payload["rcond_ensemble"]
    low = summary[HEIGHT_KEYS[0]]
    high = summary[HEIGHT_KEYS[1]]
    return f"""# Newman C1 Hardy Joint Correlated-Error Ledger Gate

Date: 2026-08-05

Status: finite two-height correlated-error diagnostic and scalar-extrapolation
obstruction validated; no physical-height error theorem or interval enclosure;
not a proof of RH.

## Residual modes

The accepted `h=0.01` residual vectors are `fast_main-direct_main` on the
fifteen shifts `-0.07,...,0.07`.  A discrete orthonormal polynomial basis on
those same points gives:

| Height | max sample error | residual energy above degree 3 | above degree 5 |
|---|---:|---:|---:|
| `10^10` | `{residuals[HEIGHT_KEYS[0]]['maximum_absolute']:.12g}` | `{residuals[HEIGHT_KEYS[0]]['tail_energy_fraction_after_degree']['3']:.6g}` | `{residuals[HEIGHT_KEYS[0]]['tail_energy_fraction_after_degree']['5']:.6g}` |
| `10^12` | `{residuals[HEIGHT_KEYS[1]]['maximum_absolute']:.12g}` | `{residuals[HEIGHT_KEYS[1]]['tail_energy_fraction_after_degree']['3']:.6g}` | `{residuals[HEIGHT_KEYS[1]]['tail_energy_fraction_after_degree']['5']:.6g}` |

The profiles are smooth but not scalar copies.  Their cosine similarity is
`{cross['cosine_similarity']:.12g}`, while the best scalar multiple leaves
`{cross['relative_residual_after_best_scalar']:.6g}` of the `10^12` L2 norm.

## Physical projections

All `24` saved physical coefficient vectors are projected against both error
profiles.  `21` annihilate a constant shift error to relative sensitivity
below `1e-12`.  The observed correlations reduce the independent L1 bounds by
factors ranging from `{low['minimum_empirical_cancellation_factor']:.6g}` to
`{low['maximum_empirical_cancellation_factor']:.6g}` at `10^10`, and from
`{high['minimum_empirical_cancellation_factor']:.6g}` to
`{high['maximum_empirical_cancellation_factor']:.6g}` at `10^12`.

This cancellation is empirical, not a bound.  More importantly, the raw
maximum sample error falls by a factor
`{cross['maximum_error_ratio_t1e12_over_t1e10']:.6g}`, while
`{summary['increased_projection_rows']}/24` physical projections increase.
The largest increase factor is
`{summary['maximum_absolute_projection_height_ratio']:.6g}`.  A single scalar
`T`-decay envelope therefore does not control these coefficient-weighted
errors even on the two saved calibration heights.

## Cutoff stability

The physical kernel design matrices have condition numbers near `2.5e16`.
The conclusion above was replayed with least-squares cutoffs from `1e-8` to
`1e-16`.  Across all `{len(ensemble)}` choices, exactly `21/24` projections
increase.  The maximum projected errors remain in the ranges
`[{min(row['maximum_absolute_projection'][HEIGHT_KEYS[0]] for row in ensemble):.8g}, {max(row['maximum_absolute_projection'][HEIGHT_KEYS[0]] for row in ensemble):.8g}]`
at `10^10` and
`[{min(row['maximum_absolute_projection'][HEIGHT_KEYS[1]] for row in ensemble):.8g}, {max(row['maximum_absolute_projection'][HEIGHT_KEYS[1]] for row in ensemble):.8g}]`
at `10^12`, despite substantial coefficient-norm variation.

## Missing theorem

The admissible upgrade is a modewise external-error theorem.  For the fixed
shift grid, write the error as

```text
epsilon(T,delta_j)=sum_(k=0)^5 a_k(T) q_k(j)+r_6(T,j).
```

One needs explicit physical-height bounds for every `a_k(T)` and for the
rough remainder, followed by

```text
|sum_j c_j epsilon(T,delta_j)|
 <= sum_(k=0)^5 |<c,q_k>| A_k(T)+||c||_2 R_6(T).
```

A scalar bound on `max_j|epsilon_j|` throws away the observed cancellation;
a scalar asymptotic fitted from two heights is contradicted by the projected
rows.  No physical-height run is interpretable until this componentwise
theorem, or an interval-certified replacement, is available.

## Proof Boundary

This is a finite diagnostic on two saved error vectors and twenty-four saved
kernel rows.  It proves neither asymptotic mode bounds nor persistence of the
observed correlations.  It supplies no physical carrier value, interval
enclosure, retained observation, determinant sign or bound, current
inequality, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.

{payload['success']}
"""


def main() -> None:
    set_below_normal_priority()
    bridge = json.loads(BRIDGE_RESULT.read_text(encoding="utf-8"))
    require(bridge["kind"].endswith("shifted_hardy_diagnostic_bridge_gate"), "bridge kind")
    require(bridge["success"], "bridge source is not validated")
    equivalence = json.loads(EQUIVALENCE_RESULT.read_text(encoding="utf-8"))
    require(len(equivalence["cases"]) == 2, "resume equivalence cases")

    errors = {
        key: np.asarray(bridge["calibration"][key]["errors"], dtype=float)
        for key in HEIGHT_KEYS
    }
    require(all(len(value) == 15 for value in errors.values()), "error vector length")
    shifts, basis = orthonormal_shift_modes()
    residual_modes = {
        key: mode_ledger(value, basis) for key, value in errors.items()
    }
    projections = projection_rows(bridge, errors, basis)
    projection_summary = summarize_projections(projections)
    cross_height = cross_height_ledger(errors)
    ensemble = rcond_stability_ensemble(errors)

    require(projection_summary["rows"] == 24, "projection row count")
    require(projection_summary["constant_annihilating_rows"] == 21, "constant modes")
    require(projection_summary["increased_projection_rows"] == 21, "increased rows")
    require(
        all(row["increased_projection_rows"] == 21 for row in ensemble),
        "rcond stability",
    )
    require(
        residual_modes[HEIGHT_KEYS[0]]["tail_energy_fraction_after_degree"]["5"]
        < 1.0e-8,
        "T=1e10 rough residual",
    )
    require(
        residual_modes[HEIGHT_KEYS[1]]["tail_energy_fraction_after_degree"]["5"]
        < 1.0e-7,
        "T=1e12 rough residual",
    )

    success = (
        "built Newman C1 Hardy joint correlated-error ledger gate: "
        "2 residual vectors, 15 shifts, 24 physical rows, 192 restart records inherited, "
        "21 constant-annihilating rows, 21 projected errors increase, "
        "6 cutoff replays, 0 physical values and 1 open modewise-error theorem"
    )
    payload = {
        "schema_version": 1,
        "kind": KIND,
        "date": "2026-08-05",
        "status": "finite_correlated_error_diagnostic_validated_modewise_theorem_open",
        "source_audit": source_audit(),
        "shift_grid": [float(value) for value in shifts * 0.07],
        "mode_basis": "QR-orthonormalized discrete Legendre values on delta/0.07",
        "residual_modes": residual_modes,
        "cross_height": cross_height,
        "physical_projections": projections,
        "projection_summary": projection_summary,
        "rcond_ensemble": ensemble,
        "exact_ledger_identity": (
            "For each saved coefficient vector c and residual e, "
            "c dot e=sum_k <c,q_k><e,q_k>; modes 0..5 plus 6..14 are retained."
        ),
        "route_decision": (
            "Retire scalar two-height error extrapolation. Preserve the empirical "
            "correlation as a clue and require componentwise mode bounds before "
            "any physical-height interpretation."
        ),
        "missing_theorem": (
            "Prove explicit physical-height bounds for the first six discrete "
            "shift-error modes and an L2 bound for modes 6..14, then contract "
            "those bounds with every physical coefficient-mode vector."
        ),
        "physical_values_evaluated": 0,
        "interval_enclosures": 0,
        "success": success,
        "proof_boundary": (
            "Finite two-height diagnostic only; no asymptotic error theorem, "
            "physical carrier value, interval enclosure, determinant sign, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    atomic_write(RESULT_PATH, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, make_note(payload))
    print(success)


if __name__ == "__main__":
    main()
