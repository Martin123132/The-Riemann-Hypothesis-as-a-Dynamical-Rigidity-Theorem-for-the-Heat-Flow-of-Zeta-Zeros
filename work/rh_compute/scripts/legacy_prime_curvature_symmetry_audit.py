#!/usr/bin/env python3
"""Audit the legacy prime-curvature plots and their RH relevance."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
import os
from pathlib import Path
from typing import Any

# Keep this exploratory audit inside the daytime one-worker policy.
for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import binary_erosion, gaussian_filter, sobel


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "legacy_prime_curvature_symmetry_audit"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DEFAULT_FIGURE = REPO_ROOT / "outputs" / f"{STEM}_controls.png"

LEGACY_ROOT = Path(
    "C:/Users/ollet/Documents/Codex/2026-07-06/"
    "so-your-job-go-through-my-2/work/corpus_v2_repos/"
    "Motion-TimeSpace-/archive/legacy-pre-formalization-2026-06"
)
LEGACY_PRIME_SOURCE = (
    LEGACY_ROOT
    / "mathematics/prime-curvature/3d-prime-curvature-landscape.md"
)
LEGACY_RH_SOURCE = (
    LEGACY_ROOT
    / "mathematics/riemann-zeta/"
    "the-riemann-hypothesis-as-a-geometric-invariance-principle.md"
)
ATTACHMENT_ROOT = Path(
    "C:/Users/ollet/.codex/codex-remote-attachments/"
    "019f28d7-de19-71a1-acad-d30d74708a99/"
    "CB6DD96B-F37E-493C-9F99-354B77C355E8"
)

CURRENT_SOURCES = {
    "planar_abel": (
        REPO_ROOT / "outputs/jensen_window_pf_mertens_planar_abel_handoff.md"
    ),
    "edge_gram": (
        REPO_ROOT
        / "outputs/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.md"
    ),
    "abel_contact": (
        REPO_ROOT
        / "outputs/jensen_window_pf_newman_polymath15_first_order_centered_"
        "abel_scalar_shear_flux_reduction.md"
    ),
}

SOURCE_MARKERS = {
    LEGACY_PRIME_SOURCE: (
        "r = int(np.sqrt(s))",
        "err = abs(s - r*r)",
        "method='cubic'",
        "grid_z_smooth = gaussian_filter(grid_z, sigma=6)",
        "mask = grad_mag > np.percentile(grad_mag, 85)",
        "eroded = binary_erosion(mask, iterations=2)",
        "root = mp.findroot(lambda t: mp.zeta(0.5 + 1j*t), g)",
        "t_vals = np.array([float(z.imag) for z in zeros])",
        "rz_x = scaled",
        "rz_y = scaled",
    ),
    LEGACY_RH_SOURCE: (
        "analyticity",
        "off-line points are attracted",
        "A fixed point of a contraction",
    ),
}

CURRENT_MARKERS = {
    "planar_abel": (
        "E_(alpha,K)",
        "|O_(alpha,K)|^2",
        "Neither estimate is proved",
    ),
    "edge_gram": (
        "G_K is positive semidefinite",
        "`A_t` is not positive semidefinite",
        "determinant `-1/4`",
        "off-diagonal edge-pair form",
    ),
    "abel_contact": (
        "|mathsf_X|<=delta_L",
        "|mathcal_C_N|>A_L+epsilon_term",
        "independent open theorems",
    ),
}

ZERO_GUESSES = [
    14.134725,
    21.02204,
    25.01085,
    30.4249,
    32.9351,
    37.5862,
    40.9187,
    43.3271,
    48.0051,
    49.7738,
    52.9703,
    56.4462,
    59.347,
    60.8317,
    65.1125,
    67.0798,
    69.5464,
    72.0672,
    75.7047,
    77.1448,
]


@dataclass(frozen=True)
class AuditRow:
    id: str
    role: str
    readiness: str
    statement: str
    boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rounded(value: float, digits: int = 15) -> float:
    return float(f"{value:.{digits}g}")


def source_audit() -> dict[str, Any]:
    missing = [
        str(path)
        for path in (*SOURCE_MARKERS, *CURRENT_SOURCES.values())
        if not path.is_file()
    ]
    if missing:
        raise FileNotFoundError("missing audit sources: " + ", ".join(missing))

    marker_lines: dict[str, dict[str, list[int]]] = {}
    for path, markers in SOURCE_MARKERS.items():
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        found: dict[str, list[int]] = {}
        for marker in markers:
            positions = [
                index + 1 for index, line in enumerate(lines) if marker in line
            ]
            if not positions:
                raise RuntimeError(f"legacy source marker missing: {marker}")
            found[marker] = positions
        marker_lines[path.name] = found

    current_hashes: dict[str, str] = {}
    for key, path in CURRENT_SOURCES.items():
        text = path.read_text(encoding="utf-8")
        for marker in CURRENT_MARKERS[key]:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
        current_hashes[key] = file_hash(path)

    return {
        "external_read_only": True,
        "github_fetch_required": False,
        "legacy_sources": {
            "prime_curvature": {
                "path": str(LEGACY_PRIME_SOURCE),
                "sha256": file_hash(LEGACY_PRIME_SOURCE),
            },
            "geometric_invariance": {
                "path": str(LEGACY_RH_SOURCE),
                "sha256": file_hash(LEGACY_RH_SOURCE),
            },
        },
        "legacy_marker_lines": marker_lines,
        "current_rh_source_sha256": current_hashes,
    }


def attachment_audit() -> dict[str, Any]:
    paths = sorted(ATTACHMENT_ROOT.glob("*.jpg"))
    records = [
        {
            "name": path.name,
            "bytes": path.stat().st_size,
            "sha256": file_hash(path),
        }
        for path in paths
    ]
    groups: dict[str, list[str]] = {}
    for record in records:
        groups.setdefault(record["sha256"], []).append(record["name"])
    duplicate_groups = [
        names for names in groups.values() if len(names) > 1
    ]
    return {
        "count": len(records),
        "unique_hashes": len(groups),
        "records": records,
        "duplicate_groups": duplicate_groups,
    }


def primes_below(limit: int) -> np.ndarray:
    sieve = np.ones(limit, dtype=bool)
    sieve[:2] = False
    for prime in range(2, math.isqrt(limit - 1) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit : prime] = False
    return np.flatnonzero(sieve)


def defect_data(values: np.ndarray) -> dict[str, np.ndarray]:
    p_grid, q_grid = np.meshgrid(values, values, indexing="ij")
    square_sum = p_grid * p_grid + q_grid * q_grid
    root = np.fromiter(
        (math.isqrt(int(value)) for value in square_sum.ravel()),
        dtype=np.int64,
        count=square_sum.size,
    ).reshape(square_sum.shape)
    defect = square_sum - root * root
    radius = np.sqrt(square_sum.astype(float))
    return {
        "p": p_grid,
        "q": q_grid,
        "square_sum": square_sum,
        "root": root,
        "defect": defect,
        "radius": radius,
    }


def eigen_record(matrix: np.ndarray) -> dict[str, Any]:
    eigenvalues = np.linalg.eigvalsh(np.asarray(matrix, dtype=float))
    tolerance = max(1.0, float(np.max(np.abs(eigenvalues)))) * 1e-10
    return {
        "minimum": rounded(float(eigenvalues[0])),
        "maximum": rounded(float(eigenvalues[-1])),
        "negative_count": int(np.count_nonzero(eigenvalues < -tolerance)),
        "zero_count": int(np.count_nonzero(np.abs(eigenvalues) <= tolerance)),
        "positive_count": int(np.count_nonzero(eigenvalues > tolerance)),
        "values": [rounded(float(value), 12) for value in eigenvalues],
    }


def transpose_record(matrix: np.ndarray) -> dict[str, float]:
    difference = np.asarray(matrix) - np.asarray(matrix).T
    field_rms = float(np.sqrt(np.mean(np.asarray(matrix) ** 2)))
    rms = float(np.sqrt(np.mean(difference**2)))
    return {
        "max_abs": rounded(float(np.max(np.abs(difference)))),
        "rms": rounded(rms),
        "relative_rms": rounded(rms / field_rms if field_rms else 0.0),
    }


def interpolate_source(
    data: dict[str, np.ndarray],
    threshold: int,
    size: int,
) -> dict[str, np.ndarray]:
    selected = data["defect"] <= threshold
    points = np.column_stack((data["p"][selected], data["q"][selected]))
    values = data["defect"][selected].astype(float)
    axis = np.linspace(
        float(data["p"].min()), float(data["p"].max()), size
    )
    grid_x, grid_y = np.meshgrid(axis, axis, indexing="ij")
    cubic = griddata(points, values, (grid_x, grid_y), method="cubic")
    nearest = griddata(points, values, (grid_x, grid_y), method="nearest")
    return {
        "axis": axis,
        "grid_x": grid_x,
        "grid_y": grid_y,
        "cubic": cubic,
        "nearest": nearest,
        "selected": selected,
    }


def interpolation_audit(
    source_data: dict[str, np.ndarray],
    threshold: int,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    plot_grid = interpolate_source(source_data, threshold, 300)
    cubic = plot_grid["cubic"]
    smooth = gaussian_filter(cubic, sigma=6)
    nearest_smooth = gaussian_filter(plot_grid["nearest"], sigma=6)

    skeleton_grid = interpolate_source(source_data, threshold, 250)
    skeleton_smooth = gaussian_filter(skeleton_grid["cubic"], sigma=4)
    sx = sobel(skeleton_smooth, axis=0)
    sy = sobel(skeleton_smooth, axis=1)
    gradient_magnitude = np.hypot(sx, sy)
    mask = gradient_magnitude > np.percentile(gradient_magnitude, 85)
    eroded = binary_erosion(mask, iterations=2)

    flow_grid = interpolate_source(source_data, threshold, 150)
    flow = flow_grid["cubic"]
    spacing = float(flow_grid["axis"][1] - flow_grid["axis"][0])
    grad_x = -(
        np.roll(flow, -1, axis=0) - np.roll(flow, 1, axis=0)
    ) / (2 * spacing)
    grad_y = -(
        np.roll(flow, -1, axis=1) - np.roll(flow, 1, axis=1)
    ) / (2 * spacing)
    divergence = (
        np.roll(grad_x, -1, axis=0) - np.roll(grad_x, 1, axis=0)
    ) / (2 * spacing) + (
        np.roll(grad_y, -1, axis=1) - np.roll(grad_y, 1, axis=1)
    ) / (2 * spacing)
    curl = (
        np.roll(grad_y, -1, axis=0) - np.roll(grad_y, 1, axis=0)
    ) / (2 * spacing) - (
        np.roll(grad_x, -1, axis=1) - np.roll(grad_x, 1, axis=1)
    ) / (2 * spacing)

    boundary_values = np.concatenate(
        (
            divergence[:2, :].ravel(),
            divergence[-2:, :].ravel(),
            divergence[:, :2].ravel(),
            divergence[:, -2:].ravel(),
        )
    )
    interior_values = divergence[2:-2, 2:-2]

    audit = {
        "grid_size": 300,
        "finite_fraction": rounded(float(np.isfinite(cubic).mean())),
        "cubic_transpose": transpose_record(cubic),
        "smooth_sigma_6_transpose": transpose_record(smooth),
        "nearest_transpose": transpose_record(plot_grid["nearest"]),
        "nearest_smooth_transpose": transpose_record(nearest_smooth),
        "skeleton": {
            "operation": (
                "top-15-percent Sobel magnitude followed by two binary "
                "erosions; not a mathematical skeletonization"
            ),
            "mask_fraction": rounded(float(mask.mean())),
            "eroded_fraction": rounded(float(eroded.mean())),
            "mask_transpose_disagreements": int(
                np.count_nonzero(mask != mask.T)
            ),
            "eroded_transpose_disagreements": int(
                np.count_nonzero(eroded != eroded.T)
            ),
        },
        "flow": {
            "gradient_is_derived_from_interpolant": True,
            "curl_max_abs": rounded(float(np.max(np.abs(curl)))),
            "curl_rms": rounded(float(np.sqrt(np.mean(curl**2)))),
            "divergence_max_abs": rounded(
                float(np.max(np.abs(divergence)))
            ),
            "divergence_rms": rounded(
                float(np.sqrt(np.mean(divergence**2)))
            ),
            "periodic_roll_boundary": True,
            "boundary_rms": rounded(
                float(np.sqrt(np.mean(boundary_values**2)))
            ),
            "interior_rms": rounded(
                float(np.sqrt(np.mean(interior_values**2)))
            ),
        },
    }
    arrays = {
        "axis": plot_grid["axis"],
        "cubic": cubic,
        "smooth": smooth,
        "transpose_residual": cubic - cubic.T,
        "gradient_magnitude": gradient_magnitude,
        "eroded": eroded,
    }
    return audit, arrays


def arithmetic_audit(
    all_primes: np.ndarray,
    sampled_primes: np.ndarray,
    source_data: dict[str, np.ndarray],
    threshold: int,
) -> dict[str, Any]:
    defect = source_data["defect"]
    root = source_data["root"]
    selected = defect <= threshold
    residue_counts = Counter(int(value) for value in (defect % 8).ravel())
    selected_residue_counts = Counter(
        int(value) for value in (defect[selected] % 8).ravel()
    )
    low_counts = Counter(
        int(value) for value in defect[defect <= 50].ravel()
    )

    p3_index = int(np.flatnonzero(sampled_primes == 3)[0])
    p3_nine = int(
        np.count_nonzero(defect[p3_index, :] == 9)
        + np.count_nonzero(defect[:, p3_index] == 9)
        - int(defect[p3_index, p3_index] == 9)
    )

    full_odd_primes = all_primes[all_primes >= 3]
    pell_hits: list[dict[str, int]] = []
    for prime in full_odd_primes:
        pell_root = math.isqrt(2 * int(prime) * int(prime))
        pell_defect = 2 * int(prime) * int(prime) - pell_root**2
        if pell_defect == 1:
            pell_hits.append({"p": int(prime), "r": pell_root})

    sampled_pell_hits = [
        hit for hit in pell_hits if hit["p"] in set(sampled_primes.tolist())
    ]
    return {
        "sampled_prime_count": int(sampled_primes.size),
        "sampled_prime_min": int(sampled_primes.min()),
        "sampled_prime_max": int(sampled_primes.max()),
        "pair_count": int(defect.size),
        "selected_pair_count": int(np.count_nonzero(selected)),
        "selected_fraction": rounded(float(selected.mean())),
        "defect_minimum": int(defect.min()),
        "defect_maximum": int(defect.max()),
        "pythagorean_equality_count": int(np.count_nonzero(defect == 0)),
        "root_is_required_prime": False,
        "exact_formula": (
            "E(p,q)=p^2+q^2-floor(sqrt(p^2+q^2))^2"
        ),
        "radial_factorization": (
            "E=(rho-floor(rho))*(rho+floor(rho)), "
            "rho=sqrt(p^2+q^2)"
        ),
        "threshold_shell_width": (
            "E<=T iff 0<=rho-floor(rho)<="
            "sqrt(floor(rho)^2+T)-floor(rho)"
        ),
        "transpose_max_abs": int(np.max(np.abs(defect - defect.T))),
        "allowed_defect_residues_mod_8": sorted(residue_counts),
        "residue_counts_mod_8": dict(sorted(residue_counts.items())),
        "selected_residue_counts_mod_8": dict(
            sorted(selected_residue_counts.items())
        ),
        "low_defect_counts_through_50": dict(sorted(low_counts.items())),
        "p_equals_3_axis_defect_9_count": p3_nine,
        "p_equals_3_axis_share_of_defect_le_10": rounded(
            p3_nine / max(1, int(np.count_nonzero(defect <= 10)))
        ),
        "diagonal_pell_equation": (
            "E(p,p)=1 iff r^2-2p^2=-1"
        ),
        "pell_prime_hits_below_2000": pell_hits,
        "sampled_pell_hits": sampled_pell_hits,
        "source_subsampling_misses_pell_hits": not sampled_pell_hits,
        "modular_reason_no_zero": (
            "odd prime squares are 1 mod 4, so p^2+q^2 is 2 mod 4 "
            "and cannot equal a square"
        ),
    }


def randomized_null_audit(
    all_primes: np.ndarray,
    sampled_primes: np.ndarray,
    source_data: dict[str, np.ndarray],
) -> dict[str, Any]:
    thresholds = np.array([10, 25, 50, 100, 200, 300, 500])

    def threshold_fractions(values: np.ndarray) -> np.ndarray:
        defect = defect_data(values)["defect"]
        return np.array(
            [np.mean(defect <= threshold) for threshold in thresholds]
        )

    observed = threshold_fractions(sampled_primes)
    odd_control = (
        np.rint(
            np.linspace(
                int(sampled_primes.min()),
                int(sampled_primes.max()),
                sampled_primes.size,
            )
        ).astype(np.int64)
        | 1
    )
    odd_data = defect_data(odd_control)

    prime_pool = all_primes[all_primes >= 3]
    prime_blocks = [
        prime_pool[index : index + 7]
        for index in range(0, prime_pool.size, 7)
    ]
    if len(prime_blocks) != sampled_primes.size:
        raise RuntimeError("prime-block null no longer matches source sample")

    rng = np.random.default_rng(20_260_728)
    repetitions = 2_000
    simulations = np.empty((repetitions, thresholds.size))
    for row in range(repetitions):
        sample = np.array(
            [rng.choice(block) for block in prime_blocks],
            dtype=np.int64,
        )
        simulations[row] = threshold_fractions(sample)

    means = simulations.mean(axis=0)
    standard_deviations = simulations.std(axis=0, ddof=1)
    z_scores = np.divide(
        observed - means,
        standard_deviations,
        out=np.zeros_like(observed),
        where=standard_deviations > 0,
    )
    upper_quantiles = np.mean(simulations <= observed, axis=0)

    return {
        "count_matched_odd_control": {
            "values": odd_control.tolist(),
            "selected_fraction_at_500": rounded(
                float(np.mean(odd_data["defect"] <= 500))
            ),
            "transpose_max_abs": int(
                np.max(
                    np.abs(
                        odd_data["defect"] - odd_data["defect"].T
                    )
                )
            ),
            "allowed_residues_mod_8": sorted(
                int(value)
                for value in np.unique(odd_data["defect"] % 8)
            ),
        },
        "prime_block_null": {
            "description": (
                "choose one prime uniformly from each consecutive block "
                "of seven primes; preserves primality, count, and rank scale"
            ),
            "seed": 20_260_728,
            "repetitions": repetitions,
            "thresholds": thresholds.tolist(),
            "observed": [rounded(float(value)) for value in observed],
            "mean": [rounded(float(value)) for value in means],
            "standard_deviation": [
                rounded(float(value)) for value in standard_deviations
            ],
            "z_score": [rounded(float(value)) for value in z_scores],
            "empirical_lower_quantile": [
                rounded(float(value)) for value in upper_quantiles
            ],
            "multiple_threshold_guard": (
                "descriptive selected-sample diagnostics only; thresholds "
                "were inspected jointly and are not theorem evidence"
            ),
        },
    }


def spectral_audit(
    source_data: dict[str, np.ndarray],
    odd_values: np.ndarray,
    threshold: int,
) -> dict[str, Any]:
    defect = source_data["defect"].astype(float)
    adjacency = (defect <= threshold).astype(float)
    closeness = np.maximum(0.0, threshold + 1.0 - defect)
    exponential = np.exp(-defect / 100.0)
    dimension = defect.shape[0]
    center = np.eye(dimension) - np.ones((dimension, dimension)) / dimension

    laplacian = np.diag(closeness.sum(axis=1)) - closeness
    odd_defect = defect_data(odd_values)["defect"].astype(float)
    odd_closeness = np.maximum(0.0, threshold + 1.0 - odd_defect)

    return {
        "raw_symmetric_kernels": {
            "defect": eigen_record(defect),
            "threshold_adjacency": eigen_record(adjacency),
            "closeness": eigen_record(closeness),
            "exp_minus_defect_over_100": eigen_record(exponential),
        },
        "centered_on_constant_orthogonal_subspace": {
            "defect": eigen_record(center @ defect @ center),
            "closeness": eigen_record(center @ closeness @ center),
        },
        "count_matched_odd_closeness": eigen_record(odd_closeness),
        "graph_laplacian_of_closeness": eigen_record(laplacian),
        "interpretation": (
            "real symmetry gives a self-adjoint matrix but not a sign; "
            "the raw and centered kernels are indefinite. The graph "
            "Laplacian is positive semidefinite by construction, which "
            "does not create an identity with Xi."
        ),
    }


def zero_overlay_audit() -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    precision_records: dict[str, Any] = {}
    scaled_residuals: dict[int, np.ndarray] = {}
    real_roots: np.ndarray | None = None

    for digits in (30, 40, 50):
        mp.mp.dps = digits
        roots = [
            mp.findroot(
                lambda height: mp.zeta(mp.mpf("0.5") + 1j * height),
                guess,
            )
            for guess in ZERO_GUESSES
        ]
        real = np.array([float(mp.re(root)) for root in roots])
        residual = np.array([float(mp.im(root)) for root in roots])
        span = float(np.ptp(residual))
        scaled = (residual - residual.min()) / span
        scaled_residuals[digits] = scaled
        if digits == 40:
            real_roots = real
        precision_records[str(digits)] = {
            "imaginary_residual_span": span,
            "maximum_abs_imaginary_residual": float(
                np.max(np.abs(residual))
            ),
            "scaled_residual_coordinates": [
                rounded(float(value), 12) for value in scaled
            ],
        }

    if real_roots is None:
        raise RuntimeError("40-digit zero roots were not computed")
    correct_scaled = (real_roots - real_roots.min()) / np.ptp(real_roots)

    correlations_to_correct = {
        str(digits): rounded(
            float(np.corrcoef(values, correct_scaled)[0, 1])
        )
        for digits, values in scaled_residuals.items()
    }
    precision_instability = {
        "30_vs_40_max_abs": rounded(
            float(np.max(np.abs(scaled_residuals[30] - scaled_residuals[40])))
        ),
        "40_vs_50_max_abs": rounded(
            float(np.max(np.abs(scaled_residuals[40] - scaled_residuals[50])))
        ),
        "30_vs_50_max_abs": rounded(
            float(np.max(np.abs(scaled_residuals[30] - scaled_residuals[50])))
        ),
    }
    plot_span = 1999 - 3
    residual_span_40 = precision_records["40"]["imaginary_residual_span"]
    audit = {
        "solver_variable": (
            "height t in zeta(1/2+i*t)=0; the desired height is Re(t)"
        ),
        "source_extracts": "Im(t), a numerical root residual",
        "source_assigns": "rz_x=scaled and rz_y=scaled",
        "diagonal_is_imposed": True,
        "correct_root_min": rounded(float(real_roots.min())),
        "correct_root_max": rounded(float(real_roots.max())),
        "correct_root_span": rounded(float(np.ptp(real_roots))),
        "correct_scaled_coordinates": [
            rounded(float(value), 12) for value in correct_scaled
        ],
        "precision_records": precision_records,
        "correlation_of_residual_scaling_with_correct_heights": (
            correlations_to_correct
        ),
        "precision_instability": precision_instability,
        "40_digit_minmax_magnification": float(
            plot_span / residual_span_40
        ),
        "verdict": (
            "the plotted cyan coordinates encode precision-dependent "
            "solver residuals on an imposed diagonal, not zero heights "
            "or an observed relation to the curvature field"
        ),
    }
    arrays = {
        "correct_scaled": correct_scaled,
        "residual_scaled_30": scaled_residuals[30],
        "residual_scaled_40": scaled_residuals[40],
        "residual_scaled_50": scaled_residuals[50],
    }
    return audit, arrays


def old_operator_audit() -> dict[str, str]:
    return {
        "laplacian_identity": (
            "For a holomorphic function f(s), "
            "(partial_sigma^2+partial_t^2)f=0 everywhere."
        ),
        "legacy_internal_conflict": (
            "The legacy note uses analyticity to state Delta xi=0 at a "
            "critical-line zero, then claims an off-line zero would have "
            "Delta xi!=0. Analyticity makes Delta xi=0 off the line too."
        ),
        "semigroup_consequence": (
            "The proposed two-dimensional Laplacian flow acts trivially "
            "on holomorphic xi wherever it is defined; it cannot attract "
            "off-line zeros to the critical line."
        ),
        "multiplier_guard": (
            "Multiplying xi(1/2+it) by a nonzero positive function of t "
            "preserves its zeros tautologically and is not spectral "
            "localization evidence."
        ),
    }


def candidate_bridge() -> dict[str, str]:
    return {
        "prime_heat_sum": (
            "A(y)=sum_(n>=1) Lambda(n) exp(-pi*n^2*y), y>0"
        ),
        "one_variable_mellin": (
            "integral_0^infinity A(y)y^(s-1)dy="
            "Gamma(s)pi^(-s)[-zeta'(2s)/zeta(2s)], Re(s)>1/2"
        ),
        "radial_pair_heat": (
            "A(y)^2=sum_(m,n>=1)Lambda(m)Lambda(n)"
            "exp(-pi*y*(m^2+n^2))"
        ),
        "radial_pair_mellin": (
            "integral_0^infinity A(y)^2 y^(s-1)dy="
            "Gamma(s)pi^(-s)sum_(m,n)Lambda(m)Lambda(n)"
            "(m^2+n^2)^(-s), initially Re(s)>1"
        ),
        "pi_provenance": (
            "pi is introduced only by the standard Gaussian Fourier/theta "
            "normalization exp(-pi*n^2*y). Replacing pi by a>0 changes "
            "the Mellin factor to a^(-s); no circle or fitted pi is read "
            "from the legacy images."
        ),
        "surviving_question": (
            "Can a completed, regularized prime-weighted radial kernel be "
            "related exactly to the Xi contact scalar or to a known Weil "
            "positivity form, with every signed prime-power correction "
            "retained?"
        ),
        "nonpromotion": (
            "The prime-weighted radial series is a legitimate zeta-aware "
            "replacement for the picture, but its symmetry alone proves "
            "neither reflection positivity nor RH."
        ),
    }


def current_rh_map() -> dict[str, str]:
    return {
        "already_realized_geometry": (
            "The current RH corpus already has an exact positive "
            "semidefinite sine-feature kernel P_(K,R), a positive "
            "edge-feature Gram kernel G_K, and a curvature-weighted energy."
        ),
        "exact_wall": (
            "After the geometric diagonal closes, the remaining term is "
            "the signed Mobius-labelled off-diagonal edge correlation. "
            "The symmetrized threshold operator is still indefinite, with "
            "an exact 2x2 determinant -1/4."
        ),
        "contact_route": (
            "The active Xi route analogously needs a contact-conditioned "
            "signed lower bound for mathcal_C_N before absolute values, "
            "plus a one-turn phase budget."
        ),
        "lesson_from_legacy": (
            "Both the old pictures and the present obstruction separate "
            "self-adjoint symmetry from coercive sign. The wall is not a "
            "lack of symmetry; it is the missing arithmetic control of "
            "the negative or signed sector."
        ),
        "admissible_use": (
            "Use the old field as a toy falsification model for proposed "
            "symmetry-to-energy conversions. Promote it only if an exact "
            "identity maps the actual Xi/Mobius coefficients into a "
            "positive form with a controlled remainder."
        ),
    }


def build_rows() -> list[AuditRow]:
    return [
        AuditRow(
            "lpc_01_provenance",
            "source_audit",
            "available_exact",
            "The supplied images have a complete local legacy source.",
            "The source is external and read-only; no MTS file is edited.",
        ),
        AuditRow(
            "lpc_02_definition",
            "exact_identity",
            "available_exact",
            "The source field is E=p^2+q^2-floor(sqrt(p^2+q^2))^2.",
            "The third coordinate is an integer floor, not a prime.",
        ),
        AuditRow(
            "lpc_03_radial_shell",
            "exact_identity",
            "available_exact",
            "E is a radial sawtooth and E<=500 selects thin circular shells.",
            "This geometry exists for arbitrary inputs, not only primes.",
        ),
        AuditRow(
            "lpc_04_transpose",
            "exact_symmetry",
            "available_exact",
            "E(p,q)=E(q,p) and the censored point cloud is transpose-symmetric.",
            "The symmetry follows algebraically from p^2+q^2.",
        ),
        AuditRow(
            "lpc_05_modular_residue",
            "exact_arithmetic",
            "available_exact",
            "Odd-prime defects lie only in residue classes 1, 2, and 6 mod 8.",
            "This does not connect those defects to zeta zeros.",
        ),
        AuditRow(
            "lpc_06_no_exact_triangle",
            "exact_arithmetic",
            "available_exact",
            "No sampled odd-prime pair has E=0 by the mod-4 obstruction.",
            "The plot studies near-squares, not prime Pythagorean triples.",
        ),
        AuditRow(
            "lpc_07_boundary_axis",
            "artifact_guard",
            "guard_active",
            "Most E<=10 source hits are the deterministic p=3/q=3 axes.",
            "Dense lower-left structure cannot be read as a bulk law.",
        ),
        AuditRow(
            "lpc_08_censor_interpolate",
            "artifact_guard",
            "guard_active",
            "The plot discards E>500 and cubic-interpolates only retained points.",
            "The displayed surface is not the full defect function.",
        ),
        AuditRow(
            "lpc_09_interpolation_symmetry",
            "finite_validation",
            "validated_finite",
            "Cubic triangulation measurably breaks exact transpose symmetry.",
            "Nearest interpolation preserves it; smoothing only reduces the break.",
        ),
        AuditRow(
            "lpc_10_skeleton",
            "artifact_guard",
            "guard_active",
            "The PGF backbone is a percentile mask plus binary erosion.",
            "No topological skeleton or persistence invariant was computed.",
        ),
        AuditRow(
            "lpc_11_flow",
            "exact_numerical_identity",
            "available_exact",
            "The reported flow is the negative gradient of the interpolant.",
            "Its curl vanishes by commuting finite differences; roll boundaries are periodic.",
        ),
        AuditRow(
            "lpc_12_zero_component",
            "bug_certificate",
            "available_exact",
            "The source plots Im(t), not the zero height Re(t).",
            "Im(t) is only a precision-dependent root residual.",
        ),
        AuditRow(
            "lpc_13_zero_diagonal",
            "artifact_guard",
            "guard_active",
            "The source explicitly sets both zero-overlay axes to one scaled array.",
            "The diagonal alignment is imposed.",
        ),
        AuditRow(
            "lpc_14_null_control",
            "finite_validation",
            "validated_finite",
            "Odd and prime-block controls retain shell and transpose structure.",
            "Selected finite threshold deviations are descriptive and multiple-tested.",
        ),
        AuditRow(
            "lpc_15_spectral_sign",
            "finite_obstruction",
            "validated_finite",
            "Natural symmetric defect kernels have mixed positive and negative spectra.",
            "Self-adjointness alone supplies no coercive energy.",
        ),
        AuditRow(
            "lpc_16_laplacian",
            "nonpromotion_guard",
            "guard_active",
            "A graph Laplacian made from nonnegative closeness weights is PSD.",
            "That positivity is constructed and has no proved Xi identity.",
        ),
        AuditRow(
            "lpc_17_pell_survivor",
            "exact_arithmetic",
            "available_exact",
            "The diagonal E=1 locus is the negative Pell equation r^2-2p^2=-1.",
            "This is genuine quadratic-form arithmetic but the source sample misses its hits.",
        ),
        AuditRow(
            "lpc_18_prime_heat_bridge",
            "hypothesis_generator",
            "open",
            "A von-Mangoldt radial heat sum gives an exact Mellin bridge to zeta'/zeta.",
            "No completed reflection-positive Xi/contact identity is proved.",
        ),
        AuditRow(
            "lpc_19_current_wall",
            "route_decision",
            "available_exact",
            "The current RH wall is signed arithmetic coercivity, not missing symmetry.",
            "Keep the legacy idea as a falsification side model; do not replace the active route.",
        ),
    ]


def render_figure(
    path: Path,
    source_data: dict[str, np.ndarray],
    interpolation_arrays: dict[str, np.ndarray],
    odd_values: np.ndarray,
    spectral: dict[str, Any],
    zero_arrays: dict[str, np.ndarray],
) -> None:
    odd_data = defect_data(odd_values)
    odd_grid = interpolate_source(odd_data, 500, 300)
    odd_smooth = gaussian_filter(odd_grid["cubic"], sigma=6)
    extent = [
        float(interpolation_arrays["axis"].min()),
        float(interpolation_arrays["axis"].max()),
        float(interpolation_arrays["axis"].min()),
        float(interpolation_arrays["axis"].max()),
    ]

    figure, axes = plt.subplots(2, 3, figsize=(15, 9))
    image = axes[0, 0].imshow(
        interpolation_arrays["smooth"].T,
        origin="lower",
        extent=extent,
        cmap="inferno",
        aspect="equal",
    )
    axes[0, 0].set_title("Legacy pipeline: prime sample")
    axes[0, 0].set_xlabel("p")
    axes[0, 0].set_ylabel("q")
    figure.colorbar(image, ax=axes[0, 0], fraction=0.046)

    residual = interpolation_arrays["transpose_residual"]
    limit = float(np.max(np.abs(residual)))
    image = axes[0, 1].imshow(
        residual.T,
        origin="lower",
        extent=extent,
        cmap="coolwarm",
        vmin=-limit,
        vmax=limit,
        aspect="equal",
    )
    axes[0, 1].set_title("Cubic interpolation: E-E^T")
    axes[0, 1].set_xlabel("p")
    axes[0, 1].set_ylabel("q")
    figure.colorbar(image, ax=axes[0, 1], fraction=0.046)

    image = axes[0, 2].imshow(
        odd_smooth.T,
        origin="lower",
        extent=extent,
        cmap="inferno",
        aspect="equal",
    )
    axes[0, 2].set_title("Count-matched odd control")
    axes[0, 2].set_xlabel("odd input")
    axes[0, 2].set_ylabel("odd input")
    figure.colorbar(image, ax=axes[0, 2], fraction=0.046)

    image = axes[1, 0].imshow(
        source_data["defect"].T,
        origin="lower",
        cmap="viridis",
        aspect="equal",
    )
    axes[1, 0].set_title("Exact 44 x 44 defect matrix")
    axes[1, 0].set_xlabel("sample index")
    axes[1, 0].set_ylabel("sample index")
    figure.colorbar(image, ax=axes[1, 0], fraction=0.046)

    colors = {
        "defect": "#d1495b",
        "threshold_adjacency": "#00798c",
        "closeness": "#edae49",
        "exp_minus_defect_over_100": "#30638e",
    }
    for name, record in spectral["raw_symmetric_kernels"].items():
        axes[1, 1].plot(
            record["values"],
            label=name.replace("_", " "),
            color=colors[name],
            linewidth=1.5,
        )
    axes[1, 1].axhline(0, color="black", linewidth=0.8)
    axes[1, 1].set_title("Symmetric kernel eigenvalues")
    axes[1, 1].set_xlabel("ordered eigenvalue")
    axes[1, 1].set_ylabel("value")
    axes[1, 1].set_yscale("symlog", linthresh=1)
    axes[1, 1].legend(fontsize=8)

    indices = np.arange(1, len(ZERO_GUESSES) + 1)
    axes[1, 2].plot(
        indices,
        zero_arrays["correct_scaled"],
        "o-",
        label="correct Re(t)",
        color="#00798c",
        markersize=4,
    )
    axes[1, 2].plot(
        indices,
        zero_arrays["residual_scaled_40"],
        "x--",
        label="source Im(t), 40 dps",
        color="#d1495b",
        markersize=5,
    )
    axes[1, 2].set_title("Zero coordinate extraction")
    axes[1, 2].set_xlabel("zero index")
    axes[1, 2].set_ylabel("min-max scaled coordinate")
    axes[1, 2].set_ylim(-0.05, 1.05)
    axes[1, 2].legend(fontsize=8)

    figure.suptitle(
        "Legacy Prime-Curvature Symmetry Controls",
        fontsize=16,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=160)
    plt.close(figure)


def make_note(artifact: dict[str, Any]) -> str:
    arithmetic = artifact["arithmetic"]
    interpolation = artifact["interpolation"]
    zero = artifact["zero_overlay"]
    spectral = artifact["spectral"]
    null = artifact["null_controls"]["prime_block_null"]
    defect_spectrum = spectral["raw_symmetric_kernels"]["defect"]
    closeness_spectrum = spectral["raw_symmetric_kernels"]["closeness"]
    laplacian_spectrum = spectral["graph_laplacian_of_closeness"]
    attachment = artifact["attachments"]

    threshold_rows = "\n".join(
        "| {threshold} | {observed:.6f} | {mean:.6f} | {z:.3f} |".format(
            threshold=threshold,
            observed=observed,
            mean=mean,
            z=z_score,
        )
        for threshold, observed, mean, z_score in zip(
            null["thresholds"],
            null["observed"],
            null["mean"],
            null["z_score"],
        )
    )

    spectrum_rows = "\n".join(
        "| {name} | {negative} | {zero} | {positive} | {minimum:.6g} | "
        "{maximum:.6g} |".format(
            name=name.replace("_", " "),
            negative=record["negative_count"],
            zero=record["zero_count"],
            positive=record["positive_count"],
            minimum=record["minimum"],
            maximum=record["maximum"],
        )
        for name, record in spectral["raw_symmetric_kernels"].items()
    )

    return f"""# Legacy Prime-Curvature Symmetry Audit

Date: 2026-07-28

Status: exploratory source reconstruction, artifact control, and
finite spectral audit. This is not a proof of an Xi contact gap,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
outputs/{STEM}_controls.png
```

## Scope And Provenance

The complete plot source was found in the local read-only legacy
Motion-TimeSpace archive. GitHub was not needed, and no external project
file was edited. The eight supplied JPEGs contain {attachment["unique_hashes"]}
unique files; Photos 1 and 2 are byte-identical.

The source uses `N=2000`, every seventh prime starting at 3, and
`threshold=500`. That gives {arithmetic["sampled_prime_count"]} sampled
primes, {arithmetic["pair_count"]} ordered pairs, and
{arithmetic["selected_pair_count"]} retained pairs.

## Exact Geometry

The plotted scalar is

```text
E(p,q)=p^2+q^2-floor(sqrt(p^2+q^2))^2.
```

Writing `rho=sqrt(p^2+q^2)` and `m=floor(rho)` gives

```text
E=(rho-m)(rho+m).
```

Thus the raw continuous formula is radial, and `E<=T` selects thin
annular pieces next to the circles `p^2+q^2=m^2`. The visible arcs are
therefore expected before any prime theorem. Swapping `p` and `q`
leaves the formula unchanged, so the exact 44 by 44 defect matrix and
the threshold mask are transpose-symmetric by construction.

The integer `r=m` is not required to be prime. Exact equality is
impossible for odd `p,q`: `p^2+q^2=2 mod 4`, whereas a square is 0 or
1 mod 4. The observed defects occupy only residue classes
{arithmetic["allowed_defect_residues_mod_8"]} modulo 8.

There is genuine arithmetic beneath the picture. On the diagonal,
`E(p,p)=1` is exactly the negative Pell equation

```text
r^2-2p^2=-1.
```

The prime hits below 2000 are `{arithmetic["pell_prime_hits_below_2000"]}`,
but the every-seventh-prime source sample misses both. This is a real
quadratic-form thread, not evidence about zeta zeros.

## Rendering Controls

Only points with `E<=500` are passed to cubic `griddata`; values above
the threshold are omitted rather than plotted. Cubic triangulation
breaks the exact transpose symmetry:

```text
raw defect max |E-E^T| = {arithmetic["transpose_max_abs"]}
cubic RMS |Z-Z^T|      = {interpolation["cubic_transpose"]["rms"]:.9g}
cubic max |Z-Z^T|      = {interpolation["cubic_transpose"]["max_abs"]:.9g}
smoothed relative RMS  = {interpolation["smooth_sigma_6_transpose"]["relative_rms"]:.9g}
nearest max |Z-Z^T|    = {interpolation["nearest_transpose"]["max_abs"]:.9g}
```

The lower-left density has a concrete boundary cause. Of the
`E<=10` ordered hits, {arithmetic["p_equals_3_axis_defect_9_count"]}
are the `p=3` or `q=3` axes with the deterministic value `E=9`,
a share of {arithmetic["p_equals_3_axis_share_of_defect_le_10"]:.3%}.

The named "topological skeleton" is the upper 15 percent of a Sobel
gradient magnitude followed by two binary erosions. It is not
skeletonization and computes no topological invariant. The divergence
is the negative discrete Laplacian of the interpolated scalar. The
curl is {interpolation["flow"]["curl_max_abs"]:.3g} at worst, as expected
when commuting the same roll-based finite differences; `np.roll`
also imposes an unphysical periodic boundary.

## Arithmetic Null

The exact radial and transpose symmetries survive count-matched odd
inputs. A more conservative finite null chooses one prime from each
successive block of seven primes:

| defect threshold | source fraction | null mean | z score |
|---:|---:|---:|---:|
{threshold_rows}

These are selected-sample, jointly inspected finite diagnostics. They
do not establish a prime anomaly. In particular, the source pictures
use threshold 500, where the shell geometry dominates.

## Spectral Sign

Every tested raw kernel is real symmetric, but all are indefinite:

| kernel | negative | zero | positive | min eigenvalue | max eigenvalue |
|---|---:|---:|---:|---:|---:|
{spectrum_rows}

For example, the defect matrix has
{defect_spectrum["negative_count"]} negative and
{defect_spectrum["positive_count"]} positive eigenvalues; the
nonnegative closeness matrix still has
{closeness_spectrum["negative_count"]} negative eigenvalues.
Centering away the constant vector does not fix the sign.

The graph Laplacian `diag(W 1)-W` is positive semidefinite, with
{laplacian_spectrum["negative_count"]} negative eigenvalues and
{laplacian_spectrum["zero_count"]} zero mode. That positivity is true
for every nonnegative symmetric graph and is created by the
transformation. It becomes relevant to RH only if the actual Xi
contact scalar is proved equal to, or coercively bounded by, that
energy with every signed remainder retained.

## Zero Overlay

The source solves

```text
zeta(1/2+i*t)=0.
```

Here the zero height is `Re(t)`, ranging from
{zero["correct_root_min"]:.12g} to {zero["correct_root_max"]:.12g}.
The code instead extracts `Im(t)`. At 40 decimal digits those values
span only {zero["precision_records"]["40"]["imaginary_residual_span"]:.3g};
they are numerical root residuals. Min-max scaling magnifies them by
about {zero["40_digit_minmax_magnification"]:.3g}, and changing the
working precision moves a scaled coordinate by as much as
{zero["precision_instability"]["40_vs_50_max_abs"]:.3f}.

Finally, the source sets `rz_x=scaled` and `rz_y=scaled`. The cyan
diagonal is therefore imposed twice: by assigning the same coordinate
to both axes, and by stretching numerical residuals rather than zero
heights. The overlay contains no measured zero/curvature alignment.

## Old Invariance Guard

For any holomorphic `f(s)`,

```text
(partial_sigma^2+partial_t^2)f=0
```

throughout its analytic domain. The legacy note correctly invokes this
at one point, but later asserts that an off-line zero would have a
nonzero Laplacian. That step is incompatible with analyticity. The
proposed two-dimensional Laplacian flow is trivial on holomorphic xi;
similarly, multiplying critical-line xi by a nonzero scalar function
preserves its zeros tautologically.

This is the same logical warning exposed by the pictures: an invariant
can be exact because it was built into the map, without localizing the
unknown zeros.

## Zeta-Aware Replacement

The defensible part of the old idea is a radial prime heat correlation.
With the von Mangoldt weight, define

```text
A(y)=sum Lambda(n) exp(-pi*n^2*y).
```

For `Re(s)>1/2`, termwise Mellin transformation gives

```text
integral A(y)y^(s-1)dy
 =Gamma(s)pi^(-s)[-zeta'(2s)/zeta(2s)].
```

Also

```text
A(y)^2
 =sum Lambda(m)Lambda(n)exp(-pi*y*(m^2+n^2)).
```

This is an exact zeta-aware version of radial prime-pair geometry. The
`pi` is from the standard Gaussian Fourier/theta convention; replacing
it by any `a>0` simply gives `a^(-s)`. It is not inferred from a circle
in the images.

The open thought experiment is precise: construct a completed,
regularized version of this radial kernel and ask whether its
reflection form is exactly a known Weil positivity form or exactly the
current Xi contact scalar. Prime restriction destroys the ordinary
lattice theta functional equation, so positivity cannot be assumed.

## What It Says About The Wall

The current corpus has already converted the useful geometric
intuition into exact mathematics:

- an explicit positive sine-feature kernel `P_(K,R)`;
- a positive edge-feature Gram kernel `G_K`;
- a curvature-weighted energy;
- exact self-adjoint symmetrization.

The diagonal is controlled. The obstruction that remains is the signed
Mobius-labelled off-diagonal correlation; the symmetrized threshold
operator itself has an exact indefinite 2 by 2 minor. On the active Xi
route, the analogous missing result is the contact-conditioned signed
lower bound for `mathcal_C_N`, followed by a one-turn phase budget.

So the old symmetry does explain the repeated wall, but not by giving
the answer: we already have symmetry and even a Gram representation.
What we do not have is a theorem forcing the actual arithmetic
coefficient vector to avoid the negative sector.

## Route Decision

Keep the active endpoint-complete Abel/contact route primary. Retain
the legacy prime-curvature field as a small falsification laboratory
for candidate symmetry-to-energy transformations. Pursue the
von-Mangoldt radial heat bridge only if it yields an exact completed
identity with the actual Xi/contact quantity; otherwise it remains a
separate quadratic-form investigation.

## Boundary

The source reconstruction, radial factorization, modular residues,
transpose identity, interpolation controls, root-coordinate bug,
finite spectra, Pell reduction, and Mellin bridge are exact or
reproducibly finite as stated. No zero/curvature correlation,
reflection-positive prime kernel, Xi contact gap, signed Mobius gain,
one-turn phase theorem, `Lambda<=0`, PF-infinity, RH proof, or
Clay-prize conclusion is obtained.
"""


def build_artifact(figure_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    sources = source_audit()
    attachments = attachment_audit()
    all_primes = primes_below(2_000)
    sampled_primes = all_primes[all_primes >= 3][::7]
    source_data = defect_data(sampled_primes)

    arithmetic = arithmetic_audit(
        all_primes, sampled_primes, source_data, threshold=500
    )
    interpolation, interpolation_arrays = interpolation_audit(
        source_data, threshold=500
    )
    null_controls = randomized_null_audit(
        all_primes, sampled_primes, source_data
    )
    odd_values = np.array(
        null_controls["count_matched_odd_control"]["values"],
        dtype=np.int64,
    )
    spectral = spectral_audit(source_data, odd_values, threshold=500)
    zero_overlay, zero_arrays = zero_overlay_audit()

    render_figure(
        figure_path,
        source_data,
        interpolation_arrays,
        odd_values,
        spectral,
        zero_arrays,
    )

    rows = build_rows()
    artifact = {
        "kind": STEM,
        "date": "2026-07-28",
        "status": "exploratory_symmetry_and_artifact_audit_complete",
        "resource_policy": (
            "serial one-process audit; BLAS/OpenMP thread caps set to one; "
            "no background worker or external-project edit"
        ),
        "sources": sources,
        "attachments": attachments,
        "arithmetic": arithmetic,
        "interpolation": interpolation,
        "null_controls": null_controls,
        "spectral": spectral,
        "zero_overlay": zero_overlay,
        "old_operator": old_operator_audit(),
        "candidate_bridge": candidate_bridge(),
        "current_rh_map": current_rh_map(),
        "control_figure": {
            "path": str(figure_path.relative_to(REPO_ROOT)),
            "sha256": file_hash(figure_path),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "row_count": len(rows),
            "exact_or_finite_rows": len(
                [
                    row
                    for row in rows
                    if row.readiness
                    in {"available_exact", "validated_finite"}
                ]
            ),
            "guard_rows": len(
                [row for row in rows if row.readiness == "guard_active"]
            ),
            "open_hypothesis_rows": len(
                [row for row in rows if row.readiness == "open"]
            ),
        },
        "route_decision": (
            "Do not use the legacy pictures or zero overlay as RH evidence. "
            "Keep the current endpoint-complete Abel/contact route primary. "
            "Retain only the exact quadratic-form/Pell residue and the "
            "von-Mangoldt radial heat/Mellin bridge as explicitly separate "
            "hypothesis generators requiring an exact Xi identity."
        ),
        "proof_boundary": (
            "This artifact proves or reproduces source provenance, the "
            "radial and transpose identities, modular restrictions, "
            "rendering controls, the zero-coordinate bug, finite spectral "
            "indefiniteness, a Pell sublocus, and a termwise Mellin bridge. "
            "It does not prove a zero/curvature relation, prime-kernel "
            "reflection positivity, an Xi contact gap, signed Mobius "
            "cancellation, a one-turn phase budget, Lambda<=0, PF-infinity, "
            "RH, or a prize-level result."
        ),
    }
    arrays = {
        "source_data": source_data,
        "interpolation": interpolation_arrays,
        "zero": zero_arrays,
    }
    return artifact, arrays


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    args = parser.parse_args()

    artifact, _ = build_artifact(args.figure)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(make_note(artifact), encoding="utf-8")
    print(
        "built legacy prime-curvature symmetry audit: "
        f"{artifact['summary']['row_count']} rows, "
        f"{artifact['arithmetic']['sampled_prime_count']} sampled primes, "
        f"{artifact['arithmetic']['selected_pair_count']} retained pairs, "
        "one serial control figure"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
