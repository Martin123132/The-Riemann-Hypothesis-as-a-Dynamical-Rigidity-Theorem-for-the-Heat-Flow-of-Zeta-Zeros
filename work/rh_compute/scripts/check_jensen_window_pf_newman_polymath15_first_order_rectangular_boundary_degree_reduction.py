#!/usr/bin/env python3
"""Validate the first-order rectangular boundary-degree reduction."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_rectangular_boundary_degree_reduction"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "target_reconciliation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_target_reconciliation_gate.json"
    ),
    "boundary_degree": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
    ),
    "first_order_boundary": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_cofinal_boundary_reduction.json"
    ),
    "contact_index": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_first_jet_winding_gate.json"
    ),
    "oriented_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_oriented_successor_winding_reduction.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_carrier_kernel_abel_prefix_reduction.json"
    ),
    "abel_shear": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_abel_scalar_shear_flux_reduction.json"
    ),
    "oscillatory_splice": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_oscillatory_spliced_outer_collar_reduction.json"
    ),
    "connector_cap": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_dominant_ray_connector_phase_cap.json"
    ),
}
EXPECTED_IDS = [
    "rbd_01_current_main",
    "rbd_02_error_rectangle",
    "rbd_03_box_gauge",
    "rbd_04_rectangular_rouche",
    "rbd_05_radial_overstrength",
    "rbd_06_error_whitening",
    "rbd_07_winding_invariance",
    "rbd_08_boundary_abel",
    "rbd_09_zero_fibre",
    "rbd_10_adjacent_box",
    "rbd_11_adjacent_homotopy",
    "rbd_12_normalizer",
    "rbd_13_contact_index",
    "rbd_14_zero_degree",
    "rbd_15_integer_trap",
    "rbd_16_ray_schedule",
    "rbd_17_hard_bottom_interval",
    "rbd_18_outer_closure",
    "rbd_19_uniform_route",
    "rbd_20_boundary_route",
    "rbd_21_boundary_guard",
    "rbd_22_handoff",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate_constants(issues: list[str]) -> None:
    value = 100_000
    derivative = 200_000
    if value**2 + derivative**2 != 50_000_000_000:
        issues.append("rectangular radial constant failed")
    if Fraction(20_000, value) != Fraction(1, 5):
        issues.append("adjacent value ratio failed")
    if Fraction(30_000, derivative) != Fraction(3, 20):
        issues.append("adjacent derivative ratio failed")

    x = Fraction(0)
    u = Fraction(21, 10)
    if max(abs(x), abs(u) / 2) <= 1:
        issues.append("box witness did not clear unit square")
    if x**2 + u**2 >= 5:
        issues.append("box witness did not separate radial target")


def validate_whitening(issues: list[str]) -> None:
    b_0, b_1, s = sp.symbols("B_0 B_1 s", positive=True)
    matrix = sp.diag(1 / b_0, 1 / b_1)
    if sp.simplify(matrix.det() - 1 / (b_0 * b_1)) != 0:
        issues.append("whitening determinant failed")
    homotopy = (1 - s) * sp.eye(2) + s * matrix
    expected = (
        (1 - s) + s / b_0
    ) * (
        (1 - s) + s / b_1
    )
    if sp.simplify(homotopy.det() - expected) != 0:
        issues.append("positive diagonal homotopy determinant failed")


def validate_bottom_interval(issues: list[str]) -> None:
    j, ell, c_epsilon = sp.symbols(
        "j L c_epsilon", nonnegative=True
    )
    t_j = sp.Rational(25, 1) / (100 + j)
    q = 2 * t_j * ell**2
    lower = sp.sqrt((100 + j) / 50)
    upper = c_epsilon * (100 + j) / 25
    if sp.simplify(q.subs(ell, lower) - 1) != 0:
        issues.append("q=1 lower endpoint failed")
    if sp.simplify((t_j * ell).subs(ell, upper) - c_epsilon) != 0:
        issues.append("low-c upper endpoint failed")


def validate_sources(
    stored: dict, sources: dict[str, dict], issues: list[str]
) -> None:
    expected_hashes = {
        key: file_hash(path)
        for key, path in SOURCES.items()
        if path.is_file()
    }
    if stored.get("source_audit", {}).get("source_sha256") != expected_hashes:
        issues.append("source hashes drifted")

    target = sources["target_reconciliation"].get("exact", {})
    if target.get("box_optimal_target", {}).get("exact_target") != (
        "G_L(X,U)>1, equivalently "
        "|X|<=delta_0(L) => |U|>gamma_0(L)"
    ):
        issues.append("source box target drifted")

    boundary = sources["first_order_boundary"].get("exact", {})
    if "1300000000*exp(-5L/2)" not in boundary.get(
        "adjacent_vector", ""
    ):
        issues.append("source adjacent budget drifted")

    contact = sources["contact_index"].get("exact", {})
    if "+floor(m/2)" not in contact.get("multiplicity_index", ""):
        issues.append("source contact index drifted")

    successor = sources["oriented_successor"].get("exact", {})
    if "wind(proxy_j)<1" not in successor.get(
        "successor", {}
    ).get("one_sided_trap", ""):
        issues.append("source integer trap drifted")

    prefix = sources["abel_prefix"].get("exact", {})
    if "mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X" not in (
        prefix.get("contact_scalar", "")
    ):
        issues.append("source Abel shear drifted")

    splice = sources["oscillatory_splice"].get("exact", {})
    if "0<tL<=c_*+epsilon" not in splice.get("open_outer_target", ""):
        issues.append("source low-c splice drifted")

    connector = sources["connector_cap"].get("exact", {})
    if "H_j<3*pi/2" not in connector.get(
        "reduced_horizontal_target", ""
    ):
        issues.append("source horizontal phase target drifted")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "reduction result", issues)
    sources = {
        key: load_json(path, f"{key} source", issues)
        for key, path in SOURCES.items()
    }
    if issues:
        return issues

    if stored.get("kind") != STEM:
        issues.append("result kind drifted")
    if stored.get("date") != "2026-07-29":
        issues.append("result date drifted")
    if [row.get("id") for row in stored.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    expected_summary = {
        "rows": 22,
        "error_box_coordinates": 2,
        "rectangular_homotopies": 2,
        "positive_index_degree_transfers": 1,
        "boundary_only_abel_targets": 1,
        "route_branches": 2,
        "nonpromotion_guards": 2,
        "pointwise_contact_exclusions": 0,
        "cofinal_descendant_theorem": False,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    exact = stored.get("exact", {})
    if exact.get("main_gauge", {}).get("complex_half") != (
        "M_L(J)=G_L(X,U)"
    ):
        issues.append("stored complex-half gauge drifted")
    if exact.get("adjacent_chart", {}).get("gauge") != (
        "M_L(Delta J)<1/5"
    ):
        issues.append("stored adjacent box drifted")
    if exact.get("degree_transfer", {}).get("successor_trap") != (
        "0<=kappa_j=wind(V_J)<1 => kappa_j=0"
    ):
        issues.append("stored integer trap drifted")
    if exact.get("ray_bottom", {}).get("hard_interval") != (
        "max(B_epsilon,sqrt((100+j)/50))<=L"
        "<=min(L_j,(c_*+epsilon)*(100+j)/25)"
    ):
        issues.append("stored hard bottom interval drifted")
    if "one-dimensional hard-bottom Abel theorem" not in (
        exact.get("route_fork", {}).get("recommended_next", "")
    ):
        issues.append("stored route decision drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "does not prove the Xi Abel gap",
        "horizontal phase bound",
        "contact exclusion",
        "RH",
    ):
        if required not in boundary:
            issues.append(f"proof boundary missing: {required}")

    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
    else:
        note = NOTE.read_text(encoding="utf-8")
        for required in (
            "# First-Order Rectangular Boundary-Degree Reduction",
            "## Error-Matched Boundary Geometry",
            "## Hard Bottom Interval",
            "## Route Fork",
            "0 contact exclusions",
        ):
            if required not in note:
                issues.append(f"note marker missing: {required}")

    validate_constants(issues)
    validate_whitening(issues)
    validate_bottom_interval(issues)
    validate_sources(stored, sources, issues)
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"ERROR: {issue}")
    if issues:
        return 1
    print(
        "validated Newman first-order rectangular boundary-degree reduction: "
        "22 rows, 2 error-box coordinates, 2 rectangular homotopies, "
        "1 positive-index degree transfer, 1 boundary-only Abel target, "
        "2 route branches, 2 nonpromotion guards, 0 contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
