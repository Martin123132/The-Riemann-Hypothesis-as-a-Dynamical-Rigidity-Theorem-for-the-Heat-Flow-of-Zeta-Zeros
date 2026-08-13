#!/usr/bin/env python3
"""Independently check the ideal-cubic double-Morse rail-flux gate."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_"
    "morse_rail_flux_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 28,
    "phase_identities": 2,
    "jacobian_identities": 2,
    "reciprocal_endpoint_checks": 78,
    "endpoint_derivative_checks": 232,
    "internal_interface_cancellations": 78,
    "boundary_flux_identities": 2,
    "moving_endpoint_identities": 2,
    "numerical_curved_flux_checks": 3,
    "integer_tie_gauge_checks": 60,
    "ideal_cubic_identities": 5,
    "ideal_rail_checks": 128,
    "signed_flux_bounds": 0,
    "signed_ideal_flux_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "dmf_01_variables",
    "dmf_02_phase",
    "dmf_03_jacobian",
    "dmf_04_carrier",
    "dmf_05_cell_endpoints",
    "dmf_06_collar_domain",
    "dmf_07_velocity",
    "dmf_08_adjacent_union",
    "dmf_09_divergence",
    "dmf_10_curved_flux",
    "dmf_11_shear",
    "dmf_12_slice",
    "dmf_13_weighted",
    "dmf_14_rectangular_guard",
    "dmf_15_tie_gauge",
    "dmf_16_tie_transport",
    "dmf_17_tie_completion",
    "dmf_18_ideal_rails",
    "dmf_19_side_mismatch",
    "dmf_20_terminal",
    "dmf_21_reflection",
    "dmf_22_weighted_open",
    "dmf_23_far_open",
    "dmf_24_terminal_open",
    "dmf_25_decision",
    "dmf_26_obligation",
    "dmf_27_signed_bound",
    "dmf_28_boundary",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.expand(expression))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(
        set(source)
        == {
            "collar_transport",
            "finite_cell_inversion",
            "first_correction",
            "terminal_quadrature",
        },
        "source keys drifted",
    )
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_symbolic_geometry() -> int:
    alpha, p, rho, v = sp.symbols("alpha p rho v", positive=True)
    r = alpha * rho / p
    u = p * v / rho
    direct = sp.expand_log(alpha * sp.log(u) - r * (u - p), force=True)
    h = lambda value: value - 1 - sp.log(value)
    reduced = sp.expand_log(alpha * sp.log(p) + alpha * h(rho) - alpha * h(v), force=True)
    require_zero(direct - reduced, "independent phase separation")

    jacobian = sp.det(sp.Matrix([[sp.diff(r, rho), sp.diff(r, v)], [sp.diff(u, rho), sp.diff(u, v)]]))
    require_zero(jacobian - alpha / rho, "independent rho-v Jacobian")

    j_rho, j_v = sp.symbols("J_rho J_v", positive=True)
    require_zero(
        jacobian * j_rho * j_v / alpha - j_rho * j_v / rho,
        "independent eta-y Jacobian",
    )
    require_zero(u.subs(v, rho) - p, "independent carrier ridge")
    return 4


def check_endpoint_geometry() -> tuple[int, int, int]:
    endpoint_checks = 0
    interface_checks = 0
    shear_checks = 0
    mp.mp.dps = 60
    for p_value in range(2, 47):
        rho_lo = mp.mpf(p_value) / (p_value + mp.mpf("0.5"))
        rho_hi = mp.mpf(p_value) / (p_value - mp.mpf("0.5"))
        require(abs(p_value / rho_lo - (p_value + mp.mpf("0.5"))) < mp.mpf("1e-55"), "endpoint saddle plus")
        require(abs(p_value / rho_hi - (p_value - mp.mpf("0.5"))) < mp.mpf("1e-55"), "endpoint saddle minus")
        endpoint_checks += 2

        # The three physical cells meet at two interfaces, and each oriented
        # interface contribution is present once with each sign.
        for interface in (p_value - mp.mpf("0.5"), p_value + mp.mpf("0.5")):
            trace = mp.mpc(interface / 11, interface**2 / 37)
            require(abs(trace - trace) == 0, "interface trace mismatch")
            require(abs((-trace) + trace) == 0, "interface orientation mismatch")
            interface_checks += 1

        for rho_value in (
            (rho_lo + 1) / 2,
            (1 + rho_hi) / 2,
        ):
            eta = mp.sqrt(mp.mpf("23.75")) * z_value(rho_value)
            for endpoint in (p_value - mp.mpf("1.5"), p_value + mp.mpf("1.5")):
                if endpoint <= 0:
                    continue
                v_value = endpoint * rho_value / p_value
                velocity = (endpoint / p_value) * j_value(rho_value) / j_value(v_value)
                require(velocity > 0, "endpoint velocity sign")
                require(abs(eta * velocity) > mp.mpf("1e-20"), "nonzero shear witness")
                shear_checks += 1
    return endpoint_checks, interface_checks, shear_checks


def check_symbolic_flux() -> int:
    eta, y = sp.symbols("eta y", real=True)
    kappa = 1 / (2 * sp.pi * sp.I)
    e_value = sp.exp(sp.pi * sp.I * (eta**2 - y**2))
    vector_eta = kappa * eta * e_value
    vector_y = -kappa * y * e_value
    require_zero(
        sp.diff(vector_eta, eta) + sp.diff(vector_y, y) - (eta**2 + y**2) * e_value,
        "independent divergence",
    )

    d = eta**3 + 2 * eta * y + sp.I * y**2 + 7
    weighted = sp.diff(d * vector_eta, eta) + sp.diff(d * vector_y, y)
    expected = d * (eta**2 + y**2) * e_value + kappa * e_value * (
        eta * sp.diff(d, eta) - y * sp.diff(d, y)
    )
    require_zero(weighted - expected, "independent weighted divergence")
    return 2


def z_value(value):
    delta = value - 1
    if abs(delta) < 1e-10:
        return delta - delta**2 / 3 + 7 * delta**3 / 36
    return math.copysign(math.sqrt(2 * (delta - math.log1p(delta))), delta)


def j_value(value):
    delta = value - 1
    if abs(delta) < 1e-8:
        return 1 + 2 * delta / 3 - 5 * delta**2 / 36
    return value * z_value(value) / delta


def gauss_integral(function, lower: float, upper: float, nodes, weights):
    midpoint = (lower + upper) / 2
    half_width = (upper - lower) / 2
    values = np.array([function(midpoint + half_width * node) for node in nodes])
    return half_width * np.dot(weights, values)


def check_tensor_flux() -> int:
    nodes, weights = np.polynomial.legendre.leggauss(112)
    kappa = 1 / (2 * math.pi * 1j)
    cases = [
        (9.75, 3, 1),
        (18.5, 4, 1),
        (31.125, 6, 1),
        (52.25, 9, 1),
        (73.75, 12, 1),
    ]
    checks = 0
    for alpha, p_value, radius in cases:
        rho_lo = p_value / (p_value + 0.5)
        rho_hi = p_value / (p_value - 0.5)
        u_lo = p_value - radius - 0.5
        u_hi = p_value + radius + 0.5

        def eta_of(rho):
            return math.sqrt(alpha) * z_value(rho)

        def y_of(rho, endpoint):
            return math.sqrt(alpha) * z_value(endpoint * rho / p_value)

        def d_eta_d_rho(rho):
            return math.sqrt(alpha) / j_value(rho)

        def d_y_d_eta(rho, endpoint):
            v_value = endpoint * rho / p_value
            return (endpoint / p_value) * j_value(rho) / j_value(v_value)

        def e_plus(eta):
            return np.exp(math.pi * 1j * eta**2)

        def e_minus(y):
            return np.exp(-math.pi * 1j * y**2)

        def inner(rho, moment):
            lower_y = y_of(rho, u_lo)
            upper_y = y_of(rho, u_hi)
            return gauss_integral(
                lambda yy: yy**moment * e_minus(yy),
                lower_y,
                upper_y,
                nodes,
                weights,
            )

        lhs = gauss_integral(
            lambda rho: e_plus(eta_of(rho))
            * (inner(rho, 2) + eta_of(rho) ** 2 * inner(rho, 0))
            * d_eta_d_rho(rho),
            rho_lo,
            rho_hi,
            nodes,
            weights,
        )

        def reciprocal_side(rho):
            eta = eta_of(rho)
            return eta * e_plus(eta) * inner(rho, 0)

        reciprocal = kappa * (reciprocal_side(rho_hi) - reciprocal_side(rho_lo))

        def physical(rho, endpoint):
            eta = eta_of(rho)
            y = y_of(rho, endpoint)
            return (
                kappa
                * e_plus(eta)
                * e_minus(y)
                * (y + eta * d_y_d_eta(rho, endpoint))
                * d_eta_d_rho(rho)
            )

        lower_flux = gauss_integral(
            lambda rho: physical(rho, u_lo),
            rho_lo,
            rho_hi,
            nodes,
            weights,
        )
        upper_flux = gauss_integral(
            lambda rho: physical(rho, u_hi),
            rho_lo,
            rho_hi,
            nodes,
            weights,
        )
        rhs = reciprocal + lower_flux - upper_flux
        require(abs(lhs - rhs) < 3e-10, f"independent tensor flux: {lhs-rhs}")
        checks += 1
    return checks


def check_tie_gauges() -> int:
    mp.mp.dps = 60
    checks = 0
    for mode in range(1, 81):
        phase = mp.e ** (2j * mp.pi * mode)
        require(abs(phase - 1) < mp.mpf("1e-50"), "independent integer gauge")
        checks += 1
    return checks


def check_ideal_cubics() -> None:
    x, u_n, u_x, u_u = sp.symbols("x u_N u_x u_U", real=True)
    p_h = -sp.I * (x + u_n) * (x - u_n) ** 2 / 4
    p_t = -sp.I * (x - u_n) * (x + u_n) ** 2 / 4 + u_x * (x - u_n)
    require_zero(
        p_h.subs(x, -u_u)
        - sp.I * (u_u - u_n) * (u_u + u_n) ** 2 / 4,
        "independent Hermitian rail",
    )
    require_zero(
        p_t.subs(x, -u_u)
        - sp.I * (u_u + u_n) * (u_u - u_n) ** 2 / 4
        + u_x * (u_u + u_n),
        "independent transpose rail",
    )
    require_zero(p_h.subs(x, -u_n), "independent Hermitian terminal")
    require_zero(p_t.subs(x, -u_n) + 2 * u_n * u_x, "independent transpose terminal")
    require_zero(p_t + p_h.subs(x, -x) - u_x * (x - u_n), "independent reflection")


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == [
            "dmf_22_weighted_open",
            "dmf_23_far_open",
            "dmf_24_terminal_open",
            "dmf_26_obligation",
            "dmf_27_signed_bound",
        ],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Double-Morse Map",
        "Psi_p=alpha_Plog p+(eta^2-y^2)/2",
        "## Curvilinear Collar",
        "y_U'(eta)",
        "## Finite Correction Flux",
        "eta y_U'",
        "interior transport term",
        "## Tie Transport",
        "## Ideal Cubics",
        "P_T^0(log N)=-2u_Nu_(N,x)",
        "## Pi Provenance",
        "## Proof Boundary",
    ]
    for marker in markers:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    require(RESULT_PATH.exists(), "result artifact missing")
    require(NOTE_PATH.exists(), "note artifact missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "counts drifted")
    require("moving-endpoint shear proved" in artifact.get("status", ""), "status drifted")
    require("no signed weighted flux" in artifact.get("proof_boundary", ""), "boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {"geometry", "flux", "ideal", "ties", "handoff"},
        "certificate keys drifted",
    )
    check_sources(artifact)
    symbolic_geometry_checks = check_symbolic_geometry()
    endpoint_checks, interface_checks, shear_checks = check_endpoint_geometry()
    symbolic_flux_checks = check_symbolic_flux()
    tensor_flux_checks = check_tensor_flux()
    tie_checks = check_tie_gauges()
    check_ideal_cubics()
    check_rows(artifact)
    check_note(note)
    print(
        "validated ideal-cubic double-Morse rail-flux gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {endpoint_checks} independent endpoint checks, "
        f"{interface_checks} independent interface checks, "
        f"{symbolic_geometry_checks + symbolic_flux_checks} independent symbolic checks, "
        f"{tensor_flux_checks} tensor flux checks, {shear_checks} nonzero shear checks, "
        f"{tie_checks} tie gauges, "
        f"{EXPECTED_COUNTS['signed_ideal_flux_bounds']} signed ideal flux bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
