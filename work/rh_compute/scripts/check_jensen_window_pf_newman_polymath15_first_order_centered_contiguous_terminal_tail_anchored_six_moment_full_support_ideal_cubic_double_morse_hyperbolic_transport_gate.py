#!/usr/bin/env python3
"""Independently check the double-Morse hyperbolic-transport gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_"
    "morse_hyperbolic_transport_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 26,
    "profile_identities": 4,
    "profile_sign_checks": 16,
    "transport_identities": 4,
    "numerical_transport_checks": 6,
    "hyperbolic_identities": 3,
    "transverse_kernel_checks": 8,
    "ideal_transport_identities": 5,
    "signed_transport_bounds": 0,
    "signed_ideal_transport_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "dht_01_profiles",
    "dht_02_jacobian",
    "dht_03_series",
    "dht_04_profile_sign",
    "dht_05_operator",
    "dht_06_closed_form",
    "dht_07_saddle_zero",
    "dht_08_ridge",
    "dht_09_hyperbolic",
    "dht_10_boost",
    "dht_11_kernel",
    "dht_12_route",
    "dht_13_ideal_derivatives",
    "dht_14_terminal",
    "dht_15_ridge_guard",
    "dht_16_constant_witness",
    "dht_17_ridge_open",
    "dht_18_central_open",
    "dht_19_outer_open",
    "dht_20_rails_open",
    "dht_21_decision",
    "dht_22_obligation",
    "dht_23_signed_join",
    "dht_24_corrections",
    "dht_25_residual",
    "dht_26_boundary",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(
        set(source) == {"rail_flux", "collar_transport", "first_correction"},
        "source keys drifted",
    )
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def r_profile(value):
    delta = value - 1
    if abs(delta) < mp.mpf("1e-25"):
        return delta - 2 * delta**2 / 3 + delta**3 / 2 - 2 * delta**4 / 5
    return 2 * (delta - mp.log1p(delta)) / delta


def b_profile(value):
    delta = value - 1
    if abs(delta) < mp.mpf("1e-25"):
        return 1 - 2 * delta / 3 + delta**2 / 2 - 2 * delta**3 / 5
    return 2 * (delta - mp.log1p(delta)) / delta**2


def z_value(value):
    delta = value - 1
    if delta == 0:
        return mp.mpf("0")
    return mp.sign(delta) * mp.sqrt(2 * (delta - mp.log1p(delta)))


def j_value(value):
    delta = value - 1
    if abs(delta) < mp.mpf("1e-25"):
        return 1 + 2 * delta / 3 - 5 * delta**2 / 36
    return value * z_value(value) / delta


def check_profiles() -> tuple[int, int]:
    mp.mp.dps = 70
    profile_checks = 0
    ridge_witnesses = 0
    values = [
        mp.mpf("0.12"),
        mp.mpf("0.3"),
        mp.mpf("0.55"),
        mp.mpf("0.78"),
        mp.mpf("0.91"),
        mp.mpf("0.985"),
        mp.mpf("1.015"),
        mp.mpf("1.09"),
        mp.mpf("1.24"),
        mp.mpf("1.7"),
        mp.mpf("2.6"),
        mp.mpf("5.0"),
    ]
    for value in values:
        r_value = r_profile(value)
        b_value = b_profile(value)
        require(mp.sign(r_value) == mp.sign(value - 1), "independent R sign")
        require(b_value > 0, "independent B sign")
        profile_checks += 2

    for value in [
        mp.mpf("0.25"),
        mp.mpf("0.6"),
        mp.mpf("0.85"),
        mp.mpf("0.97"),
        mp.mpf("0.995"),
        mp.mpf("1.005"),
        mp.mpf("1.03"),
        mp.mpf("1.15"),
        mp.mpf("1.5"),
        mp.mpf("3.0"),
    ]:
        # Constant source P=1, S=0 on v=rho has transport
        # -J(rho)^2 R(rho)/rho.
        ridge = -j_value(value) ** 2 * r_profile(value) / value
        require(mp.sign(ridge) == -mp.sign(value - 1), "ridge orientation witness")
        ridge_witnesses += 1
    return profile_checks, ridge_witnesses


def check_symbolic_profiles() -> int:
    w = sp.symbols("w", positive=True)
    delta = w - 1
    h = delta - sp.log(w)
    r_value = 2 * h / delta
    b_value = 2 * h / delta**2
    k_value = w * r_value
    log_j_prime = 1 / w + delta / (2 * w * h) - 1 / delta
    require_zero(k_value * log_j_prime - (1 - b_value), "independent profile identity")
    require_zero(k_value / w - r_value, "independent radial identity")
    return 2


def check_transport() -> int:
    mp.mp.dps = 65
    p_value = mp.mpf("9.5")
    t_value = mp.mpf("0.007")
    sigma = mp.mpf("0.5007")

    def polynomial(x):
        return (
            mp.mpf("0.8")
            - mp.mpf("0.17") * 1j * x
            + (mp.mpf("0.04") + mp.mpf("0.025") * 1j) * x**2
            - mp.mpf("0.006") * x**3
        )

    def polynomial_prime(x):
        return (
            -mp.mpf("0.17") * 1j
            + 2 * (mp.mpf("0.04") + mp.mpf("0.025") * 1j) * x
            - mp.mpf("0.018") * x**2
        )

    def source(rho, v):
        x = mp.log(p_value * v / rho)
        s_value = t_value * x**2 / 4 - sigma * x
        return mp.e**s_value * polynomial(x) * j_value(rho) * j_value(v) / rho

    cases = [
        (mp.mpf("0.68"), mp.mpf("0.74")),
        (mp.mpf("0.81"), mp.mpf("1.16")),
        (mp.mpf("0.96"), mp.mpf("0.88")),
        (mp.mpf("1.04"), mp.mpf("1.13")),
        (mp.mpf("1.21"), mp.mpf("0.79")),
        (mp.mpf("1.35"), mp.mpf("1.29")),
        (mp.mpf("0.76"), mp.mpf("1.38")),
        (mp.mpf("1.42"), mp.mpf("0.71")),
    ]
    checks = 0
    for rho, v in cases:
        direct = rho * r_profile(rho) * mp.diff(lambda rr: source(rr, v), rho)
        direct -= v * r_profile(v) * mp.diff(lambda vv: source(rho, vv), v)
        x = mp.log(p_value * v / rho)
        s_value = t_value * x**2 / 4 - sigma * x
        g_value = t_value * x / 2 - sigma
        closed = (
            mp.e**s_value
            * j_value(rho)
            * j_value(v)
            / rho
            * (
                (b_profile(v) - b_profile(rho) - r_profile(rho)) * polynomial(x)
                - (r_profile(rho) + r_profile(v))
                * (polynomial_prime(x) + g_value * polynomial(x))
            )
        )
        require(abs(direct - closed) < mp.mpf("2e-52"), "independent transport")
        checks += 1
    return checks


def check_hyperbolic() -> int:
    eta, y, s, d = sp.symbols("eta y s d", real=True)
    eta_sd = (s + d) / sp.sqrt(2)
    y_sd = (s - d) / sp.sqrt(2)
    require_zero((eta_sd**2 - y_sd**2) / 2 - s * d, "independent sd phase")
    jacobian = sp.det(
        sp.Matrix(
            [
                [sp.diff(eta_sd, s), sp.diff(eta_sd, d)],
                [sp.diff(y_sd, s), sp.diff(y_sd, d)],
            ]
        )
    )
    require(abs(jacobian) == 1, "independent sd Jacobian")
    f = sp.Function("f")(s, d)
    operator = eta_sd * (sp.diff(f, s) + sp.diff(f, d)) / sp.sqrt(2)
    operator -= y_sd * (sp.diff(f, s) - sp.diff(f, d)) / sp.sqrt(2)
    require_zero(operator - d * sp.diff(f, s) - s * sp.diff(f, d), "independent boost")
    return 3


def check_transverse_kernels() -> int:
    mp.mp.dps = 65
    kappa = 1 / (2 * mp.pi * 1j)
    checks = 0
    for s_value in [
        mp.mpf("-3.1"),
        mp.mpf("-1.4"),
        mp.mpf("-0.35"),
        mp.mpf("0.22"),
        mp.mpf("1.1"),
        mp.mpf("2.7"),
    ]:
        for lower, upper in [
            (mp.mpf("-1.3"), mp.mpf("0.6")),
            (mp.mpf("0.15"), mp.mpf("1.75")),
        ]:
            direct = mp.quad(
                lambda dd: mp.e ** (2j * mp.pi * s_value * dd),
                [lower, upper],
            )
            closed = kappa * (
                mp.e ** (2j * mp.pi * s_value * upper)
                - mp.e ** (2j * mp.pi * s_value * lower)
            ) / s_value
            require(abs(direct - closed) < mp.mpf("2e-55"), "independent transverse kernel")
            checks += 1
    return checks


def check_ideal() -> int:
    x, u_n, u_x, g = sp.symbols("x u_N u_x g", real=True)
    p_h = -sp.I * (x + u_n) * (x - u_n) ** 2 / 4
    p_t = -sp.I * (x - u_n) * (x + u_n) ** 2 / 4 + u_x * (x - u_n)
    h_prime = -sp.I * (x - u_n) * (3 * x + u_n) / 4
    t_prime = -sp.I * (x + u_n) * (3 * x - u_n) / 4 + u_x
    require_zero(sp.diff(p_h, x) - h_prime, "independent H derivative")
    require_zero(sp.diff(p_t, x) - t_prime, "independent T derivative")
    require_zero(h_prime.subs(x, -u_n) + sp.I * u_n**2, "independent H terminal derivative")
    require_zero(t_prime.subs(x, -u_n) - u_x, "independent T terminal derivative")
    require_zero(
        (t_prime + g * p_t).subs(x, -u_n) - u_x * (1 - 2 * g * u_n),
        "independent weighted T terminal",
    )
    return 5


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == [
            "dht_17_ridge_open",
            "dht_18_central_open",
            "dht_19_outer_open",
            "dht_20_rails_open",
            "dht_22_obligation",
            "dht_23_signed_join",
            "dht_24_corrections",
            "dht_25_residual",
        ],
        "open rows drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Transport Profiles",
        "mathfrak R(1)=0",
        "## Exact Interior Transport",
        "mathcal L G_P",
        "carrier ridge",
        "## Hyperbolic Coordinates",
        "(eta^2-y^2)/2=sd",
        "1/s kernel",
        "## Ideal Cubics",
        "(P_H^0)'=-i",
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
    require("carrier-ridge reduction proved" in artifact.get("status", ""), "status drifted")
    require("no signed ridge-completed ideal join" in artifact.get("proof_boundary", ""), "boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {"profiles", "transport", "hyperbolic", "ideal", "handoff"},
        "certificate keys drifted",
    )
    check_sources(artifact)
    symbolic_profile_checks = check_symbolic_profiles()
    profile_checks, ridge_witnesses = check_profiles()
    transport_checks = check_transport()
    hyperbolic_checks = check_hyperbolic()
    kernel_checks = check_transverse_kernels()
    ideal_checks = check_ideal()
    check_rows(artifact)
    check_note(note)
    require(symbolic_profile_checks == 2, "symbolic profile count")
    print(
        "validated ideal-cubic double-Morse hyperbolic-transport gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {profile_checks} independent profile checks, "
        f"{transport_checks} independent transport checks, "
        f"{hyperbolic_checks} independent hyperbolic checks, "
        f"{kernel_checks} independent transverse kernels, {ideal_checks} ideal checks, "
        f"{ridge_witnesses} ridge sign witnesses, "
        f"{EXPECTED_COUNTS['signed_ideal_transport_bounds']} signed ideal transport bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
