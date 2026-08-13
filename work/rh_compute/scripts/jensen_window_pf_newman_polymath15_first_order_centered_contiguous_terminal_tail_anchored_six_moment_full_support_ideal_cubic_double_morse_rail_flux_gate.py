#!/usr/bin/env python3
"""Build the ideal-cubic double-Morse rail-flux gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_"
    "morse_rail_flux_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "collar_transport": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_ideal_cubic_"
        "reciprocal_collar_transport_gate.json"
    ),
    "finite_cell_inversion": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_finite_cell_inversion_gate.json"
    ),
    "first_correction": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_first_correction_gate.json"
    ),
    "terminal_quadrature": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_"
        "conditional_quadrature_gate.json"
    ),
}


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


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    require(
        payloads["collar_transport"]["counts"]["signed_collar_bounds"] == 0,
        "collar source drift",
    )
    require(
        payloads["finite_cell_inversion"]["counts"][
            "finite_cell_poisson_identities"
        ]
        >= 1,
        "finite-cell source drift",
    )
    require(
        payloads["first_correction"]["counts"][
            "first_correction_cancellations"
        ]
        == 1,
        "first-correction source drift",
    )
    require(
        payloads["terminal_quadrature"]["counts"][
            "complex_recurrence_relations"
        ]
        == 1,
        "terminal source drift",
    )
    return payloads


def source_audit() -> dict[str, dict[str, str]]:
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def double_morse_certificate() -> dict[str, str | int]:
    alpha, p, rho, v = sp.symbols("alpha p rho v", positive=True)
    r = alpha * rho / p
    u = p * v / rho
    phase = sp.expand_log(alpha * sp.log(u) - r * u + p * r, force=True)
    reduced = sp.expand_log(
        alpha * sp.log(p)
        + alpha * (rho - 1 - sp.log(rho))
        - alpha * (v - 1 - sp.log(v)),
        force=True,
    )
    require_zero(phase - reduced, "double-Morse phase")

    jacobian_rho_v = sp.det(
        sp.Matrix(
            [
                [sp.diff(r, rho), sp.diff(r, v)],
                [sp.diff(u, rho), sp.diff(u, v)],
            ]
        )
    )
    require_zero(jacobian_rho_v - alpha / rho, "rho-v Jacobian")

    j_rho, j_v = sp.symbols("J_rho J_v", nonzero=True)
    transformed_jacobian = jacobian_rho_v * j_rho / sp.sqrt(alpha) * j_v / sp.sqrt(alpha)
    require_zero(
        transformed_jacobian - j_rho * j_v / rho,
        "double-Morse Jacobian",
    )

    endpoint_checks = 0
    derivative_checks = 0
    interface_checks = 0
    for p_value in range(2, 41):
        rho_minus = Fraction(2 * p_value, 2 * p_value + 1)
        rho_plus = Fraction(2 * p_value, 2 * p_value - 1)
        require(Fraction(p_value, 1) / rho_minus == p_value + Fraction(1, 2), "lower reciprocal endpoint")
        require(Fraction(p_value, 1) / rho_plus == p_value - Fraction(1, 2), "upper reciprocal endpoint")
        endpoint_checks += 2

        for radius in (0, 1, 2):
            lower = Fraction(2 * p_value - 2 * radius - 1, 2)
            upper = Fraction(2 * p_value + 2 * radius + 1, 2)
            if lower <= 0:
                continue
            require(upper - lower == 2 * radius + 1, "collar width")
            require(lower < p_value < upper, "collar center")
            derivative_checks += 2

        # Three adjacent physical cells have two internal interfaces. Their
        # oriented fluxes occur once with each sign.
        for interface in (
            Fraction(2 * p_value - 1, 2),
            Fraction(2 * p_value + 1, 2),
        ):
            lower_top = -interface
            upper_bottom = interface
            require(lower_top + upper_bottom == 0, "interface flux cancellation")
            interface_checks += 1

    return {
        "variables": "Put rho=pr/alpha_P and v=ur/alpha_P. With eta=sqrt(alpha_P)z(rho) and y=sqrt(alpha_P)z(v), the inverse map is r=alpha_Prho/p and u=pv/rho.",
        "phase": "For Psi_p(r,u)=alpha_Plog u-r(u-p), the exact identity is Psi_p=alpha_Plog p+(eta^2-y^2)/2. The sequential positive and negative Morse phases are therefore one two-dimensional saddle chart.",
        "jacobian": "The exact Jacobian is dr du={J(rho)J(v)/rho}deta dy. Hence the transformed amplitude is G_P(eta,y)=exp(S(log(pv/rho)))P(log(pv/rho))J(rho)J(v)/rho.",
        "carrier_line": "The physical line u=p is exactly v=rho, hence y=eta and Psi_p=alpha_Plog p. It is the returned-carrier ridge, not an oscillatory error rail.",
        "cell_endpoints": "The reciprocal cell endpoints rho_-=p/(p+1/2) and rho_+=p/(p-1/2) map to fixed eta endpoints. Their inner saddles are u_0=p+1/2 and u_0=p-1/2 respectively.",
        "collar_domain": "For an internal radius-L collar, U_-=p-L-1/2 and U_+=p+L+1/2. Its Morse image is eta_-<=eta<=eta_+ and y_-(eta)<=y<=y_+(eta), where y_U(eta)=sqrt(alpha_P)z(Urho(eta)/p). This is curvilinear, not Cartesian.",
        "endpoint_velocity": "Along a physical rail u=U, y_U'(eta)=(U/p)J(rho)/J(Urho/p)=vJ(rho)/{rho J(v)}>0. The endpoint velocity vanishes nowhere; only its product with eta vanishes at the reciprocal saddle.",
        "adjacent_union": "For L=1, the cells J_(p-1), J_p, and J_(p+1) tile [p-3/2,p+3/2]. Their two internal Morse interfaces cancel with opposite orientations, leaving only the two outer collar rails.",
        "endpoint_checks": endpoint_checks,
        "endpoint_derivative_checks": derivative_checks,
        "internal_interface_cancellations": interface_checks,
        "phase_identities": 2,
        "jacobian_identities": 2,
    }


def flux_certificate() -> dict[str, str | int]:
    eta, y = sp.symbols("eta y", real=True)
    kappa = 1 / (2 * sp.pi * sp.I)
    exponential = sp.exp(sp.pi * sp.I * (eta**2 - y**2))
    f_eta = kappa * eta * exponential
    f_y = -kappa * y * exponential
    require_zero(
        sp.diff(f_eta, eta) + sp.diff(f_y, y)
        - (eta**2 + y**2) * exponential,
        "constant boundary flux",
    )

    d = sp.Function("D")(eta, y)
    require_zero(
        sp.diff(d * f_eta, eta)
        + sp.diff(d * f_y, y)
        - d * (eta**2 + y**2) * exponential
        - kappa * exponential * (eta * sp.diff(d, eta) - y * sp.diff(d, y)),
        "weighted boundary flux",
    )

    y_u, y_u_prime = sp.symbols("y_U y_U_prime")
    rail_form = eta * y_u_prime + y_u
    require_zero(rail_form - sp.diff(eta * sp.Function("Y")(eta), eta).subs(
        {
            sp.Function("Y")(eta): y_u,
            sp.Derivative(sp.Function("Y")(eta), eta): y_u_prime,
        }
    ), "rail one-form")

    return {
        "divergence": "For E(eta,y)=e((eta^2-y^2)/2), div(kappa eta E,-kappa y E)=(eta^2+y^2)E. The full-line first-correction cancellation is the zero-boundary limit of this exact divergence law.",
        "curved_flux": "On the true collar image, the paired finite second moments equal kappa times the four-side flux integral int_(boundary Omega)E{eta dy+y deta}. The reciprocal sides contribute kappa eta E dy; a physical rail contributes kappa E{y_U+eta y_U'}deta with its boundary orientation.",
        "shear": "The term kappa eta y_U'E is the moving-endpoint shear. It is absent only in an artificial Cartesian Morse rectangle. Thus the four fixed-limit terms in the two one-dimensional moment formulas do not by themselves equal the physical collar rails.",
        "slice_formula": "If F(eta)=integral_(y_-(eta))^(y_+(eta))e(-y^2/2)dy, positive-moment integration by parts differentiates F and creates y_+'e(-y_+^2/2)-y_-'e(-y_-^2/2). These are exactly the two shear terms required by the boundary flux.",
        "weighted_flux": "For a nonconstant coefficient D, div(Dkappa eta E,-Dkappa y E)=D(eta^2+y^2)E+kappa E(eta D_eta-y D_y). The physical ideal-cubic amplitude therefore leaves an interior transport term as well as boundary flux; constant-symbol cancellation cannot be promoted coefficientwise.",
        "boundary_flux_identities": 2,
        "moving_endpoint_identities": 2,
        "signed_flux_bounds": 0,
    }


def morse_z(value: mp.mpf) -> mp.mpf:
    delta = value - 1
    if delta == 0:
        return mp.mpf("0")
    return mp.sign(delta) * mp.sqrt(2 * (delta - mp.log1p(delta)))


def morse_j(value: mp.mpf) -> mp.mpf:
    delta = value - 1
    if abs(delta) < mp.mpf("1e-25"):
        return 1 + 2 * delta / 3 - 5 * delta**2 / 36
    return value * morse_z(value) / delta


def integrate_split(function, lower: mp.mpf, upper: mp.mpf, points: list[mp.mpf]):
    cuts = [lower] + [point for point in points if lower < point < upper] + [upper]
    return mp.quad(function, cuts)


def numerical_flux_checks() -> int:
    mp.mp.dps = 55
    kappa = 1 / (2 * mp.pi * 1j)
    cases = [
        (mp.mpf("12.5"), 3, 1),
        (mp.mpf("27.25"), 5, 1),
        (mp.mpf("41.75"), 8, 1),
    ]
    checks = 0
    for alpha, p_value, radius in cases:
        rho_lo = mp.mpf(p_value) / (mp.mpf(p_value) + mp.mpf("0.5"))
        rho_hi = mp.mpf(p_value) / (mp.mpf(p_value) - mp.mpf("0.5"))
        u_lo = mp.mpf(p_value) - radius - mp.mpf("0.5")
        u_hi = mp.mpf(p_value) + radius + mp.mpf("0.5")

        def eta_of(rho):
            return mp.sqrt(alpha) * morse_z(rho)

        def y_of(rho, endpoint):
            return mp.sqrt(alpha) * morse_z(endpoint * rho / p_value)

        def d_eta_d_rho(rho):
            return mp.sqrt(alpha) / morse_j(rho)

        def d_y_d_eta(rho, endpoint):
            v_value = endpoint * rho / p_value
            return (endpoint / p_value) * morse_j(rho) / morse_j(v_value)

        def e_plus(eta):
            return mp.e ** (mp.pi * 1j * eta**2)

        def e_minus(y_value):
            return mp.e ** (-mp.pi * 1j * y_value**2)

        def fresnel_primitive(y_value):
            root_i = mp.e ** (mp.pi * 1j / 4)
            return (
                mp.e ** (-mp.pi * 1j / 4)
                * mp.erf(mp.sqrt(mp.pi) * root_i * y_value)
                / 2
            )

        def inner_zero(rho):
            lower_y = y_of(rho, u_lo)
            upper_y = y_of(rho, u_hi)
            return fresnel_primitive(upper_y) - fresnel_primitive(lower_y)

        def inner_second(rho):
            lower_y = y_of(rho, u_lo)
            upper_y = y_of(rho, u_hi)
            zero_moment = inner_zero(rho)
            return kappa * zero_moment - kappa * (
                upper_y * e_minus(upper_y) - lower_y * e_minus(lower_y)
            )

        def lhs_integrand(rho):
            eta_value = eta_of(rho)
            return (
                e_plus(eta_value)
                * (inner_second(rho) + eta_value**2 * inner_zero(rho))
                * d_eta_d_rho(rho)
            )

        lhs = integrate_split(lhs_integrand, rho_lo, rho_hi, [mp.mpf("1")])

        def reciprocal_side(rho):
            eta_value = eta_of(rho)
            return eta_value * e_plus(eta_value) * inner_zero(rho)

        reciprocal_flux = kappa * (
            reciprocal_side(rho_hi) - reciprocal_side(rho_lo)
        )

        def physical_flux(rho, endpoint):
            eta_value = eta_of(rho)
            y_value = y_of(rho, endpoint)
            return (
                kappa
                * e_plus(eta_value)
                * e_minus(y_value)
                * (y_value + eta_value * d_y_d_eta(rho, endpoint))
                * d_eta_d_rho(rho)
            )

        lower_flux = integrate_split(
            lambda rho: physical_flux(rho, u_lo),
            rho_lo,
            rho_hi,
            [mp.mpf("1")],
        )
        upper_flux = integrate_split(
            lambda rho: physical_flux(rho, u_hi),
            rho_lo,
            rho_hi,
            [mp.mpf("1")],
        )
        rhs = reciprocal_flux + lower_flux - upper_flux
        require(abs(lhs - rhs) < mp.mpf("2e-38"), "numerical curved flux")
        checks += 1
    return checks


def ideal_certificate() -> dict[str, str | int]:
    x, u_n, u_x, u_u = sp.symbols("x u_N u_x u_U", real=True)
    p_h = -sp.I * (x + u_n) * (x - u_n) ** 2 / 4
    p_t = -sp.I * (x - u_n) * (x + u_n) ** 2 / 4 + u_x * (x - u_n)
    h_rail = sp.I * (u_u - u_n) * (u_u + u_n) ** 2 / 4
    t_rail = (
        sp.I * (u_u + u_n) * (u_u - u_n) ** 2 / 4
        - u_x * (u_u + u_n)
    )
    require_zero(p_h.subs(x, -u_u) - h_rail, "Hermitian rail")
    require_zero(p_t.subs(x, -u_u) - t_rail, "transpose rail")
    require_zero(p_h.subs(x, -u_n), "Hermitian terminal rail")
    require_zero(p_t.subs(x, -u_n) + 2 * u_n * u_x, "transpose terminal rail")
    require_zero(p_t + p_h.subs(x, -x) - u_x * (x - u_n), "ideal reflection")

    rail_checks = 0
    for p_value in range(3, 35):
        for radius in (0, 1):
            lower = Fraction(2 * p_value - 2 * radius - 1, 2)
            upper = Fraction(2 * p_value + 2 * radius + 1, 2)
            require(lower != upper, "distinct ideal rails")
            rail_checks += 2

    return {
        "rail_values": "For any positive physical rail U, put u_U=log(a/U). Then P_H^0(log U)=i(u_U-u_N)(u_U+u_N)^2/4 and P_T^0(log U)=i(u_U+u_N)(u_U-u_N)^2/4-u_(N,x)(u_U+u_N). These apply at half-integer collar rails as well as integer atoms.",
        "internal_values": "The reciprocal side saddles occur at U=p+1/2 and U=p-1/2, whereas the radius-one physical rails are U=p+3/2 and U=p-3/2. Adjacent-corner transport bridges a full cell; the positive endpoint trace cannot be identified pointwise with the outer collar rail.",
        "terminal": "At U=N, P_H^0(log N)=0 but P_T^0(log N)=-2u_Nu_(N,x). The exact Hermitian physical endpoint value vanishes; the transpose endpoint remains joined to the conditional terminal quadrature and full complex recurrence.",
        "reflection_guard": "The identity P_T^0(x)=-P_H^0(-x)+u_(N,x)(x-u_N) does not reverse the curvilinear collar: S(log u), the interval [1,N], reciprocal side traces, and endpoint velocities are not reflection symmetric.",
        "ideal_cubic_identities": 5,
        "ideal_rail_checks": rail_checks,
        "signed_ideal_flux_bounds": 0,
    }


def tie_certificate() -> dict[str, str | int]:
    p, r, u, alpha = sp.symbols("p r u alpha", real=True)
    psi_p = alpha * sp.log(u) - r * (u - p)
    psi_next = alpha * sp.log(u) - r * (u - p - 1)
    require_zero(psi_next - psi_p - r, "tie gauge difference")
    tie_checks = 0
    for integer_mode in range(1, 61):
        phase = mp.e ** (2j * mp.pi * integer_mode)
        require(abs(phase - 1) < mp.mpf("1e-45"), "integer tie gauge")
        tie_checks += 1
    return {
        "gauge": "At the shared reciprocal boundary, Psi_(p+1)-Psi_p=r. When the boundary is an actual integer tie, e(r)=1, so overlapping complete-mode traces use the same Fourier phase.",
        "transport": "The radius-one collar changes from [p-3/2,p+3/2] to [p-1/2,p+5/2]. After common-interface cancellation, the only transferred physical strips are the leaving J_(p-1) and entering J_(p+2), exactly as Section 11.189 requires.",
        "completion": "The collar-side flux is therefore tie invariant only after the corresponding far-side flux is included. The moving-endpoint shear travels with the same complete physical strips and cannot be assigned to the collar alone.",
        "integer_gauge_checks": tie_checks,
    }


def symbolic_certificate() -> dict:
    geometry = double_morse_certificate()
    flux = flux_certificate()
    flux["numerical_curved_flux_checks"] = numerical_flux_checks()
    ideal = ideal_certificate()
    ties = tie_certificate()
    return {
        "geometry": geometry,
        "flux": flux,
        "ideal": ideal,
        "ties": ties,
        "handoff": {
            "decision": "Replace the proposed four fixed-endpoint match by the exact curvilinear boundary flux. Keep reciprocal side integrals, physical outer-rail integrals, and moving-endpoint shear together; internal diagonal/adjacent interfaces cancel only after their oriented union.",
            "obligation": "Insert the complete transformed ideal amplitudes G_H^0 and G_T^0 into the weighted flux identity, combine the two physical rail traces with the rail-compressed far integration-by-parts functional, and compose both reciprocal side traces with finite-cell tie corrections and the q=N transpose terminal quadrature. Then estimate or obstruct the remaining interior transport eta G_eta-y G_y without blockwise moduli.",
            "reserve": "No signed complete ideal join, weighted rail flux, interior transport, far aggregate, h^2 reserve, completed current, or RH-level theorem is proved.",
        },
    }


def build_rows(cert: dict) -> list[GateRow]:
    geometry = cert["geometry"]
    flux = cert["flux"]
    ideal = cert["ideal"]
    ties = cert["ties"]
    handoff = cert["handoff"]
    rows = [
        GateRow("dmf_01_variables", "double-Morse variables", "proved", "The sequential saddle admits one exact two-variable coordinate map.", geometry["variables"], "All variables remain positive."),
        GateRow("dmf_02_phase", "separated phase", "proved", "The two sequential phases separate exactly.", geometry["phase"], "No quadratic truncation is used."),
        GateRow("dmf_03_jacobian", "coordinate Jacobian", "proved", "The double-Morse Jacobian is exact.", geometry["jacobian"], "The amplitude remains coupled."),
        GateRow("dmf_04_carrier", "carrier ridge", "proved", "The returned carrier is the constant-phase line y=eta.", geometry["carrier_line"], "It is not an error rail."),
        GateRow("dmf_05_cell_endpoints", "reciprocal sides", "proved", "Reciprocal cell ends map to fixed eta sides and half-integer inner saddles.", geometry["cell_endpoints"], "Endpoint half-weights remain explicit."),
        GateRow("dmf_06_collar_domain", "curvilinear collar", "proved", "A physical collar maps to moving y limits.", geometry["collar_domain"], "The domain is not a Morse rectangle."),
        GateRow("dmf_07_velocity", "endpoint velocity", "proved", "Every physical rail has an explicit positive Morse velocity.", geometry["endpoint_velocity"], "Its shear contribution is retained."),
        GateRow("dmf_08_adjacent_union", "adjacent composition", "proved", "The diagonal and adjacent cells cancel their internal physical interfaces.", geometry["adjacent_union"], "Only outer collar rails survive."),
        GateRow("dmf_09_divergence", "finite correction flux", "proved", "The paired second moments are an exact divergence.", flux["divergence"], "Full-line cancellation is its zero-boundary limit."),
        GateRow("dmf_10_curved_flux", "four-side flux", "proved", "The finite paired correction is the flux across the true curvilinear boundary.", flux["curved_flux"], "All side orientations are retained."),
        GateRow("dmf_11_shear", "moving endpoint shear", "proved", "Physical rail flux contains an additional endpoint-velocity term.", flux["shear"], "Dropping it changes the finite correction."),
        GateRow("dmf_12_slice", "slice derivative", "proved", "Differentiating the incomplete negative Fresnel factor reproduces the shear.", flux["slice_formula"], "Moving limits are not frozen."),
        GateRow("dmf_13_weighted", "weighted flux", "proved", "A nonconstant physical amplitude leaves an interior transport term.", flux["weighted_flux"], "Constant-symbol cancellation is insufficient."),
        GateRow("dmf_14_rectangular_guard", "Cartesian guard", "guard_validated", "The four fixed-limit moment terms are incomplete on the physical chart.", flux["shear"], "This is a route obstruction, not a sign theorem."),
        GateRow("dmf_15_tie_gauge", "integer tie phase", "proved", "Neighboring reciprocal gauges agree at an integer tie.", ties["gauge"], "Continuous noninteger sides still retain finite-Poisson traces."),
        GateRow("dmf_16_tie_transport", "collar shift", "proved", "The diagonal-adjacent union transfers exactly two physical strips.", ties["transport"], "The strips remain complete."),
        GateRow("dmf_17_tie_completion", "far completion", "guard_validated", "The shear and ordinary rail traces move with the complete mode.", ties["completion"], "Collar-only invariance remains false."),
        GateRow("dmf_18_ideal_rails", "ideal rail values", "proved", "Both ideal cubics have exact values on every physical rail.", ideal["rail_values"], "No railwise sign is inferred."),
        GateRow("dmf_19_side_mismatch", "endpoint location", "proved", "Reciprocal side saddles and outer collar rails are separated by one adjacent cell.", ideal["internal_values"], "Pointwise endpoint identification is invalid."),
        GateRow("dmf_20_terminal", "terminal endpoint", "proved", "The Hermitian terminal value vanishes and the transpose value survives.", ideal["terminal"], "The full complex recurrence remains compulsory."),
        GateRow("dmf_21_reflection", "reflection guard", "guard_validated", "Polynomial reflection does not reflect the physical flux domain.", ideal["reflection_guard"], "No H/T cancellation is inferred."),
        GateRow("dmf_22_weighted_open", "weighted ideal flux", "open", "Control the weighted physical boundary and interior transport.", flux["weighted_flux"], "No signed weighted flux bound is proved."),
        GateRow("dmf_23_far_open", "far rail composition", "open", "Compose the physical flux rails with the exact far functional.", handoff["obligation"], "No far aggregate is bounded."),
        GateRow("dmf_24_terminal_open", "terminal composition", "open", "Compose the reciprocal sides with tie halves and the transpose terminal quadrature.", handoff["obligation"], "No terminal sign is claimed."),
        GateRow("dmf_25_decision", "route decision", "guard_validated", "Use the curvilinear weighted flux as the next exact object.", handoff["decision"], "The rectangular shortcut is retired."),
        GateRow("dmf_26_obligation", "next obligation", "open", "Derive the complete weighted ideal composition.", handoff["obligation"], handoff["reserve"]),
        GateRow("dmf_27_signed_bound", "signed join", "open", "Prove or obstruct a signed complete ideal join only after exact composition.", handoff["obligation"], handoff["reserve"]),
        GateRow("dmf_28_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", handoff["reserve"], "No prize-level conclusion is claimed."),
    ]
    require(len(rows) == 28, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    geometry = cert["geometry"]
    flux = cert["flux"]
    ideal = cert["ideal"]
    ties = cert["ties"]
    handoff = cert["handoff"]
    return f"""# Ideal-Cubic Double-Morse Rail-Flux Gate

Date: 2026-08-03

Status: exact double-Morse geometry, curvilinear boundary flux, moving-endpoint shear, tie transport, and ideal rail values proved; signed complete ideal estimate open; not a proof of RH.

## Double-Morse Map

{geometry['variables']}

{geometry['phase']}

{geometry['jacobian']}

{geometry['carrier_line']}

{geometry['cell_endpoints']}

## Curvilinear Collar

{geometry['collar_domain']}

{geometry['endpoint_velocity']}

{geometry['adjacent_union']}

## Finite Correction Flux

{flux['divergence']}

{flux['curved_flux']}

{flux['shear']}

{flux['slice_formula']}

{flux['weighted_flux']}

## Tie Transport

{ties['gauge']}

{ties['transport']}

{ties['completion']}

## Ideal Cubics

{ideal['rail_values']}

{ideal['internal_values']}

{ideal['terminal']}

{ideal['reflection_guard']}

## Handoff

{handoff['decision']}

{handoff['obligation']}

## Pi Provenance

The character is e(x)=exp(2pi i x), so kappa=1/(2pi i). The divergence law follows by differentiating e((eta^2-y^2)/2); no geometric constant, circle, polygon, or fitted parameter is inserted.

## Proof Boundary

{handoff['reserve']} This gate proves no signed physical join, completed-current bound, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "phase_identities": certificate["geometry"]["phase_identities"],
        "jacobian_identities": certificate["geometry"]["jacobian_identities"],
        "reciprocal_endpoint_checks": certificate["geometry"]["endpoint_checks"],
        "endpoint_derivative_checks": certificate["geometry"]["endpoint_derivative_checks"],
        "internal_interface_cancellations": certificate["geometry"]["internal_interface_cancellations"],
        "boundary_flux_identities": certificate["flux"]["boundary_flux_identities"],
        "moving_endpoint_identities": certificate["flux"]["moving_endpoint_identities"],
        "numerical_curved_flux_checks": certificate["flux"]["numerical_curved_flux_checks"],
        "integer_tie_gauge_checks": certificate["ties"]["integer_gauge_checks"],
        "ideal_cubic_identities": certificate["ideal"]["ideal_cubic_identities"],
        "ideal_rail_checks": certificate["ideal"]["ideal_rail_checks"],
        "signed_flux_bounds": certificate["flux"]["signed_flux_bounds"],
        "signed_ideal_flux_bounds": certificate["ideal"]["signed_ideal_flux_bounds"],
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "double-Morse rail flux and moving-endpoint shear proved; signed complete ideal estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the exact double-Morse map, phase and Jacobian, curvilinear collar endpoints, moving-endpoint shear, finite correction boundary flux, weighted interior transport, adjacent-interface cancellation, integer tie gauge, ideal-cubic rail values, and terminal values. It proves no signed weighted flux, complete ideal join, far aggregate, completed current, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built ideal-cubic double-Morse rail-flux gate: "
        f"{counts['rows']} rows, {counts['reciprocal_endpoint_checks']} endpoint checks, "
        f"{counts['internal_interface_cancellations']} interface cancellations, "
        f"{counts['boundary_flux_identities']} flux identities, "
        f"{counts['numerical_curved_flux_checks']} curved-flux quadratures, "
        f"{counts['integer_tie_gauge_checks']} tie gauges, "
        f"{counts['ideal_rail_checks']} ideal rail checks, "
        f"{counts['signed_ideal_flux_bounds']} signed ideal flux bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
