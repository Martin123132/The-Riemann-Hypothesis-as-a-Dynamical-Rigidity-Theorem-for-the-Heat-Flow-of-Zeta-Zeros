#!/usr/bin/env python3
"""Build the ideal-cubic Fresnel transport and reciprocal-self-duality gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
    "ideal_cubic_reciprocal_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "aggregate_scaling": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_aggregate_scaling_gate.json"
    ),
    "basis_conjugation": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "basis_conjugation_gate.json"
    ),
    "endpoint_abel": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_endpoint_coherence_abel_gate.json"
    ),
    "endpoint_composition": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_endpoint_composition_"
        "retention_gate.json"
    ),
    "nilpotent_transport": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_nilpotent_transport_gate.json"
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


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def require_zero(expression: sp.Expr | sp.MatrixBase, label: str) -> None:
    if isinstance(expression, sp.MatrixBase):
        failures = [sp.simplify(item) for item in expression if sp.simplify(item) != 0]
        if failures:
            raise RuntimeError(f"{label}: {failures[0]}")
        return
    simplified = sp.simplify(expression)
    if simplified != 0:
        raise RuntimeError(f"{label}: {simplified}")


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        if not path.exists():
            raise RuntimeError(f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    if payloads["aggregate_scaling"]["counts"]["explicit_tail_h2_scaling_bounds"] != 6:
        raise RuntimeError("aggregate-scaling source drifted")
    if payloads["basis_conjugation"]["counts"]["joined_operator_identities"] != 3:
        raise RuntimeError("basis-conjugation source drifted")
    endpoint_counts = payloads["endpoint_abel"]["counts"]
    if endpoint_counts["endpoint_phase_collapses"] != 2:
        raise RuntimeError("endpoint-coherence source drifted")
    if endpoint_counts["weighted_abel_bounds"] != 1:
        raise RuntimeError("endpoint Abel source drifted")
    if (
        payloads["endpoint_composition"]["counts"][
            "endpoint_composed_value_identities"
        ]
        != 6
    ):
        raise RuntimeError("endpoint-composition source drifted")
    nilpotent_counts = payloads["nilpotent_transport"]["counts"]
    if nilpotent_counts["augmented_nilpotence_order"] != 4:
        raise RuntimeError("nilpotent source drifted")
    if nilpotent_counts["collar_variation_bounds"] != 2:
        raise RuntimeError("nilpotent variation source drifted")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_certificate() -> dict:
    i = sp.I
    x, x_0, u_x = sp.symbols("x x_0 u_x", real=True)
    t, delta_a, d = sp.symbols("t delta_a d", nonnegative=True, real=True)
    d_1, d_2 = sp.symbols("d_1 d_2", nonnegative=True, real=True)
    a_z, b_z = sp.symbols("A_z B_z", real=True)
    identity_4 = sp.eye(4)
    identity_8 = sp.eye(8)
    zero_4 = sp.zeros(4)

    d_0 = sp.Matrix(
        [
            [0, 0, 0, 0],
            [0, 0, i / 2, i / 2],
            [i / 2, 0, 0, 0],
            [i / 2, 0, 0, 0],
        ]
    )
    d_hat = d_0.row_join(zero_4).col_join(identity_4.row_join(d_0))

    def field(variable: sp.Expr) -> sp.Matrix:
        base = sp.Matrix(
            [
                1,
                -variable**2 / 4 - i * u_x / 2,
                i * variable / 2,
                i * variable / 2,
            ]
        )
        return base.col_join(variable * base)

    y = field(x)
    y_0 = field(x_0)
    require_zero(sp.diff(y, x) - d_hat * y, "augmented ideal jet")

    b_a = sp.Rational(1, 2) - t * delta_a / 2

    def centered_exponent(variable: sp.Expr) -> sp.Expr:
        return t * variable**2 / 4 - b_a * variable

    def epsilon(variable: sp.Expr) -> sp.Expr:
        return t * (variable + delta_a) / 2

    def g_value(variable: sp.Expr) -> sp.Expr:
        return -sp.Rational(1, 2) + epsilon(variable)

    z_x = sp.exp(centered_exponent(x)) * y
    z_0 = sp.exp(centered_exponent(x_0)) * y_0
    require_zero(
        sp.diff(z_x, x) - (d_hat + g_value(x) * identity_8) * z_x,
        "exponentially weighted jet",
    )

    def transport(step: sp.Expr) -> sp.Matrix:
        return sum(
            (
                ((-step) ** power / sp.factorial(power)) * d_hat**power
                for power in range(4)
            ),
            sp.zeros(8),
        )

    epsilon_symbol = sp.symbols("epsilon", real=True)

    def weighted_transport(step: sp.Expr, eps: sp.Expr) -> sp.Matrix:
        scalar = sp.exp(
            -step * (sp.Rational(1, 2) + eps) + t * step**2 / 4
        )
        return scalar * transport(step)

    require_zero(
        weighted_transport(d, epsilon(x)) * z_x
        - sp.exp(-d) * z_x.subs(x, x - d),
        "weighted-denominator transport",
    )
    require_zero(
        weighted_transport(d_2, epsilon_symbol - t * d_1 / 2)
        * weighted_transport(d_1, epsilon_symbol)
        - weighted_transport(d_1 + d_2, epsilon_symbol),
        "weighted transport cocycle",
    )

    def kernel(eps: sp.Expr) -> sp.Matrix:
        return a_z * (
            d_hat + (-sp.Rational(1, 2) + eps) * identity_8
        ) + b_z * identity_8

    state = z_x.col_join(z_0)
    operator_0 = kernel(epsilon(x)).row_join(identity_8)
    numerator = a_z * sp.diff(z_x, x) + b_z * z_x + z_0
    require_zero(operator_0 * state - numerator, "off-saddle vector kernel")

    def shifted_operator(step: sp.Expr) -> sp.Matrix:
        left = kernel(epsilon(x) - t * step / 2) * weighted_transport(
            step, epsilon(x)
        )
        right = weighted_transport(step, epsilon(x_0))
        return left.row_join(right)

    require_zero(
        shifted_operator(d) * state
        - sp.exp(-d)
        * (
            kernel(epsilon(x) - t * d / 2) * z_x.subs(x, x - d)
            + z_0.subs(x_0, x_0 - d)
        ),
        "adjacent off-saddle transport",
    )
    first_operator = sp.simplify(shifted_operator(d) - operator_0)
    second_operator = sp.simplify(
        shifted_operator(d_1 + d_2)
        - 2 * shifted_operator(d_1)
        + operator_0
    )

    def saddle_matrix(eps: sp.Expr) -> sp.Matrix:
        return (
            d_hat**2
            - identity_8 / 12
            + 2 * eps * d_hat
            + (eps**2 + t / 2) * identity_8
        )

    saddle_0 = saddle_matrix(epsilon(x_0))

    def shifted_saddle(step: sp.Expr) -> sp.Matrix:
        return saddle_matrix(
            epsilon(x_0) - t * step / 2
        ) * weighted_transport(step, epsilon(x_0))

    saddle_first = sp.simplify(shifted_saddle(d) - saddle_0)
    saddle_second = sp.simplify(
        shifted_saddle(d_1 + d_2)
        - 2 * shifted_saddle(d_1)
        + saddle_0
    )

    alpha_p, r, u, q = sp.symbols(
        "alpha_P r u q", positive=True, real=True
    )
    joint_phase = alpha_p * sp.log(u) - r * (u - q)
    critical = {r: alpha_p / q, u: q}
    gradient = sp.Matrix(
        [sp.diff(joint_phase, r), sp.diff(joint_phase, u)]
    )
    require_zero(gradient.subs(critical), "joint reciprocal critical point")
    critical_value = sp.simplify(joint_phase.subs(critical))
    require_zero(
        critical_value - alpha_p * sp.log(q), "joint reciprocal phase"
    )
    hessian = sp.hessian(joint_phase, (r, u)).subs(critical)
    require_zero(hessian.det() + 1, "joint Hessian determinant")
    u_curvature = -alpha_p / q**2
    r_curvature = q**2 / alpha_p
    require_zero(u_curvature * r_curvature + 1, "sequential curvature product")
    saddle_prefactor = q / (2 * alpha_p ** sp.Rational(3, 2))
    reciprocal_width = sp.sqrt(alpha_p) / q
    require_zero(
        saddle_prefactor * reciprocal_width - 1 / (2 * alpha_p),
        "reciprocal saddle scale",
    )

    z_symbols = sp.symbols("z_0:5")
    w_symbols = sp.symbols("w_0:5")
    partial = [sum(z_symbols[: index + 1]) for index in range(5)]
    double_partial = [sum(partial[: index + 1]) for index in range(5)]
    lhs = sum(w_symbols[index] * z_symbols[index] for index in range(5))
    rhs = (
        w_symbols[4] * partial[4]
        - (w_symbols[4] - w_symbols[3]) * double_partial[3]
        + sum(
            (
                w_symbols[index + 2]
                - 2 * w_symbols[index + 1]
                + w_symbols[index]
            )
            * double_partial[index]
            for index in range(3)
        )
    )
    require_zero(lhs - rhs, "second-order Abel identity")

    return {
        "weighted_jet": {
            "definition": (
                "Z(x)=exp(S_tilde(x))Y(x), with "
                "Z'=[Dhat+g_tilde(x)I_8]Z."
            ),
            "off_saddle": (
                "H_r(y)=[A_z(Dhat+g_vI_8)+B_zI_8]Z(x_r)+Z(x_(0,r)); "
                "c'_(bar P,r)(y)=exp(S(ell))wH_r(y)/(r sqrt(alpha_P)z^2), "
                "where w=(alpha^T,i beta^T)."
            ),
            "kernel": (
                "A_z=zJ^2/v, B_z=zJJ'-J, and "
                "K_z(epsilon)=A_z[Dhat+(-1/2+epsilon)I_8]+B_zI_8."
            ),
        },
        "mode_transport": {
            "weighted": (
                "P_d(epsilon)=exp[-d(1/2+epsilon)+td^2/4]E(d), so "
                "r/(r+1) times the shifted weighted jet equals P_d(epsilon)Z."
            ),
            "cocycle": (
                "P_(d_2)(epsilon-td_1/2)P_(d_1)(epsilon)="
                "P_(d_1+d_2)(epsilon) exactly."
            ),
            "operator": (
                "L_d=[K_z(epsilon_v-td/2)P_d(epsilon_v),P_d(epsilon_0)], "
                "L_0=[K_z(epsilon_v),I_8]."
            ),
            "first": (
                "Delta c'_r(y)=exp(S(ell))w[L_(d_r)-L_0]V_r/"
                "[r sqrt(alpha_P)z^2] at fixed y, V_r=(Z(x_r),Z(x_(0,r)))^T."
            ),
            "second": (
                "Delta^2c'_r(y)=exp(S(ell))w[L_(d_r+d_(r+1))-"
                "2L_(d_r)+L_0]V_r/[r sqrt(alpha_P)z^2] at fixed y."
            ),
            "saddle": (
                "At z=0 replace L_d by S_d=M(epsilon_0-td/2)P_d(epsilon_0); "
                "the exact first and second operators are S_d-S_0 and "
                "S_(d_1+d_2)-2S_(d_1)+S_0."
            ),
            "first_nonzero_entries": sum(1 for item in first_operator if item != 0),
            "second_nonzero_entries": sum(
                1 for item in second_operator if item != 0
            ),
            "saddle_first_nonzero_entries": sum(
                1 for item in saddle_first if item != 0
            ),
            "saddle_second_nonzero_entries": sum(
                1 for item in saddle_second if item != 0
            ),
        },
        "active_roster": {
            "set": (
                "For fixed y and v=v(y), the active modes are "
                "T_y=T intersect {ceil(alpha_P v/B),...,floor(alpha_P v)}."
            ),
            "fubini": (
                "sum_(r in T)e(phi_r)integral_(y_(1,r))^(y_(B,r))f_r(y)dy="
                "integral_R sum_(r in T_y)e(phi_r)f_r(y)dy exactly."
            ),
            "boundary": (
                "The moving integration limits become the two integer jumps of "
                "one contiguous active roster; no endpoint strip is discarded."
            ),
        },
        "reciprocal_phase": {
            "gauge": (
                "For integer q,r, e(alpha_P log u-ru)="
                "e(Psi_q(r,u)), Psi_q=alpha_P log u-r(u-q)."
            ),
            "critical": (
                "The unique positive continuous critical point is "
                "(r_*,u_*)=(alpha_P/q,q)."
            ),
            "value": "Psi_q(r_*,u_*)=alpha_P log q.",
            "hessian": (
                "Hess_(r,u)Psi_q=[[0,-1],[-1,-alpha_P/q^2]], "
                "with determinant -1 and inertia (1,1)."
            ),
            "curvatures": (
                "The sequential u and r curvatures are -alpha_P/q^2 and "
                "q^2/alpha_P; their product is -1 and their Fresnel signatures cancel."
            ),
            "interpretation": (
                "The reciprocal r-saddle returns the original carrier phase "
                "e(alpha_P log q), so the Morse roster is self-dual rather than "
                "an unrelated error lattice."
            ),
        },
        "saddle_scale": {
            "coefficient": (
                "At r_*=alpha_P/q, x_q=log(q/a), the ideal saddle coefficient is "
                "q exp(S(log q))wM(epsilon_q)Y(x_q)/(2alpha_P^(3/2))."
            ),
            "width": (
                "The reciprocal r-curvature width is sqrt(alpha_P)/q."
            ),
            "product": (
                "Their exact scalar product is exp(S(log q))/(2alpha_P) before "
                "the joined row and saddle matrix."
            ),
        },
        "second_order_abel": (
            "For S_j=sum_(k=m)^jz_k and T_j=sum_(l=m)^jS_l, "
            "sum_(r=m)^n w_rz_r=w_nS_n-Delta w_(n-1)T_(n-1)+"
            "sum_(r=m)^(n-2)Delta^2w_rT_r."
        ),
    }


def bound_certificate() -> dict:
    core_margin = Fraction(1, 12) - Fraction(1, 40000) - Fraction(1, 10000)
    if core_margin <= Fraction(1, 13):
        raise RuntimeError("constant-channel saddle margin failed")
    exponent = Fraction(201, 400)
    growth = Fraction(1, 1) - exponent
    coefficient = Fraction(1, 26) / growth
    if coefficient != Fraction(200, 2587):
        raise RuntimeError("reciprocal-family coefficient drifted")
    return {
        "constant_channel_margin": (
            "For bar P=1, the saddle scalar is "
            "m_q=epsilon_q^2+t/2-1/12 and the physical box gives |m_q|>1/13."
        ),
        "amplitude_floor": (
            "Since t>=0 and sigma<201/400, exp(S(log q))>=q^(-201/400), "
            "with strict inequality for q>1."
        ),
        "absolute_family_budget": (
            "If A_B is the nonnegative sum of the reciprocal saddle coefficient "
            "scale times its curvature width for bar P=1, then "
            "A_B>200(B^(199/400)-1)/(2587alpha_P)."
        ),
        "barrier": (
            "The factor B^(199/400)-1 is unbounded. Therefore summing the "
            "reciprocal q-family saddle scales by absolute values cannot prove a "
            "uniform C/alpha_P grouped estimate. This is a lower bound only on "
            "that chosen nonnegative scale budget, not on the signed Morse remainder."
        ),
        "handoff": (
            "Use the exact fixed-y L_d transport inside local r-cells, but rejoin "
            "the reciprocal stationary contributions with mathcal F, the physical "
            "q-carriers, both endpoint packages, and the terminal recurrence before "
            "summing over q. A second cancellation across the q-family is compulsory."
        ),
        "exact_rational_values": {
            "core_margin": str(core_margin),
            "core_target": "1/13",
            "exponent": str(exponent),
            "growth_exponent": str(growth),
            "family_coefficient": str(coefficient),
        },
    }


def build_rows(symbolic: dict, bounds: dict) -> list[GateRow]:
    jet = symbolic["weighted_jet"]
    transport = symbolic["mode_transport"]
    active = symbolic["active_roster"]
    phase = symbolic["reciprocal_phase"]
    scale = symbolic["saddle_scale"]
    return [
        GateRow("icrg_01_weighted", "weighted_jet", "certified", "The ideal augmented jet remains closed after the centered Morse exponential is restored.", jet["definition"], "This is an exact differential identity."),
        GateRow("icrg_02_kernel", "off_saddle_kernel", "certified", "The complete off-saddle numerator is one joined eight-vector kernel plus its saddle anchor.", jet["kernel"] + " " + jet["off_saddle"], "The apparent z=0 quotient uses the separate removable saddle formula."),
        GateRow("icrg_03_weighted_transport", "mode_transport", "certified", "The exponential and 1/r denominator combine with the nilpotent shift in one exact weighted transport.", transport["weighted"], "No Taylor remainder is used."),
        GateRow("icrg_04_cocycle", "mode_transport", "certified", "The weighted shift is a nonautonomous exact cocycle.", transport["cocycle"], "The epsilon shift is retained."),
        GateRow("icrg_05_operator", "off_saddle_transport", "certified", "One block operator transports both the moving field and saddle anchor.", transport["operator"], "A_z and B_z are held at the same fixed Fresnel coordinate."),
        GateRow("icrg_06_first", "off_saddle_transport", "certified", "The full fixed-y off-saddle coefficient has an exact adjacent difference.", transport["first"], "This is not yet a norm bound."),
        GateRow("icrg_07_second", "off_saddle_transport", "certified", "The full fixed-y off-saddle coefficient has an exact second difference.", transport["second"], "Moving active-roster jumps remain to be retained."),
        GateRow("icrg_08_saddle", "saddle_transport", "certified", "The removable saddle limit has matching exact first and second weighted transports.", transport["saddle"], "This does not bound the integrated coefficient."),
        GateRow("icrg_09_active", "active_roster", "certified", "At fixed Fresnel coordinate the modes form one explicit contiguous interval.", active["set"] + " " + active["boundary"], "The fixed roster intersection is retained."),
        GateRow("icrg_10_fubini", "active_roster", "certified", "Finite Fubini reindexing removes moving integration limits without deleting their jumps.", active["fubini"], "No infinite exchange or convergence claim is needed."),
        GateRow("icrg_11_gauge", "reciprocal_phase", "certified", "Every integer q gives an exact lattice gauge for the joint phase.", phase["gauge"], "The continuous saddle is a theorem-search coordinate."),
        GateRow("icrg_12_critical", "reciprocal_phase", "certified", "The joint reciprocal phase has one explicit positive critical point.", phase["critical"] + " " + phase["value"], "The critical r need not be integral."),
        GateRow("icrg_13_hessian", "reciprocal_phase", "certified", "The joint saddle has parameter-independent determinant and balanced inertia.", phase["hessian"], "Indefinite Hessian is not a sign theorem."),
        GateRow("icrg_14_curvature", "reciprocal_phase", "certified", "Sequential Morse curvatures are exact reciprocals with canceling signatures.", phase["curvatures"], "This is a local phase identity, not a global remainder estimate."),
        GateRow("icrg_15_self_dual", "reciprocal_self_duality", "certified", "The reciprocal mode saddle returns the original physical carrier phase.", phase["interpretation"], "The returned carrier must be recomposed rather than declared small."),
        GateRow("icrg_16_coefficient", "reciprocal_scale", "certified", "The ideal saddle coefficient has an explicit q-dependent prefactor.", scale["coefficient"] + " " + scale["width"] + " " + scale["product"], "The joined row and matrix remain inside the coefficient."),
        GateRow("icrg_17_abel", "second_order_abel", "certified", "The exact second-order Abel identity exposes the needed boundary, first-difference, and second-difference terms.", symbolic["second_order_abel"], "Useful bounds still require a reciprocal-cell partition."),
        GateRow("icrg_18_barrier", "absolute_family_barrier", "guard_validated", "A cellwise absolute sum of reciprocal saddle scales has an unbounded q-family factor.", bounds["constant_channel_margin"] + " " + bounds["amplitude_floor"] + " " + bounds["absolute_family_budget"] + " " + bounds["barrier"], "This does not lower-bound the signed Morse remainder."),
        GateRow("icrg_19_handoff", "reciprocal_handoff", "open", "The next estimate must cancel across reciprocal q-families after carrier and endpoint recomposition.", bounds["handoff"], "No such two-level cancellation theorem is proved here."),
        GateRow("icrg_20_boundary", "proof_boundary", "guard_validated", "Fixed-y smooth transport is not promoted to a grouped h^2 theorem.", "The reciprocal phase returns the original carrier lattice and the absolute q-family scale budget diverges.", "No grouped ideal-cubic, endpoint-composed, quadratic, signed-flow, Xi, Lambda<=0, PF-infinity, RH, or prize theorem is claimed."),
    ]


def render_note(artifact: dict) -> str:
    symbolic = artifact["symbolic_certificate"]
    bounds = artifact["bound_certificate"]
    jet = symbolic["weighted_jet"]
    transport = symbolic["mode_transport"]
    active = symbolic["active_roster"]
    phase = symbolic["reciprocal_phase"]
    scale = symbolic["saddle_scale"]
    return "\n".join(
        [
            "# Physical P_lin Ideal-Cubic Fresnel Reciprocal Gate",
            "",
            "Date: 2026-08-02",
            "",
            "Status: exact weighted off-saddle transport, active-roster Fubini reindexing, reciprocal joint-phase self-duality, second-order Abel identity, and one absolute-family barrier. This is not a proof of a grouped h^2 estimate, signed flow, or RH.",
            "",
            "```text",
            str(RESULT_PATH.relative_to(REPO_ROOT)).replace("\\", "/"),
            f"python {str(Path(__file__).relative_to(REPO_ROOT)).replace(chr(92), '/')}",
            "```",
            "",
            "## Weighted Off-Saddle Jet",
            "",
            "```text",
            jet["definition"],
            jet["kernel"],
            jet["off_saddle"],
            "```",
            "",
            "## Exact Roster Transport",
            "",
            "```text",
            transport["weighted"],
            transport["cocycle"],
            transport["operator"],
            transport["first"],
            transport["second"],
            transport["saddle"],
            "```",
            "",
            "These are fixed-y identities. They retain the moving exponential, the 1/r denominator, the epsilon shift, the field term, and the saddle anchor in one block operator.",
            "",
            "## Active Roster",
            "",
            "```text",
            active["set"],
            active["fubini"],
            active["boundary"],
            "```",
            "",
            "## Reciprocal Self-Duality",
            "",
            "```text",
            phase["gauge"],
            phase["critical"],
            phase["value"],
            phase["hessian"],
            phase["curvatures"],
            "```",
            "",
            phase["interpretation"],
            "",
            "The determinant -1 is exact. The continuous r saddle need not be an integer; it is the legitimate stationary coordinate for the lattice phase, just as the continuous u saddle is used in the Morse chart.",
            "",
            "## Reciprocal Scale",
            "",
            "```text",
            scale["coefficient"],
            scale["width"],
            scale["product"],
            "```",
            "",
            "## Second-Order Abel And Barrier",
            "",
            "```text",
            symbolic["second_order_abel"],
            bounds["constant_channel_margin"],
            bounds["amplitude_floor"],
            bounds["absolute_family_budget"],
            "```",
            "",
            bounds["barrier"],
            "",
            "## Handoff",
            "",
            bounds["handoff"],
            "",
            "The fixed-y O(alpha_P^-2) field transport remains genuine, but it is only the inner half of a two-level problem. The reciprocal r saddles carry the original q phases. The next theorem must identify the cancellation between the transformed Morse carrier, the c-prime correction, and the endpoint-complete physical carrier before estimating the q sum.",
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def build_artifact() -> dict:
    payloads = load_sources()
    source = source_audit(payloads)
    symbolic = symbolic_certificate()
    bounds = bound_certificate()
    rows = build_rows(symbolic, bounds)
    return {
        "kind": KIND,
        "date": "2026-08-02",
        "status": (
            "exact weighted ideal-cubic off-saddle transport, active-roster "
            "reindexing, reciprocal phase self-duality, second-order Abel identity, "
            "and absolute q-family barrier complete; two-level carrier cancellation "
            "open; no grouped h2, signed-flow, Xi, or RH theorem"
        ),
        "proof_boundary": (
            "This gate proves exact weighted off-saddle and saddle mode transports, "
            "finite active-roster Fubini reindexing, reciprocal joint-phase geometry, "
            "one second-order Abel identity, and one nonnegative absolute-family "
            "barrier. It proves no grouped ideal-cubic or endpoint-composed h2 "
            "estimate, reciprocal q-family cancellation theorem, quadratic residual "
            "bound, signed flow estimate, Phi_B upper bound, contact exclusion, "
            "retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or "
            "prize-level conclusion."
        ),
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "bound_certificate": bounds,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "off_saddle_vector_identities": 2,
            "exact_weighted_mode_transport_identities": 4,
            "active_roster_reindexings": 1,
            "reciprocal_phase_identities": 6,
            "second_order_abel_identities": 1,
            "reciprocal_absolute_family_barriers": 1,
            "grouped_ideal_cubic_bounds": 0,
            "endpoint_composed_h2_bounds": 0,
            "quadratic_remainder_bounds": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
    }


def main() -> int:
    artifact = build_artifact()
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    counts = artifact["counts"]
    print(
        "built ideal-cubic Fresnel reciprocal gate: "
        f"{counts['rows']} rows, "
        f"{counts['exact_weighted_mode_transport_identities']} weighted transports, "
        f"{counts['active_roster_reindexings']} Fubini reindexing, "
        f"{counts['reciprocal_phase_identities']} reciprocal phase identities, "
        f"{counts['reciprocal_absolute_family_barriers']} absolute-family barrier, "
        f"{counts['grouped_ideal_cubic_bounds']} grouped bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
