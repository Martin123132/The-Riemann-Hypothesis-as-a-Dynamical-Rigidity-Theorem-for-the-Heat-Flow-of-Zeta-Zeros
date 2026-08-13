#!/usr/bin/env python3
"""Build the quadratic-residual two-carrier kernel closure gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_quadratic_residual_two_carrier_kernel_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "quadratic_primitive": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_current_transfer_quadratic_primitive_gate.json",
    "two_carrier_kernel": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_two_carrier_kernel_reduction.json",
    "physical_coefficients": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_amplitude_gate.json",
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_relative_lift_gate.json",
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


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def complex_symbol(prefix: str) -> sp.Expr:
    real, imag = sp.symbols(f"{prefix}_R {prefix}_I", real=True)
    return real + sp.I * imag


def is_zero(expression: sp.Expr) -> bool:
    return sp.simplify(sp.expand_complex(expression)) == 0


def load_and_audit_sources() -> dict[str, dict[str, str]]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    primitive = payloads["quadratic_primitive"]
    require(
        "E_quad^cur=dQ_quad/dxi exactly" in primitive["exact"]["quadratic_primitive"],
        "quadratic primitive drifted",
    )
    require(
        primitive["summary"]["quadratic_constituents"] == 4,
        "quadratic constituent count drifted",
    )

    carrier = payloads["two_carrier_kernel"]
    require(
        "phase differences" in carrier["exact"]["kernel_interpretation"],
        "Hermitian phase family drifted",
    )
    require(
        "phase sums" in carrier["exact"]["kernel_interpretation"],
        "transpose phase family drifted",
    )
    require(
        "H^+_(n,m)=-(u_n+u_m)^2/8" in carrier["symbolic_certificate"]["leading_q1_kernels"],
        "leading Hermitian kernel drifted",
    )

    coefficients = payloads["physical_coefficients"]
    require(
        "Q=(R+delta)C+D" in coefficients["centered_carrier_certificate"]["rate_defect"],
        "physical Q factorization drifted",
    )
    require(
        "N=R*Q+R_x*C" in coefficients["centered_carrier_certificate"]["rate_defect"],
        "physical N factorization drifted",
    )
    require(
        "R=Q=i*x/2" in coefficients["ideal_certificate"]["specialization"],
        "ideal carrier specialization drifted",
    )

    relative = payloads["relative_lift"]["symbolic_certificate"]
    require(
        "rho_X=Re{nu c_xi mathscr R_N[P_X]}" in relative["outer_remainder"]["observations"],
        "nonterminal observation functional drifted",
    )
    require(
        "mathscr R_N[i(lambda-log N)P]" in relative["relative_lifts"]["hermitian"],
        "Hermitian relative lift drifted",
    )
    require(
        "mathscr R_N[i(lambda-log a-u_N)P]" in relative["relative_lifts"]["transpose"],
        "transpose relative lift drifted",
    )

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_certificate() -> dict[str, str | int]:
    c_l = complex_symbol("c_l")
    c_m = complex_symbol("c_m")
    q_l = complex_symbol("q_l")
    q_m = complex_symbol("q_m")
    r_l = complex_symbol("r_l")
    r_m = complex_symbol("r_m")
    rx_l = complex_symbol("rx_l")
    rx_m = complex_symbol("rx_m")

    a_l = r_l * c_l
    a_m = r_m * c_m
    n_l = r_l * q_l + rx_l * c_l
    n_m = r_m * q_m + rx_m * c_m

    k_t = c_l * n_m + n_l * c_m - a_l * q_m - q_l * a_m
    k_t_factored = (r_m - r_l) * (c_l * q_m - c_m * q_l) + (rx_l + rx_m) * c_l * c_m
    require(is_zero(k_t - k_t_factored), "transpose kernel factorization")

    delta_l = complex_symbol("delta_l")
    delta_m = complex_symbol("delta_m")
    d_l = complex_symbol("d_l")
    d_m = complex_symbol("d_m")
    q_l_physical = (r_l + delta_l) * c_l + d_l
    q_m_physical = (r_m + delta_m) * c_m + d_m
    wedge = c_l * q_m_physical - c_m * q_l_physical
    wedge_factored = (
        c_l * c_m * ((r_m - r_l) + (delta_m - delta_l))
        + c_l * d_m
        - c_m * d_l
    )
    require(is_zero(wedge - wedge_factored), "physical wedge factorization")

    k_h = (
        c_l * sp.conjugate(n_m)
        + n_l * sp.conjugate(c_m)
        - a_l * sp.conjugate(q_m)
        - q_l * sp.conjugate(a_m)
    )
    k_h_factored = (
        (sp.conjugate(r_m) - r_l)
        * (c_l * sp.conjugate(q_m) - q_l * sp.conjugate(c_m))
        + (rx_l + sp.conjugate(rx_m)) * c_l * sp.conjugate(c_m)
    )
    require(is_zero(k_h - k_h_factored), "Hermitian kernel factorization")

    x_l, x_m, u_x = sp.symbols("x_l x_m u_x", real=True)
    ideal_r_l = sp.I * x_l / 2
    ideal_r_m = sp.I * x_m / 2
    ideal_rx = -sp.I * u_x / 2
    ideal_k_h = (
        (sp.conjugate(ideal_r_m) - ideal_r_l) ** 2
        + ideal_rx
        + sp.conjugate(ideal_rx)
    )
    ideal_k_t = (ideal_r_m - ideal_r_l) ** 2 + 2 * ideal_rx
    require(is_zero(ideal_k_h + (x_l + x_m) ** 2 / 4), "ideal Hermitian kernel")
    require(is_zero(ideal_k_t + (x_m - x_l) ** 2 / 4 + sp.I * u_x), "ideal transpose kernel")

    u_v = complex_symbol("u_v")
    u_n = complex_symbol("u_n")
    u_a = complex_symbol("u_a")
    u_q = complex_symbol("u_q")
    real_primitive = 2 * (
        sp.re(u_v) * sp.re(u_n) - sp.re(u_a) * sp.re(u_q)
    )
    h_projection = (
        u_v * sp.conjugate(u_n)
        + u_n * sp.conjugate(u_v)
        - u_a * sp.conjugate(u_q)
        - u_q * sp.conjugate(u_a)
    )
    t_projection = u_v * u_n + u_n * u_v - u_a * u_q - u_q * u_a
    kernel_primitive = sp.re(h_projection + t_projection) / 2
    require(is_zero(real_primitive - kernel_primitive), "double-functional polarization")

    lambda_l, lambda_m, log_a = sp.symbols("lambda_l lambda_m log_a", real=True)
    flow_l = lambda_l - log_a
    flow_m = lambda_m - log_a
    require(
        is_zero(sp.I * flow_l + sp.conjugate(sp.I * flow_m) - sp.I * (lambda_l - lambda_m)),
        "Hermitian flow multiplier",
    )
    require(
        is_zero(sp.I * flow_l + sp.I * flow_m - sp.I * (lambda_l + lambda_m - 2 * log_a)),
        "transpose flow multiplier",
    )

    return {
        "physical_rows": "P_V=C, P_mathcalN=mathcal N=RQ+R_xC, P_A=A=RC, and P_Q=Q=(R+delta)C+D.",
        "transpose_kernel": "K_T(lambda,mu)=C_lambda mathcalN_mu+mathcalN_lambda C_mu-A_lambda Q_mu-Q_lambda A_mu=(R_mu-R_lambda)(C_lambda Q_mu-C_mu Q_lambda)+(R_(lambda,x)+R_(mu,x))C_lambda C_mu.",
        "physical_wedge": "C_lambda Q_mu-C_mu Q_lambda=C_lambda C_mu[(R_mu-R_lambda)+(delta_mu-delta_lambda)]+C_lambda D_mu-C_muD_lambda.",
        "hermitian_kernel": "K_H(lambda,mu)=C_lambda conjugate(mathcalN_mu)+mathcalN_lambda conjugate(C_mu)-A_lambda conjugate(Q_mu)-Q_lambda conjugate(A_mu)=(conjugate(R_mu)-R_lambda)[C_lambda conjugate(Q_mu)-Q_lambda conjugate(C_mu)]+[R_(lambda,x)+conjugate(R_(mu,x))]C_lambda conjugate(C_mu).",
        "ideal_kernels": "For C=1, D=delta=0, R=i*x/2, R_x=-i*u_x/2, one has K_H^0=-(x_lambda+x_mu)^2/4 and K_T^0=-(x_mu-x_lambda)^2/4-i*u_x.",
        "flow_multipliers": "The fixed-coefficient auxiliary lift gives i(lambda-mu) on K_H and i(lambda+mu-2log a) on K_T.",
        "symbolic_audits": 8,
    }


def finite_functional_certificate() -> dict[str, str | int]:
    i = sp.I
    xs = (sp.Rational(-2), sp.Rational(1, 3), sp.Rational(5, 2))
    weights = (1 + i / 2, -sp.Rational(2, 3) + i / 5, sp.Rational(3, 7) - 2 * i / 9)
    c = (1 + i / 7, sp.Rational(2, 3) - i / 4, -sp.Rational(1, 5) + 2 * i / 3)
    q = (sp.Rational(4, 5) - i / 6, -sp.Rational(3, 8) + i / 2, sp.Rational(7, 9) + i / 10)
    r = (-sp.Rational(1, 3) + i / 2, sp.Rational(5, 7) - i / 9, -sp.Rational(2, 5) - 3 * i / 8)
    rx = (i / 11, -sp.Rational(1, 13) - i / 17, sp.Rational(2, 19) - i / 23)
    a = tuple(r_j * c_j for r_j, c_j in zip(r, c))
    n = tuple(r_j * q_j + rx_j * c_j for r_j, q_j, rx_j, c_j in zip(r, q, rx, c))

    def functional(values: tuple[sp.Expr, ...]) -> sp.Expr:
        return sum(weight * value for weight, value in zip(weights, values))

    def lifted(values: tuple[sp.Expr, ...]) -> sp.Expr:
        return sum(i * x * weight * value for x, weight, value in zip(xs, weights, values))

    u_c, u_n, u_a, u_q = (functional(values) for values in (c, n, a, q))
    du_c, du_n, du_a, du_q = (lifted(values) for values in (c, n, a, q))
    primitive_direct = 2 * (sp.re(u_c) * sp.re(u_n) - sp.re(u_a) * sp.re(u_q))
    current_direct = 2 * (
        sp.re(du_c) * sp.re(u_n)
        + sp.re(u_c) * sp.re(du_n)
        - sp.re(du_a) * sp.re(u_q)
        - sp.re(u_a) * sp.re(du_q)
    )

    primitive_kernel = 0
    current_kernel = 0
    hermitian_diagonal_current = 0
    transpose_diagonal_current = 0
    for left in range(len(xs)):
        for right in range(len(xs)):
            k_h = (
                c[left] * sp.conjugate(n[right])
                + n[left] * sp.conjugate(c[right])
                - a[left] * sp.conjugate(q[right])
                - q[left] * sp.conjugate(a[right])
            )
            k_t = (
                c[left] * n[right]
                + n[left] * c[right]
                - a[left] * q[right]
                - q[left] * a[right]
            )
            h_weight = weights[left] * sp.conjugate(weights[right])
            t_weight = weights[left] * weights[right]
            primitive_kernel += (k_h * h_weight + k_t * t_weight) / 2
            current_kernel += (
                i * (xs[left] - xs[right]) * k_h * h_weight
                + i * (xs[left] + xs[right]) * k_t * t_weight
            ) / 2
            if left == right:
                hermitian_diagonal_current += i * (xs[left] - xs[right]) * k_h * h_weight / 2
                transpose_diagonal_current += i * (xs[left] + xs[right]) * k_t * t_weight / 2

    primitive_difference = sp.simplify(primitive_direct - sp.re(primitive_kernel))
    current_difference = sp.simplify(current_direct - sp.re(current_kernel))
    hermitian_diagonal_current = sp.simplify(sp.re(hermitian_diagonal_current))
    transpose_diagonal_current = sp.simplify(sp.re(transpose_diagonal_current))
    require(primitive_difference == 0, "finite primitive functional")
    require(current_difference == 0, "finite current functional")
    require(hermitian_diagonal_current == 0, "Hermitian diagonal current")
    require(transpose_diagonal_current != 0, "transpose diagonal guard")

    return {
        "carrier_count": len(xs),
        "primitive_difference": str(primitive_difference),
        "current_difference": str(current_difference),
        "hermitian_diagonal_current": str(hermitian_diagonal_current),
        "transpose_diagonal_current": str(transpose_diagonal_current),
        "finite_functional_audits": 4,
    }


def exact_certificate() -> dict[str, str]:
    return {
        "domain": "Work on one fixed reciprocal-roster cell and the fixed physical coefficient chart. The terminal conditional mathcal T_N block has already been removed, and mathscr R_N is retained as one linear nonterminal functional.",
        "functional": "Put x_lambda=lambda-log a and mathcal L_xi[P]=(nu c_xi/|nu|)mathscr R_N[P]. For X in {V,mathcal N,A,Q}, U_X=mathcal L_xi[P_X], e_X=|nu|Re(U_X), and partial_xi mathcal L_xi[P]=mathcal L_xi[i x_lambda P].",
        "primitive_kernel": "Q_quad=(1/2)Re{(mathcal L_xi tensor conjugate(mathcal L_xi))[K_H]+(mathcal L_xi tensor mathcal L_xi)[K_T]}.",
        "current_kernel": "E_quad^cur=(1/2)Re{(mathcal L_xi tensor conjugate(mathcal L_xi))[i(lambda-mu)K_H]+(mathcal L_xi tensor mathcal L_xi)[i(lambda+mu-2log a)K_T]}.",
        "normalization": "The factor 1/2 is forced by 2Re(z)Re(w)=Re(z conjugate(w)+zw). It is not an adjustable convention.",
        "old_kernel_match": "Because x_lambda=-u_lambda on a retained physical atom, K_H^0/2=-(u_lambda+u_mu)^2/8 and K_T^0/2=-(u_lambda-u_mu)^2/8-i*u_x/2, exactly the kernels already certified in Section 11.160.",
        "diagonal": "The Hermitian current multiplier i(lambda-mu) vanishes identically on lambda=mu. The transpose multiplier becomes 2i(lambda-log a), so its diagonal remains unless a separate source identity kills it.",
        "phase_closure": "Hermitian products use one weight and one conjugate weight and therefore carry only phase differences. Transpose products use two unconjugated weights and carry only phase sums. The residual-residual determinant creates no third phase family.",
        "source_grouping": "Keep mathscr R_N[P]=mathcal U_N[P]-mathcal T_1[P]-mathcal U_1[P]+mathcal I_N[P] composed until K_H and K_T have been formed. Expanding its tensor square first creates sixteen coordinate terms and destroys the certified tie-invariant cancellations; those terms are not new ownership packages.",
        "signed_join": "The pointwise route must insert the displayed E_quad^cur kernel into E_main^cur=E_FAN+E_quad^cur before any Hermitian/transpose arithmetic decomposition. A separate absolute |E_quad| allowance is neither supplied nor required by this closure.",
        "live_obligation": "Prove or falsify, on every adverse actual-source witness cell, the joint signed inequality E_FAN+E_quad^cur<-(1130179/858000)rho A_T, using the two existing phase families and preserving the grouped mathscr R_N tensor square; alternatively return to the exact endpoint determinant transfer.",
        "pi_provenance": "No new pi is introduced. The u_x=h^2/(8pi) term in the ideal transpose kernel is inherited from the completed-zeta/Riemann--Siegel saddle normalization already audited in Section 11.160.",
        "proof_boundary": "This proves exact kernel closure, factorization, normalization, phase multipliers, diagonal classification, and ideal-kernel agreement. It proves no signed Type-I/II, Vaughan, reciprocal-pair, or endpoint determinant estimate; no complete-current sign, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }


def build_rows(exact: dict[str, str], symbolic: dict[str, str | int]) -> list[GateRow]:
    return [
        GateRow("qrtk_01_domain", "domain", "proved", "The quadratic residual lives on the certified fixed-cell nonterminal functional.", exact["domain"], "No moving-roster derivative is hidden."),
        GateRow("qrtk_02_functional", "normalized functional", "proved", "All four residual observations and lifts use one normalized complex functional.", exact["functional"], "The source phase remains inside mathcal L_xi."),
        GateRow("qrtk_03_rows", "physical rows", "proved", "The four coefficient rows retain their exact physical factorization.", str(symbolic["physical_rows"]), "No small correction is set to zero here."),
        GateRow("qrtk_04_transpose", "transpose kernel", "proved", "The residual-residual transpose determinant factors into one coefficient wedge and one rate derivative.", str(symbolic["transpose_kernel"]), "This is a complex identity before projection."),
        GateRow("qrtk_05_wedge", "physical wedge", "proved", "The transpose wedge exposes all rate, defect, and correction differences.", str(symbolic["physical_wedge"]), "D and delta remain explicit."),
        GateRow("qrtk_06_hermitian", "Hermitian kernel", "proved", "The residual-residual Hermitian determinant has the conjugate companion factorization.", str(symbolic["hermitian_kernel"]), "No imaginary coefficient is discarded."),
        GateRow("qrtk_07_primitive", "double-functional primitive", "proved", "The determinant primitive is exactly one Hermitian plus one transpose kernel functional.", exact["primitive_kernel"] + " " + exact["normalization"], "The factor one half is mandatory."),
        GateRow("qrtk_08_current", "frequency current", "proved", "Differentiation supplies the exact phase-difference and phase-sum multipliers.", exact["current_kernel"], "Physical coefficient rows are fixed on this auxiliary segment."),
        GateRow("qrtk_09_diagonal", "diagonal classification", "proved", "Only the Hermitian current diagonal vanishes structurally.", exact["diagonal"], "The transpose diagonal may not be deleted."),
        GateRow("qrtk_10_ideal", "ideal specialization", "proved", "The correction-free residual kernels are closed squares plus the inherited u_x term.", str(symbolic["ideal_kernels"]), "This specialization supplies no aggregate sign."),
        GateRow("qrtk_11_match", "legacy-kernel match", "proved", "The half-normalized kernels equal the already certified physical two-carrier kernels.", exact["old_kernel_match"], "This is an exact notation match, not an analogy."),
        GateRow("qrtk_12_phases", "phase-family closure", "proved", "The quadratic residual introduces no new oscillatory phase family.", exact["phase_closure"], "Both existing families remain necessary."),
        GateRow("qrtk_13_grouping", "source grouping", "guard_validated", "The nonterminal tensor square must be formed before endpoint/interior expansion.", exact["source_grouping"], "Sixteen expanded coordinates are not sixteen independent errors."),
        GateRow("qrtk_14_join", "signed-main join", "guard_validated", "The exact quadratic current belongs inside the signed main theorem.", exact["signed_join"], "No detached absolute budget is introduced."),
        GateRow("qrtk_15_target", "live signed theorem", "open", "The actual-source joint signed estimate remains the next arithmetic theorem.", exact["live_obligation"], "Neither the pointwise nor endpoint-transfer alternative is proved here."),
        GateRow("qrtk_16_boundary", "proof boundary", "guard_validated", "Kernel closure is not a current-sign or RH proof.", exact["proof_boundary"], "All prize-level claims remain open."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic"]
    finite = payload["finite_functional"]
    lines = [
        "# Newman C1 Quadratic-Residual Two-Carrier Kernel Gate",
        "",
        "Date: 2026-08-05",
        "",
        "Status: exact residual-residual kernel closure; this is not a proof of the joint signed arithmetic estimate or RH.",
        "",
        "## Functional And Kernels",
        "",
        exact["domain"],
        "",
        exact["functional"],
        "",
        "```text",
        str(symbolic["transpose_kernel"]),
        "",
        str(symbolic["physical_wedge"]),
        "",
        str(symbolic["hermitian_kernel"]),
        "```",
        "",
        "## Primitive And Current",
        "",
        "```text",
        exact["primitive_kernel"],
        exact["current_kernel"],
        "```",
        "",
        exact["normalization"],
        "",
        exact["diagonal"],
        "",
        "## Exact Kernel Match",
        "",
        "```text",
        str(symbolic["ideal_kernels"]),
        "```",
        "",
        exact["old_kernel_match"],
        "",
        exact["phase_closure"],
        "",
        f"The independent finite audit used {finite['carrier_count']} exact rational complex carriers; both functional differences and the Hermitian current diagonal were `{finite['primitive_difference']}`, `{finite['current_difference']}`, and `{finite['hermitian_diagonal_current']}`, while the transpose diagonal was nonzero.",
        "",
        "## Arithmetic Handoff",
        "",
        exact["source_grouping"],
        "",
        exact["signed_join"],
        "",
        exact["live_obligation"],
        "",
        exact["pi_provenance"],
        "",
        exact["proof_boundary"],
        "",
        payload["success"],
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    source_audit = load_and_audit_sources()
    symbolic = symbolic_certificate()
    finite = finite_functional_certificate()
    exact = exact_certificate()
    rows = build_rows(exact, symbolic)
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-05",
        "status": "quadratic nonterminal determinant closed into the existing Hermitian and transpose two-carrier kernels; joint signed estimate open",
        "source_sha256": {key: value["sha256"] for key, value in source_audit.items()},
        "source_audit": source_audit,
        "symbolic": symbolic,
        "finite_functional": finite,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "source_audits": len(source_audit),
            "symbolic_audits": symbolic["symbolic_audits"],
            "finite_functional_audits": finite["finite_functional_audits"],
            "kernel_families": 2,
            "new_phase_families": 0,
            "hermitian_diagonal_current": 0,
            "transpose_diagonal_guards": 1,
            "open_rows": sum(row.readiness == "open" for row in rows),
        },
        "next_action": exact["live_obligation"],
        "pi_provenance": exact["pi_provenance"],
        "proof_boundary": exact["proof_boundary"],
        "success": "built Newman C1 quadratic residual two-carrier kernel gate: 16 rows, 0 issues, 8 symbolic audits, 4 finite-functional audits, 4 source audits, 2 exact kernel families, 1 open signed theorem row",
    }
    atomic_write(RESULT_PATH, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(payload["success"])


if __name__ == "__main__":
    main()
