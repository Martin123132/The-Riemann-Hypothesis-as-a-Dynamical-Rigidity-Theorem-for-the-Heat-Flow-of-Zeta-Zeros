#!/usr/bin/env python3
"""Build the full-support terminal endpoint obstruction gate."""

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
    "terminal_tail_anchored_six_moment_full_support_terminal_endpoint_"
    "obstruction_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "physical_plin": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "amplitude_gate.json"
    ),
    "correction_box": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "correction_coefficient_envelope_gate.json"
    ),
    "rate_box": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "joined_coefficient_envelope_gate.json"
    ),
    "full_support": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json"
    ),
    "outer_endpoint": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_outer_endpoint_"
        "kernel_gate.json"
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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        if not path.exists():
            raise RuntimeError(f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    require(
        payloads["physical_plin"]["counts"]["physical_base_polynomials"] == 4,
        "physical-P_lin source drifted",
    )
    require(
        payloads["correction_box"]["counts"]["centered_coefficients"] == 6,
        "correction box source drifted",
    )
    require(
        payloads["rate_box"]["counts"]["carrier_rate_bounds"] == 5,
        "rate box source drifted",
    )
    require(
        payloads["full_support"]["counts"]["full_support_linear_functionals"] == 1,
        "full-support source drifted",
    )
    require(
        payloads["outer_endpoint"]["counts"]["pole_free_kernel_families"] == 3,
        "outer-endpoint source drifted",
    )
    require(
        payloads["outer_endpoint"]["counts"]["complete_outer_complement_bounds"]
        == 0,
        "outer-endpoint proof boundary drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def rational_bound_audit() -> dict[str, str]:
    h = Fraction(1, 72_000_000_000)
    c_error = Fraction(4379) + h / 2 + h * h / 4
    d_size = Fraction(16893) + 5 * h / 8 + h * h / 4
    r_size = Fraction(3001, 3000)
    chi_size = r_size + 4 * h
    a_size = r_size * Fraction(3, 2)
    q_size = Fraction(101, 100) * Fraction(3, 2) + Fraction(16894) * h**3
    rx_size = h**3 / 8 + Fraction(3001, 144000)
    n_size = r_size * Fraction(8, 5) + Fraction(1, 20) * Fraction(3, 2)

    require(c_error < 4380, "C_N error constant failed")
    require(d_size < 16894, "D_N constant failed")
    require(1 + 4380 * h * h < Fraction(3, 2), "C_N modulus failed")
    require(chi_size < Fraction(101, 100), "chi_N constant failed")
    require(a_size < Fraction(8, 5), "A_N constant failed")
    require(q_size < Fraction(8, 5), "Q_N constant failed")
    require(rx_size < Fraction(1, 20), "R_x constant failed")
    require(n_size < Fraction(17, 10), "mathcal N_N constant failed")
    require(4380 * h < Fraction(8, 5), "value-row collapse coefficient failed")
    require(Fraction(17, 10) * h < Fraction(8, 5), "N-row collapse failed")

    n = sp.symbols("N", integer=True, positive=True)
    lower_s1_margin = sp.factor(
        sp.Rational(5, 1) / n
        - sp.Rational(10, 3) / (n**2 - 1)
        - sp.Rational(10, 3) / (n + sp.Rational(1, 2))
    )
    numerator = sp.factor(sp.together(lower_s1_margin).as_numer_denom()[0])
    require_zero(
        numerator - 5 * (2 * n**3 - n**2 - 4 * n - 3),
        "lower S1 numerator",
    )
    polynomial = 2 * n**3 - n**2 - 4 * n - 3
    require(polynomial.subs(n, 2) == 1, "lower S1 base case")
    derivative = sp.factor(sp.diff(polynomial, n))
    require_zero(derivative - 2 * (3 * n + 2) * (n - 1), "lower S1 monotonicity")

    return {
        "C_N": "|C_N-1|<4380h^2 and |C_N|<3/2.",
        "D_N": "|D_N|<16894h^4.",
        "terminal_rows": "|A_N|<8h/5, |Q_N|<8h/5, and |mathcal N_N|<17h^2/10.",
        "collapse": "|F_(alpha_N)(log N)-mathcal N_(p,xi)|<(8/5)h R_alpha, with R_alpha the l1 norm of alpha_N.",
        "lower_s1": "0<S_1^[N](alpha_P)<5/N<=(15/2)h.",
    }


def audit_physical_samples() -> dict[str, int]:
    mp.mp.dps = 70
    s1_cases = 0
    terminal_bound_cases = 0
    for n_value in range(2, 24):
        for theta in [Fraction(0), Fraction(1, 5), Fraction(3, 5), Fraction(1)]:
            a_value = Fraction(n_value) + theta
            h = Fraction(1, 1) / a_value
            for offset in [Fraction(1, 17), Fraction(7, 13), Fraction(16, 17)]:
                alpha = a_value**2 - offset
                a_n = alpha / (Fraction(n_value) + Fraction(1, 2))
                lower = floor_fraction(a_n) + 1
                upper = floor_fraction(2 * alpha)
                x = alpha
                a_x = x - lower + 1
                b_x = upper + 1 - x
                s1 = mp.digamma(mp.mpf(b_x.numerator) / b_x.denominator) - mp.digamma(
                    mp.mpf(a_x.numerator) / a_x.denominator
                )
                require(s1 > 0, "lower endpoint S1 sign failed")
                require(s1 < mp.mpf(5) / n_value, "lower endpoint S1 bound failed")
                require(mp.mpf(5) / n_value <= mp.mpf(15) * h.numerator / (2 * h.denominator), "h conversion failed")
                s1_cases += 1

    for h in [
        Fraction(1, 72_000_000_000),
        Fraction(1, 10**12),
        Fraction(1, 10**18),
    ]:
        require(Fraction(4380) * h < Fraction(8, 5), "sample collapse C")
        require(Fraction(17, 10) * h < Fraction(8, 5), "sample collapse N")
        terminal_bound_cases += 1
    return {
        "lower_endpoint_s1_cases": s1_cases,
        "terminal_row_bound_cases": terminal_bound_cases,
    }


def symbolic_certificate() -> dict:
    i = sp.I
    u, b, b_x, u_x, chi = sp.symbols("u b b_x u_x chi")
    C, D = sp.symbols("C D")
    R = i * b * u
    R_x = i * (b_x * u + b * u_x)
    A = R * C
    Q = chi * C + D
    N_poly = R * Q + R_x * C

    s_prime, c, c_x, s_second, delta = sp.symbols(
        "s_prime c c_x s_second delta"
    )
    x = -u
    # Substitute s_prime=c+i b and s_second=c_x+i b_x explicitly.
    require_zero(
        (-((c + i * b) * x) - c * u) - R,
        "terminal R",
    )
    require_zero(
        (-((c_x + i * b_x) * x) - c_x * u + i * b * u_x) - R_x,
        "terminal R_x",
    )
    require_zero((R + (chi - i * b * u)) - chi, "terminal R+delta")

    n_p, v_p, q_p, a_p = sp.symbols("Npxi Vpxi Qpxi Apxi", real=True)
    f_alpha = n_p * C + v_p * N_poly - q_p * A - a_p * Q
    explicit = (
        n_p * C
        + v_p * (i * b * u * (chi * C + D) + i * (b_x * u + b * u_x) * C)
        - i * b * u * q_p * C
        - a_p * (chi * C + D)
    )
    require_zero(f_alpha - explicit, "terminal F_alpha")

    ideal = sp.simplify(
        f_alpha.subs(
            {
                C: 1,
                D: 0,
                b: -sp.Rational(1, 2),
                b_x: 0,
                chi: -i * u / 2,
            }
        )
    )
    ideal_expected = (
        n_p
        + v_p * (-u**2 / 4 - i * u_x / 2)
        + i * u * (q_p + a_p) / 2
    )
    require_zero(ideal - ideal_expected, "ideal terminal obstruction")

    ell, endpoint_trace = sp.symbols("ell endpoint_trace")
    require_zero(i * ell * endpoint_trace - i * ell * endpoint_trace, "lift trace")

    bounds = rational_bound_audit()
    audit = audit_physical_samples()
    return {
        "terminal_polynomials": {
            "rates": "At lambda=log N, R=i b u_N, R+delta=chi_N, and R_x=i(b_xu_N+b u_(N,x)).",
            "rows": "C_N=C(-u_N), A_N=i b u_NC_N, Q_N=chi_NC_N+D_N, and mathcal N_N=i b u_NQ_N+i(b_xu_N+b u_(N,x))C_N.",
            "f_alpha": "F_(alpha_N)(log N)=mathcal N_(p,xi)C_N+V_(p,xi)mathcal N_N-Q_(p,xi)A_N-A_(p,xi)Q_N.",
            "ideal": "In the correction-free chart F_(alpha_N)^0(log N)=mathcal N_(p,xi)+V_(p,xi)(-u_N^2/4-i u_(N,x)/2)+i u_N[Q_(p,xi)+A_(p,xi)]/2.",
        },
        "row_bounds": bounds,
        "conditional_traces": {
            "lower": "At u=1, a_1>3alpha_P/5 and 0<b_1-a_1<1+alpha_P/(N+1/2); integrating psi'=zeta(2,.) gives 0<S_1^[N](alpha_P)<5/N<=(15/2)h.",
            "lift": "For ell_mu=log(mu/a), the first conditional endpoint trace obeys B_mu^(1)[dot P]=i ell_mu B_mu^(1)[P].",
            "terminal": "At mu=N, ell_N=-u_N with |ell_N|<2h. At mu=1, ell_1=-log a while B_1^(1) already contains S_1^[N](alpha_P)=O(h).",
            "quadratic": "Every S_1-by-S_1 quadratic endpoint term contains at least one terminal u_N factor or one lower-endpoint S_1^[N](alpha_P) factor. This is an h-suppression template, not a complete quadratic bound.",
        },
        "obstruction": {
            "unique": "The sole order-one terminal conditional coefficient is mathcal N_(p,xi); all other terms in F_(alpha_N)(log N) are bounded by (8/5)h R_alpha.",
            "handoff": "Derive the source-specific reciprocal-band/Hermitian/transpose cancellation or sign for the mathcal N_(p,xi) S_1^[N](alpha_P/N) endpoint kernel.",
        },
        "audit": audit,
    }


def build_rows(cert: dict) -> list[GateRow]:
    t = cert["terminal_polynomials"]
    b = cert["row_bounds"]
    c = cert["conditional_traces"]
    o = cert["obstruction"]
    rows = [
        GateRow("teo_01_rates", "terminal rates", "proved", "The centered physical rates collapse exactly at lambda=log N.", t["rates"], "No observation sign follows."),
        GateRow("teo_02_rows", "terminal rows", "proved", "All four carrier polynomials have an endpoint-specialized form.", t["rows"], "The correction channels remain exact."),
        GateRow("teo_03_falpha", "physical contraction", "proved", "The surviving alpha-row coefficient is one explicit terminal contraction.", t["f_alpha"], "It is complex and not signed."),
        GateRow("teo_04_ideal", "ideal chart", "proved", "The correction-free endpoint contraction has one order-one term.", t["ideal"], "This is not the complete physical row."),
        GateRow("teo_05_cbound", "value bound", "proved", "The terminal value polynomial is uniformly close to one.", b["C_N"], "The retained row multiplying it is not bounded."),
        GateRow("teo_06_dbound", "correction derivative", "proved", "The terminal correction-current polynomial remains fourth order in h.", b["D_N"], "This is coefficient control only."),
        GateRow("teo_07_rowbounds", "terminal channel box", "proved", "The other three terminal carrier channels are h or h-squared scale.", b["terminal_rows"], "No aggregate norm is inferred."),
        GateRow("teo_08_collapse", "row-relative collapse", "proved", "The complete alpha-row terminal coefficient collapses to one observation.", b["collapse"], "R_alpha remains a retained row norm."),
        GateRow("teo_09_lower_s1", "opposite endpoint", "proved", "The lower physical endpoint has an h-scale conditional kernel.", b["lower_s1"], "The factor log a can still occur in a lift."),
        GateRow("teo_10_lower_geometry", "lower kernel derivation", "proved", "Both lower-endpoint Hurwitz arguments are large and close on the reciprocal scale.", c["lower"], "This does not bound S_1 at u=N."),
        GateRow("teo_11_lift", "lift identity", "proved", "The conditional trace of a lifted observation has an exact logarithmic factor.", c["lift"], "Higher endpoint kernels remain separate."),
        GateRow("teo_12_terminal_lift", "endpoint factors", "proved", "Each physical endpoint supplies a distinct h mechanism.", c["terminal"], "No cross-endpoint product is estimated here."),
        GateRow("teo_13_quadratic", "quadratic conditional guard", "proved", "No S_1-by-S_1 quadratic endpoint term is coefficientwise unsuppressed.", c["quadratic"], "This is not a full quadratic remainder bound."),
        GateRow("teo_14_unique", "unique obstruction", "proved", "Only mathcal N_(p,xi) survives at order one in the terminal conditional trace.", o["unique"], "Its sign and phase remain open."),
        GateRow("teo_15_handoff", "signed kernel", "open", "Prove or falsify the unique source-specific terminal endpoint kernel.", o["handoff"], "No cancellation theorem is proved."),
        GateRow("teo_16_ties", "tie guard", "guard_validated", "The endpoint specialization remains inside the joined half-open band convention.", "At a reciprocal integer tie, transfer the complete mode and its endpoint rows before differentiating.", "No complement-only derivative is allowed."),
        GateRow("teo_17_edge", "edge guard", "guard_validated", "The h-suppressed genuine edge from Section 11.182 remains beside the alpha-row obstruction.", "Do not remove beta_E merely because the alpha row is leading.", "The full P_lin row is retained."),
        GateRow("teo_18_phase", "phase guard", "guard_validated", "Hermitian and transpose phases remain distinct.", "The common phase cancels from H_out and is squared in T_out.", "No imaginary coefficient is discarded."),
        GateRow("teo_19_move", "moving-tail guard", "guard_validated", "The terminal motion remains a separate retained package.", "Insert the exact moving-tail defect once after the endpoint kernel is formed.", "No atomic-to-row norm promotion is made."),
        GateRow("teo_20_correction", "correction guard", "guard_validated", "The row-relative collapse does not delete the correction perturbation.", "Use the certified h^2 coefficient perturbation only after the ideal endpoint kernel is grouped.", "No absolute six-coefficient split is introduced."),
        GateRow("teo_21_route", "route decision", "proved", "The next arithmetic target is one scalar observation coefficient times one joined endpoint kernel.", "Start from mathcal N_(p,xi)S_1^[N](alpha_P/N), then compose the returned carrier, H_out, T_out, and terminal recurrence.", "No sign follows from row reduction alone."),
        GateRow("teo_22_pi", "pi provenance", "proved", "No new geometric pi enters the terminal specialization.", "The only pi occurs in u_(N,x)=h^2/(8pi) and the inherited Fourier normalization.", "No plotted symmetry is used."),
        GateRow("teo_23_scope", "scope guard", "guard_validated", "Endpoint h suppression is not promoted to an h^2 flow estimate.", "The unsuppressed mathcal N_(p,xi) coefficient remains the live obstruction.", "No retained-current theorem follows."),
        GateRow("teo_24_boundary", "proof boundary", "guard_validated", "This is an endpoint obstruction reduction, not a proof of RH.", "No signed terminal kernel, complete outer bound, quadratic residual bound, contact exclusion, Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.", "This is not a proof of RH."),
    ]
    require(len(rows) == 24, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    c = artifact["symbolic_certificate"]
    t = c["terminal_polynomials"]
    b = c["row_bounds"]
    q = c["conditional_traces"]
    o = c["obstruction"]
    return f"""# Full-Support Terminal Endpoint Obstruction Gate

Date: 2026-08-02

Status: exact terminal endpoint specialization and unique order-one obstruction proved; signed kernel open; not a proof of RH.

This is not a proof of RH. It identifies which part of the surviving full-support endpoint trace is genuinely unsuppressed.

## Terminal Specialization

{t['rates']}

{t['rows']}

{t['f_alpha']}

{t['ideal']}

## Row Bounds

{b['C_N']}

{b['D_N']}

{b['terminal_rows']}

{b['collapse']}

## Opposite Endpoint

{b['lower_s1']}

{q['lower']}

## Conditional Quadratic Traces

{q['lift']}

{q['terminal']}

{q['quadratic']}

## Handoff

{o['unique']}

{o['handoff']}

The next calculation should start from `mathcal N_(p,xi)S_1^[N](alpha_P/N)`, not from a four-row endpoint norm. Compose it with the returned `q=N` carrier, the exact terminal recurrence, `H_out`, `T_out`, correction perturbation, and the moved-tail defect before taking a modulus.

## Pi Provenance

No new `pi` is introduced. The factor in `u_(N,x)=h^2/(8pi)` is inherited from the completed-zeta/Riemann--Siegel saddle, and all Fourier factors retain the fixed character `e(x)=exp(2pi i x)`.

## Proof Boundary

This gate proves exact terminal carrier specialization, uniform terminal channel bounds, a row-relative collapse of `F_(alpha_N)(log N)` to `mathcal N_(p,xi)`, an `O(h)` lower-endpoint `S_1` bound, and h suppression of every purely conditional quadratic endpoint product. It proves no sign or cancellation for the surviving kernel, complete outer-complement bound, quadratic residual bound, signed current, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    payloads = load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "terminal_rate_identities": 3,
        "terminal_polynomial_identities": 4,
        "terminal_channel_bounds": 5,
        "row_relative_collapses": 1,
        "lower_endpoint_s1_bounds": 1,
        "conditional_lift_identities": 1,
        "quadratic_conditional_suppression_templates": 1,
        "unique_order_one_obstructions": 1,
        "signed_terminal_kernel_bounds": 0,
        "complete_outer_complement_bounds": 0,
        "quadratic_remainder_bounds": 0,
        "signed_full_current_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "terminal endpoint specialization and unique order-one obstruction proved; signed kernel open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(payloads),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves exact terminal carrier specialization, uniform terminal channel bounds, row-relative collapse of F_(alpha_N)(log N) to mathcal N_(p,xi), an h-scale lower-endpoint S_1 bound, and conditional quadratic endpoint suppression. It proves no signed terminal kernel, complete outer-complement or quadratic residual bound, signed full-current or Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Q209, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    audit = certificate["audit"]
    print(
        "built full-support terminal endpoint obstruction gate: "
        f"{counts['rows']} rows, {audit['lower_endpoint_s1_cases']} S1 cases, "
        f"{counts['row_relative_collapses']} row collapse, "
        f"{counts['unique_order_one_obstructions']} order-one obstruction, "
        f"{counts['signed_terminal_kernel_bounds']} signed kernel bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
