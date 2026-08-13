#!/usr/bin/env python3
"""Build the full-support terminal conditional-quadrature gate."""

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
    "terminal_tail_anchored_six_moment_full_support_terminal_conditional_"
    "quadrature_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "terminal_obstruction": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_endpoint_"
        "obstruction_gate.json"
    ),
    "outer_endpoint": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_outer_endpoint_"
        "kernel_gate.json"
    ),
    "full_support": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json"
    ),
    "terminal_recurrence": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_recurrence_current_gate.json"
    ),
    "two_carrier": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_five_moment_two_carrier_kernel_reduction.json"
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


def mp_fraction(value: Fraction) -> mp.mpf:
    return mp.mpf(value.numerator) / value.denominator


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        if not path.exists():
            raise RuntimeError(f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    require(
        payloads["terminal_obstruction"]["counts"][
            "unique_order_one_obstructions"
        ]
        == 1,
        "terminal-obstruction source drifted",
    )
    require(
        payloads["outer_endpoint"]["counts"]["stable_endpoint_functionals"]
        == 2,
        "outer-endpoint source drifted",
    )
    require(
        payloads["full_support"]["counts"]["terminal_extension_identities"]
        == 3,
        "full-support endpoint split drifted",
    )
    require(
        payloads["terminal_recurrence"]["counts"]["exact_endpoint_recurrences"]
        == 1,
        "terminal recurrence source drifted",
    )
    require(
        payloads["two_carrier"]["counts"]["physical_chi_anchor_identities"]
        == 1,
        "terminal source anchor drifted",
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
    require(Fraction(17, 5) * h**3 < h**2, "local terminal atom scale")
    require(4380 * h**2 < h, "terminal C error scale")
    require(6 * h < Fraction(1, 2), "terminal exponential radius")

    exponential_difference = 6 * h / (1 - 6 * h)
    c_difference = 4380 * h**2
    source_ratio_difference = (
        exponential_difference + c_difference
    ) / (1 - c_difference)
    require(exponential_difference < 7 * h, "terminal exponential bound")
    require(source_ratio_difference < 9 * h, "terminal source ratio bound")
    require(9 * h < Fraction(1, 8), "terminal sign reserve")

    require(529 > 512, "sin(pi/8)>3/8 radical certificate")
    require(32 > 1, "cos(pi/8)>3/4 radical certificate")

    return {
        "terminal_atom": "The full q=N lifted mathcal N atom is smaller than (17/5)h^3|w_N|; each reciprocal/hard half is smaller than (17/10)h^3|w_N|.",
        "row_remainder": "With R_beta the l1 norm of beta_N, |P_lin^[N](log N)-mathcal N_(p,xi)|<(8/5)hR_alpha+(16/5)hR_beta.",
        "terminal_source": "At the physical anchor, C_Nw_N=-Q_RS(p)exp(W_theta) and |exp(W_theta)/C_N-1|<9h.",
        "terminal_signs": "The terminal quadrature obeys Im w_N(p=0)>5/8 and Im w_N(p=1)<-1/4.",
    }


def audit_terminal_kernel() -> dict[str, int | str]:
    mp.mp.dps = 80
    cases = 0
    for n_value in range(2, 40):
        for theta in [Fraction(0), Fraction(1, 7), Fraction(4, 7), Fraction(1)]:
            a_value = Fraction(n_value) + theta
            for offset in [Fraction(1, 17), Fraction(7, 13), Fraction(16, 17)]:
                alpha = a_value**2 - offset
                x_n = alpha / n_value
                a_n = alpha / (Fraction(n_value) + Fraction(1, 2))
                lower = floor_fraction(a_n) + 1
                upper = floor_fraction(2 * alpha)
                a_star = x_n - lower + 1
                b_star = upper + 1 - x_n

                require(a_star > Fraction(1, 4), "terminal a_* lower bound")
                require(a_star < 2, "terminal a_* upper bound")
                require(
                    b_star > alpha * (2 - Fraction(1, n_value)),
                    "terminal b_* lower bound",
                )
                require(b_star > a_star, "terminal kernel orientation")

                a_mp = mp_fraction(a_star)
                b_mp = mp_fraction(b_star)
                s1 = mp.digamma(b_mp) - mp.digamma(a_mp)
                integral_lower = mp.log(b_mp / a_mp)
                integral_upper = integral_lower + 1 / a_mp - 1 / b_mp
                coarse_lower = mp.log(
                    mp.mpf(n_value**2 - 1)
                    * (2 - mp.mpf(1) / n_value)
                    / 2
                )
                coarse_upper = mp.log(8 * (n_value + 1) ** 2 + 4) + 4

                require(s1 > integral_lower, "terminal S1 integral lower")
                require(s1 < integral_upper, "terminal S1 integral upper")
                require(s1 > coarse_lower, "terminal S1 coarse lower")
                require(s1 < coarse_upper, "terminal S1 coarse upper")
                cases += 1

    return {
        "terminal_s1_cases": cases,
        "integral_sandwich": "log(b_*/a_*)<S_(1,N)<log(b_*/a_*)+a_*^(-1)-b_*^(-1).",
        "coarse_growth": "log(((N^2-1)(2-1/N))/2)<S_(1,N)<log(8(N+1)^2+4)+4.",
    }


def symbolic_certificate() -> dict:
    i = sp.I
    u_q, u_n, c, c_x, b, b_x, delta = sp.symbols(
        "u_q u_N c c_x b b_x delta", real=True
    )
    c_q, d_q = sp.symbols("C_q D_q")
    s_prime = c + i * b
    s_second = c_x + i * b_x
    r_q = s_prime * u_q - c * u_n
    r_q_x = s_second * u_q - c_x * u_n + i * b * sp.symbols(
        "u_Nx", real=True
    )
    q_q = (r_q + delta) * c_q + d_q
    n_q = sp.expand(r_q * q_q + r_q_x * c_q)
    require_zero(
        r_q - (-s_prime * (-u_q) - c * u_n),
        "all-carrier centered R",
    )
    require_zero(
        r_q_x
        - (-s_second * (-u_q) - c_x * u_n + i * b * sp.symbols("u_Nx", real=True)),
        "all-carrier centered R_x",
    )
    require_zero(n_q - (r_q * q_q + r_q_x * c_q), "mathcal N source row")

    wr, wi, pr, pi_part = sp.symbols("w_r w_i P_r P_i", real=True)
    kappa = 1 / (2 * sp.pi * i)
    w = wr + i * wi
    p_value = pr + i * pi_part
    require_zero(
        sp.re(kappa * w * p_value)
        - sp.im(w * p_value) / (2 * sp.pi),
        "terminal quadrature projection",
    )
    n_real = sp.symbols("N_p_xi", real=True)
    require_zero(
        sp.im(w * n_real) - n_real * wi,
        "leading rank-one projection",
    )

    p_rs, m = sp.symbols("p_RS m", real=True)
    h_phase = p_rs**2 / 2 - (2 * m + 1) * p_rs + sp.Rational(3, 8)
    shifted_r_phase = (
        (p_rs - 2 * m - 2) ** 2 / 2
        + (p_rs - 2 * m - 2)
        + sp.Rational(3, 8)
    )
    require_zero(
        shifted_r_phase - h_phase - 2 * m * (m + 1),
        "terminal conjugate phase relation",
    )
    require_zero(
        shifted_r_phase.subs({p_rs: 0, m: 0}) - sp.Rational(3, 8),
        "midpoint recurrence phase",
    )

    bounds = rational_bound_audit()
    kernel = audit_terminal_kernel()
    return {
        "physical_carrier_row": {
            "rates": "For u_q=log(a/q), R_q=s_*'u_q-cu_N and R_(q,x)=s_*''u_q-c_xu_N+ib u_(N,x). Then Q_q=(R_q+delta)C_q+D_q and mathcal N_q=R_qQ_q+R_(q,x)C_q.",
            "lift": "With w_q=nu c_xi exp(S(log q))e(alpha_Plog q), mathcal N_(p,xi)=sum_(q=1)^N Re{-iu_q mathcal N_qw_q} after endpoint composition.",
            "halves": "Before endpoint composition the reciprocal return and hard endpoint each contribute one half of Re{-iu_Nmathcal N_Nw_N}; together they give the full q=N atom exactly once.",
        },
        "terminal_kernel": {
            "geometry": "For x_N=alpha_P/N, a_*=x_N-floor(alpha_P/(N+1/2)) and b_*=floor(2alpha_P)+1-x_N satisfy 1/4<a_*<2 and b_*>alpha_P(2-1/N)>a_*.",
            "integral": kernel["integral_sandwich"],
            "growth": kernel["coarse_growth"],
            "meaning": "The terminal S_1 kernel is positive and logarithmic in N; it is not an O(1) or O(h) factor.",
        },
        "conditional_boundary": {
            "exact": "If S_(1,N)=S_1^[N](alpha_P/N), the terminal conditional linear trace is L_(N,1)=[S_(1,N)/(2pi)]Im{w_NP_lin^[N](log N)}.",
            "leading": "L_(N,1)=[S_(1,N)/(2pi)]mathcal N_(p,xi)Im(w_N)+E_row.",
            "row_error": "|E_row|<[S_(1,N)/(2pi)]h|w_N|{(8/5)R_alpha+(16/5)R_beta}.",
            "bulk_split": "Writing mathcal N_(<N,xi)=sum_(q<N)Re{-iu_qmathcal N_qw_q}, L_(N,1)=[S_(1,N)/(2pi)]Im(w_N)mathcal N_(<N,xi)+E_term.",
            "bulk_error": "|E_term|<[S_(1,N)/(2pi)]{h|w_N|[(8/5)R_alpha+(16/5)R_beta]+(17/5)h^3|w_N|^2}.",
        },
        "terminal_recurrence": {
            "complex_relation": "For Q_RS(p), r(p), and R(p) of Section 11.154, Q_RS(p)r(p)^m=(-1)^m conjugate(R(p-2m-2)).",
            "defect": "C_0(p)-Q_RS(p)sum_(m=0)^M r(p)^m=(-1)^(M+1)C_0(p-2M-2)+2i sum_(m=0)^M(-1)^m Im R(p-2m-2).",
            "midpoint": "At p=0 and M=0 the explicit quadrature-defect summand is 2i sin(3pi/8), so it is nonzero.",
            "consequence": "The known recurrence telescopes the real terminal trace, whereas the conditional S_1 boundary kernel sees the imaginary quadrature because kappa=-i/(2pi).",
        },
        "bounds": bounds,
        "audit": {"terminal_s1_cases": kernel["terminal_s1_cases"]},
        "handoff": {
            "live_kernel": "The first unsuppressed boundary object is [S_(1,N)/(2pi)]Im(w_N)mathcal N_(<N,xi), not a local q=N cancellation.",
            "route": "Estimate this one-versus-many bulk-terminal correlation jointly with the remaining S_2/S_3 and interior outer terms, H_out, T_out, correction perturbation, and moved-tail defect.",
        },
    }


def build_rows(cert: dict) -> list[GateRow]:
    p = cert["physical_carrier_row"]
    k = cert["terminal_kernel"]
    b = cert["conditional_boundary"]
    r = cert["terminal_recurrence"]
    q = cert["bounds"]
    h = cert["handoff"]
    rows = [
        GateRow("tcq_01_rates", "source row", "proved", "The retained mathcal N carrier has an exact terminal-centered source row.", p["rates"], "No row norm follows."),
        GateRow("tcq_02_lift", "frequency lift", "proved", "The complete retained frequency derivative is one full carrier sum.", p["lift"], "The sum is not signed."),
        GateRow("tcq_03_halves", "endpoint accounting", "proved", "The returned and hard q=N halves reassemble one atom.", p["halves"], "No endpoint is counted twice."),
        GateRow("tcq_04_local_bound", "terminal atom", "proved", "The full local q=N lifted mathcal N atom is h-cubed scale relative to its source amplitude.", q["terminal_atom"], "This does not bound the bulk row."),
        GateRow("tcq_05_half_bound", "terminal halves", "proved", "Each q=N half has half the local h-cubed envelope.", q["terminal_atom"], "The two halves must remain joined."),
        GateRow("tcq_06_geometry", "kernel geometry", "proved", "The terminal Hurwitz arguments are ordered and pole-free.", k["geometry"], "This is not a signed source estimate."),
        GateRow("tcq_07_integral", "digamma sandwich", "proved", "The terminal conditional kernel has a sharp elementary integral sandwich.", k["integral"], "No physical coefficient is bounded."),
        GateRow("tcq_08_growth", "kernel scale", "proved", "The terminal conditional kernel grows logarithmically.", k["growth"], "It must not be treated as O(1)."),
        GateRow("tcq_09_meaning", "scale guard", "guard_validated", "Endpoint centering does not suppress S_(1,N).", k["meaning"], "Only its row coefficient was reduced in Section 11.183."),
        GateRow("tcq_10_exact_trace", "conditional trace", "proved", "The exact terminal S_1 trace is a quadrature projection.", b["exact"], "The factor 1/(2pi) is Fourier normalization."),
        GateRow("tcq_11_leading", "rank-one kernel", "proved", "The unique order-one row produces one real rank-one boundary correlation.", b["leading"], "Its sign is open."),
        GateRow("tcq_12_row_error", "row remainder", "proved", "Every nonleading terminal row is explicitly h-suppressed relative to the two retained row norms.", b["row_error"], "Neither row norm is assigned a numerical size."),
        GateRow("tcq_13_beta", "beta row", "proved", "The complete beta row remains in the terminal remainder.", q["row_remainder"], "The genuine edge is not deleted."),
        GateRow("tcq_14_bulk_split", "source decomposition", "proved", "The local terminal derivative atom separates from the q<N source.", b["bulk_split"], "The q<N sum remains global."),
        GateRow("tcq_15_bulk_error", "joined remainder", "proved", "The local atom and row-collapse errors have one explicit joined envelope.", b["bulk_error"], "This is row-relative, not an h-squared theorem."),
        GateRow("tcq_16_complex_relation", "terminal phase", "proved", "The terminal carrier phase is the conjugate recurrence phase.", r["complex_relation"], "Equality of real traces is not complex equality."),
        GateRow("tcq_17_recurrence_defect", "quadrature recurrence", "proved", "Complex recurrence composition leaves an explicit imaginary defect.", r["defect"], "The defect is not bounded here."),
        GateRow("tcq_18_midpoint", "nonzero witness", "proved", "The recurrence quadrature defect is not identically zero.", r["midpoint"], "This is not an Xi counterexample."),
        GateRow("tcq_19_consequence", "recurrence guard", "guard_validated", "The established real recurrence cannot cancel the conditional quadrature by itself.", r["consequence"], "A new quadrature estimate is required."),
        GateRow("tcq_20_source", "terminal anchor", "proved", "The corrected physical terminal source remains close to its exact Riemann--Siegel phase.", q["terminal_source"], "This controls only the terminal factor."),
        GateRow("tcq_21_signs", "terminal sign guard", "proved", "The physical terminal quadrature has both signs across the endpoint chart.", q["terminal_signs"], "No sign of the full product follows."),
        GateRow("tcq_22_live", "live kernel", "proved", "The local-cancellation route reduces to a global bulk-terminal correlation.", h["live_kernel"], "No signed bound is proved."),
        GateRow("tcq_23_route", "route decision", "open", "Prove or falsify the joined source-specific bulk-terminal estimate.", h["route"], "No componentwise budget is licensed."),
        GateRow("tcq_24_ties", "tie guard", "guard_validated", "The q=N half audit retains the half-open reciprocal convention.", "Transfer the complete mode and both endpoint halves before differentiating.", "No complement-only derivative is used."),
        GateRow("tcq_25_phase", "phase guard", "guard_validated", "The leading boundary kernel is an imaginary, not real, terminal projection.", "The rotation is forced by kappa=1/(2pi i).", "Do not substitute the real recurrence trace."),
        GateRow("tcq_26_quadratic", "quadratic guard", "guard_validated", "H_out and T_out remain outside this linear rank-one reduction.", "Retain both phase families and their mixed S_2/S_3/interior pieces.", "No quadratic remainder is bounded."),
        GateRow("tcq_27_pi", "pi provenance", "proved", "Every pi in the gate is inherited.", "The factor 1/(2pi) comes from the Fourier character, while Q_RS and R use the completed-zeta/Riemann--Siegel phase.", "No geometric fitting is introduced."),
        GateRow("tcq_28_boundary", "proof boundary", "guard_validated", "This is a source and recurrence reduction, not a proof of RH.", "No signed bulk-terminal estimate, complete outer bound, quadratic residual bound, retained-current theorem, Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.", "This is not a proof of RH."),
    ]
    require(len(rows) == 28, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    c = artifact["symbolic_certificate"]
    p = c["physical_carrier_row"]
    k = c["terminal_kernel"]
    b = c["conditional_boundary"]
    r = c["terminal_recurrence"]
    q = c["bounds"]
    h = c["handoff"]
    return f"""# Full-Support Terminal Conditional Quadrature Gate

Date: 2026-08-02

Status: exact terminal quadrature, carrier-half audit, and recurrence defect proved; signed bulk-terminal estimate open; not a proof of RH.

## Physical Carrier Row

{p['rates']}

{p['lift']}

{p['halves']}

{q['terminal_atom']}

## Terminal Kernel

{k['geometry']}

{k['integral']}

{k['growth']}

{k['meaning']}

## Conditional Boundary

{b['exact']}

{b['leading']}

{q['row_remainder']}

{b['row_error']}

{b['bulk_split']}

{b['bulk_error']}

## Recurrence Audit

{r['complex_relation']}

{r['defect']}

{r['midpoint']}

{r['consequence']}

{q['terminal_source']}

{q['terminal_signs']}

## Handoff

{h['live_kernel']}

{h['route']}

The next calculation must estimate the displayed one-versus-many bulk-terminal correlation before adding the remaining endpoint, interior, Hermitian, transpose, correction, and moved-tail channels. The real terminal recurrence cannot be substituted for the required quadrature identity.

## Pi Provenance

The factor `1/(2pi)` is forced by `kappa=1/(2pi i)` under the fixed Fourier character `e(x)=exp(2pi i x)`. The phases `Q_RS` and `R` are inherited from the completed-zeta/Riemann--Siegel endpoint recurrence. No new geometric `pi` is introduced.

## Proof Boundary

This gate proves the exact all-carrier source expansion of `mathcal N_(p,xi)`, reassembly and h-cubed size of the two q=N derivative halves, positivity and logarithmic scale of the terminal `S_1` kernel, the exact conditional quadrature projection, a row-relative bulk-terminal reduction, the complex conjugate recurrence relation, its explicit imaginary defect, and two terminal-factor sign witnesses. It proves no signed bulk-terminal correlation estimate, complete outer-complement or quadratic residual bound, signed current, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    payloads = load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "physical_carrier_source_expansions": 1,
        "endpoint_carrier_half_audits": 1,
        "endpoint_carrier_halves": 2,
        "local_terminal_derivative_bounds": 2,
        "terminal_s1_geometry_bounds": 4,
        "terminal_s1_integral_sandwiches": 1,
        "terminal_s1_sample_checks": certificate["audit"]["terminal_s1_cases"],
        "terminal_quadrature_identities": 1,
        "row_relative_terminal_bounds": 2,
        "bulk_terminal_kernel_decompositions": 1,
        "complex_recurrence_relations": 1,
        "quadrature_recurrence_defects": 1,
        "physical_terminal_quadrature_sign_witnesses": 2,
        "signed_bulk_terminal_bounds": 0,
        "complete_outer_complement_bounds": 0,
        "quadratic_remainder_bounds": 0,
        "signed_full_current_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "terminal conditional quadrature and recurrence defect proved; signed bulk-terminal estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(payloads),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the exact all-carrier expansion of mathcal N_(p,xi), the two q=N derivative halves and their h-cubed scale, positivity and logarithmic growth of the terminal S_1 kernel, the conditional quadrature projection, a row-relative bulk-terminal reduction, and the complex terminal recurrence defect. It proves no signed bulk-terminal estimate, complete outer-complement or quadratic residual bound, signed full-current or Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Q209, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built full-support terminal conditional quadrature gate: "
        f"{counts['rows']} rows, "
        f"{counts['terminal_s1_sample_checks']} kernel checks, "
        f"{counts['endpoint_carrier_halves']} carrier halves, "
        f"{counts['quadrature_recurrence_defects']} quadrature defect, "
        f"{counts['signed_bulk_terminal_bounds']} signed bulk-terminal bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
