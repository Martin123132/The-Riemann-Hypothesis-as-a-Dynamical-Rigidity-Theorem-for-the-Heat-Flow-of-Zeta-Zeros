#!/usr/bin/env python3
"""Build the reciprocal-Morse first stationary-correction cancellation gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
    "ideal_cubic_reciprocal_first_correction_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "endpoint_composition": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_endpoint_composition_"
        "retention_gate.json"
    ),
    "physical_plin": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "amplitude_gate.json"
    ),
    "reciprocal": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_gate.json"
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


def require_zero(expression: sp.Expr, label: str) -> None:
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
    if (
        payloads["endpoint_composition"]["counts"][
            "endpoint_composed_value_identities"
        ]
        != 6
    ):
        raise RuntimeError("endpoint-composition source drifted")
    if payloads["physical_plin"]["counts"]["exact_morse_cprime_formulas"] != 2:
        raise RuntimeError("physical P_lin source drifted")
    reciprocal_counts = payloads["reciprocal"]["counts"]
    if reciprocal_counts["reciprocal_phase_identities"] != 6:
        raise RuntimeError("reciprocal phase source drifted")
    if reciprocal_counts["grouped_ideal_cubic_bounds"] != 0:
        raise RuntimeError("reciprocal proof boundary drifted")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_certificate() -> dict:
    delta = sp.symbols("delta", real=True)
    alpha, q = sp.symbols("alpha_P q", positive=True, real=True)
    a_0, a_1, a_2 = sp.symbols("A_0 A_1 A_2")

    h = delta - sp.log(1 + delta)
    z_series = sp.series(
        delta * sp.sqrt(sp.simplify(2 * h / delta**2)), delta, 0, 4
    ).removeO()
    j_series = sp.series((1 + delta) * z_series / delta, delta, 0, 3).removeO()
    require_zero(j_series.subs(delta, 0) - 1, "J(1)")
    require_zero(sp.diff(j_series, delta).subs(delta, 0) - sp.Rational(2, 3), "J'(1)")
    require_zero(
        sp.diff(j_series, delta, 2).subs(delta, 0) + sp.Rational(5, 18),
        "J''(1)",
    )

    j_jet = 1 + sp.Rational(2, 3) * delta - sp.Rational(5, 36) * delta**2

    def eta_derivative(expression: sp.Expr) -> sp.Expr:
        return sp.expand(j_jet * sp.diff(expression, delta) / sp.sqrt(alpha))

    amplitude_u = a_0 + q * a_1 * delta + q**2 * a_2 * delta**2 / 2
    b_u = q / sp.sqrt(alpha) * amplitude_u * j_jet
    b_u_second = sp.simplify(
        eta_derivative(eta_derivative(b_u)).subs(delta, 0)
    )
    differential_symbol = q**2 * a_2 + 2 * q * a_1 + a_0 / 6
    require_zero(
        b_u_second - q * differential_symbol / alpha ** sp.Rational(3, 2),
        "inner second derivative",
    )

    rho = 1 + delta
    shift = sp.series(q / rho - q, delta, 0, 3).removeO()
    amplitude_dual = sp.expand(a_0 + a_1 * shift + a_2 * shift**2 / 2)
    amplitude_dual = sp.series(amplitude_dual, delta, 0, 3).removeO()
    d_dual = sp.series(j_jet * amplitude_dual / rho, delta, 0, 3).removeO()
    d_dual_second = sp.simplify(
        eta_derivative(eta_derivative(d_dual)).subs(delta, 0)
    )
    require_zero(
        d_dual_second - differential_symbol / alpha,
        "dual second derivative",
    )
    transported_inner_second = sp.simplify(
        b_u_second * sp.sqrt(alpha) / q
    )
    require_zero(
        transported_inner_second - d_dual_second,
        "sequential second-jet identity",
    )

    lam, t, sigma = sp.symbols("lambda t sigma", real=True)
    p = sp.Function("P")(lam)
    s = t * lam**2 / 4 - sigma * lam
    weighted = sp.exp(s) * p
    differential_lambda = sp.diff(weighted, lam, 2) + sp.diff(weighted, lam) + weighted / 6
    g = sp.diff(s, lam)
    expected_lambda = sp.exp(s) * (
        sp.diff(p, lam, 2)
        + (2 * g + 1) * sp.diff(p, lam)
        + (g**2 + g + t / 2 + sp.Rational(1, 6)) * p
    )
    require_zero(
        differential_lambda - expected_lambda,
        "physical saddle differential operator",
    )

    kappa = sp.symbols("kappa")
    correction_sum = -kappa * d_dual_second / 2 + kappa * transported_inner_second / 2
    require_zero(correction_sum, "first sequential correction cancellation")

    rho_symbol = sp.symbols("rho", positive=True, real=True)
    reciprocal_phase = alpha * sp.log(q) + alpha * (
        rho_symbol - 1 - sp.log(rho_symbol)
    )
    r_symbol = alpha * rho_symbol / q
    original_phase = alpha * (sp.log(alpha / r_symbol) - 1) + q * r_symbol
    require_zero(original_phase - reciprocal_phase, "reciprocal Morse phase")

    leading_u = q * a_0 / sp.sqrt(alpha)
    leading_dual_jacobian = sp.sqrt(alpha) / q
    require_zero(
        leading_u * leading_dual_jacobian - a_0,
        "leading carrier amplitude",
    )

    f_minus = sp.exp(-sp.I * sp.pi / 4)
    f_plus = sp.exp(sp.I * sp.pi / 4)
    require_zero(f_minus * f_plus - 1, "Fresnel signature product")

    return {
        "cell_partition": {
            "definition": (
                "C_q(alpha_P)={r in Z_(>0): q-1/2<=alpha_P/r<q+1/2}."
            ),
            "integer_bounds": (
                "C_q(alpha_P)={floor(alpha_P/(q+1/2))+1,...,"
                "floor(alpha_P/(q-1/2))}."
            ),
            "tie_rule": (
                "The lower derivative boundary is included and the upper is "
                "excluded, so adjacent cells are disjoint and every equality "
                "alpha_P/r=q+1/2 transfers to C_(q+1)."
            ),
            "coverage": (
                "The cells q>=1 partition every positive mode with alpha_P/r>=1/2; "
                "intersection with the fixed roster handles its at-most-one top edge."
            ),
        },
        "dual_morse": {
            "phase": (
                "For rho=qr/alpha_P, Phi_q(r)=alpha_P[log(alpha_P/r)-1]+qr="
                "alpha_Plog q+alpha_P[rho-1-log rho]."
            ),
            "coordinate": (
                "zeta(rho)=sgn(rho-1)sqrt(2[rho-1-log rho]), "
                "eta=sqrt(alpha_P)zeta, so Phi_q=alpha_Plog q+eta^2/2."
            ),
            "jacobian": (
                "dr/deta=[sqrt(alpha_P)/q]J(rho), with "
                "J(1)=1, J'(1)=2/3, and J''(1)=-5/18."
            ),
            "signature": (
                "The original u chart has F_-=e(-1/8), the reciprocal r chart "
                "has F_+=e(+1/8), and F_-F_+=1."
            ),
            "leading_carrier": (
                "At r_*=alpha_P/q, [qA(q)/sqrt(alpha_P)]"
                "[sqrt(alpha_P)/q]=A(q), so the sequential leading stationary "
                "symbol is exactly e(alpha_Plog q)A(q)."
            ),
        },
        "first_correction": {
            "operator": (
                "mathcal D_qA=q^2A''(q)+2qA'(q)+A(q)/6."
            ),
            "inner_second": (
                "For b_r(y)=[sqrt(alpha_P)/r]A((alpha_P/r)v)J(v), "
                "b_(alpha_P/q)''(0)=q mathcal D_qA/alpha_P^(3/2)."
            ),
            "dual_amplitude": (
                "The reciprocal leading amplitude after its Jacobian is "
                "d_q(eta)=J(rho)A(q/rho)/rho, with "
                "d_q''(0)=mathcal D_qA/alpha_P."
            ),
            "transported_inner": (
                "[sqrt(alpha_P)/q]b_(alpha_P/q)''(0)="
                "d_q''(0)=mathcal D_qA/alpha_P."
            ),
            "physical_operator": (
                "For A(q)=exp(S(log q))P(log q), mathcal D_qA="
                "exp(S)[P''+(2g+1)P'+(g^2+g+t/2+1/6)P](log q), "
                "the same removable-saddle operator as Section 11.173."
            ),
            "cancellation": (
                "The positive reciprocal Gaussian contributes -kappa d_q''(0)/2, "
                "while the leading transform of the negative-chart c-prime "
                "correction contributes +kappa[sqrt(alpha_P)/q]b''(0)/2; "
                "their sum is exactly zero."
            ),
        },
        "series_audit": {
            "zeta": str(z_series),
            "J": str(j_series),
            "inner_second": str(b_u_second),
            "dual_second": str(d_dual_second),
            "differential_symbol": str(differential_symbol),
        },
    }


def boundary_certificate() -> dict:
    return {
        "achievement": (
            "The first full-line sequential stationary correction cancels "
            "coefficientwise for an arbitrary twice-differentiable amplitude. "
            "For the physical ideal cubic it is exactly the cancellation of the "
            "reciprocal correction to mathcal F against the saddle symbol of mathcal J."
        ),
        "finite_cell_boundary": (
            "The proof uses the full-line Gaussian moment coefficients at one "
            "continuous reciprocal saddle. The exact discrete half-open cell has "
            "integer boundary transfers which are not bounded here."
        ),
        "incomplete_fresnel_boundary": (
            "The physical u integral has incomplete Fresnel endpoints. Their tails, "
            "the lower and upper endpoint packages, and boundary q=1,B cells must "
            "be recomposed before the full-line cancellation becomes a finite theorem."
        ),
        "next_defect": (
            "In the complete-line interior symbol the order-alpha_P^(-1) correction "
            "vanishes, so the next formal interior term begins at second correction "
            "order. No O(alpha_P^(-2)) finite-cell remainder is claimed until the "
            "discrete boundaries and incomplete Fresnel tails are bounded."
        ),
        "handoff": (
            "Derive an exact finite reciprocal-cell identity with the two endpoint "
            "transfers exposed. Subtract the returned physical carrier and the "
            "cancelling first correction, then bound the remaining cell defect in a "
            "form that is summable across q without absolute family aggregation."
        ),
    }


def build_rows(symbolic: dict, boundary: dict) -> list[GateRow]:
    cells = symbolic["cell_partition"]
    dual = symbolic["dual_morse"]
    first = symbolic["first_correction"]
    return [
        GateRow("irfc_01_cells", "reciprocal_cells", "certified", "Nearest-integer phase derivatives define deterministic reciprocal cells.", cells["definition"] + " " + cells["integer_bounds"], "Only positive modes with derivative at least 1/2 are covered."),
        GateRow("irfc_02_ties", "reciprocal_cells", "certified", "The half-open convention resolves every exact derivative tie.", cells["tie_rule"], "Cell membership is combinatorial, not a phase estimate."),
        GateRow("irfc_03_coverage", "reciprocal_cells", "certified", "The reciprocal cells partition the stationary positive roster range.", cells["coverage"], "Physical endpoint cells remain separate."),
        GateRow("irfc_04_phase", "dual_morse", "certified", "Each q cell has an exact positive-sign Morse phase.", dual["phase"] + " " + dual["coordinate"], "The discrete sum has not been replaced by an integral."),
        GateRow("irfc_05_jacobian", "dual_morse", "certified", "The reciprocal chart uses the same Morse Jacobian as the original transform.", dual["jacobian"], "Only its saddle jet is used in the correction identity."),
        GateRow("irfc_06_signature", "dual_morse", "certified", "The sequential Fresnel signatures cancel exactly.", dual["signature"], "Incomplete endpoint Fresnel factors are not full signatures."),
        GateRow("irfc_07_leading", "carrier_return", "certified", "The sequential leading stationary symbol returns the physical q carrier with unit amplitude.", dual["leading_carrier"], "This is the full-line local symbol, not a discrete cell theorem."),
        GateRow("irfc_08_operator", "first_correction", "certified", "One differential operator controls both first stationary corrections.", first["operator"], "The amplitude must have the displayed derivatives."),
        GateRow("irfc_09_inner", "first_correction", "certified", "The original negative Morse chart has an exact saddle second derivative.", first["inner_second"], "This is a coefficient identity."),
        GateRow("irfc_10_dual", "first_correction", "certified", "The reciprocal leading amplitude has the matching second derivative.", first["dual_amplitude"], "This is a coefficient identity."),
        GateRow("irfc_11_transport", "first_correction", "certified", "The reciprocal Jacobian transports the inner second derivative exactly to the dual one.", first["transported_inner"], "No remainder estimate follows by itself."),
        GateRow("irfc_12_physical", "physical_operator", "certified", "For the physical amplitude the common differential symbol is exactly the removable saddle operator.", first["physical_operator"], "The joined alpha/beta row remains retained."),
        GateRow("irfc_13_cancel", "first_correction_cancellation", "certified", "The first sequential stationary correction cancels identically.", first["cancellation"], "This is a full-line stationary-symbol cancellation."),
        GateRow("irfc_14_achievement", "reciprocal_handoff", "certified", "The cancellation acts between the transformed retained main and the c-prime saddle symbol.", boundary["achievement"], "It has not yet been promoted through finite cell boundaries."),
        GateRow("irfc_15_cell_boundary", "finite_cell_boundary", "open", "Exact integer cell boundary transfers remain to be bounded.", boundary["finite_cell_boundary"], "No finite-cell remainder bound is claimed."),
        GateRow("irfc_16_fresnel_boundary", "endpoint_boundary", "open", "Incomplete Fresnel and physical endpoint packages remain to be recomposed.", boundary["incomplete_fresnel_boundary"], "No endpoint-complete h2 estimate is claimed."),
        GateRow("irfc_17_next", "next_defect", "open", "The next formal interior correction is delayed, but its finite counterpart is still open.", boundary["next_defect"] + " " + boundary["handoff"], "No alpha_P^-2 finite-cell theorem is claimed."),
        GateRow("irfc_18_boundary", "proof_boundary", "guard_validated", "The exact symbol cancellation is not promoted to the RH bridge.", "Finite cells, incomplete endpoints, q-family summation, physical normalization, and quadratic residuals remain open.", "No grouped h2, signed-flow, Xi, Lambda<=0, PF-infinity, RH, or prize theorem is claimed."),
    ]


def render_note(artifact: dict) -> str:
    symbolic = artifact["symbolic_certificate"]
    boundary = artifact["boundary_certificate"]
    cells = symbolic["cell_partition"]
    dual = symbolic["dual_morse"]
    first = symbolic["first_correction"]
    return "\n".join(
        [
            "# Ideal-Cubic Reciprocal First-Correction Gate",
            "",
            "Date: 2026-08-02",
            "",
            "Status: exact reciprocal cells, dual positive Morse chart, returned leading carrier, and cancellation of the first full-line sequential stationary correction. This is not a proof of a finite-cell h^2 estimate, signed flow, or RH.",
            "",
            "```text",
            str(RESULT_PATH.relative_to(REPO_ROOT)).replace("\\", "/"),
            f"python {str(Path(__file__).relative_to(REPO_ROOT)).replace(chr(92), '/')}",
            "```",
            "",
            "## Reciprocal Cells",
            "",
            "```text",
            cells["definition"],
            cells["integer_bounds"],
            cells["tie_rule"],
            cells["coverage"],
            "```",
            "",
            "## Dual Morse Chart",
            "",
            "```text",
            dual["phase"],
            dual["coordinate"],
            dual["jacobian"],
            dual["signature"],
            dual["leading_carrier"],
            "```",
            "",
            "## First Correction",
            "",
            "```text",
            first["operator"],
            first["inner_second"],
            first["dual_amplitude"],
            first["transported_inner"],
            first["physical_operator"],
            first["cancellation"],
            "```",
            "",
            boundary["achievement"],
            "",
            "## Remaining Boundary",
            "",
            boundary["finite_cell_boundary"],
            "",
            boundary["incomplete_fresnel_boundary"],
            "",
            boundary["next_defect"],
            "",
            "## Handoff",
            "",
            boundary["handoff"],
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
    boundary = boundary_certificate()
    rows = build_rows(symbolic, boundary)
    return {
        "kind": KIND,
        "date": "2026-08-02",
        "status": (
            "exact reciprocal q cells, dual positive Morse chart, returned leading "
            "carrier, and first full-line sequential stationary-correction "
            "cancellation complete; finite-cell and incomplete-endpoint defects open; "
            "no grouped h2, signed-flow, Xi, or RH theorem"
        ),
        "proof_boundary": (
            "This gate proves the exact half-open reciprocal cells, dual Morse phase "
            "and Jacobian, leading carrier return, equality of the two first "
            "stationary differential symbols, and their coefficientwise cancellation. "
            "It proves no finite reciprocal-cell remainder, incomplete-Fresnel or "
            "endpoint-complete h2 bound, q-family cancellation theorem, quadratic "
            "residual bound, signed flow estimate, Phi_B upper bound, contact "
            "exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, "
            "RH, or prize-level conclusion."
        ),
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "boundary_certificate": boundary,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "reciprocal_cell_identities": 4,
            "dual_morse_identities": 5,
            "first_correction_identities": 5,
            "first_correction_cancellations": 1,
            "finite_cell_remainder_bounds": 0,
            "incomplete_fresnel_bounds": 0,
            "endpoint_composed_h2_bounds": 0,
            "q_family_bounds": 0,
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
        "built reciprocal first-correction gate: "
        f"{counts['rows']} rows, "
        f"{counts['reciprocal_cell_identities']} cell identities, "
        f"{counts['dual_morse_identities']} dual-Morse identities, "
        f"{counts['first_correction_identities']} first-correction identities, "
        f"{counts['first_correction_cancellations']} cancellation, "
        f"{counts['finite_cell_remainder_bounds']} finite-cell bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
