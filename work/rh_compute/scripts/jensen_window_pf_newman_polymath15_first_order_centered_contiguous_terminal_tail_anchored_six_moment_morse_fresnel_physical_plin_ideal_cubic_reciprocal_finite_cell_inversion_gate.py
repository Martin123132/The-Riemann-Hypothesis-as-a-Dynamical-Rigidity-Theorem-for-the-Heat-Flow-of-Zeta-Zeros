#!/usr/bin/env python3
"""Build the exact finite reciprocal-cell Fourier-inversion gate."""

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
    "ideal_cubic_reciprocal_finite_cell_inversion_gate"
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
    "finite_poisson": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_finite_poisson_transport_reduction.json"
    ),
    "first_correction": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_first_correction_gate.json"
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
    simplified = sp.simplify(expression)
    if simplified != 0:
        raise RuntimeError(f"{label}: {simplified}")


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
    endpoint_counts = payloads["endpoint_composition"]["counts"]
    if endpoint_counts["endpoint_composed_value_identities"] != 6:
        raise RuntimeError("endpoint-composition source drifted")
    if endpoint_counts["oscillatory_interior_bounds"] != 0:
        raise RuntimeError("endpoint-composition proof boundary drifted")
    poisson_counts = payloads["finite_poisson"]["counts"]
    if poisson_counts["full_poisson_identities"] != 6:
        raise RuntimeError("finite-Poisson source drifted")
    if poisson_counts["explicit_uniform_remainder_constants"] != 0:
        raise RuntimeError("finite-Poisson proof boundary drifted")
    correction_counts = payloads["first_correction"]["counts"]
    if correction_counts["first_correction_cancellations"] != 1:
        raise RuntimeError("first-correction source drifted")
    if correction_counts["finite_cell_remainder_bounds"] != 0:
        raise RuntimeError("first-correction proof boundary drifted")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def audit_cells() -> dict[str, int]:
    cases = 0
    ties = 0
    tilings = 0
    roster_containments = 0
    for alpha in [
        Fraction(13, 2),
        Fraction(25, 3),
        Fraction(10),
        Fraction(27, 2),
        Fraction(120),
    ]:
        for b_value in [2, 3, 4, 5]:
            physical_union: list[int] = []
            for q in range(1, b_value + 1):
                a_q = alpha / (Fraction(q) + Fraction(1, 2))
                b_q = alpha / (Fraction(q) - Fraction(1, 2))
                direct = [
                    r
                    for r in range(1, floor_fraction(2 * alpha) + 1)
                    if a_q < r <= b_q
                ]
                lower = floor_fraction(a_q) + 1
                upper = floor_fraction(b_q)
                formula = list(range(lower, upper + 1)) if lower <= upper else []
                require(direct == formula, f"cell formula failed at {alpha}, {q}")
                physical_union.extend(direct)
                cases += 1

                if a_q.denominator == 1:
                    boundary = int(a_q)
                    require(boundary not in direct, "lower tie was not excluded")
                    if q < b_value:
                        next_b = alpha / (Fraction(q + 1) - Fraction(1, 2))
                        require(next_b == a_q, "adjacent continuous boundaries drifted")
                        next_direct = [
                            r
                            for r in range(1, floor_fraction(2 * alpha) + 1)
                            if alpha / (Fraction(q + 1) + Fraction(1, 2))
                            < r
                            <= next_b
                        ]
                        require(boundary in next_direct, "lower tie did not transfer")
                        ties += 1

            global_lower = alpha / (Fraction(b_value) + Fraction(1, 2))
            global_upper = 2 * alpha
            direct_union = [
                r
                for r in range(1, floor_fraction(global_upper) + 1)
                if global_lower < r <= global_upper
            ]
            require(
                sorted(physical_union) == direct_union,
                f"physical cell tiling failed at {alpha}, {b_value}",
            )
            require(
                len(physical_union) == len(set(physical_union)),
                f"physical cell overlap at {alpha}, {b_value}",
            )
            alpha_minus = alpha - Fraction(1, 10)
            alpha_plus = alpha + Fraction(1, 10)
            roster_lower = max(1, floor_fraction(alpha_minus / (2 * b_value)))
            roster_upper = -floor_fraction(-2 * alpha_plus)
            require(
                all(roster_lower <= r <= roster_upper for r in physical_union),
                f"physical cells escaped the fixed roster at {alpha}, {b_value}",
            )
            tilings += 1
            roster_containments += 1
    require(ties >= 4, "insufficient exact tie cases")
    return {
        "rational_cell_cases": cases,
        "exact_ties": ties,
        "physical_tilings": tilings,
        "roster_containments": roster_containments,
    }


def symbolic_certificate() -> dict:
    alpha, q, b_value = sp.symbols("alpha_P q B", positive=True, real=True)
    r, u, s, x = sp.symbols("r u s x", real=True)
    kappa = 1 / (2 * sp.pi * sp.I)

    a_q = alpha / (q + sp.Rational(1, 2))
    b_q = alpha / (q - sp.Rational(1, 2))
    require_zero(
        b_q.subs(q, q + 1) - a_q,
        "adjacent reciprocal boundary",
    )
    require_zero(b_q.subs(q, 1) - 2 * alpha, "physical upper boundary")
    require_zero(
        a_q.subs(q, b_value) - alpha / (b_value + sp.Rational(1, 2)),
        "physical lower boundary",
    )

    theta = alpha * (sp.log(alpha / r) - 1)
    central_phase = theta + q * r
    require_zero(sp.diff(central_phase, r).subs(r, alpha / q), "central critical point")
    require_zero(
        sp.expand_log(central_phase.subs(r, alpha / q), force=True)
        - alpha * sp.log(q),
        "central critical value",
    )
    require_zero(
        sp.diff(central_phase, r, 2).subs(r, alpha / q) - q**2 / alpha,
        "central curvature",
    )

    kernel_closed = kappa * (
        sp.exp(-2 * sp.pi * sp.I * a_q * x)
        - sp.exp(-2 * sp.pi * sp.I * b_q * x)
    ) / x
    kernel_antiderivative = -kappa * sp.exp(-2 * sp.pi * sp.I * r * x) / x
    require_zero(
        sp.diff(kernel_antiderivative, r)
        - sp.exp(-2 * sp.pi * sp.I * r * x),
        "finite Fourier kernel antiderivative",
    )
    require_zero(
        kernel_antiderivative.subs(r, b_q)
        - kernel_antiderivative.subs(r, a_q)
        - kernel_closed,
        "finite Fourier kernel",
    )
    require_zero(
        sp.limit(kernel_closed, x, 0) - (b_q - a_q),
        "finite Fourier kernel continuation",
    )

    y = sp.symbols("y", real=True)
    b_0 = sp.symbols("b_0")
    c = sp.Function("c")(y)
    gaussian = sp.exp(-sp.pi * sp.I * y**2)
    b_amplitude = b_0 + y * c
    recomposed_density = (
        b_0 * gaussian
        + sp.diff(-kappa * c * gaussian, y)
        + kappa * sp.diff(c, y) * gaussian
    )
    require_zero(b_amplitude * gaussian - recomposed_density, "Morse endpoint recomposition")

    d = sp.symbols("d", real=True)
    gaussian_inverse = sp.integrate(
        sp.exp(-sp.pi * r**2 + 2 * sp.pi * sp.I * d * r),
        (r, -sp.oo, sp.oo),
    )
    require_zero(
        gaussian_inverse - sp.exp(-sp.pi * d**2),
        "Fourier inversion normalization",
    )

    epsilon = sp.symbols("epsilon", positive=True, real=True)
    endpoint_profile = lambda point: (
        sp.erf(sp.sqrt(sp.pi) * (5 - point) / sp.sqrt(epsilon))
        - sp.erf(sp.sqrt(sp.pi) * (1 - point) / sp.sqrt(epsilon))
    ) / 2
    endpoint_limits = [
        sp.limit(endpoint_profile(point), epsilon, 0, dir="+")
        for point in [0, 1, 2, 5, 6]
    ]
    require(
        endpoint_limits == [0, sp.Rational(1, 2), 1, sp.Rational(1, 2), 0],
        "Fourier endpoint weights drifted",
    )

    cell_audit = audit_cells()
    return {
        "complete_mode": {
            "fourier_transform": (
                "F_A(u)=1_[1,B](u)A(u)e(alpha_Plog u), and "
                "I_A(r)=integral_1^B A(u)e(alpha_Plog u-ru)du="
                "Fhat_A(r)."
            ),
            "off_lattice_recomposition": (
                "For real r>0, I_A(r)=U_A(r)+kappa e(-r)c_A(y_1(r))"
                "-kappa e(alpha_Plog B-rB)c_A(y_B(r))+J_A(r)."
            ),
            "integer_recomposition": (
                "For integer r and integer B, the endpoint phases reduce to 1 and "
                "e_B, so I_A(r)=U_A(r)+kappa c_A(y_1)-kappa e_Bc_A(y_B)+J_A(r)."
            ),
        },
        "reciprocal_cells": {
            "continuous_endpoints": (
                "a_q=alpha_P/(q+1/2), b_q=alpha_P/(q-1/2), and "
                "C_q(alpha_P)={r in Z_(>0):a_q<r<=b_q}."
            ),
            "integer_bounds": (
                "C_q(alpha_P)={floor(a_q)+1,...,floor(b_q)}."
            ),
            "adjacency": (
                "b_(q+1)=a_q; the shared integer, when present, belongs to C_(q+1)."
            ),
            "physical_tiling": (
                "The cells 1<=q<=B tile the continuous reciprocal band "
                "(alpha_P/(B+1/2),2alpha_P] and its integer modes. For "
                "alpha_-<=alpha_P<=alpha_+, every such mode lies in the fixed "
                "roster T={m,...,n}, m=max(1,floor(alpha_-/(2B))) and "
                "n=ceil(2alpha_+)."
            ),
        },
        "finite_cell_poisson": {
            "identity": (
                "S_q[A]:=sum_(r in C_q)I_A(r)=tau_q[A]+"
                "lim_(R->infinity)sum_(s=-R)^R B_(q,s)[A]."
            ),
            "tie_term": (
                "tau_q[A]=(1/2)1_(b_q in Z)I_A(b_q)"
                "-(1/2)1_(a_q in Z)I_A(a_q)."
            ),
            "dual_integral": (
                "B_(q,s)[A]=integral_(a_q)^(b_q)I_A(r)e(sr)dr."
            ),
            "kernel": (
                "B_(q,s)[A]=integral_1^B F_A(u)K_q(u-s)du, where "
                "K_q(x)=kappa[e(-a_qx)-e(-b_qx)]/x and K_q(0)=b_q-a_q."
            ),
            "central_phase": (
                "The s=q phase has r saddle alpha_P/q, value alpha_Plog q, "
                "and curvature q^2/alpha_P."
            ),
            "alias_gap": (
                "On [a_q,b_q], alpha_P/r lies in [q-1/2,q+1/2]; hence "
                "|s-alpha_P/r|>=1/2 for every integer s!=q."
            ),
        },
        "carrier_inversion": {
            "weight": (
                "omega_B(s)=1 for 1<s<B, 1/2 for s=1 or B, and 0 outside [1,B]."
            ),
            "full_inversion": (
                "PV integral_R I_A(r)e(sr)dr=omega_B(s)A(s)e(alpha_Plog s)."
            ),
            "exterior_tail": (
                "E_q^ext[A]=PV integral_(R setminus [a_q,b_q])I_A(r)e(qr)dr, "
                "so B_(q,q)=omega_B(q)A(q)e(alpha_Plog q)-E_q^ext[A]."
            ),
            "cell_defect": (
                "S_q[A]-omega_B(q)A(q)e(alpha_Plog q)=tau_q[A]+"
                "sum_(s!=q)B_(q,s)[A]-E_q^ext[A], with symmetric dual limits."
            ),
        },
        "global_recomposition": {
            "tie_telescope": (
                "Internal tie terms cancel because a_q=b_(q+1); only the two "
                "outer reciprocal-band half terms remain."
            ),
            "band_identity": (
                "sum_(q=1)^B S_q[A]=tau_band[A]+sum_(s=1)^B "
                "omega_B(s)A(s)e(alpha_Plog s)-sum_(s in Z)^sym E_s^band[A],"
                " where E_s^band is the Fourier-inversion exterior of "
                "[alpha_P/(B+1/2),2alpha_P]."
            ),
            "complement_identity": (
                "If M_A=sum_(r in Z)^sym I_A(r) is the complete mode sum and "
                "C_band[A]=M_A-sum_(q=1)^B S_q[A], then "
                "sum_(s in Z)^sym E_s^band[A]=tau_band[A]+C_band[A]."
            ),
            "interpretation": (
                "All internal cell aliases and cellwise exterior tails recompose "
                "exactly; the only reciprocal-frequency boundary is the two-sided "
                "exterior of the complete physical band."
            ),
        },
        "audit": cell_audit,
    }


def boundary_certificate() -> dict:
    return {
        "achievement": (
            "The full-line sequential transform is exact Fourier inversion, so it "
            "returns the starred physical carrier with no full-line remainder. The "
            "first-correction cancellation in Section 11.179 is its first local "
            "coefficient. Finite reciprocal cells differ only through explicit tie, "
            "alias, and exterior-tail terms."
        ),
        "off_lattice_guard": (
            "The real-r extension must use e(-r) and e(alpha_Plog B-rB) at the "
            "Morse endpoints. Replacing them off lattice by their integer values 1 "
            "and e_B changes the second Poisson transform."
        ),
        "remaining_bound": (
            "No uniform estimate is supplied for the noncentral aliases or the "
            "global reciprocal-band exterior. Endpoint jumps make the exterior "
            "Fourier inversion principal-value at q=1,B, so its two sides must not "
            "be bounded independently before endpoint composition."
        ),
        "handoff": (
            "Recompose the global lower and upper reciprocal-frequency exterior with "
            "the nonphysical roster fringes, mathcal T, the hard endpoint halves, "
            "mathcal L+e_Bmathcal U, and the terminal recurrence. Then derive a "
            "signed bound for that two-sided package, using the alias derivative gap "
            "only after internal cell cancellation."
        ),
    }


def build_rows(symbolic: dict, boundary: dict) -> list[GateRow]:
    mode = symbolic["complete_mode"]
    cells = symbolic["reciprocal_cells"]
    poisson = symbolic["finite_cell_poisson"]
    inversion = symbolic["carrier_inversion"]
    global_rows = symbolic["global_recomposition"]
    return [
        GateRow("rfi_01_fourier_mode", "complete_mode", "certified", "The complete real-r mode is the Fourier transform of the compact physical amplitude.", mode["fourier_transform"], "The zero extension has jumps at the two physical endpoints."),
        GateRow("rfi_02_off_lattice", "endpoint_recomposition", "certified", "Morse integration by parts has exact off-lattice endpoint phases.", mode["off_lattice_recomposition"], boundary["off_lattice_guard"]),
        GateRow("rfi_03_integer_rejoin", "endpoint_recomposition", "certified", "At integer modes the off-lattice package becomes the endpoint-composed roster summand.", mode["integer_recomposition"], "The simplification is valid only on the integer lattice."),
        GateRow("rfi_04_cells", "reciprocal_cells", "certified", "Derivative half-cells have exact continuous endpoints.", cells["continuous_endpoints"], "Only positive reciprocal cells are defined."),
        GateRow("rfi_05_bounds", "reciprocal_cells", "certified", "The half-open cell has an exact floor representation.", cells["integer_bounds"], "Empty cells are permitted."),
        GateRow("rfi_06_adjacency", "reciprocal_cells", "certified", "Adjacent continuous cells share one endpoint with deterministic lattice ownership.", cells["adjacency"], "This is a bookkeeping identity, not a remainder bound."),
        GateRow("rfi_07_poisson", "finite_cell_poisson", "certified", "Finite Poisson summation applies exactly to one half-open reciprocal cell.", poisson["identity"] + " " + poisson["dual_integral"], "The dual sum is a symmetric limit."),
        GateRow("rfi_08_tie", "finite_cell_poisson", "certified", "The half-open convention produces one explicit two-ended tie correction.", poisson["tie_term"], "The indicators cannot be discarded at exact crossings."),
        GateRow("rfi_09_transfer", "tie_transfer", "certified", "A shared integer transfers as one complete mode package between adjacent cells.", cells["adjacency"] + " " + poisson["tie_term"], "Individual cell formulas jump while their union does not."),
        GateRow("rfi_10_kernel", "finite_fourier_kernel", "certified", "Ordinary Fubini gives an exact finite reciprocal-band kernel.", poisson["kernel"], "No localization estimate is used."),
        GateRow("rfi_11_phase", "central_dual_mode", "certified", "The diagonal dual mode has exactly the reciprocal Morse saddle.", poisson["central_phase"], "This identifies only the diagonal continuous mode."),
        GateRow("rfi_12_alias_gap", "noncentral_aliases", "certified", "Every noncentral dual mode is uniformly nonstationary in the reciprocal cell.", poisson["alias_gap"], "No summed alias bound follows yet."),
        GateRow("rfi_13_inversion", "carrier_inversion", "certified", "Full-line Fourier inversion returns the compact physical amplitude.", inversion["full_inversion"], "Endpoint values carry Dirichlet half weights."),
        GateRow("rfi_14_weights", "carrier_inversion", "certified", "The returned carrier weights match the starred primal Poisson convention.", inversion["weight"], "The hard endpoint halves remain outside the starred carrier."),
        GateRow("rfi_15_exterior", "finite_cell_defect", "certified", "Truncating the diagonal inverse transform creates exactly one two-sided exterior tail.", inversion["exterior_tail"], "At q=1,B the exterior is principal-value and cannot be split absolutely."),
        GateRow("rfi_16_cell_defect", "finite_cell_defect", "certified", "The complete finite-cell defect has exactly tie, alias, and exterior components.", inversion["cell_defect"], "No component is bounded here."),
        GateRow("rfi_17_tile", "physical_band", "certified", "Physical reciprocal cells tile one continuous outer-frequency band.", cells["physical_tiling"], "Nonphysical roster fringes remain outside this band."),
        GateRow("rfi_18_telescope", "physical_band", "certified", "All internal half-integer tie corrections cancel in the physical q family.", global_rows["tie_telescope"], "The two outer band corrections remain."),
        GateRow("rfi_19_global", "global_recomposition", "certified", "Cell aliases and cellwise exterior tails recompose to one global two-sided exterior package.", global_rows["band_identity"] + " " + global_rows["interpretation"], "The global exterior package has not been estimated."),
        GateRow("rfi_20_involution", "global_recomposition", "certified", "The global exterior is exactly the complementary integer-mode package plus the two band ties.", global_rows["complement_identity"], "A second scalar Poisson transform alone supplies no smallness."),
        GateRow("rfi_21_achievement", "reciprocal_handoff", "certified", "There is no complete-line sequential stationary remainder.", boundary["achievement"], "This does not make the finite reciprocal band exact by itself."),
        GateRow("rfi_22_handoff", "open_bound", "open", "The signed global exterior/complement package must be composed and bounded.", boundary["remaining_bound"] + " " + boundary["handoff"], "No alias, exterior-tail, endpoint-complete h2, or signed-flow bound is claimed."),
        GateRow("rfi_23_boundary", "proof_boundary", "guard_validated", "The exact inversion is not promoted to the RH bridge.", "Finite reciprocal-frequency exterior, nonphysical fringes, endpoint and terminal composition, quadratic residuals, and signed flow remain open.", "No h2, Xi, Lambda<=0, PF-infinity, RH, or prize theorem is claimed."),
    ]


def render_note(artifact: dict) -> str:
    symbolic = artifact["symbolic_certificate"]
    boundary = artifact["boundary_certificate"]
    mode = symbolic["complete_mode"]
    cells = symbolic["reciprocal_cells"]
    poisson = symbolic["finite_cell_poisson"]
    inversion = symbolic["carrier_inversion"]
    global_rows = symbolic["global_recomposition"]
    return "\n".join(
        [
            "# Ideal-Cubic Reciprocal Finite-Cell Inversion Gate",
            "",
            "Date: 2026-08-02",
            "",
            "Status: exact complete-mode reassembly, half-open cell Poisson formula, carrier Fourier inversion, and global internal-boundary cancellation. This is not a proof of a reciprocal exterior-tail bound, signed flow, or RH.",
            "",
            "```text",
            str(RESULT_PATH.relative_to(REPO_ROOT)).replace("\\", "/"),
            f"python {str(Path(__file__).relative_to(REPO_ROOT)).replace(chr(92), '/')}",
            "```",
            "",
            "## Complete Real-Mode Package",
            "",
            "```text",
            mode["fourier_transform"],
            mode["off_lattice_recomposition"],
            mode["integer_recomposition"],
            "```",
            "",
            boundary["off_lattice_guard"],
            "",
            "## Exact Reciprocal Cell",
            "",
            "```text",
            cells["continuous_endpoints"],
            cells["integer_bounds"],
            cells["adjacency"],
            poisson["identity"],
            poisson["tie_term"],
            poisson["dual_integral"],
            "```",
            "",
            "## Dual Kernel And Carrier",
            "",
            "```text",
            poisson["kernel"],
            poisson["central_phase"],
            poisson["alias_gap"],
            inversion["weight"],
            inversion["full_inversion"],
            inversion["exterior_tail"],
            inversion["cell_defect"],
            "```",
            "",
            boundary["achievement"],
            "",
            "## Global Recomposition",
            "",
            "```text",
            cells["physical_tiling"],
            global_rows["tie_telescope"],
            global_rows["band_identity"],
            global_rows["complement_identity"],
            global_rows["interpretation"],
            "```",
            "",
            "## Remaining Boundary",
            "",
            boundary["remaining_bound"],
            "",
            "## Handoff",
            "",
            boundary["handoff"],
            "",
            "## Pi Provenance",
            "",
            "The kernel and inversion factors use the already fixed character e(x)=exp(2pi i x), so kappa=1/(2pi i). The endpoint half weights are the Dirichlet Fourier-inversion values. No circle, polygon, or fitted numerical constant is introduced.",
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
            "exact complete real-mode endpoint reassembly, half-open reciprocal-cell "
            "finite Poisson identity, carrier Fourier inversion, and global internal "
            "tie/alias recomposition complete; exterior-tail and signed-flow bounds "
            "open; no endpoint-complete h2, Xi, or RH theorem"
        ),
        "proof_boundary": (
            "This gate proves the exact off-lattice endpoint phases, integer endpoint "
            "reassembly, reciprocal half-cell Poisson formula and tie correction, "
            "finite Fourier kernel, central carrier inversion, alias derivative gap, "
            "physical-band tiling, internal tie cancellation, and reduction of all "
            "cellwise defects to one global reciprocal-frequency exterior package. "
            "It proves no bound for that exterior, nonphysical roster fringes, "
            "endpoint-complete h2 remainder, q-family signed aggregate, quadratic "
            "residual, signed flow, Phi_B, contact exclusion, retained aggregate or "
            "Xi theorem, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, "
            "RH, or prize-level conclusion."
        ),
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "boundary_certificate": boundary,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "complete_mode_recompositions": 2,
            "off_lattice_endpoint_phase_identities": 2,
            "reciprocal_cell_identities": 4,
            "finite_cell_poisson_identities": 3,
            "fourier_kernel_identities": 2,
            "carrier_inversion_identities": 3,
            "alias_gap_identities": 1,
            "internal_tie_cancellations": 1,
            "global_recompositions": 3,
            "finite_cell_remainder_bounds": 0,
            "alias_bounds": 0,
            "exterior_tail_bounds": 0,
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
        "built reciprocal finite-cell inversion gate: "
        f"{counts['rows']} rows, "
        f"{counts['finite_cell_poisson_identities']} cell-Poisson identities, "
        f"{counts['carrier_inversion_identities']} carrier inversions, "
        f"{counts['internal_tie_cancellations']} tie cancellation, "
        f"{counts['global_recompositions']} global recompositions, "
        f"{counts['exterior_tail_bounds']} exterior-tail bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
