#!/usr/bin/env python3
"""Build the endpoint first-pivot odd-small-ball reduction and route guard."""

from __future__ import annotations

import argparse
import cmath
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "endpoint_first_pivot_odd_small_ball_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "endpoint_schur": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "endpoint_schur_cohn_first_jet_guard.json"
    ),
    "absolute_phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "joined_odd_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_dyadic_odd_prefix_first_jet_reduction.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.json"
    ),
}


@dataclass(frozen=True)
class PivotRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {key: file_hash(path) for key, path in SOURCES.items()}


def source_audit() -> dict[str, str]:
    payloads = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }
    markers = {
        "endpoint_schur": (
            "Delta_K=|r_0+B_(0,0)|^2-|B_(K,0)|^2",
            "outside-disk Schur-Cohn recursion",
        ),
        "absolute_phase_anchor": (
            "r_0=g_0/f_1=(-1)^N*beta*(T_0+i)*H_a",
            "q_n=f_n/f_1=",
        ),
        "joined_odd_prefix": (
            "M_k=floor(N/2^k)",
            "P_N(z,xi,mu)=sum_(k=0)^K",
        ),
        "adjacent_recurrence": (
            "Q_N=f_(N+1)+kappa_N*J_a",
            "J_a=H_a(p+2)+H_a(p)",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        key: str(payload.get("kind", ""))
        for key, payload in payloads.items()
    }


def norm_square(value: sp.Expr) -> sp.Expr:
    return sp.expand_complex(value * sp.conjugate(value))


def denominator_audit() -> dict[str, str]:
    o_r, o_i, r_r, r_i, c_r, c_i = sp.symbols(
        "o_r o_i r_r r_i c_r c_i", real=True
    )
    d_r, d_i = sp.symbols("d_r d_i", real=True)
    odd = o_r + sp.I * o_i
    endpoint = r_r + sp.I * r_i
    terminal = c_r + sp.I * c_i
    denominator = d_r + sp.I * d_i
    denominator_norm = norm_square(denominator)
    pivot = (
        norm_square((odd + endpoint) / denominator)
        - norm_square(terminal / denominator)
    )
    expected = (
        norm_square(odd + endpoint) - norm_square(terminal)
    ) / denominator_norm
    if sp.simplify(pivot - expected) != 0:
        raise RuntimeError("common-denominator pivot cancellation failed")

    correlation = (
        norm_square(odd + endpoint)
        - norm_square(odd)
        - norm_square(endpoint)
        - 2 * sp.re(odd * sp.conjugate(endpoint))
    )
    if sp.simplify(correlation) != 0:
        raise RuntimeError("pivot correlation expansion failed")

    return {
        "cancellation": (
            "For nonzero D, |(O+R)/D|^2-|C/D|^2="
            "(|O+R|^2-|C|^2)/|D|^2."
        ),
        "correlation": (
            "|O+R|^2=|O|^2+|R|^2+2Re(O*conj(R))."
        ),
    }


def triangle_audit() -> dict:
    a = sp.sqrt(sp.Rational(1, 3))
    b = sp.sqrt(sp.Rational(1, 5))
    # a+b>1 is equivalent to 30>7sqrt(15); both sides are positive,
    # and squaring gives the exact integer comparison 900>735.
    if not 900 > 735:
        raise RuntimeError("triangle upper-radius audit failed")
    if not sp.Rational(1, 3) > sp.Rational(1, 5):
        raise RuntimeError("triangle side ordering audit failed")
    delta = sp.N(a + b - 1, 40)
    if delta <= 0:
        raise RuntimeError("endpoint-neighborhood radius is not positive")

    omega = 62643 / 100
    witness = (
        1
        + cmath.exp(1j * omega * math.log(3)) / math.sqrt(3)
        + cmath.exp(1j * omega * math.log(5)) / math.sqrt(5)
    )
    if abs(witness) >= 0.001:
        raise RuntimeError("linked-phase numerical witness drifted")
    return {
        "a": "1/sqrt(3)",
        "b": "1/sqrt(5)",
        "upper_radius_proof": "a+b>1 because 30>7sqrt(15), as 900>735.",
        "lower_radius_proof": "|a-b|<1 because 0<b<a<1.",
        "endpoint_radius": "delta_0=1/sqrt(3)+1/sqrt(5)-1",
        "endpoint_radius_decimal": str(delta),
        "diagnostic_omega": "62643/100",
        "diagnostic_modulus": f"{abs(witness):.17e}",
    }


def build_exact() -> dict:
    algebra = denominator_audit()
    triangle = triangle_audit()
    return {
        "pi_provenance": (
            "The pi in beta and T_0 comes from the completed-zeta "
            "normalization and Riemann-Siegel saddle. The 2*pi used for "
            "phase reduction is exactly the period of exp(i*theta). The "
            "dyadic split, common-denominator cancellation, Schur pivot, "
            "and triangle geometry introduce no new pi."
        ),
        "domain": (
            "In a fixed q>=1 physical chart let h=log(2), "
            "K=floor(log_2 N), s_*=sigma_*-i*omega_*, and |d_1|<1/2."
        ),
        "terminal_singleton": (
            "Since 2^K<=N<2^(K+1), M_K=floor(N/2^K)=1. Hence the "
            "top dyadic layer contains only m=1 and "
            "B_(K,0)=C_K/(1+d_1), where "
            "C_K=exp[t*(K*h)^2/4-sigma_*K*h]*(1+d_(2^K))."
        ),
        "odd_constant_numerator": (
            "Define O_N=sum_(m<=N,m odd)"
            "exp[t*log(m)^2/4-s_*log(m)]*(1+d_m). "
            "Then B_(0,0)=O_N/(1+d_1)."
        ),
        "endpoint_numerator": (
            "Define R_N=(-1)^N*beta*(T_0+i)*H_a/M_t(s). "
            "The branch-free phase-anchor formula gives "
            "r_0=R_N/(1+d_1), with beta>0 and H_a real."
        ),
        "common_denominator": (
            "Therefore b_0=r_0+B_(0,0)=(O_N+R_N)/(1+d_1) and "
            "b_K=B_(K,0)=C_K/(1+d_1). " + algebra["cancellation"]
        ),
        "pivot_numerator": (
            "Delta_K={|O_N+R_N|^2-rho_K^2}/|1+d_1|^2, where "
            "rho_K=exp[t*(K*h)^2/4-sigma_*K*h]*|1+d_(2^K)|. "
            "Thus the explicit normalization denominator affects neither "
            "the sign nor the zero set of the first pivot. The m=1 "
            "occurrence of d_1 remains inside O_N."
        ),
        "small_ball_equivalence": (
            "Delta_K>0 iff |O_N+R_N|>rho_K, equivalently "
            "O_N is outside the closed disk Dbar(-R_N,rho_K). "
            "The first Schur obligation is an endpoint-shifted odd "
            "Dirichlet-polynomial small-ball exclusion."
        ),
        "correlation_form": (
            algebra["correlation"]
            + " Hence Delta_K>0 iff "
            "|O_N|^2+|R_N|^2+2Re(O_N*conj(R_N))>rho_K^2."
        ),
        "endpoint_phase": (
            "When H_a!=0, arg(R_N)=N*pi+arg(H_a)+arg(T_0+i)"
            "-arg(M_t(s)) modulo 2*pi; when H_a=0, R_N=0. This is "
            "branch-free phase bookkeeping, not a lower bound or an "
            "alignment theorem for O_N."
        ),
        "sufficient_routes": (
            "Any one of ||O_N|-|R_N||>rho_K or "
            "|Re(u*(O_N+R_N))|>rho_K for a proved unit complex u is "
            "sufficient. Equivalently one may prove the displayed "
            "correlation inequality. No certified current source supplies "
            "one of these uniformly."
        ),
        "linked_phase_model": (
            "Set t=0, sigma_*=1/2, d_n=0, R_N=0, and N=5. Then K=2, "
            "rho_2=1/2, and "
            "O_5(omega)=1+3^(-1/2)exp(i*omega*log(3))"
            "+5^(-1/2)exp(i*omega*log(5))."
        ),
        "density_theorem": (
            "The continuous flow omega -> "
            "(omega*log(3),omega*log(5)) mod 2*pi is dense in the "
            "two-torus. By Kronecker's theorem it is enough that no "
            "nonzero integers u,v satisfy u*log(3)+v*log(5)=0; "
            "exponentiation would give 3^u*5^v=1, impossible by unique "
            "prime factorization."
        ),
        "triangle_geometry": (
            "The two rotating vectors have lengths a=1/sqrt(3) and "
            "b=1/sqrt(5), and their sums fill the closed annulus "
            "|a-b|<=|w|<=a+b. The exact inequalities "
            "|a-b|<1<a+b hold. "
            + triangle["upper_radius_proof"]
        ),
        "linked_pivot_failure": (
            "The annulus contains -1, so torus density gives "
            "inf_omega|O_5(omega)|=0. Therefore some linked logarithmic "
            "phases satisfy |O_5(omega)|<1/2=rho_2 and Delta_2<0. "
            "The diagnostic omega=62643/100 gives modulus "
            + triangle["diagnostic_modulus"]
            + "; the exact conclusion uses density, not this decimal."
        ),
        "endpoint_neighborhood": (
            "More generally, for fixed endpoint e, "
            "O_(5,e)=1+e+a*exp(i*omega*log(3))"
            "+b*exp(i*omega*log(5)) has infimum zero whenever "
            "|a-b|<|1+e|<a+b. In particular this holds for every "
            "|e|<delta_0=a+b-1="
            + triangle["endpoint_radius_decimal"]
            + ". Thus endpoint smallness alone does not restore the pivot."
        ),
        "scope_guard": (
            "The N=5 model preserves the actual odd logarithmic phase "
            "linkage and correction-free Dirichlet amplitudes, but omega "
            "is allowed to range independently of the physical cutoff and "
            "the recurrent Xi endpoint is replaced by zero or fixed e. "
            "It is an exact route countermodel, not a counterexample to an "
            "actual Xi pivot, Newman, or RH."
        ),
        "adjacent_boundary": (
            "The exact cutoff update Q_N=f_(N+1)+kappa_N*J_a controls the "
            "physical value jump. Existing adjacent certificates bound "
            "selected real projections, not the complex displacement of "
            "O_N relative to -R_N, so they do not imply disk avoidance."
        ),
        "route_decision": (
            "Retain the first pivot as a precise falsification coordinate, "
            "but stop treating decreasing amplitudes, logarithmic phase "
            "linkage, or endpoint smallness as candidate proofs. Full "
            "unit-disk stability now requires an Xi-specific correlation "
            "or small-ball theorem coupling O_N to R_N. Unless the source "
            "formulas reveal such a coupling, prioritize actual linked-"
            "point nonvanishing or the signed three-cylinder degree route, "
            "which asks less than control of every free z on |z|=1."
        ),
        "live_q_ge_1_target": (
            "Either prove uniformly on every q>=1 fixed-N Xi chart that "
            "O_N avoids Dbar(-R_N,rho_K), using the actual relation among "
            "omega_*, N, d_n, H_a, and M_t(s), or produce a certified "
            "physical chart where the first pivot fails and retire full "
            "disk stability."
        ),
        "weaker_primary_target": (
            "For the weaker primary route, prove nonvanishing only at "
            "z_*=exp(i*omega_*log(2)), or prove the endpoint-complete "
            "signed crossing budget through the three exact transport "
            "cylinders while retaining H_2,D_0,D_1 and the recurrent "
            "endpoint first jet."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 multiplicity-compatible parabolic/Hermite chart "
            "remains separate. This first-pivot analysis supplies no "
            "small-q closure, finite connector, or chart join."
        ),
        "proof_boundary": (
            "The terminal singleton, odd and endpoint numerators, common-"
            "denominator cancellation, pivot disk equivalence, correlation "
            "form, endpoint phase bookkeeping, linked-flow density guard, "
            "triangle-annulus pivot failure, and endpoint-neighborhood "
            "extension are exact. No uniform actual Xi small-ball "
            "exclusion, physical pivot failure, later Schur pivot, full "
            "endpoint disk stability, linked-point lower bound, signed "
            "crossing theorem, strict successor flux bound, q<1 closure, "
            "contact exclusion, Lambda<=0, PF-infinity, RH proof, or "
            "Clay-prize conclusion is asserted."
        ),
        "diagnostics": triangle,
    }


def build_rows(exact: dict) -> list[PivotRow]:
    return [
        PivotRow(
            "efposbg_00_pi_provenance",
            "definition_provenance",
            "certified",
            "Every occurrence of pi has an explicit inherited source.",
            exact["pi_provenance"],
            "The small-ball reduction itself introduces no new pi.",
        ),
        PivotRow(
            "efposbg_01_terminal_singleton",
            "exact_reindexing",
            "certified",
            "The terminal dyadic layer is the singleton m=1.",
            exact["terminal_singleton"],
            "This holds for every integer N>=1.",
        ),
        PivotRow(
            "efposbg_02_odd_numerator",
            "exact_reindexing",
            "ready_to_apply",
            "The constant odd layer has one explicit unnormalized numerator.",
            exact["odd_constant_numerator"],
            "All actual d_m corrections and odd phases are retained.",
        ),
        PivotRow(
            "efposbg_03_endpoint_numerator",
            "exact_reduction",
            "ready_to_apply",
            "The recurrent endpoint has the same first-coefficient denominator.",
            exact["endpoint_numerator"],
            "The formula is branch-free but supplies no phase alignment.",
        ),
        PivotRow(
            "efposbg_04_denominator_cancellation",
            "exact_identity",
            "certified",
            "The common correction denominator cancels from the first pivot.",
            exact["common_denominator"] + " " + exact["pivot_numerator"],
            "The common denominator cancels; all numerator corrections remain.",
        ),
        PivotRow(
            "efposbg_05_small_ball",
            "exact_equivalence",
            "certified",
            "The first pivot is exactly an endpoint-shifted odd small-ball exclusion.",
            exact["small_ball_equivalence"],
            "This reformulates the open Xi inequality; it does not prove it.",
        ),
        PivotRow(
            "efposbg_06_correlation",
            "exact_equivalence",
            "ready_to_apply",
            "The missing phase information is one explicit endpoint correlation.",
            exact["correlation_form"],
            "A lower bound for the correlation remains Xi-specific.",
        ),
        PivotRow(
            "efposbg_07_endpoint_phase",
            "exact_phase_bookkeeping",
            "ready_to_apply",
            "The endpoint direction is explicit without choosing a phase branch.",
            exact["endpoint_phase"],
            "The zero-endpoint case must be handled separately.",
        ),
        PivotRow(
            "efposbg_08_sufficient_routes",
            "conditional_route",
            "not_ready_to_apply",
            "Absolute or directional separation would prove the first pivot.",
            exact["sufficient_routes"],
            "No uniform source theorem currently supplies either separation.",
        ),
        PivotRow(
            "efposbg_09_linked_model",
            "exact_countermodel_setup",
            "guard_validated",
            "A three-odd-term model retains linked logarithmic phases.",
            exact["linked_phase_model"],
            "Its free omega is not tied to a physical Xi cutoff chart.",
        ),
        PivotRow(
            "efposbg_10_density",
            "exact_number_theory",
            "certified",
            "The linked logarithmic phase flow is dense in the two-torus.",
            exact["density_theorem"],
            "This is the standard continuous Kronecker theorem.",
        ),
        PivotRow(
            "efposbg_11_triangle",
            "exact_geometry",
            "certified",
            "The two odd vectors can cancel the leading unit vector.",
            exact["triangle_geometry"],
            "The annulus statement uses both independent torus phases.",
            exact["diagnostics"],
        ),
        PivotRow(
            "efposbg_12_linked_pivot_failure",
            "countermodel",
            "guard_validated",
            "Linked logarithmic phases and decreasing amplitudes do not force the first pivot.",
            exact["linked_pivot_failure"],
            "This is a route countermodel, not an actual Xi pivot failure.",
            exact["diagnostics"],
        ),
        PivotRow(
            "efposbg_13_endpoint_neighborhood",
            "countermodel_family",
            "guard_validated",
            "A whole neighborhood of small fixed endpoints still permits pivot failure.",
            exact["endpoint_neighborhood"],
            "The physical recurrent endpoint is neither fixed nor arbitrary.",
            exact["diagnostics"],
        ),
        PivotRow(
            "efposbg_14_scope_guard",
            "nonpromotion_guard",
            "guard_validated",
            "The finite model is kept outside the actual Xi theorem boundary.",
            exact["scope_guard"],
            "No Newman or RH conclusion follows from the countermodel.",
        ),
        PivotRow(
            "efposbg_15_adjacent_boundary",
            "route_guard",
            "guard_validated",
            "Real adjacent-chart bounds do not imply complex disk avoidance.",
            exact["adjacent_boundary"],
            "A complex endpoint-relative displacement theorem would be new input.",
        ),
        PivotRow(
            "efposbg_16_route_decision",
            "route_decision",
            "guard_validated",
            "Full disk stability is downgraded unless an Xi-specific correlation appears.",
            exact["route_decision"],
            "The Schur coordinate remains useful for exact falsification.",
        ),
        PivotRow(
            "efposbg_17_open_targets",
            "open_theorem_target",
            "not_ready_to_apply",
            "The live alternatives are physical odd-disk avoidance or linked-point signed transport.",
            exact["live_q_ge_1_target"]
            + " "
            + exact["weaker_primary_target"]
            + " "
            + exact["q_lt_1_target"],
            "Both q>=1 alternatives and the q<1 chart remain open.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-27",
        "status": (
            "exact common-denominator first-pivot reduction to an endpoint-"
            "shifted odd small-ball problem, with a linked-logarithmic-phase "
            "countermodel and small-endpoint extension; the actual Xi "
            "small-ball or weaker signed-transport theorem remains open"
        ),
        "proof_boundary": exact["proof_boundary"],
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman Endpoint First-Pivot Odd Small-Ball Guard",
            "",
            "Date: 2026-07-27",
            "",
            "Status: exact first-pivot reduction and linked-phase route",
            "countermodel; not a physical Xi pivot counterexample and not a",
            "proof of `Lambda<=0`, PF-infinity, RH, or a Clay-prize result.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Pi Provenance",
            "",
            exact["pi_provenance"],
            "",
            "## Terminal Singleton",
            "",
            "```text",
            exact["terminal_singleton"],
            "```",
            "",
            "## Common Numerators",
            "",
            "```text",
            exact["odd_constant_numerator"],
            exact["endpoint_numerator"],
            exact["common_denominator"],
            "```",
            "",
            "## Exact First Pivot",
            "",
            "```text",
            exact["pivot_numerator"],
            exact["small_ball_equivalence"],
            exact["correlation_form"],
            "```",
            "",
            "## Endpoint Phase",
            "",
            "```text",
            exact["endpoint_phase"],
            exact["sufficient_routes"],
            "```",
            "",
            "## Linked-Phase Countermodel",
            "",
            "```text",
            exact["linked_phase_model"],
            exact["density_theorem"],
            exact["triangle_geometry"],
            exact["linked_pivot_failure"],
            "```",
            "",
            "## Small-Endpoint Extension",
            "",
            "```text",
            exact["endpoint_neighborhood"],
            "```",
            "",
            exact["scope_guard"],
            "",
            "## Adjacent-Chart Boundary",
            "",
            exact["adjacent_boundary"],
            "",
            "## Route Decision",
            "",
            exact["route_decision"],
            "",
            "## Live Theorems",
            "",
            "```text",
            exact["live_q_ge_1_target"],
            exact["weaker_primary_target"],
            exact["q_lt_1_target"],
            "```",
            "",
            "## Boundary",
            "",
            exact["proof_boundary"],
            "",
        ]
    )


def write_outputs(artifact: dict, output: Path, note: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_outputs(artifact, args.output, args.note)
    print(
        "built Newman endpoint first-pivot odd small-ball guard: "
        "18 rows, 1 exact denominator-cancelled pivot, "
        "1 linked-logarithmic-phase countermodel, "
        "1 small-endpoint countermodel family, 2 open Xi routes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
