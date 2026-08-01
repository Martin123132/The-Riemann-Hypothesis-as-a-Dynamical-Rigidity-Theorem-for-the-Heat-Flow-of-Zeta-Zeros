#!/usr/bin/env python3
"""Build the interior-carrier projective-current theorem and edge guard."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "interior_projective_current_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "contact_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contact_signed_transport_reduction.json"
    ),
    "absolute_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_block_composition_guard.json"
    ),
    "real_residual": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_residual_reduction.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
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
        "contact_transport": (
            "sum_(j=0)^N d_j=mathcal_C_N",
            "generally complex",
            "cumulative-real-mass",
        ),
        "absolute_phase": (
            "nu_n=partial_x arg(z_n)",
            "nu_n<=8449/x^2-u_n/2",
            "only the terminal n=N carrier",
        ),
        "real_residual": (
            "s_*'=t*D_x/4+i*(-1/2-t*C_x/4)",
            "|d_(n,x)|<4223/x^2",
        ),
        "adjacent_recurrence": (
            "Q_N=f_(N+1)+kappa_N*J_a",
            "Delta A_a=Re",
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


def carrier_current_audit() -> dict[str, str]:
    (
        x_part,
        y_part,
        radial,
        angular,
        c_s,
        b,
        k,
        u,
        c_s_x,
        b_x,
        u_x,
    ) = sp.symbols(
        "X Y varrho nu c_s b k u c_s_x b_x u_x",
        real=True,
    )
    x_part_x = radial * x_part - angular * y_part
    y_part_x = angular * x_part + radial * y_part
    c_atom = x_part
    d_atom = c_s * k * x_part - b * u * y_part
    d_atom_x = (
        c_s_x * k * x_part
        + c_s * k * x_part_x
        - (b_x * u + b * u_x) * y_part
        - b * u * y_part_x
    )
    determinant = sp.expand(c_atom * d_atom_x - d_atom * x_part_x)
    expected = sp.expand(
        c_s_x * k * x_part**2
        - (b_x * u + b * u_x) * x_part * y_part
        - b * u * angular * (x_part**2 + y_part**2)
    )
    if sp.simplify(determinant - expected) != 0:
        raise RuntimeError("division-free carrier current failed")
    if radial in determinant.free_symbols:
        raise RuntimeError("radial current did not cancel")

    return {
        "carrier_jet": (
            "For z_n=X_n+iY_n and "
            "z_(n,x)=(varrho_n+i*nu_n)z_n, "
            "X_(n,x)=varrho_nX_n-nu_nY_n and "
            "Y_(n,x)=nu_nX_n+varrho_nY_n."
        ),
        "components": (
            "Put c_s=Re(s_*'), b=Im(s_*'), "
            "k_n=u_n-u_N=log(N/n), so k_(n,x)=0. Then "
            "c_n=X_n and d_n=c_s*k_n*X_n-b*u_n*Y_n."
        ),
        "determinant": (
            "J_n:=c_n*d_(n,x)-d_n*c_(n,x)="
            "c_(s,x)*k_n*X_n^2"
            "-(b_x*u_n+b*u_x)*X_n*Y_n"
            "-b*u_n*nu_n*(X_n^2+Y_n^2)."
        ),
        "cancellation": (
            "The radial amplitude current varrho_n cancels identically "
            "from J_n."
        ),
        "ratio": (
            "When X_n!=0, h_n=d_n/c_n="
            "c_s*log(N/n)-b*u_n*tan(theta_n), and "
            "h_(n,x)=c_(s,x)*log(N/n)"
            "-(b_x*u_n+b*u_x)*tan(theta_n)"
            "-b*u_n*nu_n*sec(theta_n)^2=J_n/X_n^2."
        ),
    }


def constant_audit() -> dict[str, str | int]:
    if not 650 < 9 * 2**50:
        raise RuntimeError("curvature envelope boundary check failed")
    if not 67592 < 144 * 2**75:
        raise RuntimeError("phase envelope boundary check failed")
    b_budget = Fraction(1, 64) + Fraction(1, 24)
    if not b_budget < Fraction(1, 16):
        raise RuntimeError("mixed-current budget failed")
    row_one = (
        -Fraction(1, 8)
        + Fraction(1, 64)
        + Fraction(1, 32)
    )
    row_two = -Fraction(1, 8) + Fraction(1, 32)
    if row_one != -Fraction(5, 64):
        raise RuntimeError("first Gershgorin row failed")
    if row_two != -Fraction(3, 32):
        raise RuntimeError("second Gershgorin row failed")
    return {
        "curvature_integer_guard": "650<9*2^50",
        "phase_integer_guard": "67592<144*2^75",
        "mixed_budget_using_pi_gt_3": (
            "1/64+1/(8*pi)<1/64+1/24=11/192<1/16"
        ),
        "gershgorin_row_one": "-1/8+1/64+1/32=-5/64",
        "gershgorin_row_two": "-1/8+1/32=-3/32",
    }


def terminal_audit() -> dict[str, str]:
    x_part, y_part, b, b_x, u_x, angular = sp.symbols(
        "X Y b b_x u_x nu", real=True
    )
    u = sp.symbols("u", real=True)
    terminal_current = sp.expand(
        -(b_x * u + b * u_x) * x_part * y_part
        - b * u * angular * (x_part**2 + y_part**2)
    )
    left_edge = sp.simplify(terminal_current.subs(u, 0))
    expected = -b * u_x * x_part * y_part
    if sp.simplify(left_edge - expected) != 0:
        raise RuntimeError("terminal left-edge current failed")
    positive = left_edge.subs(
        {b: -sp.Rational(1, 2), u_x: 1, x_part: 1, y_part: 1}
    )
    negative = left_edge.subs(
        {b: -sp.Rational(1, 2), u_x: 1, x_part: 1, y_part: -1}
    )
    if positive != sp.Rational(1, 2) or negative != -sp.Rational(1, 2):
        raise RuntimeError("terminal sign guard failed")
    return {
        "formula": (
            "For n=N, k_N=0 and "
            "J_N=-(b_x*u_N+b*u_x)X_NY_N"
            "-b*u_N*nu_N*(X_N^2+Y_N^2)."
        ),
        "left_edge": (
            "At the left cutoff u_N=0, "
            "J_N=-b*u_x*X_NY_N, which has either sign under the "
            "currently proved bounds."
        ),
        "zero_vector": (
            "At u_N=0 and X_N=0, the isolated terminal component "
            "(c_N,d_N)=(X_N,0) is the zero vector, so it has no "
            "projective label."
        ),
        "guard": (
            "The choices b=-1/2,u_x=1,(X_N,Y_N)=(1,1) and (1,-1) "
            "give J_N=1/2 and -1/2. This is a local algebraic "
            "nonpromotion guard, not two asserted Xi cutoff states."
        ),
    }


def order_swap_audit() -> dict[str, str]:
    variable = sp.symbols("x", real=True)
    first = -variable
    second = -variable + sp.sin(variable) / 2
    first_x = sp.diff(first, variable)
    second_x = sp.diff(second, variable)
    if first_x != -1 or second_x != -1 + sp.cos(variable) / 2:
        raise RuntimeError("monotone order-swap derivatives failed")
    if sp.simplify(second - first - sp.sin(variable) / 2) != 0:
        raise RuntimeError("monotone order-swap difference failed")
    return {
        "example": (
            "h_1(x)=-x and h_2(x)=-x+sin(x)/2 both decrease strictly, "
            "because h_1'=-1 and -3/2<=h_2'<=-1/2."
        ),
        "swaps": (
            "Nevertheless h_2-h_1=sin(x)/2 changes sign at every "
            "successive multiple of pi."
        ),
        "boundary": (
            "This hypothetical exact pair is not an Xi carrier pair. "
            "It proves only that individual monotonicity, as a theorem "
            "type, cannot by itself bound slope-order changes."
        ),
    }


def build_exact() -> dict:
    current = carrier_current_audit()
    constants = constant_audit()
    terminal = terminal_audit()
    order_swap = order_swap_audit()
    return {
        "pi_provenance": (
            "The pi in x=4*pi*exp(L), "
            "a^2=x/(4*pi)+t/16, and u_x=1/(8*pi*a^2) comes from "
            "the completed-zeta factor "
            "pi^(-s/2)*Gamma(s/2)*zeta(s) and the Riemann-Siegel "
            "saddle a^2=T_0/(2*pi). The use of pi>3 in one rational "
            "budget is an elementary bound on that same constant. "
            "No circle fitted to the carriers and no prefix polygon "
            "introduces pi."
        ),
        "domain": (
            "Fix one N=floor(a) chart on L>=50, 0<tL<=25, "
            "q=2*t*L^2>=1. The interior theorem applies to "
            "1<=n<=N-1; when N=1 it is vacuous. It includes q=1, "
            "where t=1/(2L^2) and 0<tL=1/(2L)<=1/100."
        ),
        "carrier_jet": current["carrier_jet"],
        "components": current["components"],
        "division_free_current": (
            current["determinant"] + " " + current["cancellation"]
        ),
        "effective_slope": current["ratio"],
        "source_bounds": (
            "The source gives "
            "nu_n<=8449/x^2-u_n/2, -b>=1/2, |b|<1, and "
            "u_x=1/(8*pi*a^2). Also "
            "s_*''=-t*alpha''(s)/8. Since "
            "alpha''(s)=1/s^3+2/(s-1)^3-1/(2s^2), "
            "|s|=|s-1|=sqrt(1+x^2)/2, t<=1/2, and x>1, "
            "|alpha''|<=26/x^2 and hence "
            "|c_(s,x)|,|b_x|<=13/(8*x^2)."
        ),
        "interior_envelope": (
            "For n<=N-1, "
            "u_n>=log(N/(N-1))>1/N>=1/a and "
            "0<=k_n=log(N/n)<L. Since a^2<2*exp(L), "
            "u_n^2>exp(-L)/2. The monotone boundary checks "
            "13L<pi^2*exp(L) and 33796a<x^2 follow for L>=50 "
            "from "
            f"{constants['curvature_integer_guard']} and "
            f"{constants['phase_integer_guard']}. Therefore "
            "8449/x^2<u_n/4, "
            "|c_(s,x)|k_n<u_n^2/64, "
            "|b_x|u_n<u_n^2/64, and "
            "u_x<u_n^2/(8*pi)."
        ),
        "matrix_budget": (
            "Put B_n=b_x*u_n+b*u_x and K_n=b*u_n*nu_n. Then "
            "nu_n<-u_n/4, K_n>=u_n^2/8, and "
            f"|B_n|<u_n^2/16 because "
            f"{constants['mixed_budget_using_pi_gt_3']}. Moreover "
            "J_n=[X_n,Y_n]M_n[X_n,Y_n]^T with "
            "M_n=[[c_(s,x)k_n-K_n,-B_n/2],[-B_n/2,-K_n]]."
        ),
        "negative_definite_theorem": (
            "Gershgorin gives "
            f"{constants['gershgorin_row_one']} and "
            f"{constants['gershgorin_row_two']}. Hence every eigenvalue "
            "of M_n is at most -5*u_n^2/64, and "
            "J_n<=-(5/64)*u_n^2*|z_n|^2<0 for 1<=n<=N-1."
        ),
        "projective_flow": (
            "The interior vector v_n=(c_n,d_n) never vanishes: if "
            "X_n=0, then Y_n!=0 and d_n=-b*u_n*Y_n!=0. Its "
            "projective angle psi_n=arg(c_n+i*d_n) therefore satisfies "
            "psi_(n,x)=J_n/(c_n^2+d_n^2)<0. On X_n!=0, "
            "h_(n,x)=J_n/X_n^2<0."
        ),
        "zero_projection_join": (
            "At X_n=0 the ratio h_n is not used, but "
            "J_n=-K_n*Y_n^2<=-(1/8)*u_n^2*Y_n^2<0. Thus the "
            "division-free projective orientation crosses every "
            "interior real-projection zero with the same strict sign."
        ),
        "terminal_obstruction": " ".join(terminal.values()),
        "edge_block": (
            "Define the exceptional edge block "
            "C_edge=c_0+c_N and D_edge=d_0+d_N. Then "
            "mathsf_X=C_edge+sum_(n=1)^(N-1)c_n and "
            "mathcal_C_N=D_edge+sum_(n=1)^(N-1)d_n exactly. "
            "If C_edge=0, retain D_edge division-free. This packages "
            "the only carrier not covered by the interior theorem with "
            "the recurrent endpoint instead of assigning either one an "
            "unsupported sign."
        ),
        "endpoint_recurrence": (
            "Distinguish two source objects that were both denoted J_a. "
            "Inside one chart, J_a^(der)=H_(a,x)+mu_aH_a is generally "
            "complex, and h_0=Re(J_a^(der))/H_a"
            "-Im(J_a^(der))/(T_0H_a)-c_s*u_N only when H_a!=0. "
            "At an adjacent cutoff, "
            "J_a^(adj)=H_a(p+2)+H_a(p) and "
            "Q_N=f_(N+1)+kappa_NJ_a^(adj)"
            "=E_[1],N+1-E_[1],N. Recompute the terminal centering and "
            "transport the aggregate by "
            "mathcal_C_(N+1)-mathcal_C_N="
            "Delta mathsf_A-alpha_(N+1)*Delta mathsf_X"
            "+(alpha_N-alpha_(N+1))*mathsf_X_N, with "
            "Delta A_a=Re[Q_(N,x)-lambda_aQ_N"
            "-e_(N+1)d_(N+1,x)]. Never identify J_a^(der) with "
            "J_a^(adj), and infer no separate endpoint or terminal "
            "projective sign."
        ),
        "order_swap_guard": " ".join(order_swap.values()),
        "cumulative_mass_boundary": (
            "Strict clockwise flow of every interior atom does not sign "
            "the cumulative real masses P_k in the slope-ordered Abel "
            "identity, does not order different h_n, and does not bound "
            "their order changes. The remaining source theorem is still "
            "an Xi-specific one-sided cumulative-mass estimate, now with "
            "N-1 monotone interior projective atoms and one exact "
            "exceptional edge block."
        ),
        "route_decision": (
            "Retain the negative-definite interior current as a new "
            "analytic lemma. In the next stage, formulate the aggregate "
            "transport on oriented projective cells for n<=N-1, but "
            "carry (C_edge,D_edge) by the exact endpoint recurrence. "
            "Seek summation-by-parts or oscillatory control of the actual "
            "Xi cumulative masses; do not infer it from monotonicity "
            "alone or split the terminal and endpoint signs."
        ),
        "q_lt_1": (
            "The q<1 layer remains a separate multiplicity-compatible "
            "parabolic/Hermite or boundary-degree theorem. The interior "
            "current calculation is not promoted to t=0 endpoint "
            "simplicity."
        ),
        "proof_boundary": (
            "This proves the exact division-free carrier current, "
            "radial-current cancellation, quantitative negative "
            "definiteness and strict projective orientation for every "
            "interior carrier on the stated q>=1 domain, including q=1. "
            "It also proves only an algebraic terminal sign obstruction, "
            "an exact edge-block decomposition, and a generic order-swap "
            "nonpromotion guard. It does not prove an Xi cumulative-mass "
            "estimate, an edge-block sign, the Abel-scalar gap, successor "
            "count, q<1 closure, finite-height effectivity, contact "
            "exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize "
            "conclusion."
        ),
        "diagnostics": {
            "constant_audit": constants,
            "terminal_sign_guard": terminal["guard"],
            "order_swap_guard": order_swap,
        },
    }


def build_rows(exact: dict) -> list[GateRow]:
    rows = [
        (
            "pi_provenance",
            "definition_provenance",
            "certified",
            "Every occurrence of pi remains source-traced.",
            "No visual or fitted geometry defines pi.",
        ),
        (
            "domain",
            "exact_domain",
            "certified",
            "The theorem includes the q=1 boundary.",
            "The q<1 layer remains separate.",
        ),
        (
            "carrier_jet",
            "exact_identity",
            "ready_to_apply",
            "Each absolute carrier has explicit radial and angular currents.",
            "No bound on the radial current is assumed.",
        ),
        (
            "components",
            "exact_decomposition",
            "ready_to_apply",
            "Terminal centering freezes k_n within one N chart.",
            "This row applies before any slope division.",
        ),
        (
            "division_free_current",
            "exact_identity",
            "ready_to_apply",
            "The projective determinant cancels the radial current exactly.",
            "All d_n phase corrections remain inside nu_n.",
        ),
        (
            "effective_slope",
            "exact_identity",
            "ready_to_apply",
            "The ratio derivative agrees with the determinant off its poles.",
            "The ratio is not used at X_n=0.",
        ),
        (
            "source_bounds",
            "analytic_bound",
            "certified",
            "The first coefficient-current bounds control the moving frame.",
            "Uses L>=50 and 0<tL<=25.",
        ),
        (
            "interior_envelope",
            "analytic_bound",
            "certified",
            "Every interior logarithmic distance dominates the frame errors.",
            "The terminal distance is intentionally excluded.",
        ),
        (
            "matrix_budget",
            "analytic_reduction",
            "ready_to_apply",
            "The carrier determinant is a uniformly negative quadratic form.",
            "The constants are deliberately nonoptimized.",
        ),
        (
            "negative_definite_theorem",
            "analytic_lemma",
            "proved",
            "Every interior projective current has one strict sign.",
            "Only 1<=n<=N-1 is covered.",
        ),
        (
            "projective_flow",
            "analytic_corollary",
            "proved",
            "Every interior component moves strictly clockwise in RP1.",
            "This is an individual-atom theorem.",
        ),
        (
            "zero_projection_join",
            "division_free_join",
            "proved",
            "Interior projection zeros have a fixed oriented join.",
            "No tangent ratio is evaluated at the zero.",
        ),
        (
            "terminal_obstruction",
            "nonpromotion_guard",
            "guard_validated",
            "The isolated terminal current has no source-level uniform sign.",
            "The sign examples are not asserted Xi states.",
        ),
        (
            "edge_block",
            "exact_decomposition",
            "ready_to_apply",
            "Terminal and endpoint contributions form one exact exceptional block.",
            "No sign or nonvanishing is asserted for the edge block.",
        ),
        (
            "endpoint_recurrence",
            "exact_boundary_law",
            "ready_to_apply",
            "The exceptional edge is transported only through the recurrent real jet.",
            "J_a is generally complex.",
        ),
        (
            "order_swap_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Individual monotonicity alone does not control order swaps.",
            "The hypothetical functions are not Xi carriers.",
        ),
        (
            "cumulative_mass_boundary",
            "open_theorem_target",
            "not_ready_to_apply",
            "The remaining gap is still a signed Xi cumulative-mass estimate.",
            "The new lemma narrows but does not close that target.",
        ),
        (
            "route_decision",
            "route_decision",
            "ready_to_apply",
            "Future work should use interior projective cells plus one recurrent edge block.",
            "Do not split unsupported terminal and endpoint signs.",
        ),
        (
            "q_lt_1",
            "open_theorem_target",
            "not_ready_to_apply",
            "The multiplicity-compatible small-q theorem remains separate.",
            "No t=0 simplicity is assumed.",
        ),
        (
            "proof_boundary",
            "proof_guard",
            "guard_validated",
            "The interior theorem is not promoted to the aggregate Xi gap.",
            "RH and every prize-level conclusion remain open.",
        ),
    ]
    return [
        GateRow(
            id=f"ipcg_{index:02d}_{suffix}",
            role=role,
            readiness=readiness,
            claim=claim,
            formula=exact[suffix],
            proof_boundary=boundary,
            diagnostics=(
                exact["diagnostics"]
                if suffix in {"negative_definite_theorem", "terminal_obstruction", "order_swap_guard"}
                else None
            ),
        )
        for index, (suffix, role, readiness, claim, boundary) in enumerate(rows)
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "proved negative-definite projective current for every "
            "interior carrier; exact terminal-endpoint edge isolation; "
            "aggregate Xi cumulative-mass theorem remains open"
        ),
        "proof_boundary": exact["proof_boundary"],
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in build_rows(exact)],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Interior Carrier Projective-Current Gate",
            "",
            "Date: 2026-07-28",
            "",
            "Status: proved interior-carrier monotonicity and exact edge isolation;",
            "not a proof of the Xi cumulative-mass gap, `Lambda<=0`, RH,",
            "PF-infinity, or a Clay-prize result.",
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
            "## Carrier Current",
            "",
            "```text",
            exact["carrier_jet"],
            exact["components"],
            exact["division_free_current"],
            exact["effective_slope"],
            "```",
            "",
            "## Uniform Interior Bound",
            "",
            exact["source_bounds"],
            "",
            exact["interior_envelope"],
            "",
            "```text",
            exact["matrix_budget"],
            exact["negative_definite_theorem"],
            "```",
            "",
            "## Projective Join",
            "",
            exact["projective_flow"],
            "",
            exact["zero_projection_join"],
            "",
            "## Exceptional Edge",
            "",
            exact["terminal_obstruction"],
            "",
            "```text",
            exact["edge_block"],
            "```",
            "",
            exact["endpoint_recurrence"],
            "",
            "## Order-Swap Guard",
            "",
            exact["order_swap_guard"],
            "",
            "## Remaining Theorem",
            "",
            exact["cumulative_mass_boundary"],
            "",
            exact["route_decision"],
            "",
            exact["q_lt_1"],
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
        "built interior projective-current gate: "
        "20 rows, exact carrier determinant, proved 5/64 interior "
        "negative-definite bound, terminal-edge and order-swap guards"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
