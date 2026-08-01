#!/usr/bin/env python3
"""Build the first-order rectangular boundary-degree reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_rectangular_boundary_degree_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "target_reconciliation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_target_reconciliation_gate.json"
    ),
    "boundary_degree": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
    ),
    "first_order_boundary": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_cofinal_boundary_reduction.json"
    ),
    "contact_index": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_first_jet_winding_gate.json"
    ),
    "oriented_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_oriented_successor_winding_reduction.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_carrier_kernel_abel_prefix_reduction.json"
    ),
    "abel_shear": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_abel_scalar_shear_flux_reduction.json"
    ),
    "oscillatory_splice": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_oscillatory_spliced_outer_collar_reduction.json"
    ),
    "connector_cap": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_dominant_ray_connector_phase_cap.json"
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


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    target = payloads["target_reconciliation"].get("exact", {})
    if target.get("box_optimal_target", {}).get("exact_target") != (
        "G_L(X,U)>1, equivalently "
        "|X|<=delta_0(L) => |U|>gamma_0(L)"
    ):
        raise RuntimeError("box-optimal target drifted")
    if target.get("abel_target", {}).get("pointwise") != (
        "|mathsf_X|<=delta_L => "
        "|mathcal_C_N|>A_L+epsilon_term"
    ):
        raise RuntimeError("normalized Abel target drifted")

    first_boundary = payloads["first_order_boundary"].get("exact", {})
    if first_boundary.get("vector_remainder") != (
        "V_r=(r_[1],r_[1],x/L), "
        "||V_r||_2^2<50000000000*exp(-5L/2)"
    ):
        raise RuntimeError("first-order vector remainder drifted")
    if "1300000000*exp(-5L/2)" not in (
        first_boundary.get("adjacent_vector", "")
    ):
        raise RuntimeError("adjacent vector budget drifted")

    contact = payloads["contact_index"].get("exact", {})
    if "+floor(m/2)" not in contact.get("multiplicity_index", ""):
        raise RuntimeError("positive contact index drifted")
    if "iff Omega contains no real multiple-zero contact" not in (
        contact.get("zero_winding_equivalence", "")
    ):
        raise RuntimeError("zero-winding equivalence drifted")

    degree = payloads["boundary_degree"].get("exact", {})
    if "only its boundary winding enters" not in (
        degree.get("zero_degree_exclusion", {}).get("scope", "")
    ):
        raise RuntimeError("boundary-degree scope drifted")

    successor = payloads["oriented_successor"].get("exact", {})
    if "strictly weaker estimate wind(proxy_j)<1" not in (
        successor.get("successor", {}).get("one_sided_trap", "")
    ):
        raise RuntimeError("successor integer trap drifted")

    prefix = payloads["abel_prefix"].get("exact", {})
    if "mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X" not in (
        prefix.get("contact_scalar", "")
    ):
        raise RuntimeError("Abel terminal shear drifted")

    shear = payloads["abel_shear"].get("exact", {})
    if "|mathcal_C_N|>A_L+epsilon_term" not in (
        shear.get("q_ge_1_pointwise_target", "")
    ):
        raise RuntimeError("Abel shear target drifted")

    splice = payloads["oscillatory_splice"].get("exact", {})
    outer = splice.get("open_outer_target", "")
    if "0<tL<=c_*+epsilon" not in outer or "q=2tL^2>=1" not in outer:
        raise RuntimeError("oscillatory splice drifted")

    connector = payloads["connector_cap"].get("exact", {})
    if "<pi/2" not in connector.get("connector_phase_cap", ""):
        raise RuntimeError("connector phase cap drifted")
    if "H_j<3*pi/2" not in connector.get(
        "reduced_horizontal_target", ""
    ):
        raise RuntimeError("horizontal phase target drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "source_kinds": {
            key: payload.get("kind") for key, payload in payloads.items()
        },
    }


def build_rows() -> list[GateRow]:
    return [
        GateRow(
            "rbd_01_current_main",
            "source contract",
            "proved",
            "Use the cutoff-local first-order main and its matched remainder.",
            "Z=J_[1]+r_[1], J_[1]=2X, J_[1],x=2U.",
            "This is local to one prescribed cutoff chart.",
        ),
        GateRow(
            "rbd_02_error_rectangle",
            "exact error geometry",
            "proved",
            "The C1 remainder gives a rectangular first-jet error set.",
            "|r_[1]|<B_0=100000e^(-5L/4), "
            "|r_[1],x/L|<B_1=200000e^(-5L/4).",
            "The two coordinates must not be merged before choosing a target.",
        ),
        GateRow(
            "rbd_03_box_gauge",
            "exact target",
            "proved",
            "The error-matched main gauge is the contact-box gauge.",
            "M_L(J)=max(|J|/B_0,|J_x/L|/B_1)"
            "=G_L(X,U).",
            "Strict M_L>1 is required for a robust homotopy.",
        ),
        GateRow(
            "rbd_04_rectangular_rouche",
            "boundary homotopy",
            "proved",
            "M_L(J)>1 transfers boundary nonvanishing and winding to Z.",
            "M_L(J+s*r)>=M_L(J)-s*M_L(r)>0 for 0<=s<=1.",
            "This is a boundary theorem; it gives no interior noncontact.",
        ),
        GateRow(
            "rbd_05_radial_overstrength",
            "nonpromotion guard",
            "proved",
            "The Euclidean radial margin implies but is not equivalent to the box margin.",
            "||V_J||_2^2>B_0^2+B_1^2=50000000000e^(-5L/2)"
            " => M_L(J)>1; the converse fails.",
            "Do not retain the stronger radial target by habit.",
        ),
        GateRow(
            "rbd_06_error_whitening",
            "exact coordinate change",
            "proved",
            "A positive diagonal map sends the error rectangle to the unit square.",
            "D_L=diag(B_0^(-1),B_1^(-1)), det(D_L)>0, "
            "D_L V_J=(J/B_0,J_x/(L B_1)).",
            "The map changes phase parametrization but not degree.",
        ),
        GateRow(
            "rbd_07_winding_invariance",
            "topological transfer",
            "proved",
            "Continuous positive diagonal whitening preserves winding.",
            "D_(L,s)=(1-s)I+sD_L is invertible with positive determinant.",
            "No metric lower bound follows from orientation preservation alone.",
        ),
        GateRow(
            "rbd_08_boundary_abel",
            "arithmetic handoff",
            "open",
            "The normalized Abel gap implies M_L(J)>1 on any selected boundary arc.",
            "|mathsf_X|<=delta_L => "
            "|mathcal_C_N|>A_L+epsilon_term.",
            "The Xi-specific Abel inequality remains unproved.",
        ),
        GateRow(
            "rbd_09_zero_fibre",
            "exceptional fibre",
            "retained",
            "The boundary Abel target remains division-free at W_0=0.",
            "mathsf_A=mathcal_C_N+c*u_N*mathsf_X.",
            "No quotient by W_0 or the first-jet norm is permitted.",
        ),
        GateRow(
            "rbd_10_adjacent_box",
            "cutoff arithmetic",
            "proved",
            "The adjacent first-order chart jump occupies a strict sub-box.",
            "|Delta J|/B_0<1/5, |Delta J_x/L|/B_1<3/20.",
            "These are overlap bounds, not a global fixed-cutoff theorem.",
        ),
        GateRow(
            "rbd_11_adjacent_homotopy",
            "chart topology",
            "proved",
            "One box-separated chart joins its adjacent chart without zero.",
            "M_L(J+s*Delta J)>=M_L(J)-s/5>4/5.",
            "Every prescribed equality boundary and overlap is still included.",
        ),
        GateRow(
            "rbd_12_normalizer",
            "exact topology",
            "proved",
            "The positive Xi normalizer and its derivative shear preserve degree.",
            "V_H=A_t[[1,0],[a/L,1]]V_Z, det=A_t^2>0.",
            "The normalizer is not discarded from quantitative estimates.",
        ),
        GateRow(
            "rbd_13_contact_index",
            "exact heat theorem",
            "proved",
            "Every real heat contact has strictly positive local index in (x,t).",
            "ind_(x,t)(H,H_x/L)=floor(m/2)>0 for multiplicity m>=2.",
            "The orientation convention must remain counterclockwise in (x,t).",
        ),
        GateRow(
            "rbd_14_zero_degree",
            "global degree",
            "conditional",
            "Boundary transfer plus zero winding excludes all interior contacts.",
            "wind(V_J)=wind(V_H)=sum_p floor(m_p/2); wind(V_J)=0 => no p.",
            "The proxy may have interior zeros; only its boundary winding is used.",
        ),
        GateRow(
            "rbd_15_integer_trap",
            "successor degree",
            "conditional",
            "A sub-one proxy winding is enough on each oriented successor.",
            "0<=kappa_j=wind(V_J)<1 => kappa_j=0.",
            "Boundary nonvanishing alone does not prove kappa_j<1.",
        ),
        GateRow(
            "rbd_16_ray_schedule",
            "cofinal geometry",
            "proved",
            "Use the ray-aligned successor bottoms.",
            "t_j=25/(100+j), L_j=101+j, R_j=4*pi*exp(L_j).",
            "Finite base and shoulder ownership remain separate.",
        ),
        GateRow(
            "rbd_17_hard_bottom_interval",
            "one-dimensional reduction",
            "proved",
            "The low-c q>=1 Abel obligation on bottom j is one explicit L interval.",
            "max(B_epsilon,sqrt((100+j)/50))<=L"
            "<=min(L_j,(c_*+epsilon)(100+j)/25).",
            "An empty interval contributes no hard bottom arc.",
        ),
        GateRow(
            "rbd_18_outer_closure",
            "source composition",
            "proved",
            "Oscillatory and dominant theorems close the rest of the raised outer bottom.",
            "c_*+epsilon<t_jL<25 is oscillatory; t_jL>=25 is dominant.",
            "q<1 and L<B_epsilon remain inner or finite-shoulder obligations.",
        ),
        GateRow(
            "rbd_19_uniform_route",
            "route A",
            "open",
            "A two-dimensional Abel gap gives complete outer-collar noncontact.",
            "Uniform gap on A_(j,epsilon) => deg(Omega_j^out)=0.",
            "This stronger route removes the independent outer winding theorem.",
        ),
        GateRow(
            "rbd_20_boundary_route",
            "route B",
            "open",
            "A one-dimensional boundary Abel gap can be paired with the horizontal phase cap.",
            "M_L>1 on joined boundary and H_j<3*pi/2, C_j<pi/2"
            " => kappa_j<1.",
            "This weaker-domain route retains an independent winding obligation.",
        ),
        GateRow(
            "rbd_21_boundary_guard",
            "nonpromotion guard",
            "proved",
            "Boundary nonvanishing without a winding bound does not exclude contact.",
            "F=x^2-2t has a boundary-nonzero rectangle with one index +1 contact.",
            "This is a generic heat guard, not an Xi counterexample.",
        ),
        GateRow(
            "rbd_22_handoff",
            "route decision",
            "open",
            "Attack the error-matched boundary Abel theorem before a two-dimensional upgrade.",
            "Prove the Abel gap on every nonempty hard bottom interval, "
            "then compose the joined horizontal flux below 3*pi/2.",
            "No contact exclusion, Lambda<=0, PF-infinity, RH, or prize conclusion is proved.",
        ),
    ]


def build_payload() -> dict:
    payloads = load_sources()
    source_audit = audit_sources(payloads)
    rows = build_rows()
    proof_boundary = (
        "This artifact proves the rectangular first-order boundary homotopy, "
        "error-whitening and winding invariance, adjacent-chart sub-box join, "
        "positive-index degree transfer, exact hard-bottom interval, and the "
        "two-route comparison. It does not prove the Xi Abel gap, horizontal "
        "phase bound, q<1 or bounded-L inner theorem, complete boundary "
        "nonvanishing, contact exclusion, Lambda<=0, PF-infinity, RH, or a "
        "Clay-prize conclusion."
    )
    return {
        "kind": STEM,
        "date": "2026-07-29",
        "status": (
            "exact_rectangular_boundary_degree_reduction_with_open_"
            "abel_and_winding_targets"
        ),
        "proof_boundary": proof_boundary,
        "source_audit": source_audit,
        "exact": {
            "domain": (
                "L>=50, 0<tL<25 on prescribed first-order charts; "
                "the raised outer boundary also uses L>=B_epsilon"
            ),
            "remainder_box": {
                "value": "B_0(L)=100000*exp(-5L/4)",
                "derivative": "B_1(L)=200000*exp(-5L/4)",
                "error": (
                    "|r_[1]|<B_0(L), "
                    "|partial_x r_[1]/L|<B_1(L)"
                ),
            },
            "main_gauge": {
                "definition": (
                    "M_L(J)=max(|J_[1]|/B_0(L),"
                    "|partial_x J_[1]/L|/B_1(L))"
                ),
                "complex_half": "M_L(J)=G_L(X,U)",
                "boundary_target": "M_L(J)>1",
            },
            "rectangular_rouche": {
                "homotopy": (
                    "V_s=V_J+s*V_r, "
                    "M_L(V_s)>=M_L(V_J)-s*M_L(V_r)>0"
                ),
                "conclusion": (
                    "V_J and V_Z are boundary-nonzero and have equal winding"
                ),
            },
            "error_whitening": {
                "matrix": "D_L=diag(B_0(L)^(-1),B_1(L)^(-1))",
                "determinant": "det(D_L)=1/(B_0(L)*B_1(L))>0",
                "square": "D_L maps the certified error box into (-1,1)^2",
                "winding": (
                    "The positive diagonal homotopy from I to D_L "
                    "preserves nonvanishing and winding"
                ),
            },
            "radial_comparison": {
                "sufficient": (
                    "J_[1]^2+(partial_x J_[1]/L)^2"
                    ">50000000000*exp(-5L/2) => M_L(J)>1"
                ),
                "constant": "100000^2+200000^2=50000000000",
                "strictness": (
                    "The radial condition is strictly stronger than M_L(J)>1"
                ),
            },
            "boundary_abel": {
                "target": (
                    "|mathsf_X|<=delta_L => "
                    "|mathcal_C_N|>A_L+epsilon_term"
                ),
                "conclusion": (
                    "The target implies M_L(J)>1 on the selected boundary, "
                    "including W_0=0"
                ),
                "prefix": (
                    "U=u_N*S+sum_(k=1)^(N-1)"
                    "log((k+1)/k)*F_k"
                ),
            },
            "adjacent_chart": {
                "value_ratio": "|Delta J_[1]|/B_0(L)<1/5",
                "derivative_ratio": (
                    "|partial_x Delta J_[1]/L|/B_1(L)<3/20"
                ),
                "gauge": "M_L(Delta J)<1/5",
                "homotopy": (
                    "M_L(J+s*Delta J)>1-s/5>=4/5"
                ),
            },
            "degree_transfer": {
                "normalizer": (
                    "V_H=A_t*[[1,0],[partial_x log(A_t)/L,1]]*V_Z"
                ),
                "contact_index": (
                    "ind_(x,t)(H,H_x/L)=floor(m/2)>0"
                ),
                "zero_degree": (
                    "wind(V_J)=0 implies no interior heat contact"
                ),
                "successor_trap": "0<=kappa_j=wind(V_J)<1 => kappa_j=0",
            },
            "ray_bottom": {
                "schedule": (
                    "t_j=25/(100+j), L_j=101+j, "
                    "R_j=4*pi*exp(L_j)"
                ),
                "q": "q=50*L^2/(100+j)",
                "c": "t_j*L=25*L/(100+j)",
                "hard_interval": (
                    "max(B_epsilon,sqrt((100+j)/50))<=L"
                    "<=min(L_j,(c_*+epsilon)*(100+j)/25)"
                ),
            },
            "route_fork": {
                "uniform": (
                    "Prove the Abel gap on the full two-dimensional "
                    "raised outer wedge; excise it at zero degree"
                ),
                "boundary": (
                    "Prove the Abel gap only on joined successor boundaries "
                    "and H_j<3*pi/2; use C_j<pi/2 and the integer trap"
                ),
                "recommended_next": (
                    "Attempt the one-dimensional hard-bottom Abel theorem "
                    "first, preserving the full joined phase ledger"
                ),
            },
            "nonpromotion": {
                "radial": (
                    "Box separation does not imply the stronger Euclidean "
                    "radial threshold"
                ),
                "boundary_only": (
                    "Boundary nonvanishing without a winding bound permits "
                    "positive-index interior heat contacts"
                ),
            },
            "open_handoff": (
                "For every successor j and fixed epsilon, prove the normalized "
                "Abel gap on each nonempty hard bottom interval "
                "max(B_epsilon,sqrt((100+j)/50))<=L<="
                "min(L_j,(c_*+epsilon)(100+j)/25), retaining W_0=0, "
                "the recurrent endpoint, q=1, equality boundaries, and "
                "adjacent-cutoff homotopies. Then prove the completely joined "
                "horizontal phase contribution H_j<3*pi/2. If the same Abel "
                "estimate extends to the full raised outer wedge, use the "
                "uniform-collar excision route instead."
            ),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "error_box_coordinates": 2,
            "rectangular_homotopies": 2,
            "positive_index_degree_transfers": 1,
            "boundary_only_abel_targets": 1,
            "route_branches": 2,
            "nonpromotion_guards": 2,
            "pointwise_contact_exclusions": 0,
            "cofinal_descendant_theorem": False,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# First-Order Rectangular Boundary-Degree Reduction",
        "",
        "Date: 2026-07-29",
        "",
        "Status: exact boundary and degree reduction with open Xi Abel and",
        "winding targets. This is not a proof of contact exclusion,",
        "Lambda<=0, or RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Newman first-order rectangular boundary-degree "
            "reduction: 22 rows, 2 error-box coordinates, "
            "2 rectangular homotopies, 1 positive-index degree transfer, "
            "1 boundary-only Abel target, 2 route branches, "
            "2 nonpromotion guards, 0 contact exclusions"
        ),
        "```",
        "",
        "## Error-Matched Boundary Geometry",
        "",
        "The first-order remainder is componentwise:",
        "",
        "```text",
        exact["remainder_box"]["error"],
        "```",
        "",
        "Define",
        "",
        "```text",
        exact["main_gauge"]["definition"],
        exact["main_gauge"]["complex_half"],
        "```",
        "",
        "Then `M_L(J)>1` is exactly the box-optimal signed-band condition.",
        "On a boundary it gives the componentwise Rouche homotopy",
        "",
        "```text",
        exact["rectangular_rouche"]["homotopy"],
        "```",
        "",
        "so the first-order main and exact normalized Xi jet have the same",
        "boundary winding. The older Euclidean target remains sufficient but",
        "is strictly stronger.",
        "",
        "## Error Whitening",
        "",
        "```text",
        exact["error_whitening"]["matrix"],
        exact["error_whitening"]["determinant"],
        exact["error_whitening"]["square"],
        "```",
        "",
        "The positive diagonal homotopy preserves winding. This packages the",
        "unequal C1 value and derivative errors without replacing their",
        "rectangle by a larger ellipse.",
        "",
        "## Abel Boundary Target",
        "",
        "The division-free target is unchanged but is required only on the",
        "selected boundary arcs:",
        "",
        "```text",
        exact["boundary_abel"]["target"],
        exact["boundary_abel"]["prefix"],
        "```",
        "",
        exact["boundary_abel"]["conclusion"] + ".",
        "",
        "Adjacent charts already lie in a strict sub-box:",
        "",
        "```text",
        exact["adjacent_chart"]["value_ratio"],
        exact["adjacent_chart"]["derivative_ratio"],
        exact["adjacent_chart"]["homotopy"],
        "```",
        "",
        "## Degree Transfer",
        "",
        "The positive normalizer shear preserves degree, and every real heat",
        "contact has local index `floor(m/2)>0` in the standard `(x,t)`",
        "orientation. Consequently boundary winding zero excludes all",
        "contacts. On one successor the weaker integer trap is enough:",
        "",
        "```text",
        exact["degree_transfer"]["successor_trap"],
        "```",
        "",
        "Boundary nonvanishing alone is not enough: `F=x^2-2t` has a",
        "boundary-nonzero rectangle with one positive-index contact.",
        "",
        "## Hard Bottom Interval",
        "",
        "On the ray schedule,",
        "",
        "```text",
        exact["ray_bottom"]["schedule"],
        exact["ray_bottom"]["q"],
        exact["ray_bottom"]["c"],
        "```",
        "",
        "the remaining low-c, q>=1 bottom interval is exactly",
        "",
        "```text",
        exact["ray_bottom"]["hard_interval"],
        "```",
        "",
        "Oscillatory-zeta closes the raised middle and the dominant theorem",
        "closes t_jL>=25. q<1 and L<B_epsilon remain separate.",
        "",
        "## Route Fork",
        "",
        "Route A:",
        "",
        exact["route_fork"]["uniform"] + ".",
        "",
        "Route B:",
        "",
        exact["route_fork"]["boundary"] + ".",
        "",
        "Recommended next bounded target:",
        "",
        exact["route_fork"]["recommended_next"] + ".",
        "",
        "Open handoff:",
        "",
        exact["open_handoff"],
        "",
        "Pi provenance is unchanged: `pi` comes from the completed-zeta,",
        "Riemann-Siegel, Fourier, and phase-period normalizations already",
        "recorded in the source chain.",
        "",
        "## Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built Newman first-order rectangular boundary-degree reduction: "
        "22 rows, 2 error-box coordinates, 2 rectangular homotopies, "
        "1 positive-index degree transfer, 1 boundary-only Abel target, "
        "2 route branches, 2 nonpromotion guards, 0 contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
