#!/usr/bin/env python3
"""Build the oriented successor-winding reduction for the first-order proxy."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_oriented_successor_winding_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "first_jet_winding": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_first_jet_winding_gate.json"
    ),
    "boundary_degree": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
    ),
    "q208": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_closed_boundary_winding_certificate.json"
    ),
    "ray_schedule": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.json"
    ),
    "dominant_ray": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "dominant_saddle_global_ray_certificate.json"
    ),
    "first_order_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_global_remainder_certificate.json"
    ),
    "first_order_cofinal": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_cofinal_boundary_reduction.json"
    ),
}


@dataclass(frozen=True)
class ReductionRow:
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
        raise FileNotFoundError("missing source artifacts: " + ", ".join(missing))
    return {key: file_hash(path) for key, path in SOURCES.items()}


def source_audit() -> dict[str, str]:
    payloads = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }
    markers = {
        "first_jet_winding": (
            '"domain_orientation": "(x,t)"',
            "+floor(m/2)",
            "+sum_(p in Omega)",
        ),
        "boundary_degree": (
            "strictly positive local index",
            "wind(V_P)=0",
        ),
        "q208": (
            '"exact_winding": 0',
            "Q_1 through Q_208",
        ),
        "ray_schedule": (
            '"base_time": "1/4"',
            "t_(j+1)*L(x)>=t_(j+1)*L_j",
            "S_j=[t_(j+1),1/2]x[R_j,R_(j+1)]",
        ),
        "dominant_ray": (
            "tL>=25",
            "L_t(x)>0",
        ),
        "first_order_remainder": (
            '"eta_0": 100000',
            '"eta_1": 200000',
        ),
        "first_order_cofinal": (
            "50000000000",
            "1300000000",
            "13/500",
            "t_j=25/(100+j)",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(payload.get("kind", "")) for key, payload in payloads.items()}


def phase_audit() -> dict:
    value, value_x, value_xx = sp.symbols(
        "J J_x J_xx", real=True
    )
    value_t, value_xt = sp.symbols("J_t J_xt", real=True)
    scale = sp.symbols("ell", positive=True)
    scale_x, scale_t = sp.symbols("ell_x ell_t", real=True)
    norm_sq = value**2 + (value_x / scale) ** 2
    phase_x = sp.factor(
        (
            value
            * (value_xx / scale - value_x * scale_x / scale**2)
            - (value_x / scale) * value_x
        )
        / norm_sq
    )
    phase_t = sp.factor(
        (
            value
            * (value_xt / scale - value_x * scale_t / scale**2)
            - (value_x / scale) * value_t
        )
        / norm_sq
    )
    expected_x = (
        scale * value * value_xx
        - value * value_x * scale_x
        - scale * value_x**2
    ) / (scale**2 * norm_sq)
    expected_t = (
        scale * value * value_xt
        - value * value_x * scale_t
        - scale * value_x * value_t
    ) / (scale**2 * norm_sq)
    if sp.simplify(phase_x - expected_x) != 0:
        raise RuntimeError("scaled horizontal Pruefer identity failed")
    if sp.simplify(phase_t - expected_t) != 0:
        raise RuntimeError("scaled vertical Pruefer identity failed")
    if sp.simplify(phase_x.subs(value, 0) + scale) != 0:
        raise RuntimeError("simple-zero phase orientation failed")
    return {
        "proxy": "W_J=J+i*J_x/ell, ell(t,x)>0",
        "norm": "Q_ell=J^2+(J_x/ell)^2",
        "phase_x": (
            "partial_x arg(W_J)="
            "(ell*J*J_xx-J*J_x*ell_x-ell*J_x^2)/(ell^2*Q_ell)"
        ),
        "phase_t": (
            "partial_t arg(W_J)="
            "(ell*J*J_xt-J*J_x*ell_t-ell*J_x*J_t)/(ell^2*Q_ell)"
        ),
        "simple_zero": (
            "At J=0 and J_x!=0, partial_x arg(W_J)=-ell<0."
        ),
        "sympy_phase_x": str(phase_x),
        "sympy_phase_t": str(phase_t),
    }


def chain_audit() -> dict:
    p_j = Counter(
        {
            "bottom_old": 1,
            "right_upper": 1,
            "top_old": 1,
            "axis_upper": 1,
        }
    )
    descendant = Counter(
        {
            "bottom_new_inner": 1,
            "right_lower": 1,
            "bottom_old": -1,
            "axis_lower": 1,
        }
    )
    outer = Counter(
        {
            "bottom_new_outer": 1,
            "right_next": 1,
            "top_outer": 1,
            "right_lower": -1,
            "right_upper": -1,
        }
    )
    expected = Counter(
        {
            "bottom_new_inner": 1,
            "bottom_new_outer": 1,
            "right_next": 1,
            "top_outer": 1,
            "top_old": 1,
            "axis_upper": 1,
            "axis_lower": 1,
        }
    )
    total = p_j + descendant + outer
    if total != expected:
        raise RuntimeError(f"successor chain cancellation failed: {total}")
    return {
        "rectangles": (
            "P_j=[t_j,T]x[0,R_j], "
            "D_j=[t_(j+1),t_j]x[0,R_j], "
            "S_j=[t_(j+1),T]x[R_j,R_(j+1)]"
        ),
        "oriented_identity": "partial P_(j+1)=partial P_j+partial D_j+partial S_j",
        "p_j_segments": dict(sorted(p_j.items())),
        "descendant_segments": dict(sorted(descendant.items())),
        "outer_segments": dict(sorted(outer.items())),
        "successor_segments": dict(sorted(expected.items())),
        "cancellation_validated": True,
    }


def build_exact() -> dict:
    phase = phase_audit()
    chain = chain_audit()
    return {
        "orientation": {
            "plane": "(x,t)",
            "positive_boundary": (
                "bottom x:0->R; right t:lower->upper; "
                "top x:R->0; axis t:upper->lower"
            ),
            "contact_degree": (
                "ind_(x,t)(H,H_x/ell)=+floor(m/2)>0"
            ),
            "reversal_guard": (
                "Writing coordinates as (t,x) reverses both contact degree "
                "and boundary orientation; it must not be mixed with the "
                "standard counterclockwise edge order."
            ),
        },
        "phase": phase,
        "signed_crossing": {
            "general": (
                "For a closed regular path W=u+iv, wind(W,0) is the oriented "
                "intersection count with the positive imaginary ray: "
                "sum_(u=0,v>0) sign(-d_s u)."
            ),
            "horizontal_forward": (
                "For W_J on x increasing, each zero J=0,J_x>0 contributes -1."
            ),
            "horizontal_reverse": (
                "For the same path traversed with x decreasing, each zero "
                "J=0,J_x>0 contributes +1."
            ),
            "vertical": (
                "On x=R with t increasing, a transversal crossing "
                "J=0,J_x>0 contributes -sign(J_t)."
            ),
            "tangency_convention": (
                "A vertical tangency is interpreted by oriented intersection "
                "number, equivalently by a small regular ray rotation. "
                "Corners are assigned once by a half-open edge convention."
            ),
        },
        "chain": chain,
        "base": {
            "schedule": (
                "t_j=25/(100+j), L_j=101+j, "
                "R_j=4*pi*exp(L_j), T=1/2"
            ),
            "base_rectangle": "P_0=[1/4,1/2]x[0,R_0]",
            "base_degree": (
                "Since 1/4>1/5>=Lambda, P_0 contains no heat contact and "
                "wind(V_H(partial P_0),0)=0."
            ),
            "q208_role": (
                "Q208 is a rigorous compact calibration and core input; "
                "it is not needed to infer the P_0 degree."
            ),
        },
        "outer_strip": {
            "lower_bound": (
                "On S_j, L>=L_j and t>=t_(j+1), so "
                "tL>=t_(j+1)L_j=25 and L>=101."
            ),
            "degree": (
                "The dominant-saddle theorem gives no contact in S_j, hence "
                "wind(V_H(partial S_j),0)=0."
            ),
        },
        "successor": {
            "descendant_loop": (
                "partial D_j=B_(j+1)^[0,R_j]+V_j-reverse(B_j) "
                "+A_j, with V_j: (R_j,t_(j+1))->(R_j,t_j) "
                "and A_j on the positive symmetry axis."
            ),
            "recurrence": (
                "w_(j+1)=w_j+kappa_j, where "
                "kappa_j=wind(V_H(partial D_j),0)."
            ),
            "positive_integer": (
                "kappa_j=sum_(p in D_j) floor(m_p/2) is a nonnegative integer."
            ),
            "crossing_formula": (
                "For a derivative-compatible regular proxy, "
                "kappa_j=N_up(t_j;R_j)-N_up(t_(j+1);R_j)"
                "-sum_(t on V_j; J=0,J_x>0) sign(J_t), "
                "plus the explicitly assigned finite-shoulder/join "
                "intersection count."
            ),
            "one_sided_trap": (
                "Boundary homotopy gives wind(proxy_j)=kappa_j>=0. "
                "Therefore the strictly weaker estimate wind(proxy_j)<1 "
                "(equivalently total oriented phase change <2*pi) forces "
                "kappa_j=0; exact zero winding need not be proved directly."
            ),
        },
        "proxy_transfer": {
            "remainder_budget": (
                "||(r_[1],r_[1],x/L)||_2^2"
                "<50000000000*exp(-5L/2)"
            ),
            "margin": (
                "J_[1]^2+(J_[1],x/L)^2"
                ">50000000000*exp(-5L/2)"
            ),
            "cutoff_ratio": (
                "||Delta V_J||_2^2/50000000000*exp(-5L/2)<13/500<1"
            ),
            "conclusion": (
                "The full first jet, first-order proxy, and adjacent cutoff "
                "charts are joined by nonvanishing boundary homotopies. "
                "They have the same integer winding; cutoff joins add no "
                "independent integer theorem."
            ),
        },
        "open_arithmetic": {
            "frequency_margin": (
                "On q=2t_jL^2>=1, prove the first-order boundary margin "
                "using the saddle-centered logarithmic moment."
            ),
            "parabolic_margin": (
                "On q<1, prove a multiplicity-compatible boundary margin; "
                "its constants may deteriorate as t_j tends to zero."
            ),
            "shoulders": (
                "Certify the finite L<50, oscillatory L_epsilon, core-to-main, "
                "and vertical-connector phase cells with exact joins."
            ),
            "integer_target": (
                "For every j, prove the composite successor proxy has "
                "wind<1, or equivalently an oriented crossing count <=0. "
                "Boundary nonvanishing alone does not imply this inequality."
            ),
        },
        "countermodel": {
            "flow": "F(t,x)=x^2-2t solves F_t=-F_xx.",
            "contact": "(t,x)=(0,0) has standard local degree +1.",
            "boundary": (
                "On [-1,1]x[-1/4,1/4], F+iF_x is boundary-nonzero but "
                "has winding +1."
            ),
            "guard": (
                "A boundary norm floor, top/right nonvanishing, or cutoff "
                "continuity does not determine the winding integer."
            ),
        },
        "proof_handoff": (
            "Construct one continuous successor boundary proxy, use signed "
            "phase cells to upper-bound its winding by a number below 1, and "
            "split only the hard bottom estimates into q>=1, q<1, and finite "
            "shoulders. The positive contact degree then forces zero without "
            "an exact phase evaluation."
        ),
    }


def build_rows(exact: dict) -> list[ReductionRow]:
    return [
        ReductionRow(
            "nfosw_01_orientation",
            "exact_orientation_lemma",
            "available_exact",
            "The boundary traversal and all local contact indices use one standard (x,t) orientation.",
            exact["orientation"]["contact_degree"],
            "Changing coordinate order is allowed only if the boundary orientation is reversed too.",
            exact["orientation"],
        ),
        ReductionRow(
            "nfosw_02_scaled_pruefer_form",
            "exact_differential_identity",
            "available_exact",
            "The scaled first-order jet has an exact two-variable Pruefer one-form.",
            exact["phase"]["phase_x"] + "; " + exact["phase"]["phase_t"],
            "Valid only where the scaled first jet is nonzero.",
            exact["phase"],
        ),
        ReductionRow(
            "nfosw_03_signed_ray_count",
            "exact_topological_lemma",
            "available_exact",
            "Closed proxy winding is an oriented signed-ray intersection count.",
            exact["signed_crossing"]["general"],
            "Tangencies require the stored regularization convention.",
            exact["signed_crossing"],
        ),
        ReductionRow(
            "nfosw_04_horizontal_crossings",
            "exact_crossing_reduction",
            "available_exact",
            "Every upward simple zero on a forward horizontal edge contributes -1 and contributes +1 when reversed.",
            exact["signed_crossing"]["horizontal_forward"]
            + " "
            + exact["signed_crossing"]["horizontal_reverse"],
            "Applies on derivative-compatible proxy pieces.",
        ),
        ReductionRow(
            "nfosw_05_vertical_crossings",
            "exact_crossing_reduction",
            "available_exact",
            "The right connector contributes an oriented J_t crossing count, not automatically zero.",
            exact["signed_crossing"]["vertical"],
            "Dominant-ray nonvanishing alone does not determine this open-path phase contribution.",
        ),
        ReductionRow(
            "nfosw_06_successor_chain",
            "exact_chain_identity",
            "available_exact",
            "The old rectangle, descendant slab, and new outer strip form the next rectangle with all internal edges cancelling.",
            exact["chain"]["oriented_identity"],
            "An exact oriented-chain identity.",
            exact["chain"],
        ),
        ReductionRow(
            "nfosw_07_base_degree",
            "exact_composition",
            "available_exact",
            "The ray-aligned exhaustion starts from a rigorously zero-degree rectangle.",
            exact["base"]["base_degree"],
            "Uses only the published Lambda<=1/5 upper input at t_0=1/4.",
            exact["base"],
        ),
        ReductionRow(
            "nfosw_08_outer_strip_degree",
            "exact_composition",
            "available_exact",
            "Every newly added outer strip has degree zero.",
            exact["outer_strip"]["degree"],
            "Uses the already checked tL>=25 dominant-saddle theorem.",
            exact["outer_strip"],
        ),
        ReductionRow(
            "nfosw_09_successor_recurrence",
            "exact_degree_reduction",
            "available_exact",
            "All possible degree growth is confined to the descendant slab.",
            exact["successor"]["recurrence"]
            + " "
            + exact["successor"]["positive_integer"],
            "No arithmetic exclusion inside the descendant slab is asserted.",
            exact["successor"],
        ),
        ReductionRow(
            "nfosw_10_one_sided_integer_trap",
            "exact_integer_reduction",
            "available_exact",
            "A strict upper bound below one is enough; exact zero winding is stronger than necessary.",
            exact["successor"]["one_sided_trap"],
            "Requires the boundary homotopy and the positive contact-degree theorem.",
        ),
        ReductionRow(
            "nfosw_11_proxy_transfer",
            "exact_conditional_composition",
            "available_exact",
            "The certified first-order remainder and cutoff ratios transfer the successor integer whenever the standard margin holds.",
            exact["proxy_transfer"]["conclusion"],
            "Conditional on the still-open boundary margin.",
            exact["proxy_transfer"],
        ),
        ReductionRow(
            "nfosw_12_frequency_margin",
            "open_theorem_target",
            "not_ready_to_apply",
            "Prove the q>=1 bottom boundary margin.",
            exact["open_arithmetic"]["frequency_margin"],
            "No Xi-specific signed logarithmic-moment inequality is supplied here.",
        ),
        ReductionRow(
            "nfosw_13_parabolic_and_shoulder_margin",
            "open_theorem_target",
            "not_ready_to_apply",
            "Prove the q<1 multiplicity-compatible margin and terminate every finite shoulder.",
            exact["open_arithmetic"]["parabolic_margin"]
            + " "
            + exact["open_arithmetic"]["shoulders"],
            "A fixed floor uniform to t=0 is forbidden.",
        ),
        ReductionRow(
            "nfosw_14_one_sided_phase_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "Upper-bound every composite successor winding by a number strictly below one.",
            exact["open_arithmetic"]["integer_target"],
            "Boundary nonvanishing and dominant-ray closure do not themselves give this integer bound.",
            exact["countermodel"],
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "exact oriented successor-chain, signed-crossing, and one-sided "
            "integer reduction; three Xi boundary obligations remain open"
        ),
        "proof_boundary": (
            "This artifact proves the orientation, Pruefer, signed-crossing, "
            "successor-chain, base-degree, outer-strip, positive-integer, and "
            "conditional proxy-transfer reductions. It does not prove either "
            "bottom margin, the finite shoulders, the one-sided phase bound, "
            "all-stage contact exclusion, Lambda<=0, RH, PF-infinity, or a "
            "Clay-prize conclusion."
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    lines = [
        "# Newman First-Order Oriented Successor-Winding Reduction",
        "",
        "Date: 2026-07-26",
        "",
        "Status: exact oriented successor and one-sided winding reduction.",
        "Three Xi boundary obligations remain open. This is not a proof of",
        "`Lambda<=0`, RH, or the Clay prize.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Orientation",
        "",
        "```text",
        exact["orientation"]["positive_boundary"],
        exact["orientation"]["contact_degree"],
        exact["orientation"]["reversal_guard"],
        "```",
        "",
        "The standard `(x,t)` orientation matches the boundary order used by",
        "the finite Q208 polygon and by every ray-aligned rectangle.",
        "",
        "## Scaled Pruefer Form",
        "",
        "```text",
        exact["phase"]["proxy"],
        exact["phase"]["norm"],
        exact["phase"]["phase_x"],
        exact["phase"]["phase_t"],
        exact["phase"]["simple_zero"],
        "```",
        "",
        "## Signed Crossings",
        "",
        "```text",
        exact["signed_crossing"]["general"],
        exact["signed_crossing"]["horizontal_forward"],
        exact["signed_crossing"]["horizontal_reverse"],
        exact["signed_crossing"]["vertical"],
        exact["signed_crossing"]["tangency_convention"],
        "```",
        "",
        "The vertical connector is part of the integer. Being contact-free",
        "does not by itself make its open-path phase contribution zero.",
        "",
        "## Successor Chain",
        "",
        "```text",
        exact["chain"]["rectangles"],
        exact["chain"]["oriented_identity"],
        exact["base"]["base_degree"],
        exact["outer_strip"]["lower_bound"],
        exact["outer_strip"]["degree"],
        exact["successor"]["recurrence"],
        exact["successor"]["positive_integer"],
        "```",
        "",
        "All internal boundary segments cancel in the stored exact chain audit.",
        "The initial rectangle `P_0` is already contact-free because its lower",
        "time is `1/4>1/5`; Q208 remains a compact calibration, not the logical",
        "source of the `P_0` degree.",
        "",
        "## One-Sided Integer Trap",
        "",
        "```text",
        exact["successor"]["one_sided_trap"],
        exact["successor"]["crossing_formula"],
        "```",
        "",
        "This is the useful relaxation: a rigorous oriented phase upper bound",
        "below one turn suffices. Computing exact zero winding is unnecessary.",
        "",
        "## First-Order Transfer",
        "",
        "```text",
        exact["proxy_transfer"]["remainder_budget"],
        exact["proxy_transfer"]["margin"],
        exact["proxy_transfer"]["cutoff_ratio"],
        exact["proxy_transfer"]["conclusion"],
        "```",
        "",
        "## Remaining Xi Input",
        "",
        "```text",
        exact["open_arithmetic"]["frequency_margin"],
        exact["open_arithmetic"]["parabolic_margin"],
        exact["open_arithmetic"]["shoulders"],
        exact["open_arithmetic"]["integer_target"],
        "```",
        "",
        "## Scope Guard",
        "",
        "```text",
        exact["countermodel"]["flow"],
        exact["countermodel"]["contact"],
        exact["countermodel"]["boundary"],
        exact["countermodel"]["guard"],
        "```",
        "",
        "## Live Handoff",
        "",
        "```text",
        exact["proof_handoff"],
        "```",
        "",
        "This reduction does not prove either bottom margin, the finite",
        "shoulders, the one-sided phase estimate, `Lambda<=0`, or RH.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman oriented successor-winding reduction: "
        "14 rows, 1 exact chain audit, 2 exact crossing reductions, "
        "1 one-sided integer trap, 3 open cofinal obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
