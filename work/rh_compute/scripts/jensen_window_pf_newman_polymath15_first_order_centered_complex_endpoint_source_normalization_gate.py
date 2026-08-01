#!/usr/bin/env python3
"""Build the complex recurrent-endpoint source-normalization corrigendum."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_complex_endpoint_source_normalization_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "c0_source_extraction": (
        REPO_ROOT
        / "work/rh_compute/scripts/"
        "jensen_window_pf_newman_polymath15_"
        "critical_RS_C1_endpoint_peeling_contract.py"
    ),
    "absolute_phase_anchor_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_absolute_phase_anchor_reduction.json"
    ),
    "direct_projection_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_direct_projection_regime_reduction.json"
    ),
    "carrier_kernel_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_carrier_kernel_abel_prefix_reduction.json"
    ),
    "first_pivot_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_endpoint_first_pivot_odd_small_ball_guard.json"
    ),
    "odd_fibre_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_endpoint_odd_fibre_"
        "correlation_feasibility_guard.json"
    ),
    "contact_transport_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_contact_signed_transport_reduction.json"
    ),
    "interior_current_historical": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_interior_projective_current_gate.json"
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


def audit_sources() -> dict:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "published_source": {
            "citation": (
                "D.H.J. Polymath, Effective approximation of heat flow "
                "evolution of the Riemann xi function, and a new upper "
                "bound for the de Bruijn-Newman constant, "
                "arXiv:1904.12438"
            ),
            "location": "equation (53), C_0(p)",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def verify_symbolics() -> None:
    p = sp.symbols("p", real=True)
    pi = sp.pi
    c0 = (
        sp.exp(sp.I * pi * (p**2 / 2 + sp.Rational(3, 8)))
        - sp.I * sp.sqrt(2) * sp.cos(pi * p / 2)
    ) / (2 * sp.cos(pi * p))
    expected_real = sp.sqrt(2 - sp.sqrt(2)) / 4
    expected_imag = (
        sp.sqrt(2 + sp.sqrt(2)) / 4 - sp.sqrt(2) / 2
    )
    midpoint = sp.expand_complex(c0.subs(p, 0))
    if sp.simplify(sp.re(midpoint) - expected_real) != 0:
        raise RuntimeError("C_0(0) real part failed")
    if sp.simplify(sp.im(midpoint) - expected_imag) != 0:
        raise RuntimeError("C_0(0) imaginary part failed")
    if sp.simplify(sp.diff(c0, p, 3).subs(p, 0)) != 0:
        raise RuntimeError("C_0'''(0) parity failed")

    kappa, t0, hr, hi, jr, ji, c, b, u = sp.symbols(
        "kappa T_0 H_R H_I J_R J_I c b u_N", real=True
    )
    h = hr + sp.I * hi
    j = jr + sp.I * ji
    e = sp.expand(kappa * (t0 + sp.I) * h)
    g = sp.expand(kappa * (t0 + sp.I) * j)
    expected = {
        sp.re(e): kappa * (t0 * hr - hi),
        sp.im(e): kappa * (t0 * hi + hr),
        sp.re(g): kappa * (t0 * jr - ji),
        sp.im(g): kappa * (t0 * ji + jr),
    }
    for actual, target in expected.items():
        if sp.simplify(actual - target) != 0:
            raise RuntimeError("complex endpoint projection failed")

    hermitian = sp.expand_complex(g * sp.conjugate(e))
    target_hermitian = sp.expand_complex(
        kappa**2 * (t0**2 + 1) * j * sp.conjugate(h)
    )
    if sp.simplify(hermitian - target_hermitian) != 0:
        raise RuntimeError("endpoint Hermitian product failed")

    x, y, un = sp.symbols("X_n Y_n u_n", real=True)
    c0_atom = kappa * (t0 * hr - hi)
    d0_atom = kappa * (
        t0 * jr - ji - c * u * (t0 * hr - hi)
    )
    dn_atom = c * (un - u) * x - b * un * y
    minor = sp.expand(c0_atom * dn_atom - d0_atom * x)
    target_minor = sp.expand(
        kappa
        * (
            (t0 * hr - hi) * un * (c * x - b * y)
            - (t0 * jr - ji) * x
        )
    )
    if sp.simplify(minor - target_minor) != 0:
        raise RuntimeError("endpoint-carrier minor failed")

    sigma = sp.Matrix([1, -1, -1])
    kernel = sigma * sigma.T
    if kernel.rank() != 1 or kernel.eigenvals() != {sp.Integer(3): 1, sp.Integer(0): 2}:
        raise RuntimeError("Vaughan signed kernel rank failed")


def exact_payload() -> dict:
    return {
        "source_normalization": {
            "c0_formula": (
                "C_0(p)={exp(pi*i*(p^2/2+3/8))"
                "-i*sqrt(2)*cos(pi*p/2)}/{2*cos(pi*p)}"
            ),
            "critical_endpoint": (
                "H_a(p)=C_0(p)+C_0'''(p)/(12*pi^2*a)"
            ),
            "notation_guard": (
                "The endpoint function C_0(p)=F(p) is not the Vaughan "
                "low coefficient C_0(n) from the preceding decomposition."
            ),
            "source_semantics": (
                "C_0 and H_a are generally complex on the real p-axis; "
                "the source gives |C_0(p)|<=1/2, not C_0(p) in R."
            ),
            "pi_provenance": (
                "Every pi in this correction is inherited from the "
                "published completed-zeta/Riemann-Siegel endpoint "
                "normalization. The Cartesian repair introduces no new pi."
            ),
        },
        "midpoint_witness": {
            "parity": "C_0(-p)=C_0(p), hence C_0'''(0)=0",
            "value": (
                "C_0(0)=sqrt(2-sqrt(2))/4"
                "+i{sqrt(2+sqrt(2))/4-sqrt(2)/2}"
            ),
            "strict_sign": (
                "Im C_0(0)<0 because sqrt(2+sqrt(2))<2*sqrt(2)"
            ),
            "endpoint_consequence": (
                "For every a>0, H_a(0)=C_0(0) is nonreal."
            ),
            "physical_occurrence": (
                "p=0 occurs at every half-integer saddle a=N+1/2, "
                "including arbitrarily large cutoffs."
            ),
        },
        "complex_endpoint": {
            "coordinates": (
                "kappa in R, H_a=H_R+iH_I, "
                "J_a=J_R+iJ_I, e=kappa(T_0+i)H_a, "
                "g=kappa(T_0+i)J_a"
            ),
            "value_projection": (
                "Re e=kappa(T_0H_R-H_I), "
                "Im e=kappa(T_0H_I+H_R)"
            ),
            "jet_projection": (
                "Re g=kappa(T_0J_R-J_I), "
                "Im g=kappa(T_0J_I+J_R)"
            ),
            "hermitian_product": (
                "g*conj(e)=kappa^2(T_0^2+1)J_a*conj(H_a)"
            ),
            "factorization_survives": (
                "E_N=kappa(T_0+i)[J_a-s_*'log(a)H_a]"
            ),
        },
        "centered_endpoint_atom": {
            "definition": (
                "alpha=c*u_N, c_0=Re e, d_0=Re g-alpha*Re e"
            ),
            "value": "c_0=kappa(T_0H_R-H_I)",
            "slope": (
                "d_0=kappa{T_0J_R-J_I"
                "-c*u_N(T_0H_R-H_I)}"
            ),
            "effective_slope": (
                "h_0=(T_0J_R-J_I)/(T_0H_R-H_I)-c*u_N "
                "when T_0H_R-H_I!=0"
            ),
            "exceptional_fibre": (
                "The division-free endpoint fibre is "
                "T_0H_R-H_I=0; H_a=0 is only a subfibre."
            ),
        },
        "odd_fibre_repair": {
            "coordinates": (
                "For the physical endpoint g_0=A_N+iB_N and "
                "S_odd=X_odd+iY_odd, "
                "Re(S_odd*conj(g_0))=A_NX_odd+B_NY_odd."
            ),
            "cartesian": (
                "|S_odd+g_0|^2=(X_odd+A_N)^2+(Y_odd+B_N)^2"
            ),
            "rotation": (
                "For D=sqrt(A_N^2+B_N^2)>0, "
                "P=(A_NX_odd+B_NY_odd)/D and "
                "Q=(-B_NX_odd+A_NY_odd)/D give "
                "|S_odd+g_0|^2=(P+D)^2+Q^2."
            ),
            "zero_endpoint": (
                "For D=0 retain |S_odd|^2 without division."
            ),
        },
        "contact_minors": {
            "carrier_atom": (
                "For z_n=X_n+iY_n, "
                "d_n=c(u_n-u_N)X_n-b*u_nY_n."
            ),
            "endpoint_carrier": (
                "K_(0n)=c_0d_n-d_0X_n="
                "kappa{(T_0H_R-H_I)u_n(cX_n-bY_n)"
                "-(T_0J_R-J_I)X_n}."
            ),
            "carrier_carrier": (
                "K_(nm)=X_nd_m-d_nX_m="
                "c*log(n/m)X_nX_m"
                "-b{u_mX_nY_m-u_nY_nX_m}."
            ),
            "sign_boundary": (
                "Neither determinant has a sign from the currently "
                "proved amplitude, distance, or source-normalization data."
            ),
        },
        "rank_guard": {
            "vaughan_scalar": (
                "For r=(R_0,R_I,R_II) and sigma=(1,-1,-1), "
                "Re mathfrak B_N=sigma*r."
            ),
            "vaughan_kernel": (
                "sigma*sigma^T=[[1,-1,-1],[-1,1,1],[-1,1,1]] "
                "has rank 1 and eigenvalues 3,0,0."
            ),
            "pullback": (
                "After pullback to the endpoint/carrier atoms the square "
                "is d*d^T, the original centered-slope rank-one kernel."
            ),
            "contact_join": (
                "Adding the value row gives c*c^T+d*d^T of rank at most "
                "2; on c*x=0 only the rank-at-most-one slope square remains."
            ),
            "consequence": (
                "Exact Vaughan recombination changes coordinates but cannot "
                "create a contact margin. A new Xi restriction on the "
                "actual coefficient curve is still required."
            ),
        },
        "supersession": {
            "scope": (
                "Only endpoint rows that used H_a in R are superseded. "
                "The branch-free anchor, complex factorization, abstract "
                "component transport, Abel/Mangoldt identities, adjacent "
                "cutoff recurrence, and endpoint-composed Vaughan identity "
                "remain exact."
            ),
            "historical_artifacts": [
                "absolute_phase_anchor_reduction: real-H declaration only",
                "direct_projection_regime_reduction: endpoint projections",
                "carrier_kernel_abel_prefix_reduction: endpoint Hermitian product",
                "endpoint_first_pivot_odd_small_ball_guard: real-H sentence only",
                "endpoint_odd_fibre_correlation_feasibility_guard: fixed endpoint direction",
                "contact_signed_transport_reduction: endpoint atom and exceptional fibre",
                "interior_projective_current_gate: endpoint effective-slope text",
            ],
            "route_decision": (
                "Use the corrected Cartesian endpoint atom in every later "
                "contact or Type-I/II argument. Do not infer a fixed "
                "(T_0,1) endpoint direction, and do not divide merely "
                "under H_a!=0."
            ),
        },
    }


def build_payload() -> dict:
    verify_symbolics()
    exact = exact_payload()
    source_audit = audit_sources()
    rows = [
        GateRow(
            id="cesn_01_published_c0",
            role="primary_source_correction",
            readiness="available_published",
            claim="The published C_0 endpoint coefficient is complex-valued.",
            formula=exact["source_normalization"]["c0_formula"],
            proof_boundary="This is source normalization, not a new estimate.",
        ),
        GateRow(
            id="cesn_02_midpoint_parity",
            role="exact_identity",
            readiness="available_exact",
            claim="Evenness kills the third derivative at the midpoint.",
            formula=exact["midpoint_witness"]["parity"],
            proof_boundary="This statement is only at p=0.",
        ),
        GateRow(
            id="cesn_03_midpoint_nonreal",
            role="exact_counterexample",
            readiness="available_exact",
            claim="The actual endpoint source is nonreal at p=0.",
            formula=exact["midpoint_witness"]["value"],
            proof_boundary="This rejects H_a in R; it is not an RH counterexample.",
        ),
        GateRow(
            id="cesn_04_endpoint_occurrence",
            role="domain_audit",
            readiness="available_exact",
            claim="The nonreal midpoint recurs at unbounded physical cutoffs.",
            formula=exact["midpoint_witness"]["physical_occurrence"],
            proof_boundary="No lower bound for the endpoint projection follows.",
        ),
        GateRow(
            id="cesn_05_endpoint_factorization",
            role="surviving_identity",
            readiness="available_exact",
            claim="The complex endpoint keeps the existing common factorization.",
            formula=exact["complex_endpoint"]["coordinates"],
            proof_boundary="Factorization alone gives no sign.",
        ),
        GateRow(
            id="cesn_06_value_projection",
            role="corrected_identity",
            readiness="ready_to_apply",
            claim="Both Cartesian value projections retain H_I.",
            formula=exact["complex_endpoint"]["value_projection"],
            proof_boundary="No projection nonvanishing is asserted.",
        ),
        GateRow(
            id="cesn_07_jet_projection",
            role="corrected_identity",
            readiness="ready_to_apply",
            claim="Both Cartesian jet projections retain J_I.",
            formula=exact["complex_endpoint"]["jet_projection"],
            proof_boundary="No directional transversality is asserted.",
        ),
        GateRow(
            id="cesn_08_hermitian_product",
            role="corrected_identity",
            readiness="ready_to_apply",
            claim="The endpoint Hermitian product uses conjugate(H_a).",
            formula=exact["complex_endpoint"]["hermitian_product"],
            proof_boundary="Its real and imaginary parts remain unsigned.",
        ),
        GateRow(
            id="cesn_09_centered_atom",
            role="corrected_identity",
            readiness="ready_to_apply",
            claim="The terminal-centered endpoint atom has the full complex projection.",
            formula=(
                exact["centered_endpoint_atom"]["value"]
                + "; "
                + exact["centered_endpoint_atom"]["slope"]
            ),
            proof_boundary="This is an exact atom, not a contact gap.",
        ),
        GateRow(
            id="cesn_10_exceptional_fibre",
            role="division_guard",
            readiness="ready_to_apply",
            claim="The projective exceptional set is the zero real projection.",
            formula=exact["centered_endpoint_atom"]["exceptional_fibre"],
            proof_boundary="H_a=0 is not the complete exceptional set.",
        ),
        GateRow(
            id="cesn_11_odd_fibre_rotation",
            role="corrected_identity",
            readiness="ready_to_apply",
            claim="The odd-fibre rotation follows the actual complex endpoint vector.",
            formula=exact["odd_fibre_repair"]["rotation"],
            proof_boundary="This coordinate change gives no small-ball estimate.",
        ),
        GateRow(
            id="cesn_12_endpoint_carrier_minor",
            role="joint_source_kernel",
            readiness="available_exact",
            claim="The endpoint-carrier contact minor is explicit and unsigned.",
            formula=exact["contact_minors"]["endpoint_carrier"],
            proof_boundary="No Xi phase law currently signs this minor.",
        ),
        GateRow(
            id="cesn_13_carrier_carrier_minor",
            role="joint_source_kernel",
            readiness="available_exact",
            claim="The carrier-carrier contact minor is explicit and unsigned.",
            formula=exact["contact_minors"]["carrier_carrier"],
            proof_boundary="No aggregate cumulative-mass estimate follows.",
        ),
        GateRow(
            id="cesn_14_vaughan_rank",
            role="nonpromotion_guard",
            readiness="guard_validated",
            claim="The exact three-block Vaughan square has rank one.",
            formula=exact["rank_guard"]["vaughan_kernel"],
            proof_boundary="Separate block norms cannot create the missing sign.",
        ),
        GateRow(
            id="cesn_15_contact_rank",
            role="nonpromotion_guard",
            readiness="guard_validated",
            claim="The joined value/slope Gram has rank at most two.",
            formula=exact["rank_guard"]["contact_join"],
            proof_boundary="Positive semidefiniteness alone gives no contact margin.",
        ),
        GateRow(
            id="cesn_16_supersession",
            role="proof_corpus_corrigendum",
            readiness="ready_to_apply",
            claim="Only the real-H endpoint specializations are superseded.",
            formula=exact["supersession"]["scope"],
            proof_boundary="All surviving identities retain their existing proof boundaries.",
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact complex recurrent-endpoint source-normalization "
            "corrigendum and rank nonpromotion gate"
        ),
        "proof_boundary": (
            "This artifact corrects the recurrent endpoint's Cartesian "
            "projections, Hermitian product, centered atom, exceptional "
            "fibre, odd-fibre rotation, and contact minors. It proves that "
            "the exact Vaughan split has no new rank capable of supplying "
            "the missing contact margin. It does not prove an Xi phase law, "
            "signed endpoint/carrier correlation estimate, Abel-scalar gap, "
            "successor winding cap, contact exclusion, Q209, the cofinal "
            "descendant theorem, Lambda<=0, PF-infinity, RH, or a "
            "prize-level conclusion."
        ),
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "source_nonreal_witnesses": 1,
            "corrected_cartesian_projections": 4,
            "corrected_hermitian_products": 1,
            "corrected_centered_endpoint_atoms": 1,
            "corrected_odd_fibre_rotations": 1,
            "exact_contact_minors": 2,
            "rank_guards": 2,
            "historical_artifacts_quarantined": 7,
            "signed_lower_bounds": 0,
            "abel_gaps": 0,
            "winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    source = payload["source_audit"]["published_source"]
    historical = "\n".join(
        f"- `{item}`"
        for item in exact["supersession"]["historical_artifacts"]
    )
    return "\n".join(
        [
            "# Complex Endpoint Source-Normalization Corrigendum",
            "",
            "Date: 2026-07-31",
            "",
            "Status: exact correction gate and nonpromotion audit. This is",
            "not a proof of the Xi contact gap, `Lambda<=0`, or RH.",
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
            "validated complex endpoint source-normalization gate: "
            "16 rows, 1 nonreal source witness, "
            "4 corrected Cartesian projections, 2 contact minors, "
            "2 rank guards, 7 historical artifacts quarantined, "
            "0 signed lower bounds, 0 Abel gaps, 0 winding bounds",
            "```",
            "",
            "## Source Correction",
            "",
            f"Published source: {source['citation']}, {source['location']}.",
            "",
            source["url"],
            "",
            "```text",
            exact["source_normalization"]["c0_formula"],
            exact["source_normalization"]["critical_endpoint"],
            "```",
            "",
            exact["source_normalization"]["notation_guard"],
            "",
            exact["source_normalization"]["source_semantics"],
            "",
            "At the exact midpoint,",
            "",
            "```text",
            exact["midpoint_witness"]["parity"],
            exact["midpoint_witness"]["value"],
            exact["midpoint_witness"]["strict_sign"],
            exact["midpoint_witness"]["endpoint_consequence"],
            "```",
            "",
            exact["midpoint_witness"]["physical_occurrence"],
            "",
            "## Correct Endpoint",
            "",
            "```text",
            exact["complex_endpoint"]["coordinates"],
            exact["complex_endpoint"]["value_projection"],
            exact["complex_endpoint"]["jet_projection"],
            exact["complex_endpoint"]["hermitian_product"],
            exact["complex_endpoint"]["factorization_survives"],
            "```",
            "",
            "The centered endpoint atom is",
            "",
            "```text",
            exact["centered_endpoint_atom"]["definition"],
            exact["centered_endpoint_atom"]["value"],
            exact["centered_endpoint_atom"]["slope"],
            exact["centered_endpoint_atom"]["effective_slope"],
            exact["centered_endpoint_atom"]["exceptional_fibre"],
            "```",
            "",
            "## Odd Fibre",
            "",
            "```text",
            exact["odd_fibre_repair"]["coordinates"],
            exact["odd_fibre_repair"]["cartesian"],
            exact["odd_fibre_repair"]["rotation"],
            exact["odd_fibre_repair"]["zero_endpoint"],
            "```",
            "",
            "The rotation follows the actual complex endpoint vector. A",
            "fixed `(T_0,1)` endpoint direction is not source-derived.",
            "",
            "## Contact Minors",
            "",
            "```text",
            exact["contact_minors"]["carrier_atom"],
            exact["contact_minors"]["endpoint_carrier"],
            exact["contact_minors"]["carrier_carrier"],
            "```",
            "",
            exact["contact_minors"]["sign_boundary"],
            "",
            "## Vaughan Rank",
            "",
            "```text",
            exact["rank_guard"]["vaughan_scalar"],
            exact["rank_guard"]["vaughan_kernel"],
            exact["rank_guard"]["pullback"],
            exact["rank_guard"]["contact_join"],
            "```",
            "",
            exact["rank_guard"]["consequence"],
            "",
            "## Supersession",
            "",
            exact["supersession"]["scope"],
            "",
            historical,
            "",
            exact["supersession"]["route_decision"],
            "",
            "## Pi Provenance",
            "",
            exact["source_normalization"]["pi_provenance"],
            "",
            "## Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built complex endpoint source-normalization gate: "
        "16 rows, 1 nonreal source witness, "
        "4 corrected Cartesian projections, 2 contact minors, "
        "2 rank guards, 7 historical artifacts quarantined, "
        "0 signed lower bounds, 0 Abel gaps, 0 winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
