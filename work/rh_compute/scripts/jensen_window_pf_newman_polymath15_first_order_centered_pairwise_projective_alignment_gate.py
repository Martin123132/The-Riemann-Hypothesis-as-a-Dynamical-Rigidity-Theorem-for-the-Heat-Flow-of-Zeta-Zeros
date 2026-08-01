#!/usr/bin/env python3
"""Build the pairwise projective-alignment identity and route guard."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "pairwise_projective_alignment_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "interior_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "interior_projective_current_gate.json"
    ),
    "absolute_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_block_composition_guard.json"
    ),
    "contact_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contact_signed_transport_reduction.json"
    ),
    "legacy_symmetry": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "legacy_prime_curvature_symmetry_audit.json"
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
        "interior_current": (
            "J_n<=-(5/64)*u_n^2*|z_n|^2<0",
            "C_edge=c_0+c_N",
            "J_a^(der)=H_(a,x)+mu_aH_a",
        ),
        "absolute_phase": (
            "nu_n=partial_x arg(z_n)",
            "nu_n<=8449/x^2-u_n/2",
            "v_a+b*u_n",
        ),
        "contact_transport": (
            "sum_(j=0)^N d_j=mathcal_C_N",
            "cumulative-real-mass",
            "generally complex",
        ),
        "legacy_symmetry": (
            "transpose-symmetric",
            "raw and centered kernels are indefinite",
            "Xi contact scalar",
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


def alignment_audit() -> dict[str, str]:
    (
        b,
        c_s,
        c_s_x,
        b_x,
        u_x,
        aligned,
        u_n,
        u_m,
        k_n,
        k_m,
        u_terminal,
        e_n,
        e_m,
    ) = sp.symbols(
        "b c_s c_s_x b_x u_x H u_n u_m k_n k_m "
        "u_N e_n e_m",
        nonzero=True,
        real=True,
    )
    w_n = c_s * k_n - aligned
    w_m = c_s * k_m - aligned

    def aligned_jet(
        u_value: sp.Expr,
        k_value: sp.Expr,
        residual: sp.Expr,
        shifted: sp.Expr,
    ) -> sp.Expr:
        return sp.expand(
            -b**2 * u_value**2
            - shifted**2
            + c_s_x * k_value
            - (b_x / b + u_x / u_value) * shifted
            - b * u_value * residual
            - residual * shifted**2 / (b * u_value)
        )

    direct = sp.expand(
        aligned_jet(u_m, k_m, e_m, w_m)
        - aligned_jet(u_n, k_n, e_n, w_n)
    )
    delta = u_n - u_m
    sum_u = u_n + u_m
    a_n = b**2 * u_n**2 + w_n**2
    a_m = b**2 * u_m**2 + w_m**2
    compact = sp.expand(
        delta
        * (
            b**2 * sum_u
            + c_s * (c_s * (k_n + k_m) - 2 * aligned)
            - c_s_x
            + c_s * b_x / b
            + u_x
            * (c_s * u_terminal + aligned)
            / (u_n * u_m)
        )
        + e_n * a_n / (b * u_n)
        - e_m * a_m / (b * u_m)
    )
    relation = {
        k_n: u_n - u_terminal,
        k_m: u_m - u_terminal,
    }
    if sp.factor((direct - compact).subs(relation)) != 0:
        raise RuntimeError("compact pairwise alignment identity failed")

    return {
        "single_jet": (
            "At h_i=H put w_i=c_s*k_i-H and "
            "e_i=nu_i-b*u_i. Then "
            "h_(i,x)=-b^2*u_i^2-w_i^2+c_(s,x)*k_i"
            "-(b_x/b+u_x/u_i)w_i-b*u_i*e_i"
            "-e_i*w_i^2/(b*u_i)."
        ),
        "relative_identity": (
            "For 1<=n<m<=N-1 put delta=u_n-u_m=k_n-k_m>0, "
            "S=u_n+u_m, and A_i=b^2*u_i^2+w_i^2. At "
            "h_n=h_m=H, D_nm(H):=h_(m,x)-h_(n,x) equals "
            "delta*[b^2*S+c_s*(c_s*(k_n+k_m)-2H)-c_(s,x)"
            "+c_s*b_x/b+u_x*(c_s*u_N+H)/(u_n*u_m)]"
            "+e_n*A_n/(b*u_n)-e_m*A_m/(b*u_m)."
        ),
        "residual_definition": (
            "Here e_i=v_a+Im[d_(i,x)/(1+d_i)], so "
            "nu_i=b*u_i+e_i and |e_i|<8449/x^2."
        ),
    }


def moving_scale_audit() -> dict[str, str]:
    b, u_x, aligned, u_n, u_m = sp.symbols(
        "b u_x H u_n u_m",
        nonzero=True,
        real=True,
    )
    delta = u_n - u_m
    sum_u = u_n + u_m
    reduced = sp.factor(
        delta
        * (
            b**2 * sum_u
            + u_x * aligned / (u_n * u_m)
        )
    )
    critical = -b**2 * sum_u * u_n * u_m / u_x
    at_zero = sp.factor(reduced.subs(aligned, 0))
    below = sp.factor(reduced.subs(aligned, 2 * critical))
    if sp.simplify(at_zero - b**2 * delta * sum_u) != 0:
        raise RuntimeError("moving-scale positive guard failed")
    if sp.simplify(below + b**2 * delta * sum_u) != 0:
        raise RuntimeError("moving-scale negative guard failed")

    single = (
        -b**2 * u_n**2
        - aligned**2
        + u_x * aligned / u_n
    )
    if sp.simplify(single.subs(aligned, 0) + b**2 * u_n**2) != 0:
        raise RuntimeError("single-current zero-alignment guard failed")

    return {
        "frozen_model": (
            "If c_s=c_(s,x)=b_x=u_x=e_n=e_m=0, then "
            "D_nm=b^2*(u_n^2-u_m^2)=b^2*delta*S>0. "
            "This is the favorable frozen-scale correction-free "
            "Sturm orientation."
        ),
        "restored_scale": (
            "Restoring only the exact moving scale u_x>0 gives "
            "D_nm=delta*[b^2*S+u_x*H/(u_n*u_m)]."
        ),
        "critical_slope": (
            "Put H_mov=-b^2*S*u_n*u_m/u_x<0. Then "
            "D_nm(0)=b^2*delta*S>0, D_nm(H_mov)=0, and "
            "D_nm(2H_mov)=-b^2*delta*S<0."
        ),
        "individual_flow_guard": (
            "In that same reduced model "
            "h_(i,x)=-b^2*u_i^2-H^2+u_x*H/u_i. It is strictly "
            "negative at H=0 and at H=2H_mov<0 for both atoms. "
            "Thus strict individual clockwise motion and the actual "
            "logarithmic distance order do not fix pairwise crossing "
            "orientation once u_x is restored."
        ),
    }


def pole_audit() -> dict[str, str]:
    (
        b,
        c_s,
        c_s_x,
        b_x,
        u_x,
        aligned,
        u_n,
        u_m,
        k_n,
        k_m,
        u_terminal,
        e_n,
        e_m,
    ) = sp.symbols(
        "b c_s c_s_x b_x u_x H u_n u_m k_n k_m "
        "u_N e_n e_m",
        nonzero=True,
        real=True,
    )
    delta = u_n - u_m
    w_n = c_s * k_n - aligned
    w_m = c_s * k_m - aligned
    a_n = b**2 * u_n**2 + w_n**2
    a_m = b**2 * u_m**2 + w_m**2
    relative = (
        delta
        * (
            b**2 * (u_n + u_m)
            + c_s * (c_s * (k_n + k_m) - 2 * aligned)
            - c_s_x
            + c_s * b_x / b
            + u_x
            * (c_s * u_terminal + aligned)
            / (u_n * u_m)
        )
        + e_n * a_n / (b * u_n)
        - e_m * a_m / (b * u_m)
    )
    relation = {
        k_n: u_n - u_terminal,
        k_m: u_m - u_terminal,
    }
    limit_plus = sp.factor(
        sp.limit(relative.subs(relation) / (1 + aligned**2), aligned, sp.oo)
    )
    limit_minus = sp.factor(
        sp.limit(
            relative.subs(relation) / (1 + aligned**2),
            aligned,
            -sp.oo,
        )
    )
    direct_pole = sp.factor(
        (-(b * u_m + e_m) / (b * u_m))
        - (-(b * u_n + e_n) / (b * u_n))
    )
    if sp.factor(limit_plus - direct_pole) != 0:
        raise RuntimeError("positive pole-limit join failed")
    if sp.factor(limit_minus - direct_pole) != 0:
        raise RuntimeError("negative pole-limit join failed")

    return {
        "finite_rate": (
            "At a finite projective alignment, "
            "psi_(m,x)-psi_(n,x)=D_nm(H)/(1+H^2)."
        ),
        "exact_pole": (
            "At a common projection pole X_n=X_m=0, division-free "
            "currents give psi_(i,x)=-nu_i/(b*u_i)"
            "=-1-e_i/(b*u_i), and hence "
            "psi_(m,x)-psi_(n,x)="
            "e_n/(b*u_n)-e_m/(b*u_m)."
        ),
        "pole_limit": (
            "The two finite-chart limits H->+infinity and "
            "H->-infinity of D_nm(H)/(1+H^2) both equal the exact "
            "pole value e_n/(b*u_n)-e_m/(b*u_m)."
        ),
        "residual_order_guard": (
            "Choose any "
            "0<epsilon<min(8449/x^2,-b*delta,-b*u_m). "
            "The two hypothetical assignments "
            "(e_n,e_m)=(epsilon,0) and (0,epsilon) preserve "
            "nu_m-nu_n>0 and nu_n,nu_m<0. At the pole their "
            "relative rates are respectively epsilon/(b*u_n)<0 "
            "and -epsilon/(b*u_m)>0. Therefore the currently proved "
            "residual bounds and strict absolute angular-rate order "
            "do not determine the pole sign. These assignments are "
            "source-bound guards, not asserted Xi correction values."
        ),
    }


def occupation_audit() -> dict[str, str]:
    anchor = sp.symbols("A", real=True)
    c_1, c_2, c_3 = sp.symbols("c_1 c_2 c_3", real=True)
    h_1, h_2, h_3 = sp.symbols("h_1 h_2 h_3", real=True)
    mass = c_1 + c_2 + c_3
    moment = c_1 * h_1 + c_2 * h_2 + c_3 * h_3
    occupation = (
        c_1 * (h_1 - anchor)
        + c_2 * (h_2 - anchor)
        + c_3 * (h_3 - anchor)
    )
    if sp.expand(moment - anchor * mass - occupation) != 0:
        raise RuntimeError("occupation identity failed")

    two_direct = c_1 * h_1 + c_2 * h_2
    order_12 = h_2 * (c_1 + c_2) - (h_2 - h_1) * c_1
    order_21 = h_1 * (c_1 + c_2) - (h_1 - h_2) * c_2
    if sp.expand(two_direct - order_12) != 0:
        raise RuntimeError("first two-atom Abel identity failed")
    if sp.expand(two_direct - order_21) != 0:
        raise RuntimeError("second two-atom Abel identity failed")

    return {
        "occupation_identity": (
            "On a pole-free local chart let I be the interior atoms "
            "with c_i!=0, choose A<min_(i in I)h_i and "
            "B>max_(i in I)h_i, and put "
            "P_x(s)=sum_(i in I)c_i*1_(A<s<h_i). Then exactly "
            "sum_(i in I)d_i=A*sum_(i in I)c_i"
            "+integral_A^B P_x(s)ds."
        ),
        "swap_join": (
            "The threshold mass P_x(s) is permutation-free. At a "
            "pairwise alignment the two order-cell descriptions join "
            "without a singular term: "
            "h_2(c_1+c_2)-(h_2-h_1)c_1"
            "=h_1c_1+h_2c_2"
            "=h_1(c_1+c_2)-(h_1-h_2)c_2."
        ),
        "aggregate_decomposition": (
            "Let D_perp be the sum of d_i over interior atoms with "
            "c_i=0. With X_I=sum_(i in I)c_i, the endpoint-complete "
            "scalar is "
            "mathcal_C_N=D_edge+D_perp+A*X_I"
            "+integral_A^B P_x(s)ds, while "
            "mathsf_X=C_edge+X_I."
        ),
        "pole_chart_join": (
            "When an atom reaches c_i=0 it is retained in D_perp "
            "division-free and the next local projective chart is "
            "chosen before continuing. The exact pole-rate formula, "
            "not a tangent quotient, is the chart join."
        ),
    }


def build_exact() -> dict:
    alignment = alignment_audit()
    moving = moving_scale_audit()
    pole = pole_audit()
    occupation = occupation_audit()
    return {
        "pi_provenance": (
            "The pi in a^2=x/(4*pi)+t/16 and "
            "u_x=1/(8*pi*a^2) is inherited from the completed-zeta "
            "factor pi^(-s/2)*Gamma(s/2)*zeta(s) and its "
            "Riemann-Siegel saddle. It is not fitted from a circle, "
            "prime-curvature arc, polygon, or plot."
        ),
        "domain": (
            "Fix one N=floor(a) chart on L>=50, 0<tL<=25, and "
            "q=2*t*L^2>=1. For 1<=n<m<=N-1, put "
            "u_i=log(a/i), k_i=u_i-u_N=log(N/i), "
            "c_s=Re(s_*'), b=Im(s_*')<0, and "
            "h_i=c_s*k_i-b*u_i*tan(theta_i) whenever X_i!=0."
        ),
        "residual_definition": alignment["residual_definition"],
        "single_aligned_jet": alignment["single_jet"],
        "relative_alignment_identity": alignment["relative_identity"],
        "frozen_model_theorem": moving["frozen_model"],
        "moving_scale_restoration": moving["restored_scale"],
        "moving_scale_sign_reversal": (
            moving["critical_slope"] + " " + moving["individual_flow_guard"]
        ),
        "finite_projective_rate": pole["finite_rate"],
        "exact_pole_rate": pole["exact_pole"],
        "pole_limit_join": pole["pole_limit"],
        "residual_order_guard": pole["residual_order_guard"],
        "legacy_symmetry_boundary": (
            "The reconstructed prime-curvature matrix is exactly "
            "radial and transpose-symmetric by construction, while "
            "every tested raw symmetric kernel is indefinite. That "
            "static p<->q symmetry supplies neither the moving-scale "
            "term u_x*H/(u_n*u_m) nor an ordering of e_i/u_i. The "
            "pictures remain a useful falsification laboratory, but "
            "they cannot repair either pairwise obstruction without "
            "an exact identity to the Xi contact scalar."
        ),
        "occupation_identity": occupation["occupation_identity"],
        "swap_join": occupation["swap_join"],
        "aggregate_oriented_cell": occupation["aggregate_decomposition"],
        "pole_chart_join": occupation["pole_chart_join"],
        "exceptional_edge": (
            "The aggregate formula keeps "
            "(C_edge,D_edge)=(c_0+c_N,d_0+d_N) exact. At a cutoff, "
            "recompute terminal centering and transport this block "
            "through J_a^(adj); never replace the generally complex "
            "within-chart J_a^(der) by J_a^(adj), and never assign "
            "separate endpoint or terminal signs."
        ),
        "route_decision": (
            "Retire a uniform pairwise projective Sturm sign as a "
            "consequence of the currently proved current bounds. "
            "This does not prove that a stronger arithmetic Xi "
            "comparison is false. The next admissible target is a "
            "one-sided estimate for the signed occupation integral "
            "integral P_x(s)ds together with D_edge and D_perp, using "
            "actual Xi amplitudes and oscillation before absolute "
            "values. Pairwise alignments are harmless permutation "
            "joins; projection poles and cutoff edges remain exact "
            "separate joins."
        ),
        "q_lt_1": (
            "The q<1 layer remains a separate multiplicity-compatible "
            "parabolic/Hermite or boundary-degree theorem. No "
            "pairwise or occupation identity here proves t=0 endpoint "
            "simplicity."
        ),
        "proof_boundary": (
            "This proves the exact finite-alignment relative-current "
            "identity, the favorable frozen-scale model, the exact "
            "moving-scale sign-reversal guard, the division-free pole "
            "rate and its two-sided chart limit, the residual-order "
            "insufficiency guard, and a permutation-free local "
            "occupation identity. The guards are not actual Xi "
            "counterexamples and do not rule out a stronger "
            "source-specific theorem. This does not prove a signed Xi "
            "occupation bound, an edge-block sign, the Abel-scalar "
            "gap, successor count, q<1 closure, finite-height "
            "effectivity, contact exclusion, Lambda<=0, PF-infinity, "
            "RH, or a Clay-prize conclusion."
        ),
        "diagnostics": {
            "alignment": alignment,
            "moving_scale": moving,
            "pole": pole,
            "occupation": occupation,
        },
    }


def build_rows(exact: dict) -> list[GateRow]:
    rows = [
        (
            "pi_provenance",
            "definition_provenance",
            "certified",
            "The moving-scale constant pi remains source-traced.",
            "No plotted circle defines pi.",
        ),
        (
            "domain",
            "exact_domain",
            "certified",
            "The pair calculation is fixed-chart and interior.",
            "The terminal carrier and q<1 layer are excluded.",
        ),
        (
            "residual_definition",
            "exact_decomposition",
            "ready_to_apply",
            "Every phase correction is retained in e_i.",
            "No sign or ordering of e_i is assumed.",
        ),
        (
            "single_aligned_jet",
            "exact_identity",
            "ready_to_apply",
            "The single slope jet is exact at a finite alignment.",
            "It is not evaluated at a projection pole.",
        ),
        (
            "relative_alignment_identity",
            "exact_identity",
            "proved",
            "The actual pairwise alignment current has an exact compact form.",
            "All moving-frame and d_i corrections remain present.",
        ),
        (
            "frozen_model_theorem",
            "model_lemma",
            "proved",
            "The frozen correction-free model has favorable Sturm orientation.",
            "This is not the physical moving chart.",
        ),
        (
            "moving_scale_restoration",
            "exact_reduction",
            "ready_to_apply",
            "Restoring u_x introduces a term linear in aligned slope.",
            "No correction term is hidden.",
        ),
        (
            "moving_scale_sign_reversal",
            "nonpromotion_guard",
            "guard_validated",
            "Moving scale permits both pairwise crossing orientations.",
            "The reduced jets are not asserted Xi states.",
        ),
        (
            "finite_projective_rate",
            "exact_identity",
            "ready_to_apply",
            "Finite-slope relative projective speed has positive normalization.",
            "The formula remains chart-local.",
        ),
        (
            "exact_pole_rate",
            "division_free_join",
            "proved",
            "A common projection pole has an exact relative rate.",
            "No tangent ratio is used.",
        ),
        (
            "pole_limit_join",
            "division_free_join",
            "proved",
            "Both affine ends join the same projective pole rate.",
            "The limit does not supply a sign.",
        ),
        (
            "residual_order_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Existing residual bounds and angular order allow either pole sign.",
            "The hypothetical corrections are not actual Xi values.",
        ),
        (
            "legacy_symmetry_boundary",
            "nonpromotion_guard",
            "guard_validated",
            "Static prime-curvature symmetry does not sign dynamic pair currents.",
            "No image-to-Xi identity is asserted.",
        ),
        (
            "occupation_identity",
            "exact_identity",
            "proved",
            "The interior first moment is a signed occupation integral.",
            "It is local to a pole-free chart.",
        ),
        (
            "swap_join",
            "exact_cell_join",
            "proved",
            "Pairwise order swaps are permutation-free joins.",
            "This does not bound the signed occupation mass.",
        ),
        (
            "aggregate_oriented_cell",
            "exact_decomposition",
            "ready_to_apply",
            "The contact scalar splits into edge, pole, mass, and occupation terms.",
            "No term receives an unsupported sign.",
        ),
        (
            "pole_chart_join",
            "division_free_join",
            "ready_to_apply",
            "Projection poles are handled by chart change and vector sums.",
            "No global tangent ordering is asserted.",
        ),
        (
            "exceptional_edge",
            "exact_boundary_law",
            "ready_to_apply",
            "The terminal-endpoint block remains recurrent and exact.",
            "The two historical J_a objects remain distinct.",
        ),
        (
            "route_decision",
            "open_theorem_target",
            "not_ready_to_apply",
            "The surviving target is a signed Xi occupation estimate.",
            "No aggregate gap has been proved.",
        ),
        (
            "q_lt_1",
            "open_theorem_target",
            "not_ready_to_apply",
            "The small-q multiplicity theorem remains separate.",
            "No t=0 simplicity is assumed.",
        ),
        (
            "proof_boundary",
            "proof_guard",
            "guard_validated",
            "The pairwise route gate is not promoted to contact exclusion.",
            "RH and every prize-level conclusion remain open.",
        ),
    ]
    return [
        GateRow(
            id=f"ppag_{index:02d}_{suffix}",
            role=role,
            readiness=readiness,
            claim=claim,
            formula=exact[suffix],
            proof_boundary=boundary,
            diagnostics=(
                exact["diagnostics"]
                if suffix
                in {
                    "relative_alignment_identity",
                    "moving_scale_sign_reversal",
                    "pole_limit_join",
                    "residual_order_guard",
                    "occupation_identity",
                }
                else None
            ),
        )
        for index, (
            suffix,
            role,
            readiness,
            claim,
            boundary,
        ) in enumerate(rows)
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact pairwise projective-alignment identity proved; "
            "frozen model favorable; moving-scale and pole residual "
            "guards block promotion from current bounds; signed Xi "
            "occupation estimate remains open"
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
            "# Pairwise Projective-Alignment Gate",
            "",
            "Date: 2026-07-28",
            "",
            "Status: exact alignment and pole identities with two route guards;",
            "not a proof of the signed Xi occupation bound, contact exclusion,",
            "`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.",
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
            "## Exact Alignment",
            "",
            exact["domain"],
            "",
            exact["residual_definition"],
            "",
            "```text",
            exact["single_aligned_jet"],
            exact["relative_alignment_identity"],
            "```",
            "",
            "## Frozen Model",
            "",
            exact["frozen_model_theorem"],
            "",
            "## Moving-Scale Obstruction",
            "",
            "```text",
            exact["moving_scale_restoration"],
            exact["moving_scale_sign_reversal"],
            "```",
            "",
            "## Projection Pole",
            "",
            exact["finite_projective_rate"],
            "",
            "```text",
            exact["exact_pole_rate"],
            exact["pole_limit_join"],
            "```",
            "",
            exact["residual_order_guard"],
            "",
            "## Legacy Symmetry Boundary",
            "",
            exact["legacy_symmetry_boundary"],
            "",
            "## Aggregate Occupation Identity",
            "",
            "```text",
            exact["occupation_identity"],
            exact["aggregate_oriented_cell"],
            "```",
            "",
            exact["swap_join"],
            "",
            exact["pole_chart_join"],
            "",
            exact["exceptional_edge"],
            "",
            "## Route Decision",
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
        "built pairwise projective-alignment gate: "
        "21 rows, exact relative current and pole join, frozen-model "
        "Sturm term, moving-scale and residual-order guards, aggregate "
        "occupation identity"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
