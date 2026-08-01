#!/usr/bin/env python3
"""Build the Newman first-jet winding and boundary-flux gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_first_jet_winding_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
BOUNDARY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)
DIAGONAL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json"
)
COMPACT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_compact_transversality_"
    "interval_certificate.json"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def heat_polynomial(m: int, tau: sp.Symbol, x: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sum(
            (-tau) ** k
            * sp.diff(x**m, x, 2 * k)
            / sp.factorial(k)
            for k in range(m // 2 + 1)
        )
    )


def hermite_audits() -> list[dict]:
    x, tau = sp.symbols("x tau", real=True)
    rows: list[dict] = []
    for m in range(2, 11):
        polynomial = heat_polynomial(m, tau, x)
        if sp.simplify(
            sp.diff(polynomial, tau) + sp.diff(polynomial, x, 2)
        ) != 0:
            raise RuntimeError(f"heat identity failed at multiplicity {m}")
        degree = m * (m - 1) // 2
        expected_discriminant = (
            (2 * tau) ** degree
            * sp.prod(k**k for k in range(1, m + 1))
        )
        discriminant = sp.factor(sp.discriminant(polynomial, x))
        if sp.simplify(discriminant - expected_discriminant) != 0:
            raise RuntimeError(
                f"Hermite discriminant failed at multiplicity {m}"
            )
        roots_before = int(
            sp.Poly(polynomial.subs(tau, -1), x).count_roots(
                -sp.oo, sp.oo
            )
        )
        roots_after = int(
            sp.Poly(polynomial.subs(tau, 1), x).count_roots(
                -sp.oo, sp.oo
            )
        )
        expected_before = m % 2
        if roots_before != expected_before or roots_after != m:
            raise RuntimeError(
                f"root-count jump failed at multiplicity {m}"
            )
        rows.append(
            {
                "multiplicity": m,
                "polynomial": str(polynomial),
                "discriminant_tau_power": degree,
                "discriminant_constant": str(
                    sp.prod(k**k for k in range(1, m + 1))
                    * 2**degree
                ),
                "real_roots_tau_negative": roots_before,
                "real_roots_tau_positive": roots_after,
                "pair_births": (roots_after - roots_before) // 2,
                "local_index": (roots_after - roots_before) // 2,
            }
        )
    return rows


def symbolic_flux_audit() -> dict:
    f, fx, fxx, fxxx = sp.symbols(
        "f fx fxx fxxx", real=True
    )
    q = f**2 + fx**2
    laguerre = fx**2 - f * fxx
    phase_x = sp.simplify((f * fxx - fx**2) / q)
    phase_t = sp.simplify(
        (f * (-fxxx) - fx * (-fxx)) / q
    )
    laguerre_x = fx * fxx - f * fxxx
    if sp.simplify(phase_x + laguerre / q) != 0:
        raise RuntimeError("horizontal phase-flux identity failed")
    if sp.simplify(phase_t - laguerre_x / q) != 0:
        raise RuntimeError("vertical phase-flux identity failed")
    return {
        "Q": "Q=H^2+H_x^2",
        "L": "L=H_x^2-H*H_xx",
        "horizontal": "partial_x arg(H+iH_x)=-L/Q",
        "vertical": "partial_t arg(H+iH_x)=L_x/Q",
        "one_form": "d arg(H+iH_x)=-(L/Q)dx+(L_x/Q)dt",
    }


def symbolic_orientation_audit() -> dict:
    f_x, f_xx, f_xt = sp.symbols("f_x f_xx f_xt", real=True)
    f_t = -f_xx
    jacobian_x_t = sp.Matrix(
        [[f_x, f_t], [f_xx, f_xt]]
    ).det().subs(f_x, 0)
    jacobian_t_x = sp.Matrix(
        [[f_t, f_x], [f_xt, f_xx]]
    ).det().subs(f_x, 0)
    if sp.simplify(jacobian_x_t - f_xx**2) != 0:
        raise RuntimeError("standard (x,t) contact orientation failed")
    if sp.simplify(jacobian_t_x + f_xx**2) != 0:
        raise RuntimeError("reversed (t,x) contact orientation failed")
    if sp.simplify(jacobian_x_t + jacobian_t_x) != 0:
        raise RuntimeError("coordinate reversal did not reverse contact degree")
    return {
        "domain_orientation": "(x,t)",
        "positive_bottom_edge": "t=delta, x:0->R",
        "determinant_x_t": str(jacobian_x_t),
        "determinant_t_x": str(jacobian_t_x),
        "coordinate_reversal": (
            "Changing (x,t) to (t,x) reverses both local degree and the "
            "positive boundary orientation."
        ),
    }


def build_exact() -> dict:
    return {
        "orientation_convention": (
            "The domain has the standard (x,t) orientation: the positively "
            "oriented half-rectangle runs along its bottom edge from x=0 to "
            "x=R. Reversing the coordinate order to (t,x) reverses every "
            "local index and the boundary winding."
        ),
        "universal_heat_polynomial": (
            "P_m(tau,u)=exp(-tau*D_u^2)u^m="
            "(2tau)^(m/2)He_m(u/sqrt(2tau)) for tau>0"
        ),
        "discriminant": (
            "Disc_u P_m=(2tau)^(m(m-1)/2)"
            "*product_(k=1)^m k^k"
        ),
        "root_count_jump": (
            "For tau>0, P_m has m simple real Hermite roots. For tau<0, "
            "its coefficients have one sign in u^2, so it has zero real "
            "roots when m is even and only u=0 when m is odd. The contact "
            "therefore births floor(m/2) real pairs."
        ),
        "double_index": (
            "At a nondegenerate common zero of (F,F_x), "
            "det D_(x,t)(F,F_x)=F_x*F_xt-F_t*F_xx=F_xx^2>0; "
            "in the reversed (t,x) coordinate order the determinant is "
            "-F_xx^2<0"
        ),
        "multiplicity_index": (
            "At a real zero of spatial multiplicity m>=2 in any nonzero "
            "real-analytic solution F_t=-F_xx, the common-zero point is "
            "isolated and its local Brouwer index in the standard (x,t) "
            "orientation is +floor(m/2)"
        ),
        "global_winding": (
            "If Omega is a compact planar domain whose boundary contains no "
            "common zero, wind((F+iF_x)(partial Omega),0)="
            "+sum_(p in Omega) floor(m_p/2)>=0 in the standard (x,t) "
            "orientation"
        ),
        "zero_winding_equivalence": (
            "Under boundary nonvanishing, wind(F+iF_x)=0 iff Omega contains "
            "no real multiple-zero contact"
        ),
        "orientation_audit": symbolic_orientation_audit(),
        "flux": symbolic_flux_audit(),
        "xi_half_rectangles": (
            "Let Q_j=[1/(5j),1/4]x[0,38+j]. Then Lambda<=0 iff, for every "
            "j, Z=H+iH_x is nonzero on partial Q_j and "
            "wind(Z(partial Q_j),0)=0"
        ),
        "known_edges": (
            "On x=0, H_t(0)>0. On t=1/4, all zeros are simple because "
            "1/4>1/5>=Lambda. Thus only the bottom edge "
            "t=1/(5j), 0<=x<=38+j and the right edge "
            "x=38+j, 1/(5j)<=t<=1/4 need new boundary separation"
        ),
        "half_rectangle_flux": (
            "For Q=[delta,T]x[0,R], counterclockwise orientation gives "
            "2pi*wind(Z)=-integral_0^R L_delta/Q_delta dx"
            "+integral_delta^T L_x(t,R)/Q(t,R)dt"
            "+integral_0^R L_T/Q_T dx"
        ),
        "countermodel_guard": (
            "The universal m=2 and m=4 heat contacts have winding +1 and "
            "+2 respectively in the standard (x,t) orientation. Top-edge "
            "simplicity, axis positivity, or "
            "forward real-rootedness alone does not force zero winding."
        ),
        "proof_handoff": (
            "For every diagonal half-rectangle, prove first-jet separation "
            "on the bottom and right edges and prove that the displayed "
            "phase-flux integer is zero. This is a one-dimensional "
            "boundary theorem with sign-definite interior charges; it is "
            "not yet supplied by the compact certificate, the dominant "
            "ray, root-field balance, or generic heat-flow topology."
        ),
    }


def source_audit() -> dict:
    sources = {
        "boundary": json.loads(BOUNDARY_SOURCE.read_text(encoding="utf-8")),
        "diagonal": json.loads(DIAGONAL_SOURCE.read_text(encoding="utf-8")),
        "compact": json.loads(COMPACT_SOURCE.read_text(encoding="utf-8")),
    }
    texts = {key: json.dumps(value) for key, value in sources.items()}
    markers = {
        "boundary": (
            "Lambda<=1/5",
            "finite real multiple zero",
            "Hermite",
        ),
        "diagonal": (
            "R_j=38+j",
            "delta_j=1/(5j)",
            "R_j->infinity",
        ),
        "compact": (
            "|x|<=38",
            "0<=t<=1/5",
            "no contact",
        ),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in texts[key]:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        f"{key}_kind": value["kind"]
        for key, value in sources.items()
    }


def build_rows(exact: dict, audit: dict) -> list[GateRow]:
    flux = exact["flux"]
    return [
        GateRow(
            "njwg_01_universal_heat_polynomial",
            "exact_local_model",
            "ready_to_apply",
            "Every real multiple contact has a universal Hermite heat-polynomial tangent.",
            exact["universal_heat_polynomial"],
            "Parabolic leading-order model; higher terms do not change the local degree.",
            audit,
        ),
        GateRow(
            "njwg_02_discriminant",
            "exact_all_multiplicity_identity",
            "ready_to_apply",
            "The universal contact has no secondary multiple-zero time.",
            exact["discriminant"],
            "Exact for every integer m>=2.",
            hermite_audits(),
        ),
        GateRow(
            "njwg_03_root_count_jump",
            "exact_all_multiplicity_reduction",
            "ready_to_apply",
            "An m-fold contact creates exactly floor(m/2) real pairs.",
            exact["root_count_jump"],
            "Uses simplicity of Hermite roots and coefficient positivity for negative time.",
        ),
        GateRow(
            "njwg_04_local_index",
            "exact_topological_lemma",
            "ready_to_apply",
            "Every resolved pair birth has positive unit index in the standard (x,t) orientation, so an m-fold contact has index +floor(m/2).",
            f"{exact['double_index']}; {exact['multiplicity_index']}",
            "The index is local and independent of Xi arithmetic.",
        ),
        GateRow(
            "njwg_05_global_winding",
            "exact_topological_lemma",
            "ready_to_apply",
            "Interior contacts cannot cancel in first-jet winding.",
            f"{exact['global_winding']}; {exact['zero_winding_equivalence']}",
            "Requires first-jet nonvanishing on the chosen boundary.",
        ),
        GateRow(
            "njwg_06_phase_flux",
            "exact_differential_identity",
            "ready_to_apply",
            "The first-jet winding is an explicit Laguerre phase flux.",
            (
                f"{flux['Q']}; {flux['L']}; {flux['horizontal']}; "
                f"{flux['vertical']}; {flux['one_form']}"
            ),
            "Valid only away from common zeros, as required on the boundary.",
        ),
        GateRow(
            "njwg_07_xi_half_rectangles",
            "exact_equivalence",
            "ready_to_apply",
            "The positive Newman direction is equivalent to zero winding on a linear half-rectangle exhaustion.",
            exact["xi_half_rectangles"],
            "Uses positive-boundary attainment, evenness, and independent diagonal exhaustion.",
        ),
        GateRow(
            "njwg_08_known_edge_composition",
            "exact_composition",
            "ready_to_apply",
            "The top and symmetry-axis edges are already contact-free.",
            f"{exact['known_edges']}; {exact['half_rectangle_flux']}",
            "The bottom and right edge estimates remain open.",
        ),
        GateRow(
            "njwg_09_winding_guard",
            "nonpromotion_gate",
            "guard_validated",
            "Generic heat-flow structure does not make the winding vanish.",
            exact["countermodel_guard"],
            "Blocks promotion from top-time real-rootedness or axis positivity.",
        ),
        GateRow(
            "njwg_10_proof_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "Close two one-dimensional first-jet edges and one integer phase-flux condition at every stage.",
            exact["proof_handoff"],
            "No unbounded Xi edge or winding value is certified here.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    audit = source_audit()
    rows = build_rows(exact, audit)
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "exact all-multiplicity first-jet index, signed winding theorem, "
            "and diagonal boundary-flux reduction; the Xi edge separation "
            "and zero-winding family remain open and this is not a proof of "
            "Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "The gate proves sign-definite local contact charge and an exact "
            "two-dimensional-to-boundary reduction. It does not prove the "
            "bottom/right Xi edge bounds, zero winding, positive-time "
            "simplicity, Lambda<=0, or RH."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in (BOUNDARY_SOURCE, DIAGONAL_SOURCE, COMPACT_SOURCE)
        ],
        "source_audit": audit,
        "exact": exact,
        "hermite_audits": hermite_audits(),
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    audits = artifact["hermite_audits"]
    flux = exact["flux"]
    lines = [
        "# Newman First-Jet Winding Gate",
        "",
        "Date: 2026-07-26",
        "",
        "Status: exact all-multiplicity contact index and boundary-flux",
        "reduction. The Xi edge and zero-winding family remains open; this",
        "is not a proof of `Lambda <= 0` or RH.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_first_jet_winding_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_first_jet_winding_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_first_jet_winding_gate.py",
        "```",
        "",
        "## Universal Contact Charge",
        "",
        exact["orientation_convention"],
        "",
        "At a spatial zero of multiplicity `m>=2`, parabolic rescaling gives",
        "",
        "```text",
        exact["universal_heat_polynomial"],
        exact["discriminant"],
        exact["root_count_jump"],
        "```",
        "",
        "The identity `He_m'=m He_(m-1)` and the three-term Hermite",
        "recurrence evaluate the resultant at all roots of `He_(m-1)`,",
        "giving the displayed discriminant product.",
        "",
        "At every resolved double contact,",
        "",
        "```text",
        exact["double_index"],
        "```",
        "",
        "Degree stability therefore gives",
        "",
        "```text",
        exact["multiplicity_index"],
        exact["global_winding"],
        exact["zero_winding_equivalence"],
        "```",
        "",
        "For the local degree, perturb inside a fixed parabolic box until",
        "all contacts are double. Each has index `+1` in the standard",
        "`(x,t)` orientation; the real-root jump",
        "is `2floor(m/2)`, so there are exactly `floor(m/2)` pair births.",
        "Homotopy invariance then returns the stated multiple-contact index.",
        "",
        "Every interior charge is strictly positive. Hidden contacts cannot",
        "cancel one another in the boundary winding.",
        "",
        "| m | real roots before | real roots after | pair births | index | discriminant tau power |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in audits:
        lines.append(
            f"| {row['multiplicity']} | "
            f"{row['real_roots_tau_negative']} | "
            f"{row['real_roots_tau_positive']} | "
            f"{row['pair_births']} | {row['local_index']} | "
            f"{row['discriminant_tau_power']} |"
        )
    lines.extend(
        [
            "",
            "## Laguerre Phase Flux",
            "",
            "Put `Z=H+iH_x`. Away from common zeros,",
            "",
            "```text",
            flux["Q"],
            flux["L"],
            flux["horizontal"],
            flux["vertical"],
            flux["one_form"],
            "```",
            "",
            "Thus winding is not merely qualitative topology: it is a boundary",
            "integral of the first Laguerre expression and its spatial flux.",
            "",
            "## Xi Half-Rectangle Criterion",
            "",
            "```text",
            exact["xi_half_rectangles"],
            exact["known_edges"],
            exact["half_rectangle_flux"],
            "```",
            "",
            "The top time `1/4` is strictly above the published upper bound",
            "`Lambda<=1/5`, and the symmetry axis has `H_t(0)>0`. The only",
            "new separation estimates are therefore two one-dimensional",
            "edges. Once those are nonzero, the displayed integer is zero",
            "exactly when the two-dimensional stage has no contact.",
            "",
            "## Scope Guard",
            "",
            "```text",
            exact["countermodel_guard"],
            "```",
            "",
            "## Live Handoff",
            "",
            "```text",
            exact["proof_handoff"],
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman first-jet winding gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['hermite_audits'])} Hermite audits, "
        "5 exact local/global index reductions, 1 signed winding theorem, "
        "1 half-rectangle flux criterion, 1 open Xi edge target"
    )


if __name__ == "__main__":
    main()
