#!/usr/bin/env python3
"""Build the uniform-q>=1 outer-collar degree-excision reduction."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "uniform_q_ge_1_degree_excision_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "carrier_kernel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "abel_shear": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "abel_scalar_shear_flux_reduction.json"
    ),
    "global_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_global_remainder_certificate.json"
    ),
    "dominant_global": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "dominant_saddle_global_ray_certificate.json"
    ),
    "cofinal_boundary": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_cofinal_boundary_reduction.json"
    ),
    "oriented_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_oriented_successor_winding_reduction.json"
    ),
    "degree_transfer": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
    ),
    "connector_cap": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "dominant_ray_connector_phase_cap.json"
    ),
}


@dataclass(frozen=True)
class ExcisionRow:
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
        "carrier_kernel": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25",
            "|mathcal_C_N|>A_L+epsilon_term",
        ),
        "abel_shear": (
            "Gamma_(s,ell)=mathsf_X+i*",
            "X_m(x,t)=exp(m^2*t)*sin",
        ),
        "global_remainder": (
            "L>=50 and 0<tL<=25",
            "eta_0=100000 and eta_1=200000",
        ),
        "dominant_global": (
            "0<t<=1/2, L=log(x/(4*pi)), t*L>=25",
            "L_t(x)>0",
        ),
        "cofinal_boundary": (
            "t_j=25/(100+j), L_j=101+j, R_j=4*pi*exp(L_j)",
            "q<1 iff L<sqrt((100+j)/50)",
        ),
        "oriented_successor": (
            "partial D_j=B_(j+1)^[0,R_j]+V_j-reverse(B_j)",
            "kappa_j=sum_(p in D_j) floor(m_p/2)",
        ),
        "degree_transfer": (
            "deg(V_F,Omega,0)=wind",
            "+sum_p floor(m_p/2)",
        ),
        "connector_cap": (
            "|Delta_(right D_j)arg V|<pi/2",
            "H_j<3*pi/2",
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


def chain_audit() -> dict:
    outer = Counter(
        {
            "bottom_outer": 1,
            "right": 1,
            "top_outer": 1,
            "interface_outer": 1,
        }
    )
    inner = Counter(
        {
            "bottom_inner": 1,
            "interface_outer": -1,
            "top_inner": 1,
            "axis": 1,
        }
    )
    full = Counter(
        {
            "bottom_inner": 1,
            "bottom_outer": 1,
            "right": 1,
            "top_outer": 1,
            "top_inner": 1,
            "axis": 1,
        }
    )
    if inner + outer != full:
        raise RuntimeError("outer/inner interface cancellation failed")
    return {
        "outer_segments": dict(sorted(outer.items())),
        "inner_segments": dict(sorted(inner.items())),
        "full_segments": dict(sorted(full.items())),
        "identity": "partial D_j=partial Omega_j^in+partial Omega_j^out",
        "interface_cancels": True,
    }


def geometry_audit() -> dict:
    kink = Fraction(1, 2 * 50**2)
    if kink != Fraction(1, 5000):
        raise RuntimeError("L=50/q=1 kink time drifted")
    if not Fraction(25, 101) < Fraction(1, 4):
        raise RuntimeError("successor lower-time ordering failed")
    if not Fraction(50, 4) < 25:
        raise RuntimeError("L=50 interface left the Abel region")
    return {
        "kink_time": str(kink),
        "bottom_q_threshold": "(2*(25/L_j))^(-1/2)=sqrt(L_j/50)",
        "top_q_threshold": (
            "(2*(25/(L_j-1)))^(-1/2)=sqrt((L_j-1)/50)"
        ),
        "threshold_inside_right_edge": (
            "sqrt(L_j/50)<L_j and "
            "sqrt((L_j-1)/50)<L_j for L_j>=101"
        ),
        "interface_abel_side": (
            "On L_*=50, tL_*=50t<=25/2<25. "
            "On L_*=(2t)^(-1/2), tL_*=sqrt(t/2)<1<25."
        ),
    }


def build_exact() -> dict:
    chain = chain_audit()
    geometry = geometry_audit()
    return {
        "pi_provenance": (
            "The pi in L=log(x/(4*pi)) and X_*(t)=4*pi*exp(L_*(t)) "
            "is the same universal constant already fixed by the "
            "completed-zeta normalization and Riemann-Siegel saddle. "
            "The 2*pi in winding is the period of exp(i*theta). The "
            "excision introduces no circle, polygon, fitted pi, or new "
            "normalization."
        ),
        "successor_geometry": (
            "For L_j=101+j put t_b=25/L_j, t_t=25/(L_j-1), "
            "R_j=4*pi*exp(L_j), and "
            "D_j={(x,t):0<=x<=R_j, t_b<=t<=t_t}, with the standard "
            "(x,t) orientation."
        ),
        "outer_interface": (
            "For L(x)=log(x/(4*pi)), define "
            "L_*(t)=max(50,(2t)^(-1/2)) and "
            "X_*(t)=4*pi*exp(L_*(t)). Then "
            "Omega_j^out={(x,t) in D_j:x>=X_*(t)} is exactly the "
            "closed part of D_j with L>=50 and q=2tL^2>=1."
        ),
        "inner_region": (
            "Let Omega_j^in be the closure of D_j\\Omega_j^out. It "
            "contains the q<1 layer, every L<50 bounded/core shoulder, "
            "and their finite joins. The common interface is "
            "I_j={(X_*(t),t):t_b<=t<=t_t}."
        ),
        "endpoint_thresholds": (
            "At t_b and t_t the q=1 logarithmic thresholds are "
            "sqrt(L_j/50) and sqrt((L_j-1)/50), respectively. Both are "
            "strictly below L_j for L_j>=101, so the outer collar is "
            "nonempty on both horizontal edges."
        ),
        "interface_regularity": (
            "The interface is continuous and piecewise analytic. Its only "
            "possible corner is where (2t)^(-1/2)=50, namely t=1/5000. "
            "Thus Omega_j^in and Omega_j^out are compact Lipschitz "
            "domains even when one successor slab contains that corner."
        ),
        "two_regime_cover": (
            "Split Omega_j^out into A_j={tL<=25} and "
            "G_j={tL>=25}. On A_j the still-open uniform Abel-scalar "
            "gap lies in the certified L>=50, q>=1, 0<tL<=25 "
            "first-order remainder domain. On G_j the already proved "
            "dominant-saddle theorem gives L_t(x)>0 for the exact Xi "
            "flow. The two closed pieces meet at tL=25 and cover the "
            "entire outer collar."
        ),
        "interface_in_abel_domain": (
            "The whole moving interface lies strictly on the tL<25 "
            "side: if L_*=50 then tL_*=50t<=25/2, while on q=1 one has "
            "tL_*=sqrt(t/2)<1. Thus no dominant/Abel seam is placed on "
            "the degree-excision interface."
        ),
        "conditional_exact_noncontact": (
            "Assume the closed-domain Xi estimate "
            "|mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term in every prescribed "
            "canonical chart throughout A_j, including q=1, tL=25, "
            "W_0=0, endpoints, and cutoff overlaps. The checked "
            "shear/scale and first-order remainder transfers then keep "
            "the exact jet J_H=(H_t,partial_x H_t) nonzero on A_j. "
            "The dominant theorem keeps J_H nonzero on G_j. Hence J_H "
            "is nonzero on all of Omega_j^out."
        ),
        "chart_guard": (
            "mathcal_C_N is used only inside a canonical Abel chart to "
            "certify nonvanishing of the exact jet. It is not declared a "
            "global complex coordinate. Degree and interface cancellation "
            "are performed with the globally continuous exact Xi jet "
            "J_H; adjacent-cutoff and positive-normalizer homotopies "
            "retain their existing endpoint and overlap obligations."
        ),
        "oriented_chain": (
            "Orient the outer boundary as bottom outer left-to-right, "
            "right edge bottom-to-top, top outer right-to-left, and "
            "I_j top-to-bottom. Orient the inner boundary as bottom inner "
            "left-to-right, I_j bottom-to-top, top inner right-to-left, "
            "and the symmetry axis top-to-bottom. The two copies of I_j "
            "cancel, so partial D_j=partial Omega_j^in+"
            "partial Omega_j^out as oriented one-chains."
        ),
        "degree_excision": (
            "If the conditional outer noncontact theorem holds and the "
            "remaining inner boundary contract makes all displayed "
            "windings defined, then "
            "deg(J_H,Omega_j^out,0)="
            "wind(J_H(partial Omega_j^out),0)=0. Therefore "
            "kappa_j=wind(J_H(partial D_j),0)="
            "wind(J_H(partial Omega_j^in),0). All q>=1 and dominant-ray "
            "horizontal turns cancel through the zero-degree outer "
            "collar; no separate signed outer phase bound is needed."
        ),
        "many_turn_guard": (
            "For integer m>=1, X_m=exp(m^2t)sin(mx+pi/4) and "
            "C_m=partial_xX_m give a globally nonzero exact backward-heat "
            "jet. On [0,2*pi] each forward horizontal path winds -m, "
            "but the reversed upper path winds +m and the two vertical "
            "paths cancel. Thus the full zero-free strip has degree zero "
            "for every m. Arbitrarily many turns on one horizontal edge "
            "do not obstruct two-dimensional excision."
        ),
        "boundary_only_guard": (
            "F(t,x)=x^2-2t satisfies F_t=-F_xx. On "
            "[-1,1]x[-1/4,1/4], J_F=(F,F_x) is boundary-nonzero but "
            "has the interior zero (x,t)=(0,0), whose Jacobian "
            "determinant in (x,t) coordinates is +4. Hence its boundary "
            "winding is +1. Boundary nonvanishing alone cannot justify "
            "outer excision; nonvanishing on the whole closed collar is "
            "essential."
        ),
        "route_fork": (
            "Boundary-only route: retain the proved connector cap "
            "|C_j|<pi/2 and prove the joined complementary inequality "
            "H_j<3*pi/2. Uniform-collar route: prove the Abel gap on all "
            "of A_j, combine it with dominant noncontact on G_j, excise "
            "Omega_j^out, and prove only the inner-region successor "
            "degree theorem. The second route asks for a stronger "
            "pointwise domain statement but removes the independent "
            "many-turn outer phase theorem."
        ),
        "preferred_route": (
            "Promote the uniform-collar route as primary because the "
            "Abel target is already pointwise and its formulas are "
            "uniform in (x,t) on the source domain. Retain the connector "
            "quarter-turn theorem and H_j<3*pi/2 ledger as a valid "
            "fallback if the gap can be proved only on boundary arcs."
        ),
        "open_outer_target": (
            "Prove uniformly on the closed set "
            "L>=50, q=2tL^2>=1, 0<tL<=25 that "
            "|mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term, with equality boundaries, "
            "W_0=0, the recurrent endpoint, and every cutoff overlap "
            "included. This is the only new outer arithmetic theorem."
        ),
        "open_inner_target": (
            "On Omega_j^in, construct one multiplicity-compatible exact "
            "or rigorously dominated proxy covering q<1, L<50, the core, "
            "finite shoulders, and all joins. Prove boundary "
            "nonvanishing and an oriented winding below one (or direct "
            "zero degree) for every j. Positive local contact degree then "
            "forces the localized successor integer to vanish."
        ),
        "proof_boundary": (
            "The successor geometry, moving q=1/L=50 interface, "
            "two-regime cover, oriented chain cancellation, conditional "
            "degree excision, paired many-turn cancellation guard, "
            "boundary-only contact guard, and route comparison are exact "
            "or inherited from checked sources. This does not prove the "
            "uniform Xi Abel-scalar gap, its exact-H noncontact "
            "consequence, the inner q<1/bounded-L successor theorem, "
            "complete boundary nonvanishing, contact exclusion, "
            "Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion."
        ),
        "diagnostics": {
            "geometry": geometry,
            "chain": chain,
        },
    }


def build_rows(exact: dict) -> list[ExcisionRow]:
    diagnostics = exact["diagnostics"]
    return [
        ExcisionRow(
            "uqde_00_pi_provenance",
            "source_provenance",
            "available_exact",
            "The outer interface introduces no unexplained pi.",
            exact["pi_provenance"],
            "Only the inherited Xi normalization and winding period occur.",
        ),
        ExcisionRow(
            "uqde_01_successor_geometry",
            "exact_geometry",
            "available_exact",
            "Every successor slab has one fixed standard orientation.",
            exact["successor_geometry"],
            "Valid for every integer j>=0.",
        ),
        ExcisionRow(
            "uqde_02_outer_interface",
            "exact_geometry",
            "available_exact",
            "One moving interface cuts out exactly L>=50 and q>=1.",
            exact["outer_interface"],
            "The logarithmic coordinate is used only for x>0.",
        ),
        ExcisionRow(
            "uqde_03_endpoint_thresholds",
            "exact_geometry",
            "available_exact",
            "The outer collar reaches both successor horizontal edges.",
            exact["endpoint_thresholds"],
            "Strict for L_j>=101.",
            diagnostics["geometry"],
        ),
        ExcisionRow(
            "uqde_04_interface_regularity",
            "exact_geometry",
            "available_exact",
            "The moving interface is a valid Lipschitz degree boundary.",
            exact["interface_regularity"],
            "The t=1/5000 corner needs no smoothing.",
        ),
        ExcisionRow(
            "uqde_05_two_regime_cover",
            "exact_domain_reduction",
            "available_exact",
            "The Abel and dominant regimes cover the whole outer collar.",
            exact["two_regime_cover"],
            "The Abel estimate remains open; dominant noncontact is imported.",
        ),
        ExcisionRow(
            "uqde_06_interface_abel_side",
            "exact_domain_reduction",
            "available_exact",
            "The degree interface stays strictly inside the Abel regime.",
            exact["interface_in_abel_domain"],
            "No unproved dominant/Abel seam is hidden on I_j.",
            diagnostics["geometry"],
        ),
        ExcisionRow(
            "uqde_07_conditional_noncontact",
            "conditional_exact_reduction",
            "available_exact",
            "The uniform Abel gap would make the exact outer jet nonzero.",
            exact["conditional_exact_noncontact"],
            "Conditional on the explicitly open closed-domain Xi gap.",
        ),
        ExcisionRow(
            "uqde_08_chart_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Degree is taken in the exact global jet, not an Abel chart.",
            exact["chart_guard"],
            "Cutoff overlaps and endpoint tracks may not be discarded.",
        ),
        ExcisionRow(
            "uqde_09_oriented_chain",
            "exact_topological_lemma",
            "available_exact",
            "The inner and outer chains add to the successor boundary.",
            exact["oriented_chain"],
            "The common interface occurs with opposite orientations.",
            diagnostics["chain"],
        ),
        ExcisionRow(
            "uqde_10_degree_excision",
            "conditional_topological_reduction",
            "available_exact",
            "Uniform outer noncontact localizes the successor degree.",
            exact["degree_excision"],
            "The inner boundary theorem remains open.",
        ),
        ExcisionRow(
            "uqde_11_many_turn_guard",
            "countermodel",
            "guard_validated",
            "Paired horizontal turns cancel in a zero-free strip.",
            exact["many_turn_guard"],
            "Generic heat-flow guard, not an Xi theorem.",
        ),
        ExcisionRow(
            "uqde_12_boundary_only_guard",
            "countermodel",
            "guard_validated",
            "Boundary nonvanishing does not imply zero interior degree.",
            exact["boundary_only_guard"],
            "Generic heat-flow guard, not an Xi counterexample.",
        ),
        ExcisionRow(
            "uqde_13_route_fork",
            "exact_route_comparison",
            "available_exact",
            "Boundary-only and uniform-collar routes have different obligations.",
            exact["route_fork"],
            "Neither route is promoted to a completed Xi theorem.",
        ),
        ExcisionRow(
            "uqde_14_outer_gap_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The closed-domain Xi Abel-scalar gap remains open.",
            exact["open_outer_target"],
            "This is an RH-level arithmetic estimate.",
        ),
        ExcisionRow(
            "uqde_15_inner_degree_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The localized inner successor theorem remains open.",
            exact["open_inner_target"],
            "It includes q<1, bounded L, core, shoulders, and joins.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact conditional uniform-q>=1 outer-collar degree-excision "
            "reduction; the uniform Xi Abel gap and inner successor "
            "degree theorem remain open"
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in build_rows(exact)],
        "proof_boundary": exact["proof_boundary"],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return f"""# Newman Uniform-q>=1 Degree-Excision Reduction

Date: 2026-07-28

Status: exact conditional outer-collar degree reduction. This is not a
proof of the uniform Xi Abel gap, the inner successor theorem,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Pi Provenance

{exact["pi_provenance"]}

## Successor Collar

```text
{exact["successor_geometry"]}
{exact["outer_interface"]}
{exact["inner_region"]}
{exact["endpoint_thresholds"]}
{exact["interface_regularity"]}
```

## Two-Regime Cover

```text
{exact["two_regime_cover"]}
{exact["interface_in_abel_domain"]}
```

The first-order remainder is used only on `tL<=25`. The already proved
dominant-saddle theorem supplies exact noncontact on `tL>=25`; the two
pieces meet at equality.

## Conditional Noncontact

```text
{exact["conditional_exact_noncontact"]}
{exact["chart_guard"]}
```

## Oriented Excision

```text
{exact["oriented_chain"]}
{exact["degree_excision"]}
```

## Route Guards

```text
{exact["many_turn_guard"]}

{exact["boundary_only_guard"]}
```

The first guard shows why arbitrarily many turns on one horizontal edge do
not obstruct a zero-free two-dimensional collar. The second shows why
boundary-only nonvanishing is too weak for the same conclusion.

## Route Fork

```text
{exact["route_fork"]}
{exact["preferred_route"]}
```

## Open Theorems

Outer:

```text
{exact["open_outer_target"]}
```

Inner:

```text
{exact["open_inner_target"]}
```

## Boundary

{exact["proof_boundary"]}
"""


def write_artifact(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_artifact(args.out, artifact)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman uniform-q>=1 degree excision: "
        "16 rows, 1 exact interface chain, "
        "1 conditional outer-degree excision, "
        "1 paired many-turn cancellation guard, "
        "1 boundary-only contact guard, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
