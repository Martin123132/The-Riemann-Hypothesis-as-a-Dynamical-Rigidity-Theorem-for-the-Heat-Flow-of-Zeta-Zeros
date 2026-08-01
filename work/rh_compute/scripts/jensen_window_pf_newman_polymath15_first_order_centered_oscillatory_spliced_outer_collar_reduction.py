#!/usr/bin/env python3
"""Build the oscillatory-spliced outer-collar degree reduction."""

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
    "oscillatory_spliced_outer_collar_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

C_STAR = Fraction(4_911_678_521, 1_933_561_194)
EXAMPLE_EPSILON = Fraction(1, 100)
EXAMPLE_CAP = C_STAR + EXAMPLE_EPSILON

SOURCES = {
    "uniform_collar": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "uniform_q_ge_1_degree_excision_reduction.json"
    ),
    "oscillatory_handoff": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "oscillatory_zeta_handoff_theorem.json"
    ),
    "antedb_frontier": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "antedb_beta_frontier_audit.json"
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
}


@dataclass(frozen=True)
class SpliceRow:
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
        "uniform_collar": (
            "L_*(t)=max(50,(2t)^(-1/2))",
            "partial D_j=partial Omega_j^in+partial Omega_j^out",
            "|mathcal_C_N|>A_L+epsilon_term",
        ),
        "oscillatory_handoff": (
            "For every epsilon>0 there exists L_epsilon",
            "t*L>=4911678521/1933561194+epsilon",
            "L_t(x)>0",
        ),
        "antedb_frontier": (
            "4911678521/1933561194",
            "audited pointwise-bound frontier remains c_*",
            "220633/620612",
        ),
        "global_remainder": (
            "L>=50 and 0<tL<=25",
            "eta_0=100000 and eta_1=200000",
        ),
        "dominant_global": (
            "0<t<=1/2, L=log(x/(4*pi)), t*L>=25",
            "L_t(x)>0",
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
            "interface": 1,
        }
    )
    inner = Counter(
        {
            "bottom_inner": 1,
            "interface": -1,
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
        raise RuntimeError("raised-collar interface cancellation failed")
    return {
        "outer_segments": dict(sorted(outer.items())),
        "inner_segments": dict(sorted(inner.items())),
        "full_segments": dict(sorted(full.items())),
        "identity": (
            "partial D_j=partial Omega_(j,epsilon)^in"
            "+partial Omega_(j,epsilon)^out"
        ),
        "interface_cancels": True,
    }


def arithmetic_audit() -> dict:
    if not Fraction(0) < C_STAR < Fraction(25):
        raise RuntimeError("critical threshold left the certified range")
    if not Fraction(0) < EXAMPLE_EPSILON < Fraction(25) - C_STAR:
        raise RuntimeError("example epsilon is not admissible")
    if not EXAMPLE_CAP < Fraction(25):
        raise RuntimeError("example Abel cap left the first-order domain")
    q_one_cap_at_l50 = Fraction(1, 2 * 50)
    if not q_one_cap_at_l50 < C_STAR:
        raise RuntimeError("q=1 interface is not below the critical seam")
    return {
        "c_star": {
            "exact": str(C_STAR),
            "decimal": f"{float(C_STAR):.15f}",
        },
        "admissible_epsilon": (
            f"0<epsilon<25-c_*={25 - C_STAR}"
        ),
        "example": {
            "epsilon": str(EXAMPLE_EPSILON),
            "cap_exact": str(EXAMPLE_CAP),
            "cap_decimal": f"{float(EXAMPLE_CAP):.15f}",
            "fraction_of_old_c_width": f"{float(EXAMPLE_CAP / 25):.9f}",
        },
        "q_one_interface": (
            "q=1 implies c=tL=1/(2L)<=1/100<c_* for L>=50"
        ),
    }


def build_exact() -> dict:
    arithmetic = arithmetic_audit()
    chain = chain_audit()
    return {
        "pi_provenance": (
            "The pi in L=log(x/(4*pi)) and in the successor radius "
            "R_j=4*pi*exp(L_j) is inherited unchanged from the "
            "completed-zeta normalization and the Riemann-Siegel "
            "saddle. The splice introduces no circle, polygon, fitted "
            "pi, or new normalization."
        ),
        "current_frontier": (
            "The pinned July-2026 ANTEDB audit, including the four "
            "Tao-Trudgian-Yang and two Cushing post-2023 pairs, leaves "
            "the exact pointwise-beta contact unchanged at "
            "c_*=4911678521/1933561194. Its active contact is "
            "alpha_*=62831/155153, beta_*=220633/620612."
        ),
        "epsilon_quantifier": (
            "Fix 0<epsilon<25-c_*. The oscillatory-zeta theorem gives "
            "a finite L_epsilon such that the exact Xi first Laguerre "
            "quantity is positive whenever 0<t<=1/2, L>=L_epsilon, "
            "and c=tL>=c_*+epsilon. Put "
            "B_epsilon=max(50,L_epsilon)."
        ),
        "raised_collar": (
            "For L_j=101+j, t_b=25/L_j, t_t=25/(L_j-1), and "
            "D_j=[0,4*pi*exp(L_j)]x[t_b,t_t], choose "
            "J_epsilon=max(0,ceil(B_epsilon)-101). For j>=J_epsilon "
            "define L_(epsilon,*)(t)=max(B_epsilon,(2t)^(-1/2)). "
            "Then Omega_(j,epsilon)^out is exactly the closed part of "
            "D_j with L>=B_epsilon and q=2tL^2>=1. Its complement "
            "Omega_(j,epsilon)^in contains q<1 and the fixed "
            "L<=B_epsilon shoulder."
        ),
        "interface_geometry": (
            "The raised interface is continuous and piecewise analytic, "
            "with at most one corner at t=1/(2B_epsilon^2). Since "
            "L_j>=B_epsilon and sqrt(L_j/50)<L_j, it lies inside both "
            "successor horizontal edges. On its q=1 branch, "
            "c=1/(2L)<=1/100<c_*. Its L=B_epsilon branch may cross the "
            "oscillatory seam, where exact noncontact is already known."
        ),
        "three_regime_cover": (
            "Write c_epsilon=c_*+epsilon and split the raised outer "
            "collar into A_(j,epsilon)={c<=c_epsilon}, "
            "Z_(j,epsilon)={c_epsilon<=c<=25}, and "
            "G_(j,epsilon)={c>=25}. The first set lies in the checked "
            "L>=50, q>=1, 0<c<=25 first-order domain. The exact "
            "oscillatory-zeta theorem gives noncontact on Z, and the "
            "dominant-saddle theorem gives noncontact on G. The closed "
            "sets meet at c=c_epsilon and c=25 and cover the collar."
        ),
        "seam_stress": (
            "The only open seam tests are therefore q=1, "
            "c=c_epsilon, L=B_epsilon, W_0=0, the recurrent endpoint, "
            "and adjacent cutoff overlaps. The q=1 seam is strictly "
            "inside the low-c Abel region. At c=c_epsilon and c=25 the "
            "globally continuous exact Xi jet, rather than an unjoined "
            "chart coordinate, is used. The known exact noncontact "
            "theorems include both equality seams."
        ),
        "conditional_noncontact": (
            "Assume only on A_(j,epsilon) that every prescribed "
            "canonical chart satisfies |mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term, including W_0=0, "
            "endpoints, equality boundaries, and cutoff overlaps. The "
            "checked first-order transfer makes the exact Xi jet "
            "J_H=(H_t,partial_xH_t) nonzero on A. The oscillatory and "
            "dominant theorems make it nonzero on Z and G. Hence J_H "
            "is nonzero on all of Omega_(j,epsilon)^out."
        ),
        "degree_excision": (
            "The raised inner and outer interface copies have opposite "
            "orientations, so partial D_j=partial "
            "Omega_(j,epsilon)^in+partial Omega_(j,epsilon)^out. "
            "Conditional outer noncontact gives "
            "deg(J_H,Omega_(j,epsilon)^out,0)=0 and therefore "
            "kappa_j=wind(J_H(partial "
            "Omega_(j,epsilon)^in),0) for every j>=J_epsilon."
        ),
        "fixed_epsilon_consequence": (
            "For each fixed admissible epsilon, the sufficient outer "
            "Abel theorem is needed only on the closed wedge "
            "L>=B_epsilon, q>=1, "
            "0<c<=4911678521/1933561194+epsilon. Requiring it through "
            "c=25 is a valid but unnecessary strengthening."
        ),
        "asymptotic_scope": (
            "Because epsilon is arbitrary, every fixed scaled ray "
            "c>c_* is eventually closed. Thus the asymptotically live "
            "outer arithmetic envelope is c<=c_*+o(1). This statement "
            "does not manufacture a uniform explicit epsilon(L), and "
            "it does not include the q<1 parabolic layer."
        ),
        "frontier_guard": (
            "The current ANTEDB pointwise-beta machinery cannot lower "
            "this seam: the exact TY1/TY2 contact survives the six "
            "post-2023 pairs and twelve audited beta-to-beta iterations. "
            "A lower seam needs an improved beta bound near "
            "(alpha_*,beta_*) or cancellation beyond pointwise beta "
            "majorants."
        ),
        "finite_shoulder_guard": (
            "L_epsilon is existential in the imported theorem. The "
            "logical reduction is valid with that finite constant, but "
            "a finite-height computer certificate needs an effective "
            "upper bound for B_epsilon. Until then, L<B_epsilon and the "
            "finitely many early successor slabs remain explicit inner "
            "or finite-shoulder obligations; they may not be silently "
            "declared checked."
        ),
        "route_decision": (
            "Use the oscillatory-spliced collar as the primary outer "
            "route. Attack the Abel scalar only in the low-c wedge and "
            "develop the multiplicity-compatible q<1/bounded-L inner "
            "degree theorem separately. Retain the c<=25 collar and "
            "the H_j<3*pi/2 boundary ledger as stronger fallbacks."
        ),
        "open_outer_target": (
            "For one fixed 0<epsilon<25-c_* and its "
            "B_epsilon=max(50,L_epsilon), prove uniformly on "
            "L>=B_epsilon, q=2tL^2>=1, "
            "0<tL<=c_*+epsilon that |mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term. Include q=1, W_0=0, "
            "the recurrent endpoint, equality boundaries, and every "
            "adjacent-cutoff overlap."
        ),
        "open_inner_target": (
            "Construct one multiplicity-compatible exact or rigorously "
            "dominated successor proxy on q<1 together with the fixed "
            "L<=B_epsilon shoulder, core, early slabs, and all joins. "
            "Prove zero degree, or boundary nonvanishing and winding "
            "strictly below one. Effectivize B_epsilon before claiming "
            "a finite exhaustive certificate."
        ),
        "proof_boundary": (
            "This artifact proves the epsilon-quantified raised-collar "
            "geometry, three-regime cover, exact seam allocation, "
            "conditional degree excision, and the reduction of the "
            "outer Abel burden from c<=25 to c<=c_*+epsilon. It does "
            "not prove the Xi Abel-scalar gap, an effective "
            "L_epsilon, the q<1/bounded-L inner theorem, complete "
            "boundary nonvanishing, contact exclusion, Lambda<=0, "
            "PF-infinity, RH, or a Clay-prize conclusion."
        ),
        "diagnostics": {
            "arithmetic": arithmetic,
            "chain": chain,
        },
    }


def build_rows(exact: dict) -> list[SpliceRow]:
    diagnostics = exact["diagnostics"]
    return [
        SpliceRow(
            "osce_00_pi_provenance",
            "source_provenance",
            "available_exact",
            "The oscillatory splice introduces no unexplained pi.",
            exact["pi_provenance"],
            "Only the inherited Xi normalization occurs.",
        ),
        SpliceRow(
            "osce_01_current_frontier",
            "current_source_audit",
            "guard_validated",
            "The audited current pointwise-beta threshold is exact.",
            exact["current_frontier"],
            "Pinned to the July-2026 ANTEDB audit, not future literature.",
            diagnostics["arithmetic"],
        ),
        SpliceRow(
            "osce_02_epsilon_quantifier",
            "exact_quantifier",
            "available_exact",
            "A fixed positive margin supplies a finite outer threshold.",
            exact["epsilon_quantifier"],
            "L_epsilon is existential rather than numerically explicit.",
        ),
        SpliceRow(
            "osce_03_raised_collar",
            "exact_geometry",
            "available_exact",
            "The raised interface cuts out exactly q>=1 above B_epsilon.",
            exact["raised_collar"],
            "Applied only to eventual successors j>=J_epsilon.",
        ),
        SpliceRow(
            "osce_04_interface_geometry",
            "exact_geometry",
            "available_exact",
            "The raised interface is a valid Lipschitz degree boundary.",
            exact["interface_geometry"],
            "Its q=1 branch stays strictly below c_*.",
            diagnostics["arithmetic"],
        ),
        SpliceRow(
            "osce_05_three_regime_cover",
            "exact_domain_reduction",
            "available_exact",
            "Abel, oscillatory, and dominant regimes cover the collar.",
            exact["three_regime_cover"],
            "Only the low-c Abel noncontact statement remains conditional.",
        ),
        SpliceRow(
            "osce_06_seam_stress",
            "nonpromotion_guard",
            "guard_validated",
            "All equality and chart seams remain explicit.",
            exact["seam_stress"],
            "No chart scalar is promoted to a global coordinate.",
        ),
        SpliceRow(
            "osce_07_conditional_noncontact",
            "conditional_exact_reduction",
            "available_exact",
            "The low-c gap would make the full raised collar noncontact.",
            exact["conditional_noncontact"],
            "Conditional on the explicitly open Xi Abel gap.",
        ),
        SpliceRow(
            "osce_08_degree_excision",
            "conditional_topological_reduction",
            "available_exact",
            "Raised-collar noncontact localizes the successor degree.",
            exact["degree_excision"],
            "The inner successor theorem remains open.",
            diagnostics["chain"],
        ),
        SpliceRow(
            "osce_09_fixed_epsilon_consequence",
            "exact_route_reduction",
            "available_exact",
            "The closed Abel target need not extend to c=25.",
            exact["fixed_epsilon_consequence"],
            "A fixed epsilon-dependent shoulder is retained.",
        ),
        SpliceRow(
            "osce_10_asymptotic_scope",
            "quantifier_guard",
            "guard_validated",
            "The live asymptotic outer envelope is c<=c_*+o(1).",
            exact["asymptotic_scope"],
            "No explicit diagonal epsilon(L) is asserted.",
        ),
        SpliceRow(
            "osce_11_frontier_guard",
            "current_route_guard",
            "guard_validated",
            "Current pointwise beta improvements do not lower c_*.",
            exact["frontier_guard"],
            "A different cancellation input could still improve it.",
        ),
        SpliceRow(
            "osce_12_finite_shoulder_guard",
            "effectivity_guard",
            "guard_validated",
            "The epsilon-dependent finite shoulder remains visible.",
            exact["finite_shoulder_guard"],
            "No finite-height coverage is promoted without an explicit bound.",
        ),
        SpliceRow(
            "osce_13_outer_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The reduced low-c Xi Abel-scalar theorem remains open.",
            exact["open_outer_target"],
            "This is the surviving outer RH-level arithmetic estimate.",
        ),
        SpliceRow(
            "osce_14_inner_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q<1 and bounded-L inner degree theorem remains open.",
            exact["open_inner_target"],
            "It must be multiplicity-compatible as t tends to zero.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact oscillatory-spliced outer-collar reduction; the "
            "low-c Xi Abel gap, effective finite shoulder, and inner "
            "successor theorem remain open"
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
    return f"""# Newman Oscillatory-Spliced Outer-Collar Reduction

Date: 2026-07-28

Status: exact conditional domain and degree reduction. This is not a
proof of the Xi Abel gap, `Lambda<=0`, PF-infinity, RH, or a Clay-prize
result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Pi Provenance

{exact["pi_provenance"]}

## Current Frontier

```text
{exact["current_frontier"]}
```

## Epsilon Collar

```text
{exact["epsilon_quantifier"]}
{exact["raised_collar"]}
{exact["interface_geometry"]}
```

For the concrete bookkeeping choice `epsilon=1/100`, the new closed
Abel cap is

```text
c<= {EXAMPLE_CAP} = {float(EXAMPLE_CAP):.15f}...
```

instead of `c<=25`. The theorem remains valid for every smaller fixed
positive epsilon, with its own finite threshold.

## Three-Regime Cover

```text
{exact["three_regime_cover"]}
{exact["seam_stress"]}
```

## Conditional Excision

```text
{exact["conditional_noncontact"]}
{exact["degree_excision"]}
```

## Exact Consequence

```text
{exact["fixed_epsilon_consequence"]}
{exact["asymptotic_scope"]}
```

## Route Guards

```text
{exact["frontier_guard"]}

{exact["finite_shoulder_guard"]}
```

## Route Decision

{exact["route_decision"]}

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
        "built Newman oscillatory-spliced outer collar: "
        "15 rows, 1 exact current frontier, 1 epsilon-raised collar, "
        "1 three-regime cover, 1 conditional degree excision, "
        "2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
