#!/usr/bin/env python3
"""Build the residual-placement and cumulative-current handoff gate."""

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
    "residual_placement_c1_cumulative_handoff_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "global_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_global_remainder_certificate.json"
    ),
    "rectangular_transfer": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_rectangular_boundary_degree_reduction.json"
    ),
    "critical_ray_edge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "critical_ray_finite_height_real_edge_gate.json"
    ),
    "cutoff_edge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_cutoff_real_edge_projective_splice_gate.json"
    ),
    "interior_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "interior_projective_current_gate.json"
    ),
    "pairwise_occupation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "pairwise_projective_alignment_gate.json"
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


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    remainder = payloads["global_remainder"].get("exact", {})
    if remainder.get("global_remainder") != (
        "For every critical disk, including cutoff crossings, "
        "|r_[1](x)|<100000*exp(-5L/4) and "
        "|partial_x r_[1](x)|<200000*L*exp(-5L/4)"
    ):
        raise RuntimeError("global C1 remainder drifted")

    rectangular = payloads["rectangular_transfer"].get("exact", {})
    if rectangular.get("rectangular_rouche", {}).get("homotopy") != (
        "V_s=V_J+s*V_r, M_L(V_s)>=M_L(V_J)-s*M_L(V_r)>0"
    ):
        raise RuntimeError("rectangular homotopy drifted")
    if rectangular.get("adjacent_chart", {}).get("gauge") != (
        "M_L(Delta J)<1/5"
    ):
        raise RuntimeError("adjacent chart gauge drifted")

    edge = payloads["critical_ray_edge"].get("counts", {})
    if edge.get("uniform_q_ge_1_model_signs") != 1:
        raise RuntimeError("critical-ray edge sign drifted")
    cutoff = payloads["cutoff_edge"].get("counts", {})
    if cutoff.get("signed_affine_splices") != 1:
        raise RuntimeError("cutoff edge splice drifted")

    interior = payloads["interior_current"].get("exact", {})
    if "C_edge=c_0+c_N" not in interior.get("edge_block", ""):
        raise RuntimeError("interior edge decomposition drifted")
    if "J_n<=-(5/64)" not in interior.get(
        "negative_definite_theorem", ""
    ):
        raise RuntimeError("interior current theorem drifted")

    occupation = payloads["pairwise_occupation"].get("exact", {})
    if "mathcal_C_N=D_edge+D_perp" not in occupation.get(
        "aggregate_oriented_cell", ""
    ):
        raise RuntimeError("aggregate occupation decomposition drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind") for key, payload in payloads.items()
        },
    }


def symbolic_certificate() -> dict:
    c, d, u, w = sp.symbols("c d u w", real=True)
    c_x, d_x, u_x, w_x = sp.symbols(
        "c_x d_x u_x w_x", real=True
    )
    lam = sp.symbols("lambda", real=True)

    base = c * d_x - d * c_x
    allocated = (c + lam * u) * (d_x + lam * w_x) - (
        d + lam * w
    ) * (c_x + lam * u_x)
    difference = sp.expand(allocated - base)
    expected = sp.expand(
        lam * (c * w_x + u * d_x - d * u_x - w * c_x)
        + lam**2 * (u * w_x - w * u_x)
    )
    if sp.simplify(difference - expected) != 0:
        raise RuntimeError("allocation-current identity failed")

    witness = sp.simplify(
        difference.subs(
            {
                c: 1,
                d: 0,
                c_x: 0,
                d_x: -1,
                u: 0,
                w: 1,
                u_x: -2,
                w_x: 0,
            }
        )
        + base.subs({c: 1, d: 0, c_x: 0, d_x: -1})
    )
    if sp.simplify(witness - (-1 + 2 * lam**2)) != 0:
        raise RuntimeError("allocation sign witness failed")

    x, ell = sp.symbols("x L", positive=True)
    r_x, r_xx = sp.symbols("r_x r_xx", real=True)
    normalized_second_derivative = r_xx / ell - r_x / (x * ell**2)

    cs = sp.symbols("c0:3", real=True)
    ds = sp.symbols("d0:3", real=True)
    cxs = sp.symbols("cx0:3", real=True)
    dxs = sp.symbols("dx0:3", real=True)
    total = sp.expand(sum(cs) * sum(dxs) - sum(ds) * sum(cxs))
    diagonal = sum(cs[i] * dxs[i] - ds[i] * cxs[i] for i in range(3))
    cross = sum(
        cs[i] * dxs[j]
        + cs[j] * dxs[i]
        - ds[i] * cxs[j]
        - ds[j] * cxs[i]
        for i in range(3)
        for j in range(i + 1, 3)
    )
    if sp.simplify(total - diagonal - cross) != 0:
        raise RuntimeError("cumulative polarization identity failed")

    theta = sp.symbols("theta", real=True)
    v1 = sp.Matrix([sp.cos(theta), -sp.sin(theta)])
    v2 = sp.Matrix([sp.cos(3 * theta), -sp.sin(3 * theta)])

    def current(vector: sp.Matrix) -> sp.Expr:
        derivative = vector.diff(theta)
        return sp.expand(vector[0] * derivative[1] - vector[1] * derivative[0])

    if sp.simplify(current(v1) + 1) != 0:
        raise RuntimeError("first clockwise witness failed")
    if sp.simplify(current(v2) + 3) != 0:
        raise RuntimeError("second clockwise witness failed")
    joined = sp.simplify((v1 + v2).subs(theta, sp.pi / 2))
    if joined != sp.zeros(2, 1):
        raise RuntimeError("clockwise-sum zero witness failed")

    return {
        "allocation_current": (
            "K(v+lambda e)-K(v)="
            "lambda[c*w_x+u*d_x-d*u_x-w*c_x]"
            "+lambda^2[u*w_x-w*u_x]"
        ),
        "four_defect_specialization": (
            "K(v+e)-K(v)=c*w_x+u*d_x+u*w_x"
            "-d*u_x-w*c_x-w*u_x"
        ),
        "allocation_sign_witness": (
            "v=(1,0), v_x=(0,-1), e=(0,1), e_x=(-2,0) "
            "gives K(v+lambda e)=-1+2lambda^2"
        ),
        "normalized_residual_derivative": (
            "partial_x[(partial_x r_[1])/L]="
            "partial_x^2 r_[1]/L-partial_x r_[1]/(xL^2)"
        ),
        "cumulative_polarization": (
            "K(sum_a v_a)=sum_a K(v_a)+sum_(a<b)"
            "[det(v_a,partial_x v_b)+det(v_b,partial_x v_a)]"
        ),
        "clockwise_sum_guard": (
            "v_1=(cos(theta),-sin(theta)), "
            "v_2=(cos(3theta),-sin(3theta)); "
            "K(v_1)=-1, K(v_2)=-3, but v_1+v_2=0 at theta=pi/2"
        ),
        "normalized_second_derivative_expression": str(
            normalized_second_derivative
        ),
    }


def exact_payload() -> dict:
    return {
        "domain": "L>=50, 0<tL<=25 on prescribed first-order charts",
        "canonical_whole_jet": (
            "V_Z=(Z,partial_x Z/L)=V_J+V_r with "
            "V_J=(J_[1],partial_x J_[1]/L) and "
            "V_r=(r_[1],partial_x r_[1]/L)"
        ),
        "canonical_main_decomposition": (
            "Inside V_J only, (mathsf X,mathcal C_N)="
            "(C_edge,D_edge)+sum_(n=1)^(N-1)(c_n,d_n)"
        ),
        "allocation_gauge": (
            "For any constant lambda, assign lambda V_r to the edge and "
            "(1-lambda)V_r to the complementary main. Their sum remains "
            "V_Z, while the individual edge current varies with lambda."
        ),
        "c1_box": (
            "|r_[1]|<B_0=100000exp(-5L/4), "
            "|partial_x r_[1]/L|<B_1=200000exp(-5L/4)"
        ),
        "boundary_transfer": (
            "M_L(V_J)>1 and M_L(V_r)<1 imply "
            "M_L(V_J+sV_r)>0 for 0<=s<=1; no derivative of V_r is used"
        ),
        "chart_covariance": (
            "On an adjacent overlap, V_(J,N)+V_(r,N)="
            "V_(J,N+1)+V_(r,N+1)=V_Z, hence "
            "Delta V_r=-Delta V_J"
        ),
        "retained_edge": (
            "The edge block is strictly clockwise in every prescribed cell "
            "and across every aligned cutoff for q>=1; it remains a component "
            "theorem for V_J"
        ),
        "diagonal_reserve": (
            "R_diag=-K(v_edge)-sum_(n=1)^(N-1)K(v_n)>0; "
            "K(V_J)=-R_diag+C_cross"
        ),
        "open_invariant_target": (
            "Prove C_cross<R_diag in the required orientation, or prove the "
            "division-free Abel gap |mathsf_X|<=delta_L => "
            "|mathcal C_N|>A_L+epsilon_term. Then transfer the resulting "
            "whole-main boundary theorem to Xi with the C1 homotopy."
        ),
        "route_correction": (
            "A C2 estimate would be needed only after choosing a noncanonical "
            "residual allocation and differentiating that allocated first-jet "
            "vector. It is not required by the invariant whole-jet boundary "
            "Rouche theorem."
        ),
    }


def build_rows(exact: dict, symbolic: dict) -> list[GateRow]:
    return [
        GateRow(
            "rpc_01_whole_jet",
            "exact source identity",
            "proved",
            "The complete normalized Xi first jet is main plus one global remainder jet.",
            exact["canonical_whole_jet"],
            "This identity is made only after the retained components are summed.",
        ),
        GateRow(
            "rpc_02_c1_box",
            "uniform analytic input",
            "proved",
            "The global residual has two certified first-jet coordinates.",
            exact["c1_box"],
            "No global second x derivative is asserted here.",
        ),
        GateRow(
            "rpc_03_main_components",
            "canonical decomposition",
            "proved",
            "The endpoint-terminal edge is a component of the retained main.",
            exact["canonical_main_decomposition"],
            "The global remainder is not one of these retained atoms.",
        ),
        GateRow(
            "rpc_04_allocation_gauge",
            "nonpromotion guard",
            "proved",
            "Allocating the global remainder to one component is nonunique.",
            exact["allocation_gauge"],
            "An allocated component current is not an invariant of Xi.",
        ),
        GateRow(
            "rpc_05_four_defects",
            "exact algebra",
            "proved",
            "The former four-defect expression is the current change under one chosen allocation.",
            symbolic["four_defect_specialization"],
            "The identity is valid, but it does not canonically define a complete Xi edge.",
        ),
        GateRow(
            "rpc_06_allocation_witness",
            "sign guard",
            "proved",
            "The allocated edge-current sign can change while the whole vector is held fixed by compensation.",
            symbolic["allocation_sign_witness"],
            "The witness is algebraic and is not an asserted Xi state.",
        ),
        GateRow(
            "rpc_07_derivative_order",
            "derivative audit",
            "proved",
            "Differentiating an allocated first-jet residual introduces r_[1],xx.",
            symbolic["normalized_residual_derivative"],
            "This explains the apparent C2 demand without making it invariant.",
        ),
        GateRow(
            "rpc_08_c2_nonrequirement",
            "route correction",
            "proved",
            "The whole-jet boundary homotopy uses C1 scalar control and no C2 residual estimate.",
            exact["boundary_transfer"],
            "This is a boundary-winding transfer, not pointwise interior exclusion.",
        ),
        GateRow(
            "rpc_09_chart_covariance",
            "cutoff identity",
            "proved",
            "Matched adjacent remainders compensate the main-chart jump exactly.",
            exact["chart_covariance"],
            "Chart covariance does not assign the jump to the endpoint-terminal edge.",
        ),
        GateRow(
            "rpc_10_cutoff_roles",
            "join separation",
            "proved",
            "The signed projective edge splice and the C1 whole-jet chart transfer have different roles.",
            "Delta_cut<-h/25<0 for the retained edge; M_L(Delta J)<1/5 for the whole main chart.",
            "Neither statement is promoted into an Xi-level edge current.",
        ),
        GateRow(
            "rpc_11_retained_edge",
            "component theorem",
            "proved",
            "The completed edge theorem remains valid input to the retained main.",
            exact["retained_edge"],
            "The correction changes only the handoff, not the edge margins.",
        ),
        GateRow(
            "rpc_12_polarization",
            "exact cumulative identity",
            "proved",
            "The projective current of the retained sum contains every pair cross current.",
            symbolic["cumulative_polarization"],
            "No diagonal component theorem controls the cross aggregate by itself.",
        ),
        GateRow(
            "rpc_13_diagonal_reserve",
            "reduced target",
            "conditional",
            "The signed edge and interior currents create a positive diagonal reserve.",
            exact["diagonal_reserve"],
            "Strict cumulative clockwise current still requires C_cross<R_diag.",
        ),
        GateRow(
            "rpc_14_sum_guard",
            "generic counterexample",
            "proved",
            "Strictly clockwise component currents are not closed under addition.",
            symbolic["clockwise_sum_guard"],
            "This is a route guard, not a counterexample to the Xi estimate.",
        ),
        GateRow(
            "rpc_15_invariant_handoff",
            "open theorem target",
            "open",
            "The next admissible theorem is cumulative cross-current control or the endpoint-complete Abel gap.",
            exact["open_invariant_target"],
            "Every endpoint, projection pole, exceptional fibre, and cutoff join remains included.",
        ),
        GateRow(
            "rpc_16_boundary",
            "proof boundary",
            "proved",
            "The noncanonical Xi-edge C2 obligation is retired without claiming contact exclusion.",
            exact["route_correction"],
            "No cumulative bound, Abel gap, winding cap, contact exclusion, Lambda<=0, PF-infinity, RH, or prize conclusion is proved.",
        ),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    return f"""# Residual Placement, C1 Transfer, and Cumulative Handoff

Date: 2026-07-31

Status: exact residual-placement and C1 handoff correction. This is not a
proof of an aggregate cross-current bound, Abel gap, `Lambda<=0`, or RH.

## Result

The complete normalized first jet is

```text
{exact['canonical_whole_jet']}
```

The endpoint-terminal edge is a component of `V_J`, not a canonical
component of the global residual.  The earlier request for a complete
`Xi` edge with four residual defects is therefore retired as a proof
obligation.  The retained edge sign and cutoff margins are unchanged.

## Allocation Gauge

For an edge row `v=(c,d)` and an allocated residual row `e=(u,w)`, put
`K(v)=det(v,partial_x v)`.  Exact expansion gives

```text
{symbolic['allocation_current']}
```

At `lambda=1` this is exactly

```text
{symbolic['four_defect_specialization']}
```

However, assigning `lambda e` to the edge and `(1-lambda)e` to the
complement leaves the complete Xi jet unchanged.  The individual edge
current varies with `lambda`; the explicit sign-changing witness is

```text
{symbolic['allocation_sign_witness']}
```

The four-defect algebra is correct but its proposed Xi-edge meaning
is not invariant.

## Derivative Order

If the second residual coordinate is `w=partial_x r_[1]/L`, then

```text
{symbolic['normalized_residual_derivative']}
```

An allocated component current consequently asks for a second x
derivative.  The invariant boundary transfer does not.  It uses

```text
{exact['c1_box']}
{exact['boundary_transfer']}
```

No derivative of `V_r` along the boundary homotopy occurs.

## Cutoff Covariance

On an adjacent overlap,

```text
{exact['chart_covariance']}
```

The retained edge connector has `Delta_cut<-h/25<0`, while the whole-main
chart jump has `M_L(Delta J)<1/5`.  These are distinct joins with distinct
roles.  Neither requires an Xi-level edge allocation.

## Cumulative Current

Inside the retained main,

```text
{exact['canonical_main_decomposition']}
```

and exact polarization gives

```text
{symbolic['cumulative_polarization']}
```

The edge theorem and the interior theorem provide the diagonal reserve

```text
{exact['diagonal_reserve']}
```

but all cross currents remain.  The generic guard

```text
{symbolic['clockwise_sum_guard']}
```

shows why strict component orientation alone cannot sign the sum.

## Next Theorem

```text
{exact['open_invariant_target']}
```

No cumulative cross-current bound, Abel gap, successor
winding cap, contact exclusion, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion.
"""


def build_payload() -> dict:
    payloads = load_sources()
    exact = exact_payload()
    symbolic = symbolic_certificate()
    rows = build_rows(exact, symbolic)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact residual-placement correction, C1 whole-jet boundary "
            "transfer, and cumulative-current handoff"
        ),
        "proof_boundary": (
            "This gate proves that the global first-order remainder is "
            "canonically attached to the complete jet, that an Xi-level edge "
            "allocation and its C2 current are noninvariant, and that the "
            "certified C1 remainder is the correct whole-boundary transfer. "
            "It preserves the retained-model edge and cutoff signs and gives "
            "the exact cumulative polarization. It proves no cross-current "
            "bound, Abel gap, winding cap, contact exclusion, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
        "source_audit": source_audit(payloads),
        "exact": exact,
        "symbolic_certificate": symbolic,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "global_remainder_coordinates": 2,
            "allocation_gauge_parameters": 1,
            "allocation_dependent_edge_current_guards": 1,
            "required_c2_for_boundary_transfer": 0,
            "required_c2_for_allocated_component_current": 1,
            "whole_jet_c1_homotopies": 1,
            "retained_model_signed_edge_families": 1,
            "retained_model_signed_interior_families": 1,
            "exact_cumulative_polarizations": 1,
            "open_cross_current_or_abel_targets": 1,
            "xi_level_edge_signs": 0,
            "contact_exclusions": 0,
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "built residual-placement C1 cumulative handoff gate: "
        f"{len(payload['rows'])} rows, 1 allocation guard, "
        "0 C2 boundary requirements, 1 cumulative polarization"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
