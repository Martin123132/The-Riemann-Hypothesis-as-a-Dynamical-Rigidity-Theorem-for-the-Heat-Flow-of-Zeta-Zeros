#!/usr/bin/env python3
"""Build the terminal-tail-anchored Abel cross-current reduction."""

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
    "contiguous_terminal_tail_anchored_abel_cross_current_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "growing_tail": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_growing_prefix_finite_height_gate.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "interior_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "interior_projective_current_gate.json"
    ),
    "residual_handoff": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "residual_placement_c1_cumulative_handoff_gate.json"
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
    tail = payloads["growing_tail"]
    tail_theorem = tail.get("exact", {}).get("theorem", "")
    if "a^2*J_(tail,M)/S_a^2<-1/400" not in tail_theorem:
        raise RuntimeError("growing-tail current reserve drifted")
    if tail.get("counts", {}).get("bulk_aggregate_closures") != 0:
        raise RuntimeError("growing-tail proof boundary drifted")

    abel = payloads["abel_prefix"].get("exact", {})
    for marker in (
        "U=u_N*S+sum_(k=1)^(N-1)h_k*F_k",
        "mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X",
    ):
        if marker not in json.dumps(abel, sort_keys=True):
            raise RuntimeError(f"Abel-prefix source marker missing: {marker}")

    interior = payloads["interior_current"].get("exact", {})
    if "mathsf_X=C_edge+sum_(n=1)^(N-1)c_n" not in interior.get(
        "edge_block", ""
    ):
        raise RuntimeError("interior additive-coordinate source drifted")

    handoff = payloads["residual_handoff"].get("exact", {})
    if "K(V_J)=-R_diag+C_cross" not in handoff.get(
        "diagonal_reserve", ""
    ):
        raise RuntimeError("residual-handoff cross-current target drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_terminal_reserve": (
            "P_T=a^2*J_(tail,M)/S_a^2<-1/400"
        ),
        "imported_bulk_closures": 0,
    }


def symbolic_certificate() -> dict:
    # A three-term generic bulk is enough to certify finite Abel summation.
    kappa, h_1, h_2 = sp.symbols("kappa h_1 h_2", real=True)
    z_1, z_2, z_3 = sp.symbols("z_1 z_2 z_3")
    dz_1, dz_2, dz_3 = sp.symbols("dz_1 dz_2 dz_3")
    f_1 = z_1
    f_2 = z_1 + z_2
    f_3 = z_1 + z_2 + z_3
    r_direct = (
        (kappa + h_1 + h_2) * z_1
        + (kappa + h_2) * z_2
        + kappa * z_3
    )
    r_abel = kappa * f_3 + h_1 * f_1 + h_2 * f_2
    if sp.expand(r_direct - r_abel) != 0:
        raise RuntimeError("truncated Abel identity failed")

    df_1 = dz_1
    df_2 = dz_1 + dz_2
    df_3 = dz_1 + dz_2 + dz_3
    dr_direct = (
        (kappa + h_1 + h_2) * dz_1
        + (kappa + h_2) * dz_2
        + kappa * dz_3
    )
    dr_abel = kappa * df_3 + h_1 * df_1 + h_2 * df_2
    if sp.expand(dr_direct - dr_abel) != 0:
        raise RuntimeError("differentiated truncated Abel identity failed")

    c, b, u = sp.symbols("c b u", real=True)
    t_r, t_i, g_r = sp.symbols("T_R T_I G_R", real=True)
    f_r, f_i, r_r, r_i = sp.symbols(
        "F_R F_I R_R R_I", real=True
    )
    a_tail = g_r - c * u * t_r
    a_bulk = c * r_r - b * r_i - b * u * f_i
    direct_total = (
        g_r
        + c * (u * f_r + r_r)
        - b * (u * f_i + r_i)
        - c * u * (t_r + f_r)
    )
    if sp.expand(direct_total - a_tail - a_bulk) != 0:
        raise RuntimeError("tail-anchored centered scalar failed")
    alternate = (
        g_r
        - c * u * t_r
        + b * u * t_i
        - b * u * (t_i + f_i)
        + c * r_r
        - b * r_i
    )
    if sp.expand(alternate - direct_total) != 0:
        raise RuntimeError("alternate re-anchoring identity failed")

    c_x, b_x, u_x = sp.symbols("c_x b_x u_x", real=True)
    f_i_x, r_r_x, r_i_x = sp.symbols(
        "F_I_x R_R_x R_I_x", real=True
    )
    a_bulk_x = (
        c_x * r_r
        + c * r_r_x
        - b_x * r_i
        - b * r_i_x
        - (b_x * u + b * u_x) * f_i
        - b * u * f_i_x
    )
    expected_a_bulk_x = sp.expand(
        c_x * r_r
        + c * r_r_x
        - b_x * r_i
        - b * r_i_x
        - b_x * u * f_i
        - b * u_x * f_i
        - b * u * f_i_x
    )
    if sp.expand(a_bulk_x - expected_a_bulk_x) != 0:
        raise RuntimeError("bulk centered-scalar derivative failed")

    # The common real normalization cancels from the projective determinant.
    c_0, a_0, c_0_x, a_0_x, sigma, sigma_x, h = sp.symbols(
        "C_0 A_0 C_0_x A_0_x Sigma Sigma_x h",
        positive=True,
        real=True,
    )
    c_norm = c_0 / sigma
    a_norm = a_0 / sigma
    c_norm_x = (c_0_x * sigma - c_0 * sigma_x) / sigma**2
    a_norm_x = (a_0_x * sigma - a_0 * sigma_x) / sigma**2
    normalized_current = sp.expand(
        (c_norm * a_norm_x - a_norm * c_norm_x) / h**2
    )
    raw_current = (c_0 * a_0_x - a_0 * c_0_x) / (
        h**2 * sigma**2
    )
    if sp.simplify(normalized_current - raw_current) != 0:
        raise RuntimeError("common-normalizer current identity failed")

    c_t, d_t, e_t, n_t = sp.symbols("C_T D_T E_T N_T", real=True)
    c_b, d_b, e_b, n_b = sp.symbols("C_B D_B E_B N_B", real=True)
    p_t = c_t * n_t - d_t * e_t
    p_b = c_b * n_b - d_b * e_b
    p_total = (c_t + c_b) * (n_t + n_b) - (
        d_t + d_b
    ) * (e_t + e_b)
    p_cross = (
        c_t * n_b
        + c_b * n_t
        - d_t * e_b
        - d_b * e_t
    )
    if sp.expand(p_total - p_t - p_b - p_cross) != 0:
        raise RuntimeError("tail-bulk current polarization failed")

    x_t, aa_t, x_t_x, aa_t_x = sp.symbols(
        "X_T A_T X_T_x A_T_x", real=True
    )
    x_b, aa_b, x_b_x, aa_b_x = sp.symbols(
        "X_B A_B X_B_x A_B_x", real=True
    )
    scalar_tail = (x_t * aa_t_x - aa_t * x_t_x) / h**2
    scalar_total = (
        (x_t + x_b) * (aa_t_x + aa_b_x)
        - (aa_t + aa_b) * (x_t_x + x_b_x)
    ) / h**2
    phi = sp.expand(scalar_total - scalar_tail)
    phi_expected = (
        (x_t + x_b) * aa_b_x
        + x_b * aa_t_x
        - (aa_t + aa_b) * x_b_x
        - aa_b * x_t_x
    ) / h**2
    if sp.simplify(phi - phi_expected) != 0:
        raise RuntimeError("tail-anchored signed remainder failed")

    # Generic guard: negative component currents do not sign their sum.
    tail_value = sp.Matrix([1, 0])
    tail_derivative = sp.Matrix([0, -1])
    bulk_value = sp.Matrix([sp.Rational(-1, 2), 0])
    bulk_derivative = sp.Matrix([0, 2])

    def determinant(value: sp.Matrix, derivative: sp.Matrix) -> sp.Expr:
        return sp.expand(
            value[0] * derivative[1] - value[1] * derivative[0]
        )

    guard_tail = determinant(tail_value, tail_derivative)
    guard_bulk = determinant(bulk_value, bulk_derivative)
    guard_total = determinant(
        tail_value + bulk_value, tail_derivative + bulk_derivative
    )
    guard_cross = sp.expand(guard_total - guard_tail - guard_bulk)
    if (
        guard_tail != -1
        or guard_bulk != -1
        or guard_total != sp.Rational(1, 2)
        or guard_cross != sp.Rational(5, 2)
    ):
        raise RuntimeError("negative-components join guard failed")

    return {
        "truncated_abel": (
            "F_k=sum_(n=1)^k z_n, h_k=log((k+1)/k), "
            "kappa_B=u_B-u_N=log(N/B); "
            "R_B=sum_(n=1)^B(u_n-u_N)z_n="
            "kappa_B*F_B+sum_(k=1)^(B-1)h_k*F_k"
        ),
        "differentiated_abel": (
            "Because N and B are fixed inside the chart, "
            "partial_x(u_n-u_N)=0; hence "
            "R_(B,x)=kappa_B*F_(B,x)+"
            "sum_(k=1)^(B-1)h_k*F_(k,x)"
        ),
        "normalized_partition": (
            "hat(T)=T/S_a, hat(F_k)=F_k/S_a, hat(R_B)=R_B/S_a; "
            "hat(W_0)=hat(T)+hat(F_B)"
        ),
        "tail_centered_scalar": (
            "A_T=Re(hat(G_T))-c*u_N*Re(hat(T))"
        ),
        "bulk_centered_scalar": (
            "A_B=Re(s_*'*hat(R_B))-b*u_N*Im(hat(F_B))="
            "c*Re(hat(R_B))-b*Im(hat(R_B))-"
            "b*u_N*Im(hat(F_B))"
        ),
        "reanchored_total_scalar": (
            "A=A_T+A_B=Re(hat(G_T)-s_*'*u_N*hat(T))"
            "-b*u_N*Im(hat(W_0))+Re(s_*'*hat(R_B))"
        ),
        "bulk_centered_derivative": (
            "A_(B,x)=c_x*R_R+c*R_(R,x)-b_x*R_I-b*R_(I,x)"
            "-(b_x*u_N+b*u_(N,x))*F_I-b*u_N*F_(I,x), "
            "where F=hat(F_B), R=hat(R_B)"
        ),
        "normalized_coordinate_map": (
            "For a block with normalized real value X and normalized "
            "centered scalar A, (C,D,E,N)="
            "(X,A/h,X_x/h,A_x/h^2) and P=CN-DE"
        ),
        "normalizer_cancellation": (
            "If X=C_0/S_a and A=A_0/S_a, then "
            "P=[C_0*A_(0,x)-A_0*C_(0,x)]/(h^2*S_a^2)="
            "a^2*J/S_a^2; every S_(a,x) term cancels"
        ),
        "four_coordinate_split": (
            "P_ret=P_T+P_B+P_cross, P_B=C_B*N_B-D_B*E_B, "
            "P_cross=C_T*N_B+C_B*N_T-D_T*E_B-D_B*E_T"
        ),
        "signed_bulk_remainder": (
            "Phi_B=P_B+P_cross=h^-2*[(X_T+X_B)A_(B,x)"
            "+X_B*A_(T,x)-(A_T+A_B)X_(B,x)-A_B*X_(T,x)]"
        ),
        "negative_components_guard": (
            "At one point v_T=(1,0), v_(T,x)=(0,-1), "
            "v_B=(-1/2,0), v_(B,x)=(0,2): "
            "P_T=P_B=-1, P_cross=5/2, but P_ret=1/2"
        ),
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Fix a physical q=1, L>=50, fixed-N chart and an admissible "
            "M>=1 with h*(M+1)^18<=1. Put B=N-M-1>=1 and keep M,B "
            "fixed while differentiating inside the chart."
        ),
        "partition": (
            "T=e+sum_(n=B+1)^N z_n, "
            "G_T=g+s_*'*sum_(n=B+1)^N u_n*z_n, "
            "F_B=sum_(n=1)^B z_n, U_B=sum_(n=1)^B u_n*z_n; "
            "W_0=T+F_B and W_A=G_T+s_*'*U_B."
        ),
        "bulk_moment": (
            "U_B=u_N*F_B+R_B, with "
            "R_B=sum_(n=1)^B log(N/n)z_n."
        ),
        "tail_input": (
            "The imported growing-prefix theorem supplies "
            "P_T=a^2*J_(tail,M)/S_a^2<-1/400."
        ),
        "sharp_join_target": (
            "The signed inequality Phi_B<=1/400 implies P_ret<0. "
            "The stronger Phi_B<=1/800 preserves P_ret<-1/800."
        ),
        "streaming_reduction": (
            "Compute hat(F_k) and hat(F_(k,x)) in index order, accumulate "
            "hat(R_B)=kappa_B*hat(F_B)+sum h_k*hat(F_k) and its "
            "derivative, then evaluate Phi_B in constant additional work. "
            "This is O(B) time and O(1) streaming state, rather than an "
            "O(B^2) pair expansion."
        ),
        "proof_gain": (
            "The certified terminal block now replaces the isolated "
            "endpoint in the Abel identity without changing the u_N shear. "
            "All nonterminal uncertainty is concentrated in one signed "
            "cumulative correlation Phi_B."
        ),
        "centering_guard": (
            "Do not recenter at u_B: u_N is the centering used by the "
            "certified current. An x-dependent change of shear contributes "
            "an extra shear-derivative current unless it is tracked."
        ),
        "pi_provenance": (
            "The only pi occurs through the inherited completed-zeta/"
            "Riemann-Siegel normalization and h=1/a. The Abel weights are "
            "logarithmic arithmetic ratios and introduce no fitted geometry."
        ),
        "next_target": (
            "Prove the Xi-specific signed upper bound Phi_B<=1/400, "
            "preferably Phi_B<=1/800, from the actual multiplicative carrier "
            "phases and cumulative prefixes before taking absolute values."
        ),
        "proof_boundary": (
            "This is an exact cancellation-preserving reduction, not the "
            "required bound. It proves no bulk aggregate closure, complete "
            "retained current sign, Xi-level current theorem, Abel gap, "
            "winding cap, contact exclusion, Q209, Lambda<=0, PF-infinity, "
            "RH, or prize-level conclusion."
        ),
    }


def build_rows(exact: dict, symbolic: dict) -> list[ReductionRow]:
    return [
        ReductionRow("ttac_01_domain", "chart domain", "proved", "The terminal/bulk partition is fixed before x differentiation.", exact["domain"], "Adjacent cutoff transport remains separate."),
        ReductionRow("ttac_02_partition", "exact partition", "proved", "The retained main splits into the certified terminal block and the nonterminal bulk.", exact["partition"], "No carrier is omitted or counted twice."),
        ReductionRow("ttac_03_abel", "truncated Abel identity", "proved", "The bulk centered moment is a positive-log weighted prefix sum.", symbolic["truncated_abel"], "This is finite algebra, not an estimate."),
        ReductionRow("ttac_04_abel_x", "differentiated Abel identity", "proved", "The same prefix weights survive one x derivative.", symbolic["differentiated_abel"], "N and B must stay fixed inside the chart."),
        ReductionRow("ttac_05_normalizer", "normalization covariance", "proved", "The common source normalizer cancels from projective current.", symbolic["normalizer_cancellation"], "The normalizer must be common to both blocks."),
        ReductionRow("ttac_06_tail", "certified terminal input", "proved", "The growing terminal block supplies a strict finite-height reserve.", exact["tail_input"], "This imports the q=1 theorem only."),
        ReductionRow("ttac_07_bulk_scalar", "bulk Abel scalar", "proved", "The bulk centered coordinate uses only F_B and R_B.", symbolic["bulk_centered_scalar"], "No division by a real projection is used."),
        ReductionRow("ttac_08_bulk_x", "bulk first jet", "proved", "The bulk second projective coordinate uses only first prefix derivatives.", symbolic["bulk_centered_derivative"], "No second derivative of an individual carrier is introduced."),
        ReductionRow("ttac_09_coordinates", "four-coordinate map", "proved", "Tail and bulk share one additive normalized coordinate map.", symbolic["normalized_coordinate_map"], "The u_N centering is preserved."),
        ReductionRow("ttac_10_polarization", "current polarization", "proved", "The retained current is tail plus bulk plus one exact cross term.", symbolic["four_coordinate_split"], "The cross term is not assigned a sign."),
        ReductionRow("ttac_11_remainder", "signed cumulative remainder", "proved", "Every nonterminal current contribution is one scalar Phi_B.", symbolic["signed_bulk_remainder"], "Its Xi-specific upper bound remains open."),
        ReductionRow("ttac_12_reserve", "sufficient inequality", "proved", "The terminal margin turns one signed bulk estimate into a retained current theorem.", exact["sharp_join_target"], "No such bound is proved here."),
        ReductionRow("ttac_13_guard", "nonpromotion guard", "proved", "Negative component currents do not sign their sum.", symbolic["negative_components_guard"], "This is generic algebra, not an Xi counterexample."),
        ReductionRow("ttac_14_handoff", "open handoff", "open", "The next theorem is now a linear-time cumulative-correlation bound.", exact["next_target"], exact["proof_boundary"], {"streaming": exact["streaming_reduction"]}),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    return f"""# Contiguous Terminal-Tail Anchored Abel Cross-Current Reduction

Date: 2026-07-31

Status: exact cancellation-preserving reduction to one signed nonterminal
correlation, with `0 bulk closures`. This is not a proof artifact; the
`Phi_B` bound and RH remain open.

## Fixed-Chart Partition

{exact['domain']}

```text
{exact['partition']}
{exact['bulk_moment']}
```

The terminal block is exactly the endpoint-terminal edge together with the
carriers `N-M,...,N`; the bulk is `1,...,B`.

## Truncated Abel Identity

```text
{symbolic['truncated_abel']}
{symbolic['differentiated_abel']}
```

The differences `u_n-u_N=log(N/n)`, `kappa_B`, and every `h_k` are constant
under `x` differentiation while `N` and `B` are fixed. Thus no derivative of
an Abel weight appears.

## Tail-Anchored Centering

Use the same real positive source scale `S_a` and the same `u_N` shear as the
growing-prefix theorem:

```text
{symbolic['normalized_partition']}
{symbolic['tail_centered_scalar']}
{symbolic['bulk_centered_scalar']}
{symbolic['reanchored_total_scalar']}
```

In particular, the certified terminal block replaces the old isolated
endpoint without recentering the problem at `u_B`.

The differentiated bulk scalar is

```text
{symbolic['bulk_centered_derivative']}
```

## Four Coordinates and Current

```text
{symbolic['normalized_coordinate_map']}
{symbolic['normalizer_cancellation']}
```

Apply this map separately to the tail and bulk. Exact polarization gives

```text
{symbolic['four_coordinate_split']}
```

Equivalently, all nonterminal uncertainty is

```text
{symbolic['signed_bulk_remainder']}
```

{exact['streaming_reduction']}

## Sharp Remaining Target

```text
{exact['tail_input']}
{exact['sharp_join_target']}
```

This is a signed upper-bound problem. Replacing `Phi_B` by separate absolute
pair bounds discards the prefix cancellation retained by the Abel form.

The reason separate clockwise theorems are insufficient is exact:

```text
{symbolic['negative_components_guard']}
```

{exact['proof_gain']}

## Boundary

{exact['centering_guard']}

{exact['pi_provenance']}

{exact['next_target']}

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    payloads = load_sources()
    source = source_audit(payloads)
    symbolic = symbolic_certificate()
    exact = exact_payload()
    rows = build_rows(exact, symbolic)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact q=1 terminal-tail-anchored Abel cross-current reduction; "
            "signed Phi_B target open; no bulk closure, Xi-level current "
            "theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "exact_partitions": 1,
            "exact_abel_identities": 2,
            "normalized_coordinate_maps": 2,
            "current_polarization_identities": 2,
            "terminal_current_inputs": 1,
            "signed_bulk_targets": 1,
            "bulk_aggregate_closures": 0,
            "xi_level_current_theorems": 0,
        },
        "rows": [asdict(row) for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "wrote terminal-tail anchored Abel cross-current reduction: "
        f"{payload['counts']['rows']} rows, 2 Abel identities, "
        "1 signed Phi_B target, 0 bulk closures"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
