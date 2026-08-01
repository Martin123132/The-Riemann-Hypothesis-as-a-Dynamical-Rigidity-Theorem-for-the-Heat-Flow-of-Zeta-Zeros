#!/usr/bin/env python3
"""Build the q=1 saddle-phase quadratic transfer barrier."""

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
    "contiguous_terminal_tail_anchored_five_moment_q1_saddle_"
    "phase_quadratic_transfer_barrier"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "saddle_variance": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_q1_saddle_"
        "phase_variance_reduction.json"
    ),
    "finite_edge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "q1_finite_height_real_edge_remainder_gate.json"
    ),
    "growing_tail": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_growing_prefix_finite_height_gate.json"
    ),
    "two_carrier": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_two_carrier_"
        "kernel_reduction.json"
    ),
    "mangoldt": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_current_mangoldt_"
        "normal_form_gate.json"
    ),
}


@dataclass(frozen=True)
class QuadraticBarrierRow:
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


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.simplify(sp.expand(expression)) != 0:
        raise RuntimeError(label)


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    saddle = payloads["saddle_variance"]
    physical = saddle.get("physical_q1_certificate", {})
    variance = saddle.get("phase_variance_certificate", {})
    if "Omega=-Im(s_*)" not in physical.get("frequency", ""):
        raise RuntimeError("physical saddle frequency source drifted")
    if "D_hat" not in variance.get("branch_free", ""):
        raise RuntimeError("branch-free discriminant source drifted")
    if saddle.get("counts", {}).get("signed_physical_bounds") != 0:
        raise RuntimeError("upstream physical proof boundary drifted")

    edge = payloads["finite_edge"].get("majorant_certificate", {}).get(
        "source_block_majorants", {}
    )
    stable = edge.get("stable_terminal", {})
    if stable.get("W") != "|W|<6h":
        raise RuntimeError("terminal logarithm source drifted")

    growing = payloads["growing_tail"]
    if "K=M+1" not in growing.get("exact", {}).get("prefix", ""):
        raise RuntimeError("growing-tail cutoff source drifted")
    normalized = growing.get("symbolic_certificate", {}).get(
        "normalized_carrier", ""
    )
    if "W_theta" not in normalized or "z_(N-m)/S_a" not in normalized:
        raise RuntimeError("terminal normalized-carrier source drifted")

    terminal = payloads["two_carrier"].get("finite_sum_certificate", {})
    if "C_Nw_N=-Q(p)exp(W_theta)" not in terminal.get(
        "physical_chi_identity", ""
    ):
        raise RuntimeError("terminal correction-free carrier source drifted")
    if "|d_n|<1/2" not in terminal.get("physical_chi_bound", ""):
        raise RuntimeError("uniform correction source drifted")

    correction = payloads["mangoldt"].get("symbolic_certificate", {}).get(
        "quadratic_correction", ""
    )
    if "c_n=(1+d_n)/(1+d_1)" not in correction:
        raise RuntimeError("correction multiplier source drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_signed_physical_bounds": 0,
        "imported_phi_b_bounds": 0,
    }


def discriminant_collapse_certificate() -> dict:
    f, g, p, q, r, s, u_x = sp.symbols(
        "f g p q r s u_x", real=True
    )
    m0 = f + sp.I * g
    m1 = p + sp.I * q
    m2 = r + sp.I * s
    r2 = sp.expand(m0 * sp.conjugate(m0))
    a_coord = sp.expand(sp.re(m1 * sp.conjugate(m0)))
    c_coord = sp.expand((m2 * m0 - m1**2) * sp.conjugate(m0) ** 2)
    d_hat = sp.expand(
        f**2 * sp.re(c_coord)
        + a_coord**2 * r2
        - (2 * u_x * r2**2 + sp.im(c_coord)) * f * g
    )
    four_p = -f * r - q**2 + 2 * u_x * f * g
    require_zero(d_hat + four_p * r2**2, "D_hat collapse failed")

    return {
        "primary_current": (
            "For M_0=f+ig, M_1=p+iq, and M_2=r+is, "
            "4P_bulk^(0)=-f*r-q^2+2u_x*f*g."
        ),
        "collapse": (
            "The Section 11.162 discriminant factors identically as "
            "D_hat=(f*r+q^2-2u_x*f*g)|M_0|^4="
            "-4P_bulk^(0)|M_0|^4."
        ),
        "interpretation": (
            "D_hat is the primary quadratic current multiplied by "
            "|M_0|^4, not an independent sixth-degree positivity invariant."
        ),
        "zero_fibre": (
            "At M_0=0, D_hat=0 while "
            "P_bulk^(0)=-(Im M_1)^2/4. The discriminant therefore loses "
            "the strict exceptional-fibre information carried by P_bulk^(0)."
        ),
        "symbolic_difference": str(
            sp.simplify(d_hat + four_p * r2**2)
        ),
    }


def frequency_flow_certificate() -> dict:
    f, g, p, q, r, y, u_x = sp.symbols(
        "f g p q r y u_x", real=True
    )
    current = (-f * r - q**2 + 2 * u_x * f * g) / 4
    flow = {f: q, g: -p, q: -r, r: y}
    derivative = sum(sp.diff(current, variable) * value for variable, value in flow.items())
    derivative_expected = (
        q * r - f * y + 2 * u_x * (q * g - f * p)
    ) / 4
    require_zero(
        derivative - derivative_expected,
        "quadratic frequency derivative failed",
    )

    df, dg, dq, dr = sp.symbols("df dg dq dr", real=True)
    shifted = (
        -(f + df) * (r + dr)
        - (q + dq) ** 2
        + 2 * u_x * (f + df) * (g + dg)
    ) / 4
    endpoint_expected = (
        -df * r
        - f * dr
        - df * dr
        - 2 * q * dq
        - dq**2
        + 2 * u_x * (df * g + f * dg + df * dg)
    )
    require_zero(
        4 * (shifted - current) - endpoint_expected,
        "endpoint perturbation identity failed",
    )

    return {
        "moment_flow": (
            "For M_k(xi)=omega_a H_k(xi), "
            "partial_xi M_k=-i M_(k+1), 0<=k<=4. Hence, if "
            "M_0=f+ig, M_1=p+iq, M_2=r+is, and y=Im M_3, then "
            "f_xi=q, g_xi=-p, q_xi=-r, and r_xi=y."
        ),
        "current_derivative": (
            "4 partial_xi P_bulk^(0)="
            "q*r-f*y+2u_x(q*g-f*p)."
        ),
        "signed_line_integral": (
            "Since Omega=T_0-epsilon, "
            "4[P_bulk^(0)(Omega)-P_bulk^(0)(T_0)]="
            "-epsilon integral_0^1 K(T_0-theta*epsilon)dtheta, where "
            "K=q*r-f*y+2u_x(q*g-f*p)."
        ),
        "moment_line_integral": (
            "H_k(Omega)-H_k(T_0)="
            "i*epsilon*integral_0^1 H_(k+1)(T_0-theta*epsilon)dtheta."
        ),
        "endpoint_identity": (
            "For delta f,delta g,delta q,delta r, the exact endpoint "
            "difference is 4 delta P=-delta f*r-f*delta r-"
            "delta f*delta r-2q*delta q-(delta q)^2+2u_x["
            "delta f*g+f*delta g+delta f*delta g]."
        ),
        "absolute_line_majorant": (
            "Put A_j^+=sum_(n<=B)u_n^j A_n. Pointwise triangle "
            "inequality on the exact line integral gives "
            "|P(Omega)-P(T_0)|<=B_abs, where "
            "B_abs=epsilon[A_0^+A_3^++A_1^+A_2^+"
            "+4u_x A_0^+A_1^+]/4."
        ),
        "symbolic_differences": {
            "frequency_derivative": str(
                sp.simplify(derivative - derivative_expected)
            ),
            "endpoint_perturbation": str(
                sp.simplify(4 * (shifted - current) - endpoint_expected)
            ),
        },
    }


def terminal_amplitude_certificate() -> dict:
    h_max = sp.Rational(1, 72_000_000_000)
    log_two_lower = sp.Rational(69, 100)
    if not 8 * h_max < log_two_lower:
        raise RuntimeError("terminal exponential lower bound failed")

    return {
        "terminal_identity": (
            "At n=N, C_N w_N=-Q_RS(p)exp(W_theta), "
            "|Q_RS|=1, |W_theta|<6h, and "
            "C_N=(1+d_N)/(1+d_1)."
        ),
        "correction_bound": (
            "The uniform |d_j|<1/2 gives |C_N|<3, hence "
            "A_N=|w_N|>exp(-6h)/3."
        ),
        "profile_transfer": (
            "Since A_N=A_a exp[B_a u_N+u_N^2/(8L^2)], "
            "B_a<1/2, u_N<2h, and L>=50, the exponent is <2h. "
            "Therefore A_a>exp(-8h)/3>1/6."
        ),
        "rational_certificate": {
            "h_upper": str(h_max),
            "log_2_lower": str(log_two_lower),
            "8h_upper": str(8 * h_max),
            "8h_lt_log_2_lower": bool(8 * h_max < log_two_lower),
            "A_a_lower": "1/6",
        },
    }


def positive_mass_barrier_certificate() -> dict:
    ratio_lower = (
        sp.Integer(2) ** 50
        * sp.Rational(69, 100) ** 3
        / (
            96
            * sp.Rational(22, 7)
            * sp.Integer(50) ** 2
        )
    )
    if not ratio_lower > 490_000_000:
        raise RuntimeError("absolute-majorant gap certificate failed")

    coefficient = sp.simplify(
        sp.Rational(3, 3200) * sp.Rational(1, 6) ** 2
    )
    if coefficient != sp.Rational(1, 38_400):
        raise RuntimeError("barrier coefficient simplification failed")

    return {
        "bulk_inclusion": (
            "Let K=M+1, B=N-K, and "
            "I={n integer: ceil(a/4)<=n<=floor(a/2)}. From K^18<=a "
            "and a>72000000000 one has K<a/4 and B>a/2, so I is "
            "contained in {1,...,B}."
        ),
        "block_count": (
            "The block has #I>=a/4-1>a/5. For n in I, "
            "log(2)<=u_n=log(a/n)<=log(4) and A_n>=A_a."
        ),
        "mass_lower_bounds": (
            "Writing d=log(2), A_0^+>a*A_a/5 and "
            "A_3^+>a*A_a*d^3/5."
        ),
        "epsilon_conversion": (
            "Because a^2=x/(4pi)+1/(32L^2), x<4pi a^2. Thus "
            "epsilon>3/(8L^2x)>3/(32pi L^2a^2)."
        ),
        "majorant_lower_bound": (
            "The positive A_0^+A_3^+ term alone forces "
            "B_abs>3A_a^2 log(2)^3/(3200pi L^2)>"
            "log(2)^3/(38400pi L^2)."
        ),
        "target_gap": (
            "Since h=1/a, a^2>exp(L), exp(L)/L^2 is increasing for "
            "L>=50, exp(1)>2, log(2)>69/100, and pi<22/7, "
            "B_abs/(h^2/400)>490000000."
        ),
        "scope": (
            "This lower-bounds the stated positive-mass upper majorant. "
            "It does not lower-bound the true signed transfer error. It "
            "rules out spending an h^2/400 reserve through this canonical "
            "pointwise triangle estimate, not through cancellation-preserving "
            "analysis of the signed line integral."
        ),
        "rational_certificate": {
            "majorant_coefficient": str(coefficient),
            "ratio_lower": str(ratio_lower),
            "ratio_lower_decimal": str(sp.N(ratio_lower, 30)),
            "ratio_gt_490000000": bool(ratio_lower > 490_000_000),
        },
    }


def exact_payload() -> dict:
    return {
        "rejected_route": (
            "Do not estimate D_hat as a separate positive invariant and do "
            "not transfer P_bulk^(0) from T_0 to Omega by replacing every "
            "moment with its positive mass. The first step is algebraically "
            "redundant and the second has a certified majorant more than "
            "490000000 times the entire h^2/400 reserve."
        ),
        "surviving_route": (
            "Retain K(xi)=q*r-f*Im(M_3)+2u_x[q*g-f*p] with its signs on "
            "the whole interval T_0-epsilon<=xi<=T_0. Insert its moments "
            "into the endpoint-composed Hermitian/transpose or balanced "
            "Mangoldt/Poisson representation before taking absolute values."
        ),
        "sharp_signed_target": (
            "A P-only benchmark is "
            "|integral_0^1 K(T_0-theta*epsilon)dtheta|"
            "<=h^2/(100epsilon), which would give "
            "|P(Omega)-P(T_0)|<=h^2/400. The actual proof may instead "
            "transfer E+P_bulk^(0)+R_corr jointly and use cancellation "
            "between those composed terms."
        ),
        "proof_boundary": (
            "This proves the algebraic collapse of D_hat, exact frequency "
            "flow and signed transfer identities, A_a>1/6, a macroscopic "
            "positive-mass block, and failure of one explicit absolute-mass "
            "transfer majorant at the required scale. It does not prove the "
            "true saddle transfer is large, rule out oscillatory or composed "
            "transfer, sign P_bulk^(0) at T_0 or Omega, bound E or R_corr, "
            "upper-bound Phi_B, establish a signed Type-I/II, Vaughan, or "
            "Poisson estimate, exclude contact, transfer the retained model "
            "to Xi, prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level "
            "result."
        ),
    }


def build_rows(
    collapse: dict,
    flow: dict,
    terminal: dict,
    barrier: dict,
    exact: dict,
) -> list[QuadraticBarrierRow]:
    return [
        QuadraticBarrierRow("qpb_01_domain", "physical q=1 domain", "proved", "The transfer audit stays on the certified q=1 terminal/bulk chart.", "Use L>=50, h=1/a<1/72000000000, B=N-M-1, Omega=T_0-epsilon, and 3/(8L^2x)<epsilon<7/(8L^2x).", "No q>1 statement is made."),
        QuadraticBarrierRow("qpb_02_current", "primary quadratic current", "proved", "The ideal bulk current is a quadratic in four real moment components.", collapse["primary_current"], "No relative moment is used."),
        QuadraticBarrierRow("qpb_03_collapse", "discriminant collapse", "proved", "The branch-free discriminant has no independent positivity content.", collapse["collapse"] + " " + collapse["interpretation"], "The factorization is polynomial and division-free.", {"difference": collapse["symbolic_difference"]}),
        QuadraticBarrierRow("qpb_04_zero", "exceptional fibre", "proved", "The primary current retains information lost by D_hat at M_0=0.", collapse["zero_fibre"], "No division by M_0 is allowed."),
        QuadraticBarrierRow("qpb_05_flow", "moment frequency flow", "proved", "All five unnormalized moments obey one exact differential chain.", flow["moment_flow"] + " " + flow["moment_line_integral"], "This identity preserves oscillatory cancellation."),
        QuadraticBarrierRow("qpb_06_derivative", "quadratic frequency derivative", "proved", "The leading current has an exact signed derivative using moments only through order three.", flow["current_derivative"], "u_x is fixed along the auxiliary frequency segment.", flow["symbolic_differences"]),
        QuadraticBarrierRow("qpb_07_line", "signed transfer line integral", "proved", "The exact T_0-to-Omega transfer is one signed integral.", flow["signed_line_integral"], "No termwise absolute value has been taken."),
        QuadraticBarrierRow("qpb_08_endpoint", "endpoint perturbation audit", "proved", "Direct endpoint expansion agrees with the differential transfer.", flow["endpoint_identity"], "The line-integral form is the sharper absolute-value starting point."),
        QuadraticBarrierRow("qpb_09_terminal", "terminal amplitude normalization", "proved", "The physical terminal source gives a uniform lower bound for A_a.", terminal["terminal_identity"] + " " + terminal["correction_bound"] + " " + terminal["profile_transfer"], "This uses the correction-free amplitude and retains C_N explicitly.", terminal["rational_certificate"]),
        QuadraticBarrierRow("qpb_10_block", "macroscopic bulk block", "proved", "A fixed macroscopic block lies below every admissible bulk cutoff.", barrier["bulk_inclusion"] + " " + barrier["block_count"], "The short growing terminal prefix does not remove this block."),
        QuadraticBarrierRow("qpb_11_mass", "positive-moment lower bounds", "proved", "The positive masses in the canonical transfer bound are macroscopic.", barrier["mass_lower_bounds"], "These are lower bounds on positive masses, not signed moments."),
        QuadraticBarrierRow("qpb_12_majorant", "canonical absolute majorant", "proved", "Pointwise triangle inequality gives a valid but scale-incompatible transfer bound.", flow["absolute_line_majorant"] + " " + barrier["epsilon_conversion"] + " " + barrier["majorant_lower_bound"], "The statement concerns the size of the majorant, not the actual error."),
        QuadraticBarrierRow("qpb_13_barrier", "reserve-scale barrier", "guard_validated", "The canonical positive-mass majorant cannot fit inside the available reserve.", barrier["target_gap"], barrier["scope"], barrier["rational_certificate"]),
        QuadraticBarrierRow("qpb_14_reject", "route rejection", "guard_validated", "Two tempting saddle-transfer moves are now formally rejected.", exact["rejected_route"], "Other cancellation-preserving routes remain open."),
        QuadraticBarrierRow("qpb_15_handoff", "signed analytic handoff", "open", "The next estimate must preserve the signed frequency-flow integrand or transfer the full composed current.", exact["surviving_route"] + " " + exact["sharp_signed_target"], exact["proof_boundary"]),
    ]


def render_note(payload: dict) -> str:
    collapse = payload["discriminant_collapse_certificate"]
    flow = payload["frequency_flow_certificate"]
    terminal = payload["terminal_amplitude_certificate"]
    barrier = payload["positive_mass_barrier_certificate"]
    exact = payload["exact"]
    return f"""# Physical q=1 Quadratic Transfer Barrier

Date: 2026-08-01

Status: exact algebraic reduction and route guard. This proves no physical
current sign and is not a proof of RH.

## Discriminant Collapse

```text
{collapse['primary_current']}

{collapse['collapse']}
```

{collapse['interpretation']} {collapse['zero_fibre']}

## Exact Frequency Flow

```text
{flow['moment_flow']}

{flow['moment_line_integral']}

{flow['current_derivative']}

{flow['signed_line_integral']}
```

The independent endpoint expansion is

```text
{flow['endpoint_identity']}
```

## Terminal Scale

```text
{terminal['terminal_identity']}
{terminal['correction_bound']}
{terminal['profile_transfer']}
```

## Positive-Mass Barrier

```text
{barrier['bulk_inclusion']}
{barrier['block_count']}
{barrier['mass_lower_bounds']}

{flow['absolute_line_majorant']}
{barrier['epsilon_conversion']}
{barrier['majorant_lower_bound']}
{barrier['target_gap']}
```

{barrier['scope']}

## Route Decision

{exact['rejected_route']}

{exact['surviving_route']}

```text
{exact['sharp_signed_target']}
```

## Pi Provenance

The `pi` in `x<4pi*a^2`, `u_x=h^2/(8pi)`, and `T_0=2pi*a^2`
is inherited from the completed-zeta, Riemann--Siegel, and Poisson
normalizations already traced in Formal Core Section 11.118. No circle,
polygon, or fitted geometric constant is introduced here.

## Boundary

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
    collapse = discriminant_collapse_certificate()
    flow = frequency_flow_certificate()
    terminal = terminal_amplitude_certificate()
    barrier = positive_mass_barrier_certificate()
    exact = exact_payload()
    rows = build_rows(collapse, flow, terminal, barrier, exact)
    return {
        "kind": STEM,
        "date": "2026-08-01",
        "status": (
            "exact q=1 discriminant collapse, signed frequency transfer, "
            "and absolute-mass majorant barrier; cancellation-preserving "
            "physical estimate open; no Phi_B bound, Lambda<=0, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "discriminant_collapse_certificate": collapse,
        "frequency_flow_certificate": flow,
        "terminal_amplitude_certificate": terminal,
        "positive_mass_barrier_certificate": barrier,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "discriminant_collapses": 1,
            "moment_frequency_flow_identities": 5,
            "signed_line_integrals": 1,
            "endpoint_perturbation_identities": 1,
            "terminal_amplitude_lower_bounds": 1,
            "positive_mass_block_bounds": 2,
            "absolute_mass_transfer_barriers": 1,
            "signed_physical_bounds": 0,
            "phi_b_bounds": 0,
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
        "built q=1 quadratic transfer barrier: "
        "15 rows, 1 discriminant collapse, 1 signed line integral, "
        "1 terminal amplitude lower bound, 1 absolute-mass barrier, "
        "0 signed physical bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
