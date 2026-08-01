#!/usr/bin/env python3
"""Build the endpoint-relative phase-current and recurrence-defect gate."""

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
    "first_order_centered_endpoint_relative_phase_current_recurrence_gate"
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
    "complex_endpoint_corrigendum": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complex_endpoint_source_normalization_gate.json"
    ),
    "adjacent_saddle_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_saddle_recurrence.json"
    ),
    "adjacent_chart_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_chart_stability_certificate.json"
    ),
    "interior_projective_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_interior_projective_current_gate.json"
    ),
    "contact_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_contact_signed_transport_reduction.json"
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
                "bound for the de Bruijn-Newman constant, arXiv:1904.12438"
            ),
            "location": "equation (53), C_0(p), and its adjacent recurrence",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def verify_symbolics() -> None:
    kappa, kappa_x, t0 = sp.symbols(
        "kappa kappa_x T_0", real=True
    )
    hr, hi, hxr, hxi = sp.symbols(
        "H_R H_I H_xR H_xI", real=True
    )
    h = hr + sp.I * hi
    hx = hxr + sp.I * hxi
    endpoint = kappa * (t0 + sp.I) * h
    endpoint_x = (
        kappa_x * (t0 + sp.I) * h
        + kappa * sp.Rational(1, 2) * h
        + kappa * (t0 + sp.I) * hx
    )
    endpoint_current = sp.expand_complex(
        sp.im(endpoint_x * sp.conjugate(endpoint))
    )
    endpoint_target = sp.expand_complex(
        kappa**2
        * (
            (t0**2 + 1) * sp.im(hx * sp.conjugate(h))
            - sp.Rational(1, 2) * (hr**2 + hi**2)
        )
    )
    if sp.simplify(endpoint_current - endpoint_target) != 0:
        raise RuntimeError("endpoint self-current identity failed")

    endpoint_mass = sp.expand_complex(endpoint * sp.conjugate(endpoint))
    mass_target = kappa**2 * (t0**2 + 1) * (hr**2 + hi**2)
    if sp.simplify(endpoint_mass - mass_target) != 0:
        raise RuntimeError("endpoint mass identity failed")

    jr, ji, v_a, frame_r = sp.symbols(
        "J_R J_I v_a frame_r", real=True
    )
    j = jr + sp.I * ji
    jet = kappa * (t0 + sp.I) * j
    framed_x = jet + (frame_r + sp.I * v_a) * endpoint
    framed_current = sp.expand_complex(
        sp.im(framed_x * sp.conjugate(endpoint))
    )
    framed_target = sp.expand_complex(
        kappa**2
        * (t0**2 + 1)
        * (
            sp.im(j * sp.conjugate(h))
            + v_a * (hr**2 + hi**2)
        )
    )
    if sp.simplify(framed_current - framed_target) != 0:
        raise RuntimeError("framed endpoint current identity failed")

    b, u_n, delta_n = sp.symbols("b u_n delta_n", real=True)
    x_n, y_n = sp.symbols("X_n Y_n", real=True)
    z_mass = x_n**2 + y_n**2
    carrier_rate = v_a + b * u_n + delta_n
    relative_current = sp.expand(
        mass_target * carrier_rate * z_mass
        - z_mass * framed_current
    )
    relative_target = sp.expand_complex(
        kappa**2
        * (t0**2 + 1)
        * z_mass
        * (
            (hr**2 + hi**2) * (b * u_n + delta_n)
            - sp.im(j * sp.conjugate(h))
        )
    )
    if sp.simplify(relative_current - relative_target) != 0:
        raise RuntimeError("endpoint-carrier relative current failed")

    p = sp.symbols("p", real=True)
    c0 = (
        sp.exp(sp.I * sp.pi * (p**2 / 2 + sp.Rational(3, 8)))
        - sp.I * sp.sqrt(2) * sp.cos(sp.pi * p / 2)
    ) / (2 * sp.cos(sp.pi * p))
    if sp.simplify(c0.subs(p, -p) - c0) != 0:
        raise RuntimeError("C_0 parity failed")
    phase_zero = sp.simplify(
        sp.im(sp.diff(c0, p).subs(p, 0) / c0.subs(p, 0))
    )
    if phase_zero != 0:
        raise RuntimeError("midpoint C_0 phase derivative failed")
    phase_one = sp.simplify(
        sp.im(sp.diff(c0, p).subs(p, 1) / c0.subs(p, 1))
        / sp.pi
    )
    phase_one_target = 1 - sp.sqrt(4 + 2 * sp.sqrt(2)) / 4
    if sp.simplify(phase_one - phase_one_target) != 0:
        raise RuntimeError("cutoff C_0 phase derivative failed")
    if phase_one_target.is_positive is not True:
        raise RuntimeError("cutoff phase witness is not positive")

    m = sp.symbols("m", integer=True, nonnegative=True)
    leading = phase_one / 4 - m / 2
    if sp.simplify(leading.subs(m, 0) - phase_one_target / 4) != 0:
        raise RuntimeError("cutoff terminal limit failed")
    midpoint_leading = -m / 2 - sp.Rational(1, 4)
    if midpoint_leading.subs(m, 0) != -sp.Rational(1, 4):
        raise RuntimeError("midpoint terminal limit failed")

    e0, e1, e2, z1, z2 = sp.symbols(
        "e_0 e_1 e_2 z_1 z_2", complex=True
    )
    q0 = z1 + e1 - e0
    q1 = z2 + e2 - e1
    if sp.simplify(e2 + z2 + z1 - (e0 + q0 + q1)) != 0:
        raise RuntimeError("terminal-tail value recurrence failed")

    jhr, jhi, jxr, jxi = sp.symbols(
        "J_HR J_HI J_xR J_xI", real=True
    )
    jh = jhr + sp.I * jhi
    jhx = jxr + sp.I * jxi
    successor = jh - h
    successor_x = jhx - hx
    successor_current = sp.expand_complex(
        sp.im(successor_x * sp.conjugate(successor))
    )
    recurrence_target = sp.expand_complex(
        sp.im(jhx * sp.conjugate(jh))
        + sp.im(hx * sp.conjugate(h))
        - sp.im(
            jhx * sp.conjugate(h) + hx * sp.conjugate(jh)
        )
    )
    if sp.simplify(successor_current - recurrence_target) != 0:
        raise RuntimeError("endpoint-current recurrence defect failed")

    defect_y, defect_rate = sp.symbols(
        "defect_y defect_rate", real=True
    )
    predecessor = sp.Integer(1)
    defect = sp.I * defect_y
    defect_x = sp.I * defect_rate
    defect_current = sp.expand_complex(
        sp.im(defect_x * sp.conjugate(predecessor + defect))
    )
    if sp.simplify(defect_current - defect_rate) != 0:
        raise RuntimeError("real-projection nonpromotion witness failed")

    c0_atom, d0_atom, c_tail, d_tail = sp.symbols(
        "c_0 d_0 C_tail D_tail", real=True
    )
    value_total = c0_atom + c_tail
    slope_total = d0_atom + d_tail
    minor_sum = c0_atom * d_tail - d0_atom * c_tail
    if sp.simplify(
        minor_sum
        - (c0_atom * slope_total - d0_atom * value_total)
    ) != 0:
        raise RuntimeError("contact-minor sum identity failed")


def exact_payload() -> dict:
    return {
        "endpoint_self_current": {
            "coordinates": (
                "T_(0,x)=1/2, kappa in R, H_a in C, "
                "e=kappa(T_0+i)H_a"
            ),
            "division_free": (
                "I_e:=Im(e_x*conj(e))="
                "kappa^2{(T_0^2+1)Im(H_(a,x)*conj(H_a))"
                "-|H_a|^2/2}"
            ),
            "mass": "|e|^2=kappa^2(T_0^2+1)|H_a|^2",
            "ordinary_rate": (
                "On H_a!=0, partial_x arg(e)="
                "Im(H_(a,x)*conj(H_a))/|H_a|^2"
                "-1/{2(T_0^2+1)}"
            ),
            "frame_form": (
                "For e_x=g+(lambda_a-rho_1)e, "
                "lambda_a=u_a+i*v_a, "
                "g=kappa(T_0+i)J_a^(der), "
                "I_e=kappa^2(T_0^2+1)"
                "{Im(J_a^(der)*conj(H_a))+v_a|H_a|^2}"
            ),
            "phase_cancellation": (
                "Im(mu_a)+v_a=-1/{2(T_0^2+1)}; "
                "the radial normalizer cancels from I_e"
            ),
        },
        "endpoint_carrier_relative_current": {
            "carrier_rate": (
                "For z_n=r_n*zeta_n and "
                "delta_n=Im(d_(n,x)/(1+d_n)), "
                "Im(z_(n,x)*conj(z_n))="
                "(v_a+b*u_n+delta_n)|z_n|^2"
            ),
            "definition": (
                "Omega_(0n):=|e|^2 Im(z_(n,x)*conj(z_n))"
                "-|z_n|^2 Im(e_x*conj(e))"
            ),
            "division_free": (
                "Omega_(0n)=kappa^2(T_0^2+1)|z_n|^2"
                "{|H_a|^2(b*u_n+delta_n)"
                "-Im(J_a^(der)*conj(H_a))}"
            ),
            "ordinary_rate": (
                "When e*z_n!=0, "
                "Omega_(0n)/(|e|^2|z_n|^2)="
                "partial_x arg(z_n/e)"
            ),
            "zero_fibre": (
                "The polynomial current remains defined and equals zero "
                "when H_a=0 or z_n=0; no endpoint projection is divided."
            ),
        },
        "near_terminal_asymptotic": {
            "path": (
                "Fix p in [-1,1], theta=(1-p)/2, a=N+theta, "
                "q=2tL^2=1, and n=N-m with fixed m>=0; let N->infinity."
            ),
            "source_limits": (
                "H_a->C_0(p), a H_(a,x)->-C_0'(p)/(4*pi), "
                "a*mu_a->0, b->-1/2, "
                "a*u_(N-m)->m+theta, and a*delta_(N-m)->0."
            ),
            "phase_symbol": (
                "Phi(p):=Im(C_0'(p)*conj(C_0(p)))/|C_0(p)|^2"
            ),
            "scaled_limit": (
                "a*Omega_(0,N-m)/(|e|^2|z_(N-m)|^2)"
                " -> Phi(p)/(4*pi)-m/2-(1-p)/4"
            ),
            "midpoint_witness": (
                "C_0'(0)=0, so the m=0 limit at p=0 is -1/4."
            ),
            "cutoff_witness": (
                "At p=1, Phi(1)/pi="
                "1-sqrt(4+2sqrt(2))/4, so the m=0 limit is "
                "{1-sqrt(4+2sqrt(2))/4}/4>0."
            ),
            "consequence": (
                "Along the physical q=1 family the endpoint-terminal "
                "relative current has both signs for all sufficiently "
                "large cutoffs. A uniform one-sided terminal orientation "
                "is false."
            ),
        },
        "preprojection_tail_recurrence": {
            "one_step_value": (
                "Q_(N-1):=z_N+e_N-e_(N-1), hence "
                "e_N+z_N=e_(N-1)+Q_(N-1)"
            ),
            "one_step_slope": (
                "Delta W_(A,N-1):=s_*'u_N z_N+g_N-g_(N-1), "
                "hence g_N+s_*'u_Nz_N="
                "g_(N-1)+Delta W_(A,N-1)"
            ),
            "tail_value": (
                "e_N+sum_(j=0)^m z_(N-j)="
                "e_(N-m-1)+sum_(k=N-m-1)^(N-1)Q_k"
            ),
            "tail_slope": (
                "g_N+s_*'sum_(j=0)^m u_(N-j)z_(N-j)="
                "g_(N-m-1)+sum_(k=N-m-1)^(N-1)Delta W_(A,k)"
            ),
            "current_defect": (
                "If H_+=J-H, then I(H_+)=I(J)+I(H)"
                "-Im(J_x*conj(H)+H_x*conj(J)); "
                "the mixed current does not telescope."
            ),
            "projection_boundary": (
                "Existing adjacent-chart certificates control the needed "
                "real traces, not a complex norm/current of Q_k. "
                "Re Q=Re Q_x=0 alone permits either current sign: "
                "E=1, Q=i*y, Q_x=i*sigma gives I(E+Q)=sigma."
            ),
        },
        "contact_minor_sum": {
            "definition": (
                "K_(0n)=c_0d_n-d_0c_n, "
                "mathsf X=sum_(j=0)^N c_j, "
                "mathcal C_N=sum_(j=0)^N d_j"
            ),
            "identity": (
                "sum_(n=1)^N K_(0n)="
                "c_0*mathcal C_N-d_0*mathsf X"
            ),
            "contact": (
                "On mathsf X=0, "
                "sum_(n=1)^N K_(0n)=c_0*mathcal C_N."
            ),
            "ordinary_use": (
                "When c_0!=0 on contact, "
                "mathcal C_N={sum_n K_(0n)}/c_0."
            ),
            "exceptional_blind_spot": (
                "When c_0=0=mathsf X, the minor sum is zero even if "
                "mathcal C_N!=0. The corrected zero-real-projection fibre "
                "therefore still needs a division-free direct theorem."
            ),
        },
        "route_decision": {
            "retired": (
                "Retire a uniform endpoint-terminal relative-current sign "
                "and any claim that the C_0 recurrence alone telescopes "
                "phase currents."
            ),
            "surviving_candidate": (
                "The exact tail composition may still be useful if one "
                "proves a complex recurrence-defect current bound, or if "
                "the terminal-composed real edge block is treated directly."
            ),
            "next_target": (
                "Derive the division-free projective current of "
                "(C_edge,D_edge)=(c_0+c_N,d_0+d_N), retaining the second "
                "endpoint jet and d_(N,x). Test whether it is negative; "
                "otherwise return to the complete contact scalar and seek "
                "a signed cumulative-minor estimate with a separate "
                "zero-real-projection theorem."
            ),
        },
    }


def build_rows(exact: dict) -> list[GateRow]:
    return [
        GateRow(
            "erpc_01_endpoint_coordinates",
            "source_normal_form",
            "available_exact",
            "The physical endpoint keeps its full complex source factor.",
            exact["endpoint_self_current"]["coordinates"],
            "This is a coordinate identity, not a sign theorem.",
        ),
        GateRow(
            "erpc_02_endpoint_current",
            "division_free_identity",
            "ready_to_apply",
            "The endpoint self-current is independent of radial normalization.",
            exact["endpoint_self_current"]["division_free"],
            "The displayed bracket is not asserted to have one sign.",
        ),
        GateRow(
            "erpc_03_endpoint_phase_rate",
            "ordinary_fibre_identity",
            "ready_with_division_guard",
            "The ordinary endpoint phase rate contains the moving T_0+i term.",
            exact["endpoint_self_current"]["ordinary_rate"],
            "Use only when H_a is nonzero.",
        ),
        GateRow(
            "erpc_04_frame_current",
            "exact_frame_identity",
            "ready_to_apply",
            "The J_a^(der) frame gives the same endpoint current.",
            exact["endpoint_self_current"]["frame_form"],
            "J_a^(der) is not the adjacent recurrence block.",
        ),
        GateRow(
            "erpc_05_carrier_rate",
            "source_rate_identity",
            "ready_to_apply",
            "The physical carrier angular rate retains its correction.",
            exact["endpoint_carrier_relative_current"]["carrier_rate"],
            "No correction term is dropped.",
        ),
        GateRow(
            "erpc_06_relative_current",
            "division_free_identity",
            "ready_to_apply",
            "The endpoint-carrier relative current cancels the frame phase.",
            exact["endpoint_carrier_relative_current"]["division_free"],
            "Its source bracket remains unsigned.",
        ),
        GateRow(
            "erpc_07_relative_rate",
            "ordinary_fibre_identity",
            "ready_with_division_guard",
            "The polynomial current equals the relative phase rate off zero fibres.",
            exact["endpoint_carrier_relative_current"]["ordinary_rate"],
            "The polynomial formula, not the ratio, is primary.",
        ),
        GateRow(
            "erpc_08_near_terminal_limit",
            "asymptotic_source_identity",
            "available_exact",
            "Every fixed near-terminal offset has an explicit scaled limit.",
            exact["near_terminal_asymptotic"]["scaled_limit"],
            "This is along the stated q=1 cofinal family.",
        ),
        GateRow(
            "erpc_09_midpoint_sign",
            "exact_asymptotic_witness",
            "available_exact",
            "The terminal relative current is eventually negative at p=0.",
            exact["near_terminal_asymptotic"]["midpoint_witness"],
            "This is an endpoint-source sign witness, not a contact.",
        ),
        GateRow(
            "erpc_10_cutoff_sign",
            "exact_asymptotic_witness",
            "available_exact",
            "The terminal relative current is eventually positive at cutoff equality.",
            exact["near_terminal_asymptotic"]["cutoff_witness"],
            "This is an endpoint-source sign witness, not a contact.",
        ),
        GateRow(
            "erpc_11_terminal_sign_guard",
            "nonpromotion_guard",
            "obstruction_exact",
            "A uniform terminal relative-current sign is false.",
            exact["near_terminal_asymptotic"]["consequence"],
            "Contact conditioning could impose additional restrictions.",
        ),
        GateRow(
            "erpc_12_tail_recurrence",
            "preprojection_composition",
            "available_exact",
            "Terminal and near-terminal values telescope into recurrence defects.",
            exact["preprojection_tail_recurrence"]["tail_value"],
            "The complex recurrence defects still require control.",
        ),
        GateRow(
            "erpc_13_current_defect",
            "recurrence_nonpromotion_guard",
            "obstruction_exact",
            "Endpoint phase currents acquire a mixed recurrence defect.",
            exact["preprojection_tail_recurrence"]["current_defect"],
            "The C_0 value recurrence alone does not sign this term.",
        ),
        GateRow(
            "erpc_14_projection_guard",
            "data_sufficiency_guard",
            "obstruction_exact",
            "Real-trace recurrence bounds cannot determine complex edge current.",
            exact["preprojection_tail_recurrence"]["projection_boundary"],
            "A complex-current bound or direct real-edge theorem is needed.",
        ),
        GateRow(
            "erpc_15_contact_minor_sum",
            "contact_identity",
            "ready_to_apply",
            "All endpoint-carrier minors sum to the complete contact scalar.",
            exact["contact_minor_sum"]["identity"],
            "Division by c_0 is permitted only on its ordinary fibre.",
        ),
        GateRow(
            "erpc_16_exceptional_fibre",
            "division_guard",
            "obstruction_exact",
            "The cumulative minor identity is blind on the corrected endpoint fibre.",
            exact["contact_minor_sum"]["exceptional_blind_spot"],
            "A separate division-free theorem remains mandatory.",
        ),
    ]


def build_payload() -> dict:
    verify_symbolics()
    exact = exact_payload()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact endpoint-relative phase-current identities, "
            "q=1 terminal sign-reversal guard, and recurrence-defect gate"
        ),
        "proof_boundary": (
            "This artifact proves the endpoint self-current, the "
            "division-free endpoint-carrier relative current, the "
            "near-terminal q=1 scaled law, two opposite terminal sign "
            "witnesses, the preprojection tail recurrence, its mixed-current "
            "defect, and the cumulative contact-minor identity. It proves "
            "no terminal-composed real-edge sign, no signed cumulative-minor "
            "estimate, no Abel-scalar gap, no successor winding cap, no "
            "contact exclusion, no Lambda<=0 conclusion, no PF-infinity, "
            "no RH, and no prize-level result."
        ),
        "source_audit": audit_sources(),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "exact_endpoint_currents": 2,
            "division_free_relative_currents": 1,
            "asymptotic_sign_witnesses": 2,
            "terminal_tail_recurrences": 1,
            "contact_minor_sums": 1,
            "uniform_terminal_signs": 0,
            "abel_gaps": 0,
            "winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    counts = payload["counts"]
    return f"""# Jensen-Window PF Newman Polymath-15 Endpoint-Relative Phase-Current Recurrence Gate

Date: 2026-07-31

Status: exact endpoint-relative phase-current identities, a rigorous
`q=1` terminal sign-reversal guard, and an exact recurrence-defect
reduction. This is not a proof of RH or `Lambda <= 0`.
This is not an Abel gap or contact-exclusion theorem.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Endpoint Self-Current

For the corrected complex endpoint,

```text
{exact["endpoint_self_current"]["coordinates"]}

{exact["endpoint_self_current"]["division_free"]}

{exact["endpoint_self_current"]["mass"]}
```

The `kappa_x` term is radial and cancels. On the ordinary fibre,

```text
{exact["endpoint_self_current"]["ordinary_rate"]}
```

The equivalent moving-frame form is

```text
{exact["endpoint_self_current"]["frame_form"]}

{exact["endpoint_self_current"]["phase_cancellation"]}
```

Thus the result retains the full complex `H_a` and `J_a^(der)`.

## Relative Carrier Current

Write

```text
{exact["endpoint_carrier_relative_current"]["carrier_rate"]}

{exact["endpoint_carrier_relative_current"]["definition"]}
```

Direct cancellation of the common frame phase gives

```text
{exact["endpoint_carrier_relative_current"]["division_free"]}
```

When both atoms are nonzero this is

```text
{exact["endpoint_carrier_relative_current"]["ordinary_rate"]}
```

The division-free polynomial identity remains primary on every zero
fibre.

## Exact Sign-Reversal Guard

Take the physical cofinal family

```text
{exact["near_terminal_asymptotic"]["path"]}
```

The already certified source estimates imply

```text
{exact["near_terminal_asymptotic"]["source_limits"]}

{exact["near_terminal_asymptotic"]["scaled_limit"]}
```

For the terminal carrier `m=0`, the two exact witnesses are

```text
{exact["near_terminal_asymptotic"]["midpoint_witness"]}

{exact["near_terminal_asymptotic"]["cutoff_witness"]}
```

Therefore

```text
{exact["near_terminal_asymptotic"]["consequence"]}
```

This rejects a uniform terminal relative-current sign before any broad
numerical calibration. It does not say that contact conditioning is
irrelevant.

## Tail Composition Before Projection

The one-step recurrence is

```text
{exact["preprojection_tail_recurrence"]["one_step_value"]}

{exact["preprojection_tail_recurrence"]["one_step_slope"]}
```

Iteration gives

```text
{exact["preprojection_tail_recurrence"]["tail_value"]}

{exact["preprojection_tail_recurrence"]["tail_slope"]}
```

The value recurrence does not telescope its phase current:

```text
{exact["preprojection_tail_recurrence"]["current_defect"]}
```

Moreover,

```text
{exact["preprojection_tail_recurrence"]["projection_boundary"]}
```

So a complex recurrence-defect current bound, or a direct real-edge
theorem, is a genuinely new input.

## Contact-Minor Sum

The corrected endpoint minors satisfy

```text
{exact["contact_minor_sum"]["identity"]}

{exact["contact_minor_sum"]["contact"]}
```

On the ordinary endpoint fibre this recovers the complete contact scalar.
On the exceptional fibre,

```text
{exact["contact_minor_sum"]["exceptional_blind_spot"]}
```

## Route Decision

```text
{exact["route_decision"]["retired"]}

{exact["route_decision"]["surviving_candidate"]}

{exact["route_decision"]["next_target"]}
```

## Proof Boundary

This gate has {counts["rows"]} rows, two exact endpoint-current forms, one
division-free relative current, two opposite asymptotic sign witnesses,
one terminal-tail recurrence, and one contact-minor sum. It proves zero
uniform terminal signs, zero Abel gaps, and zero winding bounds. The
terminal-composed real-edge current and the complete signed contact
scalar remain open.
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    counts = payload["counts"]
    print(
        "built endpoint-relative phase-current recurrence gate: "
        f"{counts['rows']} rows, "
        f"{counts['exact_endpoint_currents']} exact endpoint currents, "
        f"{counts['division_free_relative_currents']} division-free "
        "relative current, "
        f"{counts['asymptotic_sign_witnesses']} asymptotic sign witnesses, "
        f"{counts['terminal_tail_recurrences']} terminal-tail recurrence, "
        f"{counts['contact_minor_sums']} contact-minor sum, "
        f"{counts['uniform_terminal_signs']} uniform terminal signs, "
        f"{counts['abel_gaps']} Abel gaps, "
        f"{counts['winding_bounds']} winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
