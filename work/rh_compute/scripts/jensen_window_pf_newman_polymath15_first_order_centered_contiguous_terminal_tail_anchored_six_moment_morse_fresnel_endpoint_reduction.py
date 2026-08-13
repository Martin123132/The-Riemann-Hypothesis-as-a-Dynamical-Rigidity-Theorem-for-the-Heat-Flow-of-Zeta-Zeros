#!/usr/bin/env python3
"""Build the logarithmic Morse-Fresnel endpoint reduction."""

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
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "endpoint_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
VANDEHEY_URL = "https://arxiv.org/abs/1205.0090"
SOURCE_PATHS = {
    "finite_poisson_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_finite_poisson_"
        "transport_reduction.json"
    ),
    "flow_matrix_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_flow_matrix_"
        "phase_reduction.json"
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


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.simplify(expression) != 0:
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
    finite = payloads["finite_poisson_transport"]
    counts = finite.get("counts", {})
    if counts.get("bare_logarithmic_moments") != 6:
        raise RuntimeError("finite-Poisson source lost six-moment closure")
    if counts.get("explicit_uniform_remainder_constants") != 0:
        raise RuntimeError("finite-Poisson proof boundary drifted")
    if counts.get("dual_cutoff_jump_laws") != 2:
        raise RuntimeError("finite-Poisson cutoff audit drifted")
    flow = payloads["flow_matrix_phase"]
    if flow.get("counts", {}).get("real_flow_observations") != 8:
        raise RuntimeError("flow source lost eight-observation handoff")
    return {
        key: {
            "path": str(SOURCE_PATHS[key].relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(SOURCE_PATHS[key]),
        }
        for key in SOURCE_PATHS
    }


def symbolic_certificate() -> dict:
    delta = sp.symbols("delta", positive=True)
    v = 1 + delta
    h = v - 1 - sp.log(v)
    z = sp.sqrt(2 * h)
    jacobian = sp.simplify(v * z / (v - 1))
    z_series = sp.series(z, delta, 0, 4).removeO()
    jacobian_series = sp.series(jacobian, delta, 0, 3).removeO()
    require_zero(
        z_series - (delta - delta**2 / 3 + 7 * delta**3 / 36),
        "signed Morse-coordinate series failed",
    )
    require_zero(
        jacobian_series - (1 + 2 * delta / 3 - 5 * delta**2 / 36),
        "Morse Jacobian series failed",
    )

    u, alpha, r = sp.symbols("u alpha r", positive=True)
    amplitude = sp.Function("A")(u)
    q = alpha / u - r
    first = sp.diff(amplitude / q, u)
    second = sp.diff(first / q, u)
    q_prime = sp.diff(q, u)
    q_second = sp.diff(q, u, 2)
    expanded = (
        sp.diff(amplitude, u, 2) / q**2
        - 3 * sp.diff(amplitude, u) * q_prime / q**3
        - amplitude * q_second / q**3
        + 3 * amplitude * q_prime**2 / q**4
    )
    require_zero(second - expanded, "twice-integrated tail operator failed")

    u0 = alpha / r
    vv = sp.symbols("v", positive=True)
    phase = alpha * sp.log(u0 * vv) - r * u0 * vv
    phase0 = alpha * (sp.log(u0) - 1)
    require_zero(
        phase - (phase0 - alpha * (vv - 1 - sp.log(vv))),
        "exact logarithmic Morse phase failed",
    )
    return {
        "h_prime": str(sp.diff(vv - 1 - sp.log(vv), vv)),
        "z_series_right": str(z_series),
        "jacobian_series_right": str(jacobian_series),
        "phase_difference": "0",
        "tail_operator_difference": "0",
        "q_prime": str(q_prime),
        "q_second": str(q_second),
        "tail_operator": str(expanded),
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use the physical q=1, L>=50, fixed-N chart, B=N-M-1, "
            "alpha=xi/(2*pi), and alpha_-<=alpha<=alpha_+ corresponding to "
            "T_0-epsilon<=xi<=T_0. The six amplitudes are A_j(u)=(log "
            "u)^j exp[t(log u)^2/4-sigma log u], 0<=j<=5."
        ),
        "source_scope": (
            "Vandehey Section 11 motivates an incomplete-Fresnel endpoint "
            "chart and explicitly notes that it is computer-evaluable near "
            "transition, but its displayed remainder and Corollary 1.4 retain "
            "implicit constants. The formulas below are instead exact direct "
            "specializations to the logarithmic phase; no sourced big-O term "
            "is promoted into a numerical bound."
        ),
        "fixed_roster": (
            "Set m=max(1,floor(alpha_-/(2B))), n=ceil(2alpha_+), and "
            "T={m,...,n}. This roster is fixed on the whole auxiliary segment. "
            "It contains every mode whose saddle can enter [1,B], as well as "
            "both endpoint-transition collars."
        ),
        "proof_boundary": (
            "This proves an exact fixed-roster Morse-Fresnel representation "
            "for every transition mode, an exact twice-integrated reassembly "
            "of the infinite nonstationary tail, stable digamma/Hurwitz-zeta "
            "endpoint sums, and explicit constant-one residual majorants. It "
            "does not evaluate the summed transition variation norms on the "
            "physical chart, bound the imaginary Hermitian coefficient, prove "
            "a signed offset-sum or signed flow estimate, upper-bound Phi_B, "
            "exclude contact, prove a retained aggregate or Xi theorem, Q209, "
            "a cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or a "
            "prize-level conclusion."
        ),
        "next_target": (
            "Build interval majorants for the six transition quantities "
            "|c(y_1)|+|c(y_B)|+integral|c'|, sum them over the fixed roster "
            "without splitting the Hermitian pairs, and evaluate the explicit "
            "Hurwitz-zeta tail envelope. Then compose both physical endpoint "
            "functionals with the terminal recurrence and insert the six value "
            "transforms into the eight observation rows."
        ),
        "exact_target": (
            "Prove Psi(T_0)-epsilon*integral_0^1 "
            "Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800."
        ),
    }


def morse_certificate() -> dict:
    return {
        "phase": (
            "For r>0 put u_0=alpha/r, v=u/u_0, "
            "h(v)=v-1-log v, and phi_r(u)=alpha log u-ru. Then "
            "phi_r(u)=phi_r(u_0)-alpha h(v), where "
            "phi_r(u_0)=alpha[log(alpha/r)-1]."
        ),
        "coordinate": (
            "Define z(v)=sgn(v-1)sqrt(2h(v)), z(1)=0. Since h>=0, "
            "h'= (v-1)/v, z'=(v-1)/(vz)>0, and z tends to -infinity "
            "and +infinity at the two ends, z is a smooth increasing "
            "bijection (0,infinity)->R."
        ),
        "jacobian": (
            "With J(v)=dv/dz=vz/(v-1) and J(1)=1, put y=sqrt(alpha)z. "
            "Then du/dy=[sqrt(alpha)/r]J(v) and "
            "phi_r(u)=phi_r(u_0)-y^2/2 exactly."
        ),
        "series": (
            "For delta=v-1, z=delta-delta^2/3+7delta^3/36+O(delta^4) "
            "and J=1+2delta/3-5delta^2/36+O(delta^3). These expansions "
            "supply the removable values at v=1; they are not asymptotic "
            "approximations in the exact transformed integral."
        ),
        "endpoint_coordinates": (
            "The transformed endpoints are y_1=sqrt(alpha)z(r/alpha) "
            "and y_B=sqrt(alpha)z(Br/alpha). The lower saddle crossing "
            "r=alpha has y_1=0; the upper crossing r=alpha/B has y_B=0."
        ),
    }


def transition_certificate() -> dict:
    return {
        "amplitude": (
            "Let v(y) be the inverse Morse coordinate and define "
            "b_(j,r)(y)=[sqrt(alpha)/r]A_j(u_0v(y))J(v(y)). Then "
            "b_(j,r)(0)=[sqrt(alpha)/r]A_j(u_0)."
        ),
        "exact_integral": (
            "For every r in the fixed roster, "
            "I_(j,r)=e(phi_r(u_0))*integral_(y_1)^(y_B) "
            "b_(j,r)(y)e(-y^2/2)dy. This identity remains valid when "
            "u_0 is just outside [1,B]."
        ),
        "uniform_main": (
            "Define F(y_1,y_B)=integral_(y_1)^(y_B)e(-y^2/2)dy and "
            "U_(j,r)=b_(j,r)(0)e(phi_r(u_0))F(y_1,y_B). Then "
            "I_(j,r)=U_(j,r)+Q_(j,r), where Q is the same integral with "
            "b(y)-b(0)."
        ),
        "fresnel_limits": (
            "F(-infinity,infinity)=e(-1/8)=exp(-i*pi/4), while each "
            "half-line integral ending at zero is e(-1/8)/2. Thus the "
            "usual full stationary main and both starred half-mains are "
            "limits of one continuous function, not separate cutoff rules."
        ),
        "residual_bound": (
            "Put c(y)=[b(y)-b(0)]/y for y!=0 and c(0)=b'(0). Exact "
            "integration by parts gives |Q_(j,r)|<={|c(y_1)|+|c(y_B)|+"
            "integral_(y_1)^(y_B)|c'(y)|dy}/(2*pi). The constant is "
            "exactly 1/(2*pi), with no hidden big-O constant."
        ),
        "removable_c": (
            "The identities c(y)=integral_0^1 b'(sy)ds and "
            "c'(y)=integral_0^1 s b''(sy)ds make the residual bound "
            "regular at y=0. Moreover c_(j,r)(0)=r^(-1){jA_(j-1)(u_0)+"
            "[2/3-sigma+(t/2)log u_0]A_j(u_0)}."
        ),
        "smoothness": (
            "Because the roster is fixed and z, J, F, and c use removable "
            "charts, U+Q is smooth through y_1=0 and y_B=0. The two sharp "
            "jump laws of Section 11.167 are replaced by one exact incomplete-"
            "Fresnel transition; no cutoff absence is assumed."
        ),
    }


def tail_certificate() -> dict:
    return {
        "gap": (
            "For r outside T, q_r(u)=alpha/u-r has no zero on [1,B]. "
            "The low tail lies below alpha_-/(2B), and the high tail lies "
            "above 2alpha_+, so both are uniformly separated from every "
            "stationary point on the full xi segment."
        ),
        "two_ibp": (
            "Let kappa=1/(2*pi*i), D_rA=(A/q_r)', and "
            "C_rA=((D_rA)/q_r)'. Then exactly "
            "I_r=kappa[e(phi_r)A/q_r]_1^B-kappa^2[e(phi_r)(D_rA)/q_r]_1^B"
            "+kappa^2 integral_1^B e(phi_r)C_rA du."
        ),
        "operator": (
            "Writing q'= -alpha/u^2 and q''=2alpha/u^3, "
            "C_rA=A''/q^2-3A'q'/q^3-Aq''/q^3+3A(q')^2/q^4."
        ),
        "stable_sums": (
            "For x in (m-1,n+1), a_x=x-m+1 and b_x=n+1-x, the symmetric "
            "outside-roster sums are S_1(x)=psi(b_x)-psi(a_x), "
            "S_2(x)=zeta(2,a_x)+zeta(2,b_x), and "
            "S_3(x)=zeta(3,a_x)-zeta(3,b_x). Equivalently "
            "S_1=pi*cot(pi*x)-sum_(r=m)^n 1/(x-r), with removable values "
            "at every integer in T."
        ),
        "boundary_functional": (
            "At an integer endpoint mu in {1,B}, x=alpha/mu and "
            "e(-rmu)=1. The complete outside-roster endpoint functional is "
            "B_j(mu)=e(alpha log mu){kappa A_j(mu)S_1(x)-kappa^2["
            "A_j'(mu)S_2(x)+(alpha/mu^2)A_j(mu)S_3(x)]}. The summed tail "
            "equals B_j(B)-B_j(1)+kappa^2 sum_(r outside T)integral e C_rA_j."
        ),
        "absolute_envelope": (
            "For k>=2 put Z_k(x)=zeta(k,a_x)+zeta(k,b_x). The remaining "
            "tail has the explicit bound 1/(4*pi^2) integral_1^B {"
            "|A_j''|Z_2+[3|A_j'||q'|+|A_j||q''|]Z_3+"
            "3|A_j||q'|^2Z_4}du, where x=alpha/u."
        ),
        "convergence": (
            "The only conditional tail term is S_1, which is evaluated in "
            "the same symmetric sense as full Poisson summation. After that "
            "endpoint term is reassembled, S_2, S_3, and the integral "
            "remainder are absolutely convergent. No logarithmic truncation "
            "loss or artificial dyadic endpoint remains."
        ),
    }


def handoff_certificate() -> dict:
    return {
        "six_values": (
            "Apply the exact representation independently to j=0,...,5. "
            "Section 11.167 already proves that the eight physical flow "
            "observations are fixed rows of these six values, so no derivative "
            "of the incomplete-Fresnel main is needed and no seventh moment "
            "is introduced."
        ),
        "route_decision": (
            "Use the fixed-roster Morse-Fresnel plus twice-integrated tail "
            "decomposition as the primary endpoint route. Do not use the "
            "inexplicit Vandehey remainder as a numerical certificate, and "
            "do not introduce a hard dyadic partition unless a later explicit "
            "comparison proves it smaller after all boundary terms are rejoined."
        ),
        "remaining_scalar_work": (
            "The new quantitative objects are explicit: six transition "
            "variation sums and six Hurwitz-zeta tail integrals. Their physical "
            "interval evaluation, Hermitian pairing, and endpoint-terminal "
            "composition remain open."
        ),
        "pi_provenance": (
            "The factor 2*pi comes from e(x)=exp(2*pi*i*x); the Fresnel limit "
            "e(-1/8)=exp(-i*pi/4) is the negative-curvature signature. The "
            "cotangent identity is the symmetric partial-fraction expansion "
            "of pi*cot(pi*x). No geometric circle or fitted pi is introduced."
        ),
    }


def build_rows(
    exact: dict,
    morse: dict,
    transition: dict,
    tail: dict,
    handoff: dict,
) -> list[ReductionRow]:
    return [
        ReductionRow("mfer_01_domain", "physical domain", "proved", "The auxiliary chart and six amplitudes are fixed.", exact["domain"], "No remainder bound is inferred from the domain."),
        ReductionRow("mfer_02_source", "source boundary", "proved", "The sourced motivation is separated from the exact specialization.", exact["source_scope"], "No implicit source constant is promoted."),
        ReductionRow("mfer_03_phase", "Morse phase", "proved", "The logarithmic phase has an exact universal defect.", morse["phase"], "Per-mode identity."),
        ReductionRow("mfer_04_coordinate", "global coordinate", "proved", "The signed coordinate is a global smooth bijection.", morse["coordinate"], "No numerical asymptotic is used."),
        ReductionRow("mfer_05_jacobian", "Jacobian", "proved", "The phase is exactly quadratic in y.", morse["jacobian"], "Exact for every positive mode."),
        ReductionRow("mfer_06_series", "removable chart", "proved", "The saddle values are removable.", morse["series"], "Series records local values only."),
        ReductionRow("mfer_07_endpoints", "endpoint coordinates", "proved", "Both dual crossings are zeros of transformed endpoints.", morse["endpoint_coordinates"], "No crossing is assumed absent."),
        ReductionRow("mfer_08_roster", "fixed roster", "proved", "One finite roster covers the full auxiliary segment.", exact["fixed_roster"], "The roster is not yet optimized."),
        ReductionRow("mfer_09_amplitude", "transformed amplitude", "proved", "The saddle coefficient is explicit.", transition["amplitude"], "Value transform only."),
        ReductionRow("mfer_10_integral", "exact transition integral", "proved", "Every roster mode has an exact Fresnel chart.", transition["exact_integral"], "No stationary approximation is claimed."),
        ReductionRow("mfer_11_main", "uniform main", "proved", "The transition main and residual are exact.", transition["uniform_main"], "Residual still requires evaluation."),
        ReductionRow("mfer_12_limits", "Fresnel limits", "proved", "Full and half stationary mains are one continuous family.", transition["fresnel_limits"], "Standard Fresnel normalization."),
        ReductionRow("mfer_13_residual", "transition residual bound", "proved", "The residual has an explicit constant-one variation bound.", transition["residual_bound"], "The variation norm is not yet bounded physically."),
        ReductionRow("mfer_14_removable_c", "residual chart", "proved", "The variation integrand is regular at the saddle.", transition["removable_c"], "No interval enclosure is claimed."),
        ReductionRow("mfer_15_smooth", "cutoff smoothness", "proved", "The exact fixed-roster representation has no sharp cutoff jump.", transition["smoothness"], "Global signed estimates remain open."),
        ReductionRow("mfer_16_gap", "nonstationary gap", "proved", "Every outside-roster mode is uniformly nonstationary.", tail["gap"], "Uses the deliberately enlarged roster."),
        ReductionRow("mfer_17_ibp", "twofold integration by parts", "proved", "The nonstationary tail has an exact second-order identity.", tail["two_ibp"], "Endpoint terms must be summed before absolute values."),
        ReductionRow("mfer_18_operator", "tail operator", "proved", "The remaining integrand has denominator powers at least two.", tail["operator"], "Exact differential identity."),
        ReductionRow("mfer_19_sums", "stable endpoint sums", "proved", "The conditional endpoint series has a stable closed form.", tail["stable_sums"], "Symmetric Poisson convention retained."),
        ReductionRow("mfer_20_boundary", "endpoint functional", "proved", "All outside-roster first and second endpoint terms are reassembled.", tail["boundary_functional"], "Must still compose with the physical terminal recurrence."),
        ReductionRow("mfer_21_envelope", "absolute tail envelope", "proved", "The post-boundary tail has an explicit convergent majorant.", tail["absolute_envelope"], "No numerical physical value is claimed."),
        ReductionRow("mfer_22_convergence", "convergence audit", "proved", "No logarithmic truncation loss remains after reassembly.", tail["convergence"], "The finite transition sum can still be large."),
        ReductionRow("mfer_23_values", "six-value handoff", "proved", "Six value transforms still suffice for eight observations.", handoff["six_values"], "No signed flow estimate follows."),
        ReductionRow("mfer_24_route", "route selection", "proved", "The constant-bearing endpoint route is selected.", handoff["route_decision"], "Selection is architectural, not the final bound."),
        ReductionRow("mfer_25_handoff", "quantitative handoff", "open", "The remaining scalar quantities are explicit.", handoff["remaining_scalar_work"], exact["next_target"]),
        ReductionRow("mfer_26_pi", "pi provenance", "proved", "Every pi has analytic provenance.", handoff["pi_provenance"], "No fitted constant."),
        ReductionRow("mfer_27_boundary", "proof boundary", "proved", "Exact reduction is separated from the missing physical estimates.", exact["proof_boundary"], "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    morse = payload["morse_certificate"]
    transition = payload["transition_certificate"]
    tail = payload["tail_certificate"]
    handoff = payload["handoff_certificate"]
    return f"""# Six-Moment Logarithmic Morse-Fresnel Endpoint Reduction

Date: 2026-08-01

Status: exact endpoint and nonstationary-tail reduction with explicit
constant-one residual formulas, `0 evaluated physical remainder constants`,
`0 signed flow bounds`, and this is not a proof of RH.

Primary analytic source used only for comparison and motivation:

```text
{VANDEHEY_URL}
```

## Domain And Source Boundary

{exact['domain']}

{exact['source_scope']}

## Exact Logarithmic Morse Coordinate

{morse['phase']}

{morse['coordinate']}

{morse['jacobian']}

{morse['series']}

{morse['endpoint_coordinates']}

## Fixed Transition Roster

{exact['fixed_roster']}

## Incomplete-Fresnel Transition

{transition['amplitude']}

{transition['exact_integral']}

{transition['uniform_main']}

{transition['fresnel_limits']}

{transition['residual_bound']}

{transition['removable_c']}

{transition['smoothness']}

## Reassembled Nonstationary Tail

{tail['gap']}

{tail['two_ibp']}

{tail['operator']}

{tail['stable_sums']}

{tail['boundary_functional']}

{tail['absolute_envelope']}

{tail['convergence']}

## Six-Value Handoff

{handoff['six_values']}

{handoff['route_decision']}

{handoff['remaining_scalar_work']}

## Pi Provenance

{handoff['pi_provenance']}

## Next Target

{exact['next_target']}

{exact['exact_target']}

## Proof Boundary

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    sources = load_sources()
    exact = exact_payload()
    morse = morse_certificate()
    transition = transition_certificate()
    tail = tail_certificate()
    handoff = handoff_certificate()
    symbolic = symbolic_certificate()
    rows = build_rows(exact, morse, transition, tail, handoff)
    return {
        "kind": STEM,
        "date": "2026-08-01",
        "status": (
            "exact fixed-roster logarithmic Morse-Fresnel transition and "
            "twice-integrated nonstationary-tail reduction complete; physical "
            "variation sums and signed flow estimate open"
        ),
        "source_audit": source_audit(sources),
        "external_source": {
            "author": "Joseph Vandehey",
            "title": "Error term improvements for van der Corput transforms",
            "url": VANDEHEY_URL,
            "usage": "comparison and Fresnel-endpoint motivation only",
        },
        "exact": exact,
        "morse_certificate": morse,
        "transition_certificate": transition,
        "tail_certificate": tail,
        "handoff_certificate": handoff,
        "symbolic_certificate": symbolic,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "bare_logarithmic_moments": 6,
            "fixed_transition_rosters": 1,
            "exact_morse_phase_identities": 1,
            "exact_transition_integrals": 6,
            "incomplete_fresnel_main_families": 6,
            "explicit_transition_residual_majorants": 6,
            "twofold_nonstationary_ibp_identities": 6,
            "stable_endpoint_sum_families": 3,
            "explicit_tail_integral_majorants": 6,
            "sharp_cutoff_jump_laws_required": 0,
            "hidden_big_o_terms_in_derived_decomposition": 0,
            "evaluated_physical_remainder_constants": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
        "proof_boundary": exact["proof_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    atomic_write(args.output, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(payload))
    counts = payload["counts"]
    print(
        "built six-moment Morse-Fresnel endpoint reduction: "
        f"{counts['rows']} rows, "
        f"{counts['exact_transition_integrals']} exact transition integrals, "
        f"{counts['explicit_transition_residual_majorants']} transition majorants, "
        f"{counts['explicit_tail_integral_majorants']} tail majorants, "
        f"{counts['evaluated_physical_remainder_constants']} evaluated physical constants, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
