#!/usr/bin/env python3
"""Build the six-moment endpoint-composition and retention gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "endpoint_composition_retention_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "endpoint_coherence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "endpoint_coherence_abel_gate.json"
    ),
    "morse_fresnel_endpoint": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "endpoint_reduction.json"
    ),
    "flow_matrix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_flow_matrix_"
        "phase_reduction.json"
    ),
}


@dataclass(frozen=True)
class CompositionRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
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
    coherence = payloads["endpoint_coherence"].get("counts", {})
    if coherence.get("endpoint_phase_collapses") != 2:
        raise RuntimeError("endpoint-coherence source lost two phase collapses")
    if coherence.get("endpoint_composed_bounds") != 0:
        raise RuntimeError("endpoint-coherence source overpromoted composition")
    morse = payloads["morse_fresnel_endpoint"].get("counts", {})
    if morse.get("exact_transition_integrals") != 6:
        raise RuntimeError("Morse-Fresnel source lost six transition integrals")
    if morse.get("stable_endpoint_sum_families") != 3:
        raise RuntimeError("Morse-Fresnel source lost stable endpoint sums")
    flow = payloads["flow_matrix"].get("counts", {})
    if flow.get("correction_free_moments") != 6:
        raise RuntimeError("flow source lost six moments")
    if flow.get("real_flow_observations") != 8:
        raise RuntimeError("flow source lost eight observations")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_audit() -> dict:
    delta, e_b, a_b = sp.symbols("delta e_B A_B")
    kappa = sp.symbols("kappa")
    cal_1, cal_b = sp.symbols("calB_1 calB_B")
    c_1, c_b = sp.symbols("c_1 c_B")
    fresnel, interior, tail = sp.symbols("F J T")

    direct = (
        delta / 2
        + e_b * a_b / 2
        + fresnel
        + kappa * c_1
        - kappa * e_b * c_b
        + interior
        + e_b * cal_b
        - cal_1
        + tail
    )
    lower = delta / 2 - cal_1 + kappa * c_1
    upper = a_b / 2 + cal_b - kappa * c_b
    composed = fresnel + lower + e_b * upper + interior + tail
    if sp.expand(direct - composed) != 0:
        raise RuntimeError("endpoint composition failed")

    p, rho, m, ell = sp.symbols("p rho m ell")
    full_quadratic = (p + rho) ** 2 * m + ell * (p + rho)
    expanded = p**2 * m + ell * p + 2 * p * m * rho + rho**2 * m + ell * rho
    if sp.expand(full_quadratic - expanded) != 0:
        raise RuntimeError("observation expansion failed")

    return {
        "composition_check": "direct endpoint pieces equal the collected form exactly",
        "observation_check": (
            "For symmetric M, (p+r)^T M(p+r)+ell^T(p+r) equals "
            "p^T M p+ell^T p+2p^TMr+r^TMr+ell^Tr."
        ),
    }


def rational_certificate() -> dict:
    lower_constant = Fraction(2603, 7200)
    if lower_constant / 4096 >= Fraction(1, 10_000):
        raise RuntimeError("lower endpoint retention reserve failed")

    upper_constant = Fraction(15_653, 1800)
    if upper_constant / 4096 >= Fraction(1, 100):
        raise RuntimeError("upper endpoint retention reserve failed")

    return {
        "lower_deviation_constant": str(lower_constant),
        "lower_retention": (
            "|Re L_0-1/2|<2603/(7200alpha)<1/10000, hence "
            "|L_0|>4999/10000."
        ),
        "upper_deviation_constant": str(upper_constant),
        "upper_retention": (
            "For every 0<=j<=5, |Re U_j-A_j(B)/2|<"
            "15653 A_j(B)/(1800alpha)<A_j(B)/100, hence "
            "|U_j|>49A_j(B)/100."
        ),
        "domain_reserve": (
            "The displayed floors use only alpha>4096, log(B)>1/10, "
            "|g'|<=51/100, pi>3, Z_2<=6u/alpha, and "
            "Z_3<=10u^2/alpha^2; the physical chart has alpha>2^49."
        ),
    }


def exact_certificate() -> dict:
    return {
        "domain": (
            "Retain q=1, j=0,...,5, the fixed roster T={m,...,n}, "
            "kappa=1/(2pi i), e_B=e(alpha log B), and all notation of "
            "Sections 11.167-11.170."
        ),
        "aggregates": (
            "Set F_j=sum_(r in T)U_(j,r), J_j=kappa sum_(r in T)"
            "e(phi_r(alpha/r)) integral c'_(j,r)e(-y^2/2), and let T_j "
            "be the twice-integrated outside-roster integral tail."
        ),
        "lower_package": (
            "L_j=delta_(j0)/2-calB_j(1)+kappa sum_(r in T)c_(j,r)(y_1)."
        ),
        "upper_package": (
            "U_j=A_j(B)/2+e(-alpha log B)calB_j(B)-"
            "kappa sum_(r in T)c_(j,r)(y_B)."
        ),
        "six_value_identity": (
            "H_j=F_j+L_j+e(alpha log B)U_j+J_j+T_j, 0<=j<=5."
        ),
        "lower_expansion": (
            "With d_j=A_j'(1)=-sigma delta_(j0)+delta_(j1), "
            "L_j=delta_(j0)/2-kappa delta_(j0)S_1(alpha)+"
            "kappa^2[d_jS_2(alpha)+alpha delta_(j0)S_3(alpha)]"
            "+kappa sum_T c_(j,r)(y_1)."
        ),
        "upper_expansion": (
            "For x_B=alpha/B, U_j=A_j(B)/2+kappa A_j(B)S_1(x_B)-"
            "kappa^2[A_j'(B)S_2(x_B)+(alpha/B^2)A_j(B)S_3(x_B)]"
            "-kappa sum_T c_(j,r)(y_B)."
        ),
        "amplitude_derivative": (
            "For lambda=log B, A_j'(B)=B^(-1)exp(t lambda^2/4-sigma lambda)"
            "{j lambda^(j-1)+(t lambda/2-sigma)lambda^j}, with the first "
            "term absent for j=0."
        ),
        "real_parts": (
            "All S_k and endpoint c values are real. Thus only the physical "
            "half endpoint and kappa^2 terms contribute to Re L_j and Re U_j."
        ),
        "retained_carrier": (
            "Define P_j^ec=F_j+L_j+e(alpha log B)U_j and "
            "rho_j^osc=J_j+T_j. Then H_j=P_j^ec+rho_j^osc exactly."
        ),
    }


def roster_certificate() -> dict:
    return {
        "endpoint_geometry": (
            "At a physical endpoint mu, q_r(mu)=alpha/mu-r and the exact "
            "Morse identity gives b_(j,r)(y_mu)/y_mu=-A_j(mu)/q_r(mu). "
            "Hence c_(j,r)(y_mu)=-A_j(mu)/q_r(mu)-b_(j,r)(0)/y_mu."
        ),
        "cotangent_audit": (
            "S_1(x)+sum_(r in T)1/(x-r)=pi cot(pi x). Substitution gives "
            "an equivalent full-cotangent form, but its cotangent pole and "
            "the b(0)/y pole cancel only after recombination. The stable "
            "digamma/Hurwitz package must be used at saddle crossings."
        ),
        "geometric_factors": (
            "1/2-kappa pi cot(pi x)=1/[1-e(x)] and "
            "1/2+kappa pi cot(pi x)=1/[1-e(-x)]. These are audit identities, "
            "not separately bounded endpoint terms."
        ),
        "one_mode_transfer": (
            "For any uniformly nonstationary mode moved across the roster "
            "boundary, its exact Morse representation U+Q equals its exact "
            "twice-integrated representation: I=U+kappa c(y_1)-kappa e_B "
            "c(y_B)+kappa e(phi_0)integral c'e(-y^2/2)="
            "kappa[e(phi)A/q]_1^B-kappa^2[e(phi)(D_rA)/q]_1^B+"
            "kappa^2 integral e(phi)C_rA."
        ),
        "invariance": (
            "Therefore the complete decomposition of H_j is roster-invariant. "
            "The separate packages L_j and U_j are bookkeeping-dependent and "
            "must not be interpreted without F_j, J_j, and T_j."
        ),
    }


def observation_certificate() -> dict:
    return {
        "normalized_split": (
            "After the fixed source normalization and common c_xi phase, let "
            "y=p+r be the twelve-real moment vector induced respectively by "
            "P^ec and rho^osc."
        ),
        "flow_expansion": (
            "Using Psi'=y^T M_xi y+ell_xi^T y and "
            "ell_xi=U_1t_T, one has exactly Psi'=p^T M_xi p+ell_xi^T p+"
            "2p^T M_xi r+r^T M_xi r+ell_xi^T r."
        ),
        "terminal_composition": (
            "The fixed terminal coordinates enter through t_T and hence "
            "through ell_xi in this identity. They do not belong inside the "
            "scalar Poisson endpoint package and no terminal cancellation is "
            "assumed before the observation map."
        ),
        "known_tail": (
            "The outside tail satisfies |T_j|<2K_jh^2 with "
            "K=(6,14,55,336,2738,27936). The grouped interior J_j remains "
            "unbounded at h^2 scale."
        ),
        "next_target": (
            "Insert the explicit carrier p into the eight observations and "
            "retain its endpoint, Hermitian, transpose, and terminal terms. "
            "Then prove a coefficient-aware bound for the five terms involving "
            "r, starting with grouped control of J_0,...,J_5; do not demand an "
            "h^2 bound for the retained endpoint carrier itself."
        ),
    }


def proof_boundary() -> str:
    return (
        "This proves six exact stable endpoint-composed identities, exact "
        "one-mode roster transfer, quantitative retention of the lower "
        "order-zero and all six factored upper endpoint packages, and the "
        "exact carrier/remainder expansion through the eight-observation flow "
        "including terminal coordinates. It does not bound the grouped "
        "oscillatory interior, prove an endpoint-composed h^2 remainder, "
        "evaluate the retained carrier's signed current, prove a signed flow "
        "or Phi_B bound, exclude contact, establish a retained aggregate or "
        "Xi theorem, prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level "
        "conclusion."
    )


def build_rows(
    exact: dict,
    roster: dict,
    rational: dict,
    observation: dict,
    boundary: str,
) -> list[CompositionRow]:
    return [
        CompositionRow("mfcr_01_domain", "physical domain", "proved", "The composition uses the certified q=1 fixed roster.", exact["domain"], "No sign follows."),
        CompositionRow("mfcr_02_aggregates", "aggregate definitions", "proved", "The roster main, interior residual, and outside tail are separated exactly.", exact["aggregates"], "Definitions only."),
        CompositionRow("mfcr_03_lower", "lower package", "proved", "All lower coherent terms form one stable package.", exact["lower_package"], "No smallness is asserted."),
        CompositionRow("mfcr_04_upper", "upper package", "proved", "All upper coherent terms form one factored stable package.", exact["upper_package"], "The common endpoint phase remains outside."),
        CompositionRow("mfcr_05_identity", "six-value identity", "proved", "Each of the six moments has one exact endpoint-composed form.", exact["six_value_identity"], "No absolute value is taken."),
        CompositionRow("mfcr_06_lower_expand", "lower expansion", "proved", "The lower package is explicit in stable endpoint sums.", exact["lower_expansion"], "Uses the fixed roster."),
        CompositionRow("mfcr_07_upper_expand", "upper expansion", "proved", "The upper package is explicit in stable endpoint sums.", exact["upper_expansion"], "Uses the fixed roster."),
        CompositionRow("mfcr_08_derivative", "amplitude derivative", "proved", "The upper derivative is explicit for all six moments.", exact["amplitude_derivative"], "No asymptotic replacement."),
        CompositionRow("mfcr_09_real", "real-part structure", "proved", "Coherent first-order traces cannot alter anchored real parts.", exact["real_parts"], "Cross-endpoint cancellation remains possible."),
        CompositionRow("mfcr_10_lower_floor", "lower retention", "proved", "The lower order-zero package retains the physical half endpoint.", rational["lower_retention"], rational["domain_reserve"]),
        CompositionRow("mfcr_11_upper_floor", "upper retention", "proved", "All factored upper packages retain their physical half endpoints.", rational["upper_retention"], rational["domain_reserve"]),
        CompositionRow("mfcr_12_geometry", "endpoint geometry", "proved", "The transformed endpoint quotient has an exact rational form.", roster["endpoint_geometry"], "Away from a crossing before removable continuation."),
        CompositionRow("mfcr_13_cotangent", "cotangent audit", "proved", "The stable and full-cotangent forms agree.", roster["cotangent_audit"], "The cotangent form is numerically unstable at crossings."),
        CompositionRow("mfcr_14_geometric", "geometric factors", "proved", "The half endpoint and cotangent combine into exact geometric factors.", roster["geometric_factors"], "Audit identity only."),
        CompositionRow("mfcr_15_transfer", "one-mode transfer", "proved", "A mode has two exactly equal nonstationary representations.", roster["one_mode_transfer"], "Uniformly nonstationary mode only."),
        CompositionRow("mfcr_16_invariance", "roster invariance", "proved", "Changing an admissible roster does not change H_j.", roster["invariance"], "Individual packages are not canonical."),
        CompositionRow("mfcr_17_carrier", "retained carrier", "proved", "The large endpoint terms are retained in the main carrier.", exact["retained_carrier"], "The oscillatory remainder remains open."),
        CompositionRow("mfcr_18_normalize", "normalization", "proved", "The carrier split passes through the fixed source and common phase.", observation["normalized_split"], "No coefficient norm is applied."),
        CompositionRow("mfcr_19_flow", "observation expansion", "proved", "The endpoint carrier and residual enter the flow polynomial exactly.", observation["flow_expansion"], "No sign follows from the expansion."),
        CompositionRow("mfcr_20_terminal", "terminal placement", "proved", "Terminal data compose at the observation map.", observation["terminal_composition"], "No scalar terminal cancellation is assumed."),
        CompositionRow("mfcr_21_tail", "known tail", "proved", "The only already bounded residual piece is the outside integral tail.", observation["known_tail"], "Large constants still need coefficient-aware propagation."),
        CompositionRow("mfcr_22_interior", "interior theorem", "open", "The grouped c-prime interior is the live analytic remainder.", observation["next_target"], "No h^2 interior bound is claimed."),
        CompositionRow("mfcr_23_pi", "pi provenance", "proved", "Every pi is inherited from Fourier normalization.", "kappa=1/(2pi i) comes from e(x)=exp(2pi i x); pi cot(pi x) is the symmetric partial-fraction sum.", "No fitted geometry."),
        CompositionRow("mfcr_24_boundary", "proof boundary", "proved", "The composition is not the missing signed estimate.", boundary, "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact_certificate"]
    roster = payload["roster_certificate"]
    rational = payload["rational_certificate"]
    observation = payload["observation_certificate"]
    return f"""# Six-Moment Morse-Fresnel Endpoint-Composition Retention Gate

Date: 2026-08-02

Status: six exact endpoint-composed value identities, exact roster transfer,
and retained endpoint floors; `0 oscillatory interior bounds`, `0 signed flow
bounds`, and this is not a proof of RH.

## Exact Composition

{exact['domain']}

{exact['aggregates']}

Define

```text
{exact['lower_package']}

{exact['upper_package']}
```

Direct substitution into finite starred Poisson summation gives, without any
absolute value,

```text
{exact['six_value_identity']}
```

The stable endpoint formulas are

```text
{exact['lower_expansion']}

{exact['upper_expansion']}
```

{exact['amplitude_derivative']}

## Endpoint Retention

{exact['real_parts']}

{rational['lower_retention']}

{rational['upper_retention']}

{rational['domain_reserve']}

Thus endpoint composition does not make the hard endpoints small. They must
be retained in the comparison carrier.

## Roster Audit

{roster['endpoint_geometry']}

{roster['cotangent_audit']}

{roster['geometric_factors']}

{roster['one_mode_transfer']}

{roster['invariance']}

## Carrier And Remainder

{exact['retained_carrier']}

{observation['known_tail']}

## Eight-Observation Composition

{observation['normalized_split']}

{observation['flow_expansion']}

{observation['terminal_composition']}

## Next Target

{observation['next_target']}

## Pi Provenance

`kappa=1/(2pi i)` is forced by `e(x)=exp(2pi i x)`, while
`pi cot(pi x)` is the symmetric partial-fraction sum over integral Poisson
modes. No circle or fitted geometric constant is introduced.

## Proof Boundary

{payload['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    sources = load_sources()
    symbolic = symbolic_audit()
    rational = rational_certificate()
    exact = exact_certificate()
    roster = roster_certificate()
    observation = observation_certificate()
    boundary = proof_boundary()
    rows = build_rows(exact, roster, rational, observation, boundary)
    return {
        "kind": STEM,
        "date": "2026-08-02",
        "status": (
            "six endpoint-composed identities and endpoint retention proved; "
            "grouped oscillatory interior and signed flow remain open"
        ),
        "source_audit": source_audit(sources),
        "symbolic_audit": symbolic,
        "rational_certificate": rational,
        "exact_certificate": exact,
        "roster_certificate": roster,
        "observation_certificate": observation,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "endpoint_composed_value_identities": 6,
            "stable_endpoint_packages": 2,
            "roster_transfer_representations": 2,
            "retained_endpoint_floors": 7,
            "observation_expansions": 1,
            "enumerated_physical_modes": 0,
            "oscillatory_interior_bounds": 0,
            "endpoint_composed_h2_bounds": 0,
            "evaluated_signed_carriers": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
        "proof_boundary": boundary,
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
        "built Morse-Fresnel endpoint-composition retention gate: "
        f"{counts['rows']} rows, "
        f"{counts['endpoint_composed_value_identities']} composed identities, "
        f"{counts['roster_transfer_representations']} roster representations, "
        f"{counts['retained_endpoint_floors']} retained endpoint floors, "
        f"{counts['observation_expansions']} observation expansion, "
        f"{counts['oscillatory_interior_bounds']} oscillatory interior bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
