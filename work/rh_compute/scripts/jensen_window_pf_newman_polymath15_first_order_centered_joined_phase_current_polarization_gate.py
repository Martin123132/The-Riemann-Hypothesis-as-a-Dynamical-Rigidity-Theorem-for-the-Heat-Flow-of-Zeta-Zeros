#!/usr/bin/env python3
"""Build the exact joined phase-current polarization and contact guard."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "joined_phase_current_polarization_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "ray_bottom": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "ray_bottom_logarithmic_flow_reduction.json"
    ),
    "phase_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_logarithmic_phase_flow_gate.json"
    ),
    "large_prime": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_length_propagation_gate.json"
    ),
    "ternary": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_ternary_length_gate.json"
    ),
    "dyadic": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_dyadic_length_gate.json"
    ),
    "short": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_short_family_counter_gate.json"
    ),
}

ComplexQ = tuple[Fraction, Fraction]


@dataclass(frozen=True)
class GateRow:
    id: str
    claim: str
    readiness: str
    exact_statement: str
    implication: str
    boundary: str


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def complex_text(value: ComplexQ) -> dict[str, str]:
    return {
        "real": fraction_text(value[0]),
        "imag": fraction_text(value[1]),
    }


def add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def scale(value: ComplexQ, scalar: Fraction) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def conjugate(value: ComplexQ) -> ComplexQ:
    return value[0], -value[1]


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def current(derivative: ComplexQ, value: ComplexQ) -> Fraction:
    return multiply(derivative, conjugate(value))[1]


def norm_squared(value: ComplexQ) -> Fraction:
    return value[0] ** 2 + value[1] ** 2


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def source_audit() -> dict:
    payloads = {}
    hashes = {}
    for key, path in SOURCES.items():
        if not path.is_file():
            raise RuntimeError(f"missing source {key}: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
        hashes[key] = file_hash(path)

    if payloads["dyadic"].get("summary", {}).get(
        "complete_dyadic_minimum_length"
    ) != 8:
        raise RuntimeError("dyadic all-length source drifted")
    if payloads["ternary"].get("summary", {}).get(
        "actual_all_length_ternary_families"
    ) != 1:
        raise RuntimeError("ternary all-length source drifted")
    if payloads["large_prime"].get("summary", {}).get(
        "analytic_all_length_large_prime_families"
    ) != 1:
        raise RuntimeError("large-prime source drifted")
    if payloads["short"].get("summary", {}).get(
        "exact_interior_negative_witnesses"
    ) != 8:
        raise RuntimeError("short-family source drifted")
    if payloads["ray_bottom"].get("summary", {}).get(
        "proved_abel_gaps"
    ) != 0:
        raise RuntimeError("ray-bottom proof boundary drifted")
    if payloads["phase_flow"].get("summary", {}).get(
        "replacement_starlikeness_targets"
    ) != 1:
        raise RuntimeError("prime-power phase-flow source drifted")
    return {"source_sha256": hashes}


def polarization_witness() -> dict:
    values = (
        (Fraction(1, 2), Fraction(-1, 3)),
        (Fraction(-2, 5), Fraction(3, 7)),
        (Fraction(4, 9), Fraction(1, 6)),
    )
    derivatives = (
        (Fraction(2, 3), Fraction(1, 5)),
        (Fraction(-1, 4), Fraction(5, 8)),
        (Fraction(7, 10), Fraction(-2, 11)),
    )
    total_value = (Fraction(0), Fraction(0))
    total_derivative = (Fraction(0), Fraction(0))
    for value, derivative in zip(values, derivatives, strict=True):
        total_value = add(total_value, value)
        total_derivative = add(total_derivative, derivative)

    diagonal = [current(derivative, value) for value, derivative in zip(
        values, derivatives, strict=True
    )]
    cross = []
    for left in range(len(values)):
        for right in range(left + 1, len(values)):
            term = current(derivatives[left], values[right])
            term += current(derivatives[right], values[left])
            cross.append(
                {
                    "left": left,
                    "right": right,
                    "current": fraction_text(term),
                }
            )
    joined = current(total_derivative, total_value)
    reconstructed = sum(diagonal, Fraction(0)) + sum(
        (Fraction(row["current"]) for row in cross),
        Fraction(0),
    )
    if joined != reconstructed:
        raise RuntimeError("joined-current polarization identity failed")
    return {
        "values": [complex_text(value) for value in values],
        "derivatives": [
            complex_text(derivative) for derivative in derivatives
        ],
        "total_value": complex_text(total_value),
        "total_derivative": complex_text(total_derivative),
        "joined_current": fraction_text(joined),
        "diagonal_currents": [
            fraction_text(value) for value in diagonal
        ],
        "cross_currents": cross,
        "reconstructed_current": fraction_text(reconstructed),
        "identity_exact": True,
    }


def chain_insertion_witness() -> dict:
    # Q=1+(1/2)i, H=(1/2)i, rho=E=0.
    q_value = (Fraction(1), Fraction(1, 2))
    h_value = (Fraction(0), Fraction(1, 2))
    q_norm = norm_squared(q_value)
    h_q_real = multiply(h_value, conjugate(q_value))[0]
    j_ray = q_norm + h_q_real
    tau = Fraction(2)
    projector_rate = Fraction(-1, 7)
    external_rate = Fraction(3, 5)
    internal_derivative = (Fraction(-1), Fraction(0))
    total_derivative = add(
        internal_derivative,
        (
            -(projector_rate + external_rate) * q_value[1],
            (projector_rate + external_rate) * q_value[0],
        ),
    )
    direct = current(total_derivative, q_value)
    inserted = tau * j_ray
    inserted += (
        projector_rate + external_rate - tau
    ) * q_norm
    if direct != inserted:
        raise RuntimeError("chain self-current insertion identity failed")
    return {
        "Q": complex_text(q_value),
        "H": complex_text(h_value),
        "Q_norm_squared": fraction_text(q_norm),
        "J_ray": fraction_text(j_ray),
        "tau": fraction_text(tau),
        "projector_rate": fraction_text(projector_rate),
        "external_rate": fraction_text(external_rate),
        "direct_self_current": fraction_text(direct),
        "inserted_self_current": fraction_text(inserted),
        "identity_exact": True,
    }


def two_frequency_guard(amplitude: Fraction) -> dict:
    # At theta=pi/2, exp(i theta)=i and exp(3i theta)=-i.
    first = (Fraction(0), Fraction(1))
    second = (Fraction(0), -amplitude)
    first_derivative = (Fraction(-1), Fraction(0))
    second_derivative = (3 * amplitude, Fraction(0))
    joined = add(first, second)
    joined_derivative = add(first_derivative, second_derivative)
    first_current = current(first_derivative, first)
    second_current = current(second_derivative, second)
    cross_current = current(first_derivative, second)
    cross_current += current(second_derivative, first)
    joined_current = current(joined_derivative, joined)
    formula = (1 - amplitude) * (1 - 3 * amplitude)
    linear_coefficient = 1 - 3 * amplitude
    cubic_coefficient = 4 * amplitude
    if joined_current != formula:
        raise RuntimeError("two-frequency joined-current formula failed")
    return {
        "amplitude": fraction_text(amplitude),
        "theta": "pi/2",
        "first_value": complex_text(first),
        "second_value": complex_text(second),
        "joined_value": complex_text(joined),
        "joined_derivative": complex_text(joined_derivative),
        "first_self_current": fraction_text(first_current),
        "second_self_current": fraction_text(second_current),
        "cross_current": fraction_text(cross_current),
        "joined_current": fraction_text(joined_current),
        "factorized_current": "(1-a)(1-3a)",
        "real_projection_chebyshev_form": (
            "(1-3a)x+4a x^3, x=cos(theta)"
        ),
        "real_projection_linear_coefficient": fraction_text(
            linear_coefficient
        ),
        "real_projection_cubic_coefficient": fraction_text(
            cubic_coefficient
        ),
        "identity_exact": True,
    }


def gate_rows() -> list[GateRow]:
    return [
        GateRow(
            "jpcp_01_polarization",
            "joined phase-current polarization",
            "proved",
            "K(W)=sum_a K(w_a)+sum_(a<b) K_(a,b).",
            "Every block theorem enters a joined current with explicit cross terms.",
            "The identity supplies no sign for the cross terms.",
        ),
        GateRow(
            "jpcp_02_contact",
            "contact-current identity",
            "proved",
            "K=X D_jY-Y D_jX=X D_jY-xY mathcal_C_N-Y E_j.",
            "At X=0, K=-Y(x mathcal_C_N+E_j).",
            "When Y=0 the current vanishes and cannot replace the Abel scalar.",
        ),
        GateRow(
            "jpcp_03_chain",
            "prime-power self-current insertion",
            "proved",
            "K_m=|B_m|^2{tau_m J_ray,m+(Omega_eta+e_m-tau_m)|Q_m|^2}.",
            "The certified J_ray margin is one diagonal summand.",
            "External rate and joined cross-current terms remain.",
        ),
        GateRow(
            "jpcp_04_long",
            "long-chain diagonal inputs",
            "proved",
            "Actual J_ray>0 is certified for p=2,M>=8; p=3,M>=4; p>=5,M>=2.",
            "All complete long families now have positive internal diagonal input.",
            "No sum or endpoint sign follows from those inputs.",
        ),
        GateRow(
            "jpcp_05_cubic",
            "positive-self-current cubic-contact guard",
            "proved",
            "cos(theta)+(1/3)cos(3theta)=(4/3)cos(theta)^3.",
            "At theta=pi/2 the joined real projection and its first derivative vanish.",
            "This is a generic exact guard, not an Xi state.",
        ),
        GateRow(
            "jpcp_06_negative",
            "negative joined-current guard",
            "proved",
            "For W_a=e^(i theta)+a e^(3i theta), K(pi/2)=(1-a)(1-3a).",
            "Both self currents are positive while K<0 for 1/3<a<1.",
            "Individual starlikeness is not closed under addition.",
        ),
        GateRow(
            "jpcp_07_zero_fibre",
            "zero-fibre nonpromotion",
            "proved",
            "W=0 implies K=Im((D_jW)conj(W))=0.",
            "A direct division-free Abel estimate is still required at W_0=0.",
            "No quotient by Y or W is permitted there.",
        ),
        GateRow(
            "jpcp_08_route",
            "joined cross-current handoff",
            "proved",
            "Retain endpoint, external rates, short blocks, singletons, and every pair cross current before estimating K or mathcal_C_N.",
            "The next phase route is an Xi-specific joined cross-current theorem.",
            "The alternative remains the direct endpoint-complete Abel gap.",
        ),
        GateRow(
            "jpcp_09_pi",
            "pi provenance",
            "proved",
            "The guard uses the ordinary periods of exp(i theta) and cosine.",
            "No new pi is introduced into the Xi saddle formulas.",
            "The guard frequency 3 is not an Xi parameter claim.",
        ),
        GateRow(
            "jpcp_10_boundary",
            "proof boundary",
            "open",
            "No joined Xi cross-current or Abel-gap estimate is proved.",
            "No successor winding or contact-exclusion theorem is promoted.",
            "Lambda<=0, PF-infinity, RH, and a prize-level conclusion remain open.",
        ),
    ]


def render_note(payload: dict) -> str:
    cubic = payload["guards"]["cubic_contact"]
    negative = payload["guards"]["negative_current"]
    boundary = payload["proof_boundary"]
    return f"""# Joined Phase-Current Polarization Gate

Date: 2026-07-30

Status: exact joined-current identity and nonpromotion gate. This is not
a proof of a joined Xi Abel-gap theorem, `Lambda<=0`, PF-infinity, or RH.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Exact Polarization

For a finite endpoint-complete decomposition `W=sum_a w_a`, define

```text
K(W)=Im((D_j W) conjugate(W)).
```

Direct expansion gives

```text
K(W)
 =sum_a Im((D_j w_a)conjugate(w_a))
  +sum_(a<b) Im((D_j w_a)conjugate(w_b)
                +(D_j w_b)conjugate(w_a)).
```

The checker reconstructs this identity over exact rational complex
arithmetic. No quotient or nonvanishing assumption is used.

## Contact Identity

Write `W=mathsf_X+i mathsf_Y` and retain the certified orientation
defect

```text
E_j=D_j mathsf_X-x mathcal_C_N.
```

Then

```text
K=mathsf_X D_j mathsf_Y
  -x mathsf_Y mathcal_C_N
  -mathsf_Y E_j.
```

At `mathsf_X=0`,

```text
K=-mathsf_Y(x mathcal_C_N+E_j).
```

If `mathsf_Y=0`, then `W=0` and `K=0`; phase current cannot replace
the division-free Abel scalar on that fibre.

## Prime-Power Insertion

For one physical chain `W_m=eta B_m Q_m`, let

```text
tau_m=-x log(p) b>0,
e_m=x(beta_m'-b log(m)),
eta_L=i Omega_eta eta.
```

The exact fixed-ray formula gives

```text
K_m
 =|B_m|^2[
    tau_m J_ray,m
    +(Omega_eta+e_m-tau_m)|Q_m|^2].
```

Thus the proved positive currents for `p=2,M>=8`,
`p=3,M>=4`, and `p>=5,M>=2` enter the joined current only through
the displayed diagonal term. External rates, the recurrent endpoint,
short chains, singletons, and every cross current remain.

## Exact Guards

For

```text
W_a(theta)=exp(i theta)+a exp(3i theta),
```

the joined current at `theta=pi/2` is exactly

```text
K=(1-a)(1-3a).
```

Both self currents are positive: `1` and `3a^2`. Nevertheless:

```text
a=1/3:
  Re(W)=0, Re(W')=0,
  cos(theta)+(1/3)cos(3theta)=(4/3)cos(theta)^3,
  K={cubic["joined_current"]};

a=1/2:
  K={negative["joined_current"]}<0.
```

The cubic row is a real-projection contact produced by exact cross-current
cancellation. The negative row shows that individual starlikeness is not
closed under addition. These are generic route guards, not Xi
counterexamples or phase-attainment claims.

## Route Decision

The long-chain theorems remain valid and useful diagonal input, but they
do not by themselves lower-bound the endpoint-complete Abel scalar or
the joined phase current. A viable continuation must either:

1. prove an Xi-specific bound for the endpoint, external-rate, short,
   singleton, and pair cross-current aggregate; or
2. prove the division-free Abel-scalar gap directly, retaining `W_0=0`.

## Pi Provenance

The `pi/2` in the guard is the ordinary quarter-period of the complex
exponential. Any `pi` in the Xi saddle scale remains inherited from
completed-zeta normalization; no polygon or prime block defines it.

## Boundary

{boundary}
"""


def main() -> int:
    polarization = polarization_witness()
    insertion = chain_insertion_witness()
    cubic = two_frequency_guard(Fraction(1, 3))
    negative = two_frequency_guard(Fraction(1, 2))
    rows = gate_rows()
    proof_boundary = (
        "This gate proves the exact joined-current polarization, contact "
        "identity, prime-power diagonal insertion, and generic cubic-contact "
        "and negative-current guards. It proves no Xi-specific cross-current "
        "bound, endpoint or singleton absorption, Abel-scalar gap, horizontal "
        "successor winding cap, contact exclusion, Q209, cofinal descendant "
        "theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
    )
    payload = {
        "kind": "joined_phase_current_polarization_gate",
        "date": "2026-07-30",
        "status": (
            "exact joined-current polarization and cubic-contact "
            "nonpromotion gate"
        ),
        "proof_boundary": proof_boundary,
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "exact": {
            "polarization": (
                "K(sum_a w_a)=sum_a K(w_a)+sum_(a<b) "
                "Im((D_jw_a)conj(w_b)+(D_jw_b)conj(w_a))"
            ),
            "contact": (
                "K=mathsf_X D_jmathsf_Y-x mathsf_Y mathcal_C_N"
                "-mathsf_Y E_j"
            ),
            "chain_insertion": (
                "K_m=|B_m|^2{tau_m J_ray,m+"
                "(Omega_eta+e_m-tau_m)|Q_m|^2}"
            ),
            "two_frequency": (
                "K_[exp(i theta)+a exp(3i theta)](pi/2)"
                "=(1-a)(1-3a)"
            ),
        },
        "witnesses": {
            "polarization": polarization,
            "chain_insertion": insertion,
        },
        "guards": {
            "cubic_contact": cubic,
            "negative_current": negative,
        },
        "summary": {
            "rows": len(rows),
            "exact_polarization_witnesses": 1,
            "exact_chain_insertion_witnesses": 1,
            "positive_self_current_cubic_contacts": 1,
            "negative_joined_current_witnesses": 1,
            "long_family_diagonal_inputs": 3,
            "proved_joined_cross_current_bounds": 0,
            "proved_joined_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
    }
    payload["content_sha256"] = hashlib.sha256(
        canonical_json(payload).encode("utf-8")
    ).hexdigest()
    atomic_write(
        RESULT,
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
    )
    atomic_write(NOTE, render_note(payload))
    print(
        "built joined phase-current polarization gate: "
        "1 exact polarization, 1 exact chain insertion, "
        "1 cubic contact, 1 negative joined-current guard"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
