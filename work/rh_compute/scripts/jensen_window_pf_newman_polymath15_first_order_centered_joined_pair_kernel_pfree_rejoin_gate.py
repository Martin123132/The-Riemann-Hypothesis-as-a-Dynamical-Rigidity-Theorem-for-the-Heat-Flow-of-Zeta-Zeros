#!/usr/bin/env python3
"""Build the physical pair kernel and p-free quadratic-rejoin audit."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "joined_pair_kernel_pfree_rejoin_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "polarization": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_phase_current_polarization_gate.json"
    ),
    "ray_bottom": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "ray_bottom_logarithmic_flow_reduction.json"
    ),
    "joined_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_dyadic_odd_prefix_first_jet_reduction.json"
    ),
    "complete_chain": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "complete_prime_power_chain_occupation_gate.json"
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


def subtract(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] - right[0], left[1] - right[1]


def scale(value: ComplexQ, scalar: Fraction) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def conjugate(value: ComplexQ) -> ComplexQ:
    return value[0], -value[1]


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def norm_squared(value: ComplexQ) -> Fraction:
    return value[0] ** 2 + value[1] ** 2


def current(derivative: ComplexQ, value: ComplexQ) -> Fraction:
    return multiply(derivative, conjugate(value))[1]


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
    if payloads["polarization"].get("summary", {}).get(
        "exact_polarization_witnesses"
    ) != 1:
        raise RuntimeError("polarization source drifted")
    if payloads["ray_bottom"].get("summary", {}).get(
        "exact_logarithmic_flow_identities"
    ) != 8:
        raise RuntimeError("ray-bottom source drifted")
    if payloads["complete_chain"].get("summary", {}).get(
        "p_free_telescoping_identities"
    ) != 1:
        raise RuntimeError("complete-chain p-free source drifted")
    return {"source_sha256": hashes}


def physical_kernel_witness() -> dict:
    x_value = Fraction(5, 3)
    c_value = Fraction(1, 7)
    b_value = Fraction(-2, 5)
    omega = Fraction(3, 11)
    logarithms = (Fraction(0), Fraction(2, 3), Fraction(5, 4))
    carriers = (
        (Fraction(1, 2), Fraction(1, 3)),
        (Fraction(-2, 5), Fraction(3, 7)),
        (Fraction(4, 9), Fraction(-1, 6)),
    )
    corrections = (
        (Fraction(0), Fraction(0)),
        (Fraction(1, 13), Fraction(-2, 17)),
        (Fraction(-3, 19), Fraction(1, 23)),
    )

    rates: list[ComplexQ] = []
    derivatives: list[ComplexQ] = []
    h_0 = (Fraction(0), Fraction(0))
    h_1 = (Fraction(0), Fraction(0))
    d_0 = (Fraction(0), Fraction(0))
    for ell, carrier, correction in zip(
        logarithms, carriers, corrections, strict=True
    ):
        rate = (
            x_value * (-c_value * ell + correction[0]),
            omega + x_value * (-b_value * ell + correction[1]),
        )
        rates.append(rate)
        derivatives.append(multiply(rate, carrier))
        h_0 = add(h_0, carrier)
        h_1 = add(h_1, scale(carrier, ell))
        d_0 = add(d_0, multiply(correction, carrier))

    total_derivative = (Fraction(0), Fraction(0))
    for derivative in derivatives:
        total_derivative = add(total_derivative, derivative)
    direct_current = current(total_derivative, h_0)

    diagonal_rows = []
    diagonal_sum = Fraction(0)
    for index, (rate, carrier) in enumerate(
        zip(rates, carriers, strict=True)
    ):
        value = rate[1] * norm_squared(carrier)
        diagonal_sum += value
        diagonal_rows.append(
            {
                "index": index,
                "current": fraction_text(value),
            }
        )

    pair_rows = []
    pair_sum = Fraction(0)
    for left in range(len(carriers)):
        for right in range(left + 1, len(carriers)):
            correlation = multiply(
                carriers[left], conjugate(carriers[right])
            )
            formula = (
                (rates[left][0] - rates[right][0]) * correlation[1]
                + (rates[left][1] + rates[right][1]) * correlation[0]
            )
            direct = current(derivatives[left], carriers[right])
            direct += current(derivatives[right], carriers[left])
            if direct != formula:
                raise RuntimeError("physical pair-kernel identity failed")
            pair_sum += formula
            pair_rows.append(
                {
                    "left": left,
                    "right": right,
                    "correlation": complex_text(correlation),
                    "direct_current": fraction_text(direct),
                    "kernel_current": fraction_text(formula),
                }
            )

    s_prime = (c_value, b_value)
    leading = multiply(scale(s_prime, -1), h_1)
    bracket = add(leading, d_0)
    moment_current = omega * norm_squared(h_0)
    moment_current += x_value * multiply(
        bracket, conjugate(h_0)
    )[1]
    expanded_moment = omega * norm_squared(h_0)
    h_correlation = multiply(h_1, conjugate(h_0))
    expanded_moment -= x_value * c_value * h_correlation[1]
    expanded_moment -= x_value * b_value * h_correlation[0]
    expanded_moment += x_value * multiply(
        d_0, conjugate(h_0)
    )[1]
    if not (
        direct_current
        == diagonal_sum + pair_sum
        == moment_current
        == expanded_moment
    ):
        raise RuntimeError("physical moment collapse failed")

    return {
        "parameters": {
            "x": fraction_text(x_value),
            "c": fraction_text(c_value),
            "b": fraction_text(b_value),
            "Omega_eta": fraction_text(omega),
        },
        "logarithmic_nodes": [
            fraction_text(value) for value in logarithms
        ],
        "carriers": [complex_text(value) for value in carriers],
        "corrections": [complex_text(value) for value in corrections],
        "rates": [complex_text(value) for value in rates],
        "H_0": complex_text(h_0),
        "H_1": complex_text(h_1),
        "D_0": complex_text(d_0),
        "diagonal_rows": diagonal_rows,
        "pair_rows": pair_rows,
        "diagonal_sum": fraction_text(diagonal_sum),
        "pair_sum": fraction_text(pair_sum),
        "direct_current": fraction_text(direct_current),
        "moment_current": fraction_text(moment_current),
        "expanded_moment_current": fraction_text(expanded_moment),
        "identity_exact": True,
    }


def heat_shift_witness() -> dict:
    rows = []
    for index in range(9):
        left_t = Fraction(index * index, 4)
        left_t += Fraction(1, 4) + Fraction(index, 2)
        right_t = Fraction((index + 1) ** 2, 4)
        left_s = -(index + 1)
        right_s = -(index + 1)
        if left_t != right_t or left_s != right_s:
            raise RuntimeError("heat-shift exponent identity failed")
        rows.append(
            {
                "k": index,
                "left_t_h2_coefficient": fraction_text(left_t),
                "right_t_h2_coefficient": fraction_text(right_t),
                "left_s_h_coefficient": fraction_text(Fraction(left_s)),
                "right_s_h_coefficient": fraction_text(Fraction(right_s)),
            }
        )
    return {
        "identity": (
            "A_t(p^k;s)A_t(p;s-t*k*h/2)=A_t(p^(k+1);s)"
        ),
        "checked_integer_k": [0, 8],
        "coefficient_rows": rows,
        "identity_exact": True,
    }


def quadratic_rejoin_witness() -> dict:
    # L_0=T_0-T_1=i and L_1=T_1=-i/3 at one exact jet.
    t_0 = (Fraction(0), Fraction(2, 3))
    t_1 = (Fraction(0), Fraction(-1, 3))
    t_2 = (Fraction(0), Fraction(0))
    dt_0 = (Fraction(0), Fraction(0))
    dt_1 = (Fraction(1), Fraction(0))
    dt_2 = (Fraction(0), Fraction(0))
    layers = (subtract(t_0, t_1), subtract(t_1, t_2))
    derivatives = (
        subtract(dt_0, dt_1),
        subtract(dt_1, dt_2),
    )
    diagonal = sum(
        (
            current(derivative, layer)
            for derivative, layer in zip(
                derivatives, layers, strict=True
            )
        ),
        Fraction(0),
    )
    cross = current(derivatives[0], layers[1])
    cross += current(derivatives[1], layers[0])
    joined_value = add(layers[0], layers[1])
    joined_derivative = add(derivatives[0], derivatives[1])
    joined = current(joined_derivative, joined_value)
    if joined_value != t_0 or joined_derivative != dt_0:
        raise RuntimeError("linear p-free telescope failed")
    if joined != diagonal + cross:
        raise RuntimeError("quadratic p-free rejoin failed")
    return {
        "T_values": [
            complex_text(value) for value in (t_0, t_1, t_2)
        ],
        "T_derivatives": [
            complex_text(value) for value in (dt_0, dt_1, dt_2)
        ],
        "layers": [complex_text(value) for value in layers],
        "layer_derivatives": [
            complex_text(value) for value in derivatives
        ],
        "linear_joined_value": complex_text(joined_value),
        "linear_joined_derivative": complex_text(joined_derivative),
        "diagonal_current_sum": fraction_text(diagonal),
        "pair_cross_current": fraction_text(cross),
        "joined_current": fraction_text(joined),
        "linear_telescope_exact": True,
        "quadratic_rejoin_exact": True,
    }


def gate_rows() -> list[GateRow]:
    return [
        GateRow(
            "jpkr_01_rate",
            "physical carrier rate",
            "proved",
            "Gamma_n=x[-s_*'log(n)+delta_n]+i Omega_eta.",
            "D_j(eta q_n)=Gamma_n eta q_n exactly.",
            "No aggregate sign is inferred.",
        ),
        GateRow(
            "jpkr_02_pair",
            "division-free pair kernel",
            "proved",
            "K_nm=(a_n-a_m)Im(q_n conjugate(q_m))+(omega_n+omega_m)Re(q_n conjugate(q_m)).",
            "Every physical pair cross current is explicit.",
            "Ordered rates multiply unrestricted oscillatory correlations.",
        ),
        GateRow(
            "jpkr_03_moment",
            "bulk moment collapse",
            "proved",
            "K_bulk=Omega_eta|H_0|^2-xc Im(H_1conj(H_0))-xb Re(H_1conj(H_0))+x Im(D_0conj(H_0)).",
            "All diagonal and pair rows rejoin to H_0,H_1,D_0.",
            "The recurrent endpoint is still separate.",
        ),
        GateRow(
            "jpkr_04_heat_shift",
            "heat-shifted p-free telescope",
            "proved",
            "L_k=T_k-T_(k+1), so sum_k L_k=T_0=S_(N,t)(s).",
            "The p-free decomposition remains exact under the heat shift.",
            "Actual d_n corrections retain the unique partition but not this simple factorization.",
        ),
        GateRow(
            "jpkr_05_quadratic",
            "quadratic current rejoin",
            "proved",
            "K(sum_k L_k)=sum_k K(L_k)+sum_(k<l)K_(k,l)=K(T_0).",
            "All p-free cross currents reconstruct the original full-prefix current.",
            "Linear telescoping is not a positivity transformation.",
        ),
        GateRow(
            "jpkr_06_guard",
            "strict cross-cancellation guard",
            "proved",
            "Two exact telescoping layers have diagonal current 4/3, cross current -4/3, and joined current 0.",
            "The full positive diagonal reserve can disappear after rejoining.",
            "This is a generic algebraic guard, not an Xi state.",
        ),
        GateRow(
            "jpkr_07_correction",
            "correction-current localization",
            "proved",
            "The complete coefficient correction is x Im(D_0 conjugate(H_0)).",
            "No pair correction is omitted or double-counted.",
            "Its crude absolute bound is not promoted as an Abel-gap estimate.",
        ),
        GateRow(
            "jpkr_08_route",
            "p-free route decision",
            "proved",
            "P-free reindexing is an exact coordinate change for the joined current, not an algebraic sign theorem.",
            "Any gain must be an Xi-specific oscillatory-correlation estimate.",
            "Otherwise return to the direct division-free Abel scalar.",
        ),
        GateRow(
            "jpkr_09_pi",
            "pi provenance",
            "proved",
            "No pi enters the carrier-pair or p-free telescope identities.",
            "Any Xi saddle pi remains inherited from completed zeta.",
            "The generic guard uses only ordinary complex phases.",
        ),
        GateRow(
            "jpkr_10_boundary",
            "proof boundary",
            "open",
            "No Xi pair-correlation sign or Abel gap is proved.",
            "No endpoint-complete winding cap or contact exclusion is promoted.",
            "Lambda<=0, PF-infinity, RH, and a prize-level conclusion remain open.",
        ),
    ]


def render_note(payload: dict) -> str:
    witness = payload["witnesses"]["quadratic_rejoin"]
    boundary = payload["proof_boundary"]
    return f"""# Joined Pair Kernel And P-Free Rejoin Gate

Date: 2026-07-30

Status: exact carrier-pair and p-free rejoin audit. This is not a proof
of an Xi cross-current sign, an Abel gap, `Lambda<=0`, PF-infinity, or RH.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Physical Pair Kernel

On one fixed ray cell put

```text
w_n=eta q_n,
Gamma_n=x[-s_*'log(n)+delta_n]+i Omega_eta
       =a_n+i omega_n.
```

Then `D_j w_n=Gamma_n w_n`. For `n<m`, exact rational
polarization gives the division-free pair current

```text
K_(n,m)
 =Im[(D_jw_n)conjugate(w_m)+(D_jw_m)conjugate(w_n)]

 =(a_n-a_m)Im[q_n conjugate(q_m)]
  +(omega_n+omega_m)Re[q_n conjugate(q_m)].
```

With `s_*'=c+ib` and `delta_n=xi_n+i chi_n`,

```text
a_n-a_m
 =x[c log(m/n)+xi_n-xi_m],

omega_n+omega_m
 =2 Omega_eta+x[-b log(nm)+chi_n+chi_m].
```

The ordered logarithmic rates therefore multiply unrestricted sine and
cosine correlations. Rate ordering alone supplies no sign.

## Moment Collapse

Let

```text
H_0=sum_n q_n,
H_1=sum_n log(n)q_n,
D_0=sum_n delta_nq_n.
```

Summing every diagonal and pair row gives exactly

```text
K_bulk
 =Omega_eta|H_0|^2
  -xc Im[H_1 conjugate(H_0)]
  -xb Re[H_1 conjugate(H_0)]
  +x Im[D_0 conjugate(H_0)].
```

The checker reconstructs both sides over exact rational complex
arithmetic. The recurrent endpoint and its cross terms still have to be
added before an endpoint-complete estimate.

## Heat-Shifted P-Free Telescope

For

```text
A_t(n;s)=exp[t log(n)^2/4-s log(n)],
h=log(p),
M_k=floor(N/p^k),
sigma_k=s-t k h/2,
```

put

```text
T_k=A_t(p^k;s)S_(M_k,t)(sigma_k),
L_k=A_t(p^k;s)O_(M_k,t)^(p)(sigma_k).
```

The exact heat-shift identity

```text
A_t(p^k;s)A_t(p;s-t k h/2)=A_t(p^(k+1);s)
```

gives

```text
L_k=T_k-T_(k+1),
sum_k L_k=T_0=S_(N,t)(s).
```

This is a genuine linear telescope. The phase current is quadratic:

```text
K(T_0)
 =sum_k K(L_k)+sum_(k<l)K_(k,l).
```

Thus retaining every cross term reconstructs the original full-prefix
moment current. P-free reindexing is exact, but it does not create a
new algebraic sign.

## Exact Cancellation Guard

One two-layer exact jet has

```text
L_0=i,       D L_0=-1,
L_1=-i/3,    D L_1=1.
```

It telescopes linearly to `T_0=2i/3`, `D T_0=0`, while

```text
sum_k K(L_k)={witness["diagonal_current_sum"]},
K_(0,1)={witness["pair_cross_current"]},
K(T_0)={witness["joined_current"]}.
```

The cross current cancels the entire positive diagonal reserve. This is
the same cubic-contact algebra seen by the preceding polarization gate,
now written as an exact telescoping-layer audit. It is not an Xi
counterexample.

## Route Decision

The old prime symmetry remains useful as an exact organization of the
coefficients, but its linear telescope does not reduce the quadratic
joined current below the original `H_0,H_1,D_0` correlation. Any further
gain must come from a new Xi-specific estimate for those oscillatory
correlations together with the recurrent endpoint. Without such a
theorem, the more robust target is the direct endpoint-complete,
division-free Abel-scalar gap.

## Pi Provenance

No `pi` enters the pair-kernel or p-free telescope identities. Any `pi`
in the Xi saddle scale remains the ordinary constant inherited from
completed-zeta normalization.

## Boundary

{boundary}
"""


def main() -> int:
    rows = gate_rows()
    proof_boundary = (
        "This gate proves the physical division-free pair kernel, exact "
        "bulk moment collapse, correction-free heat-shifted p-free "
        "telescope, quadratic rejoin, and strict cross-cancellation guard. "
        "It proves no Xi-specific pair-correlation sign, recurrent-endpoint "
        "absorption, Abel-scalar gap, horizontal successor winding cap, "
        "contact exclusion, Q209, cofinal descendant theorem, Lambda<=0, "
        "PF-infinity, RH, or prize-level conclusion."
    )
    payload = {
        "kind": "joined_pair_kernel_pfree_rejoin_gate",
        "date": "2026-07-30",
        "status": (
            "exact physical pair kernel and p-free quadratic-rejoin audit"
        ),
        "proof_boundary": proof_boundary,
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "exact": {
            "pair_kernel": (
                "K_(n,m)=(a_n-a_m)Im(q_nconj(q_m))"
                "+(omega_n+omega_m)Re(q_nconj(q_m))"
            ),
            "moment_collapse": (
                "K_bulk=Omega_eta|H_0|^2-xc Im(H_1conj(H_0))"
                "-xb Re(H_1conj(H_0))+x Im(D_0conj(H_0))"
            ),
            "heat_telescope": "L_k=T_k-T_(k+1); sum_kL_k=T_0",
            "quadratic_rejoin": (
                "K(T_0)=sum_kK(L_k)+sum_(k<l)K_(k,l)"
            ),
        },
        "witnesses": {
            "physical_kernel": physical_kernel_witness(),
            "heat_shift": heat_shift_witness(),
            "quadratic_rejoin": quadratic_rejoin_witness(),
        },
        "summary": {
            "rows": len(rows),
            "exact_pair_kernel_witnesses": 1,
            "exact_moment_collapses": 1,
            "heat_shifted_pfree_telescopes": 1,
            "quadratic_rejoin_witnesses": 1,
            "strict_cross_cancellation_guards": 1,
            "proved_cross_current_signs": 0,
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
        "built joined pair-kernel p-free rejoin gate: "
        "1 exact pair kernel, 1 moment collapse, 1 heat telescope, "
        "1 quadratic rejoin, 1 strict cross-cancellation guard"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
