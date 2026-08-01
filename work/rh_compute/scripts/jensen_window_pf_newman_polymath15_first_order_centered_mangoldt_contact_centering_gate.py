#!/usr/bin/env python3
"""Build the contact-centered Mangoldt/Abel equivalence gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "mangoldt_contact_centering_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "mangoldt_abel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_abel_contact_gate.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "absolute_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
}

ComplexQ = tuple[Fraction, Fraction]
LogVector = dict[int, ComplexQ]


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


def negate(value: ComplexQ) -> ComplexQ:
    return -value[0], -value[1]


def scale(value: ComplexQ, scalar: Fraction) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def sum_complex(values: list[ComplexQ]) -> ComplexQ:
    total = (Fraction(0), Fraction(0))
    for value in values:
        total = add(total, value)
    return total


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def source_audit() -> dict:
    payloads: dict[str, dict] = {}
    hashes: dict[str, str] = {}
    for key, path in SOURCES.items():
        if not path.is_file():
            raise RuntimeError(f"missing source {key}: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
        hashes[key] = file_hash(path)

    mangoldt_text = json.dumps(payloads["mangoldt_abel"], sort_keys=True)
    abel_text = json.dumps(payloads["abel_prefix"], sort_keys=True)
    anchor_text = json.dumps(payloads["absolute_anchor"], sort_keys=True)
    markers = {
        "mangoldt": (
            "varrho_1",
            "H_j=sum_(dm<=N)Lambda(d)",
            "mathcal_C_N=Re[g+s_*'U]",
        ),
        "abel": (
            "U=u_N*S+sum_(k=1)^(N-1)h_k*F_k",
            "V_N=g-s_*'*u_N*e",
        ),
        "anchor": (
            "eta=f_1/|f_1|",
            "q_n=f_n/f_1",
        ),
    }
    for marker in markers["mangoldt"]:
        if marker not in mangoldt_text:
            raise RuntimeError(f"Mangoldt source marker missing: {marker}")
    for marker in markers["abel"]:
        if marker not in abel_text:
            raise RuntimeError(f"Abel source marker missing: {marker}")
    for marker in markers["anchor"]:
        if marker not in anchor_text:
            raise RuntimeError(f"anchor source marker missing: {marker}")
    return {"source_sha256": hashes}


def symbolic_audit() -> dict[str, str]:
    h_0, h_1, h_2, h_3 = sp.symbols("H_0 H_1 H_2 H_3")
    varrho_1, varrho_2 = sp.symbols("varrho_1 varrho_2")
    log_a = sp.symbols("log_a", real=True)
    eta = sp.symbols("eta")
    bulk = eta * (h_0 + varrho_1 * h_1 + varrho_2 * h_2)
    abel = eta * (
        log_a * h_0
        + (varrho_1 * log_a - 1) * h_1
        + (varrho_2 * log_a - varrho_1) * h_2
        - varrho_2 * h_3
    )
    mangoldt = eta * (h_1 + varrho_1 * h_2 + varrho_2 * h_3)
    if sp.expand(abel - log_a * bulk + mangoldt) != 0:
        raise RuntimeError("contact-centering polynomial failed")

    c, b, x_value, y_value, log_n = sp.symbols(
        "c b mathsf_X mathsf_Y log_N", real=True
    )
    s_prime = c + sp.I * b
    w_0 = x_value + sp.I * y_value
    u_n = log_a - log_n
    shear = sp.re(sp.expand(s_prime * log_a * w_0)) - c * u_n * x_value
    expanded = c * log_n * x_value - b * log_a * y_value
    if sp.simplify(shear - expanded) != 0:
        raise RuntimeError("Cartesian contact centering failed")

    return {
        "physical_log_moment": (
            "Put Z_Lambda=eta[H_1+varrho_1 H_2+varrho_2 H_3]. "
            "Then Z_Lambda=sum_(n<=N)log(n)z_n"
            "=sum_(dm<=N)Lambda(d)z_(dm)"
            "=(1/2)sum_(dm<=N)[Lambda(d)+Lambda(m)]z_(dm)."
        ),
        "abel_centering": (
            "With S=sum_(n<=N)z_n and L_a=log(a), "
            "U=L_a S-Z_Lambda=L_a(W_0-e)-Z_Lambda."
        ),
        "centered_jet": (
            "Define E_N=g-s_*'L_a e and "
            "mathfrak B_N=E_N-s_*'Z_Lambda+s_*'L_a W_0. "
            "Then mathfrak B_N=g+s_*'U exactly."
        ),
        "centered_contact_scalar": (
            "The endpoint-complete scalar is "
            "mathcal_C_N=Re(mathfrak B_N)-c u_N mathsf_X. "
            "Equivalently, for W_0=mathsf_X+i mathsf_Y and "
            "u_N=L_a-log(N), mathcal_C_N="
            "Re[E_N-s_*'Z_Lambda]+c log(N)mathsf_X"
            "-b L_a mathsf_Y."
        ),
        "zero_fibre": (
            "At W_0=0, mathfrak B_N=E_N-s_*'Z_Lambda and "
            "mathcal_C_N=Re[E_N-s_*'Z_Lambda]. "
            "No division by W_0, H_a, or a carrier is used."
        ),
        "real_crossing": (
            "At mathsf_X=0, mathcal_C_N="
            "Re[E_N-s_*'Z_Lambda]-b L_a mathsf_Y. "
            "The unrestricted imaginary fibre remains explicit."
        ),
        "abel_duality": (
            "For F_k=sum_(n<=k)z_n and h_k=log((k+1)/k), "
            "Z_Lambda=log(N)S-sum_(k=1)^(N-1)h_kF_k. "
            "Hence at W_0=0, E_N-s_*'Z_Lambda="
            "V_N+s_*'sum_k h_kF_k, where V_N=g-s_*'u_Ne."
        ),
        "endpoint_factor": (
            "With kappa=(-1)^N B_0 real, "
            "e=kappa(T_0+i)H_a and g=kappa(T_0+i)J_a, "
            "E_N=kappa(T_0+i)[J_a-s_*'L_a H_a]. "
            "This identity remains division-free when H_a=0."
        ),
        "balanced_hyperbola": (
            "Let D=floor(sqrt(N)). Since no pair dm<=N has d>D and "
            "m>D, Z_Lambda equals "
            "(1/2)sum_(d,m<=D)[Lambda(d)+Lambda(m)]z_(dm)"
            "+sum_(d<=D<m<=N/d)[Lambda(d)+Lambda(m)]z_(dm)."
        ),
        "band_handoff": (
            "On |mathsf_X|<=delta_L, the existing shear bound "
            "|c u_N mathsf_X|<epsilon_term shows that "
            "|Re(mathfrak B_N)|>A_L+2epsilon_term is sufficient for the "
            "Abel gap. This estimate is not proved here."
        ),
        "route_decision": (
            "Treat Z_Lambda, the Abel-prefix sum, and the centered jet as "
            "three exact coordinates of one object. Apply a Type-I/II or "
            "Vaughan decomposition to Z_Lambda only after composing E_N "
            "and s_*'L_aW_0, and retain the balanced-wing boundary. "
            "Do not count the coordinates as independent evidence."
        ),
    }


def deterministic_carrier(number: int) -> ComplexQ:
    return (
        Fraction((7 * number) % 19 - 9, number + 5),
        Fraction((11 * number) % 23 - 11, number + 7),
    )


def add_log_vector(
    target: LogVector,
    source: dict[int, int],
    value: ComplexQ,
    scalar: Fraction = Fraction(1),
) -> None:
    for prime, exponent in source.items():
        contribution = scale(value, scalar * exponent)
        target[prime] = add(
            target.get(prime, (Fraction(0), Fraction(0))),
            contribution,
        )
        if target[prime] == (Fraction(0), Fraction(0)):
            del target[prime]


def factor_vector(number: int) -> dict[int, int]:
    return {
        int(prime): int(exponent)
        for prime, exponent in sp.factorint(number).items()
    }


def lambda_vector(number: int) -> dict[int, int]:
    factors = factor_vector(number)
    if number <= 1 or len(factors) != 1:
        return {}
    return {next(iter(factors)): 1}


def hyperbola_vector(cutoff: int) -> LogVector:
    result: LogVector = {}
    for divisor in range(1, cutoff + 1):
        lam = lambda_vector(divisor)
        if not lam:
            continue
        for other in range(1, cutoff // divisor + 1):
            add_log_vector(
                result,
                lam,
                deterministic_carrier(divisor * other),
            )
    return result


def direct_vector(cutoff: int) -> LogVector:
    result: LogVector = {}
    for number in range(1, cutoff + 1):
        add_log_vector(
            result,
            factor_vector(number),
            deterministic_carrier(number),
        )
    return result


def balanced_vector(cutoff: int) -> tuple[LogVector, int, int]:
    split = math.isqrt(cutoff)
    result: LogVector = {}
    central_pairs = 0
    wing_pairs = 0
    for divisor in range(1, split + 1):
        for other in range(1, split + 1):
            value = deterministic_carrier(divisor * other)
            add_log_vector(
                result,
                lambda_vector(divisor),
                value,
                Fraction(1, 2),
            )
            add_log_vector(
                result,
                lambda_vector(other),
                value,
                Fraction(1, 2),
            )
            central_pairs += 1
        for other in range(split + 1, cutoff // divisor + 1):
            value = deterministic_carrier(divisor * other)
            add_log_vector(result, lambda_vector(divisor), value)
            add_log_vector(result, lambda_vector(other), value)
            wing_pairs += 1
    return result, central_pairs, wing_pairs


def hyperbola_audit() -> list[dict]:
    rows: list[dict] = []
    for cutoff in (2, 3, 4, 5, 10, 17, 31, 64):
        direct = direct_vector(cutoff)
        one_sided = hyperbola_vector(cutoff)
        balanced, central_pairs, wing_pairs = balanced_vector(cutoff)
        if direct != one_sided or direct != balanced:
            raise RuntimeError(f"balanced hyperbola audit failed at N={cutoff}")
        rows.append(
            {
                "N": cutoff,
                "D": math.isqrt(cutoff),
                "log_prime_coordinates": len(direct),
                "central_pairs": central_pairs,
                "wing_pairs": wing_pairs,
                "mismatches": 0,
            }
        )
    return rows


def rational_centering_witness() -> dict:
    h_values = (
        (Fraction(5, 7), Fraction(-1, 4)),
        (Fraction(-2, 9), Fraction(3, 8)),
        (Fraction(4, 11), Fraction(1, 6)),
        (Fraction(-3, 13), Fraction(-2, 5)),
    )
    varrho_1 = (Fraction(1, 13), Fraction(-2, 17))
    varrho_2 = (Fraction(-1, 19), Fraction(1, 23))
    eta = (Fraction(3, 5), Fraction(4, 5))
    e = (Fraction(-2, 7), Fraction(1, 8))
    g = (Fraction(5, 12), Fraction(-3, 10))
    s_prime = (Fraction(1, 17), Fraction(-1, 2))
    log_a = Fraction(7, 6)
    log_n = Fraction(5, 6)
    u_n = log_a - log_n

    h_0, h_1, h_2, h_3 = h_values
    s_bulk = multiply(
        eta,
        add(
            h_0,
            add(
                multiply(varrho_1, h_1),
                multiply(varrho_2, h_2),
            ),
        ),
    )
    u_bulk = multiply(
        eta,
        add(
            scale(h_0, log_a),
            add(
                multiply(
                    add(
                        scale(varrho_1, log_a),
                        (Fraction(-1), Fraction(0)),
                    ),
                    h_1,
                ),
                add(
                    multiply(
                        add(
                            scale(varrho_2, log_a),
                            negate(varrho_1),
                        ),
                        h_2,
                    ),
                    multiply(negate(varrho_2), h_3),
                ),
            ),
        ),
    )
    z_lambda = multiply(
        eta,
        add(
            h_1,
            add(
                multiply(varrho_1, h_2),
                multiply(varrho_2, h_3),
            ),
        ),
    )
    if u_bulk != add(scale(s_bulk, log_a), negate(z_lambda)):
        raise RuntimeError("rational U centering failed")

    w_0 = add(e, s_bulk)
    endpoint_center = add(g, negate(multiply(s_prime, scale(e, log_a))))
    centered_jet = add(
        endpoint_center,
        add(
            negate(multiply(s_prime, z_lambda)),
            multiply(s_prime, scale(w_0, log_a)),
        ),
    )
    direct_jet = add(g, multiply(s_prime, u_bulk))
    if centered_jet != direct_jet:
        raise RuntimeError("rational centered jet failed")

    contact_direct = direct_jet[0] - s_prime[0] * u_n * w_0[0]
    contact_expanded = (
        add(endpoint_center, negate(multiply(s_prime, z_lambda)))[0]
        + s_prime[0] * log_n * w_0[0]
        - s_prime[1] * log_a * w_0[1]
    )
    if contact_direct != contact_expanded:
        raise RuntimeError("rational contact scalar centering failed")

    zero_e = negate(s_bulk)
    zero_endpoint_center = add(
        g,
        negate(multiply(s_prime, scale(zero_e, log_a))),
    )
    zero_jet = add(
        zero_endpoint_center,
        negate(multiply(s_prime, z_lambda)),
    )
    if zero_jet != direct_jet:
        raise RuntimeError("zero-fibre centered jet failed")

    return {
        "H": [complex_text(value) for value in h_values],
        "varrho_1": complex_text(varrho_1),
        "varrho_2": complex_text(varrho_2),
        "eta": complex_text(eta),
        "e": complex_text(e),
        "g": complex_text(g),
        "s_prime": complex_text(s_prime),
        "log_a": fraction_text(log_a),
        "log_N": fraction_text(log_n),
        "u_N": fraction_text(u_n),
        "S": complex_text(s_bulk),
        "U": complex_text(u_bulk),
        "Z_Lambda": complex_text(z_lambda),
        "W_0": complex_text(w_0),
        "direct_jet": complex_text(direct_jet),
        "centered_jet": complex_text(centered_jet),
        "mathcal_C_N": fraction_text(contact_direct),
        "zero_fibre_e": complex_text(zero_e),
        "zero_fibre_jet": complex_text(zero_jet),
        "interpretation": (
            "Exact rational-complex identity witness only; it is not an "
            "Xi state, Abel-gap witness, or contact counterexample."
        ),
    }


def abel_duality_witness() -> dict:
    z_values = [
        (Fraction(1, 3), Fraction(-1, 7)),
        (Fraction(-2, 5), Fraction(1, 11)),
        (Fraction(3, 8), Fraction(2, 13)),
        (Fraction(-1, 4), Fraction(-3, 17)),
        (Fraction(2, 9), Fraction(1, 6)),
    ]
    increments = [
        Fraction(1, 10),
        Fraction(1, 12),
        Fraction(1, 14),
        Fraction(1, 16),
    ]
    log_n = Fraction(7, 6)
    log_a = Fraction(3, 2)
    u_n = log_a - log_n
    prefixes: list[ComplexQ] = []
    running = (Fraction(0), Fraction(0))
    for value in z_values:
        running = add(running, value)
        prefixes.append(running)

    logarithms: list[Fraction] = []
    for index in range(len(z_values)):
        logarithms.append(log_n - sum(increments[index:], Fraction(0)))
    z_lambda = sum_complex(
        [
            scale(value, logarithm)
            for value, logarithm in zip(
                z_values, logarithms, strict=True
            )
        ]
    )
    prefix_sum = sum_complex(
        [
            scale(prefixes[index], increments[index])
            for index in range(len(increments))
        ]
    )
    s_bulk = prefixes[-1]
    dual = add(scale(s_bulk, log_n), negate(prefix_sum))
    u_bulk = sum_complex(
        [
            scale(value, log_a - logarithm)
            for value, logarithm in zip(
                z_values, logarithms, strict=True
            )
        ]
    )
    abel = add(scale(s_bulk, u_n), prefix_sum)
    if z_lambda != dual or u_bulk != abel:
        raise RuntimeError("rational Abel duality failed")
    return {
        "terms": len(z_values),
        "increments": [fraction_text(value) for value in increments],
        "log_N": fraction_text(log_n),
        "log_a": fraction_text(log_a),
        "S": complex_text(s_bulk),
        "Z_Lambda_direct": complex_text(z_lambda),
        "Z_Lambda_prefix": complex_text(dual),
        "U_direct": complex_text(u_bulk),
        "U_abel": complex_text(abel),
        "mismatches": 0,
    }


def build_rows(exact: dict[str, str]) -> list[GateRow]:
    return [
        GateRow(
            "mccg_01_log_moment",
            "The three corrected nonconstant moments combine into one physical logarithmic moment.",
            "proved",
            exact["physical_log_moment"],
            "The correction polynomial stays inside the physical carrier.",
            "This identity alone supplies no sign.",
        ),
        GateRow(
            "mccg_02_abel_centering",
            "The Abel bulk is exactly the endpoint-centered value minus one logarithmic moment.",
            "proved",
            exact["abel_centering"],
            "H_0 and every bulk log(a) coefficient cancel in one step.",
            "The recurrent endpoint is not absorbed.",
        ),
        GateRow(
            "mccg_03_centered_jet",
            "The full directional jet has one centered endpoint, one Mangoldt moment, and one W_0 term.",
            "proved",
            exact["centered_jet"],
            "The four-moment and Abel-prefix routes become one object.",
            "No lower bound is asserted.",
        ),
        GateRow(
            "mccg_04_contact_scalar",
            "The endpoint-complete contact scalar has an exact centered Cartesian form.",
            "proved",
            exact["centered_contact_scalar"],
            "The X and Y fibres remain explicit without division.",
            "The signed real part is not controlled here.",
        ),
        GateRow(
            "mccg_05_zero_fibre",
            "At W_0=0 the H_0 and fibre terms vanish exactly.",
            "proved",
            exact["zero_fibre"],
            "The exceptional fibre reduces to endpoint minus one Mangoldt moment.",
            "This does not prove that the real part is nonzero.",
        ),
        GateRow(
            "mccg_06_real_crossing",
            "At a real crossing the imaginary fibre remains a written term.",
            "proved",
            exact["real_crossing"],
            "Any signed estimate must control or exploit mathsf_Y.",
            "No fixed-sector assumption is introduced.",
        ),
        GateRow(
            "mccg_07_abel_duality",
            "The logarithmic moment and the Abel-prefix sum are exact dual coordinates.",
            "proved",
            exact["abel_duality"],
            "They cannot be counted as independent evidence.",
            "The Xi-specific lower bound remains open.",
        ),
        GateRow(
            "mccg_08_endpoint_factor",
            "The centered endpoint retains the exact recurrent Riemann-Siegel factor.",
            "proved",
            exact["endpoint_factor"],
            "The H_a=0 fibre remains division-free.",
            "No endpoint dominance is claimed.",
        ),
        GateRow(
            "mccg_09_hyperbola",
            "The physical logarithmic moment has an exact balanced square-plus-wing split.",
            "proved",
            exact["balanced_hyperbola"],
            "Every retained pair has at least one index at most sqrt(N).",
            "The wing and floor boundaries remain part of any estimate.",
        ),
        GateRow(
            "mccg_10_band",
            "A centered-jet real-part estimate is sufficient for the existing Abel gap.",
            "conditional",
            exact["band_handoff"],
            "This is the narrowed theorem-search handoff.",
            "The sufficient estimate is not proved.",
        ),
        GateRow(
            "mccg_11_notation",
            "The correction coefficients are distinct from the established derivative rho_1.",
            "proved",
            (
                "Use varrho_1=a_1/a_0 and varrho_2=a_2/a_0 for the "
                "correction polynomial; retain rho_1=partial_x log|f_1| "
                "for the anchored derivative."
            ),
            "The two unrelated quantities cannot be silently identified.",
            "Notation hygiene only.",
        ),
        GateRow(
            "mccg_12_route",
            "The next analytic step is one endpoint-composed bilinear estimate, not parallel coordinate searches.",
            "open",
            exact["route_decision"],
            "A balanced Type-I/II split can now be written without losing W_0=0.",
            (
                "No Abel gap, winding cap, contact exclusion, Lambda<=0, "
                "PF-infinity, RH, or prize-level conclusion is proved."
            ),
        ),
    ]


def build_note(payload: dict) -> str:
    exact = payload["exact"]
    witnesses = payload["witnesses"]
    return f"""# Contact-Centered Mangoldt/Abel Equivalence Gate

Date: 2026-07-30

Status: exact contact centering, Abel duality, and balanced hyperbola
decomposition. The required signed lower bound remains open. This is not
a proof of contact exclusion, `Lambda<=0`, PF-infinity, RH, or a
Clay-prize result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Pi Provenance

The `pi` in `a^2=x/(4*pi)+t/16` is still the completed-zeta and
Riemann-Siegel normalization. Contact centering, Abel summation, the
Mangoldt identity, and the balanced hyperbola introduce no new `pi`.

## One Physical Mangoldt Moment

```text
{exact["physical_log_moment"]}
```

The new coefficients are called `varrho_1,varrho_2`; the established
symbol `rho_1=partial_x log|f_1|` keeps its earlier meaning.

```text
{exact["abel_centering"]}

{exact["centered_jet"]}
```

Thus the four correction-free moments do not remain four independent
objects at contact. They combine into one physical logarithmic derivative
plus the already present value `W_0`.

## Contact Fibres

```text
{exact["centered_contact_scalar"]}

{exact["zero_fibre"]}

{exact["real_crossing"]}
```

The zero fibre is genuinely simpler, but the real-crossing fibre still
contains `mathsf_Y`; no sector or sign has been assumed.

## Abel Duality

```text
{exact["abel_duality"]}
```

The independent rational audit uses
`{witnesses["abel_duality"]["terms"]}` carriers and has zero mismatches.
This identity proves that the prefix and Mangoldt coordinates are two
representations of the same jet, not two sources of evidence.

## Endpoint

```text
{exact["endpoint_factor"]}
```

The formula does not divide by `H_a`, so the endpoint-zero fibre remains
inside the theorem statement.

## Balanced Hyperbola

```text
{exact["balanced_hyperbola"]}
```

The checker reconstructs this split at
`{len(witnesses["hyperbola"])}` cutoffs with exact rational-complex
carriers and zero coefficient mismatches. The small square has no hidden
cutoff because `floor(sqrt(N))^2<=N`; every remaining pair lies in the
written wing. This is an exact organization, not a cancellation bound.

## Live Handoff

```text
{exact["band_handoff"]}

{exact["route_decision"]}
```

The next candidate must estimate the combined endpoint, Mangoldt moment,
and `W_0` term before absolute values. It must retain the wing boundary,
`W_0=0`, `mathsf_X=0`, `H_a=0`, `q=1`, cutoff equality, and adjacent
charts.

## Boundary

This gate proves the contact-centering identity, one physical Mangoldt
moment, exact Abel-prefix duality, recurrent endpoint factorization,
zero-fibre and real-crossing formulas, notation separation, and balanced
hyperbola decomposition. It proves no signed lower bound, Abel-scalar gap,
horizontal successor winding cap, contact exclusion, Q209, cofinal
descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level result.
"""


def build_payload() -> dict:
    exact = symbolic_audit()
    rows = build_rows(exact)
    hyperbola = hyperbola_audit()
    payload = {
        "kind": "mangoldt_contact_centering_gate",
        "date": "2026-07-30",
        "status": (
            "exact contact centering and Mangoldt/Abel equivalence; "
            "signed lower bound remains open"
        ),
        "source_audit": source_audit(),
        "exact": exact,
        "witnesses": {
            "rational_centering": rational_centering_witness(),
            "abel_duality": abel_duality_witness(),
            "hyperbola": hyperbola,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "physical_mangoldt_moments": 1,
            "contact_fibres": 2,
            "abel_duality_witnesses": 1,
            "balanced_hyperbola_audits": len(hyperbola),
            "hyperbola_mismatches": sum(
                int(row["mismatches"]) for row in hyperbola
            ),
            "proved_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
        "proof_boundary": (
            "This gate proves exact contact centering, one physical "
            "Mangoldt moment, Abel-prefix duality, endpoint factorization, "
            "zero-fibre and real-crossing formulas, notation separation, "
            "and a balanced hyperbola split. It proves no signed lower "
            "bound, Abel gap, winding cap, contact exclusion, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
    }
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    args = parser.parse_args()

    payload = build_payload()
    atomic_write(
        args.result,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    atomic_write(args.note, build_note(payload))
    summary = payload["summary"]
    print(
        "built contact-centered Mangoldt gate: "
        f"{summary['rows']} rows, "
        f"{summary['physical_mangoldt_moments']} physical moment, "
        f"{summary['balanced_hyperbola_audits']} hyperbola audits"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
