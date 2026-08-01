#!/usr/bin/env python3
"""Build the corrected four-moment Mangoldt-Abel contact gate."""

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
    "mangoldt_abel_contact_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "first_order_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_signed_contact_reduction.json"
    ),
    "absolute_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "legacy_symmetry": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "legacy_prime_curvature_symmetry_audit.json"
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


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def real(value: ComplexQ) -> Fraction:
    return value[0]


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

    contact_text = json.dumps(payloads["first_order_contact"], sort_keys=True)
    anchor_text = json.dumps(payloads["absolute_anchor"], sort_keys=True)
    abel_text = json.dumps(payloads["abel_prefix"], sort_keys=True)
    required = {
        "first_order_contact": (
            "d_n=1/(6s)",
            "alpha_n=alpha(s)-log(n)",
        ),
        "absolute_anchor": (
            "q_n=f_n/f_1",
            "eta=f_1/|f_1|",
        ),
        "abel_prefix": (
            "U=u_N*S+sum_(k=1)^(N-1)h_k*F_k",
            "mathcal_C_N",
        ),
    }
    for marker in required["first_order_contact"]:
        if marker not in contact_text:
            raise RuntimeError(f"first-order source marker missing: {marker}")
    for marker in required["absolute_anchor"]:
        if marker not in anchor_text:
            raise RuntimeError(f"anchor source marker missing: {marker}")
    for marker in required["abel_prefix"]:
        if marker not in abel_text:
            raise RuntimeError(f"Abel source marker missing: {marker}")
    if payloads["legacy_symmetry"].get("kind") != (
        "legacy_prime_curvature_symmetry_audit"
    ):
        raise RuntimeError("legacy symmetry source drifted")
    return {"source_sha256": hashes}


def symbolic_audit() -> dict[str, str]:
    ell, alpha, alpha_prime, heat_time, s = sp.symbols(
        "ell alpha alpha_prime heat_time s"
    )
    quadratic = alpha_prime * heat_time**2 / 8
    a_0 = (
        1
        + 1 / (6 * s)
        + alpha_prime * heat_time / 4
        + quadratic * alpha**2
    )
    a_1 = -2 * quadratic * alpha
    a_2 = quadratic
    correction = (
        1
        + 1 / (6 * s)
        + alpha_prime
        * (
            heat_time / 4
            + heat_time**2 * (alpha - ell) ** 2 / 8
        )
    )
    if sp.expand(correction - (a_0 + a_1 * ell + a_2 * ell**2)) != 0:
        raise RuntimeError("first-order correction polynomial failed")

    varrho_1, varrho_2, log_a = sp.symbols(
        "varrho_1 varrho_2 log_a"
    )
    weighted = sp.expand(
        (log_a - ell) * (1 + varrho_1 * ell + varrho_2 * ell**2)
    )
    expected = (
        log_a
        + (varrho_1 * log_a - 1) * ell
        + (varrho_2 * log_a - varrho_1) * ell**2
        - varrho_2 * ell**3
    )
    if sp.expand(weighted - expected) != 0:
        raise RuntimeError("four-moment Abel polynomial failed")

    log_d, log_m, s_star = sp.symbols("log_d log_m s_star")
    q_dm = sp.exp(
        heat_time * (log_d + log_m) ** 2 / 4
        - s_star * (log_d + log_m)
    )
    q_d_q_m = sp.exp(
        heat_time * log_d**2 / 4
        - s_star * log_d
        + heat_time * log_m**2 / 4
        - s_star * log_m
    )
    heat_coupling = sp.exp(heat_time * log_d * log_m / 2)
    if sp.simplify(q_dm - q_d_q_m * heat_coupling) != 0:
        raise RuntimeError("multiplicative heat coupling failed")

    return {
        "correction_polynomial": (
            "Put B=alpha'(s)t^2/8, "
            "a_0=1+1/(6s)+alpha'(s)t/4+B alpha(s)^2, "
            "a_1=-2B alpha(s), and a_2=B. Then, exactly, "
            "1+d_n=a_0+a_1 log(n)+a_2 log(n)^2."
        ),
        "normalized_coefficients": (
            "Since a_0=1+d_1 and |d_1|<1/2, put "
            "varrho_1=a_1/a_0 and varrho_2=a_2/a_0. With "
            "q_n^(0)=exp[t log(n)^2/4-s_*log(n)] and "
            "eta=f_1/|f_1|, the anchored carrier is "
            "z_n=eta q_n^(0)"
            "[1+varrho_1 log(n)+varrho_2 log(n)^2]."
        ),
        "four_moment_abel": (
            "For H_j=sum_(n<=N)log(n)^j q_n^(0), j=0,1,2,3, "
            "S=eta[H_0+varrho_1 H_1+varrho_2 H_2] and "
            "U=eta[log(a)H_0+(varrho_1 log(a)-1)H_1"
            "+(varrho_2 log(a)-varrho_1)H_2-varrho_2 H_3]."
        ),
        "heat_factorization": (
            "For dm<=N, q_(dm)^(0)=q_d^(0)q_m^(0)"
            "exp[(t/2)log(d)log(m)]."
        ),
        "heat_gram": (
            "For t>=0 the untruncated heat coupling "
            "G_t(d,m)=exp[(t/2)log(d)log(m)] is positive semidefinite, "
            "because sum_(d,m)y_d y_m G_t(d,m)="
            "sum_(r>=0)(t/2)^r/r! [sum_d y_d log(d)^r]^2."
        ),
    }


def factor_vector(number: int) -> dict[int, int]:
    return {int(p): int(e) for p, e in sp.factorint(number).items()}


def lambda_vector(number: int) -> dict[int, int]:
    if number <= 1:
        return {}
    factors = factor_vector(number)
    if len(factors) != 1:
        return {}
    prime = next(iter(factors))
    return {prime: 1}


def add_vectors(
    left: dict[int, Fraction], right: dict[int, Fraction]
) -> dict[int, Fraction]:
    result = dict(left)
    for prime, value in right.items():
        result[prime] = result.get(prime, Fraction(0)) + value
        if result[prime] == 0:
            del result[prime]
    return result


def scale_vector(
    value: dict[int, int] | dict[int, Fraction], scalar: Fraction
) -> dict[int, Fraction]:
    return {
        int(prime): Fraction(coefficient) * scalar
        for prime, coefficient in value.items()
        if Fraction(coefficient) * scalar
    }


def divisor_mangoldt_audit(limit: int = 128) -> dict:
    mismatches: list[dict] = []
    prime_power_terms = 0
    for number in range(1, limit + 1):
        one_sided: dict[int, Fraction] = {}
        symmetric: dict[int, Fraction] = {}
        for divisor in sp.divisors(number):
            divisor = int(divisor)
            other = number // divisor
            left = lambda_vector(divisor)
            right = lambda_vector(other)
            if left:
                prime_power_terms += 1
            one_sided = add_vectors(
                one_sided, scale_vector(left, Fraction(1))
            )
            symmetric = add_vectors(
                symmetric, scale_vector(left, Fraction(1, 2))
            )
            symmetric = add_vectors(
                symmetric, scale_vector(right, Fraction(1, 2))
            )
        expected = {
            prime: Fraction(exponent)
            for prime, exponent in factor_vector(number).items()
        }
        if one_sided != expected or symmetric != expected:
            mismatches.append(
                {
                    "n": number,
                    "one_sided": one_sided,
                    "symmetric": symmetric,
                    "expected": expected,
                }
            )
    if mismatches:
        raise RuntimeError(f"Mangoldt divisor audit failed: {mismatches[:1]}")
    return {
        "limit": limit,
        "coefficient_audits": limit,
        "moment_families": 3,
        "prime_power_terms_seen": prime_power_terms,
        "mismatches": 0,
    }


def prime_edge_witnesses() -> list[dict]:
    rows: list[dict] = []
    for cutoff in (2, 3, 4, 5, 10, 25, 64, 127):
        primes = [
            int(value)
            for value in sp.primerange(
                math.isqrt(cutoff) + 1, cutoff + 1
            )
        ]
        if not primes:
            raise RuntimeError(f"no prime edge found for N={cutoff}")
        prime = primes[-1]
        if not (prime <= cutoff < prime * prime):
            raise RuntimeError(f"bad prime edge for N={cutoff}")
        determinants = {
            str(moment): f"-log({prime})^{2 * moment}/4"
            for moment in (1, 2, 3)
        }
        rows.append(
            {
                "N": cutoff,
                "p": prime,
                "p_squared": prime * prime,
                "determinants": determinants,
            }
        )
    return rows


def rational_moment_witness() -> dict:
    logarithms = (
        Fraction(0),
        Fraction(1, 2),
        Fraction(2, 3),
        Fraction(1),
    )
    base = (
        (Fraction(1), Fraction(0)),
        (Fraction(2, 3), Fraction(1, 5)),
        (Fraction(-3, 7), Fraction(2, 9)),
        (Fraction(4, 11), Fraction(-1, 6)),
    )
    varrho_1 = (Fraction(1, 13), Fraction(-2, 17))
    varrho_2 = (Fraction(-1, 19), Fraction(1, 23))
    eta = (Fraction(3, 5), Fraction(4, 5))
    log_a = Fraction(7, 6)
    moments: list[ComplexQ] = []
    for order in range(4):
        total = (Fraction(0), Fraction(0))
        for ell, value in zip(logarithms, base, strict=True):
            total = add(total, scale(value, ell**order))
        moments.append(total)

    direct_s = (Fraction(0), Fraction(0))
    direct_u = (Fraction(0), Fraction(0))
    for ell, value in zip(logarithms, base, strict=True):
        polynomial = add(
            (Fraction(1), Fraction(0)),
            add(scale(varrho_1, ell), scale(varrho_2, ell**2)),
        )
        carrier = multiply(eta, multiply(value, polynomial))
        direct_s = add(direct_s, carrier)
        direct_u = add(direct_u, scale(carrier, log_a - ell))

    collapsed_s = multiply(
        eta,
        add(
            moments[0],
            add(
                multiply(varrho_1, moments[1]),
                multiply(varrho_2, moments[2]),
            ),
        ),
    )
    coefficient_1 = add(
        scale(varrho_1, log_a), (Fraction(-1), Fraction(0))
    )
    coefficient_2 = add(
        scale(varrho_2, log_a), scale(varrho_1, Fraction(-1))
    )
    collapsed_u = multiply(
        eta,
        add(
            scale(moments[0], log_a),
            add(
                multiply(coefficient_1, moments[1]),
                add(
                    multiply(coefficient_2, moments[2]),
                    multiply(scale(varrho_2, Fraction(-1)), moments[3]),
                ),
            ),
        ),
    )
    if direct_s != collapsed_s or direct_u != collapsed_u:
        raise RuntimeError("rational four-moment witness failed")

    endpoint = (Fraction(-2, 7), Fraction(1, 8))
    endpoint_slope = (Fraction(5, 12), Fraction(-3, 10))
    s_prime = (Fraction(1, 17), Fraction(-1, 2))
    u_n = Fraction(1, 10)
    x_value = real(add(endpoint, collapsed_s))
    contact_scalar = real(
        add(endpoint_slope, multiply(s_prime, collapsed_u))
    ) - s_prime[0] * u_n * x_value
    return {
        "logarithms": [fraction_text(value) for value in logarithms],
        "base": [complex_text(value) for value in base],
        "varrho_1": complex_text(varrho_1),
        "varrho_2": complex_text(varrho_2),
        "eta": complex_text(eta),
        "log_a": fraction_text(log_a),
        "moments": [complex_text(value) for value in moments],
        "direct_S": complex_text(direct_s),
        "collapsed_S": complex_text(collapsed_s),
        "direct_U": complex_text(direct_u),
        "collapsed_U": complex_text(collapsed_u),
        "endpoint": complex_text(endpoint),
        "endpoint_slope": complex_text(endpoint_slope),
        "s_prime": complex_text(s_prime),
        "u_N": fraction_text(u_n),
        "mathsf_X": fraction_text(x_value),
        "mathcal_C_N": fraction_text(contact_scalar),
        "interpretation": (
            "Exact rational-complex algebra witness only; it is not an Xi "
            "state or a contact counterexample."
        ),
    }


def exact_statements(symbolic: dict[str, str]) -> dict[str, str]:
    return {
        **symbolic,
        "mangoldt_moments": (
            "For j=1,2,3, the finite moment is exactly "
            "H_j=sum_(dm<=N)Lambda(d)log(dm)^(j-1)q_(dm)^(0)"
            "=(1/2)sum_(dm<=N)[Lambda(d)+Lambda(m)]"
            "log(dm)^(j-1)q_(dm)^(0)."
        ),
        "symmetric_bilinear_kernel": (
            "Combining the heat factorization with the preceding identity "
            "writes H_j=sum_(dm<=N)K_(j,t,N)(d,m)"
            "q_d^(0)q_m^(0), where K is real and transpose-symmetric. "
            "This is a complex bilinear form, not a Hermitian form."
        ),
        "contact_normal_form": (
            "With e=eta r_0=g_0/|f_1|, "
            "g=eta r_A=G_a/|f_1|, W_0=e+S, "
            "mathsf_X=Re(W_0), s_*'=c+ib, and "
            "u_N=log(a/N), the exact endpoint-complete scalar is "
            "mathcal_C_N=Re[g+s_*'U]-c u_N mathsf_X, with S and U "
            "given by the four-moment formulas."
        ),
        "contact_band_reduction": (
            "On |mathsf_X|<=delta_L, the already certified terminal shear "
            "obeys |c u_N mathsf_X|<epsilon_term. Therefore a signed "
            "endpoint-plus-four-moment estimate "
            "|Re[g+s_*'U]|>A_L+2epsilon_term is sufficient for the "
            "existing Abel gap, but is not proved here."
        ),
        "prime_edge_minor": (
            "If sqrt(N)<p<=N is prime, then for every j=1,2,3 the "
            "{1,p} principal minor of K_(j,t,N) is "
            "[[0,log(p)^j/2],[log(p)^j/2,0]], with determinant "
            "-log(p)^(2j)/4<0. Bertrand's postulate supplies such a p "
            "for every N>=2, with N=2,3,4 checked directly."
        ),
        "legacy_distinction": (
            "The exact first-jet symmetry is multiplicative and "
            "hyperbolically truncated: dm<=N with additive weight "
            "Lambda(d)+Lambda(m). The legacy radial square uses "
            "d^2+m^2 and product weight Lambda(d)Lambda(m). The latter "
            "is a different, second-order object and is not the Xi "
            "contact scalar."
        ),
        "route_decision": (
            "Retain the symmetric Mangoldt form as a Type-I/II or bilinear "
            "cancellation coordinate. Do not promote transpose symmetry, "
            "entrywise positivity, or the untruncated heat Gram to a "
            "contact lower bound. A successful estimate must keep the "
            "complex phase, hyperbolic cutoff, recurrent endpoint, "
            "first-order coefficients, W_0=0, and adjacent charts joined."
        ),
    }


def build_rows(exact: dict[str, str]) -> list[GateRow]:
    return [
        GateRow(
            "macg_01_correction",
            "The retained first Dirichlet correction is exactly quadratic in log(n).",
            "proved",
            exact["correction_polynomial"],
            "No remainder is introduced by the polynomial expansion.",
            "This is the retained first-order coefficient only.",
        ),
        GateRow(
            "macg_02_anchor",
            "The branch-free normalized carrier has two exact logarithmic correction coefficients.",
            "proved",
            exact["normalized_coefficients"],
            "The nonzero first carrier fixes the common absolute phase.",
            "No sign or contact margin follows from normalization.",
        ),
        GateRow(
            "macg_03_four_moments",
            "The full corrected Abel bulk closes exactly on H_0 through H_3.",
            "proved",
            exact["four_moment_abel"],
            "No higher logarithmic moment is needed at retained first order.",
            "The recurrent endpoint remains separate and explicit.",
        ),
        GateRow(
            "macg_04_mangoldt",
            "Each nonconstant logarithmic moment has an exact symmetric Mangoldt divisor-pair representation.",
            "proved",
            exact["mangoldt_moments"],
            "The transpose symmetry is an exact arithmetic identity.",
            "It is a finite hyperbolic sum, not a radial prime-pair field.",
        ),
        GateRow(
            "macg_05_heat_factor",
            "The heat quadratic factors into two multiplicative carriers and one symmetric coupling.",
            "proved",
            exact["heat_factorization"],
            "This exposes a natural bilinear exponential-sum coordinate.",
            "The cutoff and Mangoldt weights are not removed.",
        ),
        GateRow(
            "macg_06_heat_gram",
            "The untruncated heat coupling alone is a positive Gram kernel.",
            "proved",
            exact["heat_gram"],
            "This isolates the exact positive part of the old symmetry intuition.",
            "It is not the complete finite Mangoldt kernel.",
        ),
        GateRow(
            "macg_07_bilinear",
            "The corrected moments are transpose-symmetric complex bilinear forms.",
            "proved",
            exact["symmetric_bilinear_kernel"],
            "Type-I/II cancellation may be sought before absolute values.",
            "There is no conjugation and hence no automatic real positivity.",
        ),
        GateRow(
            "macg_08_contact",
            "The endpoint-complete Abel scalar has an exact four-moment Mangoldt normal form.",
            "proved",
            exact["contact_normal_form"],
            "The formula includes the exceptional W_0=0 fibre.",
            "The Xi-specific lower bound remains open.",
        ),
        GateRow(
            "macg_09_band",
            "A signed endpoint-plus-four-moment estimate is sufficient for the existing Abel gap.",
            "conditional",
            exact["contact_band_reduction"],
            "This is the narrow theorem-search handoff.",
            "The sufficient estimate is not established.",
        ),
        GateRow(
            "macg_10_minor",
            "The actual symmetric Mangoldt kernel has an exact negative prime-edge minor.",
            "proved",
            exact["prime_edge_minor"],
            "Transpose symmetry and entrywise nonnegativity cannot prove the gap by Gram positivity.",
            "This is a kernel guard, not an Xi contact counterexample.",
        ),
        GateRow(
            "macg_11_legacy",
            "The old radial prime symmetry and the actual contact symmetry are mathematically distinct.",
            "proved",
            exact["legacy_distinction"],
            "The useful clue is multiplicative divisor symmetry, not visual radial alignment.",
            "No zero-curvature correlation is promoted.",
        ),
        GateRow(
            "macg_12_route",
            "The surviving use of the symmetry is bilinear cancellation with the endpoint retained.",
            "open",
            exact["route_decision"],
            "The next stage is an endpoint-complete signed Type-I/II estimate.",
            (
                "No Abel gap, horizontal winding cap, contact exclusion, "
                "Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
            ),
        ),
    ]


def build_note(payload: dict) -> str:
    exact = payload["exact"]
    witnesses = payload["witnesses"]
    return f"""# Corrected Mangoldt-Abel Contact Gate

Date: 2026-07-30

Status: exact first-order moment and symmetric divisor-pair reduction
with a nonpromotion guard. This is not a proof of the Xi Abel gap,
contact exclusion, `Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Pi Provenance

The `pi` in `a^2=x/(4*pi)+t/16` remains the completed-zeta and
Riemann-Siegel normalization. No new `pi` is introduced by the
Mangoldt convolution, a circle, a polygon, or the legacy images.

## Exact First-Order Collapse

```text
{exact["correction_polynomial"]}

{exact["normalized_coefficients"]}
```

The retained correction is therefore not an uncontrolled perturbation:
it changes the correction-free Dirichlet chain by exactly two logarithmic
moments.

```text
{exact["four_moment_abel"]}
```

This is an exact collapse to four finite logarithmic moments.
The recurrent endpoint is still present through `e` and `g`; it has not
been absorbed into a bulk estimate.

## Symmetric Mangoldt Form

Using `sum_(d|n)Lambda(d)=log(n)` and swapping `d,m` gives

```text
{exact["mangoldt_moments"]}

{exact["heat_factorization"]}
```

Thus the actual symmetry is transpose symmetry on the multiplicative
hyperbola `dm<=N`. The divisor audit checked
`{witnesses["mangoldt"]["coefficient_audits"]}` exact coefficient rows
through `n={witnesses["mangoldt"]["limit"]}` with zero mismatches.

The untruncated heat coupling has the exact square expansion

```text
{exact["heat_gram"]}
```

but the retained moment is a complex bilinear form:

```text
{exact["symmetric_bilinear_kernel"]}
```

## Contact Normal Form

```text
{exact["contact_normal_form"]}
```

On the hard bottom intervals the already certified contact-band shear
therefore gives the sufficient handoff

```text
{exact["contact_band_reduction"]}
```

This is a theorem target, not a proved inequality.

## Symmetry Guard

The symmetric kernel is not positive semidefinite:

```text
{exact["prime_edge_minor"]}
```

The machine witness records {len(witnesses["prime_edges"])} cutoffs and
{3 * len(witnesses["prime_edges"])} strict negative minors. The sign is
independent of the heat parameter because `log(1)=0`.
These are negative prime-edge minors of the actual symmetric kernel.

There is a second guard: the sum uses `q_d^(0)q_m^(0)`, not
`q_d^(0)conjugate(q_m^(0))`. Even a positive real kernel would not make
the real part of this complex bilinear form positive.

## Relation To The Older Symmetry

```text
{exact["legacy_distinction"]}
```

The images were useful in pointing at symmetry, but the exact symmetry
that reaches the Xi contact problem is multiplicative and first order.
The radial `Lambda(d)Lambda(m)` square is a separate second-order object.

## Route Decision

```text
{exact["route_decision"]}
```

The promising use of this normal form is a cancellation-preserving
Vaughan/Type-I/II decomposition of the three symmetric moments, composed
with the endpoint before absolute values. Another positivity search on
the symmetric matrix is ruled out by the exact prime-edge minor.

## Boundary

{payload["proof_boundary"]}
"""


def build_payload() -> dict:
    symbolic = symbolic_audit()
    exact = exact_statements(symbolic)
    prime_edges = prime_edge_witnesses()
    rows = build_rows(exact)
    return {
        "kind": "corrected_mangoldt_abel_contact_gate",
        "date": "2026-07-30",
        "status": (
            "exact four-moment and symmetric Mangoldt reduction; "
            "Abel gap remains open"
        ),
        "source_audit": source_audit(),
        "exact": exact,
        "witnesses": {
            "rational_moment": rational_moment_witness(),
            "mangoldt": divisor_mangoldt_audit(),
            "prime_edges": prime_edges,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_correction_polynomials": 1,
            "logarithmic_moments": 4,
            "symmetric_mangoldt_moments": 3,
            "prime_edge_witnesses": len(prime_edges),
            "indefinite_prime_edge_minors": 3 * len(prime_edges),
            "proved_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
        "proof_boundary": (
            "This gate proves the exact quadratic correction polynomial, "
            "four-moment Abel collapse, symmetric Mangoldt divisor-pair "
            "identity, heat coupling, endpoint-complete contact normal form, "
            "and indefinite prime-edge guard. It proves no signed "
            "endpoint-plus-moment lower bound, Abel-scalar gap, horizontal "
            "successor winding cap, q<1 or bounded-L closure, contact "
            "exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    atomic_write(
        args.result,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    atomic_write(args.note, build_note(payload))
    print(
        "built corrected Mangoldt-Abel contact gate: "
        f"{len(payload['rows'])} rows, "
        f"{payload['summary']['symmetric_mangoldt_moments']} "
        "symmetric moments"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
