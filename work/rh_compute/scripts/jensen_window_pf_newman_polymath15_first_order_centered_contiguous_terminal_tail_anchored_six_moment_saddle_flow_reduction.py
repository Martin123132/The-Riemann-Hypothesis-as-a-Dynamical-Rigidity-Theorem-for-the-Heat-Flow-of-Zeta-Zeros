#!/usr/bin/env python3
"""Build the endpoint-composed six-moment saddle-flow reduction."""

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
    "contiguous_terminal_tail_anchored_six_moment_saddle_flow_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "quadratic_barrier": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_q1_saddle_"
        "phase_quadratic_transfer_barrier.json"
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
class SaddleFlowRow:
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
    barrier = payloads["quadratic_barrier"]
    flow = barrier.get("frequency_flow_certificate", {})
    if "one signed integral" not in next(
        (
            row.get("claim", "")
            for row in barrier.get("rows", [])
            if row.get("id") == "qpb_07_line"
        ),
        "",
    ):
        raise RuntimeError("quadratic-flow source drifted")
    if "K=q*r-f*y" not in flow.get("signed_line_integral", ""):
        raise RuntimeError("signed integrand source drifted")
    if barrier.get("counts", {}).get("signed_physical_bounds") != 0:
        raise RuntimeError("quadratic-barrier proof boundary drifted")

    kernel = payloads["two_carrier"].get("symbolic_certificate", {})
    for key, marker in (
        ("full_two_carrier_identity", "h^2*Phi_B"),
        ("hermitian_symmetrization", "H^+_(n,m)"),
        ("transpose_symmetrization", "T^+_(n,m)"),
    ):
        if marker not in kernel.get(key, ""):
            raise RuntimeError(f"two-carrier source drifted: {key}")

    mangoldt = payloads["mangoldt"].get("symbolic_certificate", {})
    if "0<=r<=4" not in mangoldt.get("correction_free_moments", ""):
        raise RuntimeError("five-moment source drifted")
    if "Lambda(d)+Lambda(m)" not in mangoldt.get("mangoldt_moments", ""):
        raise RuntimeError("Mangoldt source drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_signed_flow_bounds": 0,
        "imported_phi_b_bounds": 0,
    }


def full_current_flow_certificate() -> dict:
    count = 3
    u = sp.symbols(f"u0:{count}", real=True)
    z = sp.symbols(f"z0:{count}")
    zb = sp.symbols(f"zb0:{count}")
    linear = sp.symbols(f"L0:{count}")
    hermitian = sp.symbols(f"H0:{count * count}")
    transpose = sp.symbols(f"T0:{count * count}")

    expression = sum(linear[n] * z[n] for n in range(count))
    expression += sp.Rational(1, 2) * sum(
        hermitian[n * count + m] * z[n] * zb[m]
        + transpose[n * count + m] * z[n] * z[m]
        for n in range(count)
        for m in range(count)
    )
    variables = (*z, *zb)
    derivatives = (
        *(-sp.I * u[n] * z[n] for n in range(count)),
        *(sp.I * u[n] * zb[n] for n in range(count)),
    )
    direct = sum(
        sp.diff(expression, variable) * derivative
        for variable, derivative in zip(variables, derivatives)
    )
    expected = sum(
        -sp.I * u[n] * linear[n] * z[n] for n in range(count)
    )
    expected += sp.Rational(1, 2) * sum(
        -sp.I
        * (u[n] - u[m])
        * hermitian[n * count + m]
        * z[n]
        * zb[m]
        - sp.I
        * (u[n] + u[m])
        * transpose[n * count + m]
        * z[n]
        * z[m]
        for n in range(count)
        for m in range(count)
    )
    require_zero(direct - expected, "full-current frequency flow failed")

    return {
        "auxiliary_carriers": (
            "For T_0-epsilon<=xi<=T_0 put "
            "z_n(xi)=omega_a A_n exp(-i*xi*u_n). Then "
            "z_(n,xi)=-i*u_n*z_n and "
            "conj(z_n)_(xi)=i*u_n*conj(z_n)."
        ),
        "current_family": (
            "Define Psi(xi)=Re sum_(n<=B)L_nz_n(xi)+"
            "(1/2)Re sum_(n,m<=B)[H^+_(n,m)z_n(xi)conj(z_m(xi))"
            "+T^+_(n,m)z_n(xi)z_m(xi)]. At xi=Omega, "
            "Psi(Omega)=h^2*Phi_B exactly."
        ),
        "derivative": (
            "Psi'(xi)=Re sum_n[-i*u_n*L_n]z_n+"
            "(1/2)Re sum_(n,m){[-i(u_n-u_m)H^+_(n,m)]"
            "z_nconj(z_m)+[-i(u_n+u_m)T^+_(n,m)]z_nz_m}."
        ),
        "transfer": (
            "h^2*Phi_B-Psi(T_0)="
            "-epsilon*integral_0^1 Psi'(T_0-theta*epsilon)dtheta."
        ),
        "symmetry": (
            "If H^+_(m,n)=conj(H^+_(n,m)), then "
            "-i(u_n-u_m)H^+_(n,m) is Hermitian. If "
            "T^+_(m,n)=T^+_(n,m), then "
            "-i(u_n+u_m)T^+_(n,m) is transpose-symmetric."
        ),
        "interpretation": (
            "The complete saddle transfer remains one endpoint-linear, one "
            "phase-difference, and one phase-sum family. No new correlation "
            "family or four-variable expansion is created."
        ),
        "symbolic_difference": str(sp.simplify(direct - expected)),
        "finite_carriers": count,
    }


def leading_pair_flow_certificate() -> dict:
    count = 3
    distances = sp.symbols(f"u0:{count}", real=True)
    x = sp.symbols(f"x0:{count}", real=True)
    y = sp.symbols(f"y0:{count}", real=True)
    u_x = sp.symbols("u_x", real=True)

    f = sum(x)
    g = sum(y)
    p_1 = sum(distances[n] * x[n] for n in range(count))
    q_1 = sum(distances[n] * y[n] for n in range(count))
    p_2 = sum(distances[n] ** 2 * x[n] for n in range(count))
    q_3 = sum(distances[n] ** 3 * y[n] for n in range(count))
    k_direct = sp.expand(q_1 * p_2 - f * q_3 + 2 * u_x * (q_1 * g - f * p_1))

    k_pair = sp.Integer(0)
    diagonal = sp.Integer(0)
    for n in range(count):
        diagonal += -2 * u_x * distances[n] * (x[n] ** 2 - y[n] ** 2)
        for m in range(count):
            delta = distances[n] - distances[m]
            sigma = distances[n] + distances[m]
            im_sum = x[n] * y[m] + y[n] * x[m]
            im_difference = y[n] * x[m] - x[n] * y[m]
            re_sum = x[n] * x[m] - y[n] * y[m]
            k_pair += -sp.Rational(1, 4) * (
                delta**2 * sigma * im_sum
                + delta * sigma**2 * im_difference
            ) - u_x * sigma * re_sum
    require_zero(k_direct - k_pair, "leading pair-flow identity failed")

    pair_diagonal = sp.expand(
        sum(
            -u_x
            * (2 * distances[n])
            * (x[n] ** 2 - y[n] ** 2)
            for n in range(count)
        )
    )
    require_zero(pair_diagonal - diagonal, "leading diagonal flow failed")

    return {
        "coordinate_integrand": (
            "For M_0=f+ig, M_1=p_1+iq_1, M_2=p_2+iq_2, "
            "and q_3=Im(M_3), put "
            "K=q_1p_2-fq_3+2u_x(q_1g-fp_1)=4P_bulk^(0)'(xi)."
        ),
        "pair_identity": (
            "With Delta_nm=u_n-u_m and Sigma_nm=u_n+u_m, "
            "K=-(1/4)sum_(n,m)[Delta_nm^2 Sigma_nm Im(z_nz_m)"
            "+Delta_nm Sigma_nm^2 Im(z_nconj(z_m))]"
            "-u_x sum_(n,m)Sigma_nm Re(z_nz_m)."
        ),
        "diagonal": (
            "The Hermitian and cubic transpose pieces vanish on n=m. "
            "The exact diagonal is "
            "K_diag=-2u_x sum_n u_n Re(z_n^2), which has no fixed sign."
        ),
        "kernel_origin": (
            "This pair identity is exactly the xi derivative of the ideal "
            "kernels H_0^+=-Sigma_nm^2/8 and "
            "T_0^+=-Delta_nm^2/8-i*u_x/2."
        ),
        "scope": (
            "The identity preserves both physical phase families. Its "
            "diagonal coefficient is small, but no diagonal or aggregate "
            "sign is inferred without the actual phases."
        ),
        "symbolic_difference": str(sp.simplify(k_direct - k_pair)),
        "finite_carriers": count,
    }


def degree_certificate() -> dict:
    ell_n, ell_m, log_a = sp.symbols("ell_n ell_m log_a", real=True)
    c0, c1, c2 = sp.symbols("c0 c1 c2")
    d0, d1, d2 = sp.symbols("d0 d1 d2")
    r0, r1 = sp.symbols("r0 r1")
    q0, q1 = sp.symbols("q0 q1")
    x0, x1 = sp.symbols("x0 x1")
    x_t, a_tx, a_t, x_tx = sp.symbols("x_t a_tx a_t x_tx")

    def c(value: sp.Expr) -> sp.Expr:
        return c0 + c1 * value + c2 * value**2

    def d(value: sp.Expr) -> sp.Expr:
        return d0 + d1 * value + d2 * value**2

    def r(value: sp.Expr) -> sp.Expr:
        return r0 + r1 * value

    def q(value: sp.Expr) -> sp.Expr:
        return (q0 + q1 * value) * c(value) + d(value)

    def r_x(value: sp.Expr) -> sp.Expr:
        return x0 + x1 * value

    def n_poly(value: sp.Expr) -> sp.Expr:
        return sp.expand(r(value) * q(value) + r_x(value) * c(value))

    linear = sp.expand(
        x_t * n_poly(ell_n)
        + a_tx * c(ell_n)
        - a_t * q(ell_n)
        - x_tx * r(ell_n) * c(ell_n)
    )

    def pair(first: sp.Expr, second: sp.Expr) -> sp.Expr:
        return sp.expand(
            c(first) * n_poly(second)
            - r(first) * c(first) * q(second)
        )

    h_plus = sp.expand(
        (pair(ell_n, ell_m) + pair(ell_m, ell_n)) / 2
    )
    t_plus = h_plus
    flow_linear = sp.expand((log_a - ell_n) * linear)
    flow_h = sp.expand((ell_m - ell_n) * h_plus)
    flow_t = sp.expand((2 * log_a - ell_n - ell_m) * t_plus)

    degrees = {
        "L_n": int(sp.degree(linear, ell_n)),
        "u_n_L_n": int(sp.degree(flow_linear, ell_n)),
        "H_plus_ell_n": int(sp.degree(h_plus, ell_n)),
        "H_plus_ell_m": int(sp.degree(h_plus, ell_m)),
        "flow_H_ell_n": int(sp.degree(flow_h, ell_n)),
        "flow_H_ell_m": int(sp.degree(flow_h, ell_m)),
        "flow_T_ell_n": int(sp.degree(flow_t, ell_n)),
        "flow_T_ell_m": int(sp.degree(flow_t, ell_m)),
    }
    if max(degrees.values()) != 5:
        raise RuntimeError("saddle-flow degree audit failed")
    if any(value > 5 for value in degrees.values()):
        raise RuntimeError("saddle-flow degree exceeds five")

    return {
        "degree_statement": (
            "In ell_n=log(n), C,D have degree at most 2, R,R_x degree "
            "at most 1, Q degree at most 3, N and L degree at most 4. "
            "Multiplication by u_n, u_n-u_m, or u_n+u_m raises each "
            "one-variable degree to at most 5."
        ),
        "six_moment_closure": (
            "Therefore Psi'(xi) closes exactly on the six correction-free "
            "moments G_0(xi),...,G_5(xi), with no seventh moment and no new "
            "correction denominator."
        ),
        "degrees": degrees,
        "max_degree": max(degrees.values()),
        "moment_count": 6,
    }


def order_five_mangoldt_certificate() -> dict:
    log_d, log_m, t_symbol, s_symbol = sp.symbols(
        "log_d log_m t_symbol s_symbol"
    )
    product_exponent = (
        t_symbol * log_d**2 / 4
        - s_symbol * log_d
        + t_symbol * log_m**2 / 4
        - s_symbol * log_m
        + t_symbol * log_d * log_m / 2
    )
    combined_exponent = (
        t_symbol * (log_d + log_m) ** 2 / 4
        - s_symbol * (log_d + log_m)
    )
    require_zero(
        combined_exponent - product_exponent,
        "auxiliary heat-factor exponent failed",
    )

    checked_cutoffs = (17, 31, 64, 97)
    rows = 0
    symmetric_rows = 0
    balanced_rows = 0
    for cutoff in checked_cutoffs:
        for integer in range(2, cutoff + 1):
            factors = sp.factorint(integer)
            primes = tuple(sorted(factors))
            symbols = {prime: sp.Symbol(f"L_{prime}") for prime in primes}
            log_n = sum(factors[prime] * symbols[prime] for prime in primes)

            def mangoldt(divisor: int) -> sp.Expr:
                factorization = sp.factorint(divisor)
                if len(factorization) != 1:
                    return sp.Integer(0)
                prime = next(iter(factorization))
                return symbols.get(prime, sp.Symbol(f"L_{prime}"))

            divisors = sp.divisors(integer)
            one_sided = sum(mangoldt(divisor) for divisor in divisors)
            if sp.expand(log_n**4 * one_sided - log_n**5) != 0:
                raise RuntimeError(f"order-five Mangoldt identity failed at {integer}")
            rows += 1
            symmetric = sum(
                (mangoldt(divisor) + mangoldt(integer // divisor)) / 2
                for divisor in divisors
            )
            if sp.expand(log_n**4 * symmetric - log_n**5) != 0:
                raise RuntimeError(
                    f"symmetric order-five Mangoldt identity failed at {integer}"
                )
            symmetric_rows += 1

        root = int(cutoff**0.5)
        all_pairs = {
            (d_value, m_value)
            for d_value in range(1, cutoff + 1)
            for m_value in range(1, cutoff // d_value + 1)
        }
        square = {
            pair for pair in all_pairs if pair[0] <= root and pair[1] <= root
        }
        wings = {
            pair
            for pair in all_pairs
            if (pair[0] <= root < pair[1])
            or (pair[1] <= root < pair[0])
        }
        if square | wings != all_pairs or square & wings:
            raise RuntimeError(f"order-five balanced split failed at {cutoff}")
        balanced_rows += 1

    return {
        "phase_factor": (
            "Put c_xi=exp[i(Omega-xi)log(a)], s_xi=sigma-i*xi, and "
            "q_n^(0)(xi)=exp[t log(n)^2/4-s_xi log(n)]. Then "
            "z_n(xi)=(eta/S_a)c_xi q_n^(0)(xi). The factor c_xi is "
            "common to every carrier and stays outside each divisor sum."
        ),
        "order_five": (
            "G_5(xi)=(eta/S_a)c_xi sum_(dm<=B)Lambda(d)log(dm)^4"
            "q_(dm)^(0)(xi)=(eta*c_xi/(2S_a))sum_(dm<=B)"
            "[Lambda(d)+Lambda(m)]log(dm)^4q_(dm)^(0)(xi)."
        ),
        "heat_factor": (
            "The exact factor q_(dm)^(0)=q_d^(0)q_m^(0)"
            "exp[(t/2)log(d)log(m)] is unchanged; xi only changes the "
            "common logarithmic phase."
        ),
        "balanced_split": (
            "For D=floor(sqrt(B)), the order-five sum has the same exact "
            "D-by-D square and two hyperbolic wings as orders one through four."
        ),
        "prime_edge_guard": (
            "For every prime sqrt(B)<p<=B, the order-five symmetric "
            "{1,p} minor is [[0,log(p)^5/2],[log(p)^5/2,0]] and has "
            "determinant -log(p)^10/4<0."
        ),
        "audit": {
            "cutoffs": list(checked_cutoffs),
            "order_five_rows": rows,
            "symmetric_order_five_rows": symmetric_rows,
            "balanced_hyperbola_rows": balanced_rows,
            "max_moment_order": 5,
            "heat_factor_exponent_difference": str(
                sp.expand(combined_exponent - product_exponent)
            ),
        },
    }


def exact_payload() -> dict:
    return {
        "composed_target": (
            "Prove directly "
            "Psi(T_0)-epsilon integral_0^1"
            "Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably "
            "h^2/800, for the actual endpoint and six-moment Xi vector."
        ),
        "analytic_handoff": (
            "Apply a uniform signed Type-I/II, Vaughan, or reciprocal-Poisson "
            "estimate to the endpoint-linear, Hermitian flow, and transpose "
            "flow families together. The new arithmetic input is only the "
            "order-five logarithmic moment; no new correlation geometry is "
            "needed."
        ),
        "intermediate_boundary": (
            "Only xi=Omega is the physical carrier state. Intermediate xi "
            "and xi=T_0 form an exact auxiliary homotopy and are not asserted "
            "to satisfy the terminal source recurrence separately."
        ),
        "proof_boundary": (
            "This proves the exact endpoint-composed frequency-flow and "
            "transfer identities, preservation of Hermitian and transpose "
            "symmetry, the leading pair-flow formula, degree-five closure, "
            "and the order-five Mangoldt/balanced-hyperbola extension. It "
            "does not prove a signed flow estimate, saddle-proxy upper bound, "
            "Phi_B upper bound, contact exclusion, retained aggregate or "
            "Xi-level current theorem, Q209, Lambda<=0, PF-infinity, RH, or "
            "a prize-level result."
        ),
    }


def build_rows(
    flow: dict,
    leading: dict,
    degree: dict,
    mangoldt: dict,
    exact: dict,
) -> list[SaddleFlowRow]:
    return [
        SaddleFlowRow("sfr_01_domain", "auxiliary saddle segment", "proved", "The flow stays on the certified q=1 frequency segment.", flow["auxiliary_carriers"], exact["intermediate_boundary"]),
        SaddleFlowRow("sfr_02_family", "endpoint-composed current family", "proved", "The full joined current extends to one exact auxiliary family.", flow["current_family"], "All physical coefficients are held fixed along xi."),
        SaddleFlowRow("sfr_03_derivative", "full frequency derivative", "proved", "Termwise differentiation preserves the three-family decomposition.", flow["derivative"], flow["interpretation"], {"difference": flow["symbolic_difference"], "finite_carriers": flow["finite_carriers"]}),
        SaddleFlowRow("sfr_04_symmetry", "flow-kernel symmetries", "proved", "The differentiated pair kernels retain Hermitian and transpose symmetry.", flow["symmetry"], "Symmetry is an organization theorem, not a sign theorem."),
        SaddleFlowRow("sfr_05_transfer", "full signed saddle transfer", "proved", "The actual Phi_B differs from its saddle proxy by one signed integral.", flow["transfer"], "No positive-mass transfer bound is applied."),
        SaddleFlowRow("sfr_06_leading", "leading coordinate flow", "proved", "The Section 11.163 signed integrand is the derivative of the ideal bulk current.", leading["coordinate_integrand"], "Moments through order three suffice only for the ideal part."),
        SaddleFlowRow("sfr_07_pair", "leading pair-flow identity", "proved", "The leading signed integrand splits into the same two phase families.", leading["pair_identity"] + " " + leading["kernel_origin"], "No four-variable sum is introduced.", {"difference": leading["symbolic_difference"], "finite_carriers": leading["finite_carriers"]}),
        SaddleFlowRow("sfr_08_diagonal", "flow diagonal", "proved", "Only the shear transpose term survives on the pair diagonal.", leading["diagonal"], leading["scope"]),
        SaddleFlowRow("sfr_09_degree", "degree-five coefficient audit", "proved", "Frequency differentiation raises every one-variable logarithmic degree by at most one.", degree["degree_statement"], "The audit uses the exact C,D,R,Q,R_x,N,L degree hierarchy.", degree["degrees"]),
        SaddleFlowRow("sfr_10_six", "six-moment closure", "proved", "The full saddle-flow current closes on G_0 through G_5.", degree["six_moment_closure"], "Sufficiency is proved; no claim of irreducible minimality is made."),
        SaddleFlowRow("sfr_11_mangoldt", "order-five Mangoldt extension", "proved", "The one new moment has the same exact symmetric divisor-pair form.", mangoldt["phase_factor"] + " " + mangoldt["order_five"] + " " + mangoldt["heat_factor"], "The common auxiliary phase stays outside the divisor sum; the form remains complex bilinear, not Hermitian.", mangoldt["audit"]),
        SaddleFlowRow("sfr_12_hyperbola", "order-five balanced split", "proved", "The new moment has the same square-and-wing Type-I/II organization.", mangoldt["balanced_split"], "Organization alone supplies no cancellation.", mangoldt["audit"]),
        SaddleFlowRow("sfr_13_guard", "order-five positivity guard", "guard_validated", "The extended symmetric kernel is still indefinite on every prime edge.", mangoldt["prime_edge_guard"], "Transpose symmetry cannot be promoted to positivity."),
        SaddleFlowRow("sfr_14_target", "composed signed target", "open", "The exact transfer and saddle proxy must now be estimated together.", exact["composed_target"] + " " + exact["analytic_handoff"], "No signed constant is proved."),
        SaddleFlowRow("sfr_15_boundary", "proof boundary", "open", "The six-moment handoff narrows but does not close the arithmetic theorem.", exact["proof_boundary"], "RH remains open."),
    ]


def render_note(payload: dict) -> str:
    flow = payload["full_current_flow_certificate"]
    leading = payload["leading_pair_flow_certificate"]
    degree = payload["degree_certificate"]
    mangoldt = payload["order_five_mangoldt_certificate"]
    exact = payload["exact"]
    return f"""# Endpoint-Composed Six-Moment Saddle Flow

Date: 2026-08-01

Status: exact saddle-flow reduction with `0 signed flow bounds` and `0 Phi_B
bounds`. This is not a proof of RH.

## Auxiliary Current

```text
{flow['auxiliary_carriers']}

{flow['current_family']}
```

{exact['intermediate_boundary']}

## Exact Flow

```text
{flow['derivative']}

{flow['transfer']}
```

{flow['symmetry']} {flow['interpretation']}

## Leading Pair Kernel

```text
{leading['coordinate_integrand']}

{leading['pair_identity']}

{leading['diagonal']}
```

{leading['kernel_origin']} {leading['scope']}

## Six-Moment Closure

```text
{degree['degree_statement']}

{degree['six_moment_closure']}
```

The exact degree audit is

```json
{json.dumps(degree['degrees'], indent=2)}
```

## Order-Five Mangoldt Form

```text
{mangoldt['phase_factor']}

{mangoldt['order_five']}

{mangoldt['heat_factor']}

{mangoldt['balanced_split']}

{mangoldt['prime_edge_guard']}
```

The formal audit checks `{mangoldt['audit']['order_five_rows']}` one-sided
and `{mangoldt['audit']['symmetric_order_five_rows']}` symmetric rows.

## Exact Target

```text
{exact['composed_target']}
```

{exact['analytic_handoff']}

## Pi Provenance

No new `pi` is introduced. The existing `u_x=h^2/(8pi)` and saddle phase
inherit the completed-zeta, Riemann--Siegel, and Poisson normalizations
recorded in Formal Core Section 11.118.

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
    flow = full_current_flow_certificate()
    leading = leading_pair_flow_certificate()
    degree = degree_certificate()
    mangoldt = order_five_mangoldt_certificate()
    exact = exact_payload()
    rows = build_rows(flow, leading, degree, mangoldt, exact)
    return {
        "kind": STEM,
        "date": "2026-08-01",
        "status": (
            "exact endpoint-composed six-moment saddle-flow reduction and "
            "order-five Mangoldt handoff; signed flow estimate open; no "
            "Phi_B bound, Lambda<=0, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "full_current_flow_certificate": flow,
        "leading_pair_flow_certificate": leading,
        "degree_certificate": degree,
        "order_five_mangoldt_certificate": mangoldt,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "full_current_flow_identities": 1,
            "signed_transfer_integrals": 1,
            "flow_kernel_symmetries": 2,
            "leading_pair_flow_identities": 1,
            "correction_free_moments": 6,
            "max_moment_order": 5,
            "new_mangoldt_orders": 1,
            "order_five_prime_edge_guards": 1,
            "signed_flow_bounds": 0,
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
        "built six-moment saddle-flow reduction: 15 rows, "
        "1 full-current flow, 2 flow symmetries, 1 leading pair identity, "
        "6 moments, 1 new Mangoldt order, 0 signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
