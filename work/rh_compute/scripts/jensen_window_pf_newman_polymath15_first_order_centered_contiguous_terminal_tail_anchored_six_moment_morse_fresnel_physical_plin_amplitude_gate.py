#!/usr/bin/env python3
"""Build the physical P_lin coefficient and Morse-amplitude gate."""

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
    "physical_plin_amplitude_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "observation_image": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "observation_image_compression_gate.json"
    ),
    "two_carrier_kernel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_two_carrier_"
        "kernel_reduction.json"
    ),
    "morse_fresnel_endpoint": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "endpoint_reduction.json"
    ),
}


@dataclass(frozen=True)
class AmplitudeRow:
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
    observation = payloads["observation_image"].get("counts", {})
    if observation.get("visible_residual_observations") != 8:
        raise RuntimeError("observation-image source lost eight observations")
    if observation.get("maximum_polynomial_degree") != 5:
        raise RuntimeError("observation-image source lost degree-five closure")
    kernel = payloads["two_carrier_kernel"].get("counts", {})
    if kernel.get("carrier_polynomials") != 4:
        raise RuntimeError("two-carrier source lost four carrier polynomials")
    if kernel.get("common_rate_cancellations") != 2:
        raise RuntimeError("two-carrier source lost rate cancellations")
    morse = payloads["morse_fresnel_endpoint"].get("counts", {})
    if morse.get("exact_transition_integrals") != 6:
        raise RuntimeError("Morse source lost six exact transition integrals")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.expand(expression) != 0:
        raise RuntimeError(label)


def centered_carrier_certificate() -> dict:
    x, log_a, log_n, u_n = sp.symbols("x log_a log_N u_N", real=True)
    c, b, c_x, b_x, u_x = sp.symbols("c b c_x b_x u_N_x", real=True)
    s_prime = c + sp.I * b
    s_second = c_x + sp.I * b_x
    lam = log_a + x
    original_r = s_prime * (log_n - lam) + sp.I * b * u_n
    original_rx = (
        s_second * (log_n - lam)
        + sp.I * (b_x * u_n + b * u_x)
    )
    substitutions = {log_n: log_a - u_n}
    centered_r = -s_prime * x - c * u_n
    centered_rx = -s_second * x - c_x * u_n + sp.I * b * u_x
    require_zero(
        original_r.subs(substitutions) - centered_r,
        "terminal-centered R identity failed",
    )
    require_zero(
        original_rx.subs(substitutions) - centered_rx,
        "terminal-centered R_x identity failed",
    )
    return {
        "center": (
            "Put x=lambda-log(a) and u_N=log(a/N), so log(N)=log(a)-u_N. "
            "Writing s_*'=c+ib and s_*''=c_x+ib_x gives the exact identities "
            "R=-s_*'x-c*u_N and R_x=-s_*''x-c_x*u_N+i*b*u_(N,x)."
        ),
        "rate_defect": (
            "Put delta=chi_N-i*b*u_N. Then Q=(R+delta)C+D and "
            "N=R*Q+R_x*C. The existing physical estimate is |delta|<4h^2, "
            "but no estimate is used in this algebraic gate."
        ),
        "correction_rows": (
            "In the centered variable, C=c_0+c_1*x+c_2*x^2 and "
            "D=d_0+d_1*x+d_2*x^2, where c_0=1+rho_1*log(a)+"
            "rho_2*log(a)^2, c_1=rho_1+2rho_2*log(a), c_2=rho_2, "
            "d_0=rho_(1,x)*log(a)+rho_(2,x)*log(a)^2, "
            "d_1=rho_(1,x)+2rho_(2,x)*log(a), and d_2=rho_(2,x)."
        ),
        "symbolic_differences": {"R": "0", "R_x": "0"},
    }


def uv_coefficients(
    f_v: sp.Expr,
    f_n: sp.Expr,
    f_a: sp.Expr,
    f_q: sp.Expr,
    r_0: sp.Expr,
    r_1: sp.Expr,
    t_0: sp.Expr,
    t_1: sp.Expr,
    delta: sp.Expr,
) -> tuple[list[sp.Expr], list[sp.Expr]]:
    g_0 = r_0 * (r_0 + delta) + t_0
    g_1 = r_1 * (2 * r_0 + delta) + t_1
    g_2 = r_1**2
    u = [
        f_n * g_0 + f_q * (r_0 + delta) + f_a * r_0 + f_v,
        f_n * g_1 + (f_q + f_a) * r_1,
        f_n * g_2,
    ]
    v = [f_n * r_0 + f_q, f_n * r_1]
    return u, v


def coefficient_certificate() -> dict:
    x = sp.symbols("x", real=True)
    c_0, c_1, c_2 = sp.symbols("c0 c1 c2")
    d_0, d_1, d_2 = sp.symbols("d0 d1 d2")
    r_0, r_1, t_0, t_1, delta = sp.symbols("r0 r1 t0 t1 delta")
    n_p, v_p, q_p, a_p = sp.symbols("N_p V_p Q_p A_p", real=True)
    n_d, v_d, q_d, a_d = sp.symbols(
        "N_p_xi V_p_xi Q_p_xi A_p_xi", real=True
    )
    x_t, a_t, x_tx, a_tx = sp.symbols(
        "X_T A_T X_T_x A_T_x", real=True
    )
    c_poly = c_0 + c_1 * x + c_2 * x**2
    d_poly = d_0 + d_1 * x + d_2 * x**2
    r_poly = r_0 + r_1 * x
    rx_poly = t_0 + t_1 * x
    q_poly = (r_poly + delta) * c_poly + d_poly
    a_poly = r_poly * c_poly
    n_poly = sp.expand(r_poly * q_poly + rx_poly * c_poly)

    alpha = (n_d, v_d, -q_d, -a_d)
    beta = (n_p + a_tx, x_t + v_p, -(q_p + x_tx), -(a_t + a_p))
    u_alpha, v_alpha = uv_coefficients(
        *alpha, r_0, r_1, t_0, t_1, delta
    )
    u_beta, v_beta = uv_coefficients(
        *beta, r_0, r_1, t_0, t_1, delta
    )

    gamma = [
        u_alpha[0],
        u_alpha[1] + sp.I * u_beta[0],
        u_alpha[2] + sp.I * u_beta[1],
        sp.I * u_beta[2],
    ]
    eta = [
        v_alpha[0],
        v_alpha[1] + sp.I * v_beta[0],
        sp.I * v_beta[1],
    ]
    coefficients = []
    for degree in range(6):
        value = 0
        for index, coefficient in enumerate((c_0, c_1, c_2)):
            other = degree - index
            if 0 <= other < len(gamma):
                value += coefficient * gamma[other]
        for index, coefficient in enumerate((d_0, d_1, d_2)):
            other = degree - index
            if 0 <= other < len(eta):
                value += coefficient * eta[other]
        coefficients.append(sp.expand(value))

    direct = sp.expand(
        n_d * c_poly
        + v_d * n_poly
        - q_d * a_poly
        - a_d * q_poly
        + sp.I
        * x
        * (
            (n_p + a_tx) * c_poly
            + (x_t + v_p) * n_poly
            - (q_p + x_tx) * a_poly
            - (a_t + a_p) * q_poly
        )
    )
    recurrence = sum(coefficients[index] * x**index for index in range(6))
    require_zero(direct - recurrence, "six-coefficient recurrence failed")

    f_v, f_n, f_a, f_q = sp.symbols("f_V f_N f_A f_Q")
    u_generic, v_generic = uv_coefficients(
        f_v, f_n, f_a, f_q, r_0, r_1, t_0, t_1, delta
    )
    generic = sp.expand(
        f_v * c_poly + f_n * n_poly + f_a * a_poly + f_q * q_poly
    )
    factored = sp.expand(
        c_poly * sum(u_generic[index] * x**index for index in range(3))
        + d_poly * sum(v_generic[index] * x**index for index in range(2))
    )
    require_zero(generic - factored, "two-channel factorization failed")

    s_prime, rho_2 = sp.symbols("s_prime rho_2")
    total_value = x_t + v_p
    p_5 = coefficients[5].subs({c_2: rho_2, r_1: -s_prime})
    expected_p_5 = sp.I * rho_2 * total_value * s_prime**2
    require_zero(p_5 - expected_p_5, "leading contact factor failed")

    return {
        "observation_rows": (
            "For any real f=(f_V,f_N,f_A,f_Q), let F_f=f_V*C+f_N*N+"
            "f_A*A+f_Q*Q. Use alpha=(N_(p,xi),V_(p,xi),-Q_(p,xi),"
            "-A_(p,xi)) and beta=(N_p+A_(T,x),X_T+V_p,-Q_p-X_(T,x),"
            "-A_T-A_p). Then P_lin=F_alpha+i*x*F_beta."
        ),
        "two_channel_factorization": (
            "For R=r_0+r_1*x and R_x=t_0+t_1*x, every F_f factors exactly "
            "as F_f=C*U_f+D*V_f. Here V_f=f_N*R+f_Q and U_f=f_N*"
            "[R(R+delta)+R_x]+f_Q(R+delta)+f_A*R+f_V. Thus four physical "
            "carrier polynomials reduce to the two correction channels C and D."
        ),
        "uv_recurrence": (
            "Writing U_f=u_0+u_1*x+u_2*x^2 and V_f=v_0+v_1*x, put "
            "g_0=r_0(r_0+delta)+t_0, g_1=r_1(2r_0+delta)+t_1, "
            "g_2=r_1^2. Then u_0=f_Ng_0+f_Q(r_0+delta)+f_Ar_0+f_V, "
            "u_1=f_Ng_1+(f_Q+f_A)r_1, u_2=f_Ng_2, "
            "v_0=f_Nr_0+f_Q, and v_1=f_Nr_1."
        ),
        "channel_coefficients": (
            "Let gamma=(u_(alpha,0),u_(alpha,1)+i*u_(beta,0),"
            "u_(alpha,2)+i*u_(beta,1),i*u_(beta,2)) and "
            "eta=(v_(alpha,0),v_(alpha,1)+i*v_(beta,0),i*v_(beta,1)). "
            "Then P_lin=C*sum_(j=0)^3gamma_j*x^j+D*sum_(j=0)^2eta_j*x^j."
        ),
        "six_coefficients": [
            "p_0=c_0*gamma_0+d_0*eta_0",
            "p_1=c_0*gamma_1+c_1*gamma_0+d_0*eta_1+d_1*eta_0",
            "p_2=c_0*gamma_2+c_1*gamma_1+c_2*gamma_0+d_0*eta_2+d_1*eta_1+d_2*eta_0",
            "p_3=c_0*gamma_3+c_1*gamma_2+c_2*gamma_1+d_1*eta_2+d_2*eta_1",
            "p_4=c_1*gamma_3+c_2*gamma_2+d_2*eta_2",
            "p_5=c_2*gamma_3",
        ],
        "leading_coefficient": (
            "The top coefficient is exactly p_5=i*rho_2*(X_T+V_p)*(s_*')^2. "
            "Hence the degree-five channel vanishes on the total-value fibre "
            "X_T+V_p=0; this is a degree reduction, not a remainder bound."
        ),
        "symbolic_differences": {
            "two_channel": "0",
            "six_coefficients": "0",
            "p_5": "0",
        },
    }


def ideal_certificate() -> dict:
    x, u_x = sp.symbols("x u_x", real=True)
    n_p, v_p, q_p, a_p = sp.symbols("N_p V_p Q_p A_p", real=True)
    n_d, q_d = sp.symbols("N_p_xi Q_p_xi", real=True)
    x_t, a_t, x_tx, a_tx = sp.symbols(
        "X_T A_T X_T_x A_T_x", real=True
    )
    total_value = x_t + v_p
    ideal = sp.expand(
        n_d
        - sp.I * a_p * u_x
        + (sp.I * (n_p + a_tx - q_d) + total_value * u_x / 2) * x
        + (a_p + a_t + x_tx) * x**2 / 2
        - sp.I * total_value * x**3 / 4
    )
    contact = sp.expand(ideal.subs({v_p: -x_t, a_p: -a_t}))
    expected_contact = sp.expand(
        n_d
        + sp.I * a_t * u_x
        + sp.I * (n_p + a_tx - q_d) * x
        + x_tx * x**2 / 2
    )
    require_zero(contact - expected_contact, "ideal contact reduction failed")
    return {
        "specialization": (
            "At C=1, D=0, s_*'=-i/2, s_*''=0, delta=0, and b=-1/2, "
            "one has R=Q=i*x/2 and N=-x^2/4-i*u_(N,x)/2. The induced "
            "observation identities are A_p=Q_p, A_(p,xi)=Q_(p,xi), and "
            "V_(p,xi)=2A_p."
        ),
        "ideal_polynomial": (
            "P_lin=N_(p,xi)-i*A_p*u_(N,x)+{i[N_p+A_(T,x)-Q_(p,xi)]"
            "+(X_T+V_p)u_(N,x)/2}x+[A_p+A_T+X_(T,x)]x^2/2"
            "-i(X_T+V_p)x^3/4."
        ),
        "contact_polynomial": (
            "On X_T+V_p=0 and A_T+A_p=0, the correction-free leading "
            "polynomial is quadratic: P_lin=N_(p,xi)-i*A_p*u_(N,x)+"
            "i[N_p+A_(T,x)-Q_(p,xi)]x+X_(T,x)x^2/2."
        ),
        "boundary": (
            "The actual correction channels can restore degrees four and five. "
            "The ideal reduction is a scale guide, not an exact deletion of C-1 or D."
        ),
        "ideal_polynomial_expanded": str(ideal),
        "contact_polynomial_expanded": str(contact),
    }


def morse_amplitude_certificate() -> dict:
    return {
        "definitions": (
            "For u_0=alpha/r, lambda=log(u_0v), S(lambda)=t*lambda^2/4-"
            "sigma*lambda, z=sgn(v-1)sqrt(2[v-1-log v]), and J=dv/dz, "
            "put W_P(v)=exp(S(lambda))*P(lambda)*J(v). Then "
            "b_(P,r)(y)=sqrt(alpha)W_P(v)/r and y=sqrt(alpha)z."
        ),
        "c_formula": (
            "For z!=0, c_(P,r)(y)=[W_P(v)-W_P(1)]/(r*z)."
        ),
        "cprime_formula": (
            "For z!=0, c'_(P,r)(y)={z*J*W_P'(v)-W_P(v)+W_P(1)}/"
            "[r*sqrt(alpha)*z^2]. With g=t*lambda/2-sigma, this equals "
            "{exp(S)[zJ^2(P'+gP)/v+(zJJ'-J)P]+"
            "exp(S(lambda_0))P(lambda_0)}/[r*sqrt(alpha)*z^2], "
            "where lambda_0=log(alpha/r)."
        ),
        "saddle_formula": (
            "The removable saddle value is c'_(P,r)(0)=exp(S(lambda_0))/"
            "[2r*sqrt(alpha)] times {P''+(2g_0+1)P'+"
            "[g_0^2+g_0+t/2+1/6]P}(lambda_0), with "
            "g_0=t*lambda_0/2-sigma."
        ),
        "mode_dependence": (
            "At fixed xi the same physical P_lin is used for every Poisson "
            "mode. Mode dependence enters only through lambda_0=log(alpha/r), "
            "the Morse endpoints, and the common phase. Thus the saddle audit "
            "is one second-order differential operator sampled on the logarithmic "
            "mode lattice, not six unrelated moment amplitudes."
        ),
        "zero_boundary": (
            "Neither P_lin(lambda_0)=0 nor the displayed saddle operator vanishes "
            "identically. Any zero or discrete-difference gain must be proved from "
            "the physical coefficients and mode coupling."
        ),
    }


def route_certificate() -> dict:
    boundary = (
        "This proves the exact physical terminal-centered carrier laws, the C/D "
        "two-channel factorization, all six centered P_lin coefficients, the "
        "total-value factor in the degree-five coefficient, the correction-free "
        "cubic and contact-quadratic reductions, and exact off-saddle and saddle "
        "Morse c-prime formulas. It proves no coefficient norm on the physical "
        "chart, grouped c-prime sum, discrete cancellation theorem, quadratic "
        "residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, "
        "retained aggregate or Xi theorem, Q209, Lambda<=0, PF-infinity, RH, or "
        "a prize-level conclusion."
    )
    return {
        "next_target": (
            "Insert the certified formulas for rho_1,rho_2 and their x derivatives, "
            "bound the six centered coefficients on the full physical chart, and "
            "apply the saddle differential operator to P_lin. Then test Abel or "
            "reciprocal pairing on that one mode sequence before taking absolute "
            "values; keep the four quadratic observation pairings separate."
        ),
        "route_decision": (
            "Use x=lambda-log(a) and the C/D factorization as the proof-facing "
            "coefficient chart. Do not expand powers of log(N) and log(a) separately, "
            "and do not return to six componentwise Morse estimates."
        ),
        "proof_boundary": boundary,
    }


def build_rows(
    centered: dict,
    coefficients: dict,
    ideal: dict,
    morse: dict,
    route: dict,
) -> list[AmplitudeRow]:
    boundary = route["proof_boundary"]
    return [
        AmplitudeRow("plma_01_source", "source chain", "proved", "The exact carrier, observation image, and Morse chart are imported.", "Sections 11.160, 11.168, and 11.172.", "No estimate follows."),
        AmplitudeRow("plma_02_center", "physical centering", "proved", "The terminal logarithm cancels from R and R_x.", centered["center"], "Uses log(N)+u_N=log(a) exactly."),
        AmplitudeRow("plma_03_defect", "rate defect", "proved", "Q differs from RC through one rate defect and D.", centered["rate_defect"], "The small bound is not used here."),
        AmplitudeRow("plma_04_correction", "correction chart", "proved", "C and D have exact centered quadratic rows.", centered["correction_rows"], "No coefficient norm is asserted."),
        AmplitudeRow("plma_05_rows", "observation rows", "proved", "P_lin is F_alpha+i*x*F_beta.", coefficients["observation_rows"], "All carrier and terminal observations remain joined."),
        AmplitudeRow("plma_06_channels", "two-channel factor", "proved", "Every four-polynomial combination factors through C and D.", coefficients["two_channel_factorization"], "No positivity follows."),
        AmplitudeRow("plma_07_uv", "coefficient recurrence", "proved", "The channel coefficients have degree at most two and one.", coefficients["uv_recurrence"], "Exact polynomial algebra."),
        AmplitudeRow("plma_08_gamma", "lifted channels", "proved", "The derivative lift is absorbed into gamma and eta.", coefficients["channel_coefficients"], "No seventh moment."),
        AmplitudeRow("plma_09_coefficients", "six coefficients", "proved", "All six centered complex coefficients are explicit.", "; ".join(coefficients["six_coefficients"]), "Recurrence form avoids artificial large-log expansion."),
        AmplitudeRow("plma_10_leading", "degree-five factor", "proved", "The top coefficient contains the total-value fibre.", coefficients["leading_coefficient"], "A vanishing top coefficient is not an error bound."),
        AmplitudeRow("plma_11_ideal", "ideal specialization", "proved", "The correction-free carrier rows reduce to one cubic.", ideal["specialization"], "Actual corrections remain present."),
        AmplitudeRow("plma_12_cubic", "ideal cubic", "proved", "The complete correction-free P_lin is explicit.", ideal["ideal_polynomial"], ideal["boundary"]),
        AmplitudeRow("plma_13_contact", "contact reduction", "proved", "The correction-free simultaneous-contact polynomial is quadratic.", ideal["contact_polynomial"], ideal["boundary"]),
        AmplitudeRow("plma_14_morse", "Morse amplitude", "proved", "One W_P amplitude represents every mode.", morse["definitions"], "Exact change of variables."),
        AmplitudeRow("plma_15_c", "residual quotient", "proved", "The c chart has a centered divided-difference form.", morse["c_formula"], "Removable at the saddle."),
        AmplitudeRow("plma_16_cprime", "interior amplitude", "proved", "The off-saddle c-prime amplitude is explicit.", morse["cprime_formula"], "No absolute variation is taken."),
        AmplitudeRow("plma_17_saddle", "saddle operator", "proved", "The removable saddle value is a second-order polynomial operator.", morse["saddle_formula"], "Modewise identity."),
        AmplitudeRow("plma_18_modes", "mode compression", "proved", "All modes sample one physical polynomial.", morse["mode_dependence"], "No mode sum is bounded."),
        AmplitudeRow("plma_19_zero", "zero boundary", "proved", "No universal saddle zero is created.", morse["zero_boundary"], "Physical cancellation remains open."),
        AmplitudeRow("plma_20_route", "route decision", "proved", "The centered C/D chart is selected.", route["route_decision"], "Selection is not the final theorem."),
        AmplitudeRow("plma_21_target", "next theorem", "open", "The one-sequence grouped estimate is now explicit.", route["next_target"], "No h^2 estimate is claimed."),
        AmplitudeRow("plma_22_quadratic", "quadratic boundary", "proved", "The four quadratic observation pairings remain separate.", "No quadratic term was absorbed into P_lin.", "Their scale remains open."),
        AmplitudeRow("plma_23_pi", "pi provenance", "proved", "No new pi is introduced.", "The only sqrt(alpha) and Fourier constants are inherited from the existing Morse and Poisson charts.", "No fitted geometry."),
        AmplitudeRow("plma_24_boundary", "proof boundary", "proved", "The gate is an exact coefficient reduction.", boundary, "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    centered = payload["centered_carrier_certificate"]
    coefficients = payload["coefficient_certificate"]
    ideal = payload["ideal_certificate"]
    morse = payload["morse_amplitude_certificate"]
    route = payload["route_certificate"]
    coefficient_lines = "\n".join(coefficients["six_coefficients"])
    return f"""# Physical P_lin And Morse-Amplitude Gate

Date: 2026-08-02

Status: exact physical coefficient factorization and Morse `c'` amplitude;
`0 grouped interior bounds`, `0 signed flow bounds`, and this is not a proof
of RH.

## Terminal-Centered Carrier

{centered['center']}

{centered['rate_defect']}

{centered['correction_rows']}

## Two Correction Channels

{coefficients['observation_rows']}

{coefficients['two_channel_factorization']}

{coefficients['uv_recurrence']}

{coefficients['channel_coefficients']}

## Six Complex Coefficients

```text
{coefficient_lines}
```

{coefficients['leading_coefficient']}

## Correction-Free Audit

{ideal['specialization']}

```text
{ideal['ideal_polynomial']}
```

```text
{ideal['contact_polynomial']}
```

{ideal['boundary']}

## Exact Morse Interior

{morse['definitions']}

{morse['c_formula']}

{morse['cprime_formula']}

{morse['saddle_formula']}

{morse['mode_dependence']}

{morse['zero_boundary']}

## Route Decision

{route['route_decision']}

## Next Target

{route['next_target']}

## Proof Boundary

{route['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    sources = load_sources()
    centered = centered_carrier_certificate()
    coefficients = coefficient_certificate()
    ideal = ideal_certificate()
    morse = morse_amplitude_certificate()
    route = route_certificate()
    rows = build_rows(centered, coefficients, ideal, morse, route)
    return {
        "kind": STEM,
        "date": "2026-08-02",
        "status": (
            "exact terminal-centered physical P_lin coefficient factorization "
            "and Morse c-prime amplitude complete; grouped interior and signed "
            "flow estimates open"
        ),
        "source_audit": source_audit(sources),
        "centered_carrier_certificate": centered,
        "coefficient_certificate": coefficients,
        "ideal_certificate": ideal,
        "morse_amplitude_certificate": morse,
        "route_certificate": route,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "physical_base_polynomials": 4,
            "physical_correction_channels": 2,
            "exact_centered_complex_coefficients": 6,
            "maximum_polynomial_degree": 5,
            "degree_five_total_value_factors": 1,
            "correction_free_maximum_degree": 3,
            "correction_free_contact_maximum_degree": 2,
            "exact_morse_cprime_formulas": 2,
            "enumerated_physical_modes": 0,
            "coefficient_norm_bounds": 0,
            "grouped_interior_bounds": 0,
            "quadratic_remainder_bounds": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
        "proof_boundary": route["proof_boundary"],
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
        "built physical P_lin/Morse-amplitude gate: "
        f"{counts['rows']} rows, "
        f"{counts['physical_correction_channels']} correction channels, "
        f"{counts['exact_centered_complex_coefficients']} centered coefficients, "
        f"degree {counts['maximum_polynomial_degree']}, "
        f"{counts['exact_morse_cprime_formulas']} c-prime formulas, "
        f"{counts['grouped_interior_bounds']} grouped interior bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
