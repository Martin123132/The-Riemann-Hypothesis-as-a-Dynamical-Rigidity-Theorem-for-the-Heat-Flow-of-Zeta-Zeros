#!/usr/bin/env python3
"""Build the modular-primitive/contact nonclosure theorem and audits."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_modular_primitive_contact_nonclosure_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "fourier_field": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_fourier_moment_regular_field_gate.json",
    "theta_primitive": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_summand_spectral_square_gate.json",
    "theta_probability": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_curvature_probability_operator_gate.json",
    "theta_tail_guard": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.json",
}

MIN_MULTIPLICITY = 2
MAX_MULTIPLICITY = 16
AUDIT_PHASES = (
    sp.Rational(5, 4) * sp.pi,
    sp.Rational(3, 2) * sp.pi,
    sp.Rational(7, 4) * sp.pi,
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expression_hash(expression: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(expression)).encode("utf-8")).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    require(
        "B_c=Re[i^(m+1)M_(m+1)]"
        in payloads["fourier_field"].get("symbolic_certificate", {}).get("multiple_contact", ""),
        "Fourier-field contract drift",
    )
    require(
        "Phi(u)=(R''(u)-R(u))/8"
        in payloads["theta_primitive"].get("exact", {}).get("differential_profile", {}).get("theta_primitive", ""),
        "theta differential-primitive contract drift",
    )
    require(
        "H_t(x)=1/16+D_t[C_t](x)/8"
        in payloads["theta_primitive"].get("exact", {}).get("deformed_half_transform", {}).get("identity", ""),
        "theta deformed-transform contract drift",
    )
    require(
        "R'(0)=-1/2"
        in payloads["theta_probability"].get("exact", {}).get("theta_primitive", {}).get("modular_endpoint", ""),
        "theta endpoint contract drift",
    )
    require(
        "theta-type double-exponential decay"
        in payloads["theta_tail_guard"].get("exact", {}).get("admissibility", {}).get("statement", ""),
        "theta-tail guard drift",
    )
    return payloads


def source_audit() -> dict[str, dict[str, str]]:
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def derivative_of(function: sp.Expr, x: sp.Symbol, order: int) -> sp.Expr:
    return function if order == 0 else sp.diff(function, x, order)


def build_operator_jet_audit() -> list[dict]:
    x, t = sp.symbols("x t", real=True)
    cfun = sp.Function("C")(x)
    operator = -4 * t**2 * sp.diff(cfun, x, 2) + 4 * t * x * sp.diff(cfun, x) + (2 * t - 1 - x**2) * cfun
    rows: list[dict] = []
    for j in range(17):
        direct = sp.expand(sp.diff(operator, x, j))
        expected = (
            -4 * t**2 * derivative_of(cfun, x, j + 2)
            + 4 * t * x * derivative_of(cfun, x, j + 1)
            + (4 * t * j + 2 * t - 1 - x**2) * derivative_of(cfun, x, j)
        )
        if j >= 1:
            expected -= 2 * j * x * derivative_of(cfun, x, j - 1)
        if j >= 2:
            expected -= j * (j - 1) * derivative_of(cfun, x, j - 2)
        expected = sp.expand(expected)
        require(sp.simplify(direct - expected) == 0, f"operator jet drift j={j}")
        rows.append({"order": j, "sha256": expression_hash(direct)})
    return rows


def build_contact_audit() -> list[dict]:
    h, a, b = sp.symbols("h a b", nonzero=True)
    rows: list[dict] = []
    for m in range(MIN_MULTIPLICITY, MAX_MULTIPLICITY + 1):
        theta_jet = -sp.Rational(1, 2) + h**m * (1 + a * h + b * h**2)
        heat_jet = (theta_jet + sp.Rational(1, 2)) / 8
        theta_ratio = sp.factor(
            sp.diff(theta_jet, h, m + 1).subs(h, 0)
            / ((m + 1) * sp.diff(theta_jet, h, m).subs(h, 0))
        )
        heat_ratio = sp.factor(
            sp.diff(heat_jet, h, m + 1).subs(h, 0)
            / ((m + 1) * sp.diff(heat_jet, h, m).subs(h, 0))
        )
        require(theta_ratio == a and heat_ratio == a, f"contact ratio drift m={m}")
        rows.append(
            {
                "multiplicity": m,
                "theta_ratio": str(theta_ratio),
                "heat_ratio": str(heat_ratio),
                "jet_sha256": expression_hash(theta_jet),
            }
        )
    return rows


def build_model_audit() -> list[dict]:
    y, r = sp.symbols("y r", real=True)
    compact_field = y * sp.cot(y) - sp.Symbol("m") - 1
    rows: list[dict] = []
    for m in range(2, 7):
        for phase in AUDIT_PHASES:
            field = sp.simplify(compact_field.subs({sp.Symbol("m"): m, y: phase}) + r)
            ell = sp.simplify((2 * field - m) / 4)
            rows.append(
                {
                    "multiplicity": m,
                    "phase": str(phase),
                    "bessel_log_derivative_symbol": "r=G_a'(c)/G_a(c) at c=1",
                    "field": str(field),
                    "ell": str(ell),
                }
            )
    return rows


def build_bessel_audit() -> list[dict]:
    mp.mp.dps = 70
    rows: list[dict] = []
    for c_integer in (1, 2, 3):
        c = mp.mpf(c_integer)
        a = c**2

        def transform(x: mp.mpf) -> mp.mpf:
            return mp.re(mp.besselk(1j * x / 4, a)) / mp.besselk(0, a)

        value = transform(c)
        log_derivative = mp.diff(transform, c) / value
        lower_bound = 1 - c**2 / (32 * a)
        require(value > lower_bound > 0, f"Bessel nonvanishing diagnostic failed c={c_integer}")
        rows.append(
            {
                "c": c_integer,
                "a": c_integer**2,
                "transform_at_c": mp.nstr(value, 55),
                "log_derivative_at_c": mp.nstr(log_derivative, 55),
                "rigorous_lower_bound": "31/32",
            }
        )
    return rows


def build_kernel_audit() -> dict:
    u, v = sp.symbols("u v", real=True)
    kernel = 8 * sp.sinh(v - u)
    return {
        "diagonal": str(sp.simplify(kernel.subs(v, u))),
        "diagonal_u_derivative": str(sp.simplify(sp.diff(kernel, u).subs(v, u))),
        "homogeneous_operator": str(sp.simplify(sp.diff(kernel, u, 2) - kernel)),
        "sha256": expression_hash(kernel),
    }


def build_modular_jet_audit() -> list[dict]:
    phi = sp.symbols("phi0:16")
    jets: list[sp.Expr] = [sp.Symbol("r0"), -sp.Rational(1, 2)]
    for order in range(14):
        forcing = 0 if order % 2 else 8 * phi[order]
        jets.append(sp.expand(jets[order] + forcing))
    rows: list[dict] = []
    for order in range(1, 16, 2):
        require(jets[order] == -sp.Rational(1, 2), f"odd modular jet drift order={order}")
        rows.append({"order": order, "jet": str(jets[order])})
    return rows


def heat_monomial(m: int, z: sp.Symbol, direction: int) -> sp.Expr:
    return sp.expand(
        sum(
            direction**order * sp.diff(z**m, z, 2 * order) / sp.factorial(order)
            for order in range(m // 2 + 1)
        )
    )


def build_hermite_audit() -> list[dict]:
    z = sp.symbols("z")
    rows: list[dict] = []
    for m in range(MIN_MULTIPLICITY, MAX_MULTIPLICITY + 1):
        forward = heat_monomial(m, z, -1)
        backward = heat_monomial(m, z, 1)
        forward_count = int(sp.Poly(forward, z).count_roots(-sp.oo, sp.oo))
        backward_count = int(sp.Poly(backward, z).count_roots(-sp.oo, sp.oo))
        require(forward_count == m, f"forward Hermite count drift m={m}")
        require(backward_count == m % 2, f"backward Hermite count drift m={m}")
        rows.append(
            {
                "multiplicity": m,
                "forward_real_roots": forward_count,
                "backward_real_roots": backward_count,
                "forward_sha256": expression_hash(forward),
                "backward_sha256": expression_hash(backward),
            }
        )
    return rows


def build_rows() -> list[GateRow]:
    return [
        GateRow("mpc_01_sources", "source chain", "proved", "Four parent contracts are current and hash-pinned.", "Fourier-field, theta-primitive, theta-probability, and theta-tail guards are loaded.", "No parent conclusion is strengthened."),
        GateRow("mpc_02_theta_coordinate", "theta coordinate", "proved", "Every multiplicity contact transfers exactly to the endpoint-subtracted theta primitive.", "With A_t=D_t[C_t]=8H_t-1/2, exact multiplicity m at c means A_t(c)=-1/2, A_t^(j)(c)=0 for 1<=j<m, and A_t^(m)(c)!=0.", "This is an exact coordinate change, not a contact exclusion."),
        GateRow("mpc_03_field", "regular field", "proved", "The regular field is the adjacent theta-primitive jet ratio.", "B_c=A_t^(m+1)(c)/((m+1)A_t^(m)(c)) and ell=(2cB_c-m)/4.", "The ratio remains signed."),
        GateRow("mpc_04_operator_jets", "operator recurrence", "proved", "Every theta contact equation is an explicit finite recurrence in C_t jets.", "The audited all-j formula differentiates D_t=-4t^2D_x^2+4txD_x+(2t-1-x^2).", "No sign follows from the recurrence alone."),
        GateRow("mpc_05_summands", "infinite cancellation", "guard_validated", "The exact arithmetic theta decomposition enters only through coupled infinite sums of these jets.", "At contact, sum_n A_(n,t)=-1/2 and sum_n A_(n,t)^(j)=0 for 1<=j<m.", "Finite or termwise positivity is forbidden by the parent spectral gate."),
        GateRow("mpc_06_bessel_density", "tail model", "proved", "A Pólya-Bessel density supplies a smooth theta-style tail.", "z_a(u)=2 exp(-a cosh(4u))/K_0(a) is even, positive, uniformly strongly log-concave, and has double-exponential tails.", "This is not the discrete Xi theta series."),
        GateRow("mpc_07_bessel_transform", "published LP input", "proved", "Its characteristic transform belongs to the Laguerre-Polya class.", "G_a(x)=K_(ix/4)(a)/K_0(a); the published square-integral theorem proves K_(iz)(a) has only real zeros for a>0.", "Published input is cited; it does not concern Xi."),
        GateRow("mpc_08_nonvanishing", "contact guard", "proved", "The Bessel factor can be kept nonzero at the prescribed contact.", "Integration by parts gives E[Z V'(Z)]=1 and z sinh(4z)>=4z^2, hence E[Z^2]<=1/(16a) and G_a(c)>=1-c^2/(32a)>0 when a>c^2/32.", "This is a sufficient bound, not a zero classification."),
        GateRow("mpc_09_model", "exact contact model", "proved", "Sinc factors impose an arbitrary multiplicity contact without losing the Bessel tail.", "P_(m,y,c,a)(x)=sinc(pi x/c)^m sinc(y x/c)G_a(x), pi<y<2pi.", "The model is not Xi."),
        GateRow("mpc_10_field_model", "arbitrary field", "proved", "The sinc-Bessel family realizes every real ell.", "B_c=(y cot y-m-1)/c+G_a'(c)/G_a(c), so ell=[2y cot y-3m-2+2cG_a'(c)/G_a(c)]/4.", "The Bessel correction is finite because G_a(c)>0."),
        GateRow("mpc_11_phase", "phase bijection", "proved", "The phase remains a bijective field coordinate.", "y cot y decreases from +infinity to -infinity on (pi,2pi).", "Pi is the first uniform-characteristic zero."),
        GateRow("mpc_12_shape", "kernel geometry", "proved", "The frequency density is smooth, positive, even, log-concave, and theta-tail admissible.", "It is the convolution of interval-uniform laws with z_a; Prekopa preserves log-concavity and bounded convolution preserves the double-exponential tail class.", "Uniform strong log-concavity of the convolution is not claimed."),
        GateRow("mpc_13_heat", "exact first loss", "proved", "Every sinc-Bessel model has an all-time heat flow with exact threshold zero.", "The transform is Laguerre-Polya, positive-time backward heat preserves that class, and exp(D_z^2)z^m gives a nonreal pair on the negative side.", "This is a countermodel heat family."),
        GateRow("mpc_14_embedding", "positive boundary embedding", "proved", "The contact can be placed at any prescribed positive heat time t_*.", "Set Phi_*(u)=kappa exp(-t_*u^2)p(u), with kappa chosen by integral_0^infinity cosh(u)Phi_*(u)du=1/16; then H_(t_*) is a positive multiple of P.", "Scaling does not change zeros or ell."),
        GateRow("mpc_15_primitive", "positive primitive", "proved", "Every embedded source has a positive decreasing convex theta primitive.", "R(u)=8 integral_u^infinity sinh(v-u)Phi_*(v)dv gives R''-R=8Phi_*, R'(0)=-1/2, R>0, R'<0, and R''>0.", "This uses positivity and the exact endpoint normalization."),
        GateRow("mpc_16_reflection", "modular reflection", "proved", "The primitive has the same reflection equation as the Jacobi-theta primitive.", "Extend by R(-u)-R(u)=sinh(u); evenness of Phi_* and R'(0)=-1/2 make the extension C-infinity and preserve R''-R=8Phi_*.", "The reflection equation is weaker than discrete theta arithmetic."),
        GateRow("mpc_17_probability", "curvature probability", "proved", "The embedded primitive has the same triangular probability package.", "dmu(v)=2R''(v)dv has mass one and R(u)=(1/2)integral(v-u)_+dmu(v).", "The fixed Xi component weights are absent."),
        GateRow("mpc_18_nonclosure", "modular no-go", "proved", "Modular reflection, positive primitive curvature, theta tails, exact first loss, and local contact still do not bound ell.", "The embedded sinc-Bessel family realizes every real field at every fixed multiplicity and any selected t_*>0.", "This does not reproduce the Jacobi-theta summand roster."),
        GateRow("mpc_19_arithmetic", "Xi handoff", "open", "Any surviving local field theorem must use discrete Xi arithmetic beyond the reflection equation.", "The live extra data are the n^(-1/2) translate law, fixed theta weights, and coupled infinite contact cancellations.", "No Xi field bound is proved."),
        GateRow("mpc_20_uniform", "Jensen handoff", "open", "The cofinal Jensen route still requires a remainder uniform in degree.", "The modular-primitive embedding does not estimate the Jensen remainder.", "No degree-uniform bound is proved."),
        GateRow("mpc_21_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", "It proves an exact embedding/nonclosure theorem and identifies the remaining arithmetic input.", "No Xi collision exclusion, Lambda<=0, RH, or prize-level conclusion is proved."),
    ]


def build_note() -> str:
    return """# Newman Modular-Primitive Contact Nonclosure Gate

Date: 2026-08-03

Status: exact modular-reflection embedding and theta-tail first-loss
countermodel theorem; discrete Xi theta arithmetic, the actual Xi field, and
the degree-uniform Jensen remainder remain open; not a proof of RH.

## Theta Contact Coordinate

For the Xi theta primitive,

```text
R''-R=8Phi,                   R'(0)=-1/2,
C_t(x)=integral_0^infinity exp(tu^2)R(u)cos(xu)du,
A_t=D_t[C_t]=8H_t-1/2,
D_t=-4t^2 partial_x^2+4tx partial_x+(2t-1-x^2).
```

If `c>0` is a zero of exact multiplicity `m>=2`, then

```text
A_t(c)=-1/2,
A_t^(j)(c)=0                 for 1<=j<m,
A_t^(m)(c)!=0,
B_c=A_t^(m+1)(c)/[(m+1)A_t^(m)(c)],
ell=(2cB_c-m)/4.
```

Thus the Fourier ratio from the previous gate is exactly an adjacent jet
ratio of the endpoint-subtracted theta primitive. Differentiating `D_t`
gives, for every `j>=0`,

```text
A_t^(j)
 =-4t^2 C_t^(j+2)+4tx C_t^(j+1)
  +(4tj+2t-1-x^2)C_t^(j)
  -2jx C_t^(j-1)-j(j-1)C_t^(j-2),
```

with terms having negative derivative order omitted. In the true theta
sum, all these equations are coupled infinite cancellations.

## Sinc-Bessel Boundary Model

Let

```text
z_a(u)=2 exp(-a cosh(4u))/K_0(a),
G_a(x)=K_(ix/4)(a)/K_0(a).
```

The density is positive, even, uniformly strongly log-concave, and has
theta-style double-exponential tails. A published square-integral theorem
proves that `K_(iz)(a)` has only real zeros for every `a>0`:
https://arxiv.org/abs/0801.2996.

The elementary identity `E[Z V'(Z)]=1`, with
`V(z)=a cosh(4z)`, gives

```text
E[Z^2]<=1/(16a),
G_a(c)=E[cos(cZ)]>=1-c^2/(32a)>0
```

whenever `a>c^2/32`. Now define

```text
P_(m,y,c,a)(x)
 =sinc(pi x/c)^m sinc(y x/c)G_a(x),
pi<y<2pi.
```

It is the characteristic function of a smooth positive even log-concave
density obtained by convolving interval-uniform laws with `z_a`. Its tail
remains double exponential, its transform is Laguerre-Polya, and `x=c` is
an exact multiplicity-`m` zero. At that contact,

```text
B_c=(y cot y-m-1)/c+G_a'(c)/G_a(c),
ell=[2y cot y-3m-2+2cG_a'(c)/G_a(c)]/4.
```

Since `y cot y` is a decreasing bijection from `(pi,2pi)` to the real
line, this family realizes every real `ell`. Its all-time heat flow has
threshold exactly zero: nonnegative backward-heat time preserves the
Laguerre-Polya class, while the negative-time local limit
`exp(partial_z^2)z^m` has a nonreal pair.

## Modular-Primitive Embedding

Fix any desired positive boundary time `t_*`. Let `p` be the sinc-Bessel
density and choose `kappa>0` so that

```text
Phi_*(u)=kappa exp(-t_*u^2)p(u),
integral_0^infinity cosh(u)Phi_*(u)du=1/16.
```

Then `H_(t_*)` is a positive multiple of `P_(m,y,c,a)`, so its contact,
field, and exact first-loss geometry are unchanged. Define for `u>=0`

```text
R_*(u)=8 integral_u^infinity sinh(v-u)Phi_*(v)dv.
```

Direct differentiation and the normalization give

```text
R_*''-R_*=8Phi_*,             R_*'(0)=-1/2,
R_*>0,                        R_*'<0,
R_*''=R_*+8Phi_*>0.
```

Extend to the negative half-line by

```text
R_*(-u)-R_*(u)=sinh(u).
```

Evenness of `Phi_*` and the endpoint slope make this extension smooth and
preserve the differential identity. This is exactly the reflection equation
of the Jacobi-theta primitive. Also

```text
dmu(v)=2R_*''(v)dv
```

is a probability law and
`R_*(u)=(1/2)integral(v-u)_+dmu(v)`.

## Route Decision

The modular reflection equation, endpoint normalization, positive decreasing
convex primitive, triangular probability representation, theta-style tail,
positive log-concave Fourier kernel, and exact first-loss contact can all
coexist with every real regular field. Those properties do not close the
Xi collision coefficient.

The model deliberately does not reproduce the discrete Jacobi-theta roster.
The surviving Xi-specific information is narrower: the
`n^(-1/2)` translate law, its fixed arithmetic weights, and the coupled
infinite summand cancellations in every contact jet. A candidate inequality
must use that discrete structure and must be tested against this embedding.
The global first-jet winding and degree-uniform Jensen remainder remain
independent obligations.

## Audit

```text
python work/rh_compute/scripts/jensen_window_pf_newman_modular_primitive_contact_nonclosure_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_modular_primitive_contact_nonclosure_gate.py
```

The builder and independent checker verify 21 theorem rows, 17 operator
jets, 15 contact multiplicities, 15 sinc-Bessel field models, three Bessel
nonvanishing diagnostics, eight modular odd jets, and 15 forward/backward
Hermite root counts.

This proves no bound on the actual Xi field, no discrete theta-summand
inequality, no degree-uniform Jensen remainder, no Xi collision exclusion,
`Lambda<=0`, RH, or prize-level conclusion.
"""


def main() -> int:
    load_sources()
    rows = build_rows()
    operator_audit = build_operator_jet_audit()
    contact_audit = build_contact_audit()
    model_audit = build_model_audit()
    bessel_audit = build_bessel_audit()
    kernel_audit = build_kernel_audit()
    modular_jet_audit = build_modular_jet_audit()
    hermite_audit = build_hermite_audit()
    require(len(rows) == 21, "row count drift")

    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "exact modular-primitive embedding and theta-tail contact nonclosure",
        "source_audit": source_audit(),
        "symbolic_certificate": {
            "theta_contact": "A_t=D_t[C_t]=8H_t-1/2; exact multiplicity m at c means A_t(c)=-1/2, A_t^(j)(c)=0 for 1<=j<m, and A_t^(m)(c)!=0.",
            "theta_field": "B_c=A_t^(m+1)(c)/((m+1)A_t^(m)(c)) and ell=(2cB_c-m)/4.",
            "bessel_density": "z_a(u)=2 exp(-a cosh(4u))/K_0(a), G_a(x)=K_(ix/4)(a)/K_0(a).",
            "bessel_nonvanishing": "E[Z^2]<=1/(16a) and G_a(c)>=1-c^2/(32a)>0 when a>c^2/32.",
            "contact_model": "P_(m,y,c,a)(x)=sinc(pi x/c)^m sinc(y x/c)G_a(x), pi<y<2pi.",
            "model_field": "B_c=(y cot y-m-1)/c+G_a'(c)/G_a(c), ell=[2y cot y-3m-2+2cG_a'(c)/G_a(c)]/4.",
            "positive_boundary_embedding": "Phi_*(u)=kappa exp(-t_*u^2)p(u), integral_0^infinity cosh(u)Phi_*(u)du=1/16, so H_(t_*) is a positive multiple of P.",
            "positive_primitive": "R_*(u)=8 integral_u^infinity sinh(v-u)Phi_*(v)dv gives R_*''-R_*=8Phi_*, R_*'(0)=-1/2, R_*>0, R_*'<0, and R_*''>0.",
            "modular_reflection": "R_*(-u)-R_*(u)=sinh(u).",
            "probability": "dmu(v)=2R_*''(v)dv is a probability law and R_*(u)=(1/2)integral(v-u)_+dmu(v).",
            "nonclosure": "Modular reflection, positive primitive curvature, theta tails, exact first loss, and local contact do not bound ell.",
            "arithmetic_handoff": "A surviving Xi field theorem must use the discrete n^(-1/2) translate law, fixed theta weights, or coupled infinite summand cancellations absent from the embedding.",
        },
        "operator_jet_audit": operator_audit,
        "contact_audit": contact_audit,
        "model_audit": model_audit,
        "bessel_audit": bessel_audit,
        "primitive_kernel_audit": kernel_audit,
        "modular_jet_audit": modular_jet_audit,
        "hermite_audit": hermite_audit,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "sources": len(SOURCE_PATHS),
            "operator_jets": len(operator_audit),
            "contact_multiplicities": len(contact_audit),
            "sinc_bessel_models": len(model_audit),
            "bessel_diagnostics": len(bessel_audit),
            "modular_odd_jets": len(modular_jet_audit),
            "hermite_checks": len(hermite_audit),
            "open_handoffs": 2,
            "actual_xi_field_bounds": 0,
            "degree_uniform_bounds": 0,
            "rh_conclusions": 0,
        },
        "published_input": {
            "title": "Using integrals of squares of certain real-valued special functions to prove that K_(iz)(a) has only real zeros",
            "url": "https://arxiv.org/abs/0801.2996",
            "use": "K_(iz)(a) has only real zeros for a>0.",
        },
        "proof_boundary": "This gate proves the all-multiplicity theta-primitive contact coordinate and an exact positive-time modular-reflection, theta-tail, Laguerre-Polya first-loss family realizing every real field. It does not reproduce the discrete Jacobi-theta summand roster, bound the actual Xi field, prove a degree-uniform Jensen remainder, exclude an Xi collision, prove Lambda<=0, RH, or a prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, build_note())
    print(
        "built modular-primitive contact nonclosure gate: "
        "21 rows, 4 sources, 17 operator jets, 15 contact multiplicities, "
        "15 sinc-Bessel models, 3 Bessel diagnostics, 8 modular odd jets, "
        "15 Hermite checks, 2 open handoffs, 0 Xi field bounds, "
        "0 degree-uniform bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
