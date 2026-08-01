#!/usr/bin/env python3
"""Build the Jordan-Muntz causal-energy and natural-mollifier bridge."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_jordan_muntz_causal_energy_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md"
)


def row(
    row_id: str,
    role: str,
    status: str,
    statement: str,
    boundary: str,
) -> dict[str, str]:
    return {
        "id": row_id,
        "role": role,
        "status": status,
        "statement": statement,
        "proof_boundary": boundary,
    }


def build_payload() -> dict:
    rows = [
        row(
            "jmce_01_dirichlet_series",
            "exact_identity",
            "available_exact",
            "For c_omega(n)=sum_(dm=n)mu(d)d^(-omega)m^omega, "
            "sum c_omega(n)n^(-s)=zeta(s-omega)/zeta(s+omega) "
            "when Re(s)>1+omega.",
            "Absolute Dirichlet-series identity only.",
        ),
        row(
            "jmce_02_jordan_error",
            "definition",
            "available_exact",
            "C_omega(x)=sum_(n<=x)c_omega(n), "
            "A_omega=1/((1+omega)zeta(1+2omega)), and "
            "E_omega(x)=C_omega(x)-A_omega*x^(1+omega).",
            "Definition and residue normalization only.",
        ),
        row(
            "jmce_03_power_remainder",
            "definition",
            "available_exact",
            "R_omega(y)=sum_(m<=y)m^omega-y^(1+omega)/(1+omega), "
            "with R_omega(y)=-y^(1+omega)/(1+omega) for 0<y<1.",
            "Definition of the generalized fractional-part remainder.",
        ),
        row(
            "jmce_04_full_muntz_identity",
            "exact_identity",
            "available_exact",
            "E_omega(x)=sum_(d>=1)mu(d)d^(-omega)R_omega(x/d).",
            "The d>x tail is part of the identity and must not be dropped.",
        ),
        row(
            "jmce_05_tail_absolute",
            "exact_convergence",
            "available_exact",
            "For d>x the summand is "
            "-mu(d)x^(1+omega)/((1+omega)d^(1+2omega)); "
            "hence the tail is absolutely convergent for omega>0.",
            "Pointwise absolute tail convergence is not Hilbert convergence.",
        ),
        row(
            "jmce_06_remainder_hilbert_space",
            "exact_bound",
            "available_exact",
            "R_omega(y)=O_omega(y^omega) at infinity and is "
            "-y^(1+omega)/(1+omega) below one, so "
            "R_omega belongs to H=L2((0,infinity),dx/x^2) exactly in "
            "the needed range 0<omega<1/2.",
            "This concerns the base remainder, not the infinite Mobius sum.",
        ),
        row(
            "jmce_07_remainder_mellin",
            "exact_transform",
            "available_exact",
            "Integral_0^infinity R_omega(x)x^(-s)dx/x "
            "=zeta(s-omega)/s for omega<Re(s)<1+omega.",
            "Generalized Muntz transform in its honest convergence strip.",
        ),
        row(
            "jmce_08_unitary_dilations",
            "exact_identity",
            "available_exact",
            "U_d f(x)=sqrt(d)f(x/d) is unitary on H, and "
            "d^(-omega)R_omega(x/d)="
            "d^(-omega-1/2)U_dR_omega(x).",
            "Unitary orbit identity only.",
        ),
        row(
            "jmce_09_partial_mobius_sum",
            "definition",
            "available_exact",
            "E_(omega,N)=sum_(d<=N)mu(d)d^(-omega)R_omega(x/d).",
            "Every E_(omega,N) is a finite H vector.",
        ),
        row(
            "jmce_10_partial_mellin",
            "exact_transform",
            "available_exact",
            "On Re(s)=1/2, Mellin(E_(omega,N))(s)="
            "zeta(s-omega)M_N(s+omega)/s, where "
            "M_N(u)=sum_(d<=N)mu(d)d^(-u).",
            "Finite Dirichlet polynomial; no limiting interchange is used.",
        ),
        row(
            "jmce_11_partial_plancherel",
            "exact_identity",
            "available_exact",
            "||E_(omega,N)||_H^2=(1/(2pi))integral_R "
            "|zeta(1/2-omega+it)/(1/2+it)|^2"
            "|M_N(1/2+omega+it)|^2dt.",
            "Exact finite-N norm, not a bound uniform in N.",
        ),
        row(
            "jmce_12_mollified_right_line",
            "exact_identity",
            "available_exact",
            "The functional equation rewrites the integrand as "
            "W_omega(t)|zeta(sigma+it)M_N(sigma+it)|^2, "
            "sigma=1/2+omega and "
            "W_omega(t)=|chi(1/2-omega+it)|^2/|1/2+it|^2.",
            "An exact weighted natural-mollifier norm on a line right of "
            "one half.",
        ),
        row(
            "jmce_13_pointwise_limit",
            "exact_convergence",
            "available_exact",
            "E_(omega,N)(x) tends pointwise to E_omega(x) because the "
            "tail after d>x is absolutely convergent.",
            "Pointwise convergence does not imply H convergence.",
        ),
        row(
            "jmce_14_fatou_handoff",
            "exact_implication",
            "available_exact",
            "sup_N ||E_(omega,N)||_H<infinity implies E_omega in H by "
            "Fatou and the pointwise limit.",
            "One-way handoff; no Bessel property is assumed.",
        ),
        row(
            "jmce_15_positive_time_energy",
            "exact_identity",
            "available_exact",
            "E_omega in H iff J_omega="
            "integral_1^infinity |E_omega(x)|^2dx/x^2<infinity, "
            "since the known x<1 contribution is "
            "A_omega^2/(1+2omega).",
            "Exact removal of the harmless negative-log-time tail.",
        ),
        row(
            "jmce_16_positive_laplace_transform",
            "exact_transform",
            "available_exact",
            "For s=z+1/2, the positive-time transform is "
            "G_omega(s)=zeta(s-omega)/(s*zeta(s+omega))"
            "-A_omega/(s-1-omega), initially for Re(s)>1-omega.",
            "The second term removes the residue-main-term pole at "
            "s=1+omega.",
        ),
        row(
            "jmce_17_boundary_energy_automatic",
            "route_guard",
            "guard_validated",
            "G_omega(1/2+it) is in L2(dt) after removable boundary "
            "zeros: the zeta-ratio term is O(|t|^(omega-1)) and "
            "omega<1/2.",
            "Finite boundary energy does not prove positive-time Hardy "
            "support.",
        ),
        row(
            "jmce_18_l2_implies_fixed_shift",
            "theorem_candidate",
            "theorem_candidate",
            "J_omega<infinity makes G_omega an H2 function in "
            "Re(s)>1/2; therefore every denominator zero in that "
            "half-plane is canceled and Suzuki's D(omega) holds.",
            "Paley-Wiener implication; independent expert review required.",
        ),
        row(
            "jmce_19_fixed_shift_implies_l2",
            "theorem_candidate",
            "theorem_candidate",
            "If D(omega) holds, the completed quotient Q_omega is inner. "
            "The completed gamma quotient, its single removed pole, and "
            "the O(|t|^(omega-1)) boundary estimate put G_omega in H2, "
            "so transform uniqueness gives J_omega<infinity.",
            "Uses the audited fixed-shift innerness theorem and the "
            "standard gamma-quotient Hardy estimate.",
        ),
        row(
            "jmce_20_fixed_shift_equivalence",
            "theorem_candidate",
            "theorem_candidate",
            "For every fixed 0<omega<1/2, "
            "J_omega<infinity iff D(omega).",
            "New corpus theorem candidate, not an unconditional J bound.",
        ),
        row(
            "jmce_21_cofinal_error_criterion",
            "theorem_candidate",
            "theorem_candidate",
            "RH iff there is a sequence omega_j->0 with "
            "integral_1^infinity |E_(omega_j)(x)|^2dx/x^2 finite for "
            "every j.",
            "Combines the fixed-shift candidate with the audited cofinal "
            "phase theorem.",
        ),
        row(
            "jmce_22_rh_natural_convergence",
            "conditional_estimate",
            "conditional_on_rh",
            "Under RH, Balazard-Saias natural-Mobius approximation and "
            "Burnol's shifted-ratio estimate give "
            "E_(omega,N)->E_omega in H for every fixed "
            "0<omega<1/2.",
            "Conditional theorem technology; it cannot supply the missing "
            "unconditional bound.",
        ),
        row(
            "jmce_23_cofinal_partial_norm",
            "theorem_candidate",
            "theorem_candidate",
            "RH iff there exist omega_j->0 such that "
            "sup_N ||E_(omega_j,N)||_H<infinity for every j.",
            "RH gives convergence; the reverse direction uses Fatou, the "
            "error-energy criterion, and cofinal accumulation.",
        ),
        row(
            "jmce_24_finite_energy_formula",
            "exact_identity",
            "available_exact",
            "For X>=1, J_omega(X)=integral_1^X E_omega(x)^2dx/x^2 "
            "equals sum_(n,m<=X)c_n*c_m*(1/max(n,m)-1/X)"
            "-(2A_omega/omega)sum_(n<=X)c_n*(X^omega-n^omega)"
            "+A_omega^2*(X^(1+2omega)-1)/(1+2omega).",
            "All three terms must remain together.",
        ),
        row(
            "jmce_25_signed_measure_energy",
            "exact_identity",
            "available_exact",
            "With eta_(omega,X)=sum_(n<=X)c_n delta_n"
            "-A_omega(1+omega)u^omega 1_(0,X](u)du, "
            "J_omega(X)=double_integral K_X(u,v)deta(u)deta(v), where "
            "K_X=1/max(1,u,v)-1/X on [0,X]^2.",
            "Exact weighted H-minus-one discrepancy norm.",
        ),
        row(
            "jmce_26_max_kernel_tn",
            "exact_total_positivity",
            "available_exact",
            "K_X(u,v)=integral_1^X 1_(u<=x)1_(v<=x)dx/x^2 is positive "
            "semidefinite; in reciprocal coordinates it is a Brownian "
            "min kernel and every ordered finite matrix is totally "
            "nonnegative by a cumulative-matrix factorization.",
            "Kernel total nonnegativity gives J>=0, not a uniform upper "
            "bound.",
        ),
        row(
            "jmce_27_tn_nonclosure_guard",
            "route_guard",
            "guard_validated",
            "The positive-positive, cross, and continuum-continuum pieces "
            "of the finite formula each have X^(1+2omega)-scale; only "
            "their signed combination can stay bounded.",
            "A proof using positivity of the max kernel alone cannot close "
            "the criterion.",
        ),
        row(
            "jmce_28_causality_countermodel",
            "countermodel_gate",
            "guard_validated",
            "1/(a+z) and 1/(a-z) have the same L2 boundary norm 1/(2a). "
            "The first is the transform of e^(-at)1_(t>0); the second has "
            "anti-causal boundary inverse e^(at)1_(t<0), while its actual "
            "positive-time Bromwich inverse is -e^(at).",
            "Equal finite spectral energy does not distinguish the pole "
            "side; Hardy support does.",
        ),
        row(
            "jmce_29_l2_coefficients_guard",
            "countermodel_gate",
            "guard_validated",
            "The coefficients d^(-omega-1/2) are square summable, but "
            "square-summable coefficients do not force convergence of a "
            "series of unit vectors without a Bessel bound.",
            "Generic model v_d=v and coefficients 1/d is l2 while its "
            "vector partial sums diverge.",
        ),
        row(
            "jmce_30_burnol_fit",
            "literature_fit",
            "source_fit",
            "Burnol's f_epsilon theorem has the same shifted-zeta-ratio, "
            "square-integrability, and causality architecture. The present "
            "R_omega identity is its fractional-power Jordan analogue.",
            "Literature fit and proof template, not a transferred "
            "unconditional estimate.",
        ),
        row(
            "jmce_31_nyman_nontransfer",
            "literature_guard",
            "guard_validated",
            "Baez-Duarte's general strong Nyman-Beurling sufficiency "
            "requires compact support in the audited theorem, while "
            "R_omega is not compactly supported; arbitrary closure also "
            "does not prove convergence of this fixed natural Mobius "
            "series.",
            "The general closure theorem cannot be cited as the missing "
            "natural-mollifier bound.",
        ),
        row(
            "jmce_32_open_mollifier_gate",
            "open_arithmetic_gate",
            "open_target",
            "For one explicit omega_j->0, prove uniformly in N the exact "
            "finite inequality integral_R W_(omega_j)(t)"
            "|zeta(1/2+omega_j+it)M_N(1/2+omega_j+it)|^2dt<C_j.",
            "This all-height natural-mollifier estimate is not proved; "
            "finite t grids and finite X energies do not close it.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_jordan_muntz_causal_energy_bridge",
        "date": "2026-07-23",
        "status": (
            "exact Jordan-Muntz and finite-energy reductions with two "
            "cofinal theorem candidates and one open mollifier gate"
        ),
        "proof_boundary": (
            "The exact identities expose a finite natural-mollifier norm "
            "and a totally nonnegative max-kernel discrepancy energy. "
            "They do not prove the required uniform bound, RH, or "
            "Lambda<=0. Fixed-shift and cofinal equivalences are internally "
            "audited theorem candidates requiring independent review."
        ),
        "rows": rows,
        "exact": {
            "hilbert_space": "H=L2((0,infinity),dx/x^2)",
            "remainder": (
                "R_omega(y)=sum_(m<=y)m^omega-y^(1+omega)/(1+omega)"
            ),
            "full_muntz": (
                "E_omega=sum_(d>=1)mu(d)d^(-omega)R_omega(x/d)"
            ),
            "partial_mellin": (
                "Mellin(E_(omega,N))(s)="
                "zeta(s-omega)M_N(s+omega)/s"
            ),
            "positive_transform": (
                "G_omega(s)=zeta(s-omega)/(s*zeta(s+omega))"
                "-A_omega/(s-1-omega)"
            ),
            "finite_energy_kernel": (
                "K_X(u,v)=1/max(1,u,v)-1/X"
            ),
        },
        "sources": [
            {
                "title": (
                    "On monotonicity of certain weighted summatory "
                    "functions associated with L-functions"
                ),
                "url": "https://arxiv.org/abs/1204.1823",
                "use": "Jordan coefficients and shifted quotient",
            },
            {
                "title": (
                    "On an analytic estimate in the theory of the Riemann "
                    "Zeta function and a Theorem of Baez-Duarte"
                ),
                "url": "https://arxiv.org/abs/math/0202166",
                "use": (
                    "shifted-ratio L2 theorem, natural Mobius "
                    "approximation, and causality template"
                ),
            },
            {
                "title": (
                    "A general strong Nyman-Beurling Criterion for the "
                    "Riemann Hypothesis"
                ),
                "url": "https://arxiv.org/abs/math/0505453",
                "use": "Muntz dilation framework and compact-support guard",
            },
            {
                "title": "Suzuki cofinal L2 hierarchy",
                "url": (
                    "outputs/"
                    "jensen_window_pf_suzuki_cofinal_l2_hierarchy.md"
                ),
                "use": "fixed-shift and cofinal Hardy handoff",
            },
        ],
        "audit": {
            "row_count": len(rows),
            "theorem_candidate_count": 5,
            "conditional_estimate_count": 1,
            "guard_count": 5,
            "open_arithmetic_gate_count": 1,
            "finite_energy_formula_checked": True,
            "muntz_identity_checked": True,
            "max_kernel_tn_checked": True,
            "uniform_mollifier_bound_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note(payload: dict) -> str:
    return """# Jensen-Window PF Jordan-Muntz Causal-Energy Bridge

Date: 2026-07-23

Status: exact Jordan-Muntz and finite-energy reductions with internally
audited fixed-shift/cofinal theorem candidates and one open natural-mollifier
gate. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_jordan_muntz_causal_energy_bridge.json
python work/rh_compute/scripts/jensen_window_pf_jordan_muntz_causal_energy_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_causal_energy_bridge.py
```

## Generalized Muntz Identity

For `0<omega<1/2`, retain

```text
c_omega(n)
 =sum_(d*m=n) mu(d)d^(-omega)m^omega,

C_omega(x)=sum_(n<=x)c_omega(n),

A_omega=1/((1+omega)*zeta(1+2omega)),

E_omega(x)=C_omega(x)-A_omega*x^(1+omega).
```

Define the fractional-power remainder

```text
R_omega(y)
 =sum_(m<=y)m^omega-y^(1+omega)/(1+omega).           (JMCE.1)
```

For `0<y<1`, this is exactly
`-y^(1+omega)/(1+omega)`. Expanding `C_omega` by its Dirichlet
convolution and using

```text
sum_(d>=1)mu(d)d^(-(1+2omega))=1/zeta(1+2omega)
```

gives the full identity

```text
E_omega(x)
 =sum_(d>=1)mu(d)d^(-omega)R_omega(x/d).             (JMCE.2)
```

This version includes the earlier `d>x` tail rather than writing it
separately. That tail is absolutely convergent because its terms are

```text
-mu(d)*x^(1+omega)/((1+omega)*d^(1+2omega)).
```

The elementary power-sum estimate gives
`R_omega(y)=O_omega(y^omega)` at infinity. Together with its exact
behavior below one, this proves

```text
R_omega in H,
H=L2((0,infinity),dx/x^2),                           (JMCE.3)
```

precisely for `0<omega<1/2`. Its Mellin transform in the honest
convergence strip is

```text
integral_0^infinity R_omega(x)x^(-s)dx/x
 =zeta(s-omega)/s,
omega<Re(s)<1+omega.                                 (JMCE.4)
```

## Natural Mobius Dilations

The operators

```text
(U_d f)(x)=sqrt(d)f(x/d)
```

are unitary on `H`. Hence the finite natural approximants are

```text
E_(omega,N)(x)
 =sum_(d<=N)mu(d)d^(-omega)R_omega(x/d)
 =sum_(d<=N)mu(d)d^(-omega-1/2)U_dR_omega(x).        (JMCE.5)
```

Put

```text
M_N(u)=sum_(d<=N)mu(d)d^(-u).
```

Mellin-Plancherel gives the exact finite identity

```text
||E_(omega,N)||_H^2
 =1/(2*pi)*integral_R
  |zeta(1/2-omega+i*t)/(1/2+i*t)|^2
  |M_N(1/2+omega+i*t)|^2 dt.                         (JMCE.6)
```

No infinite series has been moved through this integral. The zeta functional
equation rewrites the same integrand as

```text
W_omega(t)
*|zeta(sigma+i*t)M_N(sigma+i*t)|^2,

sigma=1/2+omega,

W_omega(t)
 =|chi(1/2-omega+i*t)|^2/|1/2+i*t|^2,               (JMCE.7)
```

where `zeta(v)=chi(v)zeta(1-v)`. Thus the live arithmetic object is a
weighted all-height norm of the ordinary natural mollifier on a line
strictly to the right of one half.

For every fixed `x`, (JMCE.2) gives

```text
E_(omega,N)(x)->E_omega(x).
```

Consequently,

```text
sup_N ||E_(omega,N)||_H < infinity
=> E_omega in H                                      (JMCE.8)
```

by Fatou. Square summability of the scalar coefficients
`d^(-omega-1/2)` is not enough: the dilation orbit has not been proved to
be a Bessel sequence.

## Causal Error Criterion

The part below `x=1` is explicit:

```text
integral_0^1 |E_omega(x)|^2 dx/x^2
 =A_omega^2/(1+2omega).
```

Therefore define the only open part by

```text
J_omega
 =integral_1^infinity |E_omega(x)|^2 dx/x^2.         (JMCE.9)
```

With `s=z+1/2`, its positive-log-time Laplace transform is

```text
G_omega(s)
 =zeta(s-omega)/(s*zeta(s+omega))
  -A_omega/(s-1-omega).                              (JMCE.10)
```

The subtraction removes the residue at `s=1+omega`. Initially this is an
ordinary transform for `Re(s)>1-omega`, using the elementary error bound.

If `J_omega<infinity`, Paley-Wiener puts (JMCE.10) in `H2` of
`Re(s)>1/2`. Every possible pole from a zero of `zeta(s+omega)` must then
be removable, so Suzuki's reduced fixed-shift condition `D(omega)` holds.

Conversely, if `D(omega)` holds, the completed quotient

```text
Q_omega(s-1/2)
 =xi(s-omega)/xi(s+omega)
```

is inner. Factoring (JMCE.10) through `Q_omega`, removing its single
archimedean pole at `1+omega`, and using the standard gamma quotient
estimate `O(|t|^(omega-1))` gives `G_omega in H2`. Transform uniqueness
then recovers the actual arithmetic error. The resulting fixed-shift
theorem candidate is

```text
J_omega<infinity
iff D(omega).                                        (JMCE.11)
```

Combining this with the audited cofinal phase theorem gives

```text
RH
iff there are omega_j->0 such that
    integral_1^infinity |E_(omega_j)(x)|^2 dx/x^2
    is finite for every j.                           (JMCE.12)
```

Burnol's shifted-ratio estimate and the Balazard-Saias natural-Mobius
approximation show under RH that `E_(omega,N)->E_omega` in `H` for every
fixed `omega`. In the reverse direction, (JMCE.8) and (JMCE.12) require
only a uniform norm bound. This yields the finite natural criterion

```text
RH
iff there are omega_j->0 such that
    sup_N ||E_(omega_j,N)||_H < infinity
    for every j.                                     (JMCE.13)
```

Equations (JMCE.11)-(JMCE.13) are new corpus theorem candidates requiring
independent review. They do not assert the required bound.

## Finite Brownian Energy

For every finite `X>=1`, direct expansion inside the finite integral gives

```text
J_omega(X)
 =integral_1^X E_omega(x)^2 dx/x^2

 =sum_(n,m<=X)c_omega(n)c_omega(m)
   *(1/max(n,m)-1/X)

  -(2*A_omega/omega)
   *sum_(n<=X)c_omega(n)*(X^omega-n^omega)

  +A_omega^2*(X^(1+2omega)-1)/(1+2omega).            (JMCE.14)
```

There is also a clean signed-measure form. Put

```text
d eta_(omega,X)(u)
 =sum_(n<=X)c_omega(n)delta_n
  -A_omega*(1+omega)u^omega*1_(0,X](u)du.
```

Then

```text
J_omega(X)
 =double_integral_[0,X]^2
  K_X(u,v)deta(u)deta(v),

K_X(u,v)
 =integral_1^X 1_(u<=x)1_(v<=x)dx/x^2
 =1/max(1,u,v)-1/X.                                  (JMCE.15)
```

The kernel is positive semidefinite. Under reciprocal coordinates it is
a Brownian `min` kernel, and every ordered finite matrix is totally
nonnegative by a cumulative-matrix factorization. This is exact contact
with total positivity, but it supplies only `J_omega(X)>=0`. The
positive-positive, cross, and continuum-continuum terms in (JMCE.14) each
have the same potentially divergent scale; a uniform upper bound requires
their arithmetic cancellation.

## Countermodel Gates

Boundary energy is automatic. For `a>0`,

```text
1/(a+z),  1/(a-z)
```

have the same boundary squared norm `1/(2a)`. The first is the Laplace
transform of the causal `exp(-a*t)1_(t>0)`. The second has anti-causal
boundary inverse `exp(a*t)1_(t<0)`, while its actual positive-time
Bromwich inverse is `-exp(a*t)` and is not square integrable. Thus
Plancherel without Hardy support cannot prove (JMCE.11).

Likewise, scalar coefficients in `l2` do not make an arbitrary unit-vector
series convergent: `v_d=v` and coefficients `1/d` are the elementary
countermodel. A Bessel or cancellation estimate for the actual dilation
orbit is a theorem obligation, not a Hilbert-space formality.

## Literature Fit

Burnol's `f_epsilon` theorem has the same shifted-zeta-ratio,
square-integrability, natural-Mobius, and causality architecture. The
present `R_omega` is the fractional-power Jordan remainder naturally
selected by Suzuki's coefficients.

Baez-Duarte's general strong Nyman-Beurling theorem does not close
(JMCE.13): its audited sufficiency direction requires compact support,
whereas `R_omega` is not compactly supported, and arbitrary closure does
not prove convergence of this fixed natural Mobius sequence.

Primary sources:

- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- Jean-Francois Burnol, `On an analytic estimate in the theory of the Riemann Zeta function and a Theorem of Baez-Duarte`: https://arxiv.org/abs/math/0202166
- Luis Baez-Duarte, `A general strong Nyman-Beurling Criterion for the Riemann Hypothesis`: https://arxiv.org/abs/math/0505453

## Open Gate

For one explicit cofinal sequence, prove the all-height finite inequality

```text
sup_N integral_R W_omega(t)
 |zeta(1/2+omega+i*t)M_N(1/2+omega+i*t)|^2 dt
 <infinity.                                          (JMCE.16)
```

This is a sharper target than another finite `x` or finite `t` grid. It is
also exactly RH-strength through (JMCE.13), so no claim of an easy
mean-square estimate is made.

## Proof Boundary

Equations (JMCE.1)-(JMCE.10), (JMCE.14), and (JMCE.15) are exact
identities, convergence statements, or finite Hilbert-space formulas.
Equations (JMCE.11)-(JMCE.13) are internally audited theorem candidates
using standard Hardy-space and conditional natural-approximation
machinery. The uniform estimate (JMCE.16) is open. This artifact does not
prove that estimate, RH, or `Lambda <= 0`.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jordan-Muntz causal-energy bridge: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['theorem_candidate_count']} theorem candidates, "
        f"{payload['audit']['guard_count']} guards, "
        f"{payload['audit']['open_arithmetic_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
