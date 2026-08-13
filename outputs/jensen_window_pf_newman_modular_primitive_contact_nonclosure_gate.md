# Newman Modular-Primitive Contact Nonclosure Gate

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
