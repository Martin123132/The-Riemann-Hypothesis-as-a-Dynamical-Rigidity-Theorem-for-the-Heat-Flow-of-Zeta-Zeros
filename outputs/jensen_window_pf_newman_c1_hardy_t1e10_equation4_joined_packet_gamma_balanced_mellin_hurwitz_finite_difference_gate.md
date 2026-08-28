# Gamma-balanced Mellin-Hurwitz coordinate for the joined packet

Date: 2026-08-27

Status: exact scaled coordinate and explicit absolute-majorant barrier certified; saddle-aligned enclosure open.

Put

```text
omega=exp(i*pi/4),
q_-=79788+1/2,       q_+=2561211+1/2,
Delta_q(u)=[exp(-q_-u)-exp(-q_+u)]/[1-exp(-u)].
```

The denominator is removable at zero and `Delta_q(0)=2481423`.  In the
decaying half-line source use `u=2*pi*omega*x`.  The Gaussian becomes
`exp(i*u^2/(4*pi))`.  Rotating the resulting ray from `arg(u)=pi/4` to the
positive real axis is legal: the Gaussian does not grow in that sector, the
finite difference is regular at zero, and its real part decays at infinity.
Therefore

```text
M_s(q_-,q_+)
 =integral_0^infinity u^(-s)exp(i*u^2/(4*pi))Delta_q(u)du,

S_W=E_s M_s(q_-,q_+),
E_s=exp(pi*t/2+i*pi/4)(2*pi)^(s-1).                 (MH1)
```

Removing only the Gaussian chirp gives the exact finite Hurwitz anchor

```text
integral_0^infinity u^(-s)Delta_q(u)du
 =Gamma(1-s)[zeta(1-s,q_-)-zeta(1-s,q_+)].          (MH2)
```

The apparently large factor in (MH1) is exactly gamma-balanced:

```text
|E_s Gamma(1-s)|^2=1/[1+exp(-2*pi*t)].              (MH3)
```

Thus no exponential loss is intrinsic to this coordinate.

For every integer `K>=1`, Taylor's formula gives the exact finite expansion

```text
M_s=sum_(k=0)^(K-1) (i/(4*pi))^k/k!
       Gamma(1-s+2k)Delta_zeta(1-s+2k)+R_K,         (MH4)

Delta_zeta(w)=zeta(w,q_-)-zeta(w,q_+),

R_K=integral_0^infinity u^(-s)Delta_q(u)
    [exp(ix)-sum_(k<K)(ix)^k/k!]du,
x=u^2/(4*pi).                                       (MH5)
```

The elementary real-axis remainder estimate is explicit:

```text
|E_s R_K| <= |E_s| Gamma(2K+1/2)Delta_zeta(2K+1/2)
              /[(4*pi)^K K!].                       (MH6)
```

It is also unusable at the production height.  If `B_K` denotes the right
side of (MH6), positivity of the finite Hurwitz sum gives

```text
B_(K+1)/B_K
 <= [K+3/(16(K+1))]/(pi*q_-^2).                    (MH7)
```

Exact rational arithmetic with a 99-decimal enclosing interval for `pi`
proves that the right side of (MH7) is below one through

```text
K=20,000,022,018
```

and exceeds one at the next integer.  Consequently this particular absolute
Taylor majorant keeps decreasing for more than twenty billion orders before
it can even reach a turning point.  It is not a viable evaluator.  The useful
part is (MH1)--(MH4): the remainder must remain oscillatory and be moved onto
a saddle-aligned contour or an endpoint-complete Mordell representation.

An altered `t=3.25`, five-label row compares the original decaying half-line,
the scaled Mellin integral, the Hurwitz anchor, and the order-three finite
expansion.  The independent checker changes the height, roster, order,
quadrature partition, and precision.

Primary-source provenance: the entire auxiliary-kernel contour move follows
the no-residue mechanism in Arias de Reyna, https://arxiv.org/abs/2407.02016.
The unequal-truncation/Mordell route context is O'Sullivan,
https://arxiv.org/abs/1811.01130.  No implicit asymptotic constant is used.

Pi provenance: every `pi` comes from the inherited Riemann-Siegel Gaussian,
the fixed quarter-turn, the Mellin scaling, or the standard gamma identity.

## Proof boundary

Exact gamma-balanced Mellin-Hurwitz finite-difference coordinate, finite expansion with explicit remainder, and a rigorous non-scalability result for its positive-axis absolute Taylor majorant only. No saddle-aligned remainder bound, quantitative joined-packet, J_Z, or D_K enclosure, and no non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
