# Jensen-Window PF Edrei Hankel Boundary-Flux Gate

Date: 2026-07-24

Status: exact finite-rank Hankel boundary flux, infinite-support
canonical-order escape, and high-shift collision sensor with one open
Xi/Phi determinant-growth gate. This is not a proof of PF-infinity,
Jensen hyperbolicity for zeta, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_edrei_hankel_boundary_flux_gate`.

```text
work/rh_compute/results/jensen_window_pf_edrei_hankel_boundary_flux_gate.json
python work/rh_compute/scripts/jensen_window_pf_edrei_hankel_boundary_flux_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_edrei_hankel_boundary_flux_gate.py
```

## Polynomial Flux Identity

Retain the exact moment hierarchy from Lemma 11.22C:

```text
a_n'=-2(n+1)(2n+3)a_(n+1)+4(n+1)sum_(k=0)^n a_k*a_(n-k)
```

At a time where the logarithmic derivative has a Stieltjes measure
`nu`, write `P(t)=sum_i c_i t^i`, put `E=t*d/dt`, and define

```text
Q_s(P)=sum_(i,j)c_i*c_j*a_(s+i+j)=integral x^s*P(x)^2 dnu(x)
f_s(t)=t^s*P(t)^2
g_s(t)=t*f_s(t)
h_s(t)=t*g_s'(t)
```

Coefficient extraction and

```text
(n+1)sum_(k=0)^n x^k*y^(n-k)
 =[x*d/dx(x^(n+1))-y*d/dy(y^(n+1))]/(x-y)
```

give the exact divided-difference identity

```text
Q_s'(P)=-2 integral x(E+1)(2E+3)[x^s P^2]dnu+4 double_integral [h_s(x)-h_s(y)]/(x-y)dnu(x)dnu(y), h_s(t)=t*d/dt[t^(s+1)P(t)^2]
```

with the diagonal value interpreted as `h_s'(x)`. This is an identity
for every polynomial `P`; no sign has been imposed.

## Finite-Rank Edrei Boundary

Sokal's logarithmic-derivative criterion supplies the quantized measure

```text
For F in LP+: dnu=gamma*delta_0+sum_j rho_j*delta_(beta_j), rho_j=m_j*beta_j, m_j positive integers
```

Suppose `s>=1` and `P` vanishes on every positive atom. Off-diagonal
divided differences vanish. At one atom `beta`, direct differentiation
gives

```text
-2*beta*(E+1)*(2E+3)[t^s P(t)^2]_(t=beta)
 =-8*beta^(s+3)*P'(beta)^2
h_s'(beta)=2*beta^(s+2)*P'(beta)^2.
```

The linear and diagonal quadratic pieces therefore combine to

```text
If s>=1 and P(beta_j)=0 on every positive atom, then Q_s'(P)=8 sum_j rho_j*beta_j^(s+2)*(rho_j-beta_j)*P'(beta_j)^2
For Edrei residues rho_j=m_j*beta_j: Q_s'(P)=8 sum_j m_j*(m_j-1)*beta_j^(s+4)*P'(beta_j)^2>=0
```

This proves forward nonnegativity at every finite-rank Edrei wall.
It is strict precisely when a factor is repeated; simple factors are
tangent at first order.

For `r` positive atoms and the monic annihilator
`P(t)=product_j(t-beta_j)`, the shifted moment matrix has rank `r`.
Its adjugate is the outer product of the coefficient vector of `P`
times its leading Vandermonde cofactor. Hence

```text
For r positive atoms, monic P=product_j(t-beta_j), and D_(r,s)=det(a_(s+i+j))_(i,j=0)^r: D_(r,s)'=tau_(r,s)*Q_s'(P), tau_(r,s)=product_j(rho_j*beta_j^s)*product_(i<j)(beta_i-beta_j)^2
```

At `r=1` this recovers

```text
D_(1,s)'=8*m^2*(m-1)*beta^(2s+5).
```

The builder checks 16 exact rational null-form
instances and 12 exact rational determinant
instances.

## Residue-Quantization Guard

The sign is not a generic consequence of Stieltjes positivity:

```text
For nu=(1/2)delta_1, P(t)=t-1, and s=1, Q_s'(P)=-2: Stieltjes positivity without the integer Edrei residue condition does not orient the wall forward
```

The formal logarithmic primitive has a fractional power and is not
entire. The integer Edrei residue, not positivity of `nu` alone, is the
structural input behind the forward orientation.

## Infinite-Support Boundary And Escaping Order

For an LP+ endpoint with distinct atoms
`beta_1>beta_2>...>0` and `s>=1`, Cauchy-Binet gives

```text
D_(d,s)=sum_(|J|=d+1) product_(j in J)(rho_j*beta_j^s)
          *V(beta_J)^2>0.
```

For `s=0`, the same positive-atom sum is strictly positive and the
optional `gamma*delta_0` atom contributes only additional nonnegative
Cauchy-Binet terms. Thus every fixed block in both canonical shifts
`s=0,1` is positive definite.

Thus every fixed finite shifted block is positive definite when the
support is infinite. If a continuous family approaches such a boundary
from outside LP+, each fixed canonical block at shifts zero and one
stays positive in a neighborhood of the boundary. Sokal's criterion
then proves

```text
If a continuous family exits LP+ through an infinite-support boundary, every outside sequence approaching the boundary requires canonical Stieltjes witness order d->infinity in one of the shifts s=0 or s=1
```

This is stronger than saying that a large finite computation has not
yet found the right row: no fixed canonical order can be the boundary
detector.

## High-Shift Collision Sensor

Assume the first repeated atom is `beta_r`, so
`m_1=...=m_(r-1)=1` and `m_r>=2`. Let `pi_(r,s)` be the monic
degree-`r` orthogonal polynomial for `x^s*dnu`, and put
`h_(r,s)=Q_s(pi_(r,s))`. Schur complementation and differentiation
at fixed leading coefficient give the exact identities

```text
For the monic x^s*dnu-orthogonal polynomial pi_(r,s): h_(r,s)=Q_s(pi_(r,s))=D_(r,s)/D_(r-1,s) and partial_lambda log D_(r,s)=partial_lambda log D_(r-1,s)+Q_s'(pi_(r,s))/Q_s(pi_(r,s))
```

Indeed, `partial_lambda pi_(r,s)` has degree at most `r-1`, so its
cross term with `pi_(r,s)` vanishes by orthogonality. Cauchy-Binet
then gives

```text
D_(r,s)=C_D*(beta_1*...*beta_r*beta_(r+1))^s*(1+o(1)).
D_(r-1,s)=C_0*(beta_1*...*beta_r)^s*(1+o(1)).
h_(r,s)=C_H*beta_(r+1)^s*(1+o(1)),
```

To control the numerator, use the Lagrange polynomial `L_i` on the
first `r` atoms in the orthogonality equation. Since `a_0<infinity`,
the tail is uniformly summable and

```text
If P_r(t)=product_(j=1)^r(t-beta_j), then pi_(r,s)(beta_i)=O((beta_(r+1)/beta_i)^s) for i<=r and pi_(r,s)=P_r+O((beta_(r+1)/beta_r)^s) coefficientwise
```

Inserting this into the exact divided-difference flux shows that every
simple leading atom and every tail term is
`O(s^2*beta_(r+1)^s)`, while the repeated atom contributes
`8*m_r*(m_r-1)*beta_r^(s+4)*P_r'(beta_r)^2*(1+o(1))`.
The lower logarithmic derivative is only `O(s^2)`, because its
leading `r` distinct atoms already give full rank. Consequently

```text
partial_lambda D_(r,s)
 =C_N*(beta_1*...*beta_(r-1)*beta_r^2)^s*(1+o(1)),
If beta_1>beta_2>... and r is the first index with m_r>=2, then lim_(s->infinity)(partial_lambda log D_(r,s))^(1/s)=beta_r/beta_(r+1)>1
```

where `C_D,C_N>0`. In the exact two-atom benchmark
`nu=m*delta_1+n*q*delta_q, 0<q<1`, the simple leading atom gives

```text
-2*(q + 1)*(q**2*s + 3*q**2 - 2*q*s - 14*q + s + 3)/(q - 1)**2
```

whose `s`th-root growth is one. Replacing the leading atom by a double
factor and taking `q=1/2` gives

```text
128*2**s + 4*s**2 + 29*s + 131
```

whose `s`th-root growth is two, exactly `beta_1/beta_2`.

## Newman Reduction

If `Lambda<=0`, Lemma 11.9 says every zero is simple for every
`0<lambda<=1/5`. For each fixed `r`, the same orthogonal-quotient
argument now has no repeated atom among the annihilated leading support:
both `Q_s'(pi_(r,s))` and `Q_s(pi_(r,s))` have the
`beta_(r+1)^s` scale up to polynomial factors. Recursing through the
lower quotients gives `partial_lambda log D_(r,s)=O(s^2)`, so the
subexponential condition below holds.

Conversely, the audited positive-boundary theorem says that `Lambda>0`
produces a finite multiple real zero of `H_Lambda` with
`0<Lambda<=1/5`. The transform has infinitely many distinct zero
locations, so its Edrei measure has a next atom. The first repeated
descending reciprocal-square atom then violates the condition by the
collision sensor. Therefore the following criterion is exactly
equivalent to `Lambda<=0`:

```text
Prove directly from the Phi transform that for every fixed r and every lambda in (0,1/5] for which H_lambda has only real zeros, limsup_(s->infinity)max(1,partial_lambda log D_(r,s))^(1/s)<=1
```

This equivalence does not prove its open direction. A proof must use the
actual Phi family. The repeated quadratic from
Lemma 11.22C violates this bound, so it cannot follow from LP+,
Stieltjes positivity, radial heat, or positive scale mixing alone.

## Literature Audit

Sokal's Proposition 6 supplies the discrete integer-residue measure.
Standard Gram/Cauchy-Binet identities supply the positive Vandermonde
expansion. The searched Hankel/Toda and orthogonal-polynomial machinery
organizes moment determinants but does not provide the required
Phi-specific subexponential heat-flux estimate. The Polymath positive-time
asymptotics supply infinite high zeros and boundary localization, not
the missing determinant bound.

## Proof Boundary

This artifact proves the polynomial Hankel-flux identity, the finite-rank integer-residue orientation and determinant formula, strictness and canonical-order escape at an infinite-support LP+ boundary, and the high-shift determinant signature forced by a repeated atom. It does not prove the required Xi/Phi subexponential determinant-growth estimate, all-order Stieltjes positivity at lambda zero, PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda<=0.
