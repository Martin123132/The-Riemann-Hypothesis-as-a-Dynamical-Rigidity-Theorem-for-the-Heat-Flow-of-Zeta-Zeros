# W4 contour recombination and complete endpoint Hurwitz majorants

Date: 2026-08-07

Status: rigorous symbolic selector-gap-uniform local endpoint theorem with an all-374 finite-roster application; not a proof of the recurrence budget or RH

## W4 is an endpoint half-ray difference

For `Phi1>0`, deform the finite zero-mode integral to the two lower
half-rays on which `exp(-2*pi*i*F)` decays.  Cauchy's theorem gives

```text
zero_mode=H_exact(0,0)-H_exact(N,0).                  (1)
```

The quadratic half-ray in W4 is `H_quad(0,0)`, while

```text
H_linear(N,0)=-i*exp(-2*pi*i*F(N))/(2*pi*F'(N)).     (2)
```

Hence

```text
zero_mode-W4
 =[H_exact(0,0)-H_quad(0,0)]
  -[H_exact(N,0)-H_linear(N,0)].                     (3)
```

The explicit W4 endpoint term is therefore the missing linear `n=0` ray,
not an independent error.  Arb verifies (1)--(3) at 70 and 110 digits on all
209 positive-`Phi1` calls.  The largest
identity gap is below `5.90081052919642030898517202213198110104970123942090E-41`.

## Complete endpoint majorant

For an exact-minus-linear mode with starting gap `q>0`, local quadratic
coefficient `B>0`, and `C=|Phi3|`, a common-ray homotopy and the lower bound
`sin(theta)>=1/2` give

```text
sum_(k>=0)|D_(q+k)|
 <= 4*B*zeta(3,q)/pi^2 + 12*C*zeta(4,q)/pi^3.        (4)
```

The upper family uses (4) at `q=m-Phi1` at zero and `q=m+1-xi` at `N`,
plus the W2 bound of Section 11.271.  The lower W3 family starts at
`q=2+Phi1` and `q=1+xi`; the W4 family starts at `q=1+Phi1` and `q=xi`,
plus the common quadratic zero-mode bound.  Thus every Hurwitz argument is
bounded away from the exceptional selector gap that was removed before
estimation.

Arb certifies the resulting upper, lower, and complete local bounds on all
374 calls:

```text
maximum upper-family bound   <= 8.50672090967257588797318287010218076920024411023214E-3
maximum lower-family bound   <= 4.54917938885938990371468578204369843342550290050494E-3
maximum complete local bound <= 1.06817811702884647221191642133974934954252917944494E-2
worst complete bound/actual ratio <= 2.78424234420720423802933313625611960816939856354114E+3
```

Pi in (1)--(4) comes from the Fourier phase `exp(2*pi*i*n*u)` and the exact
Gaussian/exponential moments; it is not fitted.  The theorem is uniform in
the selector gaps but conservative.  It does not yet prove that its bound
fits the complete recurrence budget, control coefficient transport, Legendre
tails, cross-block scaling, the outer Hardy remainder, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
