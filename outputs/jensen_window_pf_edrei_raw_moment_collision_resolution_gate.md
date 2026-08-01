# Jensen-Window PF Edrei Raw-Moment Collision-Resolution Gate

Date: 2026-07-24

Status: exact Phi-moment transfer and LP+ collision-resolution
obstruction with one open Xi/Phi precision handoff. This is not a
proof of PF-infinity, Jensen hyperbolicity for zeta, RH, or
`Lambda <= 0`.

Artifact kind:
`jensen_window_pf_edrei_raw_moment_collision_resolution_gate`.

```text
work/rh_compute/results/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.json
python work/rh_compute/scripts/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py
```

## Exact Phi-Moment Transfer

Put

```text
M_n(lambda)=integral_R exp(lambda*u^2)*Phi(u)*u^(2n)du
F_lambda(z)=sum_(n>=0)M_n(lambda)*z^n/(2n)!
H_lambda(z)=sum_(n>=0)d_n z^n, d_n=M_n/((2n)!*M_0)
```

For
`H_lambda'/H_lambda=sum_(n>=0)(-1)^n*a_n*z^n`,
coefficient comparison in `H'=H*(H'/H)` gives

```text
(n+1)d_(n+1)=sum_(k=0)^n(-1)^k*a_k*d_(n-k); a_n=(-1)^n[(n+1)d_(n+1)-sum_(k=0)^(n-1)(-1)^k*a_k*d_(n-k)]
```

Thus `D_(r,s)` depends on the raw moments through at least order
`s+2r+1`; it is not a local fixed-order functional of a saddle
expansion. The builder performs
104 exact rational triangular checks.

## Symmetric Collision Family

Fix `0<q<beta` and `0<epsilon<beta-q`. Consider

```text
H_0(z)=(1+beta*z)^2(1+q*z), H_epsilon(z)=(1+(beta+epsilon)z)*(1+(beta-epsilon)z)*(1+q*z)
```

Both functions are normalized LP+ polynomials. Direct multiplication
gives the exact compact-open gap

```text
H_epsilon-H_0=-epsilon^2*z^2*(1+q*z)
```

Multiplying both functions by one common normalized LP+ tail whose
atoms lie below `q` gives the same identity times that tail, so the
obstruction is not a finite-support artifact. Their Edrei moments
differ by

```text
a_n(epsilon)-a_n(0)=(beta+epsilon)^(n+1)+(beta-epsilon)^(n+1)-2*beta^(n+1)
 =2*sum_(j>=1)binom(n+1,2j)
    *beta^(n+1-2j)*epsilon^(2j).
```

The gap is quadratic for each fixed `n`, but not uniformly small at
the high orders needed by the shifted determinant.

## Exponential Cancellation Condition Number

At `epsilon=0`, the distinct Edrei atoms have residues `2*beta` and
`q`. Cauchy-Binet gives

```text
D_(1,s)(0)=2*(beta*q)^(s+1)*(beta-q)^2
```

The two products subtracted in
`D_(1,s)=a_s*a_(s+2)-a_(s+1)^2`
are each of order `beta^(2s)`. With the natural componentwise
subtraction condition number

```text
kappa_s=(a_s*a_(s+2)+a_(s+1)^2)/D_(1,s)~[4*beta^2/(beta-q)^2]*(beta/q)^(s+1)
```

a worst-case entrywise relative remainder must therefore be
`o((q/beta)^s)` to resolve the collision determinant uniformly.
Any fixed algebraic remainder `O(s^-N)` loses this comparison.
For `beta=1,q=1/2`, the exact audit confirms
`kappa_s/2^s -> 32 for beta=1, q=1/2`.

This is a black-box reconstruction barrier. A correlated integral
identity, sum-of-squares formula, or operator inequality that preserves
the cancellation can bypass it; mere independent error bars cannot.

## Split-Pair Crossover

For `epsilon>0`, Cauchy-Binet gives

```text
D_(1,s)(epsilon)=4*epsilon^2*(beta^2-epsilon^2)^(s+1)+[((beta+epsilon)*q)^(s+1)*(beta+epsilon-q)^2]+[((beta-epsilon)*q)^(s+1)*(beta-epsilon-q)^2]
```

The first term has the largest exponential base but coefficient
`4*epsilon^2`. Relative to the collision determinant,

```text
T_split/D_collision=[2*epsilon^2/(beta-q)^2]*[(beta^2-epsilon^2)/(beta*q)]^(s+1)
```

so its benchmark crossover is

```text
s_cross=2*log(beta/epsilon)/log(beta/q)+O(1)
```

The builder checks 4 decreasing gaps.
Every finite shift range can therefore be made collision-like by
taking a sufficiently small but nonzero split.

### Arbitrary Fixed Rank

Let `B={b_1,...,b_(r-1)}` be any fixed string of larger simple
atoms, `b_1>...>b_(r-1)>beta>q`, and put
`P_B=product_(b in B)b`. Multiplying both collision polynomials
by `product_(b in B)(1+b*z)` moves the repeated factor to rank
`r`. Cauchy-Binet gives

```text
For B={b_1,...,b_(r-1)} above beta>q and P_B=product B: D_(r,s)(0)=2*(P_B*beta*q)^(s+1)*V(B,beta,q)^2
T_top(epsilon)=[P_B*(beta^2-epsilon^2)]^(s+1)*V(B,beta+epsilon,beta-epsilon)^2
```

The top split term divided by the collision determinant is

```text
(1/2)*[(beta^2-epsilon^2)/(beta*q)]^(s+1)
 *V(B,beta+epsilon,beta-epsilon)^2/V(B,beta,q)^2.
```

As `epsilon->0`, the Vandermonde ratio is

```text
4*epsilon^2/(beta-q)^2
 *product_(b in B)[(b-beta)/(b-q)]^2*(1+O(epsilon^2)).
```

The prefix therefore changes only the bounded prefactor. The
crossover remains
`2*log(beta/epsilon)/log(beta/q)+O(1)` at every fixed rank.
The builder checks 3 ranks,
24 exact collision determinants,
24 split Cauchy-Binet
identities, and 24 exact top-term
ratios.

## Noncommuting Flux Limits

The exact repeated-model logarithmic flux is

```text
With x=q/beta, partial_lambda log D_(1,s)(0)=beta*[16*x^(-s-1)-P_s(x)]/(1-x)^2, P_s(x)=2*x^3*s+6*x^3-4*x^2*s^2-18*x^2*s-50*x^2+8*x*s^2+38*x*s+18*x-4*s^2-22*s-30
```

For `beta=1,q=1/2` this reduces to

```text
partial_lambda log D_(1,s)
 =128*2^s+4*s^2+29*s+131.
```

For every fixed `epsilon>0`, the atoms are distinct. Differentiating
their positive Cauchy-Binet pair sum under radial heat shows that the
largest pair contributes an affine function of `s`, while all other
pairs are exponentially smaller. Hence

```text
partial_lambda log D_(1,s)(epsilon)
 =A_epsilon*s+B_epsilon+o(1),
I_r(H)=limsup_(s->infinity)max(1,partial_lambda log D_(r,s)(H))^(1/s); after any fixed larger simple prefix, I_r(H_epsilon)=1 for epsilon>0 and I_r(H_0)=beta/q>1
```

For each fixed `s`, the coefficient and heat-flow formulas are
polynomial in `epsilon^2`, so the split flux tends to the repeated
flux. The limits `epsilon->0` and `s->infinity` do not commute.
The builder checks 36 exact split
determinants, 75 exact velocity/flux
identities, and 36 repeated benchmarks.

## Proof Consequence

This artifact proves:

- the exact `Phi raw moments -> Taylor coefficients -> Edrei moments`
  transfer;
- an exponentially ill-conditioned collision determinant;
- a logarithmically escaping split-resolution shift;
- a compact-open LP+ family with noncommuting collision limits.

It does not prove:

- the Xi/Phi subexponential determinant-flux estimate;
- all-order Stieltjes positivity at `lambda=0`;
- PF-infinity, Jensen hyperbolicity for zeta, RH, or `Lambda<=0`.

The surviving handoff is precise:

```text
Control the Xi/Phi determinant flux by a structural identity that preserves the rank-one cancellation, or obtain gap-adapted exponentially small remainders; fixed algebraic raw-moment asymptotics alone do not resolve the collision.
```

The existing cancellation-preserving Phi formulation is the strict
Laguerre correlation target

```text
L_t(x)=H_t'(x)^2-H_t(x)H_t''(x)
      =Fourier[K_(1,t)](2x).
```

It is recorded in
`outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md`.
The present gate explains why returning to that correlated kernel is
structurally preferable to extending entrywise raw-moment expansions;
it does not prove `L_t(x)>0`.

A finite Poincare saddle expansion with an algebraic remainder is
not, by itself, strong enough at a possible collision.
