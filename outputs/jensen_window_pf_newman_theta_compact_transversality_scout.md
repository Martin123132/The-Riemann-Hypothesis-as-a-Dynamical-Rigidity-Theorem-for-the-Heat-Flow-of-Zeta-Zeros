# Newman Theta Compact-Transversality Scout

Date: 2026-07-24

Status: one exact near-origin theorem, one compact numerical
scout, and an exact proof-partition handoff. This is not a proof
of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_compact_transversality_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_compact_transversality_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_compact_transversality_scout.py
```

Current result:

```text
validated Newman theta compact-transversality scout: 12 rows, 0 issues, 3 exact moment/kernel inequalities, 1 uniform near-origin no-contact theorem, 61951 compact diagnostic points, 0 compact grid failures, 1 independent quadrature crosscheck, 1 outer nonpromotion guard, 2 high-frequency composition theorems, 2 open certification targets
```

## Exact Near-Origin Theorem

phi_t(u)=exp(t*u^2)*Phi(u), H_t(x)=integral_0^infinity phi_t(u)*cos(xu)du

phi_n(u)>0 for n>=1 and u>=0, hence phi_t(u)>0 for 0<=t<=1/5.

Put `M_k(t)=integral u^k*phi_t(u)du`. Rational exponential
bounds give

```text
M_0(t)=integral phi_t(u)du>9/2900 for 0<=t<=1/5
M_2(t)<40/513 for 0<=t<=1/5
H_t(x)>=M_0(t)-x^2*M_2(t)/2>9/2900-5/2052
```

Since the final rational margin is `248/371925>0`,

```text
0<=t<=1/5 and |x|<=1/4
H_t(x)>0 on the stated rectangle, so H_t and H_t' cannot vanish together there.
```

This removes the artificial degeneracy caused by multiplying
the contact functional by `x^4`.

## Compact First-Block Scout

The independently checked pointwise implication is

```text
|J_1|>B_0 or |J_1'|>B_1 implies (J,J')!=(0,0).
```

On the sampled rectangle

```text
0<=t<=1/5 in steps of 0.005, 1/4<=x<=38 in steps of 0.025
points=61951
minimum max(|J_1|/B_0,|J_1'|/B_1)=1.849415879849908
minimum location: t=0.2, x=38.0
grid failures=0
maximum 600/900-node ratio difference=4.129106e-10
```

These are point values. They do not control the continuum between
points and therefore are not promoted to a compact theorem.

## Scope Guard

The raw first-block disjunction loses its margin near x=39.6 and cannot be promoted to an all-frequency theorem.

```text
outer diagnostic failures=2042/2205
t=0.0: first failure x=39.65, failure points=408
t=0.05: first failure x=39.65, failure points=408
t=0.1: first failure x=39.65, failure points=408
t=0.15: first failure x=39.6, failure points=409
t=0.2: first failure x=39.6, failure points=409
```

This confirms that raw moment dominance is a compact tool only.

## Global Proof Partition

The existing Polymath-15 certificates already supply the correct
cancellation-preserving high-frequency framework:

```text
L=log(x/(4*pi)), c=t*L
L>=50 and c>=25: exact dominant-saddle ray closed
for every epsilon>0, c>=c_*+epsilon: asymptotically closed
c_*=4911678521/1933561194
L>=50 and 0<c<=c_*+o(1): corrected phase transversality open
bounded/intervening L: compact certification still open
```

Endpoint subtraction should not replace this phase-aware machinery:
absolute decay of both `H_t` and `H_t'` cannot itself exclude their
simultaneous vanishing.

## Live Handoff

First promote the robust sampled rectangle to a rigorous
two-parameter interval certificate using Arb quadrature or a
Taylor/Cauchy enclosure. In parallel, retain the corrected
Riemann-Siegel phase target on the residual scaled layer. Neither
task may assume zero simplicity, `Lambda<=0`, or RH.
