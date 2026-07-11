# Jensen-Window PF Newman Gasper Residual Two-Block Gate

Date: 2026-07-11

Status: exact and interval-certified positive two-block Gasper
obstruction. This is not a proof or disproof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_newman_gasper_residual_two_block_gate`.

```text
work/rh_compute/results/jensen_window_pf_newman_gasper_residual_two_block_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_gasper_residual_two_block_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_gasper_residual_two_block_gate.py
```

Current result:

```text
validated Jensen-window PF Newman Gasper residual two-block gate: 8 rows, 0 issues, 2 exact kernel theorems, 1 exact Laguerre quadratic, 2 Acb derivative certificates, 2 beta intervals covered, 1 exhaustive positive-residual obstruction, 1 multiplier guard, 1 signed handoff
```

## Positive First Residual

For `q=exp(-4u)` and `a=3/(2pi)`,

```text
With q=exp(-4u) and a=3/(2*pi), phi_1(u)/Psi_9(u)=exp(pi*q)*(1-a*q)/(1+q^(9/2)).
exp(pi*q)*(1-a*q)>(1+pi*q)*(1-a*q)>1+q>=1+q^(9/2) for 0<q<=1, since pi>3 gives pi-a-3/2>1.
Phi(u)>phi_1(u)>Psi_9(u) for every finite u>=0
```

Thus subtracting the original fake-Xi block leaves a strictly positive
kernel. Positivity alone is not enough: the interval theorem below
includes `beta=0` and proves that this residual transform violates the
first Laguerre inequality.

## Two-Block Family

```text
Psi_beta(u)=4*pi^2*exp(-2*pi*cosh(4u))*(cosh(9u)+beta*cosh(5u)), beta>=0
P_beta=P_(9/4)+beta*P_(5/4), R_beta=H_0-P_beta
Phi(u)-Psi_beta(u)~2*pi^2*(pi-3/(2*pi)-beta)*exp(5u-pi*exp(4u))
```

A globally nonnegative residual must therefore have
`0<=beta<=pi-3/(2pi)`. On this complete parameter interval,

```text
L[R_beta]=L[E_0]-beta*B[E_0,P_(5/4)]+beta^2*L[P_(5/4)]
```

## Interval Theorem

For real entire f, central differences with step h enclose f' and f''. A square Acb evaluation bounds abs(f) by M on every Cauchy circle of radius r around the finite-difference segment; the added errors are h^2*M/r^3 and 2*h^2*M/r^4.
The resulting certified quadratics are:

```text
x=66, beta in [0, [2.3000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 4.66e-111]]
  c0=[-1.483e-21 +/- 5.70e-25]
  c1=[-6.213e-21 +/- 4.82e-25]
  c2=[2.8475e-21 +/- 3.10e-26]
  left=[-1.483e-21 +/- 5.70e-25]
  right=[-7.1e-22 +/- 1.78e-24]
x=50, beta in [[2.3000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 4.66e-111], [2.6641278243141072311559920931619597980937904621537364747319425601311260133835293933576934170433578101139295583 +/- 3.89e-110]]
  c0=[1.411e-15 +/- 1.31e-19]
  c1=[-1.4587e-15 +/- 6.49e-20]
  c2=[3.4559e-16 +/- 4.59e-21]
  left=[-1.16e-16 +/- 3.95e-19]
  right=[-2.2e-17 +/- 5.16e-19]
```

Each quadratic is strictly convex and negative at both endpoints,
hence negative throughout its closed beta interval. The intervals
meet at `23/10` and cover every beta allowed by tail positivity.
Consequently no member of this positive two-block family leaves a
positive residual whose transform is independently Laguerre-Polya.

## Multiplier Guard

The tail-matched multiplier also fails the standard universal-factor
test. Writing `U_beta(z)=cosh(9z)+beta*cosh(5z)` gives

```text
U_beta(z)=cosh(z)*Q_beta(cosh(z)^2), where Q_beta(y)=5*beta + 256*y**4 - 576*y**3 + y**2*(16*beta + 432) - y*(20*beta + 120) + 9
disc_y Q_beta=-16777216*(20*beta**5 - 64*beta**4 + 184*beta**3 - 432*beta**2 + 972*beta - 729)
```

For beta>=11/10 the quartic discriminant is negative, so U_beta has off-imaginary zeros and fails the standard Laguerre-Polya universal-multiplier test. In particular the tail-matched beta is not certified by that theorem.

## Live Handoff

Any surviving Gasper comparison must now use at least one signed block or a coupled square identity; decomposing Phi into a positive 9/5 Gasper block plus a separately Laguerre-Polya positive residual is impossible. The live target is a sign-aware multi-block identity whose mixed terms are controlled together after modular cancellation.

The comparison programme has therefore moved from positive block
domination to a genuinely signed, coupled-square problem.

References: https://arxiv.org/abs/0801.2996 and
https://arxiv.org/abs/1502.06844.
