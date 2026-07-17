# Jensen-Window PF Newman Strict-Laguerre Correlation Target

Date: 2026-07-17

Status: exact strict-Laguerre/Wiener equivalence with a generic-kernel
guard. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_correlation_target.json
python work/rh_compute/scripts/jensen_window_pf_newman_strict_laguerre_correlation_target.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_correlation_target.py
```

Current result:

```text
validated Jensen-window PF Newman strict-Laguerre correlation target: 10 rows, 0 issues, 1 strict-Laguerre equivalence, 1 exact correlation identity, 1 Wiener-density equivalence, 1 RH-equivalent density target, 1 exact strict-log-concavity/positive-definiteness countermodel, 2 non-promotion gates, 1 open Xi handoff
```

## Strict Laguerre Target

Set

```text
phi_t(u)=exp(t*u^2)*Phi(u), u in R
H_t(x)=integral_0^infinity phi_t(u)*cos(x*u)du
L_t(x)=H_t'(x)^2-H_t(x)*H_t''(x)
K_(1,t)(v)=integral_R phi_t(s+v)*phi_t(s-v)*s^2 ds
```

The positive-boundary attainment theorem turns the first Laguerre
inequality into an exact endgame:

```text
Lambda<=0 if and only if L_t(x)>0 for every real x and every 0<t<=1/5. The forward implication uses the simple real-zero factorization for t>Lambda; the reverse implication uses the finite multiple zero of H_Lambda forced when Lambda>0 together with the published bound Lambda<=1/5.
```

This is weaker in form than proving every `H_t` is Laguerre-Polya:
only one strict real-axis inequality is requested, but it must hold
for the complete positive-time continuum.

## Correlation Identity

The midpoint/difference change of variables gives exactly

```text
L_t(x)=integral_R K_(1,t)(v)*cos(2*x*v)dv; equivalently, Fourier[K_(1,t)](xi)=L_t(xi/2).
```

There is no missing normalization factor. Wiener's theorem now yields

```text
Translations of K_(1,t) are dense iff L_t has no real zero. Since L_t(0)>0, zero-freeness is equivalent to L_t(x)>0 for all x.
Lambda<=0 if and only if, for every 0<t<=1/5, the translations of K_(1,t) are dense in L1(R).
```

Primary sources: https://arxiv.org/abs/1309.0055 and
https://arxiv.org/abs/1606.05011.

At a hypothetical positive boundary the correlation kernel is already
positive definite, but its transform touches zero at the multiple root.
The missing property is therefore zero-freeness or translate density,
not ordinary nonnegativity.

## Exact Shape Guard

Let

```text
T_a(y)=(a-|y|)_+ = 1_[-a/2,a/2]*1_[-a/2,a/2]
G_sigma(x)=exp(-x^2/(2*sigma^2))
K_(a,sigma)(x)=integral_(-a)^a (a-|y|)*exp(-(x-y)^2/(2*sigma^2))dy
```

For `a=1`, `sigma=2`, conditional differentiation gives

```text
With p_x(y) proportional to T_a(y)G_sigma(x-y), (log K)''=Var_x(Y)/sigma^4-1/sigma^2
Since Y lies in [-a,a], Var_x(Y)<=a^2; for sigma>a, (log K)''<=(a^2-sigma^2)/sigma^4<0.
(log K)''<=-3/16
```

so `K` is smooth, positive, even, and strictly log-concave. Its
Fourier transform is nonnegative, hence `K` is positive definite, but

```text
8*sqrt(2)*sqrt(pi)*exp(-2*xi**2)*sin(xi/2)**2/xi**2
xi=2*pi*n/a for every nonzero integer n
```

The zeros are double, so the translates of `K` are not dense. This
blocks promotion from strict log-concavity plus positive definiteness.
This Gaussian model alone does not reproduce the Xi-specific tail.
The separate theta-tail weighted countermodel closes promotion from
uniform strong and root-variable log-concavity plus theta-type decay;
it still leaves Xi arithmetic or modular correlation structure open.

## Live Handoff

Prove uniformly for every 0<t<=1/5 that Fourier[K_(1,t)] has no real zero, equivalently that translations of K_(1,t) are dense in L1(R). A viable proof must use Xi-specific structure beyond generic strict log-concavity and positive definiteness. The exact higher-shift expansion must be grouped through the theta modular identity before the spectral transform: its ordinary termwise Bessel-transform sum diverges already at zero frequency. Continuum-cell subtraction repairs the endpoint transform and gives a normally convergent coupled matrix at t=0, but each cell block retains a nonzero exp(-5u) tail and cannot be deformed to any t>0. The surviving route needs a t-compatible modular grouping or a direct curvature estimate for the endpoint-subtracted theta primitive. The former sharper sufficient subtarget M_t(x)=-L_t'(x)>0 has now been rejected for Xi: Arb certifies M_0(1401016343/100000)<0, and continuity preserves that sign for sufficiently small positive t. The surviving routes are direct zero-freeness of Fourier[K_(1,t)] or corrected C1 double-zero transversality; do not impose global monotonicity of L_t.

The exact route refinement is recorded in
`outputs/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.md`.
The convergent endpoint matrix and its positive-time obstruction are in
`outputs/jensen_window_pf_newman_theta_cell_renormalization_gate.md`.
The rejected monotonicity subtarget, its misleading finite diagnostics,
and the rigorous Xi counterexample are in
`outputs/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.md`.
