# Jensen-Window PF Weighted Fractional Autocorrelation Gram Bridge

Date: 2026-07-23

Status: exact weighted autocorrelation Gram/spectral reduction with one open
signed off-diagonal gate. This is not a proof of RH, PF-infinity, or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.json
python work/rh_compute/scripts/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py
```

## Weighted Autocorrelation

Fix

```text
0 < omega < 1/2,
alpha = 2*omega,
a = (1-alpha)/2 = 1/2-omega,
beta = (1+alpha)/2 = 1/2+omega.
```

For `lambda>0`, define

```text
A_alpha(lambda)
  = integral_(0,infinity) {x}{lambda*x} x^(alpha-2) dx.       (WFAGB.1)
```

The integrand is `O(x^alpha)` at zero and `O(x^(alpha-2))` at
infinity, so (WFAGB.1) converges absolutely. Substitution in
`A_alpha(1/lambda)` gives

```text
A_alpha(lambda)
  = lambda^(1-alpha) A_alpha(1/lambda).                       (WFAGB.2)
```

The logarithmically stationary normalization is

```text
C_alpha(u) = exp(-a*u) A_alpha(exp(u)).                       (WFAGB.3)
```

Equation (WFAGB.2) says exactly that `C_alpha(-u)=C_alpha(u)`.

## Mellin And Fourier Forms

In `-1<Re(s)<0`,

```text
integral_(0,infinity) {x} x^(s-1) dx = zeta(-s)/s.            (WFAGB.4)
```

Using (WFAGB.4) twice and first working in the common absolute
convergence strip `alpha-1<Re(s)<0` gives

```text
integral_(0,infinity) A_alpha(lambda) lambda^(s-1) d(lambda)
  = -zeta(-s) zeta(1-alpha+s) / [s(s+1-alpha)].               (WFAGB.5)
```

On the central line `s=-a+it`, the two zeta factors are conjugates:

```text
Phi_alpha(t) = |zeta(a+it)|^2 / (a^2+t^2),                   (WFAGB.6)

C_alpha(u)
  = (1/(2*pi)) integral_R Phi_alpha(t) exp(-itu) dt.          (WFAGB.7)
```

The standard convexity bound for zeta makes `Phi_alpha` integrable
because `a>0`. In particular, `Phi_alpha>=0`, so `C_alpha` is positive
definite and

```text
|C_alpha(u)| <= C_alpha(0).                                  (WFAGB.8)
```

## Exact Burnol Gram Form

Set

```text
f_(alpha,N)(t)
  = sum_(d<=N) mu(d)d^(-alpha) {1/(d*t)}.
```

The reciprocal substitution `x=1/t` gives

```text
integral_(0,infinity) t^(-alpha)|f_(alpha,N)(t)|^2 dt
  = integral_(0,infinity) x^(alpha-2)
      |sum_(d<=N) mu(d)d^(-alpha){x/d}|^2 dx.                 (WFAGB.9)
```

For the cross-kernel,

```text
I_alpha(d,e)
  = integral_(0,infinity) x^(alpha-2){x/d}{x/e} dx
  = d^(alpha-1) A_alpha(d/e)
  = (de)^((alpha-1)/2) C_alpha(log(d/e)).                    (WFAGB.10)
```

Consequently, the common value in (WFAGB.9) is

```text
Gamma_(alpha,N)
  = sum_(d,e<=N) mu(d)mu(e)(de)^(-beta)
      C_alpha(log(d/e))                                      (WFAGB.11)

  = (1/(2*pi)) integral_R Phi_alpha(t)
      |sum_(d<=N) mu(d)d^(-beta-it)|^2 dt.                   (WFAGB.12)
```

This is exactly the Burnol-Hardy Mellin-Plancherel norm already used in
the fixed-shift bridge: `beta=1/2+omega` and `a=1/2-omega`. Thus the
reciprocal-cell, Mellin, and stationary-Gram formulations have now been
normalization-checked against one another.

## The Real Remaining Gate

Split (WFAGB.11) into

```text
Diag_(alpha,N)
  = C_alpha(0) sum_(d<=N) mu(d)^2 d^(-(1+alpha)),             (WFAGB.13)

Off_(alpha,N)
  = 2 sum_(d<e<=N) mu(d)mu(e)(de)^(-beta)
      C_alpha(log(d/e)).                                     (WFAGB.14)
```

The diagonal is already harmless:

```text
0 <= Diag_(alpha,N) <= C_alpha(0) zeta(1+alpha).              (WFAGB.15)
```

All growth is therefore a signed off-diagonal arithmetic question.
Positive definiteness gives a lower bound for the full Gram form, not the
uniform upper bound that the proof programme needs. Indeed, applying
(WFAGB.8) absolutely gives only

```text
Gamma_(alpha,N)
  <= C_alpha(0)(sum_(d<=N)d^(-beta))^2
  = O(N^(1-alpha)),                                          (WFAGB.16)
```

which diverges. Kernel positivity, pointwise positivity, the bounded
diagonal, and the unweighted autocorrelation theorem therefore cannot be
promoted to the missing estimate.

The exact open target is:

> **Weighted signed off-diagonal gate.** Prove
> `sup_N Gamma_(alpha,N)<infinity` along one explicit cofinal sequence
> `alpha->0+`, using a weighted multiplicative large-sieve, bilinear
> Mobius-cancellation theorem, or an equivalent estimate strong enough to
> control (WFAGB.14).

Equation (WFAGB.16) is only a nonuniform barrier estimate.
The weighted signed off-diagonal gate is open.

## Source Boundary

- [Báez-Duarte et al., *A general strong Nyman-Beurling criterion for the
  Riemann hypothesis*](https://arxiv.org/abs/math/0306251) supplies the
  `alpha=0` fractional-part autocorrelation/Mellin prototype.
- [Burnol, *A note on Nyman's equivalent formulation of the Riemann
  hypothesis*](https://arxiv.org/abs/math/0202166) supplies the weighted
  approximants and Hardy-space setting.
- [Báez-Duarte et al., later autocorrelation analysis](https://arxiv.org/abs/1305.4395)
  studies the unweighted autocorrelation further.

The `0<alpha<1` stationary kernel, its finite weighted Mobius Gram form,
and the diagonal/off-diagonal reduction are derived here. None of the
sources is cited as proving the open weighted uniform bound.
