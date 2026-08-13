# Continuous-height prototype full-line join

Date: 2026-08-12

Status: one prototype crossing mode is certified throughout its full
`pi/8` height cell; this is not a proof of 399-event propagation or RH

Let `tau=pi*m*(C-2m)`, `t*=pi*C^2/8`, `beta=(t*)^(1/3)`, and
`h=t-tau`.  Since `eta=t*/t`, the two phases in the common full-line
chart simplify exactly to

```text
Theta_exact(t)=t*z/beta-t* tanh(z/beta)+normal terms,
Theta_Airy(t)=z^3/3+(t-t*)z/beta+normal terms.       (CH1)
```

Every remaining amplitude and normal Fresnel factor is independent of `t`.
Consequently their signed difference satisfies the exact transport identity

```text
D_t(z,y)=exp(i*h*z/beta) D_tau(z,y).                 (CH2)
```

Write `h=(pi/16) theta`, `|theta|<=1`, and
`kappa=pi/(16 beta)`.  On the finite deformed contour we certified the
moments `M_n=int z^n D_tau(z) dz` through degree
12 and used

```text
Delta I(theta)=sum_(n=0)^N (i*kappa*theta)^n M_n/n!+R_N(theta).
```

The exponential Taylor remainder is bounded against a direct Arb enclosure
of the finite-contour absolute mass.  The already certified tails on
`Im z=1` acquire at most the exact factor `exp(kappa)`.  This gives uniformly
for every real `t` in `[tau-pi/16,tau+pi/16]`

```text
|Delta I(t)| <= [1.724523524699450268019997313235305614036585717485405793944120092782962311981291897163e-5 +/- 5.03e-90],
(sqrt(2)/pi)|Delta I(t)| <= [7.763083334418385926749928496183770915805861924707203271091478460611901357420196129833e-6 +/- 9.61e-91].
```

Both are below `0.000019` and `0.0000086`, respectively.  The interval
reconstruction also overlaps the separately certified lower face, event,
and upper face balls.

## Boundary

This closes continuous height only for the present one-mode prototype and
its present event parameters.  It does not yet make the constants uniform
in the event index, treat all 399 event cells, control the remaining modes,
establish complete `T_upper`, prove `Lambda<=0`, prove RH, or establish a
prize-level conclusion.
