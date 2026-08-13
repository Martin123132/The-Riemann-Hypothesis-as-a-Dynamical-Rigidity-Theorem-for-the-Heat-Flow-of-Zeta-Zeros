# Finite-Poisson portcullis saddle reduction

Date: 2026-08-10

Status: exact finite-Poisson and joint-saddle reduction validated; not a proof of the uniform remainder

Insert the exact incomplete-theta derivative from the equation-(10) current
gate into the equation-(9) Kummer integral.  Finite Poisson summation on the
contiguous odd-alpha roster gives dual mode `m` with joint phase

```text
Phi(alpha,x;m,t)
 =pi*alpha^2*x/4-pi*m*alpha+(t/2)log((1-x)/x).
```

The two stationary equations have the unique positive solution

```text
x_m     =2*pi*m^2/(t+2*pi*m^2),
alpha_m =2*m+t/(pi*m).                                  (P1)
```

Equation (P1) is exactly paper equation (124) with `m=N`.  This derives the
appearance of `pi`: it comes from the quadratic Kummer phase and the Fourier
dual frequency, not from a fitted circle or polygon.

The full Hessian is nondegenerate:

```text
det Hess(Phi)= -(2*pi*m**2 + t)**3/(8*m**2*t) < 0.           (P2)
```

Thus the turning point of the one-dimensional map `alpha(N)` is a projection
effect.  The global two-variable saddle itself does not degenerate.

At `t=10^10`, interval arithmetic certifies

```text
N_-(159577) in (39852,39853),
N_+(159577) in (39936,39937),
N_t            in (39894,39895),
N_-(5122421) in (621,622),
N_+(5122421) in (2560588,2560589).
```

The interior lower branch therefore contains modes
`622..39852`.  The lower-alpha boundary gap has 84
modes: `39853..39894` below `N_t` and
`39895..39936` above it, exactly 42 on each side.  The
classical source upper range `622..39894` consists of 39231
interior modes plus those 42 lower endpoint-Fresnel modes.

For each Poisson mode, the alpha integral has the exact grouped form

```text
J_m=(-1)^m/x*{[E_m(B)-E_m(A)]/(i*pi)+m*F_m}.
```

The boundary current and Fresnel integral must remain paired for each `m`.
Their two series need not converge separately.  Likewise, the continuous
reflection `N*=t/(2*pi*N)` satisfies `x_(N*)=1-x_N`, but it does not preserve
the integer mode lattice in general.  The reflected `x` branches must be
recombined by the exact finite-Poisson/theta formula, not by rounding a
reciprocal partner.

The next analytic target is now precise: prove an endpoint-uniform finite
Poisson/stationary-phase theorem for the weighted Kummer integral, retain the
two reflected branches and all 84 turning modes, and bound symmetric
nonstationary mode tails.  Then compare the resulting exact dual expression
with `T_upper` and only afterward insert the source cubic-saddle truncation and
Gaussian evaluator defects.

Proof boundary: exact algebra, symmetric finite-Poisson identity under its
stated smoothness assumptions, and finite endpoint classification only.  No
interchange theorem for the Kummer endpoint singularities, explicit stationary
remainder, source-error bound, `Lambda<=0`, RH, or prize-level conclusion is
proved.
