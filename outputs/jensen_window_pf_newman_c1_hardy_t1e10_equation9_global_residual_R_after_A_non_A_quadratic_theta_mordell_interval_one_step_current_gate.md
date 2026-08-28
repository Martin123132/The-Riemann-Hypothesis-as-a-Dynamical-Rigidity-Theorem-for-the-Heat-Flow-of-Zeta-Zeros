# Certified Mordell interval one-step current

Date: 2026-08-24

Status: uniform central-strip tail bounds and two full-roster pointwise
Mordell current enclosures certified; parameter-box recursion open

## Exact tail bounds

For `theta=exp(i*pi/4)`, Kuznetsov's rotated representation is

```text
h(z,tau)=2 theta integral_0^infinity exp(-pi tau y^2)
          cosh(2pi z theta y)/cosh(pi theta y) dy.     (MI1)
```

For real `|z|<=1/2`, set `u=pi*y/sqrt(2)`.  The elementary identities

```text
|cosh(a+ib)|^2=sinh(a)^2+cos(b)^2,
|sinh(a+ib)|^2=sinh(a)^2+sin(b)^2                   (MI2)
```

give, for `y>=Y>0`,

```text
|cosh(2pi z theta y)/cosh(pi theta y)|<=coth(u),
|sinh(2pi z theta y)/cosh(pi theta y)|<=coth(u).     (MI3)
```

Since `coth(u)` decreases, truncating (MI1) and its `z` derivative at `Y`
has the uniform modulus bounds

```text
T_h(Y,tau)
 <=coth(pi Y/sqrt(2))*erfc(sqrt(pi tau)Y)/sqrt(tau),

T_hz(Y,tau)
 <=2*coth(pi Y/sqrt(2))*exp(-pi tau Y^2)/tau.        (MI4)
```

No asymptotic constant is hidden in (MI4).  The finite segments are enclosed
by Arb/ACB complex ball integration.  Each complex tail disk is added as a
containing real-imaginary rectangle.  Arguments outside the central strip
are transported by the exact Mordell unit recurrence and negative `tau` by
conjugation.

## Full one-step certificates

Production uses 320 bits,
`Y=10`, and absolute integration
tolerance `1e-55`.  At
`L=2481422`, both `x=1/2,s=0` and `x=2/5,s=1/3` complete-current balls overlap
their independently enclosed exact rational-period source currents, and the
ball for their difference contains zero.

```text
largest complete radius / source mass:
  2.77496430578298726880530328172e-60,

largest h tail bound:
  2.1087487593897386755726628813e-56,

largest h_z tail bound:
  1.33019665275897183720721926882e-54,

largest main-plus-endpoint condition ratio:
  3734933.28525266563519835472107.        (MI5)
```

The transformed periods are `[2, 6]` and the
source periods are `[4, 15]`.  Thus neither million-term
sum is trusted to floating recurrence in this certificate.

## Decision

Equation (MI4) closes the analytic truncation tail uniformly on the central
strip, and the two physical test rows now have end-to-end rigorous one-step
enclosures.  This removes the nonrigorous quadrature qualification from
those two rows.  It does not yet give a parameter-box evaluator over all
`(x,s)`, a recursive interval algorithm, or the small-`tau`
Euler--Maclaurin branch.  Those remain necessary before physical quadrature
or a non-A bound.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_interval_one_step_current_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_interval_one_step_current_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_interval_one_step_current_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_interval_one_step_current_gate.py
```

Primary-source boundary: (MI1) and the exact unit recurrence are taken from
Alexey Kuznetsov, *Computing the truncated theta function via Mordell
integral*, arXiv:1306.4081v2.  The tail estimates (MI2)--(MI4) are derived
here.  The paper's practical Gauss--Laguerre rule is not used.

Pi provenance: every `pi` comes from the inherited quadratic Fourier phase
or the stated rotated Mordell representation.  No fitted or geometric
constant is inserted.

This gate proves (MI2)--(MI4) and certifies two full-roster pointwise interval
rows.  It does not prove a uniform parameter-box Mordell evaluator, the
small-`tau` branch, recursive interval control, physical quadrature, the
non-A bound, joined `R_after_A`, `R_Dir`, `Q_K-T`, an all-height theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
