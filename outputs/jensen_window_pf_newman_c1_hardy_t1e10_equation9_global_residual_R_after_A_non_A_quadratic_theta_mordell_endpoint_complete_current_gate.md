# Endpoint-complete Mordell one-step theta current

Date: 2026-08-24

Status: exact endpoint-complete one-step current identity certified;
numerical Mordell validation diagnostic; recursive interval evaluator open

## Scope

Let

```text
F_n(a,tau)=sum_(k=0)^n exp(2*pi*i*(a*k+tau*k^2)),
a=tau*c,                         tau=x/2>0,
J_n=exp(i*pi*a^2/(2*tau))*[c F_n+(pi*i)^(-1)partial_a F_n].   (ME1)
```

Kuznetsov's exact Theorem 1 in
[Computing the truncated theta function via Mordell integral](https://arxiv.org/abs/1306.4081v2)
writes `F_n` as one shorter theta sum plus two explicit Mordell-integral
endpoint terms.  This gate differentiates that identity in the one joined
direction required by the non-A source.  No unnamed constants or code from
the paper are imported.

## Exact joined transform

Choose the unique integer `r` and reduced phase `z` with

```text
a=z+r,                         -1/2<=z<1/2,
m=floor(2*n*tau),              w=z/(2*tau),
sigma=-1/(4*tau),
C=exp(i*pi/4-i*pi*z^2/(2*tau))/sqrt(2*tau).          (ME2)
```

Differentiating the exact Mordell identity and keeping the physical current
tied gives

```text
[a/tau+(pi*i)^(-1)partial_z] [C F_m(w,sigma)]
 =C/tau [r F_m(w,sigma)+(2*pi*i)^(-1)partial_w F_m(w,sigma)]
 =C/tau sum_(k=0)^m (r+k)exp(2*pi*i*(w*k+sigma*k^2)). (ME3)
```

Thus the apparent zeroth and first transformed moments are one absolute-
frequency first moment.  With `v=r+k`, the phase identity

```text
(a^2-z^2)/(2*tau)+z*k/tau-k^2/(2*tau)
 =(a/tau)*v-v^2/(2*tau)                              (ME4)
```

recovers the pure dual current of Section 11.462 without evaluating a huge
unreduced Mordell argument.

Write

```text
u_-=z-tau+1/2,
u_+=z+(2n+1)tau-m-1/2,
E_-=exp[-pi*i*(z-tau/2)],
E_+=exp[2*pi*i*(n+1/2)*(z+tau*(n+1/2))].             (ME5)
```

The endpoint-complete joined current is exactly

```text
J_endpoint=exp(i*pi*a^2/(2*tau))*(-i/2)*{
 E_-[(a/tau-1)h(u_-,-2tau)+h_z(u_-,-2tau)/(pi*i)]
 +(-1)^m E_+[(a/tau+2n+1)h(u_+,-2tau)
              +h_z(u_+,-2tau)/(pi*i)]}.            (ME6)
```

Equations (ME3) and (ME6) are one exact modular step for the full current;
they do not split `R_0` from `R_1`.

## Validation

Two short rows compare (ME3)+(ME6) with independent 90-digit direct sums.
Their largest absolute discrepancy is
`3.79161294924389047227623641004e-88`.  Two full
`L=2481422` rows compare the complete one-step expression with the validated
O(`L`) source oracle.  The transformed weighted sums are also checked
against independent exact rational-period compression, with periods
2 through
6.  The largest discrepancies relative
to absolute term mass are

```text
transformed weighted sum: 6.121841e-17,
complete source current:  1.025821e-16.  (ME7)
```

The Mordell values and derivatives are evaluated twice, at 60 and 90 digits,
from the rotated defining integral after exact unit-recurrence reduction.
The largest cross-precision changes are

```text
h:   2.8046750840405164019e-61,
h_z: 3.5122251476440650858e-61.                 (ME8)
```

These numerical rows validate the specialization and implementation only.
They are not interval enclosures.  Kuznetsov explicitly distinguishes his
rigorously analyzable algorithm from the faster practical Gauss-Laguerre
version and notes that he did not obtain rigorous error estimates for that
practical quadrature.  This gate therefore does not promote (ME7)--(ME8) to
a theorem.

## Decision

The previously open endpoint ownership problem for one modular step is now
closed at the exact-identity level: the remainder is the joined pair (ME6),
not separately estimated `R_0,R_1`.  The next obligation is a certified
evaluator for `h,h_z` and a small-`tau` Euler--Maclaurin branch, followed by a
recursive interval implementation.  Only then can this become a uniform
fast evaluator for the physical `(x,s)` quadrature.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_endpoint_complete_current_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_endpoint_complete_current_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_endpoint_complete_current_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_endpoint_complete_current_gate.py
```

Pi provenance: every `pi` in (ME1)--(ME8) comes from the inherited quadratic
Fourier phase or Kuznetsov's stated Mordell identity.  No fitted or geometric
constant is inserted.

This gate proves the exact endpoint-complete one-step current identity and
the rational transformed-period formulas only.  The four floating rows are
diagnostics.  It does not prove a certified Mordell evaluator, the small-
`tau` branch, recursive interval control, physical quadrature, the non-A
bound, joined `R_after_A`, `R_Dir`, `Q_K-T`, an all-height theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
