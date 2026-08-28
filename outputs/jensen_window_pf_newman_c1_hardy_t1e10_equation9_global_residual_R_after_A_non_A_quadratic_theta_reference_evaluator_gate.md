# Non-A quadratic-theta reference evaluator

Date: 2026-08-24

Status: one-worker floating reference evaluator validated against independent
high-precision and exact-periodic oracles; interval error theorem and physical
quadrature open

The exact non-A current from Section 11.460 is

```text
S_0(x,s)=sum_(n=0)^(L-1) z_n,
S_1(x,s)=sum_(n=0)^(L-1) n z_n,
z_n=exp(i*pi*x[n^2+n(A+2s)]),

F_x(s)=exp(i*pi*x(A+2s)^2/4)[(A+2s)S_0+2S_1],
L=2481422.                                               (RE1)
```

The reference evaluator uses

```text
z_(n+1)=z_n r_n,
r_(n+1)=r_n exp(2*pi*i*x).                          (RE2)
```

It never sends a large unreduced phase to the platform trigonometric
functions.  At the start of every 256-term block, `z_n` and `r_n`
are reseeded from the exact rational input phases reduced modulo two.  Real
and imaginary parts of `S_0` and `S_1` are accumulated by blockwise
`math.fsum`.  The implementation is O(L), uses one below-normal worker, and
is the audit oracle for a later fast incomplete-Gauss evaluator.

The short-roster gate compares four nonperiodic decimal/rational phase rows,
through length 4096, against direct 90-digit term evaluation.  The
full-roster gate uses an independent exact compression.  For rational `x,s`,

```text
z_(n+P)/z_n
 =exp(i*pi[2xPn+xP(P+A+2s)]).                       (RE3)
```

Consequently `P` is a period whenever

```text
xP is an integer,
xP(P+A+2s) is an even integer.                      (RE4)
```

Writing `L=kP+R`, one cycle gives exactly

```text
S_0=k Z_0+Z_(0,R),
S_1=P k(k-1)Z_0/2+k Z_1+kP Z_(0,R)+Z_(1,R),        (RE5)
```

with the analogous affine-cycle formula applied directly to the original
weighted current `F_x(s)`.  Thus the seven full-roster rows are checked
without replaying 2,481,422 terms in the independent oracle.  Their periods
range from 1 to 4096.

All 11 rows pass the declared mass-normalized floating
allowance `2.0e-10`.  The largest observed normalized
errors are

```text
S_0: 2.575944e-13
S_1: 3.896109e-13
F:   2.577557e-13.       (RE6)
```

These denominators are the corresponding absolute term masses, not a claim
that cancellation is absent.  Agreement at rational test phases validates
the implementation and catches long-run phase drift; it does not bound the
roundoff error uniformly for every physical `x,s`.

The next stage is to build a genuinely sublinear incomplete quadratic-Gauss
evaluator, cross-check it against this O(L) oracle on deterministic points,
and only then compose it with the pole-free `C_U,epsilon` cell kernel.  A
rigorous physical computation additionally needs interval phase reduction,
summation-error enclosures, `x,s` quadrature remainders, and the ordered
`M`-then-`epsilon` limit.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.py
```

Pi provenance: `pi` in (RE1)--(RE3) is inherited from the original Kummer
quadratic phase.  Rational periodicity only reduces that existing phase
modulo its ordinary `2*pi` period; it introduces no fitted or geometric
constant.

Proof boundary: floating implementation validation and exact rational-phase
periodicity/cycle identities only.  No uniform floating-error theorem,
interval-certified theta evaluator, fast incomplete-Gauss algorithm,
physical `x,s` quadrature, non-A bound, joined `R_after_A`, `R_Dir`, `Q_K-T`,
all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
