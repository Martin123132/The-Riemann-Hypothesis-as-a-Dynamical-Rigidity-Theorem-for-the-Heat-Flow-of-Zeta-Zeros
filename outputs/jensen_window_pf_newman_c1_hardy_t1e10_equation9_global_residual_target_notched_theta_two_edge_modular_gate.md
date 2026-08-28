# Two-edge modular transform of the target-notched theta kernel

Date: 2026-08-23

Status: exact regulated complement transform and cubic dual-tail certificate;
not an evaluation or bound for `R_after_A`

Put

```text
T={622,...,39894},       c=621.5,       d=39894.5,
q_epsilon(y)=exp(-pi*epsilon*y^2),       epsilon>0,

h_epsilon(y)=q_epsilon(y)[1_(y<c)+1_(y>d)].          (TN1)
```

Because `c` and `d` are half-integers, the integer samples of `h_epsilon`
are exactly the target complement.  After the lawful `M->infinity` passage
at fixed positive `epsilon` in Section 11.443,

```text
C_epsilon(u)
 =sum_(m in Z)[1-chi_T(m)]q_epsilon(m)e^(-2pi*i*m*u)
 =sum_(k in Z)hat h_epsilon(k+u).                    (TN2)
```

The second equality is Poisson summation for the integrable bounded-variation
function (TN1), in the inherited symmetric prescription.  Its exact Fourier
transform is

```text
hat h_epsilon(xi)
 =e^(-pi*xi^2/epsilon)/(2sqrt(epsilon))
  {erfc[-sqrt(pi*epsilon)c-i sqrt(pi/epsilon)xi]
   +erfc[ sqrt(pi*epsilon)d+i sqrt(pi/epsilon)xi]}.  (TN3)
```

For `j=0,1`, define

```text
b_j(xi)
 =[q_epsilon^(j)(d)e^(-2pi*i*xi*d)
   -q_epsilon^(j)(c)e^(-2pi*i*xi*c)]
   /(2pi*i*xi)^(j+1),                                (TN4)

E_j(u)=q_epsilon^(j)(d)e^(-2pi*i*u*d)
       -q_epsilon^(j)(c)e^(-2pi*i*u*c).
```

The half-integer phases give `exp(-2pi*i*k*c)=exp(-2pi*i*k*d)=(-1)^k`.
The two standard Mittag-Leffler sums therefore close the boundary currents:

```text
sum_k b_0(k+u)= E_0(u)/(2i sin(pi*u)),
sum_k b_1(k+u)=-E_1(u)cos(pi*u)/(4sin(pi*u)^2).       (TN5)
```

Three integrations by parts, performed on the two exterior intervals before
an absolute value, give for `rho_2=hat h-b_0-b_1`

```text
|rho_2(xi)|
 <=[|q_epsilon''(c)|+|q_epsilon''(d)|+20pi*epsilon]
    /(2pi|xi|)^3.                                    (TN6)
```

Here the explicit constant follows from

```text
q_epsilon'''=(12pi^2 epsilon^2 y-8pi^3 epsilon^3 y^3)
              q_epsilon,
integral_R |q_epsilon'''|<=20pi*epsilon.             (TN7)
```

Consequently, for noninteger `u`,

```text
C_epsilon(u)
 =E_0(u)/(2i sin(pi*u))
  -E_1(u)cos(pi*u)/(4sin(pi*u)^2)
  +sum_(k in Z)rho_2(k+u),                           (TN8)
```

with continuous extension of the complete right side at integer `u`.  The
closed terms and residual must not be evaluated separately near those
removable points; the stable evaluator uses the direct Gaussian sum there
and (TN8) away from the integer lattice.

For `|k|<=64` and `0.1<=u<=0.9`, the certified dual-tail enclosures are

```text
epsilon=1e-6: [6.528937450666721695060600177892742540446103991715442238236289336715962e-11 +/- 1.75e-81]
epsilon=1e-8: [6.867111002727063288042984668013149188701141530822527977244566166608852e-13 +/- 1.42e-83]
epsilon=1e-10: [6.889460236398084202921179195653848520652434069487249017629400230613491e-15 +/- 1.88e-85].     (TN9)
```

The separate floating pilot checks nine real-edge rows down to
`epsilon=10^-10`: its direct cutoff grows to `328014`, while the dual cutoff
stays at `64`, and the largest discrepancy is below `4.19e-11`.  Those
floating comparisons select the route but are not used as proof of (TN1)--
(TN9) and do not estimate `R_after_A`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.py
```

Pi provenance: every `pi` in (TN1)--(TN9) comes from the inherited Gaussian
Abel weight, integer Fourier character, Fourier transform convention, and
standard Jacobi/Poisson normalization.  No geometric or fitted occurrence is
introduced.

Proof boundary: exact fixed-`epsilon` target-complement transform, explicit
boundary-current closure, and cubic dual-tail bound only.  No integration
against the Kummer source, A/B extraction evaluation, numerical or analytic
bound for `R_after_A`, complete `R_Dir` or `Q_K-T`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
