# Pole-free target-notch cell and one-cell source fold

Date: 2026-08-23

Status: exact punctured-cell recombination, uniform cubic tail, and one-cell
fold certified; Kummer-weighted quadrature and `R_after_A` bound open

Retain the two-edge objects `hat h`, `b_0`, `b_1`, and `rho_2` from Section
11.444.  On the centred cell `|s|<=1/2`, put

```text
g_0(s)=[csc(pi*s)-1/(pi*s)]/(2i),
g_1(s)={1/(pi*s)^2-cos(pi*s)/sin(pi*s)^2}/4.         (PC1)
```

Their apparent singularities are removable:

```text
g_0(s)=(1/(2i))[pi*s/6+7(pi*s)^3/360+...],
g_1(s)=1/24+7(pi*s)^2/480+...,
g_0(0)=0,  g_1(0)=1/24.                              (PC2)
```

The `k=0` residual identity is

```text
rho_2(s)=hat h_epsilon(s)-b_0(s)-b_1(s).
```

Combining it with the two closed boundary sums before evaluation gives the
single pole-free formula

```text
C_epsilon(s)
 =hat h_epsilon(s)+E_0(s)g_0(s)+E_1(s)g_1(s)
  +sum_(k in Z, k!=0)rho_2(k+s),       |s|<=1/2,     (PC3)

C_epsilon(0)
 =hat h_epsilon(0)+E_1(0)/24+sum_(k!=0)rho_2(k).     (PC4)
```

No direct high-cutoff switch is needed at the integer lattice.  Uniformly on
the whole centred cell,

```text
sum_(|k|>K)|rho_2(k+s)|
 <=K3/[(2pi)^3(K-1/2)^2],
K3=|q_epsilon''(c)|+|q_epsilon''(d)|+20pi*epsilon.   (PC5)
```

For `K=64`, the certified tails are

```text
epsilon=1e-6: [6.548302940258891799542952621948798162414800652753522509695882672469211e-11 +/- 4.54e-81]
epsilon=1e-8: [6.887479549317750506554295966713475103808344402866605972845544097886508e-13 +/- 4.19e-83]
epsilon=1e-10: [6.909895073078900089481325342875753204122569115857216372961103461136594e-15 +/- 3.74e-85]. (PC6)
```

The common complement is one-periodic.  Since `L=2481422` is an integer, define

```text
Theta_L(x,s)
 =sum_(n=0)^(L-1)e^[i*pi*x*(A+2n+2s)^2/4],

F_x(s)
 =sum_(n=0)^(L-1)(A+2n+2s)e^[i*pi*x*(A+2n+2s)^2/4]
 =(i*pi*x)^(-1)partial_s Theta_L(x,s),       x>0.     (PC7)
```

Cellwise change of variables gives the exact fold

```text
integral_0^L f_x(u)C_epsilon(u)du
 =integral_0^1 F_x(s)C_epsilon(s)ds.                 (PC8)
```

The finite current shifts from labels `A,...,B-2` at `s=0` to
`A+2,...,B` at `s=1`, so

```text
Theta_L(x,1)-Theta_L(x,0)
 =e^(i*pi*x*B^2/4)-e^(i*pi*x*A^2/4).                (PC9)
```

This is the same endpoint difference that drives the certified Poisson and
double-Weber cancellations.  Equation (PC8) is now an executable one-cell
target for a finite quadratic-Gauss current evaluator; it is not yet a
quadrature theorem.

The separate floating pilot passes 21 direct/dual rows, including the exact
integer and points within `10^-5`, down to `epsilon=10^-10`.  Its worst
relative discrepancy is below `3.91e-11`.  Those samples are not used as
proof of (PC1)--(PC9).

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.py
```

Pi provenance: every `pi` in (PC1)--(PC9) is inherited from the Gaussian
Abel/Fourier kernel and the original Kummer quadratic phase.  No fitted or
geometric occurrence is introduced.

Proof boundary: exact fixed-regulator pole cancellation, one-cell folding,
finite theta-current identity, and cubic dual-tail enclosures only.
No finite-Gauss-current error theorem, Kummer `x` integration, A/B extraction
evaluation, quantitative `R_after_A` or `R_Dir` bound, complete `Q_K-T`,
all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
