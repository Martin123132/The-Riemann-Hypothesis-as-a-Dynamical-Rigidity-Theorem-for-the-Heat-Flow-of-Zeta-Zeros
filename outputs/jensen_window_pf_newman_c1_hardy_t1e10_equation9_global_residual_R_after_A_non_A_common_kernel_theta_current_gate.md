# Non-A common kernel and quadratic-theta current

Date: 2026-08-24

Status: exact executable saved-height reduction; quantitative non-A bound open

At one common cutoff `M>=B=5122421` and regulator `epsilon>0`, Section
11.459 defines

```text
mathfrak K_nonA,(M,epsilon)
 =mathfrak R_A,(M,epsilon)-mathfrak E_(A,42),epsilon. (NK1)
```

The coefficient-mask evaluator is

```text
H_x+I_0
 +sum_(m=1)^621       w_m(P_m+P_-m+A_m+A_-m+B_m+B_-m)
 +sum_(m=622)^39852   w_m(P_-m+A_m+A_-m+B_m+B_-m)
 +sum_(m=39853)^39936 w_m(P_-m+B_m+B_-m)
 +sum_(m=39937)^M     w_m(P_m+P_-m+A_m+A_-m+B_m+B_-m)
 -chi_W Btr_(M,epsilon)-O_(M,epsilon)
 -mathfrak E_(A,42),epsilon.                         (NK2)
```

The two adjacent 42-mode rows have been combined in (NK2) because their
surviving coefficient vectors are identical.  No coefficient changed when
the analytic A-face channel was extracted.

Let `U={622,...,39936}`.  The exactly equivalent extended-notch form is

```text
H_x+integral_0^L f_x(u)[D_(M,epsilon)-G_(U,epsilon)]du
 +sum_(m=622)^39852 w_m(A_m+B_m)
 +sum_(m=39853)^39936 w_m(B_m-A_-m)
 -chi_W Btr_(M,epsilon)-O_(M,epsilon)
 -mathfrak E_(A,42),epsilon.                         (NK3)
```

Equivalently, after the exact two-jet reassembly,

```text
H_x+E_x+2*pi*i sum_(m=1)^M m w_m
       [hat Phi_x(m)-hat Phi_x(-m)]
 -sum_(m=622)^39936 w_m P_m
 -sum_(m=39853)^39936 w_m(A_m+A_-m)
 -chi_W Btr_(M,epsilon)-O_(M,epsilon)
 -mathfrak E_(A,42),epsilon.                         (NK4)
```

The pair series in (NK4) is absolutely and uniformly convergent after the
two endpoint jets are removed.  Its inherited raw derivative majorant is
still far too large for the target, so this identity is a convergence-safe
cross-check rather than the selected numerical enclosure.

The source term in (NK3) folds exactly to one cell:

```text
integral_0^L f_x(u)C_(U,M,epsilon)(u)du
 =integral_0^1 F_x(s)C_(U,M,epsilon)(s)ds,           (NK5)

F_x(s)=sum_(n=0)^(L-1)(A+2n+2s)
       exp(i*pi*x(A+2n+2s)^2/4),
L=2481422.                                               (NK6)
```

After `M` tends to infinity at fixed `epsilon`, the certified pole-free
centred-cell formula evaluates `C_U,epsilon` with a cubic dual tail.  The
remaining large object in (NK5) is not arbitrary: put

```text
S_0=sum_(n=0)^(L-1) exp(i*pi*x[n^2+n(A+2s)]),
S_1=sum_(n=0)^(L-1) n exp(i*pi*x[n^2+n(A+2s)]).      (NK7)
```

Exact quadratic completion gives

```text
F_x(s)=exp(i*pi*x(A+2s)^2/4)
       [(A+2s)S_0+2S_1].                            (NK8)
```

For direct reference evaluation, if `z_n` is the summand in `S_0`, then

```text
z_(n+1)/z_n=exp(i*pi*x[A+2n+2s+1]),

(z_(n+2)/z_(n+1))/(z_(n+1)/z_n)=exp(2*pi*i*x).      (NK9)
```

Thus one initial phase and two multiplicative recurrences evaluate both
`S_0` and `S_1` without repeated large-argument trigonometric calls.  At
`x=0`, the removable value is exactly

```text
F_0(s)=L(A+2s+L-1).                                  (NK10)
```

Equations (NK7)--(NK10) define the reference oracle for a fast incomplete
quadratic-Gauss recursion.  The next numerical stage must cross-check that
fast evaluator against this recurrence on short and full rosters, then
compose it with the pole-free `C_U,epsilon` kernel on a bounded lattice.
Neither a finite-Gauss error theorem nor the physical `x,s` quadrature has
yet been proved.

The inherited physical limit is

```text
R_nonA=lim_(epsilon down 0)lim_(M to infinity)
       mathcal P_t[mathfrak K_nonA,(M,epsilon)],      (NK11)
```

and the sufficient target remains

```text
|R_nonA|<0.0331747039947.                            (NK12)
```

No row of (NK2) may be normed independently: its zero, half-current,
negative, B-completion, and remote-positive cancellations are precisely why
(NK3)--(NK5) are retained as common-kernel evaluators.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_common_kernel_theta_current_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_common_kernel_theta_current_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_common_kernel_theta_current_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_common_kernel_theta_current_gate.py
```

Pi provenance: every `pi` in (NK3)--(NK10) is inherited from the original
Kummer quadratic phase, integer Fourier character, Gaussian Abel regulator,
and Poisson normalization.  The algebraic phase factorization introduces no
new geometric or fitted occurrence.

Proof boundary: (NK2)--(NK10) are exact saved-height common-kernel and finite
quadratic-theta reductions.  No fast-evaluator error theorem, physical
quadrature, non-A bound, joined `R_after_A`, `R_Dir`, `Q_K-T`, all-height
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
