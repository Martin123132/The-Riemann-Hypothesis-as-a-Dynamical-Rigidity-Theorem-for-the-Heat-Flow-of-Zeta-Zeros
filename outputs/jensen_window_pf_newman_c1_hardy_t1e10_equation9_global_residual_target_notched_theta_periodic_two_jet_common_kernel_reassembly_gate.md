# Two-jet reassembly of the common projector kernel

Date: 2026-08-23

Status: exact finite-regulator target-current cancellation and absolutely
convergent symmetric common-kernel representation certified; bound open

Let `hat Phi_x(m)` be the Fourier coefficient of the periodic `C^1` defect
from the two-jet extraction gate.  For every nonzero integer `m`, integration
by parts gives

```text
I_m=2*pi*i*m*hat Phi_x(m)+i*K_x/(2*pi*m),
I_0=E_x.                                                (TR1)
```

The second term in (TR1) is the extracted first-derivative endpoint jet.  In
the target-notched source integral it contributes

```text
-i*K_x*H_T(epsilon)/(2*pi).                            (TR2)
```

For every positive target mode, however,

```text
A_m+B_m=I_m-P_m,                                      (TR3)
```

where `P_m` is the exact full-line Gamma bulk.  Substitution of (TR1) into
(TR3) contributes the opposite current

```text
+i*K_x*H_T(epsilon)/(2*pi).                            (TR4)
```

Thus (TR2) and (TR4) cancel coefficientwise before any norm.  The target
Fourier coefficients of `Phi_x` are restored at the same time.  With
`w_m=exp(-pi*epsilon*m^2)` and every common cutoff `M>=B=5122421`, the complete
Gamma-normalized projector kernel is exactly

```text
Delta_(M,epsilon)(x)
 =H_x+E_x
  +2*pi*i*sum_(m=1)^M m*w_m
     [hat Phi_x(m)-hat Phi_x(-m)]
  -sum_(m=622)^39894 w_m*P_m(x).                      (TR5)
```

The apparent target harmonic current is therefore not an additional error
term and must not be counted after (TR5).  Since `Phi_x` is periodic through
its first derivative,

```text
|hat Phi_x(m)|<=K_Phi(x)/(2*pi*|m|)^3,                (TR6)
```

so the symmetric series in (TR5) is absolutely convergent.  The finite source
roster makes `K_Phi(x)` continuous on `0<=x<=1`; hence its maximum is finite
and the convergence is uniform on the Kummer interval.  This licenses the
zero-regulator representation structurally, but the raw derivative maximum is
not asserted to be quantitatively useful.

The post-A object remains

```text
mathfrak R_A=Delta-chi_W*Btr-O-mathcal A,              (TR7)
```

with exactly the ownership and physical limit order already certified.  No A
or B subtraction has been discarded or counted twice.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate.py
```

Pi provenance: every `pi` in (TR1)--(TR7) is inherited from the Kummer
quadratic phase, integer Fourier character, and Gaussian Abel regulator.  No
fitted or geometric occurrence is introduced.

Proof boundary: exact finite-regulator coefficient reassembly, cancellation
of the target harmonic endpoint jet, and absolute/uniform convergence of the
periodic two-jet pair series only.  No useful numerical tail constant,
physical Kummer quadrature, A/B-subtracted `R_after_A` or `R_Dir` bound,
complete `Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
