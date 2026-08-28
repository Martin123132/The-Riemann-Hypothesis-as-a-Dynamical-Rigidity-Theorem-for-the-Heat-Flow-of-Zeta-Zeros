# Periodic two-jet extraction of the target-notched theta current

Date: 2026-08-23

Status: exact finite-regulator endpoint-current extraction and absolutely
convergent differentiated leakage certified; quantitative `R_after_A` bound open

Let `e_D(x)=exp(i*pi*x*D^2/4)` for `D` equal to `A=159577` or `B=5122421`, and
retain the one-cell objects from Section 11.445.  Put

```text
D_x=e_B(x)-e_A(x),
E_x=D_x/(i*pi*x),              E_0=(B^2-A^2)/4,
K_x=B*e_B(x)-A*e_A(x),        K_0=B-A=2L.             (TJ1)
```

The endpoint divided difference has the stable entire form

```text
E_x=((B^2-A^2)/4)e^[i*pi*x*(A^2+B^2)/8]
    sinc[pi*x*(B^2-A^2)/8],                           (TJ2)
```

where `sinc(z)=sin(z)/z` and `sinc(0)=1`.  For `x>0`, define

```text
Psi_x(s)
 =[Theta_L(x,s)-Theta_L(x,0)-s*D_x]/(i*pi*x),

Phi_x(s)=Psi_x(s)-(K_x/2)(s^2-s).                    (TJ3)
```

Both quotients extend continuously to `x=0`, and the roster sums collapse to

```text
Psi_0(s)=L(s^2-s),             Phi_0(s)=0.            (TJ4)
```

The value and first-derivative jumps are now removed exactly:

```text
Phi_x(0)=Phi_x(1)=0,
partial_s Phi_x(0)=partial_s Phi_x(1).                (TJ5)
```

For `T={622,...,39894}`, set

```text
H_T(epsilon)=sum_(m in T)e^(-pi*epsilon*m^2)/m.
```

At every common finite cutoff `M>=B`, Fourier-pair cancellation gives

```text
integral_0^1 (s^2-s)C_(M,epsilon)'(s)ds
 =i*H_T(epsilon)/pi.                                  (TJ6)
```

Since `integral_0^1 C_(M,epsilon)=1`, integration by parts with (TJ5) gives
the exact joined source identity

```text
integral_0^1 F_x(s)C_(M,epsilon)(s)ds
 =E_x-i*K_x*H_T(epsilon)/(2*pi)
  -integral_0^1 Phi_x(s)C_(M,epsilon)'(s)ds.          (TJ7)
```

Thus the coherent endpoint value and derivative currents are extracted
before any norm.  At zero regulator,

```text
H_T(0)=[4.16185772996900304025207600348705642984596876815154322119292404133349293634 +/- 7.20e-76],
H_T(0)/pi=[1.32476046033956276831194354001237996844173426659028391392901539427103246713 +/- 4.26e-75].              (TJ8)
```

If `hat Phi_x(m)=integral_0^1 Phi_x(s)e^(-2*pi*i*m*s)ds`, three integrations
by parts use (TJ5) to give

```text
|hat Phi_x(m)|
 <=K_Phi(x)/(2*pi*|m|)^3,
K_Phi(x)=|Delta Phi_x''|+integral_0^1|Phi_x'''(s)|ds,

sum_(|m|>K)2*pi*|m|*|hat Phi_x(m)|
 <=K_Phi(x)/(2*pi^2*K).                               (TJ9)
```

The differentiated leakage in (TJ7) is therefore absolutely convergent after
the two endpoint jets are removed, including when the Abel regulator tends to
zero.  The raw derivative majorant in (TJ9) is only a convergence certificate;
it is not claimed to be quantitatively useful.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_extraction_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_extraction_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_extraction_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_extraction_gate.py
```

Pi provenance: every `pi` in (TJ1)--(TJ9) is inherited from the Kummer
quadratic phase, integer Fourier character, and Gaussian Abel regulator.  No
fitted or geometric occurrence is introduced.

Proof boundary: exact finite-regulator two-jet periodization, closed target
harmonic current, and absolute convergence of the differentiated leakage only.
No quantitative finite-current evaluator, cancellation with the separately
owned A/B extraction terms, Kummer `x` bound, `R_after_A` or `R_Dir` bound,
complete `Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
