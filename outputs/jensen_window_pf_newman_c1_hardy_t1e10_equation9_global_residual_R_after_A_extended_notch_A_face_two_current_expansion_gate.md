# Translated A-face two-current expansion

Date: 2026-08-23

Status: exact certificate for the first two normal currents and uniform
Fresnel-tail remainder; signed edge/block Morse assembly remains open

Put

```text
d_A=39936.5,
delta_A(x)=d_A-A*x/2,
q_A(x)=-delta_A(x)sqrt(2/x),
a_A(x)=d_A sqrt(2/x),              0<x<=1/2.          (TC1)
```

Since `d_A-A/4=169/4`, `|q_A|` decreases to its endpoint minimum

```text
q_A(x)<=-169/2=-84.5.                                (TC2)
```

For `q<0`, retain the exact lower Fresnel tail

```text
J_-(q)=F(q)+(1+i)/2
      =integral_(-infinity)^q exp(i*pi*v^2/2)dv.      (TC3)
```

Repeated integration by parts gives

```text
J_-(q)
 =exp(i*pi*q^2/2)
   [1/(i*pi*q)+1/((i*pi)^2*q^3)]+r_2(q),

|r_2(q)|<=6/(pi^3|q|^5).                             (TC4)
```

The bound follows by exposing the next `3/((i*pi)^3q^5)` boundary term and
bounding the remaining `v^-6` integral absolutely.  It is not a formal
asymptotic equality.

Use the continuous lower-endpoint current inherited from the exact finite
endpoint decomposition,

```text
L_A(d,x)=-exp(i*pi*q_A^2/2)/(i*pi)-a_A J_-(q_A).      (TC5)
```

The bare endpoint exponential and the first Fresnel-tail term cancel
algebraically.  With `E_A=exp(i*pi*q_A^2/2)`, the exact two-current form is

```text
L_A(d_A,x)=E_A[c_0(x)+c_1(x)]+R_2(x),

c_0(x)= A*x/[2*i*pi*delta_A(x)],
c_1(x)=-d_A*x/[2*pi^2*delta_A(x)^3],

|R_2(x)|
 <=3*d_A*x^2/[2*pi^3*delta_A(x)^5].                  (TC6)
```

All three moduli in (TC6) increase on the half-domain.  Their certified
endpoint maxima are

```text
|c_0| <= [300.56175566598085202356407464728667870263866140029371268410073092173460681747891 +/- 1.25e-78]
|c_1| <= [0.013413129711267953330719296106310526399048064109263654055414797153951257968760117 +/- 8.76e-83]
|R_2| <= [3.5877162215010314125233977262322975103185990957651822905521245538522008135148908e-6 +/- 8.00e-87]. (TC7)
```

Quadratic completion retains the exact A-face phase:

```text
exp(-i*pi*d_A^2/x)E_A
 =exp(i*pi*A^2*x/4-i*pi*d_A*A).
```

Here `A` and the numerator `79873` are both `1 mod 4`, so

```text
exp(-i*pi*d_A*A)=-i.                                  (TC8)
```

Thus both displayed currents and the remainder enter the same tangential
phase `-i exp(i Phi_A(x))`; no fitted phase alignment is being used.

The finite-block boundary is exact:

```text
q_A(39936,1/2)=-83.5,
q_A(d_A,1/2)=-84.5.                                  (TC9)
```

The translated tail starts one normal-coordinate unit beyond the final
exact mode `39936`.  Equations (TC6)--(TC9) must therefore be joined to the
already-certified modes `39895..39936` before norms.  The small remainder in
(TC7) does not make the leading current small; its endpoint modulus is about
`[300.56176 +/- 4.34e-6]`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_two_current_expansion_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_two_current_expansion_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_two_current_expansion_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_two_current_expansion_gate.py
```

Pi provenance: every `pi` in (TC1)--(TC9) is inherited from the exact
Kummer/Fourier phase and canonical Fresnel primitive.  No fitted or
geometric occurrence is introduced.

Proof boundary: exact continuous-edge Fresnel expansion, first two normal
currents, uniform pointwise remainder, common A-face phase, and roster
adjacency only.  No signed edge/42-mode cancellation bound, tangential Morse
integral, quantitative `R_after_A`, `R_Dir`, or `Q_K-T` estimate, all-height
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
