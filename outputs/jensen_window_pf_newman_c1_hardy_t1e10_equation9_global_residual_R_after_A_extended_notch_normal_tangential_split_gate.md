# Post-A extended-notch normal/tangential split

Date: 2026-08-23

Status: exact edge-translation identity and mixed endpoint/Morse geometry
certified; quantitative signed assembly remains open

Let `T={622,...,39894}`, `U={622,...,39936}`, and retain the common
Gaussian regulator.  Since `U minus T={39895,...,39936}`, exactly

```text
C_U,epsilon(s)-C_T,epsilon(s)
 =-sum_(m=39895)^39936
    exp(-pi*epsilon*m^2)exp(-2*pi*i*m*s).             (NT1)
```

On the continuous interpolants this is the 42-mode edge translation

```text
h_U,epsilon(y)-h_T,epsilon(y)
 =-q_epsilon(y)1_(39894.5<y<39936.5),                (NT2)
```

up to measure-zero endpoint conventions.  Poisson summation therefore gives

```text
sum_(k in Z)[hat h_U,epsilon(k+s)-hat h_T,epsilon(k+s)]
 =-sum_(m=39895)^39936 w_m exp(-2*pi*i*m*s).          (NT3)
```

At zero regulator the right side is the finite polynomial

```text
-exp(-pi*i*79831*s)sin(42*pi*s)/sin(pi*s),            (NT4)
```

with continuous value `-42` at integer `s`.  Thus the apparent analytical
gain from moving the upper edge is exactly balanced by the explicit 42-mode
outside-target half of the A window.  The edge and that block must stay in
one signed assembly.

For the exact joint phase

```text
Phi(alpha,x;m)
 =pi*alpha^2*x/4-pi*m*alpha+(t/2)log((1-x)/x),

partial_alpha Phi(A,x;m)=pi(A*x/2-m),
partial_x Phi(A,x;m)=pi*A^2/4-t/[2*x(1-x)].           (NT5)
```

On `0<=x<=1/2`, the old and new upper edges obey

```text
|partial_alpha Phi(A,x;39894.5)| >= pi/4,
|partial_alpha Phi(A,x;39936.5)| >= 169*pi/4.         (NT6)
```

The edge shift therefore improves the uniform A-face normal denominator by
the exact factor `169`.  At the tangential stationary point

```text
x_*=[1-sqrt(1-8t/(pi*A^2))]/2
   =[0.49947538036472974127653905997513611415517135950293767675876951565355152221205792 +/- 2.86e-81],                        (NT7)
```

the mode-unit normal gaps are

```text
old: [42.108613768761038156863213173852155730109982299857177932918500276604369983216422 +/- 4.85e-79]
new: [84.108613768761038156863213173852155730109982299857177932918500276604369983216422 +/- 4.85e-79].       (NT8)
```

But `partial_x Phi(A,x_*;m)=0` for every `m`, and the positive tangential
curvature is `[83939326.461260164873455052806884843660556268898598797452281582630954798389243676 +/- 2.27e-73]`.  The new edge is
therefore uniformly nonstationary in the normal direction, not in the
tangential direction.  The licensed order is a signed normal endpoint
reduction followed by a tangential Morse/Fresnel estimate.  Integrating by
parts through `x_*` is invalid.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate.py
```

Pi provenance: every `pi` in (NT1)--(NT8) is inherited from the Gaussian
Abel/Fourier kernel and the exact Kummer phase.  No fitted or geometric
occurrence is introduced.

Proof boundary: exact fixed-regulator edge translation, finite-block
ownership, uniform normal gap, and tangential Morse classification only.
No endpoint-expansion remainder, signed 42-mode cancellation bound,
quantitative `R_after_A`, `R_Dir`, or `Q_K-T` estimate, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
