# Post-A extended-notch pole-free kernel

Date: 2026-08-23

Status: exact instantiated two-edge transform, pole-free centred-cell formula,
and cubic dual-tail certificate; the joined stationary estimate remains open

After the exact A-window subtraction, put

```text
U={622,...,39936},       c=621.5,       d_A=39936.5,
q_epsilon(y)=exp(-pi*epsilon*y^2),
h_A,epsilon(y)=q_epsilon(y)[1_(y<c)+1_(y>d_A)].       (EA1)
```

The edge width is exactly `d_A-c=39315=|U|`.  Both edge numerators are odd,
so for every integer `k`

```text
exp(-2*pi*i*k*c)=exp(-2*pi*i*k*d_A)=(-1)^k.          (EA2)
```

Thus the extended complement has the exact fixed-regulator Poisson form

```text
C_A,epsilon(s)
 =sum_(m in Z)[1-chi_U(m)]q_epsilon(m)e^(-2*pi*i*m*s)
 =sum_(k in Z)hat h_A,epsilon(k+s),                  (EA3)

hat h_A,epsilon(xi)
 =e^(-pi*xi^2/epsilon)/(2sqrt(epsilon))
  {erfc[-sqrt(pi*epsilon)c-i sqrt(pi/epsilon)xi]
   +erfc[ sqrt(pi*epsilon)d_A+i sqrt(pi/epsilon)xi]}. (EA4)
```

Define, for `j=0,1`,

```text
b_A,j(xi)
 =[q_epsilon^(j)(d_A)e^(-2*pi*i*xi*d_A)
   -q_epsilon^(j)(c)e^(-2*pi*i*xi*c)]
   /(2*pi*i*xi)^(j+1),

E_A,j(s)=q_epsilon^(j)(d_A)e^(-2*pi*i*s*d_A)
         -q_epsilon^(j)(c)e^(-2*pi*i*s*c).           (EA5)
```

With `rho_A,2=hat h_A,epsilon-b_A,0-b_A,1`, the same `k=0`
recombination as in Section 11.445 gives one formula on the closed centred
cell:

```text
C_A,epsilon(s)
 =hat h_A,epsilon(s)+E_A,0(s)g_0(s)+E_A,1(s)g_1(s)
  +sum_(k in Z,k!=0)rho_A,2(k+s),       |s|<=1/2,    (EA6)

g_0(s)=[csc(pi*s)-1/(pi*s)]/(2i),
g_1(s)={1/(pi*s)^2-cos(pi*s)/sin(pi*s)^2}/4,
g_0(0)=0,       g_1(0)=1/24.                        (EA7)
```

Three integrations by parts on the two exterior intervals give

```text
|rho_A,2(xi)|
 <=K3_A/(2*pi*|xi|)^3,
K3_A=|q_epsilon''(c)|+|q_epsilon''(d_A)|+20*pi*epsilon,

sum_(|k|>K)|rho_A,2(k+s)|
 <=K3_A/[(2*pi)^3(K-1/2)^2],       |s|<=1/2.        (EA8)
```

At `K=64`, the fresh Arb enclosures are

```text
epsilon=1e-6: [6.548302940258891799542952621948798162414800652753522509695882672469211e-11 +/- 4.54e-81]
epsilon=1e-8: [6.887479549317750506553120495590115032220669962120633970286253021548241e-13 +/- 1.58e-83]
epsilon=1e-10: [6.910696904212682533221616965764740177706638690339591064842848614155363e-15 +/- 2.09e-85]. (EA9)
```

The 42.25-mode edge displacement is useful, but it is not a nonstationary
certificate.  The surviving exact A common phase is

```text
Phi_A(x)=pi*A^2*x/4+(t/2)log((1-x)/x),
x_*=[1-sqrt(1-8*t/(pi*A^2))]/2
   =[0.49947538036472974127653905997513611415517135950293767675876951565355152221205792 +/- 2.86e-81].                          (EA10)
```

Its positive curvature is `[83939326.461260164873455052806884843660556268898598797452281582630954798389243676 +/- 2.27e-73]`.  Neither `Phi_A`
nor `x_*` contains `d_A`.  Moving the projector edge out of the A
half-boundary collar therefore does not remove this interior stationary
point.  The upper-edge current and the already-deleted A window must be
estimated in one phase-adapted Morse/Fresnel assembly before norms; blanket
nonstationary integration by parts is invalid.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_pole_free_kernel_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_pole_free_kernel_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_pole_free_kernel_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_pole_free_kernel_gate.py
```

Pi provenance: every `pi` in (EA1)--(EA10) is inherited from the Gaussian
Abel/Fourier kernel and the exact Kummer common phase.  No fitted or geometric
occurrence is introduced.

Proof boundary: exact fixed-regulator extended-notch transform, pole-free
cell recombination, cubic dual tails, and common-phase stationary guard only.
No Kummer-weighted cell integral, joined phase-adapted remainder estimate,
quantitative `R_after_A`, `R_Dir`, or `Q_K-T` bound, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
