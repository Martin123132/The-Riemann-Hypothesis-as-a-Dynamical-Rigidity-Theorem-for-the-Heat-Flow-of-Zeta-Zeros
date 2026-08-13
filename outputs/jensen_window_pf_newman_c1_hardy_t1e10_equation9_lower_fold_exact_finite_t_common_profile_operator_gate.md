# Exact finite-t common-profile operator

Date: 2026-08-13

Status: exact transformed-strip Fourier operator proved; uniform profile norm
remains open

Keep the common source carrier suppressed and put

```text
t*=beta^3=pi C^2/8,       t=t*-beta lambda,
sigma=4beta/(pi C),       Y=2/sigma,
h=4beta/C,                hY=2pi,
x(z)=[1-tanh(z/beta)]/2.                              (FP1)
```

For `m=39895+n`, exact arithmetic gives

```text
d_m=h(n+3/4),
exp(-i d_mYs)=exp(-2pi i ns)exp(-3pi i s/2).          (FP2)
```

Before the mode factor is inserted, define the exact reduced fold profile

```text
J_ex(lambda,s)=Integral_Gamma cosh(z/beta)^(-3/2)
 (1+sigma Ys/C)
 exp(i{beta^3[z/beta-tanh(z/beta)]-lambda z
        -beta Ys tanh(z/beta)+x(z)Y^2s^2/(2beta)})dz. (FP3)
```

`Gamma` denotes the same admissible exact fold contour/Abel continuation as
the completed-strip theorem.  The beta-minus-four Airy contraction is

```text
J_4(lambda,s)=2pi exp(iY^2s^2/(4beta))
 [U(lambda,Ys,beta)Ai(-lambda-Ys)
  +V(lambda,Ys,beta)d_X Ai(-lambda-Ys)].              (FP4)
```

The inverse Airy normalization and `dy=Y ds` force

```text
g_ex,lambda(s)=Y/(2pi)e^(-3pi i s/2)J_ex(lambda,s),
g_4,lambda(s) =Y/(2pi)e^(-3pi i s/2)J_4(lambda,s),
Y/(2pi)=1/h.                                          (FP5)
```

No fitted scale appears.  Substitution of (FP2) into the exact transformed
mode integral proves

```text
G_ex,(39895+n)=Integral_0^1 g_ex,lambda(s)e^(-2pi i ns)ds,
G_4,(39895+n) =Integral_0^1 g_4,lambda(s)e^(-2pi i ns)ds. (FP6)
```

Therefore, with `delta g=g_ex-g_4` and `w_39894=0`, the entire weighted
exact-minus-beta-minus-four correction is exactly

```text
Delta W=sum_(m=39696)^40094 w_m(G_ex,m-G_4,m)
       =Integral_0^1 delta g_lambda(s)P_w(s)ds,
P_w(s)=sum_(n=-199)^199 w_(39895+n)e^(-2pi i ns).     (FP7)
```

This proves the operator identity left open in Section 11.397.  It neither
separates modes nor applies absolute values before the finite Fourier pairing.
At `s=0`, (FP3)--(FP5) reduce to the saved completed-strip difference and

```text
|delta g_lambda(0)|
 <1.75842251125259831742031899764141519403286918141964047104483729464309891537402844719769180E-10<1.759e-10.       (FP8)
```

The existing whole-kernel enclosure now gives the unconditional implication

```text
sup_(lambda,s)|delta g_lambda(s)|<=epsilon_profile
  ==> |Delta W|<3.840283e-5 epsilon_profile.           (FP9)
```

What remains is analytic rather than combinatorial: certify a uniform bound
for (FP3) minus (FP4) on the top corridor and `0<=s<=1`, preferably by
integrating their difference on one common contour.  The scalar value (FP8)
does not supply that norm.

Pi provenance: `pi C^2/8` is inherited from the exact Kummer phase; `2pi` in
(FP2), (FP5), and (FP6) follows algebraically from `hY=2pi` and the inverse
Airy Fourier normalization.  No circle, polygon, fitted period, or inserted
geometric constant is used.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.py
```

No uniform exact-minus-beta-minus-four profile estimate, numerical bound for
`Delta W`, complete physical source/initial-data splice, complete `Q_K-T` or
`T_upper`, all-corridor or height-uniform theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
