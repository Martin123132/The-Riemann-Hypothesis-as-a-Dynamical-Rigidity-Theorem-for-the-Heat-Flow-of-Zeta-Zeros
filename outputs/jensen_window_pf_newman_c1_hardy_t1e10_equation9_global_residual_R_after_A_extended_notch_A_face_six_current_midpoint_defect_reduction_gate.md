# Six-current elementary A-face midpoint defect

Date: 2026-08-23

Status: exact certificate for the six-current rational-log reduction and its
uniform tail error on `0<x<=x_16`; tangential integration remains open

Write `delta(y,x)=y-A*x/2`.  Repeated integration by parts in the lower
Fresnel tail gives, for `q<0`,

```text
H(q)=exp(-i*pi*q^2/2)J_-(q)
 =sum_(n=0)^5 (2n-1)!!/[(i*pi)^(n+1)q^(2n+1)]
  +R_(H,6)(q),

|R_(H,6)(q)|<=2*11!!/[pi^7|q|^13].                 (SC1)
```

The final bound exposes the next boundary current and bounds the remaining
`v^-14` integral absolutely; it is a finite remainder theorem, not a formal
asymptotic series.

Substitution into the phase-stripped A current gives the elementary six-term
amplitude

```text
G_6(y,x)=A*x/[2*i*pi*delta]
 +sum_(n=1)^5
  y(2n-1)!!(x/2)^n/[(i*pi)^(n+1)delta^(2n+1)].       (SC2)
```

On the tail region of Section 11.453,

```text
0<x<=x_16,       39894.5<=y<=39936.5,
```

the exact pointwise error obeys

```text
|G_A-G_6|
 <=2*11!!*y*(x/2)^6/[pi^7 delta^13].                (SC3)
```

The right side increases with `x` and decreases with `y`, so its maximum is
at `(x_16,39894.5)`.  Arb certifies

```text
sup |G_A-G_6| <= [1.2196376196424647064542115065855113548487968664640424386524084086007424206468328e-10 +/- 1.09e-90]. (SC4)
```

Define the un-oriented midpoint defect

```text
D_6(x)=sum_(m=39895)^39936 G_6(m,x)
       -integral_(39894.5)^39936.5 G_6(y,x)dy.       (SC5)
```

Every strip integral is explicit.  With `a=A*x/2`,
`delta_l=39894.5-a`, `delta_r=39936.5-a`, and
`b_n=(2n-1)!!(x/2)^n/(i*pi)^(n+1)`,

```text
I_0=A*x*log(delta_r/delta_l)/(2*i*pi),

I_n=b_n
 [delta^(1-2n)/(1-2n)-a*delta^(-2n)/(2n)]
       |_(delta_l)^(delta_r),       1<=n<=5.         (SC6)
```

Thus `D_6=sum_m G_6(m,x)-sum_(n=0)^5 I_n` is a finite rational-log
function.  The exact oriented channel from Section 11.453 is
`E_42=-D_42`; consequently

```text
E_42(x)=-D_6(x)+R_(D,6)(x),
|R_(D,6)(x)|
 <=[1.0244956004996703534215376655318295380729893678297956484680230632246236333433395e-8 +/- 4.30e-88].          (SC7)
```

The factor `84` in (SC7) is the 42 point masses plus the strip length 42.
No sampled cancellation is used in this bound.

At the two route landmarks, direct Arb evaluation of the elementary formula
gives

```text
|D_6(x_*)|  =[0.22334938235023055820720098515552603199135573719079839399667939752172129681191217 +/- 3.22e-81],
|D_6(x_16)| =[8.0332514109606178739320304394292574051736813671351340739952946894536599070486274 +/- 2.83e-80]. (SC8)
```

The apparent `1/x` in the inherited Kummer mode response is harmless at the
lower endpoint.  The exact regularized limit is

```text
lim_(x down 0)D_6(x)/x
 =A[sum_m 1/m-log(r/l)]/(2*i*pi)
  +[sum_m 1/m^2-(1/l-1/r)]/[2(i*pi)^2]
 =[2.0955289415791334259241934506916681691983685475424340046961794801727923949295590e-19 +/- 2.97e-99]
  +i*[1.3977644513426248958019629906367633225674813829102786653503993968114233272913631e-9 +/- 2.63e-89]. (SC9)
```

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.py
```

Pi provenance: every `pi` in (SC1)--(SC9) is inherited from the canonical
Fresnel primitive and exact Kummer/Fourier phase.  No fitted occurrence is
introduced.

Proof boundary: exact six-current expansion, elementary rational-log defect,
uniform tail error below `1.025e-8`, and regular lower-endpoint limit on
`0<x<=x_16` only.  No Morse-core integral, endpoint-layer integration by
parts, full signed A-face bound, quantitative `R_after_A`, `R_Dir`, or
`Q_K-T` estimate, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
