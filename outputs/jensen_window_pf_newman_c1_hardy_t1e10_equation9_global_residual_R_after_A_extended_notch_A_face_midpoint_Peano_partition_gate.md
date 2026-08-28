# A-face midpoint-Peano reduction and saddle/endpoint partition

Date: 2026-08-23

Status: exact certificate for the fixed-regulator signed defect, Peano
reduction, and analytic tail/endpoint partition; quantitative derivative and
x-integral bounds remain open

After stripping the exact common A-face phase, put

```text
c(x)=sqrt(2/x),
q_A(y,x)=-(y-A*x/2)c(x),
H(q)=exp(-i*pi*q^2/2)J_-(q),
G_A(y,x)=-1/(i*pi)-y*c(x)H(q_A(y,x)).                (MP1)
```

For every positive Abel regulator define

```text
f_(epsilon,x)(y)=exp(-pi*epsilon*y^2)G_A(y,x).       (MP2)
```

The integers `39895,...,39936` are exactly the midpoints of the 42 unit
cells partitioning `[39894.5,39936.5]`.  Compact-strip Fubini therefore
reassembles the oriented primal-block/minus-translated-strip A-face channel as

```text
E_(42,epsilon)(x)
 =[-sum_(m=39895)^39936 f_(epsilon,x)(m)]
   -[-integral_(39894.5)^39936.5 f_(epsilon,x)(y)dy]. (MP3)
```

For `|y-m|<=1/2`, let

```text
K_m(y)=(1/2)(1/2-|y-m|)^2.                           (MP4)
```

Twice integrating the cell error, with the value and slope terms cancelling
by midpoint symmetry, gives the exact Peano identity

```text
E_(42,epsilon)(x)
 =sum_(m=39895)^39936
   integral_(m-1/2)^(m+1/2)K_m(y)
   partial_y^2 f_(epsilon,x)(y)dy.                   (MP5)
```

The kernel is positive and has the exact masses

```text
integral_cell K_m(y)dy=1/24,
sum_(42 cells) integral_cell K_m(y)dy=7/4.           (MP6)
```

Since the strip is finite and `G_A` is smooth there for `0<x<=1/2`, the
Abel-zero limit may be taken inside (MP5):

```text
E_42(x)=sum_m integral_cell K_m(y)G_(A,yy)(y,x)dy,
|E_42(x)|<=(1/24)sum_m sup_cell |G_(A,yy)(.,x)|.     (MP7)
```

This is the signed cancellation object.  Bounding the 42 modes or the strip
current separately discards (MP5).

The phase-stripped Fresnel ratio obeys

```text
H_q=1-i*pi*qH,
H_qq=-i*pi*q-(i*pi+pi^2*q^2)H,                       (MP8)

G_(A,y)=-cH+y*c^2 H_q,
G_(A,yy)=2*c^2 H_q-y*c^3 H_qq.                       (MP9)
```

At positive regulator the exact differentiated amplitude is

```text
partial_y^2 f_(epsilon,x)
 =exp(-pi*epsilon*y^2)
  [G_yy-4*pi*epsilon*y*G_y
   +(4*pi^2*epsilon^2*y^2-2*pi*epsilon)G].           (MP10)
```

No phase was fitted in this reduction.  Direct quadratic completion gives,
for real `y`,

```text
exp(i*pi*A*y)exp(-i*pi*y^2/x)exp(i*pi*q_A^2/2)
 =exp(i*pi*A^2*x/4).                                 (MP11)
```

The endpoint growth seen in floating telemetry requires a two-region
argument, not a false uniform-smallness claim.  Fix the rational normal
threshold `Q=16` and let

```text
x_16=
 [(sqrt(2*16^2+8*A*39894.5)-16sqrt(2))/(2A)]^2
 =[0.49990287794624744742232537112847999336114843365113812702103708034289544120397343 +/- 4.37e-81].                           (MP12)
```

Monotonicity of `(y-A*x/2)sqrt(2/x)` in both variables gives

```text
0<x<=x_16, 39894.5<=y<=39936.5
    ==> -q_A(y,x)>=16.                               (MP13)
```

The unique tangential saddle on the half-domain is

```text
x_*=[0.49947538036472974127653905997513611415517135950293767675876951565355152221205792 +/- 2.86e-81],
```

and `x_*<x_16`.  On the remaining endpoint layer, whose width is
`[9.7122053752552577674628871520006638851566348861872978962919657104558796026565631e-5 +/- 3.94e-85]`, the tangential phase is uniformly
nonstationary:

```text
Phi_A'(x)>=Phi_A'(x_16)
 =[21263.469682183676802956899864693847369499352382443478492379135748490719738268651 +/- 3.59e-76] >21000.     (MP14)
```

Thus the next rigorous estimate has a forced architecture:

1. On `0<x<=x_16`, use the uniform `|q_A|>=16` normal-tail hierarchy and
   retain the Morse treatment through `x_*`.
2. On `x_16<=x<=1/2`, keep the exact Peano amplitude and use tangential
   integration by parts, where (MP14) is valid.

The route pilot is intentionally nonrigorous.  It found
`|sum G_A(m)-integral G_A|` about
`0.223349382387` at the saddle but
`2413.07730137` at `x=1/2`.  These
numbers select (MP12)--(MP14); neither is used in the proof.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate.py
```

Pi provenance: every `pi` in (MP1)--(MP14) is inherited from the Gaussian
Abel regulator, Fourier character, canonical Fresnel primitive, and exact
Kummer phase.  The rational threshold `16` is an explicit analytic partition
choice, not a fitted occurrence of `pi` or a numerical theorem constant.

Proof boundary: exact fixed-regulator A-face block/strip orientation, Peano
identity, Abel-zero compact-strip passage, normal derivative ODE, continuous
phase cancellation, and tail/endpoint phase partition only.  No certified
`G_yy` envelope, endpoint-layer integration-by-parts remainder, tangential
Morse integral, quantitative `R_after_A`, `R_Dir`, or `Q_K-T` estimate,
all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
