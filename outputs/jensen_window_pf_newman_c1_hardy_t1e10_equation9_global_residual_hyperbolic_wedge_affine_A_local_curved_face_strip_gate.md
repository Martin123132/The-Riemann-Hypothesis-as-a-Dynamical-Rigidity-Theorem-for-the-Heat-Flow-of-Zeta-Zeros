# Exact local curved-face strip on the A Airy box

Date: 2026-08-13

Status: rigorous finite-box affine curved-face correction on all 84 A modes;
not a complete curved-face, amplitude, A-endpoint, or global residual bound

Put `y=v-1`.  The logistic Morse coordinate is analytic at `y=0` in the
form

```text
sigma=y sqrt(q(y)),
q(y)=2[y-log(1+y)]/y^2
    =2 sum_(k>=0)(-1)^k y^k/(k+2),
S=sqrt(t/2)sigma.                                     (CF1)
```

The evaluator uses `24` terms.  For `|y|<=r<1`, the omitted
tails are bounded directly by

```text
|q-q_N| <=2r^N/[(N+2)(1-r)],
|q'-q_N'|<=2r^(N-1)/(1-r).                            (CF2)
```

At `r=0.0037` their saved bounds are respectively

```text
[3.34586100195286744804739351626593387541596214e-60 +/- 3.18e-105],
[2.35114556893985280133060084926795353407608150e-56 +/- 1.56e-101].            (CF3)
```

This removes the apparent `0/0` at the null point without a floating-point
branch choice.  In the same variable the exact A face is

```text
P_A(y)=sqrt(pi/2){A/[1+r_m(1+y)]^(1/2)
                   -2m[1+r_m(1+y)]^(1/2)}.           (CF4)
```

For `P_lin=a+rho S`, the exact affine strip integral on
`y_half<=y<=0.0037` is evaluated through the analytic primitives

```text
J_0(P)=integral_0^P exp(iu^2/2)du,
J_1(P)=-i exp(iP^2/2),                                (CF5)

R_face,loc=(2pi)^(-1) integral e^(-iS^2/2)
 {(1+mu S)[J_0(P_A)-J_0(P_lin)]
   +lambda[J_1(P_A)-J_1(P_lin)]}dS.                  (CF6)
```

All 84 half-boundaries lie inside this box; its upper edge corresponds to
`S` a little above `261`.  Restoring the exact A sign, common carrier,
paper rotation, and normalization gives

```text
target-side local strip=[0.0038773379822857392867178339380679923575896015411274111704999872516 +/- 2.89e-13],
outer-side local strip =[-0.014095331213899910877810847387545630478321779977922139070879982727 +/- 1.33e-13],
complete signed strip  =[-0.010217993231614171591093013449477638120732178436794727900380000945 +/- 4.21e-13],
sum of modewise moduli =[0.037680742616479709888957651046532718265226196027022780903179899666 +/- 4.21e-13],
signed/modulus ratio   =[0.2711728199100121288064001168205408021094626747071743011 +/- 1.42e-11]. (CF7)
```

The exact local face correction is therefore not a perturbation at the
`1.4e-4` target scale.  It has to be retained as part of the A Airy carrier,
not charged as an error to the affine tangent wedge.  Equation (CF7) does
not include `y>0.0037`, the exact-minus-affine transformed amplitude,
or any of the global projector companions, so it is not an `R_Dir` estimate.

The production run is resumable mode by mode in an fsynced JSONL cache.  An
independent higher-precision replay recomputes every mode with more series
terms and tighter quadrature tolerance.

Pi provenance: every `pi` in (CF1)--(CF7) comes from equation (9), the exact
bi-Morse transformation, Fresnel primitives, odd-endpoint phase carrier,
and paper normalization.  No fitted constant is introduced.

Proof boundary: the exact affine curved-face strip on the finite real box
`y_half<=y<=0.0037` for modes 39853..39936 only.  No exterior face
strip, transformed-amplitude remainder, complete A endpoint theorem,
`R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
