# Compact exact-minus-affine transformed amplitude on A

Date: 2026-08-14

Status: rigorous compact nonlinear-amplitude value on all 84 exact A domains;
not yet a complete A endpoint block or global residual bound

The exact transformed amplitude factors as

```text
A_exact(P,S)=E(P)O(S),
E=2u(pi A m u-i)/[(u+1)(pi A m-i)],
O=v^(3/4)sigma/(v-1).                                (NA1)
```

For `x=1/[1+r(1+y)]`, the endpoint-face variable obeys

```text
P_A=sqrt(pi/2)[A sqrt(x)-2m/sqrt(x)],
I_E(P_A)=integral_(-infinity)^P_A E(P)exp(iP^2/2)dP
        =(-1)^m H_m(x)/C_E,
C_E=2sqrt(2)i(pi A m-i)/(sqrt(pi)A).                 (NA2)
```

Equation (NA2) follows from the exact identities
`C_E E(P_A)dP_A/dx=x^(-1/2)(2+i pi A^2 x)` and
`P_A^2/2=pi m^2/x+pi A^2x/4-pi A m`; `A=159577` is odd.
The previously certified inverse-Gaussian `erfc` primitive supplies `H_m`.

Thus the compact two-dimensional amplitude defect reduces exactly to

```text
R_A,m^amp=(2pi)^(-1) integral_(y_half)^(0.0037) exp(-iS(y)^2/2)
 [O(S(y))I_E(P_A(y))-I_aff(P_A(y),S(y))] S'(y)dy.    (NA3)
```

For rigorous complex-disk callbacks, each Fresnel current is evaluated at
the disk center and its variation is enclosed by the supremum of its exact
elementary ODE derivative times the disk radius.  This avoids a numerically
pathological interval-`erfc` cancellation without changing the function or
weakening the enclosure.

After restoring the A orientation, common carrier, paper rotation, and
equation-(9) normalization, deterministic summation before norms gives

```text
target-side compact nonlinear = [2.5696412105123670413267277021248301385408441634547664560989986973e-9 +/- 1.16e-14],
outer-side compact nonlinear  = [8.4926550509974915531656641478982480707752766845599557033959838571e-9 +/- 1.03e-14],
complete compact nonlinear    = [1.1062296261509858594492391850023078209316120848014722159494977918e-8 +/- 2.19e-14],
sum of modewise moduli        = [1.1426783870400914659757019233280496312673763748326855294218978464e-8 +/- 2.19e-14],
signed/modulus ratio          = [0.9681023457087576389312744140625000000000000000000000000 +/- 3.77e-6]. (NA4)
```

The production calculation is resumable through a modewise fsynced JSONL
cache.  The independent checker recomputes every mode at higher precision,
with more logistic-series terms and a tighter quadrature tolerance.

Pi provenance: every `pi` in (NA1)--(NA4) comes from the original equation-(9)
triangle, exact bi-Morse/Fresnel maps, odd-endpoint carrier, and paper
normalization.  No circle, polygon, fitted constant, or plotted symmetry is
inserted.

Proof boundary: exact-minus-affine transformed-amplitude value on the compact
exact A domain `y_half<=y<=0.0037` for modes 39853..39936 at `t=10^10`
only.  The affine compact value and full exact exterior are separate certified
inputs.  Their complete A assembly, the positive full-line projector jump,
`R_Dir`, complete `Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, and any prize-level conclusion remain unproved.
