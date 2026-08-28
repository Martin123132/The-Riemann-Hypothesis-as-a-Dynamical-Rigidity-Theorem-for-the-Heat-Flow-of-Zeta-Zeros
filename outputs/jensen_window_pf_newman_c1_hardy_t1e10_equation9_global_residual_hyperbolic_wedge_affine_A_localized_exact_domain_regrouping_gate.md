# Localized exact affine A-domain and exterior saddle rosters

Date: 2026-08-13

Status: exact regrouping and rigorous exterior-stationarity guard; not an
exterior integral bound or complete A carrier

Let `chi_E` and `chi_T` denote the exact and tangent A domains, and split
`S>h` at `S_0` into `chi_L+chi_X=1`.  Before norms,

```text
chi_T+(chi_E-chi_T)=chi_E=chi_L chi_E+chi_X chi_E.    (LR1)
```

However, retaining only the already certified local face strip gives

```text
chi_T+chi_L(chi_E-chi_T)=chi_L chi_E+chi_X chi_T.     (LR2)
```

Thus global tangent plus local face is not the local exact carrier: it still
contains the entire tangent exterior.  The cancellation-safe next object is

```text
C_aff,exact=C_aff,exact,loc+C_aff,exact,ext,           (LR3)
```

with no tangent exterior in the second term.

This matters because the tangent boundary phase

```text
F_T(S)=((a+rho S)^2-S^2)/2,
F_T'(S)=rho a+(rho^2-1)S                              (LR4)
```

has an artificial exterior stationary point for exactly

```text
39884..39894 (11 modes).   (LR5)
```

The adjacent guards are

```text
mode 39883: S_T*=[256.18469326220139877077991157425505168676834783235825909528499445 +/- 4.31e-63] < S_0,
mode 39884: S_T*=[284.89005353413202879627612989906793171835037326053336398811672213 +/- 1.04e-63] > S_0,
mode 39894: rho^2-1=[6.2665739085848621125804364673098811428099067711555368331915967423e-6 +/- 3.20e-71] >0,
mode 39895: rho^2-1=[-1.8799722628127119247836956436522493028139012880850372956003044308e-5 +/- 1.91e-70] <0.        (LR6)
```

These saddles belong to the unbounded linearization, not the exact triangle.
They must be removed algebraically by (LR3), not bounded as though the
exterior were nonstationary.

After the common carrier is restored, the exact boundary phase on `z=x` is

```text
Phi_face,A(x)=pi A^2 x/4+(t/2)log((1-x)/x),           (LR7)
```

which is independent of `m`.  Its unique half-domain stationary point is

```text
kappa=sqrt(1-8t/(pi A^2))=[0.00104923927054051744692188004972777168965728099412464648246096869289695557588 +/- 4.16e-78],
x_*=(1-kappa)/2=[0.499475380364729741276539059975136114155171359502937676758769515653551522212 +/- 5.80e-77]. (LR8)
```

At the common cutoff `y_0=0.0037`, this true stationary point remains in
the exterior for exactly

```text
39927..39936 (10 modes).          (LR9)
```

The adjacent exact guards are

```text
mode 39926: y_*=[0.0036974740994109080437239055366726113716723866412071427990950150721 +/- 1.18e-68] < y_0,
mode 39927: y_*=[0.0037477526168467323315050806315950745621191843354422418421615966123 +/- 3.44e-68] > y_0. (LR10)
```

Finally, every exact exterior inner endpoint satisfies

```text
P_A(y_0)<-261.2,                                      (LR11)
```

so the inner Fresnel tail can be expanded with an explicit inverse-`P`
remainder.  Equations (LR7)--(LR10) show that its leading boundary currents
must be summed in the common `x` phase and supplied with a ten-mode incomplete
stationary transition; a blanket exterior integration-by-parts estimate is
not valid.

Pi provenance: every `pi` in (LR4)--(LR11) comes from the exact equation-(9)
bi-Morse phase, the triangle face `z=x`, and the saved endpoint normalization.
No fitted constant is introduced.

Proof boundary: exact localization algebra, the complete 84-mode tangent and
exact exterior stationary rosters, common exact-face phase, and cutoff
Fresnel-tail margin only.  No local exact-affine value, exterior Fresnel
current estimate, transformed-amplitude remainder, complete A endpoint
theorem, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
