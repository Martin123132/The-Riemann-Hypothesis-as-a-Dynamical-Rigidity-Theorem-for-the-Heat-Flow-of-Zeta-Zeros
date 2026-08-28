# Null-coordinate unfolding of the A wedge

Date: 2026-08-13

Status: exact characteristic-coordinate reduction and saved-height unfolding
certificate; not a numerical wedge or residual bound

The flat-face derivative formula is nonuniform when `c=rho^2-1` crosses
zero.  Near the A corner use instead

```text
R=P+S,       delta=1+rho.                              (NU1)
```

Then

```text
(P^2-S^2)/2=R^2/2-RS,
P<a+rho S  iff  R<a+delta S,
r0=a+delta h.                                         (NU2)
```

Swapping the two half-plane integrations and performing the `S` integral
gives, for `delta<0`,

```text
C=(2pi i)^(-1) int_(-infinity)^r0 exp(iR^2/2)
 [exp(-iRh)-exp(-iR(R-a)/delta)] dR/R.                (NU3)
```

The numerator in (NU3) vanishes at `R=0`.  For `delta>0`, the same common
Abel boundary value is

```text
C=1/2+(2pi i)^(-1)PV{
 int_(-infinity)^r0 exp(iR^2/2-iRh)dR/R
 +int_r0^infinity exp(iR^2/2-iR(R-a)/delta)dR/R}.    (NU4)
```

At `delta=0`,

```text
C=H(a)/2+(2pi i)^(-1)PV
 int_(-infinity)^a exp(iR^2/2-iRh)dR/R.               (NU5)
```

Equations (NU3)--(NU5) are one-dimensional characteristic formulas.  They do
not divide by `c` and expose the half-jump rather than hiding it in a
singular Fresnel scaling.

The continuous half-boundary mode and its dynamic null endpoint are

```text
nu=sqrt(t/(2pi)),       D_t=4nu=sqrt(8t/pi).           (NU6)
```

For endpoint `D`, the codimension-two unfolding is exactly

```text
eta_D=sqrt(pi)(D-D_t)/2
     =4(pi D^2/8-t)/[sqrt(pi)(D+D_t)].                 (NU7)
```

Thus `eta_D=0` exactly at `t=pi D^2/8`.  For the saved A endpoint,

```text
eta_A=[0.07784566526066391096382359515456067508724212627806181942904160576418438086731859649704241921376306123 +/- 7.09e-103].       (NU8)
```

The two adjacent integer modes straddle the null slope:

```text
m=39894: delta=[-3.13328204556424252866822147293402176482672263935920007949562872932369163060000000000000000e-6 +/- 3.79e-80],
         r0-eta_A=[-7.80986384653262648402704839777748652242993340038109351314719176745345060021499138266698244e-12 +/- 3.10e-80],
m=39895: delta=[9.39990549317519993656784885098273694712611441631386115569173700781843121987000000000000000e-6 +/- 2.13e-81],
         r0-eta_A=[2.97628158663382047055126502392450238548738549103509566680847759242686602548684633365482762e-10 +/- 5.90e-81]. (NU9)
```

At the exact continuous null point, `rho=-1`,
`P_A''(0)=-2/(3sqrt(pi)A)`, and the boundary phase has cubic coefficient
`1/(3A sqrt(pi))`.  This recovers the Airy cubic from the triangle geometry
itself and explains the provenance of its `pi`.

The result changes the quantitative route: evaluate the affine A wedge with
(NU3)--(NU5), then compare the exact curved boundary to this null-coordinate
carrier.  It remains essential to retain the paired-triangle amplitude and
the global projector orientation; the older fold atlas cannot be inserted
by matching phases alone.

Proof boundary: exact null-coordinate algebra, one-dimensional Abel formulas,
continuous cubic coefficient, and saved-height corner localization only.  No
canonical wedge value, curved-face or amplitude bound, `R_Dir` estimate,
complete `Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
