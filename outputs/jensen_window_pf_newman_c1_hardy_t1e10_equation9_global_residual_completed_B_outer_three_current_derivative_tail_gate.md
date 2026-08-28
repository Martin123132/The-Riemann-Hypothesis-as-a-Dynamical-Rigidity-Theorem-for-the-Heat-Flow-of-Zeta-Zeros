# Outer three-current derivative tail for the completed B trace

Date: 2026-08-13

Status: uniform analytic rational-tail derivative certificate; not a proof
of the complete completed-B estimate

Let `c=Bx/2`, `D=c^2-m^2`, and retain the first three paired rational
currents `C_0,C_1,C_2` from the sign-adapted B-endpoint Fresnel expansion.
The exact formulas are

```text
C_0=-(1+i*pi*B*c)/(pi^2 D),

C_1=c/(pi^3 B){(-4*pi*B*c^3+4*i*c^2)D^-3
                  +(5*pi*B*c-3*i)D^-2},

C_2=-i*c^2/(pi^4 B^2){(-48*pi*B*c^5+48*i*c^4)D^-5
                         +(84*pi*B*c^3-60*i*c^2)D^-4
                         +(-35*pi*B*c+15*i)D^-3}.       (OT1)
```

For the analytic outer tail choose `m>=B=5122421`.  On the entire extracted
domain `delta=1e-4<=x<=1/2`,

```text
0<=c<=B/4,       c/m<=1/4,
|D|=m^2-c^2 >= (15/16)m^2.                              (OT2)
```

Each term in (OT1) is a monomial `alpha B^q c^p D^-k`.  It is differentiated
symbolically through order two before absolute values.  Every resulting
lattice tail is then enclosed by

```text
sum_(m=M)^infinity m^(-2k)
 <= M^(-2k)+M^(1-2k)/(2k-1).                            (OT3)
```

Consequently every finite cutoff and every positive Abel weight bounded by
one obeys the following uniform bounds for
`S=sum_(m>=B)(C_0+C_1+C_2)`:

```text
|S|    <= [434804.6836814489523901306578426608984464570483212588771091035388685979543718723043536214626367495692 +/- 9.66e-95]
|S_x|  <= [908258.6876692366761304035452011803398450411947010787036020128401842620227254481919702999535726587044 +/- 1.10e-94]
|S_xx| <= [244263.7091651926897864433397310038733277378464058639595233653698059922863052157368307001616989879347 +/- 6.52e-95]                                    (OT4)
```

After multiplication by `w=[x(1-x)]^(-1/4)` and conversion to
`xi=sqrt(H_B)(x-x_B)`, the corresponding bounds are

```text
|w S|          <= [4348155.544779742646804112504715805885715758606408442383088345961292490854566605480445437547506605638 +/- 9.93e-94]
|d_xi(w S)|    <= [37.33484983165941542016850288662853249998908206656058813678319058896503905927253234601415142108665805 +/- 8.39e-97]
|d_xixi(w S)|  <= [0.001600968749856890216854167359496763059942147428576471262338257433923060153268383356021505369408026206 +/- 7.16e-101].                    (OT5)
```

The large unnormalized derivatives in (OT4) are expected near the artificial
corner `x=delta`; (OT5) records the actual normalized scale needed by the
tangential argument.  No mode enumeration is used.

What this closes: the first three *rational* boundary currents have a
cutoff-uniform, Abel-uniform analytic outer tail through two derivatives.

What remains open: differentiated bounds for the exact Fresnel expansion
remainder, the finite block `m<B`, both cutoff/window edge terms, and the
common-regulator combination with `R_join`.  Thus this gate proves no complete
B estimate, `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH,
or prize-level conclusion.

Pi provenance: every `pi` in (OT1) is inherited from the equation-(9) Fresnel
phase and its endpoint integrations by parts.  No geometric or fitted value
of `pi` is introduced in this certificate.
