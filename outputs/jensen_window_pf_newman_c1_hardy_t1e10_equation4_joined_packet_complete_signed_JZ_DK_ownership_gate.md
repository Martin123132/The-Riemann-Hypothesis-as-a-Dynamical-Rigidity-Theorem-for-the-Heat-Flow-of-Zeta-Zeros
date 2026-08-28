# Complete signed production-height J_Z and D_K ownership gate

Date: 2026-08-28

Status: the complete saved-height packet, `J_Z`, and `D_K` are rigorously
enclosed.  The previously selected sufficient corridor is rigorously missed;
this is a route result, not a proof of RH.

The exact preprojection ownership identity is

```text
H(t)J_Z=Hardy_t[V_L+C_U+E_s*T_U+O_join+T_A^join],
Hardy_t[X]=2 Re[exp(i*theta(t))X].                    (FJ1)
```

Here `theta(t)=Im log Gamma(1/4+it/2)-t log(pi)/2`, evaluated directly by
Arb and overlapped with the independent saved Hardy calibration.  The natural
A lift is not a sixth top-level packet: it is subtracted exactly once inside

```text
T_A^join=transition source-transition Gamma lift-mathcal A_A^nat. (FJ2)
```

The ordinary packet is independently reconstructed as
`-T_lower-T_upper+Gamma_defect`, and the transition packet is independently
reconstructed from all three terms in (FJ2).  Every non-negligible top-level
single-sign flip is disjoint from the certified result; adding the natural A
lift again or subtracting it twice is also rejected.

The retained Hardy projections are

```text
V_L       [2.26043634183139660885107246158579069780200026464929934629735827630911101e-5 +/- 6.47e-16]
C_U       [-0.0941220401927439895121289775020126153982520385191340545390389174441034423 +/- 7.24e-26]
E_s_T_U   [+/- 7.96e-1873216240]
O_join    [0.000459381166333759201451968371424248350612509327217104202355670835494171569 +/- 9.68e-7]
T_A_join  [0.127084238678725938583621064863264533056792859443750353830162619680957083 +/- 5.99e-12]
```

Their signed sum gives

```text
J_arg = [0.0287068752306179038632993047580633888729410624413696899125614678871038433 +/- 3.45e-7]
       +i[-0.00904238270971239170145056348074849376027109386128693320715099719148613238 +/- 3.45e-7],

H J_Z = [0.0334441840157340222390325664572920239161313502544798964869423466551109038 +/- 9.68e-7],
J_Z   = [0.0334441839950159192085266113281250000000000000000000000000000000000000000 +/- 9.68e-7].                    (FJ3)
```

Using the independently certified physical `A_transition`, exact `H`,
`rho_RS=U-T`, and the rigorous absolute Gamma-to-classical correction gives

```text
R_KGamma = [-0.00322994080275772750943947977079551201240542597941180774312003822088929339 +/- 9.68e-7],
Q_K-T    = [-0.00322994080275772750943947977079551201240542597941180774312003822088929339 +/- 9.68e-7],
Delta_KU = [-0.00322993949548862253501230458193635938153939727198319914708435390350382801 +/- 9.68e-7],
D_K      = [0.00322993949548862253435783132407150023158920987945732675036270185246516652 +/- 9.68e-7].                 (FJ4)
```

The two algebraically exact evaluations
`D_K=(1-H)U-H*Delta_KU` and `D_K=U-H[T+(Q_K-T)]` overlap.

The stronger sufficient route does not close:

```text
|J_Z|<0.02939975                         false,
-0.0663408<R_KGamma<-0.00727437          false,
0.0072743687<D_K<0.0663407986            false.     (FJ5)
```

The rigorous misses of the relevant walls are respectively
`[0.00404346658369153738021850585937500000000000000000000000000000000000000000 +/- 1e-79]`,
`[0.00404346178288742589423558746771414298637387089558819225687996177911070661 +/- 3.71e-75]`, and
`[0.00404346179013699094408383315932669884983168855804267324963729814753483348 +/- 2.91e-75]`.  The dominant interval
uncertainty is the `O_join` row,
but its width is far smaller than the miss.  This is not a failure due to
enclosure width.  Sharpening the same numerical
enclosures cannot reverse (FJ5); the next route must use additional signed
structure or a different sufficient theorem.

Pi provenance: every `pi` comes from the inherited Fresnel/Fourier phase,
the standard Riemann--Siegel theta normalization, the Mellin quarter turn, or
`p=t/(2*pi)`.  No circle, polygon, fitted constant, or visual pattern supplies
`pi`.

Proof boundary: rigorous final signed assembly and scalar enclosures at the
single height `t=10^10`.  This rejects only the displayed sufficient corridor.
It does not disprove the original one-sided route, prove the non-A bound,
certify equation (4) as a global infinite-series identity, give an all-height
theorem, prove `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
