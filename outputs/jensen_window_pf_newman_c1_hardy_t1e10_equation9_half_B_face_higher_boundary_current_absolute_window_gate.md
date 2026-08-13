# Higher pole-subtracted B face currents on the Gaussian window

Date: 2026-08-13

Status: exact boundary-current algebra plus saved-height absolute enclosure;
not a proof of the complete B-face expansion

For the exact B endpoint triangle put

```text
phi(z)=pi*m^2/z+pi*B^2*z/4,
h(z)=z^(-1/2)(2+i*pi*B^2*z).
```

Repeated integration by parts with
`g_0=h/(i*phi')`, `g_(j+1)=g_j'/(i*phi')` produces alternating face
currents `+g_0,-g_1,+g_2`.  After `c=Bx/2` and summation over all positive
modes except 621 and 622, define

```text
R_k(c)=sum_(m>=1,m notin {621,622})(c^2-m^2)^(-k).  (HC1)
```

The next two exact currents after Section 11.373 are

```text
C_1=c/(pi^3 B){(-4pi Bc^3+4i c^2)R_3
                 +(5pi Bc-3i)R_2},                  (HC2)

C_2=-i c^2/(pi^4 B^2){(-48pi Bc^5+48i c^4)R_5
                 +(84pi Bc^3-60i c^2)R_4
                 +(-35pi Bc+15i)R_3}.               (HC3)
```

These formulas retain the signed pair before lattice summation.  Since the
dangerous first current has already been integrated with its trace phase,
the two higher currents may now be bounded absolutely.  On `|xi|<=70`, a
finite Arb sum through `100000` plus the analytic tail

```text
sum_(m>M)|m^2-c^2|^(-k)
 <=rho^(-k) M^(1-2k)/(2k-1),
rho=1-[c_high/(M+1)]^2,                               (HC4)
```

gives

```text
sup |C_1| <= [108.791743189643624417807046156222555142925463094883408685469386337494124669230297588638226 +/- 8.25e-83],
sup |C_2| <= [0.0147921235929150805379404839787615261291101718944964638875742865475547584764039636320681972 +/- 2.03e-86]. (HC5)
```

Restoring the complete Kummer weight, window width, real projection factor,
and `(pi/(32t))^(1/4)` normalization yields

```text
E_1 <= [1.48312944326356763788887383542121368739711210375126858851942741469182667463787421697127152e-6 +/- 1.13e-90],
E_2 <= [2.01657160606416907872098317688713309893620030452649621808249063339158522508621348453446089e-10 +/- 2.76e-94],
E_1+E_2 <= [1.48333110042417405479674593373890240070700572378172123814123566375516583316038283831972496e-6 +/- 1.13e-90]
          <1.491e-6.                                  (HC6)
```

The `C_2` cost is negligible; the deliberately coarse absolute enclosure of
`C_1` is below `1.49e-6`.  No oscillatory gain is claimed for either.

This does not certify the remainder after the third face current, the exact
local replacements, or either outside-window tail.  It also must not be
added to the first-current complex modulus: the certified first-current
quantity is its signed physical projection.

Pi provenance: all pi factors are differentiated from the equation-(9)
triangle phase or inherited from its paper normalization.  No fitted
constant is used.

Proof boundary: the second and third pole-subtracted B face currents on the
saved-height `|xi|<=70` window only.  No complete boundary expansion,
complete B estimate, A-fold splice, complete paired residual, complete
`T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
