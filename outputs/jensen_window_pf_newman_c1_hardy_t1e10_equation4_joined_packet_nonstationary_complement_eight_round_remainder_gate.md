# Eight-round enclosures of the nonstationary complementary tails

Date: 2026-08-27

Status: `R_B`, `R_L`, and the lower complementary tail are each enclosed as
rigorous complex balls after eight integration-by-parts rounds, but the balls
are numerically useless at the `D_K` scale.  This certifies rejection of the
unsplit absolute eight-round route.  This is not a proof of the complete
packet, `D_K`, or RH.

For each label put

```text
D_q(y)=q-y-p/y,             phi_q'(y)=2*pi*D_q(y),
h_0(y)=y^(-1/2),            h_(n+1)=d/dy[h_n/D_q].  (ER1)
```

Then the exact finite-cutoff identity is

```text
I_q=sum_(n=0)^7 (-1)^n(i2*pi)^(-n-1)
    [exp(i*phi_q)h_n/D_q]_left^right
    +(i2*pi)^(-8) integral exp(i*phi_q)h_8.         (ER2)
```

Every `h_n` is generated as an exact rational monomial table in
`p`, half-integer powers of `y`, and powers of `D_q^(-1)`.  The level-eight
table has 45
exact monomials over denominator powers 8 through 16.  Thus its summed
integral is absolutely convergent even though the first endpoint sum uses the
root-of-unity Abel convention.

The endpoint terms use the certified Section 11.500 balls before any norm.
For the integrated remainder, 8192 rational slabs are fourth-power clustered
at the smallest-gap endpoint.  On each slab, negative powers of `y` are
maximized at its left edge and the denominator roster is bounded by

```text
sum_(m>=0)(delta+m)^(-k)
 <=delta^(-k)+delta^(1-k)/(k-1),                   (ER3)
```

with the better finite-count bound used for `R_B`.  This gives

```text
R_B endpoint partial = [-4.11660623453488067569578899602897450034401601181896526495947792254178360 +/- 7.98e-72]
                     + i [-4.92471501890956475011118666948927239470383536169649622045165896951825707 +/- 1.66e-71],
R_B remainder <= [30.6901190340873834322525840853840722782824480054156256473544973486669373 +/- 4.25e-72],
R_B ball = [-4.11660623453488067569578899602897450034401601181896526495947792254178360 +/- 30.7]
         + i [-4.92471501890956475011118666948927239470383536169649622045165896951825707 +/- 30.7];

R_L endpoint partial = [21793298.1711958647149052691317208587814160198746629495450753745422006482 +/- 6.23e-65]
                     + i [-5616365.17165100409684276700675302376041783084721386349391637139020699196 +/- 5.35e-66],
R_L remainder <= [50174039.9270250201720705340827501599681734571528805799444189322723891003 +/- 1.14e-65],
R_L ball = [21793298.1711958647149052691317208587814160198746629495450753745422006482 +/- 5.02e+7]
         + i [-5616365.17165100409684276700675302376041783084721386349391637139020699196 +/- 5.02e+7];

R_lower endpoint partial = [-0.000373272062088314702969740925245201897226372473397722497201758732433057794 +/- 5.00e-76]
                          + i [-0.000407687675674933890463366256825595132084247433389667730681352292889794539 +/- 4.77e-76],
R_lower remainder <= [25.6198115131043205624420429553941813562250692480943562746261504120688319 +/- 4.16e-71],
R_lower ball = [-0.000373272062088314702969740925245201897226372473397722497201758732433057794 +/- 25.7]
             + i [-0.000407687675674933890463366256825595132084247433389667730681352292889794538 +/- 25.7].
```

The obstruction is structural, not a failed sign check.  With
`chi=delta/sqrt(abs(Q'(y)))`, the three endpoint detunings are

```text
chi_RB(B)    = [0.432498723422339150403428103492868814866634270562778970887143458919095667 +/- 3.99e-73],
chi_RL(L)    = [0.160796682591000508614098969767614400871039796005576578416068794194230573 +/- 6.15e-74],
chi_lower(a) = [21.8416864195495674185578030560095043246436482760102056536434628970383403 +/- 4.06e-71].
```

Thus `R_B` and especially `R_L` begin within the exterior quadratic
endpoint-transition scale.  Their enormous high-order endpoint terms are the
expected asymptotic failure of ordinary integration by parts near an excluded
saddle.  They require a grouped Fresnel/erfc endpoint layer before the remote
tail is integrated by parts.  The lower tail has a different geometry and
should be sharpened separately rather than hidden inside the same failed bound.

Pi provenance: every `pi` comes from the inherited Fresnel/Fourier phase,
the derivative `phi'=2*pi*D`, or `p=t/(2*pi)`.  No fitted circle constant is
used.

Proof boundary: rigorous complex enclosures of the three isolated
nonstationary complementary tails at `t=10^10` only.  They have not yet been
joined with `P_240`, the lower cell, transition packet, upper arc, natural A
lift, or Gamma defect.  No complete ordinary or joined packet, `J_Z`, `D_K`,
non-A, all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
