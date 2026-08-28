# Joined transition K on the first height subcell

Date: 2026-08-28

Status: rigorous transition-component interval certificate; the complete
`K_T` is not yet enclosed

On

```text
I_1=[10^10-10^-4,10^10+10^-4],
```

write

```text
J_tr=I_tr-C_G D_(39853..39894),
K_tr=J_tr'+(i theta'-H'/H)J_tr.                    (TK1)
```

The finite transition source is differentiated before its alternating
boundary-jet collapse.  Since

```text
d_t F(y,t)=-i log(y)F(y,t),
```

the complete Hardy transport operator has the cancellation-stable source
weight

```text
(d_t+i theta'-H'/H)F
  =[i(theta'-log(y))-H'/H]F.                       (TK2)
```

The 42-mode Gamma block is likewise summed with its operator weight per mode:

```text
C_G sum_m [i(theta'-log(m))+C_G'/C_G-H'/H]m^(-1/2-it).
                                                               (TK3)
```

No separately widened `D'` and `i theta'D` balls enter the direct result.
The source and Gamma operator packets are

```text
source K: real [5.705656817908623508368822281209214431728540332200977863e-5 +/- 2.48e-7]
          imag [4.529518035964684206652038761786063134682230400186764772e-5 +/- 2.55e-7],

Gamma K:  real [3.536019129435206792641993278849844502350919726150049031e-6 +/- 8.69e-8]
          imag [5.153874803093435437454248389480590774084511477870679484e-5 +/- 6.33e-8].
```

Their direct joined enclosure is

```text
K_tr: real [5.352054904965102829104622953324229981493448359585972960e-5 +/- 3.34e-7]
      imag [-6.243567671287512308022096276945276394022810776839147118e-6 +/- 3.18e-7].             (TK4)
```

It overlaps the deliberately dependency-heavy assembly from `J_tr'` and
`(i theta'-H'/H)J_tr`.  After the common Hardy projection and division by
positive `H`, this component contributes

```text
Hardy_t[K_tr]/H=[7.564537099824519827961921691894531250000000000000000000e-5 +/- 1.06e-6] > 0.             (TK5)
```

The independent checker raises precision from 384 to 448 bits, changes jet
order `14 -> 16`, derivative slabs `256 -> 384`, integral panels `48 -> 64`,
and logarithm-series terms `34 -> 40`.  It also repeats the altered five-label
direct finite-source comparison.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
integer Fourier spacing, Gamma normalization, Riemann--Siegel phase, or
`p=t/(2*pi)`.  No fitted geometric constant supplies `pi`.

Proof boundary: (TK1)--(TK5) certify only the transition contribution to
`K_T` on `I_1`.  They do not enclose the lower-cell, upper-arc, ordinary-packet,
or tiny target-correction contributions, hence do not prove a complete
nonzero-radius `K_T` bound, a wider `Q_K-T` sign interval, the full event-cell
sign theorem, an event-wall handoff, an all-height theorem, `Lambda<=0`, RH,
or a prize-level conclusion.
