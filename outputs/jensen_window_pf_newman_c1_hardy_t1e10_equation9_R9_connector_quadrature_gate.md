# R=9 characteristic connector quadrature

Date: 2026-08-11

Status: finite R=9 connector quadrature validated; not a proof of a
height-uniform contour-tail theorem or the `y>64` saddle join

Section 11.316 shows that the prior `|z|>4` contribution contains migrating
Airy saddles once `y>16-lambda`.  The stationary core must therefore include

```text
-9<=z<=9,       0<=y<=64.                               (R9Q1)
```

Forty quarter-unit connector panels on `[-9,-4] union [4,9]` give

```text
I_connector=[3.0850800042833600664426731371205705 +/- 1.95e-12]
            +i*[0.017323409590500173273959762950833513 +/- 1.95e-12].           (R9Q2)
```

Combining (R9Q2) with the certified `R=4` core gives

```text
I_canonical(9,64)=[-43.188350903648067576266875197562766 +/- 2.97e-12]
                  +i*[1.4321090565520764095697441099321242 +/- 2.97e-12].  (R9Q3)
```

The exact finite-`t` connector again uses the degree-17 stable tanh form and
closed Fresnel first moment.  Its omitted Taylor contribution to the complete
connector is below

```text
[1.0942975720272002501369707419191559e-28 +/- 6.53e-62] <1e-27.     (R9Q4)
```

The exact finite-`t` expanded rectangle divided by `sqrt(2)/pi` is

```text
I_exact(9,64)=[-43.188859078644004368995839202028855 +/- 2.92e-12]
              +i*[1.4321004698155396763530253889522353 +/- 2.92e-12],       (R9Q5)
```

and its relative correction to (R9Q3) is

```text
[1.1761559886313065640480075846330266e-5 +/- 1.39e-13].             (R9Q6)
```

Subtracting (R9Q3) from the independent full-line Airy representation leaves
the true outgoing `R=9` contour tail

```text
I_tail=[-0.14458228191428833422544111922501544 +/- 2.97e-12]
       +i*[-0.0039190273487292765251016547929787113 +/- 2.97e-12],               (R9Q7)

|I_tail|/|I_canonical(9,64)|
 =[0.0033471047473600350384070819002510500 +/- 7.06e-14].                           (R9Q8)
```

Unlike the old `|z|>4` quantity, (R9Q7) contains no canonical stationary
point for `0<=y<=64`.  It is the term that must now be certified directly on
the outgoing `pi/6` and `5pi/6` rays.

Pi provenance is unchanged: the rays are forced by the cubic Airy phase and
all scale factors retain the Kummer/Fourier--Poisson normalization.

Proof boundary: rigorous finite connector quadrature and a saved-height true
tail measurement only.  No direct ray-integral bound, height-uniform theorem,
`y>64` partition, complete `T_upper`, `Lambda<=0`, RH, or prize-level
conclusion is proved.
