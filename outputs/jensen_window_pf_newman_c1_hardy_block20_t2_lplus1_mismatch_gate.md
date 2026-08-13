# Hardy block-20 `t2` L+1 source/paper mismatch gate

Date: 2026-08-06

Status: exact source algebra plus rigorous finite correction test; not a proof of the remaining nonsaddle estimates or RH

## Exact source contribution

The published W2 formula, equation (89), contains the endpoint denominator

```text
delta = ceil(xi)-xi = 1-fracL
```

but no denominator `L+1`.  The accepted source instead forms

```text
cr3 = special_term - 1/delta - 1/(L+1),
t2  = i*conj(endpoint)*cr3/tpp,
qq  = conj(t1+t2+t4)+t3+t5.                            (1)
```

Therefore the source-only `-1/(L+1)` term contributes exactly

```text
t2_extra = -i*conj(endpoint)/[tpp(L+1)],
q_extra  = conj(t2_extra)
         =  i*endpoint/[tpp(L+1)].                     (2)
```

The two constructions in (2) agree as complex intervals on all 374 calls.
Removing that source term applies `Delta_drop=-q_extra` to `qq`.

## Exact equation-(69) comparison

Import the residual `R69` left after the independently certified exact
equation-(69) saddle-family replacement.  The corrected residual is

```text
R_drop=R69-Delta_drop=R69+q_extra.                     (3)
```

The complete roster gives

```text
minimum |q_extra|                         >= 5.30516476972984452562945877908381208761524385154073E-2
maximum |q_extra|                         <= 7.95774715459476678844418816862571878686315931944912E-2
minimum |R_drop|                          >= 1.31747794172241320377153049743248173173567252284294E-6
maximum |R_drop|                          <= 1.53021049034020299628325364129787469972147321786431E-3
maximum |R_drop|/|R69|                    <= 1.90162305166434140964595186451264973143083675054272E-2
median upper ratio |R_drop|/|R69|            1.36101443880888613746063608537428159269526236495565E-3
rigorously improved calls                    374 / 374
R_drop balls excluding zero                  374 / 374
```

Thus the single published-source mismatch explains at least 98.09 percent of
the prior residual in the worst call by magnitude and substantially more in
the median call.  This is a source correction, not a fitted compensation: its
phase, scale, and `L` dependence were fixed before comparison by (1)--(2).

## Pi provenance

Equation (2) uses the source's exact binary128 `tpp=8*atan(1)`.  Replacing it
by mathematical `2*pi` changes any roster correction by at most
`4.39319033352878311616778359399023300960258702138028E-36`.  Pi is solely the Fourier
normalization in `exp(2*pi*i*x)`.

## Boundary

The gate proves the source algebra and a finite exact-point improvement after
removing a term absent from published equation (89).  It does not prove that
every unavailable implementation omitted the term, does not yet enclose the
remaining W2--W5 nonsaddle approximations uniformly in height, and proves no
determinant/current sign, `Lambda<=0`, PF-infinity, RH, or prize-level result.
