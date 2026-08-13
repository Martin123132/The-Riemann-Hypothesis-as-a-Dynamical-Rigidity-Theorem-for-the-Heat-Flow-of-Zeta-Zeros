# Hardy block-20 exact equation-(69) saddle-family gate

Date: 2026-08-06

Status: rigorous finite source-normalized integration of the exact equation-(69) saddle family; not a proof of a height-uniform recurrence estimate or RH

## Exact finite object

For every recursive block-20 call, the parent has length `N=104`, the
transformed length is `L=1` or `2`, and the exact starting index is
`jbot=ceil(Phi1)` in `{0,1}`.  With

```text
g_n(y)=n y-Phi1 y-Phi2 y^2-Phi3 y^3,
```

this gate encloses every entire-function integral

```text
I_n=integral_0^104 exp(i*tpp*g_n(y)) dy,
jbot <= n <= L.                                        (1)
```

There are 544 integrals over the 374 calls.
Each is evaluated independently at 70 and 110 decimal digits by Arb's
certified complex integrator, with one active FLINT thread.  The two precision
levels overlap for every integral and every assembled quantity.

Let `M*C` be the saved recurrence multiplier times the independently summed
and exactly adapted child, and let `t5_src` be the emitted source component.
The source-relative exact saddle-family correction and replacement are

```text
W69_exact = sum_n I_n-M*C,
Delta69   = W69_exact-t5_src.                           (2)
```

No fitted compensation is used.  If `P` is the independently summed parent
and `q_src` the emitted full correction, the residual after replacing `t5`
has the two exactly equivalent constructions

```text
R_after = [P-M*C-q_src]-Delta69
        = P-sum_n I_n-(q_src-t5_src).                   (3)
```

The second line cancels the transformed saddle model completely.  The
checker requires the interval difference between the two constructions to
contain zero on all 374 calls.

## Result

```text
maximum |Delta69|                         <= 1.25603104855648916744327506343761262107625955214361E-1
maximum |Delta69|/|source native target|  <= 1.62235878853359039347207443374292348139654368867647E+1
minimum |R_after|                         >= 5.29499689694274065561570293758332289230513919863099E-2
maximum |R_after|                         <= 8.04686548683205039781523627556905785958526359206402E-2
improved / worsened / unresolved calls       229 / 145 / 0
R_after balls excluding zero                 374 / 374
```

Thus the complete finite equation-(69) saddle family, not merely the displayed
finite-`ip` approximation, can be replaced exactly on this roster.  It changes
the discrepancy materially on some calls but does not close any of the 374
recurrences.  The surviving local defect belongs to the other Euler--Maclaurin
pieces represented by `t1`--`t4` (paper W2--W5 and the endpoint half-sum), or
to their source-model identification, rather than to W1 alone.

## Pi provenance

The source-normalized integral uses the exact binary128 constant produced by
`tpp=8*atan(1)`.  The paper integral uses Arb's mathematical `2*pi`.  The
elementary inequality `|exp(iu)-exp(iv)|<=|u-v|`, integrated against an
explicit polynomial majorant for `|g_n|`, proves

```text
maximum source-tpp to mathematical-2*pi integral shift
  <= 3.00657677567342392853813594173836315668006185953958E-30.
```

Every observed interval shift lies inside that independent bound.  Pi here is
the Fourier normalization in `exp(2*pi*i*x)`; no circle or polygon is inserted
into the recurrence.

## Boundary

This is a rigorous finite exact-point replacement of the saddle-bearing
equation-(69) family on one saved block.  It is not a uniform asymptotic bound,
does not yet integrate the nonsaddle families (67b)--(67c), does not certify
W2--W5 or the outer Hardy representation, and proves no determinant/current
sign, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
