# Hardy block-20 W2--W5 source/paper correspondence gate

Date: 2026-08-07

Status: exact orientation ledger and rigorous finite published-component evaluation; not a proof of an exact nonsaddle-contour theorem or RH

## Term ledger

The fixed paper is Lewis--Brereton, equations (89), (93), (94), (96), and
(97).  The source evaluates the printed W terms internally with the conjugate
of its saved endpoint, then the final `qq` assembly conjugates `t1`, `t2`, and
`t4` back into the source's `exp(-2*pi*i*f)` recurrence orientation:

```text
paper endpoint half-sum        (1+E)/2       source t3
paper W2                       equation 89   source conj(t2), after deleting the source-only L+1 term
paper W3 if Phi1<0             equation 93   source conj(t4), negative branch
paper W4 if Phi1>0             equation 94   source conj(t4), positive branch
paper W5                       equation 97   source conj(t1)
```

The source emits `qq=conjg(t1+t2+t4)+t3+t5`; hence the conjugations above are
part of the implementation orientation, not an inferred sign adjustment.
The previously certified `t2` deletion removes only
`i*E/[tpp*(L+1)]` from final `qq`.

## Published components

For each exact binary128 cubic parent, this gate reconstructs

```text
xi       = Phi1 + 2 Phi2 N + 3 Phi3 N^2,
delta    = ceil(xi)-xi,
E_src    = exp(-2*pi*i*f(N)),
Q_paper  = source-oriented conjugate of [(1+conj(E_src))/2 + W2 + (W3 or W4) + W5].
```

Every `erfc` and `digamma` in W2--W5 is evaluated by Arb at 90 and 150
decimal digits.  The selector is derived from exact dyadic `Phi1`.  Exact
dyadic evaluation of `xi=f'(N)` reproduces the saved child-length selector;
its last-bit difference from the Fortran-rounded `fracL` is retained as an
explicit transport term.  The roster contains 165
W3 calls and 209 W4 calls.

Let `P_pi` be the directly summed parent with mathematical `2*pi`, and let
`I69_pi` be the independently certified mathematical-pi equation-(69)
integral family.  The complete published finite residual is

```text
R_paper = P_pi - I69_pi - Q_paper.                     (1)
```

It is independently reconstructed from the prior source-normalized residual
after the `L+1` deletion by transporting the parent, integral family, and
nonsaddle aggregate to mathematical pi.  Both constructions overlap on all
374 calls.

## Finite result

```text
minimum |R_paper|                         >= 1.31747795381777048182240847325351739387979120944877E-6
maximum |R_paper|                         <= 1.53020780143258097321647209044896474499324699793254E-3
improved versus post-L+1 source residual     52 / 374
worsened                                     322 / 374
indeterminate                                0 / 374
R_paper balls excluding zero                 374 / 374
maximum |W2 source-to-paper change|        <= 6.50574187159186715516947816820722907187202393213236E-4
maximum |W3/W4 source-to-paper change|     <= 7.15601558672073567191984226014517346683052089320002E-4
maximum |W5 source-to-paper change|        <= 1.09892857654902906055288227041366128495228754602460E-29
maximum |half-sum source-to-paper change|  <= 7.80548554807693347335500988549631143792600657107606E-32
maximum |complete nonsaddle transport|     <= 7.15521887106601530132103123361733529685432509201032E-4
```

Thus the exact published W2--W5 special-function formulas are now separated
from the source's finite PSI/ERF implementations and from the earlier
source-only denominator.  Any residual in (1) is not licensed as an error
bound: it is the finite discrepancy left by the contour approximations used
to obtain equations (84), (85), (87), (90), and (91), together with any
remaining normalization issue not already excluded by the two-construction
identity.

## Pi provenance

Pi enters only through the Fourier character `e(x)=exp(2*pi*i*x)` used in
the paper.  The source instead stores binary128 `p=4*atan(1)` and
`tpp=2*p`.  This gate does not insert pi from a circle, polygon, or unrelated
geometric assumption: it evaluates both normalizations and explicitly
transports the parent sum, endpoint phase, equation-(69) integrals, and all
W2--W5 factors.

## Boundary

This is a rigorous finite exact-input audit of the published approximation
on 374 block-20 calls.  It does not yet integrate the exact infinite
nonsaddle families (67b)--(67c), prove a height-uniform contour remainder,
control all recurrence levels or the outer Hardy representation, or prove a
determinant/current sign, `Lambda<=0`, PF-infinity, RH, or a prize-level
conclusion.
