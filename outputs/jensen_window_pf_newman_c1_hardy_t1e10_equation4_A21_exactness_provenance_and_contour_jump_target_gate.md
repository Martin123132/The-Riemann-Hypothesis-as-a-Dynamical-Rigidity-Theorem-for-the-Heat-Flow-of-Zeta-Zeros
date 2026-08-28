# Equation-(4) A21 provenance and exact contour-defect target

Date: 2026-08-27

Status: provenance gap isolated; exact split-contour theorem open

The finite odd-label transform `Q_K` remains an exact finite mathematical
object, and the repaired Euler-Kummer identity still maps the complete finite
source roster to it.  What is not currently licensed is the further statement
that the infinite Kummer construction is an exact representation of Hardy's
`Z(t)`.

The 2026 paper calls equation (4) an exact formulation and sends the reader to
the 2015 Appendix A for its derivation.  That appendix starts from the exact
Riemann-Siegel contour

```text
z=1/2+q exp(-i*pi/4),     -infinity<q<infinity,       (CP1)
```

but uses different geometric expansions on the two open half-contours.  Put

```text
w=pi q exp(-i*pi/4),     1/sin(pi z)=1/cos(w).
```

Then the convergent identities are

```text
q<0:  1/cos(w)=2 sum_(k>=0)(-1)^k exp(i(2k+1)w),    (CP2)
q>0:  1/cos(w)=2 sum_(k>=0)(-1)^k exp(-i(2k+1)w).   (CP3)
```

These are the source paper's A13 and A18 expansions.  A13 is not convergent
on `q>0`, A18 is not convergent on `q<0`, and neither ordinary series converges
at `q=0`.

For `N` terms the exact residuals are

```text
q<0: R_N^-(q)=sec(w)[-exp(2iw)]^N,
q>0: R_N^+(q)=sec(w)[-exp(-2iw)]^N.                  (CP4)
```

The Arb rows in the machine artifact verify (CP4) directly at six points and
two precisions.  More importantly, convergence is not uniform at the contour
crossing.  On `q=-u/N`,

```text
|[-exp(2iw)]^N|=exp(-sqrt(2) pi u),                 (CP5)
```

which is independent of `N`.  For `u=1` the certified common value is
`[0.01176198053138912168914475002296359047769707872905960066692824630593644201133705395720180104299821895 +/- 3.55e-102]`.  Pointwise decay away from
`q=0` therefore does not by itself justify exchanging the whole-contour
integral, infinite sum, and limit.

Appendix A acknowledges this issue.  It says A13 is only convergent on its own
half-contour, reports that matching the two termwise asymptotic constructions
fails, introduces A21 as a heuristic continuation of A13 across the entire
`q` range, and sets the convergence question aside.  A33 is subsequently
conditional on suitable convergence criteria.  Later numerical and
asymptotic agreement is evidence for the proposed construction, but is not a
proof of the missing interchange or boundary term.

The rigorous replacement is now precise.  Starting from the exact contour,
split at `q=0`, retain A13 only on `q<0` and A18 only on `q>0`, introduce finite
cutoffs before interchanging sum and integral, and analyze the `q=O(1/N)`
boundary layer.  The chart change at `q=-1/sqrt(2)`, where `Re(z)=0`, must also
retain its branch convention.  After matching the finite roster, the resulting
projected complement must be proved to equal

```text
C_K=Z-H Q_K,       D_K=C_K-L=U-H Q_K.               (CP6)
```

At the saved height the sufficient scalar target is

```text
0.0072743687<D_K<0.0663407986.       (CP7)
```

The source of every `pi` in (CP1)--(CP5) is explicit: it comes from the
Riemann-Siegel denominator `sin(pi z)` and the fixed contour rotation
`exp(-i*pi/4)`, not from an inserted geometric circle constant.

Proof boundary: this gate certifies the two valid geometric half-contour
identities, their exact finite residuals, the nonuniform crossing layer, the
source-provenance nonpromotion, and the inherited sufficient `D_K` target.
It does not derive the contour defect, prove that a single boundary layer is
its only contribution, enclose `C_K`, `D_K`, `Delta_KU`, `Q_K`, or `J_Z`, or
prove a non-A bound, an all-height theorem, `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
