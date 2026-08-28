# Exact H(t) and continuation-defect target

Date: 2026-08-27

Status: exact normalization and target certified; continuation theorem open

Starting from the Kummer prefactor stated in equation (4) and Euler's
integral, the normalization multiplying the correctly scaled finite Kummer
main simplifies exactly to

```text
H(t)=2^(5/4)t^(1/4)pi^(3/2)exp(-3pi t/4)
     /[(1+exp(-2pi t)) |Gamma(1/4+it/2)|
       |Gamma(3/4+it/2)|^2].                          (HC1)
```

The duplication formula gives the independent equivalent expression

```text
H(t)=2^(3/4)t^(1/4)sqrt(pi)exp(-3pi t/4)sqrt(cosh(pi t))
     /[(1+exp(-2pi t)) |Gamma(3/4+it/2)|].            (HC2)
```

Arb evaluates both logarithmic forms at two precisions.  At `t=10^10`,

```text
log H=[3.12500000000000000001953125000000000000039713541666666666668357340494791666666790008951822916666680413640340216861423673e-22 +/- 1.97e-130],
H-1=[3.12500000000000000002001953125000000000040328979492187500001705087423324584961061450731381773948683223310896999761468520e-22 +/- 1.97e-130].       (HC3)
```

The exactness in (HC1)--(HC3) is algebraic exactness of the displayed
special-function prefactor.  It does not certify that the proposed
equation-(4) infinite Kummer series equals Hardy's `Z(t)`; that contour and
interchange theorem remains a separate obligation.

In particular, the entire exact `H-1` correction to the complementary upper
component is below `[6.554826139571993444923852038951940633996691562513911810017914398705962484052407592918304693695534439e-22 +/- 1.04e-122]`.
It is about nineteen orders of magnitude too small to supply the bridge
target.

Define the exact finite-roster continuation complement and its lower-adjusted
defect by

```text
C_K=Z-H Q_K,
D_K=C_K-L=U-H Q_K.                                    (HC4)
```

Since `Q_K=U+Delta_KU`, this gives the exact reversible identities

```text
D_K=(1-H)U-H Delta_KU,
Delta_KU=((1-H)U-D_K)/H.                              (HC5)
```

Transporting the certified `Delta_KU` target through (HC5) yields the simple
conservative target

```text
0.0072743687<D_K<0.0663407986.         (HC6)
```

The complete interval endpoints and the corresponding target for `C_K` are
stored in the machine artifact.  Thus a successful continuation theorem must
show that the exact complement `Z-HQ_K` differs from the independently
certified lower Hardy component by a **positive** amount between about
`0.00727` and `0.06634`.

The paper's Appendix-A Euler-Maclaurin route cannot currently certify (HC6):

```text
A64: asymptotic Kummer term with 1+O(Lambda_alpha^-1);
A65: transition approximation with O(g^2) and a stated restricted range;
A71: derivative formula explicitly only to leading order;
A73: Z(t) is written with an approximation sign;
A74: the tail integral is a Laplace estimate with 1+O(1/t).              (HC7)
```

The surrounding text also labels part of the analytic-continuation
cancellation as an unproven assertion and speculative.  Consequently A73 is
a route hint, not an exact identity or a constant-bearing remainder theorem.
The next derivation must begin from the exact equation-(4)/Riemann-Siegel
contour and produce `C_K` with every omitted lower label, upper continuation,
and transition contribution explicitly owned.

Proof boundary: exact `H(t)` normalization, exact defect identities, and a
saved-height sufficient defect target only.  No enclosure of `D_K`, `C_K`,
`Delta_KU`, `Q_K`, or `J_Z`, no non-A bound, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
