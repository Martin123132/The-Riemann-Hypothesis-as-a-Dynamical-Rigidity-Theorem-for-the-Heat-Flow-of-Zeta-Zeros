# Hardy Cubic Child Legendre Adapter Gate

Date: 2026-08-06
Status: exact cubic transformed-child adapter; not a proof and analytic q/tail enclosure open

## Exact Transform

For `F(x)=a1*x+a2*x^2+a3*x^3`, put `y=2*a2` and solve `u=y*X+3*a3*X^2` on the branch `X(0)=0`. Its Legendre dual satisfies

`H(u)=u^2/(2*y)-a3*u^3/y^3+[9*a3^2/(2*y^5)]u^4+...`.

The source's `scheme` constructs the displayed cubic truncation, and `phifactors` is exactly its binomial translation `H3(k-a1)`. The subsequent integer, parity, and orientation operations preserve the mathematical integer-index phase. The untranslated constant is carried by `exp(i*tpp*phi0)` in the recurrence multiplier.

## Denominator

The source kernel denominator is exactly `D(k)=1+h1*u+(6*h2-3*h1^2/2)*u^2`. The gate performs `131` saved checks and `168` synthetic raised-order checks. The minimum saved absolute denominator is `9.996754486533344014627189293769687437028934147E-1`.

## Saved Chains

All `64` chains have `MIT=2` and degree pair `3 -> 3`. The original parent has `[104, 104]` terms under the source's inclusive convention, while the transformed child length is only `[1, 2]`. The raise-order test therefore uses `L(1)`, not `L(0)`.

The largest exact raw-coefficient source gap is `6.745580893151429043377809152436420699849702234E-33`. The largest first-omitted-term indicator is `8.417345991236922180305484176349741525391016322E-6`, below the source threshold `0.001`.

## Rigorous Point Audit

The maximum exact Legendre-tail ball over all saved child indices is `[8.42281119291988804337324269552774205384748176633e-6 +/- 4.61e-54]`. The maximum weighted-kernel/source gap is `[5.15116150894070182771202768142436005749122805483e-34 +/- 3.59e-83]`. These are finite exact-input point statements, not interval theorems.

## What Remains

The 0.001 first-omitted-term test is a heuristic, not a tail theorem. Source formulas above quartic, selector-stable interval transport, analytic W1-W5/qq, and the outer Hardy remainder remain open.

Turn the exact cubic stationary branch into an outward-rounded interval cell enclosure, then bound the full Legendre tail and analytic qq on that same selector-stable domain.

## Proof Boundary

This gate proves the exact mathematical cubic level-one-to-level-two Legendre/parity/orientation adapter, the exact centered denominator identity, and finite exact-input point reconstructions for sixty-four saved chains. It does not turn the source's first-omitted-term criterion into a rigorous tail bound, audit all source formulas above quartic, enclose analytic W1-W5 or qq on a real cell or disk, control selector transitions or the outer Hardy representation, evaluate a physical carrier, prove a determinant sign, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.
