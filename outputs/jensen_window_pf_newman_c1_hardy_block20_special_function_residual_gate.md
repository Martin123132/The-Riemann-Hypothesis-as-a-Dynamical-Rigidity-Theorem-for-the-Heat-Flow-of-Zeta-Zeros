# Hardy block-20 rigorous special-function residual gate

Date: 2026-08-06

Status: rigorous finite-point PSI/ERF residuals and correlated q replacement; not a proof or uniform q theorem

## Admitted source observations

The source-derived one-CPU probe copied the accepted initialization and the
actual `PSI` and `ERF` routines.  At all 374 recursive block-20 calls it replayed
the saved binary128 `t1`, `t2`, and `t4` components exactly:

```text
374 recursive calls,
2244 real PSI evaluations,
 748 complex ERF evaluations,
1122 / 1122 saved component pairs reproduced bit-for-bit.
```

The rigorous references are Arb `digamma(x)` and
`erf((1-i)x/sqrt(2))`.  The latter is exactly the source target
`erf(exp(-i*pi/4)x)` because `exp(-i*pi/4)=(1-i)/sqrt(2)`.  The evaluator's
ordinary constant is explicitly `p=4*atan(1)`, followed by `tpp=2*p`,
`sp=sqrt(p)`, and `EPI4=(1+i)/sqrt(2)`; no unexplained pi is inserted.

## Correlated PSI cancellation

Number the six source calls in execution order by `P1,...,P6`, let `E` be the
saved endpoint and `T=tpp`.  Direct substitution into the source formulas gives

```text
conj(t1)_PSI + (t5)_PSI
  = (i/T) * ( E*(P6-P1) + (P3-P5) ).
```

The `P2=PSI(1-fracL)` and `P4=PSI(L+1-a1)` residuals cancel exactly.  The gate
therefore propagates the four surviving PSI residuals and the two conjugated
ERF residuals as one complex ball for each call.

## Finite physical result

```text
maximum PSI residual             <= 1.98855351952964895837508918355876914694131421988386E-9
maximum residual difference P6-P1 <= 1.99599455505460966467287296144218592858639951227871E-26
maximum residual difference P3-P5 <= 3.70263062189466654805881746449025632476422897008803E-32
maximum complex ERF residual     <= 1.71706753132667998469767895474729937998208094340512E-4
maximum correlated PSI q shift  <= 3.17672401963412474910334025042163526090576033362147E-27
maximum correlated ERF q shift  <= 7.15521887106601530132103123337381519254144071900596E-4
maximum total special q shift   <= 7.15521887106601530132103123337363958369685514054674E-4
median total special q shift    <= 3.04040101775151785307800603377420549440132875776056E-9
calls with total shift > 1e-4  = 24
```

Every bound is evaluated at the exact saved binary128 physical argument and
contains the corresponding rigorous special-function value.  Because every
recursive chain has one q step and its final source transformations are affine
conjugation/subtraction, the same magnitude bounds the local parent-state change
from this special-function replacement.

This is not yet the full q error.  The real intrinsic `erfc` calls in `t5`, the
Euler-Maclaurin/saddle truncation, binary128 roundoff outside the observed
special-function outputs, recurrence and block accumulation, the Legendre tail,
and the outer Hardy representation remain separate obligations.  The result is
finite block-20 evidence, not a height-uniform theorem or RH.  No prize-level conclusion follows.
