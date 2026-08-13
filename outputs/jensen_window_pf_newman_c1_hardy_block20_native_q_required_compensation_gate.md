# Hardy block-20 native-q required-compensation gate

Date: 2026-08-06

Status: rigorous finite native-q inversion and required-compensation target; not a proof of the missing remainder

## Exact affine inversion

Let `P_raw` be the independently summed parent, `C` the source-adapted child,
`M` the saved multiplier, and `O` the saved conjugation/subtract-one map.  The
source recurrence has the form

```text
P_adapted = O(P_raw),
model(q)  = O(M C + q).
```

Therefore the unique additive q value required for exact finite equality is

```text
q_need = P_raw - M C.                                  (1)
```

The same value is independently recovered as
`O^-1(P_adapted)-M C`.  If

```text
q_full = q_source + delta_q_PSI/complex-ERF
                  + delta_q_intrinsic-real-erfc,
```

then the native correction still required by the finite sums is

```text
R_q = q_need - q_full.                                 (2)
```

Subtract-one cancels from the difference.  Consequently the already admitted
fully corrected recurrence defect is exactly `R_q` when the parent orientation
is direct and `conj(R_q)` when it is conjugated.

## Rigorous result

Every construction in (1)--(2) is evaluated at 180 and 260 decimal digits from
the exact saved binary128 source normalization.  All 374 precision pairs
overlap; both q-need constructions overlap; all orientation identities contain
zero; and all reconstructed defects overlap the prior fully corrected balls.

```text
minimum |R_q|                           >= 4.66665048362406825963434180168989386723463716839527E-3
maximum |R_q|                           <= 2.05015438240021097313750686731496617932304852838501E-1
median  |R_q|                           <= 8.02846043009659936679984490819105861054668326153596E-2
maximum |R_q|/|q_full|                  <= 2.58197579373967616663003920881484451121755203728339E+0
maximum q-need cross-construction gap   <= 4.28697880939147175028515614721148442591181734205048E-258
maximum oriented-defect identity gap    <= 1.72564779004958218308191201568228192600425818814959E-177
nonzero native compensation balls        = 374 / 374
```

This converts the observed recurrence discrepancy into the exact finite
complex target that any omitted endpoint-saddle, finite-`ip`,
Euler--Maclaurin, or compensation term must supply.  It does not show that any
particular omitted source term equals `R_q`, does not assign a sign or formula
to the remainder, and does not provide a height-uniform bound.  It is not a
proof of an outer Hardy estimate, `Lambda<=0`, PF-infinity, RH, or a prize-level
theorem.
