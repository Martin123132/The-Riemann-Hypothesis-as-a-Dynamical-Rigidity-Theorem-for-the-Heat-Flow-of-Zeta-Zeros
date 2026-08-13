# Hardy block-20 t5 intrinsic real-erfc residual gate

Date: 2026-08-06

Status: rigorous finite exact-argument intrinsic-erfc residuals and correlated t5/q replacement; not a proof or saddle-error theorem

## Admitted source observation

The source-derived one-CPU probe reproduces all 374 saved binary128 `t5`
values bit-for-bit.  Its admitted call roster is

```text
374 recursive calls,
539 lower-endpoint intrinsic real-erfc evaluations,
374 upper-endpoint intrinsic real-erfc evaluations,
913 total intrinsic real-erfc evaluations.
```

For every call, this gate starts from the exact emitted binary128 argument `x`
and encloses `Arb.erfc(x) - source_erfc(x)`.  It does not reconstruct `x` from
decimal data.  In the source, pi enters the argument through `p=4*atan(1)` and
`sp=sqrt(p)`; the rigorous reference inserts no additional pi because the
already-generated binary128 argument is the exact comparison point.

## Linear correlated replacement

For each emitted source weight `W_k`, define

```text
delta_k = W_k * (Arb.erfc(x_k) - source_erfc(x_k)),
delta_q_real_erfc = sum_k delta_k.
```

The sum is taken jointly within each recursive call, preserving cancellation.
The source-term association gap `source_term - W_k*source_erfc(x_k)` is measured
separately and is not relabelled as intrinsic-erfc error.

## Result

```text
maximum intrinsic-erfc residual             <= 4.57549191136680692950631607546020951218639195427451E-35
maximum emitted weight magnitude            <= 4.59161739479680543333751670369362214475171438925164E+0
maximum one-call weighted replacement       <= 1.97694940650258454908778340217706180033093304502771E-34
maximum correlated t5/q replacement         <= 1.97694940650258454908778340217708987471483470054228E-34
median correlated t5/q replacement          <= 9.14274357445012316826743355390387067803280737398938E-71
maximum source-term association gap         <= 4.65280702911089086178078076077377321817957979370104E-34
source outputs enclosed by rigorous erfc ball = 0 / 913
source outputs exactly zero                   = 0
```

This isolates the accuracy of the intrinsic real `erfc` calls at the 913 finite
source arguments.  It does not validate the saddle locations, endpoint formula,
finite `ip` endpoint sums, omitted Euler--Maclaurin terms, other binary128
arithmetic, recurrence accumulation, outer Hardy remainder, height uniformity,
RH, or a prize-level theorem.
