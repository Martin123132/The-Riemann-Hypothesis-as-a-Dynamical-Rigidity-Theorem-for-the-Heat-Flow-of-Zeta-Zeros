# Newman C1 Hardy Exact Low-Mode Projected-Kernel Gate

Date: 2026-08-05

Status: `exact_low_mode_projection_with_diagnostic_kernel_replay`; not a proof of the external error theorem, RH, or a prize-level conclusion.

## Exact projection

The twenty-four saved physical coefficient vectors already showed a sparse low-mode pattern, but only to floating precision.  For each row, this gate defines the replacement vector

```text
c_tilde=c-sum_(k in Z) <c,p_k> p_k/||p_k||_2^2,
```

where the `p_k` are the primitive integer Gram vectors and `Z` is the row's listed inactive-mode set.  The saved decimals are finite rationals, so every entry of `c_tilde` and every identity `<c_tilde,p_k>=0` is exact.

Across `24` rows this enforces `84` exact annihilations.  The largest coefficient-vector L2 change is `8.26409961206e-08`.

## Kernel replay

The complete midpoint, Chebyshev, and endpoint held-out grids were replayed at 180-digit phase precision.  The maximum projected-kernel relative error is `9.32193472449e-11`, while the largest projected coefficient L1 norm remains `2311004.4995`. This is a floating diagnostic, not an interval certificate.

## Derivative sensitivity

Substituting the exact moment zeros into Section 11.237 gives the following largest coefficient multipliers after normalization by each target scale:

| derivative order | maximum normalized multiplier | row |
|---:|---:|---|
| 0 | 0.484616617238 | root 1 `base_rows:P_V` |
| 1 | 0.0193846646895 | root 0 `base_rows:P_A` |
| 2 | 0.000928848516373 | root 0 `base_rows:P_V` |
| 3 | 3.10154636108e-05 | root 0 `tangent_rows:P_N` |
| 4 | 6.49847806816e-07 | root 2 `base_rows:P_N` |
| 5 | 2.03323579822e-08 | root 1 `tangent_rows:P_N` |
| 6 | 5.20555700473e-08 | root 0 `tangent_rows:P_N` |

Thus the million-scale independent L1 amplification is not the relevant structured derivative budget.  The exact low-mode design reduces the worst normalized multipliers for orders zero through six to the values above.  It still does not supply the missing external `D_r(T)` bounds.

## Falsification retained

At the two calibration heights, `21/24` absolute projected residuals still increase.  The largest increase factor is `3.48540571632`.  Exact moment cleanup therefore does not revive the rejected scalar-extrapolation route.

No new `pi` is introduced.  The grid is the evaluator's literal `0.01` shift grid, and the projection uses rational finite-grid inner products only.

This gate proves eighty-four exact rational low-mode annihilations and replays finite kernel and two-height diagnostics. It does not prove the external evaluator derivative theorem, an interval kernel approximation, a physical carrier value, retained observation, determinant sign or bound, complete-current inequality, all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.

```text
built Newman C1 Hardy exact low-mode projected-kernel gate: 24 rows, 84 exact annihilations, max held-out error 9.32193472449e-11, 21/24 projected errors still increase, 0 physical values and 1 open derivative theorem
```
