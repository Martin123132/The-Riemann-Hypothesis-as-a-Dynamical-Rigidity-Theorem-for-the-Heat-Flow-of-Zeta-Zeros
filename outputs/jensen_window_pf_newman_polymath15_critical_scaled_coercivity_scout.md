# Jensen-Window PF Newman Polymath-15 Critical Scaled Coercivity Scout

Date: 2026-07-17

Status: high-precision corrected-main diagnostics. This is not a proof
of a uniform sign, `Lambda <= 0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.py
```

The table evaluates the uncorrected finite curvature, the explicit
Polymath-15 endpoint correction, and the corrected curvature for six
scaled times at four frequencies:

| x | c=tL | N | p | C[P] | Q=Ct/A | C[P-Q] | exact xi C (c=0) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000 | 0 | 8 | -0.841241 | 6.558404 | -0.251397 | 6.505134 | 6.5053979 |
| 1000 | 0.5 | 8 | -0.842042 | 5.9802861 | -0.219485 | 6.0052784 | - |
| 1000 | 1 | 8 | -0.842842 | 5.5122024 | -0.191624 | 5.5932732 | - |
| 1000 | 2 | 8 | -0.844442 | 4.8414656 | -0.146063 | 4.9852987 | - |
| 1000 | 3 | 8 | -0.846043 | 4.4373256 | -0.111335 | 4.5999261 | - |
| 1000 | 4 | 8 | -0.847643 | 4.2178639 | -0.0848634 | 4.3740351 | - |
| 10000 | 0 | 28 | 0.581042 | 3.5981149 | -0.102571 | 3.792411 | 3.7925524 |
| 10000 | 0.5 | 28 | 0.580876 | 3.94015 | -0.0832328 | 4.0705623 | - |
| 10000 | 1 | 28 | 0.58071 | 4.2953453 | -0.0675405 | 4.3766137 | - |
| 10000 | 2 | 28 | 0.580378 | 5.0175063 | -0.0444738 | 5.0348579 | - |
| 10000 | 3 | 28 | 0.580047 | 5.7243468 | -0.0292849 | 5.7097357 | - |
| 10000 | 4 | 28 | 0.579715 | 6.3932938 | -0.0192834 | 6.3659651 | - |
| 100000 | 0 | 89 | 0.587588 | 33.225674 | 0.0580922 | 33.885483 | 33.885666 |
| 100000 | 0.5 | 89 | 0.587549 | 31.430371 | 0.0438732 | 31.929807 | - |
| 100000 | 1 | 89 | 0.58751 | 30.035818 | 0.0331345 | 30.41297 | - |
| 100000 | 2 | 89 | 0.587432 | 28.086912 | 0.0188992 | 28.299491 | - |
| 100000 | 3 | 89 | 0.587354 | 26.857202 | 0.0107797 | 26.974455 | - |
| 100000 | 4 | 89 | 0.587276 | 26.040851 | 0.00614849 | 26.103736 | - |
| 1000000 | 0 | 282 | 0.810416 | 262.75474 | -0.0428923 | 262.42238 | 262.42234 |
| 1000000 | 0.5 | 282 | 0.810407 | 206.71341 | -0.0301456 | 206.49439 | - |
| 1000000 | 1 | 282 | 0.810397 | 168.32894 | -0.021187 | 168.17816 | - |
| 1000000 | 2 | 282 | 0.810377 | 121.39495 | -0.0104655 | 121.31572 | - |
| 1000000 | 3 | 282 | 0.810358 | 95.01822 | -0.00516949 | 94.9732 | - |
| 1000000 | 4 | 282 | 0.810338 | 78.577155 | -0.00255351 | 78.551014 | - |

All 24 corrected point values are positive. Coarse/fine precision
reruns have maximum relative corrected-curvature delta
`0.0`.

At `c=0`, the corrected finite expression is independently compared
with `H_0(x)=xi((1+i*x)/2)/8`, including derivatives of the real
normalization. The largest relative curvature discrepancy is
`4.05691788e-5` and the smallest is
`1.383212199e-7`.

The correction is not cosmetic: it changes the curvature visibly and
systematically improves agreement with exact xi. The sampled margins
are substantial, but point positivity is not an interval or asymptotic
theorem. The live obligation remains a uniform corrected coercivity
bound throughout the critical scaled layer.
