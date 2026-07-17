# Jensen-Window PF Newman Polymath-15 Endpoint C0 Recurrence Transition Gate

Date: 2026-07-17

Status: exact adjacent-cutoff recurrence and corrected-jump diagnostics.
This is not a proof of `Lambda <= 0` or RH; the complex-collar
transition theorem remains open.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.py
```

## Exact Recurrence

Direct simplification of the explicit C0 formula gives

```text
C0(p+2)+C0(p)=exp(pi*i*(p^2/2+p+3/8))
p_N=1-2*(a-N); replacing N by N+1 sends p_N to p_N+2
(-1)^(N+1)C0(p+2)-(-1)^N C0(p)=(-1)^(N+1)*exp(pi*i*(p^2/2+p+3/8))
```

Therefore

```text
J_(N+1)-J_N equals one added two-saddle Dirichlet block minus the explicit signed endpoint-recurrence block
```

The large adjacent corrected sums never need to be subtracted.

## Boundary Scout

At exact cutoff boundaries, independent precision runs give

| n | c=tL | L | Delta J | Delta J'/L | scaled value | scaled first |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 0 | 4.60517 | 0.00031271044 | -2.5312969e-5 | 0.09888772 | 0.008004664 |
| 10 | 0.1 | 4.605157 | 0.00030386881 | -2.4594912e-5 | 0.098896 | 0.008004568 |
| 10 | 1 | 4.605034 | 0.00023520037 | -1.8982991e-5 | 0.09916551 | 0.008003635 |
| 100 | 0 | 9.21034 | 9.8971341e-7 | -3.991278e-8 | 0.09897134 | 0.003991278 |
| 100 | 0.1 | 9.21034 | 9.343553e-7 | -3.7680099e-8 | 0.09897192 | 0.003991278 |
| 100 | 1 | 9.21034 | 5.5687021e-7 | -2.2444639e-8 | 0.09902699 | 0.00399128 |
| 100 | 4 | 9.210338 | 9.9858539e-8 | -3.9913034e-9 | 0.09985813 | 0.003991287 |
| 1000 | 0 | 13.81551 | 3.129775e-9 | -8.4141194e-11 | 0.09897218 | 0.002660778 |
| 1000 | 0.1 | 13.81551 | 2.8708745e-9 | -7.7180692e-11 | 0.09897242 | 0.002660778 |
| 1000 | 1 | 13.81551 | 1.3201433e-9 | -3.5482048e-11 | 0.09899678 | 0.002660778 |
| 1000 | 4 | 13.81551 | 9.9365859e-11 | -2.6607784e-12 | 0.09936586 | 0.002660778 |
| 1000000 | 0 | 27.63102 | 9.8972185e-17 | -1.3303887e-18 | 0.09897218 | 0.001330389 |
| 1000000 | 0.1 | 27.63102 | 8.3274767e-17 | -1.1193826e-18 | 0.09897225 | 0.001330389 |
| 1000000 | 1 | 27.63102 | 1.7601114e-17 | -2.3658029e-19 | 0.09897834 | 0.001330389 |
| 1000000 | 4 | 27.63102 | 9.9070603e-20 | -1.3303887e-21 | 0.0990706 | 0.001330389 |
| 100000000000 | 0 | 50.65687 | 3.1297753e-29 | -2.2947592e-31 | 0.09897218 | 0.0007256666 |
| 100000000000 | 0.1 | 50.65687 | 2.2804116e-29 | -1.6720032e-31 | 0.0989722 | 0.0007256666 |
| 100000000000 | 1 | 50.65687 | 1.3198397e-30 | -9.6769192e-33 | 0.09897401 | 0.0007256666 |
| 100000000000 | 4 | 50.65687 | 9.9001466e-35 | -7.2566657e-37 | 0.09900147 | 0.0007256666 |
| 100000000000 | 25 | 50.65687 | 1.3350682e-63 | -9.6769192e-66 | 0.100116 | 0.0007256666 |

Here the final columns multiply by `exp(5L/4+cL/16)`. Their maxima
are `0.100115992091320294618488` for the value and
`0.008004663775951158518726554` for the first jet. This strongly supports
the next-saddle scale but remains finite point evidence.

## Live Estimate

The cutoff problem is now the explicit statement

```text
Prove |Delta J|+|Delta J'|/L<=exp(-5L/4) on every L>=50 radius-1/L cutoff collar with 0<tL<=25
```

If established,

```text
The live bound is absorbed by the existing 2500*exp(-3L/4) fixed-cell C1 collar budget
```

The remaining difficulty is a uniform complex-ratio estimate for one
block. Neither the recurrence nor the point scout proves that estimate,
corrected first-jet transversality, `Lambda <= 0`, or RH.
