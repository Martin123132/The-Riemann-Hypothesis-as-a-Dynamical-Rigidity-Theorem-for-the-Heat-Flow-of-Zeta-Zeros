# Jensen-Window PF Newman Polymath-15 Cancellation/Zero-Free Wall Gate

Date: 2026-07-17

Status: exact theorem-search gate separating the cancellation and
zero-free walls. This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.py
```

## Exact Frontier

```text
L=log(x/(4*pi)), c=t*L, N=exp(L/2+o(1)), sigma=Re(s_*)=1/2+c/4+o(1), tau=|Im(s_*)|=2*pi*N^2*(1+o(1))
If M=N^r and sup_I|sum_(n in I)n^(-i*tau)|<<N^(Phi(r)+eta), then the weighted dyadic block is <<N^(Phi(r)-r/2-c*r*(2-r)/8+eta)
Decay at scaled time c requires Phi(r)<r/2+c*r*(2-r)/8 on every untreated radius r
```

The eleven-pair published envelope gives

```text
The current published exponent-pair hull gives c_*=4911678521/1933561194=2.540223984760009...
```

Its maximum occurs at the pair-2/pair-1 transition:

```text
At r_*=125662/155153, Phi_known=220633/310306 but the c=2 frontier is 15549101725/24072453409; the exact phase-exponent excess is 3133668399/48144906818
At r_*, a symmetric pair (theta,1/2+theta) can meet the local c=2 frontier only if theta<=2900341791/24072453409=0.120483846898448..., versus the published theta=13/84
```

The symmetric-pair number is only a necessary local benchmark. A clean
sufficient condition for `c=2` is the complete phase curve
`Phi(r)<=r-r^2/4-delta*r` for some `delta>0`. The radius-proportional
margin becomes an `M^(-delta)` dyadic gain even as `r` tends to zero.

## Frontier Table

| active transition | r | required c | excess above c=2 curve |
|---:|---:|---:|---:|
| 10 -> 9 | 0.461538461538462 | 2.426666666666667 | 0.037869822485207 |
| 9 -> 8 | 0.502948240762826 | 2.459744358354176 | 0.043269963708578 |
| 8 -> 7 | 0.565028002489110 | 2.486656146960047 | 0.049322561645942 |
| 7 -> 6 | 0.584722760759985 | 2.498315468187929 | 0.051547297714817 |
| 6 -> 5 | 0.597504888587748 | 2.503774097755819 | 0.052770190978506 |
| 5 -> 4 | 0.622271377649349 | 2.506330852075120 | 0.054261014614580 |
| 4 -> 3 | 0.755407653910150 | 2.520153596080002 | 0.061129398866559 |
| 3 -> 2 | 0.788496732026144 | 2.535377425991032 | 0.063928506130121 |
| 2 -> 1 | 0.809923108157754 | 2.540223984760009 | 0.065088263870695 |
| 1 -> 0 | 0.857142857142857 | 2.527777777777778 | 0.064625850340136 |
| pair 0 endpoint | 1.000000000000000 | 2.476190476190476 | 0.059523809523810 |

## Two Different Walls

```text
For every fixed c>2, sigma>1 and the Euler product supplies a zeta floor; therefore 2<c<=c_* is a cancellation gap, not a zero-free gap
```

Thus `2<c<=c_*` is a concrete analytic-number-theory target: improve
the weighted log-phase cancellation while retaining the existing Euler
product floor. It is not yet closed.

At the endpoint `c=2`, the leading saddle lies on the zeta 1-line.
The relevant published results are:

- [Mossinghoff-Trudgian-Yang zero-free region](https://arxiv.org/abs/2212.06867)
- [Cully-Hugill-Leong 1-line estimates](https://arxiv.org/abs/2312.09412)
- [Leong reciprocal estimates](https://arxiv.org/abs/2405.04869)

```text
At c=2, sigma=1+o(1). Published zero-free estimates reach sigma>=1-1/(13*log(tau)); in that strip |1/zeta|<=4904*log(tau)^(11/12), while published 1-line logarithmic-derivative bounds are O(log(tau)/loglog(tau))
A cancellation envelope satisfying Phi(r)<=r-r^2/4-delta*r for some delta>0 transfers D_0,D_1,D_2 to zeta at c=2; the resulting M^(-delta) dyadic majorant handles growing blocks, while each fixed block vanishes because t=2/L; the published zero-free/1-line bounds then leave phase speed -L/4+o(L) and a polylogarithmic nonzero amplitude floor
```

The `o(L)` logarithmic phase correction cannot cancel the normalizer's
`-L/4` phase speed, and the polylogarithmic zeta floor is still vastly
larger than the exponentially small Polymath-15 remainder. This closes
the phase-amplitude step only conditionally on the missing strict
`c=2` cancellation frontier.

For every fixed `c<2`, however, the saddle enters a fixed vertical line
inside the critical strip:

```text
For fixed c<2, sigma=1/2+c/4<1. Existing zero-free regions approach 1 and do not supply a uniform zeta floor on that line; improved exponential-sum cancellation alone does not close the phase-critical-value target
```

This is a route boundary, not an impossibility theorem. A direct
Wronskian, zero-dynamical, or other Xi-specific argument could still
cross it; ordinary exponent-pair improvement by itself cannot provide
the missing nonzero phase amplitude.

## Proof Handoff

```text
The outer improvement would shrink the live thin-collar Wronskian problem from 0<c<=c_*+o(1) to fixed c<2 (plus a vanishing transition strip); it would not prove the inner phase theorem
```

So the next analytic-number-theory milestone is unambiguous: lower the
known cancellation frontier from `c_*` toward `2`. The genuinely inner
problem remains quantitative phase-critical-value avoidance for the
corrected Riemann-Siegel main.
