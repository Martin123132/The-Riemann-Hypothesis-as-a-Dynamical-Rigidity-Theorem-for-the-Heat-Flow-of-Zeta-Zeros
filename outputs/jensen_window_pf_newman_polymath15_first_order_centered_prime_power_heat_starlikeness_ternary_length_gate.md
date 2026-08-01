# Ternary All-Length Prime-Power Heat-Starlikeness Gate

Date: 2026-07-30

Status: exact all-length complete-chain theorem for p=3,M>=5.
This does not promote RH or Lambda<=0.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_ternary_length_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_ternary_length_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_ternary_length_gate.py
```

## Outcome

```text
p=3,M>=5: J_0>17/1000
p=3,M>=5: J_ray>2/125
```

This closes the complete-chain ternary family from length five
onward. It does not close short chains or any global join.

## Unique Defect

With `D_d=A_d-2A_(d+1)+A_(d+2)`, every common group for
`d>=1` is positive because `1/sqrt(3)<3/5`. Pairing the last
common group with the boundary pair gives

```text
C(B,d)=B-2(B+1)aq^3+(B+2)a^2q^4+a^2q^(2d+6)[B+2-2(B+3)aq+(B+4)a^2q^(2d+2)]
(5-8/sqrt(3))/3>0
3C(d+2,d)>=5d+14-(8d+28)/sqrt(3)>0, d>=6
```

The exact small-`d` certificates are:

| d | claimed lower | shifted Bernstein minimum | degree |
|---:|---:|---:|---:|
| 1 | 17/100 | 4.81072000000000000000e-3 | 12 x 4 |
| 2 | 17/100 | 2.40551680000000000000e-1 | 16 x 4 |
| 3 | 17/100 | 4.76292640000000000000e-1 | 20 x 4 |
| 4 | 17/100 | 7.12033600000000000000e-1 | 24 x 4 |
| 5 | 17/100 | 9.41192090285714285714e-1 | 28 x 4 |

Therefore only `D_(M-2)` may be negative.

## Positive Reserve

```text
D_0>11/250 for M>=6; D_0>1/25 for M=5
D_0+(M-1)^2D_(M-2)>1/25, 5<=M<=16
D_(M-2)<0 implies q^(2M-2)<R_M=3[2(M+1)/sqrt(3)-M]/(M+2)
M>=17: (M-1)^2(-D_(M-2))<1/100
```

The exact finite absorption certificates are:

| M | claimed lower | shifted Bernstein minimum | degree |
|---:|---:|---:|---:|
| 5 | 1/25 | 3.72195894238048133793e-1 | 32 x 8 |
| 6 | 1/25 | 4.33485644883324633254e-1 | 50 x 10 |
| 7 | 1/25 | 4.68794473579823615651e-1 | 72 x 12 |
| 8 | 1/25 | 4.80463288164528754686e-1 | 98 x 14 |
| 9 | 1/25 | 4.72081175253224778246e-1 | 128 x 16 |
| 10 | 1/25 | 4.50020115991936912761e-1 | 162 x 18 |
| 11 | 1/25 | 4.18980196893853341437e-1 | 200 x 20 |
| 12 | 1/25 | 3.64659090005845094104e-1 | 242 x 22 |
| 13 | 1/25 | 3.13138611738198154447e-1 | 288 x 24 |
| 14 | 1/25 | 2.69384681205723419588e-1 | 338 x 26 |
| 15 | 1/25 | 2.30081992089100598089e-1 | 392 x 28 |
| 16 | 1/25 | 1.96568186141325412741e-1 | 450 x 30 |

Since `0<=F_(M-2)<=M-1`, the finite certificates give
`J_0>1/50`; the analytic tail gives `J_0>17/1000`.

## Actual Fixed-Ray Transfer

```text
|R_k-1|<r_x=18000/x
sum c_k<5/2; sum k c_k<15/4; absolute current mass<16
|J_ang-J_0|<48r_x=864000/x
|J_ray-J_0|<864020/x+633450/x^2<1/1000
p=3,M>=5: J_ray>2/125
```

No new `pi` is inserted here. Any `pi` in `x=4*pi*exp(L)`
is inherited from the completed-zeta saddle normalization;
the Fejer kernels use the ordinary `2*pi` angular period.

## Boundary

The exact Fourier grouping, nonterminal cluster theorem, finite rational Bernstein absorption, analytic terminal tail, and length-uniform actual transfer prove J_0>17/1000 and J_ray>2/125 for every p=3,M>=5 complete chain. Open: p=2,M>=10, p=2 lengths 2..7, p=3 lengths 2..3, singleton chains, p-free bases, recurrent endpoints, cutoff equality, adjacent charts, the Xi Abel gap, joined winding, contact exclusion, Lambda<=0, PF-infinity, RH, and a prize-level proof.
