# Dyadic All-Length Prime-Power Heat-Starlikeness Gate

Date: 2026-07-30

Status: exact all-length complete-chain theorem for p=2,M>=10.
This does not promote RH or Lambda<=0.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_length_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_length_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_length_gate.py
```

## Outcome

```text
p=2,M>=10: J_0>3/1000
p=2,M>=10: J_ray>1/500
with preceding M=8 and M=9 gates: p=2,M>=8
```

## Two-Defect Absorption

With `K_d=|1+z+...+z^d|^2`,

```text
2J_0=sum_(d=0)^(M-1)D_dK_d
P|Y+u|^2-N|Y|^2=(P-N)|Y+P*u/(P-N)|^2-PN/(P-N)
L_M=(M+1)A^(alpha+1)/(M+2)^alpha*t^alpha(1-t)/[M-(M+1)r+At]
```

The last two kernels are kept together. Completing the square
replaces the old worst-case quadratic Fejer loss by one scalar
`t^alpha(1-t)` loss.

## Finite Rows

| M | target | shifted Bernstein minimum | degree |
|---:|---:|---:|---:|
| 10 | 1/10 | 4.179359192621830930228904698737409331408e-2 | 72 x 18 x 7 |
| 11 | 1/10 | 6.853046984702472631609457095841557896616e-2 | 90 x 20 x 7 |
| 12 | 1/10 | 8.441219989466944779603490648525077223886e-2 | 110 x 22 x 7 |
| 13 | 1/10 | 9.360971201758083204588845994191218820455e-2 | 132 x 24 x 7 |
| 14 | 1/10 | 9.883679240631589144419273283113298925423e-2 | 156 x 26 x 7 |

```text
10<=M<=14: D_0+sum_(d=1)^7D_dK_d>1/10
10<=M<=14: terminal pair>-1/25
10<=M<=14: J_0>3/100
```

## All-Length Tail

For `M>=15`, retain seven initial groups in each of
`D_0,...,D_4`. With `R=r_12`,

```text
R=r_12; c_k=R^kq^[k(25-k)]
sum_(d=0)^4K_d sum_(k=0)^7T_(k,d)>1/100
M>=15: terminal pair>-1/250
M>=15: J_0>3/1000
```

The length-free core certificate has:

| target | shifted Bernstein minimum | degree |
|---:|---:|---:|
| 1/100 | 1.198385086563061904882218746952313374068e-3 | 141 x 20 x 4 |

## Actual Fixed-Ray Transfer

```text
|R_k-1|<r_x=18000/x
sum c_k<4; sum k c_k<12; absolute current mass<64
|J_ray-J_0|<3456096/x+1621632/x^2<1/1000
p=2,M>=10: J_ray>1/500
```

No new `pi` is inserted here. The Fejer numerators use the
ordinary `2*pi` period of `z=exp(i theta)`; any `pi` in the
saddle scale is inherited from completed-zeta normalization.

## Boundary

The exact seven-kernel finite certificates, all-length four-kernel initial-group certificate, completed-square terminal estimate, and length-uniform actual transfer prove J_0>3/1000 and J_ray>1/500 for every p=2,M>=10 complete chain. Together with preceding gates, complete dyadic chains are covered for M>=8. Open: p=2 lengths 2..7, p=3 lengths 2..3, singleton chains, p-free bases, recurrent endpoints, cutoff equality, adjacent charts, the Xi Abel gap, joined winding, contact exclusion, Lambda<=0, PF-infinity, RH, and a prize-level proof.
