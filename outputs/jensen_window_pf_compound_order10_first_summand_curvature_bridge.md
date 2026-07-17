# Jensen-Window PF Order-Ten First/Full Curvature Bridge

Date: 2026-07-16

Status: exact full-kernel transfer with one open continuous
first-summand theorem. This is not a proof of PF-infinity, RH, or
`Lambda <= 0`, and it is not an order-ten entry theorem.

```text
work/rh_compute/results/jensen_window_pf_compound_order10_first_summand_curvature_bridge.json
python work/rh_compute/scripts/jensen_window_pf_compound_order10_first_summand_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_first_summand_curvature_bridge.py
```

## Stable Floor

```text
W_j^(1)>=8/(2*j+1)-4201/j^2>1/(5*j), j>=1251
W_j>=1004/(125*(2*j+1))-4751/j^2>1/(5*j), j>=1251
min(W_j,W_j^(1))>1/(5*j), j>=1251
phi(min(W_j,W_j^(1)))<=1/min(W_j,W_j^(1))<5*j
```

Both floor inequalities are certified by shifted polynomials with
strictly positive coefficients.

## Rational Error Envelope

Starting from the proved inverse-twelfth-power moment defect, the exact
positive-stencil lemma propagates `C/j^p` bounds through the seven
stable logarithms. The final scaled transfer envelope is

```text
k^2*|Z_k-Z_k^(1)| <= 8.9855798484483117128975379053613939076528584196934E+0
|Z_k-Z_k^(1)|<10/k^2 for every integer k>=1252
```

The power-ten estimate used at order nine is not asserted to be
impossible here; the power-twelve input is the effective estimate that
makes this particular rigorous envelope close cleanly.

## Remaining Premise

```text
z_1''(t)<=4200/t^2 for every real t>=1251
z_1''(t)<=4200/t^2 => Z_k^(1)<=4200*[-log(1-1/k^2)]<4201/k^2, k>=1252
z_1''(t)<=4200/t^2 on t>=1251 => Z_k<4201/k^2+10/k^2=4211/k^2<5500/k^2, k>=1252
z_1''(t)<=4200/t^2 on t>=1251 => Q_(10,n)(-100)>0 for every n>=1243
```

Together with the finite block, that premise would prove positivity only
for `n>=4`; the four negative endpoint rows remain negative.
