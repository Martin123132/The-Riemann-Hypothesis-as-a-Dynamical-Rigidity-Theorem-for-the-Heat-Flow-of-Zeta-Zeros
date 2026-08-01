# Order-Twelve Lambda=-100 Endpoint Completion

Date: 2026-07-22

Status: rigorous finite endpoint handoff through `n=1492`. This is not
all-shift endpoint positivity or PF-infinity, and it is not a proof of RH or `Lambda<=0`.

```text
python work/rh_compute/scripts/jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.py
```

## Exact Reduction

```text
Q_(12,n)*Q_(10,n+2)=Q_(11,n+1)^2-Q_(11,n)*Q_(11,n+2)
R_n=Q_(11,n+1)^2/(Q_(11,n)*Q_(11,n+2))-1
```

Four retained-integral chunks add exactly 252 coefficients, extending
coverage from `A_1262` through `A_1514`. The enlarged stable chain
overlaps every one of the 1,243 inherited order-eleven balls.

## Finite Theorem

```text
Q_(12,n)(-100)>0 for every 1241<=n<=1492
Q_(12,n)(-100)<0 for n=0,1,2,3 and Q_(12,n)(-100)>0 for every 4<=n<=1492
minimum relative Q11 margin at n=1492: [0.0047350624466398245050529275417750663551151537009613069090810308496548089062461104 +/- 3.79E-83]
```

All 252 new endpoint rows have positive raw and factored numerators,
positive denominators, positive Q12 balls, and zero inconclusive rows.
The first four endpoint rows remain rigorously negative. The sole
remaining endpoint obligation is the analytic tail `n>=1493`, supplied
conditionally by the separate ninth-coordinate curvature target.
