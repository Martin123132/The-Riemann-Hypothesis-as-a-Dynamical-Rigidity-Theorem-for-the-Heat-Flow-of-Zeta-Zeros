# Newman Theta Modular-Tail Derivative-Budget Scout

Date: 2026-07-24

Status: finite floating-point stress scout. This is not an
interval certificate and not a proof of `Lambda<=0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_modular_tail_derivative_budget_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_modular_tail_derivative_budget_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_tail_derivative_budget_scout.py
```

## Derivative Envelopes

At `T=1/5`, the exact positive-integrand formulas were evaluated
for `m=5,...,9`. Selected rows are:

| N | m | d0 envelope | d1 envelope |
|---:|---:|---:|---:|
| 4 | 5 | 4.41702731795306139e-09 | 1.31456996519174594e-09 |
| 4 | 9 | 4.29203433268627310e-01 | 1.38343398792887401e-01 |
| 5 | 5 | 3.95082652456632926e-14 | 1.26984573011130012e-14 |
| 5 | 9 | 6.25725171603453266e-06 | 2.14442785771739149e-06 |
| 6 | 5 | 1.50345943633798747e-19 | 5.14151045332616098e-20 |
| 6 | 9 | 3.57918340483601246e-11 | 1.29239256840862763e-11 |
| 7 | 5 | 4.05137147171933968e-25 | 1.49656915328871326e-25 |
| 7 | 9 | 1.73288832538941739e-16 | 6.90660969264984362e-17 |
| 8 | 5 | 3.50391277708966872e-30 | 1.34874007778487959e-30 |
| 8 | 9 | 5.69297995828021036e-21 | 2.31531622477655840e-21 |
| 9 | 5 | 1.41934277240627921e-35 | 5.45350082107192712e-36 |
| 9 | 9 | 5.98657049275710999e-26 | 2.39675830292615292e-26 |
| 10 | 5 | 1.48756812227439006e-41 | 5.74558379954753645e-42 |
| 10 | 9 | 1.37044486640016513e-31 | 5.45957185481750024e-32 |

The node ladder, arithmetic-cap comparison, and final-panel mass
are recorded in the JSON artifact. They are convergence checks,
not directed-rounding enclosures.

## Direct-C1 Stress

For each sampled row the scout composes

```text
max(|S|/epsilon0,|4S+xS'|/(4epsilon0+xepsilon1))>1
N=ceil(max(n_*5,n_*9))+kappa.
```

At `x=200`, `kappa=2` first passes at orders `[9, 8]` for the two
sampled times. Immediately below the next arithmetic-count jump,
`x=300`, every `kappa=2,3,4` row fails all sampled orders
`m=5,...,9`, whereas `kappa=5` passes at both sampled times.

This is a finite non-promotion guard. It does not prove that every
fixed collar fails, and it does not certify the passing rows.

## Revised Handoff

The interval programme should test

```text
N=max(N_sad,ceil(K*(1+x)^(3/4)))
```

with rigorous derivative integrals and retained first-jet boxes.
No strict Laguerre theorem, `Lambda<=0`, RH, or Clay-prize result
is claimed.
