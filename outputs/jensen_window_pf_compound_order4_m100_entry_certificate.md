# Jensen-Window PF Compound Order-Four Lambda=-100 Entry Certificate

Date: 2026-07-13

Status: all-shift contiguous order-four entry theorem at `lambda=-100`.
This proves one complete compound layer at one heat parameter. It is not
a proof of forward order-four invariance, arbitrary-column order four,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_m100_entry_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_m100_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_m100_entry_certificate.py
```

## Sign Coordinate

Desnanot-Jacobi condensation and the completed order-two and order-three
signs give

```text
H_(4,n)>0 iff T_(n+1)^2>T_n*T_(n+2).
```

The exact order-three gap factorization then gives the equivalent stable
penalty condition

```text
P_n=log(G_n*G_(n+2)/G_(n+1)^2)
H_(4,n)>0 iff P_n<-3*log(x_(n+3)).
```

## Repaired Prefix

The 1024-bit repaired Arb source proves

```text
H_(4,n)(-100)>0, 0<=n<=316.
```

All `317` margins are strict. The minimum lower enclosure is
`4.71789538689370302040913771438131967598144836682129515711914E-3` at `n=None`.

## Global Curvature

The compact interval certificate covers `319<=t<=V'(2)`. The complete
finite exact-corridor cover and the analytic `u>=20` ray prove the
localized ceiling on every mode `u>=2`; the stable localization converts
that ceiling to the first-summand curvature. Therefore

```text
K_1(t)<=7/(2*t^2) for every real t>=319.
```

## Analytic Tail

For `n>=317`, put `k=n+3>=320`. Then `k+s>=319` throughout the
tent interval, and

```text
P_n^(1)=integral_[-1,1](1-|s|)*K_1(k+s) ds
       <=(7/2)*[-log(1-1/k^2)]
       <=7/(2*(k^2-1))
       <=18/(5*k^2).
```

The already-proved complete-kernel perturbation gives

```text
|P_n-P_n^(1)|<=2/(5*k^2),
P_n<=4/k^2.
```

The scaled-defect anchor supplies the strict remaining buffer:

```text
-3*log(x_k)>3*d_k>=753/(250*(2*k+1))>4/k^2.
```

After `k=320+m`, the numerator of the last strict comparison has
positive coefficients

```text
[753, 479920, 76466200]
```

and endpoint value `76466200`. Hence
`P_n<-3*log(x_k)` and therefore `H_(4,n)>0` for every `n>=317`.

## Entry Theorem

Combining the finite prefix and analytic tail proves

```text
H_(4,n)(-100)>0 for every integer n>=0.
```

This closes the contiguous order-four entry problem at `lambda=-100`.
The next live theorem is forward propagation: derive the order-four
compound flow and prevent a first loss of positivity from escaping to
infinite index before carrying the layer through `lambda=0`.

```text
outputs/jensen_window_pf_compound_order4_condensation_gate.md
outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/formal_core.md
```
