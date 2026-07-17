# Arb Xi Lambda-Zero Order-Four Prefix Certificate

Date: 2026-07-13

Status: rigorous direct-Xi-series finite certificate through `n=500`.
This is not a proof of an all-shift order-four theorem, PF-infinity, RH, or
`Lambda <= 0`.

```text
work/rh_compute/results/arb_xi_lambda0_order4_prefix_certificate.json
work/rh_compute/results/arb_xi_lambda0_order4_prefix_coefficients_n0_n506_bits24576.jsonl
python work/rh_compute/scripts/arb_xi_lambda0_order4_prefix_certificate.py
python work/rh_compute/scripts/check_arb_xi_lambda0_order4_prefix_certificate.py
```

## Direct Series

Python-flint/Arb expands the exact identity

```text
xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2
A_k(0)=k!*[z^(2k)]xi(1/2+z)/4^(k+1).
```

The expansion uses `24576` bits and series order
`1016`. Each serialized decimal Arb ball is rounded
outward and checked to contain the original high-precision ball.

## Certified Prefix

Recomputing the raw `4x4` Hankel determinant and the stable gap
factorization from the cached coefficient balls proves

```text
A_k(0)>0 for every 0<=k<=506
H_(4,n)(0)>0 for every 0<=n<=500
F_n(0)>0 for every 0<=n<=500
```

The smallest recorded stable margin occurs at `n=500`:

```text
F_500=[5.481268356518127684426961803646365408360017160550209264948367373383654753206792239922810852195036203e-20 +/- 2.95e-120]
F_500/G_501^2=[0.003716298961440233294098765075005020047401588139987551949881424960342387775139012190166540704409788816 +/- 1.73e-103]
```

This extends the previous rigorous endpoint prefix from `n=58` to
`n=500`. It remains finite evidence; it composes with the separate
eventual Xi-asymptotic theorem only after that theorem's threshold is
made explicit.

```text
outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md
outputs/formal_core.md
```
