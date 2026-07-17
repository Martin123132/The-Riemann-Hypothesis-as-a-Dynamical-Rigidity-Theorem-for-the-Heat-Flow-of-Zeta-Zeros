# Jensen-Window PF Compound Order-Five Lambda=-100 Prefix Certificate

Date: 2026-07-13

Status: rigorous finite `lambda=-100` contiguous order-five prefix
through `n=316`. This is not a proof of the analytic tail, all-shift
order five, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order5_m100_prefix_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_m100_prefix_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_m100_prefix_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_prefix_certificate.py
```

## Stable Coordinate

Use

```text
x_k=A_(k-1)*A_(k+1)/A_k^2, d_k=1-x_k
G_n=d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)
F_n=G_(n+1)^2-x_(n+3)^3*G_n*G_(n+2)
J_n=d_(n+3)*d_(n+5)*F_(n+1)^2-x_(n+4)^4*d_(n+4)^2*F_n*F_(n+2)
```

Exact symbolic determinant algebra gives

```text
H_(5,n)=W_n*J_n
W_n=A_n^5*rho_n^20*x_(n+1)^15*x_(n+2)^10*x_(n+3)^5/(d_(n+3)*d_(n+4)^2*d_(n+5)*G_(n+2))
```

Every factor in `W_n` is positive inside the already completed lower
cones. Thus `H_(5,n)` has exactly the sign of `J_n`.

## Arb Prefix

Seven source files, merged in recorded precedence order, give 1024-bit
outward-rounded balls for

```text
A_k(-100)>0 for every 0<=k<=324.
```

Direct Arb evaluation of the stable coordinate proves

```text
J_n(-100)>0 for every 0<=n<=316,
relative_margin_n(-100)>0 for every 0<=n<=316,
H_(5,n)(-100)>0 for every 0<=n<=316.
```

The weakest bounds occur at the final row:

```text
minimum J row: n=316
J_316=[1.195506752987560996592413174E-45 +/- 8.91E-73]
minimum relative row: n=316
relative_316=[0.00626902754512557813924681177 +/- 7.83E-30]
```

The relative lower bound is about `0.006269`, leaving a visible strict
buffer rather than a zero-containing endpoint enclosure.

## Remaining Tail

The sole endpoint-entry target is now

```text
J_n(-100)>0 for every n>=317,
equivalently H_(5,n)(-100)>0 for every n>=317.
```

The compact-uniform eventual theorem proves positivity after some
non-effective threshold, but it does not splice this explicit prefix to
all shifts. That analytic tail remains open.

```text
outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
