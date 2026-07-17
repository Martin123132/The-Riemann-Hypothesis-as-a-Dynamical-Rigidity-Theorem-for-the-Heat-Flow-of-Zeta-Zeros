# Jensen-Window PF Order-Four Formal Cumulant Second-Next Certificate

Date: 2026-07-13

Status: exact formal epsilon-ten and second-next parity theorem.
This is not a proof of the exact cumulant corridors, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.py
```

## Exact Extension

The tilted-Gaussian recurrence now includes normalized potential jets
`L_3,...,L_12` through `epsilon^10`. All 56 coefficients at epsilon powers
one through eight agree exactly with the epsilon-eight expansion.

After corridor scaling, the new layer is

```text
orders 2,3,4: scaled(kappa_r^[10]-kappa_r^[8])=q^-4*D_r
orders 5,6:   scaled(kappa_r^[10]-kappa_r^[8])=q^-3*D_r
orders 7,8:   scaled(kappa_r^[10]-kappa_r^[8])=q^-2*D_r
```

The seven exact rational polynomials contain 30 monomials for odd orders
and 42 for even orders. Their complete formulas are stored in the JSON.

## Proof Boundary

This is a finite subtraction layer for the central theorem, not a claim that
the asymptotic expansion converges. Global coefficient bounds and the actual
beyond-epsilon-ten central and two-tail estimates remain open.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md
outputs/formal_core.md
```
