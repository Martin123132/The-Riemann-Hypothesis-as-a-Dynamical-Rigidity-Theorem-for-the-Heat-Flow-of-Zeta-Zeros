# Jensen-Window PF Order-Four Exact Cumulant Central-Residual Certificate

Date: 2026-07-13

Status: global exact central residual and unit-disk partition theorem.
This is not a proof of order-four entry, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.py
```

## Central Decomposition

The central exact-minus-epsilon-ten residual is decomposed before any
cumulant differentiation into

```text
the full Gaussian partition coefficients Z_11,...,Z_14,
their two added formal Gaussian tails,
the epsilon-fifteen Bell-polynomial exponential remainder,
the exact seventeenth-order potential Taylor remainder.
```

The finite interval uses the cancellation-preserving coefficient caps
`2, 2, 21/10, 12/5`. The ray uses explicit, deliberately crude high-jet
majorants; exponential q growth leaves ample room.

## Closed Budgets

The exact scalar compositions prove

```text
central residual <1/(500000*q^3), 2<=u<=20
central residual <1/(300000*u*q^3), u>=20
```

The finite central budget ratio is below `3.63963906179695801080133493347414003233043040008906489983026E-1`.
The ray central budget ratio at its worst endpoint is below `3.28701745816252305106076317148152296713745359650908799785395E-63`.

Together with the two formal and two exact tail theorems, this gives

```text
|A_u-P_u^[10]|<1/(100000*q^3),       2<=u<=20,
|A_u-P_u^[10]|<1/(20000*u*q^3),      u>=20.
```

The complex-disk contract therefore yields the simultaneous scaled
exact-minus-epsilon-ten cumulant budgets `9/1000` and `1/(100u)` for
orders two through eight.

## Remaining Boundary

The exact-density remainder is no longer open. The next job is to compose
these exact cumulant corridors with the order-four determinant algebra and
then audit the first-summand/full-kernel handoff. No claim of PF-infinity,
RH, or `Lambda <= 0` is made here.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/formal_core.md
```
