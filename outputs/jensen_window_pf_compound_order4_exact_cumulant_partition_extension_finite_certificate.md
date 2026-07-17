# Jensen-Window PF Order-Four Partition-Extension Finite Certificate

Date: 2026-07-13

Status: rigorous finite partition extension and high-jet theorem.
This is not a proof of the exact cumulant corridors, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.py
```

## Partition Extension

Exact graded exponential and Gaussian tilted-moment recurrences extend the
formal partition through epsilon fourteen. Centered Arb Taylor enclosures
preserve the cancellations in all 78 scalar coefficient functions and prove
on `2<=u<=20`, uniformly for `|z|<=1`,

```text
||Z_11||_1<2,
||Z_12||_1<2,
||Z_13||_1<21/10,
||Z_14||_1<12/5.
```

## High Jets

The finite covers also prove

```text
L_13<200, L_14<400, L_15<800, L_16<1600,
L_17(v)<4000 for 39/20<=v<=401/20.
```

The partition cover contains `5400` blocks;
the shifted collar cover contains `5430` blocks.
Both caches are deterministic and source-hashed in the JSON artifact.

## Remaining Boundary

The remaining finite step is scalar composition: combine these four small
partition coefficients with the epsilon-fifteen Bell remainder, the added
formal tails, and the exact potential Taylor remainder. No exact cumulant
corridor is promoted by this finite input theorem.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md
outputs/formal_core.md
```
