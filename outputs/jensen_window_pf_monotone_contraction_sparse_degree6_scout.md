# Jensen-Window PF Monotone-Contraction Sparse Degree-6 Scout

Date: 2026-07-06

Status: exact sparse degree-6 diagnostic. This is not a proof of
Jensen-window PF-infinity, all-shape Schur positivity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_monotone_contraction_sparse_degree6_scout`.

Proof boundary: this artifact proves finite exact Bernstein sign
certificates for bounded degree-6 column rows under monotone
contractions. It does not prove that the zeta coefficients satisfy
monotone contractions, and it does not prove all column rows or all
Schur shapes.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree6_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_sparse_degree6_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree6_scout.py
```

Current result:

```text
validated Jensen-window PF monotone-contraction sparse degree-6 scout: 10 degree-6 rows, 63347 Bernstein coefficients, m<=10, 0 negative Bernstein rows, 0 zero Bernstein rows, 0 issues
```

## Certified Region

```text
degree 6: 0 <= x1 <= x2 <= x3 <= x4 <= x5 <= 1
x_i=1-prod_{r<=i}(1-s_r), all s_r in [0,1]
```

## Sparse Recurrence

After removing the positive monomial `A^m*rho^m`, the column recurrence is:

```text
Q_0=1 and Q_m=sum_{j=1..min(6,m)} (-1)^(j-1) binom(6,j) prod_{i=1..j-1} x_i^(j-i) Q_{m-j}, after factoring A^m*rho^m.
```

## Degree-6 Rows

```text
mcs6_d6_m1: m=1, recurrence terms=1, power terms=1, Bernstein multidegree=[], Bernstein count=1, min=6, negative=0, zero=0
mcs6_d6_m2: m=2, recurrence terms=2, power terms=2, Bernstein multidegree=[1], Bernstein count=2, min=21, negative=0, zero=0
mcs6_d6_m3: m=3, recurrence terms=3, power terms=5, Bernstein multidegree=[3, 1], Bernstein count=8, min=56, negative=0, zero=0
mcs6_d6_m4: m=4, recurrence terms=5, power terms=29, Bernstein multidegree=[6, 3, 1], Bernstein count=56, min=126, negative=0, zero=0
mcs6_d6_m5: m=5, recurrence terms=7, power terms=332, Bernstein multidegree=[10, 6, 3, 1], Bernstein count=616, min=252, negative=0, zero=0
mcs6_d6_m6: m=6, recurrence terms=11, power terms=5916, Bernstein multidegree=[15, 10, 6, 3, 1], Bernstein count=9856, min=462, negative=0, zero=0
mcs6_d6_m7: m=7, recurrence terms=14, power terms=5937, Bernstein multidegree=[15, 10, 6, 3, 1], Bernstein count=9856, min=792, negative=0, zero=0
mcs6_d6_m8: m=8, recurrence terms=20, power terms=6571, Bernstein multidegree=[16, 10, 6, 3, 1], Bernstein count=10472, min=1287, negative=0, zero=0
mcs6_d6_m9: m=9, recurrence terms=26, power terms=8472, Bernstein multidegree=[18, 11, 6, 3, 1], Bernstein count=12768, min=2002, negative=0, zero=0
mcs6_d6_m10: m=10, recurrence terms=35, power terms=13458, Bernstein multidegree=[21, 13, 7, 3, 1], Bernstein count=19712, min=14794/7, negative=0, zero=0
```

## Consequence

A sparse exact common-denominator Bernstein transform certifies degree 6 monotone-contraction column rows through m=10, with 63,347 strict positive Bernstein coefficients and no zero or negative Bernstein rows. This extends the bounded column-family evidence beyond the degree-5 band, but it remains finite theorem-search algebra rather than an all-m or zeta cone-entry theorem.

Integration:

```text
outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

A sparse exact common-denominator Bernstein transform certifies degree 6 monotone-contraction column rows through m=10, with 63,347 strict positive Bernstein coefficients and no zero or negative Bernstein rows. This extends the bounded column-family evidence beyond the degree-5 band, but it remains finite theorem-search algebra rather than an all-m or zeta cone-entry theorem.
