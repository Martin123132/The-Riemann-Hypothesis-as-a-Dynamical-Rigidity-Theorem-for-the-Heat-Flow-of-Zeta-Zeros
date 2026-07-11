# Jensen-Window PF Monotone-Contraction Sparse Degree-7 Frontier Scout

Date: 2026-07-06

Status: exact sparse degree-7 frontier diagnostic. This is not a proof
of Jensen-window PF-infinity, all-shape Schur positivity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout`.

Proof boundary: this artifact proves finite exact Bernstein sign
certificates for bounded degree-7 column rows `m<=9` under monotone
contractions. It also records that the same one-shot global Bernstein
certificate fails at `m=10`. That failure is not a proof that the
polynomial is negative, and it does not prove all column rows or all
Schur shapes.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.py
```

Current result:

```text
validated Jensen-window PF monotone-contraction sparse degree-7 frontier scout: 9 positive rows, 1 certificate-obstruction row, 932691 Bernstein coefficients, first obstruction m=10, 126 negative Bernstein coefficients, 0 zero Bernstein coefficients, 0 issues
```

## Certified Region

```text
degree 7: 0 <= x1 <= x2 <= x3 <= x4 <= x5 <= x6 <= 1
x_i=1-prod_{r<=i}(1-s_r), all s_r in [0,1]
```

## Sparse Recurrence

After removing the positive monomial `A^m*rho^m`, the column recurrence is:

```text
Q_0=1 and Q_m=sum_{j=1..min(7,m)} (-1)^(j-1) binom(7,j) prod_{i=1..j-1} x_i^(j-i) Q_{m-j}, after factoring A^m*rho^m.
```

## Degree-7 Rows

```text
mcs7_d7_m1: m=1, class=strict_positive_global_bernstein_certificate, recurrence terms=1, power terms=1, Bernstein multidegree=[], Bernstein count=1, min=7, min_index=[], negative=0, zero=0
mcs7_d7_m2: m=2, class=strict_positive_global_bernstein_certificate, recurrence terms=2, power terms=2, Bernstein multidegree=[1], Bernstein count=2, min=28, min_index=[1], negative=0, zero=0
mcs7_d7_m3: m=3, class=strict_positive_global_bernstein_certificate, recurrence terms=3, power terms=5, Bernstein multidegree=[3, 1], Bernstein count=8, min=84, min_index=[3, 0], negative=0, zero=0
mcs7_d7_m4: m=4, class=strict_positive_global_bernstein_certificate, recurrence terms=5, power terms=29, Bernstein multidegree=[6, 3, 1], Bernstein count=56, min=210, min_index=[6, 0, 0], negative=0, zero=0
mcs7_d7_m5: m=5, class=strict_positive_global_bernstein_certificate, recurrence terms=7, power terms=333, Bernstein multidegree=[10, 6, 3, 1], Bernstein count=616, min=462, min_index=[10, 0, 0, 0], negative=0, zero=0
mcs7_d7_m6: m=6, class=strict_positive_global_bernstein_certificate, recurrence terms=11, power terms=5916, Bernstein multidegree=[15, 10, 6, 3, 1], Bernstein count=9856, min=924, min_index=[15, 0, 0, 0, 0], negative=0, zero=0
mcs7_d7_m7: m=7, class=strict_positive_global_bernstein_certificate, recurrence terms=15, power terms=142578, Bernstein multidegree=[21, 15, 10, 6, 3, 1], Bernstein count=216832, min=1716, min_index=[21, 0, 0, 0, 0, 0], negative=0, zero=0
mcs7_d7_m8: m=8, class=strict_positive_global_bernstein_certificate, recurrence terms=21, power terms=142650, Bernstein multidegree=[21, 15, 10, 6, 3, 1], Bernstein count=216832, min=3003, min_index=[21, 0, 0, 0, 0, 0], negative=0, zero=0
mcs7_d7_m9: m=9, class=strict_positive_global_bernstein_certificate, recurrence terms=28, power terms=152563, Bernstein multidegree=[22, 15, 10, 6, 3, 1], Bernstein count=226688, min=72835/22, min_index=[21, 0, 10, 0, 0, 0], negative=0, zero=0
mcs7_d7_m10: m=10, class=global_bernstein_certificate_obstruction, recurrence terms=38, power terms=182753, Bernstein multidegree=[24, 16, 10, 6, 3, 1], Bernstein count=261800, min=-4928, min_index=[22, 0, 10, 0, 0, 0], negative=126, zero=0
```

The first obstruction row has representative negative Bernstein coefficients:

```text
index=[22, 0, 9, 0, 0, 0], coefficient=-36267/230
index=[22, 0, 9, 0, 0, 1], coefficient=-35651/230
index=[22, 0, 9, 0, 1, 0], coefficient=-155617/828
index=[22, 0, 9, 0, 1, 1], coefficient=-770693/4140
index=[22, 0, 9, 0, 2, 0], coefficient=-225841/1035
index=[22, 0, 9, 0, 2, 1], coefficient=-9779/45
index=[22, 0, 9, 0, 3, 0], coefficient=-342881/1380
index=[22, 0, 9, 0, 3, 1], coefficient=-342881/1380
index=[22, 0, 9, 1, 1, 0], coefficient=-174251/12420
index=[22, 0, 9, 1, 1, 1], coefficient=-155771/12420
index=[22, 0, 9, 1, 2, 0], coefficient=-974897/24840
index=[22, 0, 9, 1, 2, 1], coefficient=-956417/24840
```

## Consequence

The sparse exact axis-Bernstein transform certifies degree 7 monotone-contraction column rows through m=9, covering 670,891 strictly positive Bernstein coefficients. At degree 7 m=10 the same global Bernstein certificate first fails, with 126 negative Bernstein coefficients and minimum -4928. This marks a certificate frontier, not a positivity counterexample or proof obstruction to a subdivided/stronger certificate.

Integration:

```text
outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md
outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The sparse exact axis-Bernstein transform certifies degree 7 monotone-contraction column rows through m=9, covering 670,891 strictly positive Bernstein coefficients. At degree 7 m=10 the same global Bernstein certificate first fails, with 126 negative Bernstein coefficients and minimum -4928. This marks a certificate frontier, not a positivity counterexample or proof obstruction to a subdivided/stronger certificate.
