# Jensen-Window PF Monotone-Contraction Sparse Degree-7 Subdivision Scout

Date: 2026-07-06

Status: exact sparse degree-7 subdivision certificate. This is not a
proof of Jensen-window PF-infinity, all-shape Schur positivity, RH,
or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout`.

Proof boundary: this artifact repairs the finite degree-7 `m=10`
one-shot global Bernstein certificate obstruction by dyadic subdivision
in `s0`. It does not prove all column rows, all degrees, all Schur
shapes, zeta cone entry, RH, or `Lambda <= 0`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.py
```

Current result:

```text
validated Jensen-window PF monotone-contraction sparse degree-7 subdivision scout: 3 dyadic slabs, 785400 slab Bernstein coefficients, 0 negative slab coefficients, 0 zero slab coefficients, repaired m=10 obstruction, 0 issues
```

## Global Obstruction

The one-shot global Bernstein net for degree 7, `m=10`, has:

```text
shape=[25, 17, 11, 7, 4, 2]
min=-4928
negative coefficients=126
zero coefficients=0
```

## Subdivision Certificate

Splitting only the first monotone parameter `s0` gives:

```text
s0 in [0,1/2]: count=261800, min=297689420329/16384, min_index=[24, 0, 10, 0, 0, 0], negative=0, zero=0
s0 in [1/2,3/4]: count=261800, min=345475759715905/268435456, min_index=[24, 0, 10, 0, 0, 0], negative=0, zero=0
s0 in [3/4,1]: count=261800, min=225767871/129536, min_index=[20, 0, 10, 0, 0, 0], negative=0, zero=0
```

## Consequence

The degree-7 m=10 one-shot global Bernstein obstruction is repaired by splitting only s0 into the three dyadic slabs [0,1/2], [1/2,3/4], and [3/4,1]. Across those slabs, 785,400 exact Bernstein coefficients are strictly positive. Combined with the global degree-7 m<=9 certificates, this gives bounded degree-7 column positivity through m=10 under monotone contractions, while leaving all-m/all-degree and zeta cone-entry theorems open.

Integration:

```text
outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md
outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The degree-7 m=10 one-shot global Bernstein obstruction is repaired by splitting only s0 into the three dyadic slabs [0,1/2], [1/2,3/4], and [3/4,1]. Across those slabs, 785,400 exact Bernstein coefficients are strictly positive. Combined with the global degree-7 m<=9 certificates, this gives bounded degree-7 column positivity through m=10 under monotone contractions, while leaving all-m/all-degree and zeta cone-entry theorems open.
