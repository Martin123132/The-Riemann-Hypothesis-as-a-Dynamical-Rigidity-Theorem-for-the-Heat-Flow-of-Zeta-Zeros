# Jensen-Window PF Cofinal-Degree Polar-Closure Lemma

Date: 2026-07-10

Status: exact cofinal-degree reduction with the terminal theorem open.
This is not a proof of PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_cofinal_degree_polar_closure_lemma.json
python work/rh_compute/scripts/jensen_window_pf_cofinal_degree_polar_closure_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_cofinal_degree_polar_closure_lemma.py
```

Current result:

```text
validated Jensen-window PF cofinal-degree polar-closure lemma: 10 rows, 0 issues, 3 exact polar identities, 1 interlacing theorem, 1 multiplicity theorem, 1 finite-tower closure, 1 cofinal-degree closure, 1050 finite Sturm rows, 2875 contraction-only rows, 1 open terminal handoff
```

## Polar Descent

At a fixed shift `n`, normalized Jensen windows satisfy

```text
P_(d,n)=P_(d+1,n)-w*P_(d+1,n)'/(d+1)
```

If `P_D=product_i(1+alpha_i*w)` with `alpha_i>0`, then

```text
If P_D=product_i(1+alpha_i*w), then P_(D-1)/P_D=(1/D)*sum_i 1/(1+alpha_i*w).
```

The reciprocal sum is strictly decreasing between its poles. Its zeros
are real, negative, and interlace the roots of `P_D`; repeated terminal
roots lose one unit of multiplicity. Therefore

```text
H_(D,n) => H_(d,n) for every 0<=d<=D
```

## Cofinal Closure

For any fixed shift and finite target degree `d`, one may choose a
hyperbolic terminal degree `D>=d` and descend. Consequently

```text
If {D:H_(D,n)} is unbounded for a fixed n, then H_(d,n) holds for every finite d.
```

Thus the all-degree target does not require certifying every degree
directly. An unbounded terminal subsequence at every shift is sufficient.

## Evidence Audit

Current rigorous finite Sturm data contain

```text
degrees=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
rows=1050
failed_or_inconclusive=0
shifts=0..20 on five cached nonnegative heat parameters.
```

The 2,875 contraction rows remain a logically separate diagnostic; the
new root counts come from Sturm certificates, not from promotion of
contraction inequalities.

## Remaining Handoff

The next theorem must construct, for every shift `n`, an unbounded
sequence `D_j` of hyperbolic terminal degrees. The reduction is exact,
but no large-degree asymptotic or compactness theorem currently supplies
those terminals for the zeta coefficients.
