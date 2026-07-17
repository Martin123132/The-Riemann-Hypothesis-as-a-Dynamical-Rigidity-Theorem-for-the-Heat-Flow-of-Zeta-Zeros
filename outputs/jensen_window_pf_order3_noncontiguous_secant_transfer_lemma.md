# Jensen-Window PF Order-Three Noncontiguous Secant-Transfer Lemma

Date: 2026-07-12

Status: exact all-column reshaped-Hankel order-two and order-three
theorem at lambda=0 with order four open. This is not a proof of
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_order3_noncontiguous_secant_transfer_lemma`.

```text
work/rh_compute/results/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.json
python work/rh_compute/scripts/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.py
```

Current result:

```text
validated Jensen-window PF order-three noncontiguous secant-transfer lemma: 8 rows, 0 issues, 2 exact identities, 1 secant-averaging lemma, 1 arbitrary-column order-two theorem, 1 arbitrary-column order-three theorem, 1 open order-four handoff
```

## Planar Normalization

Fix `n` and divide the column with offset `j` by its positive first
entry. The three-row Hankel column becomes

```text
A_(n+j)*(1, r_(n+j), r_(n+j)*r_(n+j+1))^T
P_j=(u_j,v_j)=(r_(n+j),r_(n+j)*r_(n+j+1)).
```

Strict moment log-concavity gives `x_k<1`, so

```text
u_(j+1)<u_j.
```

## Local Orientation

Define the edge slope

```text
sigma_j=(v_(j+1)-v_j)/(u_(j+1)-u_j)
```

Direct determinant algebra gives

```text
det[(1,u_j,v_j)^T]_(j=l,l+1,l+2)=(u_(l+1)-u_l)*(u_(l+2)-u_(l+1))*(sigma_(l+1)-sigma_l)
```

The two abscissa gaps have the same negative sign. Therefore the proved
contiguous theorem `D_(3,n+j)<0` for every shift is exactly

```text
sigma_(j+1)<sigma_j for every j>=0.
```

## Arbitrary Columns

For `a<b`, the longer secant slope

```text
S_(a,b)=(v_b-v_a)/(u_b-u_a)
```

is the weighted average of `sigma_a,...,sigma_(b-1)` with positive
weights `u_j-u_(j+1)`. If `a<b<c`, every edge in the first interval
has larger slope than every edge in the second. Hence

```text
S_(a,b)>S_(b,c),
det(P_a,P_b,P_c)=(u_b-u_a)*(u_c-u_b)*(S_(b,c)-S_(a,b))<0
```

and restoring the positive column factors proves

```text
R_(3,n)(j_1,j_2,j_3)<0
for every n>=0 and 0<=j_1<j_2<j_3 at lambda=0.
```

For two columns, the normalized determinant is simply the difference
of two decreasing ratios, so

```text
R_(2,n)(j_1,j_2)<0 for every j_1<j_2.
```

Thus the complete arbitrary-column reshaped-Hankel layers of orders two
and three are now closed. Order four has no automatic planar secant
reduction and remains the next compound target.

```text
outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
