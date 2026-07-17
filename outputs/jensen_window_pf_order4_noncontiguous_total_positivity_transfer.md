# Jensen-Window PF Order-Four Noncontiguous Total-Positivity Transfer

Date: 2026-07-13

Status: exact arbitrary-column reshaped-Hankel order-four theorem at
lambda zero. This is not a proof of contiguous order five, PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_order4_noncontiguous_total_positivity_transfer`.

```text
work/rh_compute/results/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json
python work/rh_compute/scripts/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.py
python work/rh_compute/scripts/check_jensen_window_pf_order4_noncontiguous_total_positivity_transfer.py
```

## Published Criterion

Gasca and Pena, `Total positivity and Neville elimination`, Linear
Algebra and its Applications 165 (1992), 25-44, prove the rectangular
initial-minor criterion for strict total positivity:

```text
an m by p real matrix is strictly totally positive if and only if
all its initial minors are strictly positive.
```

An initial minor uses consecutive rows and consecutive columns and
touches the first row or first column. A primary-source restatement is
given in Launois and Lenagan, arXiv:1207.3613, lines 15-17 of the
introduction.

## Reversed Finite Block

Fix `n>=0` and `N>=0`, and define the finite reversed Hankel block

```text
B_(i,q)^(n,N)=A_(n+i+N-q)(0)
0<=i<=3, 0<=q<=N.
```

For a solid minor of order `1<=k<=4`, direct index substitution and
column reversal give

```text
det B[r:r+k,c:c+k]=epsilon_k*H_(k,n+r+N-c-k+1)(0)
```

where the mapped shift is nonnegative. The completed lambda-zero
theorems give

```text
A_s(0)>0, -H_(2,s)(0)>0, -H_(3,s)(0)>0, H_(4,s)(0)>0
for every integer s>=0.
```

Thus every solid minor, hence every initial minor, of `B^(n,N)` is
positive. The published criterion makes `B^(n,N)` strictly totally
positive.

## Arbitrary Columns

Given `0<=j_1<...<j_k<=N`, use the increasing columns
`N-j_k<...<N-j_1` of `B`. Exact reversal gives

```text
det B[0:k|N-j_k,...,N-j_1]=epsilon_k*R_(k,n)(j_1,...,j_k)
```

and strict total positivity makes the left side positive. Since `N`
can be the largest requested offset,

```text
R_(4,n)(j_1,j_2,j_3,j_4)>0
for every n>=0 and 0<=j_1<j_2<j_3<j_4 at lambda=0.
```

More generally, for every fixed `m`,

```text
epsilon_k=(-1)^binom(k,2)
[epsilon_k H_(k,s)>0 for 1<=k<=m and every s]
  => [epsilon_k R_(k,n)(j_1,...,j_k)>0 for 1<=k<=m].
```

Therefore arbitrary columns require no new analytic theorem once the
contiguous layers are complete. The first new layer is contiguous
order five.

## Exact Audits

The checker verifies the reversal identity through order 4 and 240 solid-block index maps.
The rational benchmark `a_k=1/k!` has 1020 strictly signed arbitrary minors across four shifts and orders one through four.

The lower layers are essential. The exact positive sequence

```text
[10, 9, 29, 18, 21, 25, 3, 16]
```

has

```text
H_(4,0)=288076>0, H_(4,1)=264875>0,
R_(4,0)(0,1,3,4)=-231169<0,
H_(2,0)=209>0 (wrong signed layer).
```

So contiguous order four alone cannot be promoted; the proof genuinely
uses every initial-minor sign through order four.

```text
outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md
outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
