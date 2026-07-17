# Jensen-Window PF Endpoint Deep-Schur Coordinate

Date: 2026-07-16

Status: exact normalized deep-Schur coordinate and rigorous endpoint
PF separation. The proposed all-order deep rectangle hierarchy is
rejected by an order-ten endpoint counterexample. This does not settle
the separate Jensen-window PF bridge. This is not a proof of RH or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_endpoint_deep_schur_coordinate.json
python work/rh_compute/scripts/jensen_window_pf_endpoint_deep_schur_coordinate.py
python work/rh_compute/scripts/check_jensen_window_pf_endpoint_deep_schur_coordinate.py
```

## Normalized Endpoint Specialization

Use the genuine complete-homogeneous specialization

```text
h_k=A_k(-100)/A_0(-100) for k>=0, h_k=0 for k<0, h_0=1
```

The normalization is essential: ordinary Schur theory has `h_0=1`.
Every size-`m` determinant in the unnormalized `A_k` acquires only the
strictly positive common factor `A_0(-100)^m`.

## Rectangular Identity

Reverse the `m` Hankel columns. Since the reversal has sign
`epsilon_m=(-1)^binom(m,2)`, transpose the resulting matrix, and put
`N=n+m-1`. The Jacobi-Trudi matrix is obtained entry by entry:

```text
Q_(m,n)(-100)=det[A_(n+m-1+i-j)(-100)]_(0<=i,j<m)
Q_(m,n)(-100)=A_0(-100)^m*s_((n+m-1)^m)(h)
```

Thus the candidate static endpoint statement was precisely

```text
s_((N^m))(h)>0 for every m>=10 and N>=m-1
```

This is an exact change of coordinates; the counterexample below
shows that the universal positivity statement is false.

## Arbitrary Columns And The Deep Cone

For increasing offsets `0<=j_0<...<j_(m-1)`, define `R_(m,n)` as in
the existing arbitrary-column transfer theorem. Reversal and transpose
give

```text
lambda_(q+1)=n+j_(m-1-q)+q for 0<=q<m
epsilon_m*R_(m,n)(j_0,...,j_(m-1))=A_0(-100)^m*s_lambda(h)
```

The adjacent partition differences are
`lambda_r-lambda_(r+1)=j_(m-r)-j_(m-r-1)-1>=0`, and
`lambda_m=n+j_0+m-1>=m-1`. Conversely, for every
`lambda_1>=...>=lambda_m>=m-1`, the canonical inverse is

```text
n=lambda_m-(m-1), j_l=lambda_(m-l)-lambda_m+l for 0<=l<m
```

It has `j_0=0` and strictly increasing columns. Hence the arbitrary
signed minors are in bijection with

```text
D_m={lambda_1>=...>=lambda_m>=m-1}; s_lambda(h)=det[h_(lambda_i-i+j)]_(1<=i,j<=m)
```

The threshold is not decorative:

```text
min_(i,j)(lambda_i-i+j)=lambda_m-m+1>=0 exactly on D_m
```

So the deep cone is exactly the Jacobi-Trudi region that never touches
the artificial values `h_k=0` for `k<0`.

## Exact Endpoint Equivalence

The completed orders through nine provide every lower initial minor.
At each fixed higher order, the Gasca-Pena initial-minor theorem then
transfers rectangular positivity to every arbitrary column set. Since
rectangles are themselves deep shapes, this proves

```text
[Q_(m,n)(-100)>0 for every m>=10,n>=0] iff [s_lambda(h)>0 for every m>=10,lambda in D_m]
```

Within this candidate hierarchy, only the rectangles are independent;
the rest of the deep cone would be a structural consequence once all
lower rectangles were available.

## Rigorous PF Separation

The rigorous acb endpoint enclosures for `A_0,...,A_3` were propagated
with exact rational interval arithmetic. They give

```text
s_(1,1,1)(h)=h_1^3-2*h_1*h_2+h_3
interval = [-4.84848642182060969711098547098349E-11,
            -4.84848642182060969711098547098349E-11 ]
```

The upper endpoint is strictly negative. Therefore the actual normalized
endpoint sequence is not `PF_3`, hence not `PF-infinity`. But
`(1,1,1)` lies outside `D_3` because its smallest part is `1<2`.
Thus full Schur positivity imposes a genuinely additional inequality
that is already false here. The ordinary Edrei route is unavailable,
while this failed inequality is not part of the rectangular endpoint
family.

## Deep Rectangle Counterexample

The first-open-order certificate now gives required deep failures:

```text
s_((N^10))(h)<0 for N=9,10,11,12
s_((N^m))(h)>0 for every m>=10 and N>=m-1 is false
```

These are not shallow boundary shapes. At order ten the depth condition
is `N>=9`, so `(9^10)` through `(12^10)` belong to the required cone.

## Primary-Literature Fit Gate

### gasca_pena_initial_minors

Classification: `applicable_transfer_only`.

Positive initial minors characterize strict total positivity of each finite reversed block. This is exactly the fixed-order rectangle-to-arbitrary-shape transfer already used here.

It does not supply: new positivity for any open rectangle.

Source: https://doi.org/10.1016/0024-3795(92)90226-Z

### edrei_pf_classification

Classification: `strictly_too_strong`.

The theorem classifies full Toeplitz total positivity, whereas the endpoint specialization has a rigorously negative order-three Toeplitz minor outside the deep cone.

It does not supply: restricted deep-cone positivity.

Source: https://doi.org/10.4153/CJM-1953-010-3

### pena_hankel_toeplitz_orientation

Classification: `structural_orientation_support`.

Finite Hankel column reversal produces the signature (-1)^binom(k,2) in the associated Toeplitz matrix. It supports the orientation coordinate but assumes total positivity on the Hankel side and does not prove the Xi endpoint signs.

It does not supply: the open all-order endpoint hierarchy.

Source: https://doi.org/10.3390/math13142278

### kushel_eventual_total_positivity

Classification: `coordinate_mismatch`.

Eventual total positivity there means total positivity of all sufficiently high matrix powers, not positivity of Toeplitz minors whose Jacobi-Trudi shapes lie beyond a moving depth boundary.

It does not supply: an eventual-depth Toeplitz theorem.

Source: https://doi.org/10.7153/oam-09-56

No direct closing theorem was found in this audited primary-source
set. That is a bounded route audit, not a claim about all literature.
The phrase `eventual total positivity` is especially unsafe here:
in the audited matrix literature it refers to high matrix powers,
whereas our coordinate concerns minors beyond a moving Toeplitz depth.

## Conditional Heat Composition

The abstract implication into the heat theorem remains exact:

```text
deep endpoint positivity => all-order endpoint hierarchy => all-order signed-Hankel positivity on -100<=lambda<=0
```

For the actual endpoint sequence its antecedent is false, so this
composition is unavailable as an RH route.

The separate Jensen-window obligation remains:

```text
deep Schur positivity of h_k=A_k(-100)/A_0(-100) is not Schur positivity of the binomially weighted Jensen-window sequences
```

## Machine Audit

The generator records `4096` rectangle-map checks,
`984` arbitrary-column checks,
`1023` bounded inverse-partition checks, and
symbolic determinant residuals through order `5`.
The finite checks audit implementation. General validity comes from
the explicit entrywise maps and finite determinant identities above.
