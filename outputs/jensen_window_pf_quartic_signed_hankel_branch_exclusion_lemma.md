# Jensen-Window PF Quartic Signed-Hankel Branch-Exclusion Lemma

Date: 2026-07-25

Status: exact Xi quartic branch-exclusion lemma. The remaining outer
threshold is open; this is not a proof of degree-four Xi
hyperbolicity, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.json
python work/rh_compute/scripts/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.py
```

Current result:

```text
validated Jensen-window PF quartic signed-Hankel branch-exclusion lemma: 10 rows, 0 issues, 4 exact identities, 1 imported Xi Hankel theorem, 2 excluded boundary strata, 1 reduced outer threshold, 1 countermodel separation, 1 open outer threshold
```

## Shifted Hankel Normalization

For a fixed shift `n`, write

```text
A_(n+j)=A_n*r_n^j*B_j, r_n=A_(n+1)/A_n,
B_0=B_1=1, B_2=x, B_3=x^2*y, B_4=x^3*y^2*z.
```

Then

```text
D_(3,n)=det[A_(n+i+j)]_(i,j=0..2)
       =A_n^3*r_n^6*det[B_(i+j)]_(i,j=0..2),
det[B_(i+j)]=x**3*(x*y**2*z - x*y**2 - y**2*z + 2*y - 1).
```

The prefactor is strictly positive.

## Quartic Boundary Factorization

At a nondegenerate quartic boundary point, write

```text
P(w)=(1+a*w)^2*(1+b*w)*(1+c*w),
2*a+b+c=4, p=b*c.
```

Substitution of the exact contraction coordinates gives

```text
C=3*a^2-4*a+p=(a-b)*(a-c),
det[B_(i+j)]=-C^3/216.
```

The completed compound-order-three theorem proves

```text
D_(3,n)(lambda)<0
```

for every shift and every finite `lambda>=-100` in the propagated
Xi ratio cone. Therefore `C>0`. The repeated root lies outside the
two simple roots; the middle-root branch `C<0` and triple-root
stratum `C=0` are impossible at a hypothetical Xi quartic contact.

## Reduced Quartic Target

The exact branch-aware condition from the threshold lemma was

```text
C*(u-U(a,p))<=0,
U(a,p)=(-a**2 + 2*a + p)*(3*a**2 - 5*a + 5*p)/(6*p**2).
```

Since the Xi signed-Hankel theorem forces `C>0`, the complete live
first-crossing condition is now only

```text
u<=U(a,p).
```

The tangent triple-root analysis no longer belongs to the live Xi
boundary set. Proving the remaining outer threshold is still a
genuine open step.

## Countermodel Separation

The strong-log-concave local quartic countermodel has

```text
D_3=[9.515769718203299328698745521778012352202713614459947107E-7 +/- 1.86E-62]>0.
```

It is therefore excluded by structure already proved for Xi. This
explains exactly why that model blocks a generic local theorem while
not blocking the signed-Hankel route.

## Proof Boundary

The lemma selects the only possible Xi quartic boundary branch. It
does not prove `u<=U`, quartic heat invariance, an all-degree bridge,
PF-infinity, `Lambda<=0`, or RH.

```text
outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md
outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md
outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md
outputs/jensen_window_pf_strong_logconcave_local_quartic_countermodel.md
```

Summary:

At a quartic double-root boundary, the normalized contiguous order-three Hankel determinant is -((a-b)(a-c))^3/216. The strict Xi sign therefore excludes the middle-root and triple-root strata throughout the target heat interval, leaving only the outer branch and the single unproved threshold u<=U(a,p).
