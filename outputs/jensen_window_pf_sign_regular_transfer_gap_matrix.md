# Jensen-Window PF Sign-Regular Transfer Gap Matrix

Date: 2026-07-16

Status: finite theorem-search diagnostic. This is not a proof of
Jensen-window PF-infinity, Jensen hyperbolicity, Laguerre-Polya
membership, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_sign_regular_transfer_gap_matrix`.

Proof boundary: this artifact combines exact low-degree algebra and
finite countermodel gates. It does not prove the missing weaker Xi/Phi
to Jensen-window PF transfer theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_sign_regular_transfer_gap_matrix.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_sign_regular_transfer_gap_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_sign_regular_transfer_gap_matrix.py
```

Current result:

```text
validated Jensen-window PF sign-regular transfer gap matrix: 9 transfer rows, 2 countermodel gates, 3 open requirements, 3 rejected shortcuts, 0 ready-to-apply rows, 0 issues
```

## Countermodel Gates

Shared positive rational sequence:

```text
A_0..A_4 = 1, 33/40, 429/640, 4719/12800, 4719/35840
```

It passes the finite reshaped-Hankel signs used by the bridge algebra gate:

```text
finite reshaped signs pass: True
degree-3 Jensen discriminant: -2476526481/3276800000
degree-3 Jensen hyperbolicity breaks: True
```

It also passes selected low-order Jensen-window Toeplitz tests while
later obligations fail:

```text
d=3 selected Toeplitz minors positive: True
d=3 first negative contiguous size: 8
d=3 first negative contiguous determinant: -435846079534239/104857600000000
d=4 selected Toeplitz minors positive: True
d=4 quartic discriminant: -668519580649275927/359661568000000000
d=4 first negative contiguous size: 6
d=4 first negative contiguous determinant: -229760849637/28672000000
```

## Transfer Contract

A usable theorem must supply all of:

```text
1. a weaker Xi/Phi-specific antecedent satisfied despite the negative order-ten endpoint minors
2. noncircular Xi/Phi-specific structure absent from arbitrary positive sequences
3. binomial-weight and shift uniformity for every Jensen-window Toeplitz minor
4. no endpoint PF, Jensen hyperbolicity, Laguerre-Polya, RH, or Lambda <= 0 assumptions
```

Rows:

```text
srgt_01_degree2_exact_contact: exact_contact - Degree 2 is exact: the m=1 signed-Hankel condition is the Jensen-window quadratic discriminant threshold.
srgt_02_degree3_new_obligation: exact_obstruction - Degree 3 introduces a cubic discriminant and Toeplitz minors that are not signed reshaped-Hankel minors.
srgt_03_finite_reshaped_hankel_countermodel: rejected_shortcut - A positive rational sequence passes finite reshaped-Hankel signs for k=2,3 at N=3 while its degree-3 Jensen discriminant is negative.
srgt_04_selected_toeplitz_countermodel: rejected_shortcut - The same positive rational sequence has selected low-order degree-3 and degree-4 Toeplitz minors positive while later contiguous Toeplitz minors and discriminants are negative.
srgt_05_weaker_xi_phi_antecedent_requirement: open_requirement - The bridge needs a weaker Xi/Phi-specific antecedent that is actually satisfied despite the negative order-ten signed-Hankel endpoint minors.
srgt_06_xi_phi_specific_transfer_requirement: open_requirement - A viable bridge must add noncircular Xi/Phi-specific structure that is absent from arbitrary positive sequences.
srgt_07_binomial_shift_uniformity_requirement: open_requirement - The transfer must output every binomially weighted Jensen-window Toeplitz matrix for all degrees and shifts.
srgt_08_acceptable_theorem_shape: live_contract - An acceptable theorem may prove a sign-regular-to-Toeplitz transfer, a positive determinant integral, or an Xi/Phi positive-kernel identity that directly gives the Jensen-window PF target.
srgt_09_forbidden_transfer_shortcuts: rejected_shortcut - Finite grids, degree-2 analogy, selected minors, endpoint PF/LP assumptions, and asymptotic Jensen statements are forbidden as bridge replacements.
```

Integration:

```text
outputs/jensen_hankel_bridge_algebra.md
outputs/jensen_window_pf_obligation_algebra.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/jensen_window_pf_bridge_obligations.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
```

Summary:

The signed-Hankel/Jensen bridge is exact only at degree 2. The degree-3/4 countermodel gates show that finite reshaped-Hankel signs and selected low-order Toeplitz positivity do not imply Jensen-window PF-infinity. The all-order signed-Hankel antecedent is itself false at order ten, so a proof needs a weaker Xi/Phi-specific kernel, determinant, or direct Jensen theorem with binomial and shift uniformity.
