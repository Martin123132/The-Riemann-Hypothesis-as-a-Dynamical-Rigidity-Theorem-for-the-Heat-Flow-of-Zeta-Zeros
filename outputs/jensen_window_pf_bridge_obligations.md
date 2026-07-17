# Jensen-Window PF Bridge Obligations

Date: 2026-07-16

Status: theorem-obligation ledger. This is not a proof of PF-infinity, Jensen
hyperbolicity, Laguerre-Polya membership, RH, or `Lambda <= 0`; it decomposes
the open Jensen-window PF bridge into exact reformulations, finite evidence,
open obligations, and rejection tests.

## Purpose

The target `target_jensen_window_pf_bridge` is concrete but still too broad to
attack as one sentence. This ledger splits it into checkable obligations so a
future proof attempt can be rejected or advanced at the right layer.

Machine-readable ledger:

```text
work/rh_compute/results/jensen_window_pf_bridge_obligations.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py
```

Current result:

```text
validated Jensen-window PF bridge obligations: 11 obligations, 0 issues, 3 open obligations
```

## Exact And Finite Layers

```text
jwpf_01_window_pf_jensen_equivalence:
  finite PF-infinity of B^{d,n,0} is the Toeplitz-total-positivity form of
  the real-nonpositive-root condition for one fixed Jensen window

jwpf_02_degree2_signed_hankel_contact:
  degree 2 matches the m=1 signed-Hankel condition exactly

jwpf_03_low_degree_extra_toeplitz_obligations:
  degrees 3 and 4 introduce extra banded Toeplitz obligations

jwpf_04_current_finite_pf_sturm_evidence:
  current Arb/Sturm/PF manifests validate finite evidence only
```

These rows organize what we already know. None has `would_close_target=true`.

## Rejected Antecedent And Open Bridge Obligations

```text
jwpf_05_all_order_shifted_sign_consistency:
  rejected: the actual endpoint sequence has Q_(10,n)(-100)<0 for
  n=0,1,2,3

jwpf_05b_weaker_xi_specific_antecedent:
  identify a weaker Xi/Phi-specific condition that survives those order-ten
  failures and is provably satisfied for every required degree and shift

jwpf_06_sign_regular_to_jensen_pf_conversion:
  legacy id: convert the weaker jwpf_05b structure, without all-shift
  signed-Hankel positivity, into every binomially weighted Jensen-window
  Toeplitz conclusion or directly into Jensen hyperbolicity

jwpf_07_binomial_weight_and_shift_uniformity:
  handle binomial weights binom(d,j) and all shifts n uniformly
```

The central open theorem is `jwpf_06_sign_regular_to_jensen_pf_conversion`.
Its legacy identifier is retained for dependency stability; its admissible
antecedent is now `jwpf_05b_weaker_xi_specific_antecedent`, not the rejected
all-order signed-Hankel hierarchy.
The checker permits `would_close_target=true` only on open theorem-obligation
rows, never on finite evidence or countermodel rows.

The theorem machinery audit for this row is:

```text
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
```

Current audit result:

```text
validated Jensen-window PF theorem machinery fit matrix: 7 rows, 0 issues, 0 ready-to-apply rows
```

The structural ansatz workbench for this row is:

```text
outputs/jensen_window_pf_structural_ansatz_matrix.md
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
```

Current ansatz result:

```text
validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows
```

## Downstream Conditional Layer

```text
jwpf_08_jensen_to_laguerre_polya_limit:
  once all-degree/all-shift Jensen hyperbolicity is proved noncircularly,
  document the limiting theorem and the normal-family/growth hypotheses for
  this heat-flow normalization
```

This row is conditional. It cannot be invoked before the all-degree Jensen
theorem is actually proved.

## Rejection Tests

```text
jwpf_09_finite_rectangle_promotion_rejected:
  finite Jensen-window PF/Sturm rectangles cannot be promoted to all-shift
  Jensen hyperbolicity by wording alone

jwpf_10_ordinary_coefficient_pf_route_separated:
  PF evidence for c_k=A_k/k! is a separate coefficient route and does not by
itself prove every binomially weighted Jensen window
```

The rigorous route counterexample is recorded in:

```text
outputs/jensen_window_pf_endpoint_order10_counterexample.md
work/rh_compute/results/jensen_window_pf_endpoint_order10_counterexample.json
python work/rh_compute/scripts/check_jensen_window_pf_endpoint_order10_counterexample.py
```

## Integration Points

This obligation ledger is tied to:

```text
outputs/jensen_window_pf_bridge_target.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
outputs/jensen_window_pf_structural_ansatz_matrix.md
outputs/sign_regularity_theorem_fit_matrix.md
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_sign_regularity_theorem_fit_matrix.py
```

## Boundary

Passing this checker means the Jensen-window bridge has a reproducible
obligation decomposition. It does not prove any open obligation, and it does
not promote finite signed-Hankel, Jensen-window, Sturm, PF, or coefficient-PF
evidence into a proof of `Lambda <= 0`.
