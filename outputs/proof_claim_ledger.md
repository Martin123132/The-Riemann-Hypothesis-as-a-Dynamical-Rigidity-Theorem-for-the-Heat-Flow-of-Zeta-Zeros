# Proof Claim Ledger

Date: 2026-07-05

Status: claim-classification ledger. This is not a proof of RH or `Lambda <= 0`; it records which claims are exact, finite, diagnostic, open, rejected, or hygiene gates.

## Purpose

The proof programme now has many artifacts. This ledger keeps the main claims
separated by type so finite evidence and diagnostics cannot be silently
promoted into an all-order theorem.

Machine-readable ledger:

```text
work/rh_compute/results/proof_claim_ledger.json
```

Validator:

```text
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

Current result:

```text
validated proof-claim ledger: 27 claims, 0 issues, 6 open theorem targets
```

## Categories

```text
exact_lemma:
  exact noncircular identities already available for proof development

finite_certificate:
  promoted finite interval or manifest-backed certificates

diagnostic:
  finite necessary-condition or theorem-search diagnostics

algebraic_reindexing:
  exact algebraic translations such as Toeplitz/Jacobi-Trudi

countermodel_gate:
  executable proof-safety obstructions

theorem_target:
  open bridge theorem needed before any proof promotion

forbidden_promotion:
  invalid proof step rejected by countermodel or logic

hygiene_gate:
  reproducibility, language, status, and reference integrity checks
```

## Current Open Targets

```text
target_direct_coefficient_pf:
  prove all Toeplitz minors of c_k(0)=mu_{2k}(0)/(2k)! are nonnegative

target_signed_hankel_jensen_bridge:
  prove that all-order signed-Hankel/sign-consistency structure for A_k,
  with Xi/Phi-specific analytic hypotheses, implies all-degree/all-shift
  Jensen hyperbolicity, Laguerre-Polya membership, or Lambda <= 0
  specification: outputs/signed_hankel_jensen_bridge_target.md

target_jensen_window_pf_bridge:
  prove that every binomially weighted Jensen window
  B^{d,n,0}_j=binom(d,j) A_{n+j}(0) is finite PF-infinity, or prove this
  from all-order shifted reshaped-Hankel sign consistency
  specification: outputs/jensen_window_pf_bridge_target.md

target_schur_positive_specialization:
  construct a positive Schur/Edrei-Thoma specialization h_k -> d_k(0)
  specification: outputs/positive_schur_specialization_target.md

target_positive_determinant_integral:
  derive a positive determinant integral formula for every structurally
  nonzero Toeplitz minor

target_edrei_log_power_representation:
  prove an all-order positive Edrei log-power representation for H_0
```

Each is explicitly recorded as `open_target`. The checker rejects theorem
targets that lack a current blocker or required upgrade.

## New Finite Diagnostic

```text
hankel_sign_consistency_reduction_point_audit:
  exact-rationalized cache point audit for the Grussler-Damm reshaped-Hankel
  finite condition; validates five lambdas, k=2..5, N=18; not an interval
  certificate or all-order sign-consistency proof

hankel_sign_consistency_reduction_finite_certificate:
  Arb/enclosure-backed finite certificate for the same reshaped-Hankel
  frontier; validates 689,795 finite minors for k=2..7, N=20; not an all-order
  sign-consistency theorem

shifted_hankel_sign_consistency_finite_certificate:
  Arb/enclosure-backed finite certificate for shifted reshaped-Hankel blocks;
  staircase manifest validates 3,154,515 finite minors for shifts n=0..20:
  k=2..5 at N=18, k=6 at N=16, k=7 at N=15, and k=8 at N=14; not an
  all-shift or all-order sign-consistency theorem

jensen_hankel_bridge_algebra_gate:
  exact degree-2 signed-Hankel/Jensen identity plus a positive rational
  degree-3 countermodel; blocks promotion from finite low-order reshaped signs
  to Jensen hyperbolicity

jensen_window_pf_obligation_algebra_gate:
  exact low-degree Jensen-window PF obligation algebra; degree 2 matches the
  signed-Hankel threshold, while degree 3 and degree 4 introduce additional
  banded Toeplitz obligations and finite low-order countermodel failures

arb_jensen_window_pf_obligation_diagnostic:
  Arb/enclosure-backed finite diagnostic for selected degree-3/4
  Jensen-window contiguous Toeplitz determinants; validates 1,470/1,470
  finite determinants for the five-lambda grid and shifts n=0..20

arb_jensen_window_sturm_hyperbolicity_diagnostic:
  Arb/Sturm finite diagnostic for selected degree-3/4 Jensen-window
  positive-root counts; validates 210/210 finite root-count rows for the
  five-lambda grid and shifts n=0..20

arb_jensen_window_sturm_d5_hyperbolicity_diagnostic:
  Arb/Sturm finite diagnostic for selected degree-5 Jensen-window
  positive-root counts; validates 105/105 finite root-count rows for the
  five-lambda grid and shifts n=0..20

jensen_window_sturm_pf_finite_consequence:
  finite window-by-window PF consequence of the Arb/Sturm root-count
  manifests; validates 315/315 checked Jensen windows across degrees
  d=3,4,5, five lambdas, and shifts n=0..20

jensen_window_pf_bridge_obligation_ledger:
  theorem-obligation hygiene gate for target_jensen_window_pf_bridge;
  validates 10 exact, finite, open, conditional, rejected, and route-separated
  obligations, with 3 open obligations and finite rows blocked from closing
  the target

signed_hankel_jensen_dependency_graph:
  dependency hygiene gate connecting finite evidence, countermodel gates, open
  theorem targets, and the `lambda_le_0_goal` node; validates that finite and
  diagnostic nodes only support open targets and have no direct proving edge
  to the not-proved conclusion

target_jensen_window_pf_bridge:
  theorem target that reformulates all-degree/all-shift Jensen hyperbolicity
  as finite PF-infinity of every binomially weighted window
  B^{d,n,0}_j=binom(d,j) A_{n+j}(0); open and not proved

countermodel_gates:
  validates 11 proof-safety examples, including local heat birth, finite
  coefficient-prefix promotion, finite Jensen-window rectangle promotion,
  finite Schur-prefix promotion, finite signed-Hankel grid promotion, finite
  moment-recurrence promotion, and Stieltjes/Hankel-to-Toeplitz promotion traps
```

## Boundary

Passing the ledger checker means the corpus has not confused finite
certificates, diagnostics, or algebraic translations with the missing bridge
theorem. It does not prove any open target.
