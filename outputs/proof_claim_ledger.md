# Proof Claim Ledger

Date: 2026-07-04

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
validated proof-claim ledger: 17 claims, 0 issues, 5 open theorem targets
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

target_schur_positive_specialization:
  construct a positive Schur/Edrei-Thoma specialization h_k -> d_k(0)

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
  frontier; validates 62,985 finite minors; not an all-order
  sign-consistency theorem
```

## Boundary

Passing the ledger checker means the corpus has not confused finite
certificates, diagnostics, or algebraic translations with the missing bridge
theorem. It does not prove any open target.
