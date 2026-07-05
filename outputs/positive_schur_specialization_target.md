# Positive Schur Specialization Target

Date: 2026-07-05

Status: theorem target note. This is not a proof of PF-infinity, Laguerre-Polya membership, RH, or `Lambda <= 0`.

## Purpose

This note isolates the currently best algebraic coefficient-PF bridge target:

```text
construct a positive Schur/Edrei-Thoma specialization h_k -> d_k(0)
```

without using finite Toeplitz evidence, ASW/Edrei, Jensen hyperbolicity, RH, or
`Lambda <= 0` as hidden hypotheses.

## Objects

The coefficient-PF route uses:

```text
c_k(lambda) = A_k(lambda) / k! = mu_{2k}(lambda) / (2k)!
d_k(lambda) = c_k(lambda) / c_0(lambda)
H_lambda(z) = sum_{k >= 0} d_k(lambda) z^k
```

The endpoint statement is:

```text
all upper-triangular Toeplitz minors of c_k(0) are nonnegative
```

which is the same as PF-infinity for the one-sided coefficient sequence.

## Exact Algebra Gate

For row and column sets:

```text
R = (r_1 < ... < r_m)
Q = (q_1 < ... < q_m)
```

the exact reindexing gate proves that every structurally nonzero Toeplitz
minor:

```text
det(c_{q_j-r_i})_{i,j=1}^m
```

can be written as:

```text
c_0^m s_{lambda/mu}(d)
```

where `s_{lambda/mu}(d)` is the Jacobi-Trudi skew-Schur determinant under the
specialization:

```text
h_k -> d_k.
```

The executable algebra gate is:

```text
python work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py
outputs/toeplitz_jacobi_trudi_bridge_note.md
```

Current result:

```text
validated Toeplitz/Jacobi-Trudi reindexing: N=10, orders<=5, 124129 minors, 39094 nonzero maps, 85035 structural zeros
```

This validates the translation. It does not prove positivity of any infinite
family of skew-Schur evaluations.

## Target S

Prove one of the following noncircular sufficient statements.

### S1: Positive Specialization

Construct a nonnegative Edrei-Thoma or Schur-positive specialization `alpha`
such that:

```text
h_k(alpha) = d_k(0) for every k >= 0.
```

This would imply:

```text
s_{lambda/mu}(d) >= 0
```

for every skew shape and hence all Toeplitz minors of `c_k(0)` are
nonnegative.

### S2: Planar Network

Construct a locally finite planar network with nonnegative edge weights whose
path matrix satisfies:

```text
P_{r,q} = d_{q-r}(0)
```

for the Toeplitz support convention. Lindstrom-Gessel-Viennot would then give
nonnegative minors.

### S3: Production Matrix Or Continued Fraction

Construct a production matrix, continued fraction, or equivalent total
positivity certificate whose hypotheses are verified directly for `d_k(0)` and
whose conclusion is Schur positivity.

### S4: Positive Determinant Integral

For every structurally nonzero Toeplitz minor, derive:

```text
det(c_{q_j-r_i}(0)) = integral W_{R,Q}(x) dnu(x)
W_{R,Q}(x) >= 0
```

with the measure and integrand defined independently of the desired
PF-infinity conclusion.

## Evidence That May Guide The Search

Finite certificates may guide conjecture selection:

```text
python work/rh_compute/scripts/check_toeplitz_certificate_manifest.py
python work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py
python work/rh_compute/scripts/check_edrei_log_sign_manifest.py
python work/rh_compute/scripts/check_edrei_power_hankel_manifest.py
```

Relevant notes:

```text
outputs/finite_toeplitz_certificate_status.md
outputs/toeplitz_jacobi_trudi_bridge_note.md
outputs/edrei_log_sign_diagnostic.md
outputs/missing_coefficient_pf_theorem.md
outputs/sign_regularity_theorem_fit_matrix.md
```

These artifacts are evidence and falsification pressure. They are not a proof
of Target S.

## Kill Gates

Reject a proposed Target S proof if it:

```text
1. uses only finitely many Toeplitz, Schur, Hankel, or Edrei-log checks;
2. invokes ASW/Edrei after assuming PF-infinity or real-negative zeros;
3. proves only Stieltjes or Hankel moment positivity;
4. assumes Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0;
5. identifies the direct de Bruijn-Newman kernel PF problem with the coefficient PF problem;
6. leaves a signed determinant or signed integral with no independent positivity theorem;
7. proves positivity only for a bounded shape/order/degree family.
```

Executable countermodel support:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
outputs/countermodel_library.md
```

These gates block finite-prefix promotion, moment-positivity promotion,
finite Schur-prefix promotion, local-repulsion promotion, and closure arguments
whose hypotheses do not imply Toeplitz total positivity.

The exact Schur-prefix trap preserves `h_0..h_6` from the PF-infinity model
`h_k = 1/k!`, validates `2,940` finite Toeplitz/Schur grid tests, and then
chooses a positive `h_7` that makes the Jacobi-Trudi determinant
`s_(6,6) = h_6^2 - h_5 h_7` negative. This shows that even a finite prefix
copied from a genuinely positive specialization cannot replace an all-order
construction.

## Current Blocker

No positive specialization, planar network, production matrix, continued
fraction, or determinant integral formula is currently known for:

```text
d_k(0) = c_k(0) / c_0(0).
```

The missing theorem is exactly the construction of such an object, or an
equivalent all-order total-positivity proof for the coefficient sequence.

## Acceptance Condition

This target can be promoted only when a proposed construction provides:

```text
1. exact definitions for all parameters, weights, measures, or kernels;
2. a proof that those objects are nonnegative in the required domain;
3. an identity recovering d_k(0), every Toeplitz minor, or every skew-Schur
   evaluation;
4. convergence and limiting arguments sufficient for all orders;
5. no appeal to RH, Lambda <= 0, PF-infinity, Laguerre-Polya membership, or
   Jensen hyperbolicity as an input.
```

Until then, Target S remains an open theorem target.
