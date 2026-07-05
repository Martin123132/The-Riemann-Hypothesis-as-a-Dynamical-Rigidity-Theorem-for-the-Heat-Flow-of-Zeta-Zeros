# RH Dynamical Rigidity Proof Programme

Date: 2026-07-03

Status: proof-programme roadmap. This is not a proof of RH or `Lambda <= 0`; it organizes exact results, evidence, gates, and missing bridge theorems.

Goal:

Move the current bundle from exploratory notes and numerical evidence toward a referee-legible proof programme for RH via the de Bruijn-Newman family.

This roadmap does not claim the proof is complete. It identifies the exact bridge theorem that would be needed and organizes work around proving, weakening, or replacing it.

## North Star

In the de Bruijn-Newman framework, RH is equivalent to

```text
Lambda <= 0.
```

Rodgers-Tao proves the complementary statement

```text
Lambda >= 0.
```

Therefore this programme must prove

```text
Lambda <= 0
```

without assuming RH at `lambda = 0`.

Equivalently, it must exclude every positive Newman boundary:

```text
lambda_* > 0,
Xi_{lambda_*}(t_*) = 0,
partial_t Xi_{lambda_*}(t_*) = 0.
```

## Current Status

### Usable exact material

- Simple real zeros track by the implicit function theorem.
- The zero velocity identity

```text
t_i' = -F_lambda/F_t
```

is exact wherever the zero is simple.

- For the Newman sign,

```text
F_lambda = -F_tt,
```

the real-zero velocity has local singular structure

```text
t_i' = 2 sum_{j != i} 1/(t_i - t_j) + regular external field.
```

- Weierstrass preparation gives square-root local birth/death normal form at a simple double zero.
- Already-real finite isolated clusters do not collapse in the increasing-lambda repulsive direction, provided the external field is locally Lipschitz.

### Strong numerical material

- Finite-window reduced repulsive flows show robust low-dimensional geometry.
- The leading mode aligns with Hadamard repulsion.
- A secondary density-balancing mode appears stable across many windows.
- Jensen/Hankel diagnostics suggest a signed regularity pattern rather than ordinary positivity.

### Not yet proved

- Reduced finite-window dynamics are not proved equivalent to the full Xi_lambda zero dynamics.
- Rank-2 closure is not proved for the full infinite zero system.
- Local repulsion does not exclude a positive Newman boundary.
- Tail positivity in the finite pairwise model is false as stated; the finite pairwise tail is compressive.
- Candidate arguments assuming RH at lambda = 0 are circular for proving RH.

## Main Bridge Problem

Prove one of the following.

### Bridge A: No Positive Newman Boundary

Show there is no `lambda_* > 0` and real `t_*` satisfying

```text
Xi_{lambda_*}(t_*) = 0,
partial_t Xi_{lambda_*}(t_*) = 0.
```

This is the most direct route to `Lambda <= 0`.

### Bridge B: Laguerre-Polya Membership

Show `Xi_lambda` belongs to the Laguerre-Polya class for every `lambda > 0`.

This is equivalent in spirit to proving all zeros are real in the positive-lambda regime.

### Bridge C: Signed Total Positivity / Sign-Regularity

Turn the observed signed Hankel and Jensen polynomial patterns into a uniform theorem strong enough to imply real-rootedness of the limiting entire function.

This is currently the most promising continuation of the later notebook evidence.

### Bridge D: Xi-Specific Comparison Principle

Prove a global comparison theorem for neighboring gaps in the actual Xi_lambda flow. It must use the specific structure of Xi_lambda and cannot rely only on the local `4/g` repulsion law.

## Workstreams

### Workstream 1: Formal Core Cleanup

Deliverable:

```text
formal_core.md
```

Tasks:

- Fix all sign conventions: heat parameter, Newman lambda, reduced repulsive parameter.
- Define Xi_lambda, Phi, Lambda, and the zero indexing convention.
- State exact lemmas only where assumptions are explicit.
- Separate theorem, conditional theorem, numerical observation, and conjecture.
- Remove or mark all circular claims.

Success criterion:

A skeptical reader can tell exactly what has been proved and what remains conjectural.

### Workstream 2: Bridge Theorem Selection

Deliverable:

```text
bridge_theorem_options.md
```

Tasks:

- State Bridges A-D precisely.
- For each bridge, list known supporting evidence, obstacles, and possible existing theorem families.
- Choose one route as the first serious attack.

Recommended first route:

```text
Bridge C: signed total positivity / sign-regularity.
```

Reason:

The local zero-flow route has a known obstruction: square-root birth is locally compatible with repulsion. The signed Hankel/Jensen evidence is less exhausted and may connect to established real-rootedness machinery.

### Workstream 3: Signed Hankel/Jensen Audit

Deliverable:

```text
signed_hankel_jensen_audit.md
```

Tasks:

- Extract the exact coefficient sequence being tested.
- Define all ordinary and shifted Hankel matrices used in the notebooks.
- Re-run finite tests with reproducible scripts.
- Replace floating checks with interval or rational-enclosure checks where possible.
- Identify the exact theorem that would turn these signs into real-rootedness.

2026-07-04 status update:

- The reproduced signed-Hankel grid has been extended and is now enclosure-backed Arb certified for `lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}`, `m = 0..19`, `s = 0..24`: `2,500/2,500` sign-separated positive.
- The certificate is finite only; the missing workstream theorem is still an all-order bridge from signed-Hankel/sign-regularity data to Jensen hyperbolicity or Laguerre-Polya membership.
- Certificate details are in `outputs/signed_hankel_certificate_status.md`.
- The theorem ecosystem fit/misfit matrix is now recorded in `outputs/sign_regularity_theorem_fit_matrix.md`. Its current conclusion is that ASW/Edrei, Schoenberg, Jensen, and Hankel moment theorems are endpoint or necessary-condition tools unless we independently construct a positive Schur/Edrei-Thoma specialization, an Edrei log-power representation, or a positive determinant integral formula for the zeta coefficient sequence.
- The Grussler-Damm sign-consistency route has been sharpened: shifted-principal signed-Hankel minors are not enough, so a reshaped-Hankel audit checks all k-column minors of the k-row matrix `(A_{i+j})` for `k = 2..5`, `N = 18`, across the five-lambda grid. The exact-rationalized cache audit is now accompanied by an Arb/enclosure-backed finite certificate validating `62,985/62,985` finite minors. This is still not an all-order sign-consistency theorem.
- The Edrei log-power representation now has an Arb moment-recurrence scout in `outputs/edrei_moment_quadrature_scout.md`: orders `2..12` are positive across the five-lambda grid, while a broader order `2..20` scout is inconclusive from order `13` onward and has no negative rows. This is a finite theorem-search diagnostic, not an all-order representation proof.

Success criterion:

We either identify a known theorem that nearly applies, or we isolate the missing hypothesis.

Failure criterion:

The signed pattern fails at higher degree or under rigorous enclosures.

### Workstream 4: Toy Model and Countermodel Library

Deliverable:

```text
countermodel_library.md
```

Tasks:

- Keep the minimal square-root birth heat-equation model.
- Build additional toy models satisfying the same local zero-flow laws but failing RH-like conclusions.
- Use these to prevent false bridge lemmas.

2026-07-04 status update:

- The countermodel library now exists at `outputs/countermodel_library.md`.
- The executable gate `work/rh_compute/scripts/countermodel_gate_examples.py` validates the local heat-birth model and finite-prefix extension traps for the current five-lambda coefficient cache.
- The new finite-prefix traps show that preserving all current coefficients through `k = 32` still allows positive one-term extensions that break order-2 Toeplitz/PF, `m = 1` signed-Hankel, and degree-2 Jensen hyperbolicity. This blocks any manuscript step that tries to promote finite-prefix evidence into an all-order theorem without additional structure.
- The gate also contains an exact finite shifted-principal signed-Hankel grid trap: `a_k = 1/k!` passes `m = 0..4`, `s = 0..8`, then a positive `a_17` extension breaks the next shifted `m = 1` signed-Hankel/Jensen test. This blocks finite-grid promotion unless a separate Hankel sign-consistency reduction and all-order proof are supplied.

Success criterion:

Every proposed proof step is tested against toy models before entering the manuscript.

### Workstream 5: Rigorous Numerics

Deliverable:

```text
rigorous_numerics_plan.md
```

Tasks:

- Convert key finite Jensen/Hankel tests to interval arithmetic.
- Record precision, truncation, tail bounds, and root-certification method.
- Produce machine-checkable certificates for finite claims.

Success criterion:

Finite computational claims become reliable supporting lemmas, not informal plots.

### Workstream 6: Manuscript Reconstruction

Deliverable:

```text
paper_draft.md
```

Suggested title:

```text
The Riemann Hypothesis as a Dynamical Rigidity Programme:
Local Zero Motion, Reduced Repulsive Flow, and the Newman Boundary Gap
```

Structure:

1. Introduction and status
2. de Bruijn-Newman setup
3. Exact zero-motion lemmas
4. Local square-root boundary normal form
5. Why local repulsion is insufficient
6. Reduced-flow numerics
7. Signed Hankel/Jensen evidence
8. Bridge theorem / open problem
9. Appendices with reproducible computations

Success criterion:

The manuscript becomes something a serious analyst could engage with, even before the bridge theorem is proved.

## First Sprint

The first sprint should be narrow:

1. Create `formal_core.md`.
2. Create `bridge_theorem_options.md`.
3. Extract the signed Hankel/Jensen definitions from the notebooks.
4. Re-run the finite tests in a clean script.
5. Search for the closest known theorem connecting signed Hankel patterns, Jensen polynomials, total positivity, and Laguerre-Polya class.

## What We Should Not Do Next

- Do not claim a proof of RH yet.
- Do not add more rank-2 numerical evidence until the bridge theorem is clearer.
- Do not rely on a local no-collision lemma to exclude positive Newman birth.
- Do not use any argument that assumes RH at lambda = 0.

## Practical Decision Rule

Every future step should answer one of these:

1. Does it prove a non-circular bridge theorem?
2. Does it rigorously support a hypothesis inside a bridge theorem?
3. Does it rule out a tempting but false route?
4. Does it make the manuscript clearer and harder to misread?

If not, it is probably a distraction.
