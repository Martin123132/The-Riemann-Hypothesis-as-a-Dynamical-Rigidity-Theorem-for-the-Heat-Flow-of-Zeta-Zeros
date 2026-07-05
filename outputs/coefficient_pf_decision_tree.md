# Coefficient PF Decision Tree

Date: 2026-07-04

Status: proof-programme triage note. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.

## Starting Point

The live coefficient object is:

```text
c_k(lambda) = A_k(lambda) / k! = mu_{2k}(lambda) / (2k)!
F_lambda(s) = sum_{k >= 0} c_k(lambda) s^k.
```

The finite Arb evidence supports selected Toeplitz/PF inequalities for this sequence. The missing bridge is:

```text
all upper-triangular Toeplitz minors of {c_k(0)} are nonnegative.
```

By ASW/Edrei, that would put `F_0` in the restricted Laguerre-Polya class, after the normalization checks recorded elsewhere.

## Decision Tree

### Route 1: Direct All-Order PF-Infinity

Try to prove:

```text
{c_k(0)} is PF-infinity.
```

Acceptable evidence:

```text
an exact theorem applying to the whole infinite sequence,
not a finite minor pattern.
```

Kill gate:

```text
find one certified negative Toeplitz minor for c_k(0).
```

### Route 2: Positive Determinantal Integral Formula

Try to prove, for every structurally nonzero row/column pair:

```text
det(c_{q_b-r_a}) = integral W_{R,Q}(x) dnu(x)
W_{R,Q}(x) >= 0.
```

Acceptable evidence:

```text
Cauchy-Binet, Andreief, Karlin-style, or Schur-positive identity
with a manifestly nonnegative integrand/measure.
```

Kill gate:

```text
the candidate identity introduces a signed determinant with no independent
positivity theorem, or assumes real zeros of F_0.
```

### Route 3: Schur-Positive Specialization

Use Jacobi-Trudi:

```text
Toeplitz minors = c_0^m * skew-Schur evaluations under h_k -> d_k,
d_k = c_k/c_0.
```

Try to prove that the specialization:

```text
h_k |-> d_k(lambda)
```

is Schur-positive.

Acceptable evidence:

```text
positive specialization, production matrix, planar network, or measure on partitions
whose complete homogeneous coordinates are exactly d_k(lambda) = c_k(lambda)/c_0(lambda).
```

Kill gate:

```text
the construction proves only Hankel/Stieltjes positivity, not Schur/Toeplitz positivity.
the construction invokes ASW/Edrei after assuming the same PF-infinity conclusion.
```

### Route 4: Closure Or Multiplier Factorization

Try to factor:

```text
c_k(lambda) = a_k(lambda) b_k
```

with:

```text
{a_k(lambda)} PF-infinity
{b_k} PF-preserving as a multiplier
```

or factor `F_lambda` into restricted Laguerre-Polya pieces.

Acceptable evidence:

```text
a named closure theorem whose hypotheses are verified for this exact sequence.
```

Kill gate:

```text
using ordinary moment positivity of mu_{2k} as if it implied coefficient PF-infinity.
```

### Route 5: Heat-Interval Preservation

Try to prove:

```text
{c_k(lambda)} is PF-infinity for lambda in [0, epsilon].
```

Acceptable evidence:

```text
a noncircular variation-diminishing or total-positivity preservation theorem
for the coefficient sequence under the heat deformation.
```

Kill gate:

```text
the proof assumes real-rootedness at lambda = 0, or assumes Lambda <= 0.
```

## Current Priority Order

1. Route 2: seek a positive determinant integral identity.
2. Route 3: translate the Toeplitz minors into skew-Schur positivity and search for a positive specialization.
3. Route 4: audit multiplier/closure options carefully.
4. Route 5: attempt only if the heat-flow structure gives a genuine preservation theorem.
5. Route 1: keep as the umbrella statement, but do not treat it as a proof method by itself.

## Evidence Policy

Finite computations may:

```text
falsify a route;
identify fragile minors;
calibrate conjectures;
prioritize theorem search.
```

Finite computations may not:

```text
prove PF-infinity by extrapolation;
replace ASW/Edrei hypotheses;
justify any RH or Lambda <= 0 claim without the missing all-order bridge.
```

## Next Theorem-Search Action

For each structurally nonzero Toeplitz minor:

```text
det(c_{q_b-r_a}) = det(mu_{2(q_b-r_a)} / (2(q_b-r_a))!)
```

first map it to its exact skew shape using:

```text
work/rh_compute/scripts/toeplitz_jacobi_trudi_map.py
```

then attempt to rewrite the corresponding skew-Schur evaluation as a positive combination of:

```text
Hankel moment determinants,
Vandermonde-square integrals,
skew-Schur evaluations,
or planar-network path weights.
```

Any rewrite that leaves a signed sum with no independent positivity theorem is not yet a bridge.

The exact reindexing and source anchors for this step are recorded in:

```text
outputs/toeplitz_jacobi_trudi_bridge_note.md
```

The broader theorem ecosystem fit/misfit audit is recorded in:

```text
outputs/sign_regularity_theorem_fit_matrix.md
```
