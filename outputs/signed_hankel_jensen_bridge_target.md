# Signed-Hankel/Jensen Bridge Target

Date: 2026-07-05

Status: theorem target. This is not a proof of RH or `Lambda <= 0`; it is the
precise open bridge statement suggested by the signed-Hankel, reshaped-Hankel,
and Jensen evidence.

## Purpose

The current finite evidence supports a signed-Hankel/sign-consistency pattern
for the exponential/Jensen coefficients:

```text
A_k(lambda) = mu_{2k}(lambda) k! / (2k)!
```

The missing proof step is not another finite grid. It is an all-order theorem
that turns the right sign-regularity structure into all-degree/all-shift
Jensen hyperbolicity, Laguerre-Polya membership, or directly `Lambda <= 0`.

## Objects

For each `lambda`, define:

```text
P_{d,n,lambda}(x)
  = sum_{j=0}^d binom(d,j) A_{n+j}(lambda) x^j
```

The Jensen target is:

```text
for every d >= 1 and n >= 0,
P_{d,n,0}(x) has only real nonpositive zeros.
```

Equivalently, `Q_{d,n,0}(y)=P_{d,n,0}(-y)` has `d` positive real zeros.

The total-positivity reformulation is recorded in:

```text
outputs/jensen_window_pf_bridge_target.md
python work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py
```

It asks for every binomially weighted Jensen window:

```text
B^{d,n,0}_j = binom(d,j) A_{n+j}(0)
```

to be a finite PF-infinity sequence. This is an equivalent Jensen-window
target, not a proof that the target holds.

For the reshaped-Hankel sign-consistency route, define the shifted row block:

```text
R_{k,n}(j_1,...,j_k)
  = det(A_{n+i+j_l})_{i=0..k-1, l=1..k}
```

where:

```text
0 <= j_1 < ... < j_k.
```

The all-order sign-consistency hypothesis would be:

```text
(-1)^(k(k-1)/2) R_{k,n}(j_1,...,j_k) > 0
```

for every `k >= 1`, `n >= 0`, and every strictly increasing column set, with
the correct weak/zero clauses if structural degeneracies are proved.

The first Arb certificate verifies only a finite unshifted frontier:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0
k = 2..7
N = 20 columns
689,795/689,795 finite minors positive
```

The shifted Arb certificate extends the finite evidence to a bounded shift
range:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0..20
k = 2..5
N = 18 columns
1,322,685/1,322,685 finite minors positive
```

The order-6 shifted Arb certificate adds the next-order finite slice:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0..20
k = 6
N = 16 columns
840,840/840,840 finite minors positive
```

The order-7 shifted Arb certificate adds another bounded next-order slice:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0..20
k = 7
N = 15 columns
675,675/675,675 finite minors positive
```

The order-8 shifted Arb certificate adds the next bounded slice:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0..20
k = 8
N = 14 columns
315,315/315,315 finite minors positive
```

The consolidated shifted staircase checker validates:

```text
3,154,515/3,154,515 finite shifted minors positive
```

All certificates remain finite. They are evidence for the theorem target, not
the all-shift/all-order theorem itself.

## Candidate Theorem B-Star

Prove a theorem of the following shape:

```text
Assume:
  1. A_k(lambda) are the Xi/Phi heat-flow coefficients, not an arbitrary
     positive sequence.
  2. A_k(lambda) satisfy the all-order shifted reshaped-Hankel
     sign-consistency condition above.
  3. The associated entire functions have the required convergence, parity,
     growth, and heat-flow normalization.

Then:
  for lambda = 0, all Jensen polynomials P_{d,n,0} are hyperbolic;
  equivalently H_0 belongs to the required Laguerre-Polya class;
  hence Lambda <= 0.
```

Stronger acceptable version:

```text
the same conclusion holds uniformly for lambda in [0, epsilon]
```

or:

```text
the sign-consistency structure excludes positive Newman boundary birth
directly.
```

## Exact Low-Degree Gates

Degree 1 is automatic from positivity.

Degree 2 is exact:

```text
P_{2,n}(x) = A_n + 2 A_{n+1} x + A_{n+2} x^2
Delta_2 = 4(A_{n+1}^2 - A_n A_{n+2})
Delta_2 = -4 det [[A_n, A_{n+1}], [A_{n+1}, A_{n+2}]]
```

Thus the `m = 1` signed-Hankel condition is exactly degree-2 Jensen
hyperbolicity for positive coefficients.

Degree 3 is the first nontrivial bridge obstruction. The exact algebra gate in:

```text
outputs/jensen_hankel_bridge_algebra.md
work/rh_compute/results/jensen_hankel_bridge_algebra.json
```

gives a positive rational sequence that passes finite reshaped-Hankel signs for
`k = 2,3`, `N = 3`, but has negative degree-3 Jensen discriminant. Therefore
low-order finite sign checks cannot be promoted into Jensen hyperbolicity.

## Proof Obligations

1. Prove the all-order shifted reshaped-Hankel sign-consistency condition for
   the actual `A_k(lambda)` sequence, or replace it with a stronger condition
   that is provably satisfied by the Xi/Phi coefficients.
2. Prove a sign-regularity theorem that converts that all-order structure into
   hyperbolicity of every Jensen polynomial, not just degree 2.
3. Equivalently, prove that every binomially weighted Jensen window
   `B^{d,n,0}_j = binom(d,j) A_{n+j}(0)` is PF-infinity, or prove this from
   the signed-Hankel/sign-consistency structure.
4. Identify the exact theorem ecosystem: signed total positivity,
   sign-regular kernels, compound matrices, variation-diminishing transforms,
   multiplier sequences, or a new Xi-specific determinant identity.
5. Prove the limiting passage from all Jensen polynomials to the required
   Laguerre-Polya or Newman conclusion without assuming RH.
6. Preserve the heat-flow parameter discipline: finite positive-lambda
   certificates can guide the theorem, but they cannot replace a lambda-zero
   proof or a uniform interval theorem.

## Kill Gates

Reject a proposed proof if it:

```text
uses only finitely many A_k(lambda);
uses only the unshifted n=0 reshaped-Hankel block;
treats degree-2 Jensen as representative of all degrees;
assumes Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0;
uses asymptotic fixed-degree Jensen results as all-shift Jensen results;
depends only on local zero repulsion;
confuses ordinary coefficient PF with signed-Hankel sign consistency.
```

## Current Status

This theorem target is open. The finite Arb certificate and the exact
low-degree algebra gate make the target sharper, but they do not prove the
bridge.
