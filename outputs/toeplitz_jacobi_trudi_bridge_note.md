# Toeplitz Jacobi-Trudi Bridge Note

Date: 2026-07-04

Status: exact reindexing note for theorem search. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.

## Purpose

The active coefficient-PF route studies the upper-triangular Toeplitz minors:

```text
det(c_{q_j-r_i})_{i,j=1}^m
```

for:

```text
c_k(lambda) = A_k(lambda) / k! = mu_{2k}(lambda) / (2k)!.
```

The finite Arb certificates test these determinants directly. The theorem-search question is whether every structurally nonzero determinant has a positive all-order explanation.

## Exact Reindexing

For strictly increasing row and column sets:

```text
R = (r_1 < ... < r_m)
Q = (q_1 < ... < q_m)
```

the minor is structurally nonzero exactly when:

```text
q_i >= r_i for every i.
```

When this holds, choose:

```text
K = max(0, max_i(r_i - i), max_i(q_i - i))
lambda_i = K + i - r_i
mu_i     = K + i - q_i
```

with `i = 1..m`. Then `lambda` and `mu` are partitions and `lambda_i >= mu_i`. Moreover:

```text
lambda_i - mu_j - i + j = q_j - r_i.
```

Thus, after normalizing:

```text
d_k = c_k / c_0
```

the Toeplitz minor is:

```text
det(c_{q_j-r_i}) = c_0^m det(d_{q_j-r_i})
                 = c_0^m s_{lambda/mu}(d),
```

where the last determinant is the Jacobi-Trudi skew-Schur determinant under the specialization `h_k -> d_k`.

The factor `c_0^m` is harmless for signs as long as `c_0 > 0`; the actual positivity burden is Schur positivity of the specialization `h_k -> d_k`.

## Workbench Script

The exact map is implemented in:

```text
work/rh_compute/scripts/toeplitz_jacobi_trudi_map.py
work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py
```

Example:

```text
python work/rh_compute/scripts/toeplitz_jacobi_trudi_map.py --rows 0,1 --cols 1,2
```

returns the skew shape:

```text
lambda = (1, 1)
mu = (0, 0)
Jacobi-Trudi indices = [[1, 2], [0, 1]]
```

corresponding to:

```text
det [[d_1, d_2],
     [d_0, d_1]].
```

Executable validation:

```text
python work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py
```

Current result:

```text
OK Toeplitz/Jacobi-Trudi example maps: one nonzero skew shape and one structural zero
validated Toeplitz/Jacobi-Trudi reindexing: N=10, orders<=5, 124129 minors, 39094 nonzero maps, 85035 structural zeros
```

This is an exact algebraic reindexing gate. It checks the nonzero example
above, a structural-zero example, and the full `N = 10`, orders `<= 5`
combinatorial grid against fixed counts. It does not assert positivity of the
resulting skew-Schur evaluations.

## Why This Helps

This turns the missing bridge into a precise Schur-positivity problem:

```text
prove s_{lambda/mu}(d_0,d_1,d_2,...) >= 0
for every skew shape arising from the Toeplitz support condition.
```

In fact, the arising shapes are exactly the finite Toeplitz minor shapes under the chosen indexing. So a positive specialization, planar network, production matrix, or measure on partitions with complete homogeneous coordinates `d_k` would imply the full coefficient PF-infinity target.

## Source Cross-Check

The ASW/Edrei theorem family is the final coefficient-PF bridge:

```text
PF-infinity of d_k
  => restricted Laguerre-Polya generating function
```

but using that theorem directly as a hypothesis would be circular. The noncircular target is stronger structural data, such as:

```text
1. d_k is realized as complete homogeneous symmetric functions h_k(alpha)
   for a nonnegative specialization alpha;
2. a planar-network or Lindstrom-Gessel-Viennot model has path weights whose
   path matrix gives exactly d_{q-r};
3. a production-matrix or continued-fraction model verifies the same Schur
   positivity without assuming real zeros of F_0.
```

The useful literature anchors are:

```text
Aissen-Schoenberg-Whitney / Edrei generating-function theorem:
https://www.pnas.org/doi/10.1073/pnas.37.5.303
https://link.springer.com/article/10.1007/BF02786970

Khare, "Polya frequency sequences: analysis meets algebra":
https://www.imsc.res.in/~amri/algcomb/201112_khare.pdf

Toeplitz minors and skew-Schur specializations:
https://arxiv.org/abs/1706.02574
```

## Boundary

This note does not prove such a positive specialization exists. It only removes ambiguity in the algebraic translation.

Finite certificates remain evidence and falsification pressure. A referee-grade proof still needs an all-order theorem proving Schur positivity for the zeta heat-flow coefficient specialization, or a different bridge of equal strength.
