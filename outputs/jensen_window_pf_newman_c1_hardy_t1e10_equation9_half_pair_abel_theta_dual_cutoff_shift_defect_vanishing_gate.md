# Vanishing of the modular dual cutoff-shift defect

Date: 2026-08-13

Status: fixed-height weighted-L1 tail certificate; not a proof of a
quantitative bound for the surviving common Gaussian-tail series

The exact finite-cutoff reindexing leaves the A-endpoint block

```text
R_N=sum_(n=N-L+1)^N A_n,   L=2481422.                    (SD1)
```

To control it without losing the `x=0` geometry, scale `z=xq`.  The A dual
kernel is exactly

```text
K_(A,n)(x)=e^(i*pi/4)x int_0^1 (1-q)^(-1/2)
 [2+i*pi*A^2*x*q] exp(i*phi_(A,n)(q,x)) dq,

phi_(A,n)=pi*x[A^2*q/4-n^2*q/(1-q)].                 (SD2)
```

The prefactor `x` in (SD2) is essential.  A direct absolute bound gives

```text
|K_(A,n)(x)| <= 4x+2*pi*A^2*x^2.                    (SD3)
```

For `n>=A`, integrate once in `q`.  Since

```text
|phi_q|=pi*x[n^2/(1-q)^2-A^2/4],                    (SD4)
```

the `q=1` boundary vanishes and an explicit total-variation estimate yields

```text
|K_(A,n)(x)|
 <= [368/(63*pi)+(668/315)A^2*x]/n^2.               (SD5)
```

Let `a_t(x)=x^(-7/4)(1-x)^(-1/4)` and split at `x_0=n^(-2)`.
Use (SD3) below the split and (SD5) above it.  Since
`(1-x)^(-1/4)<=2^(1/4)` on the half interval,

```text
M_n=int_0^(1/2)|a_t(x)K_(A,n)(x)|dx
 <= C_half*n^(-1/2)
    +4*(668/315)A^2*n^(-2)
    +2^(1/4)*(8*pi*A^2/5)n^(-5/2),                  (SD6)

C_half=[21.975491371568892328131411901384911237176208934419417604341223702055270048795308 +/- 2.95e-79].
```

For `N-L+1>=A`, every term of (SD1) satisfies (SD6), hence

```text
||R_N||_weighted-L1
 <= L times the right side of (SD6) at n=N-L+1
 -> 0 as N->infinity.                               (SD7)
```

This proves that the fixed-width shift defect vanishes and licenses the
matched continuation reindexing in the weighted integral at this height.
The dual `+/-n` multiplicity only multiplies (SD7) by two.

Pi provenance: `pi` in (SD2)-(SD6) is the equation-(9)/Jacobi phase
normalization.  No fitted constant is introduced.

Proof boundary: (SD7) is fixed-height because `A` and `L` are held fixed.
It proves cutoff-shift vanishing, not a useful finite-`N` source-error bound,
not a height-uniform reindexing theorem, and not a bound for the common
Gaussian-tail series, exceptional `k=0` term, complete source-minus-target
residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or a prize-level result.
