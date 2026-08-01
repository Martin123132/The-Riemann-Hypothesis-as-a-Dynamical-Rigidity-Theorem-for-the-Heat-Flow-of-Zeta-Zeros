# Jensen-Window PF Burnol Tail-Discrepancy/Dyadic Reduction

Date: 2026-07-23

Status: exact finite-tail and dyadic square-function reduction with one open
arithmetic gate. This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.json
python work/rh_compute/scripts/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py
```

## Finite Versus Limiting Discrepancy

Fix `0<omega<1/2`, set `alpha=2omega`, and write

```text
a_omega(d)=mu(d)d^(-alpha),

A_(omega,N)=sum_(d<=N)a_omega(d)/d,

A_(omega,infinity)=1/zeta(1+alpha),

r_(omega,N)=A_(omega,infinity)-A_(omega,N).         (BTDR.1)
```

For integer `k>=1`, define

```text
B_(omega,N)(k)
 =sum_(d<=N)a_omega(d)floor(k/d),

S_(omega,N)(k)=A_(omega,N)k-B_(omega,N)(k).
```

The limiting divisor sum is finite at every `k`:

```text
B_(omega,infinity)(k)
 =sum_(d<=k)a_omega(d)floor(k/d),

S_(omega,infinity)(k)
 =A_(omega,infinity)k-B_(omega,infinity)(k).
```

Put

```text
H_(omega,N)(k)
 =B_(omega,infinity)(k)-B_(omega,N)(k)
 =sum_(N<d<=k)a_omega(d)floor(k/d).                 (BTDR.2)
```

Direct subtraction gives the exact finite-tail discrepancy

```text
D_(omega,N)(k)
 =S_(omega,N)(k)-S_(omega,infinity)(k)
 =H_(omega,N)(k)-k*r_(omega,N).                    (BTDR.3)
```

Equivalently, the omitted fractional-part tail is

```text
T_(omega,N)(k)
 =sum_(d>N)a_omega(d){k/d}
 =k*r_(omega,N)-H_(omega,N)(k)
 =-D_(omega,N)(k).                                 (BTDR.4)
```

The sum in (BTDR.4) is finite up to `d=k`; beyond `k`,
`{k/d}=k/d`, so the remainder is absolutely convergent.

There is also an exact summation-by-parts form. Let

```text
g_k(d)=d*{k/d},
r_(omega,d)=sum_(m>d)a_omega(m)/m.
```

Since `g_k(d)=k` for `d>k`,

```text
T_(omega,N)(k)
 =r_(omega,N)g_k(N+1)
  +sum_(d=N+1)^k r_(omega,d)
   [g_k(d+1)-g_k(d)].                              (BTDR.5)
```

Thus the remaining operator is driven by jumps of the divisor sawtooth,
not by a generic bounded Hardy operator.

## Exact Three-Gate Split

Define the discrete energies

```text
Q_(omega,N)
 =sum_(k>=1)k^(alpha-2)S_(omega,N)(k)^2,

Q_(omega,infinity)
 =sum_(k>=1)k^(alpha-2)S_(omega,infinity)(k)^2,

R_(omega,N)
 =sum_(k>=1)k^(alpha-2)D_(omega,N)(k)^2.            (BTDR.6)
```

For each fixed `k`, `S_(omega,N)(k)` tends to
`S_(omega,infinity)(k)`. Fatou and the Hilbert triangle inequality give

```text
sup_N Q_(omega,N)<infinity

iff

Q_(omega,infinity)<infinity
and
sup_N R_(omega,N)<infinity.                         (BTDR.7)
```

On the stable prefix `k<=N`, `H_(omega,N)(k)=0`. Hence

```text
P_(omega,N)
 :=sum_(k<=N)k^(alpha-2)D_(omega,N)(k)^2
 =r_(omega,N)^2 sum_(k<=N)k^alpha.                 (BTDR.8)
```

The elementary bounds

```text
N^(1+alpha)/(1+alpha)
 <=sum_(k<=N)k^alpha
 <=N^(1+alpha)
```

show the exact scalar equivalence

```text
sup_N P_(omega,N)<infinity

iff

r_(omega,N)=O(N^(-(1+alpha)/2))
            =O(N^(-1/2-omega)).                    (BTDR.9)
```

The post-prefix energy is

```text
U_(omega,N)
 =sum_(k>N)k^(alpha-2)
  |H_(omega,N)(k)-k*r_(omega,N)|^2,

R_(omega,N)=P_(omega,N)+U_(omega,N).                (BTDR.10)
```

Combining (BTDR.7)-(BTDR.10) yields the exact three-gate criterion

```text
sup_N Q_(omega,N)<infinity

iff

  Q_(omega,infinity)<infinity,
  r_(omega,N)=O(N^(-1/2-omega)),
  sup_N U_(omega,N)<infinity.                       (BTDR.11)
```

The first line is the limiting Jordan-error energy gate. The second is the
forced reciprocal-zeta tail rate. The third is the genuinely finite-tail
square-function gate. None is discarded or hidden inside a generic norm
statement.

## Multiplicative-Interval Form

Let

```text
M_alpha(x)=sum_(d<=floor(x))a_omega(d).
```

Counting `floor(k/d)` by multiples gives

```text
H_(omega,N)(k)
 =sum_(m<=floor(k/(N+1)))
  [M_alpha(floor(k/m))-M_alpha(N)].                 (BTDR.12)
```

Alternatively, partitioning by `q=floor(k/d)` gives

```text
H_(omega,N)(k)
 =sum_(q<=floor(k/(N+1))) q*
  [M_alpha(floor(k/q))
   -M_alpha(max(N,floor(k/(q+1))))].                (BTDR.13)
```

Each bracket in (BTDR.13) is a weighted Mobius increment on the
multiplicative interval

```text
(k/(q+1), k/q],
```

intersected with `(N,infinity)`.

## Dyadic Square Function

For `K>=N`, define

```text
V_(omega,N)(K)
 =sum_(K<k<=2K)
  |H_(omega,N)(k)-k*r_(omega,N)|^2,

K_j=2^j*N.
```

Since `alpha-2<0`, dyadic weight comparison gives

```text
2^(alpha-2) sum_(j>=0) K_j^(alpha-2)V_(omega,N)(K_j)
 <=U_(omega,N)
 <=sum_(j>=0) K_j^(alpha-2)V_(omega,N)(K_j).        (BTDR.14)
```

Therefore the remaining post-prefix condition is exactly a uniform
summable dyadic square function. For example, the sufficient envelope

```text
V_(omega,N)(2^j*N)
 <=C*(2^j*N)^(2-alpha)*2^(-eta*j)                  (BTDR.15)
```

for some `eta>0`, uniformly in `N` and `j`, would close this component.
The critical estimate without the summable factor `2^(-eta*j)` is not
enough by itself.

## Literature Boundary

Burnol records the weighted infinite function and RH-conditional convergence
of its finite natural approximants:

```text
https://arxiv.org/abs/math/0202166
```

Báez-Duarte proves analogous scalar lower-bound obstructions and divergence
results for unweighted natural approximants:

```text
https://arxiv.org/abs/math/0011254
```

The multiplicative fractional-part autocorrelation gives an exact Gram-kernel
framework:

```text
https://arxiv.org/abs/math/0306251
```

Those results motivate (BTDR.5) and the dyadic route, but they do not state or
prove the weighted uniform estimate (BTDR.14). In particular, an unweighted
natural-approximation theorem cannot be transferred by deleting the
`d^(-2omega)` and `k^(2omega-2)` weights.

## Open Gate

For every shift in one explicit cofinal sequence, prove all three lines of
(BTDR.11). The new finite-tail obligation is

```text
sup_N sum_(j>=0)(2^j*N)^(2omega-2)
 V_(omega,N)(2^j*N)<infinity.                      (BTDR.16)
```

Equation (BTDR.16) is open. The exact decomposition does not establish its
hypotheses, RH, PF-infinity, or `Lambda <= 0`.
