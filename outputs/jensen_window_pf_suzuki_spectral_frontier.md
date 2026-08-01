# Jensen-Window PF Suzuki Spectral Frontier

Date: 2026-07-23

Status: exact truncation-path spectral reformulation, finite-ceiling
route guards, one open signed quadratic-form target, and an
internally audited cofinal terminal reduction. This is not
a proof of global Fredholm nonvanishing, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_suzuki_spectral_frontier.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_spectral_frontier.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_spectral_frontier.py
```

## Fixed-Space Truncation

Start from Suzuki's supported Hankel operator

```text
(K[t]f)(x)=1_(x<t)*integral_(-infinity)^t K(x+y)f(y)dy, with K(s)=0 for s<0
```

and translate `L2(-infinity,t)` to one fixed copy of
`L2(0,infinity)`:

```text
For (U_t f)(u)=f(t-u), u>0: A_t=U_t K[t] U_t^(-1), (A_t g)(u)=integral_0^infinity K(2t-u-v)g(v)dv
```

K[t] vanishes on inputs supported in (-infinity,-t) and its output vanishes on (-infinity,-t); it is therefore represented by the continuous kernel K(x+y) on the finite square (-t,t)^2.

The translated kernel is supported on the triangle `u+v<=2t`.
Consequently,

```text
||K[t]||_HS^2=integral_0^(2t) (2t-s)*|K(s)|^2 ds
```

For the check kernel `K(s)=s for s>=0` at `t=1`,
both sides equal `4/3`.

## Hilbert-Schmidt Guard

If the kernel is not identically zero, choose `R` with
`c_R=integral_0^R |K(s)|^2 ds>0`. Then

```text
If c_R=integral_0^R |K(s)|^2 ds>0, then for t>R/2: ||K[t]||_HS^2 >= (2t-R)*c_R -> infinity
```

Thus the generic estimate `||K[t]||<=||K[t]||_HS<1` can prove
only a bounded initial interval. This is a limitation of the
Hilbert-Schmidt certificate, not evidence that the operator norm
itself exceeds one.

## Exact Spectral Reformulation

Suzuki's continuity and support assumptions imply

```text
For continuous K with K(s)=0 on s<0, t -> A_t is continuous in Hilbert-Schmidt norm and A_0=0
For 0<=t_1<=t_2, ||K[t_1]||<=||K[t_2]||, because K[t_1] is the quadratic-form compression of the self-adjoint K[t_2]
```

Since a compact self-adjoint operator attains its norm at an
eigenvalue, continuity from the zero operator gives

```text
For every T>=0: det(I+/-K[t])!=0 for all 0<=t<=T iff ||K[T]||<1 iff ||K[t]||<1 for all 0<=t<=T
```

Equivalently,

```text
If T_*=inf{t>=0: ||K[t]||=1}<infinity, compact self-adjointness forces +1 or -1 to be an eigenvalue of K[T_*], so det(I-K[T_*])=0 or det(I+K[T_*])=0
```

The all-times qualifier cannot be dropped:

```text
At one isolated time determinant nonvanishing is weaker: A=2P for a rank-one orthogonal projection has ||A||=2, det(I-A)=-1, and det(I+A)=3
```

So Suzuki condition (3) is neither an unspecified determinant
problem nor a request for one t-uniform margin. It is exactly the
problem of preventing a finite first crossing of either spectral
edge.

## High-Contour Ceiling

Suzuki's unconditional local proof supplies the sufficient test

```text
Suzuki's Proposition 4.2 eigenfunction contradiction is certified when M_v^2*exp(4v*t)<1 for some v>1/2+omega, where M_v=sup_(u in R)|Theta_(omega,nu)(u+i*v)|
```

For zeta, Stirling's formula gives the lower witness

```text
For zeta and a=nu*omega: Theta_(omega,nu)(i*v)=[xi(v+1/2-omega)/xi(v+1/2+omega)]^nu ~(2*pi/v)^a as v->infinity
```

while Suzuki's upper estimate has the same polynomial scale.
Therefore

```text
T_HC(omega,nu)=sup_(v>1/2+omega) [log(1/M_v)]_+/(2v)<infinity; the certificate M_v^2*exp(4v*t)<1 cannot prove any t>T_HC
```

The exponential support loss eventually defeats every polynomial
gain in contour height. Sharpening only the polynomial constant or
power cannot turn Proposition 4.2's mechanism into an all-`t`
proof for a fixed pair.

## Correct Global Target

Under the desired Hermite-Biehler hypothesis, Suzuki proves

```text
If E_(omega,nu) belongs to HB, Suzuki Proposition 5.1 gives a full L2 isometry K, K[t]=P_t K P_t, and ||K[t]||<1 for every finite t
```

but this conditional argument also shows why a uniform margin is
the wrong target:

```text
Under the same HB hypothesis, P_t->I strongly and K is an isometry, so lim_(t->infinity)||K[t]||=1; no epsilon>0 can satisfy ||K[t]||<=1-epsilon for every t
```

The surviving noncircular obligation is therefore

```text
The noncircular live target is, for every finite t and every nonzero f: |<K[t]f,f>|<||f||^2, equivalently I+K[t]>0 and I-K[t]>0
```

A useful proof would factor or otherwise control both signed forms
without identifying the high-contour kernel with an inner real-axis
multiplier. For one fixed pair this spectral task does not
directly prove the terminal limit. The cofinal logic is sharper:

```text
For one fixed pair, pointwise-in-t strict contractivity does not directly prove J(t;z,z)->0. On a strictly decreasing cofinal omega_n->0 sequence, however, all-time contractivity gives a bounded causal multiplier, forces RH by shifted-zero accumulation, and then Suzuki's Theorem 2.3 supplies the terminal limit
```

The causal-multiplier and shifted-zero proof is recorded in
`outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`.

## Source

- Masatoshi Suzuki, `Hamiltonians arising from L-functions in the Selberg class`, Propositions 4.1-4.3 and 5.1, Theorem 2.4: https://arxiv.org/abs/1606.05726
- `outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md`
- `outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`

## Proof Boundary

This artifact proves the fixed-space translation, exact Hilbert-Schmidt identity, path continuity and monotonicity, the all-history determinant/contractivity equivalence, the first-crossing alternative, and finite-ceiling guards for Hilbert-Schmidt and Suzuki sup-contour certificates. It records Suzuki's conditional HB contraction theorem and proves that its norms tend to one. Together with the separately checked causal-multiplier reduction, it records that the terminal premise is redundant for a cofinal determinant family. It does not prove the signed quadratic-form target, global Fredholm nonvanishing for zeta without HB, the Phi Pick sign, PF-infinity, RH, or Lambda<=0.
