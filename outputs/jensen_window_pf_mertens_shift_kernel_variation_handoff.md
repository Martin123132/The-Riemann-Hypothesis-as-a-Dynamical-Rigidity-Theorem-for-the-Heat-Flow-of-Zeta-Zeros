# Jensen-Window PF Mertens Shift-Kernel Variation Handoff

Date: 2026-07-24

Status: exact shift decomposition, two-axis bounded variation, and
joint-cancellation handoff with one open arithmetic gate.
This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_shift_kernel_variation_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_shift_kernel_variation_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_shift_kernel_variation_handoff.py
```

## Exact Shift Decomposition

Retain the notation of Corollary 11.22Z.3:

```text
N=2K+1,
R=ceil(K^(1-alpha/2)),
P=P_(K,R),
G=G_(alpha,K,R),
rho_K:=pi^2*K/N^2<pi^2/(4K).               (SKV.1)
```

For `1<=h<=3K-2`, put

```text
I_(K,h):=[K+1,4K-1-h] intersect Z,
W_(K,h)(n):=G(n,n+h),

Q_(K,h):=sum_(n in I_(K,h))
 mu(n)mu(n+h)W_(K,h)(n).                    (SKV.2)
```

Every pair `K<n<m<4K` occurs exactly once with `h=m-n`, so

```text
O_(alpha,K)=sum_(h=1)^(3K-2)Q_(K,h).        (SKV.3)
```

This identity preserves the sign across shifts.

## Three Consecutive Blocks

The shift support is the disjoint consecutive union

```text
I_h^(cc)=[K+1,2K-1-h],

I_h^(cf)=[
 max(K+1,2K-h),
 min(2K-1,4K-1-h)
],

I_h^(ff)=[2K,4K-1-h],                       (SKV.4)
```

with an interval omitted when empty. With `i=n-K`, `m=n+h`, and
`p_i=P(i,K)`,

```text
W_h(n)=
  P(i,i+h),                    n in I_h^(cc),
  p_i*b_K(m),                  n in I_h^(cf),
  P(K,K)b_K(n)b_K(n+h),        n in I_h^(ff). (SKV.5)
```

The future weight

```text
b_K(x)=(2K/x)^(1+alpha)*(4K-x)/(2K)
```

is strictly decreasing on `[2K,4K)`. PSD Cauchy-Schwarz and the
diagonal estimate from Corollary 11.22Z.3 give

```text
0<W_h(n)<=rho_K<pi^2/(4K).                   (SKV.6)
```

## Fixed-Shift Variation

On the current-current block, the sine-integral formula becomes

```text
P(i,i+h)
 =2/N integral_(h*pi/N)^((2i+h)*pi/N)S_R(u)du. (SKV.7)
```

Since `S_R>0`, this is strictly increasing in `i`, and its variation
is less than `rho_K`. On the current-future block, `p_i` increases
while `b_K(K+i+h)` decreases. For positive monotone factors,

```text
Var(p_i*b_K(K+i+h))<2rho_K.                  (SKV.8)
```

The future-future block is decreasing and has variation less than
`rho_K`. Each of the at most two block joins costs at most `rho_K`.
Consequently,

```text
Var_(n in I_(K,h))W_(K,h)(n)<6rho_K.         (SKV.9)
```

Define the fixed-origin maximal correlation

```text
C_(K,h)(x):=sum_(n=K+1)^x mu(n)mu(n+h),

A_(K,h):=max_(K+1<=x<=4K-1-h)|C_(K,h)(x)|.  (SKV.10)
```

Finite Abel summation gives, with `b_h=4K-1-h`,

```text
Q_(K,h)
 =C_(K,h)(b_h)W_h(b_h)
  +sum_(n=K+1)^(b_h-1)
    C_(K,h)(n)[W_h(n)-W_h(n+1)],             (SKV.11)

|Q_(K,h)|
 <7rho_K*A_(K,h)
 <7pi^2/(4K)*A_(K,h).                        (SKV.12)
```

Thus each shift weight is a clean bounded-variation test of a local
two-point Mobius correlation.

## Fixed-Base Variation

Fix a current base point `n=K+i`. As `h` grows, the second position
first traces

```text
P(i,i+1),...,P(i,K),
```

then traces `p_i*b_K(n+h)` down the future interval. Successive
sine-integral boundary strips are disjoint, and

```text
integral_0^pi S_R(u)du
 =2sum_(r<=R)1/(2r-1)^2<pi^2/4.
```

Therefore

```text
Var_(i<j<=K)P(i,j)<pi^2/N.                  (SKV.13)
```

The future tail costs at most `p_i<=rho_K`. A future base point has a
decreasing `h`-row. Uniformly,

```text
Var_h W_(K,h)(n)<pi^2/N+rho_K.               (SKV.14)
```

Put

```text
B_(K,n)(H):=mu(n)sum_(h=1)^H mu(n+h).
```

A second Abel summation gives the exact fixed-base expansion and the
bound

```text
|sum_(h=1)^(4K-1-n)
  mu(n)mu(n+h)W_(K,h)(n)|

 <pi^2/K*
  max_(H<=4K-1-n)|sum_(h<=H)mu(n+h)|.        (SKV.15)
```

Thus the other coordinate is an `O(1/K)` bounded-variation test of
local Mertens intervals.

## Mass Calibration

The future weights satisfy

```text
B_1:=sum_(2K<=n<4K)b_K(n)<=(2K+1)/2,

B_2:=sum_(2K<=n<4K)b_K(n)^2
 <=(2K+1)(4K+1)/(12K).                      (SKV.16)
```

The total current-current, current-future, and future-future masses
are respectively less than

```text
pi^2*K/8,  pi^2*K/4,  pi^2*K/8.
```

Hence

```text
sum_(K<n<m<4K)G(n,m)<pi^2*K/2,

|O_(alpha,K)|<pi^2*K/2.                     (SKV.17)
```

For a single shift, PSD Cauchy-Schwarz, (SKV.16), and the trace bound
give

```text
sum_(n in I_(K,h))W_(K,h)(n)<13pi^2/24.     (SKV.18)
```

The absolute normalized scale remains `K^(1-alpha)`, one full power
above the summable target.

## Why Axiswise Cancellation Is Insufficient

Suppose, conditionally, every maximal fixed-shift correlation in
(SKV.10), or every local Mertens maximum in (SKV.15), had the
square-root-plus-epsilon size `O(K^(1/2+epsilon))`. Summing the other
coordinate absolutely still gives

```text
O_(alpha,K)=O(K^(1/2+epsilon)).              (SKV.19)
```

After multiplication by `K^(-alpha)`, this closes only
`alpha>1/2+epsilon`, not a cofinal sequence tending to zero.

Matomaki, Radziwill, and Tao prove an averaged Chowla theorem whose
natural terminal two-point scale is `o(K^2)` when both averaging
ranges have size `K`. Their theorem does not provide the varying
maximal endpoints in (SKV.10). Even granting that promotion, (SKV.12)
would give only `o(K)`, with no fixed power saving sufficient for the
cofinal criterion.

Similarly, scalar all-interval logarithmic cancellation applied one
base point at a time yields at best `K*log^(-A)K` after the `1/K`
normalization. Logarithms do not replace the missing power.

The stronger condition

```text
sum_(h<=3K-2)A_(K,h)
 =O_epsilon(K^(1+epsilon))                   (SKV.20)
```

would suffice, but it is stronger than random square-root behavior
after absolute summation and is not proposed as the likely route.

## Joint Gate

The exact remaining criterion is

```text
sup_J sum_(K dyadic<=2^J)
 K^(-alpha)sum_(h=1)^(3K-2)Q_(K,h)
 <infinity                                   (SKV.21)
```

for every member of one fixed cofinal sequence `alpha_j->0`.
Equivalently, before any axiswise absolute value, the scale is

```text
<u_I,G u_I>+<u_II,G u_II>-2<u_I,G u_II>.    (SKV.22)
```

The new bounded-variation formulas are theorem interfaces for a
genuinely joint estimate in base point, shift, Vaughan type, and
possibly dyadic scale. They do not themselves provide that estimate.

The weighted joint Mobius gain, full Burnol bound, RH, PF-infinity,
and `Lambda<=0` remain open.

Primary sources:

- Matomaki, Radziwill, and Tao, *An Averaged Form of Chowla's
  Conjecture*:
  https://doi.org/10.2140/ant.2015.9.2167
- Matomaki, Shao, Tao, and Teravainen, *Higher Uniformity of
  Arithmetic Functions in Short Intervals I. All Intervals*:
  https://doi.org/10.1017/fmp.2023.28
- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
