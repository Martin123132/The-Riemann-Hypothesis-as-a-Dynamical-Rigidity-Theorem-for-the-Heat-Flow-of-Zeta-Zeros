# Jensen-Window PF Mertens Centered-Bridge Vaughan Handoff

Date: 2026-07-23

Status: exact trace/off-diagonal and centered Brownian-bridge Vaughan
reduction with one open two-component gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_centered_bridge_vaughan_handoff.py
```

## One Global Coefficient Kernel

Fix `0<alpha<1` and put

```text
q_n:=mu(n)n^(-1-alpha),
r_N:=sum_(n>N)q_n,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2.                    (MCBVH.1)
```

For dyadic `K`, let

```text
y_j:=r_(K+j-1),                       1<=j<=K,
c_(K,r):=<y,phi_r>,                   0<=r<K.
```

Define coefficient kernels by

```text
h_(K,0)(n)
 :=K^(-1/2)min((n-K)_+,K),                            (MCBVH.2)

h_(K,r)(K+m)
 :=sqrt(2/K)
   sin(pi*r*m/K)/[2sin(pi*r/(2K))],
                                      1<=m<K, 1<=r<K, (MCBVH.3)
```

and set the nonconstant kernels to zero elsewhere. Counting how many
tails contain each `q_n`, together with Lemma 11.22V, gives

```text
c_(K,r)
 =sum_(n>=1)mu(n)n^(-1-alpha)h_(K,r)(n),
                                      0<=r<K.          (MCBVH.4)
```

Thus the dyadic cosine energy

```text
E_alpha
 :=sum_(K dyadic)K^alpha
   sum_(r=0)^(K-1)|c_(K,r)|^2                         (MCBVH.5)
```

is the squared norm of one analysis operator `A_alpha` applied to the
sequence `mu(n)`. The dyadic weight comparison from Lemma 11.22V gives

```text
R_alpha<infinity iff E_alpha<infinity.                (MCBVH.6)
```

Let

```text
F_(K,r)(n):=n^(-1-alpha)h_(K,r)(n).
```

The mean row splits at `2K`:

```text
||F_(K,0)||_2^2
 <=[1+2^(-1-2alpha)/(1+2alpha)]K^(-2alpha).           (MCBVH.7)
```

For `r>0`, finite sine orthogonality gives

```text
sum_(m=1)^(K-1)|h_(K,r)(K+m)|^2
 =1/[4sin^2(pi*r/(2K))].
```

Since `sin(pi*r/(2K))>=r/K`,

```text
||F_(K,r)||_2^2
 <=K^(-2alpha)/(4r^2).                                (MCBVH.8)
```

Consequently

```text
sum_(K dyadic)K^alpha
 sum_(r=0)^(K-1)||F_(K,r)||_2^2

 <=[1+2^(-1-2alpha)/(1+2alpha)+pi^2/24]
   sum_(K dyadic)K^(-alpha)
 <infinity.                                           (MCBVH.9)
```

This proves that `A_alpha` is Hilbert-Schmidt from `l2(N)` to the
dyadic coefficient space. It does not close the problem because
`mu` is not in `l2`.

For a finite dyadic cutoff `J`, define

```text
G_(alpha,J)(m,n)
 :=sum_(K dyadic<=2^J)K^alpha
   sum_(r=0)^(K-1)F_(K,r)(m)F_(K,r)(n).
```

Every row entering a fixed `J` is in `l1`, so expansion before passage
to the limit is legitimate:

```text
E_alpha(J)
 =D_alpha(J)+2O_alpha(J),                              (MCBVH.10)

D_alpha(J)
 :=sum_n mu(n)^2 G_(alpha,J)(n,n),

O_alpha(J)
 :=sum_(m<n)mu(m)mu(n)G_(alpha,J)(m,n).
```

The diagonal converges by (MCBVH.9), while `E_alpha(J)` is monotone and
nonnegative. Hence

```text
E_alpha<infinity
 iff
sup_J O_alpha(J)<infinity.                             (MCBVH.11)
```

The RH-strength content is therefore purely signed and off-diagonal.
There is no generic `l_infinity -> l2` operator theorem available:
the scalar construction in Lemma 11.22W has

```text
x*_n:=q*_n n^(1+alpha),      |x*_n|<=1,
```

but its first-mode output energy diverges. That sequence is not
Mobius; it is a proof-method guard.

## Centered Abel Quotient

For `0<=m<K`, define the ordinary and weighted local paths

```text
S_(K,m):=M(K+m)-M(K),

T_(K,m):=sum_(ell=1)^m
          mu(K+ell)(K+ell)^(-1-alpha).                 (MCBVH.12)
```

Both start at zero. On K-vectors use the quotient norm modulo constants

```text
||[x]||_K
 :=min_a||x-a1||_2,

||[x]||_K^2
 =sum_(m=0)^(K-1)
   |x_m-K^(-1)sum_j x_j|^2.                            (MCBVH.13)
```

With `a_m=(K+m)^(-1-alpha)`, extend the Abel map to every K-vector:

```text
(A_Kx)_0:=a_1x_0,

(A_Kx)_m
 :=a_1x_0+sum_(j=1)^m a_j(x_j-x_(j-1)),
                                      1<=m<K.          (MCBVH.14)
```

For the anchored Mertens path, `A_K S_K=T_K`. More importantly,

```text
A_K(c1)=a_1c1.                                        (MCBVH.15)
```

So `A_K` descends to an invertible map on the quotient by constants.
The diagonal-plus-prefix formula gives

```text
||A_Kx||_2
 <=(2+alpha)K^(-1-alpha)||x||_2.                       (MCBVH.16)
```

Writing `b_m=1/a_m`, the inverse is

```text
(A_K^(-1)y)_0=b_1y_0,

(A_K^(-1)y)_m
 =b_my_m+sum_(j=1)^(m-1)(b_j-b_(j+1))y_j,             (MCBVH.17)
```

and

```text
||A_K^(-1)y||_2
 <=2^alpha(3+alpha)K^(1+alpha)||y||_2.                 (MCBVH.18)
```

Because constants map to constants, the same bounds descend to the
quotient:

```text
[2^alpha(3+alpha)]^(-1)
 K^(-1-alpha)||[S_K]||_K

 <=||[T_K]||_K

 <=(2+alpha)K^(-1-alpha)||[S_K]||_K.                  (MCBVH.19)
```

This resolves the commutator concern from Lemma 11.22W.
No commutator estimate is needed: the Abel map itself is uniformly invertible after
passing to the correct quotient.

The centered DCT identity is

```text
||[T_K]||_K^2
 =sum_(r=1)^(K-1)|tau_(K,r)|^2.                        (MCBVH.20)
```

All modes `r>=ceil(sqrt(K))` are already dyadically summable. Therefore
(MCBVH.19)-(MCBVH.20) give

```text
sum_(K dyadic)K^alpha
 sum_(1<=r<ceil(sqrt(K)))|tau_(K,r)|^2
 <infinity

iff

L_alpha
 :=sum_(K dyadic)K^(-2-alpha)||[S_K]||_K^2
 <infinity.                                            (MCBVH.21)
```

Retain the affine mean series

```text
M_alpha
 :=sum_(K dyadic)K^alpha
   |sqrt(K)r_K-tau_(K,0)|^2.                           (MCBVH.22)
```

The exact two-component criterion is now

```text
R_alpha<infinity
 iff
M_alpha+L_alpha<infinity.                              (MCBVH.23)
```

For every member of any fixed cofinal sequence `alpha_j->0`,
(MCBVH.23) is equivalent to RH. Neither term has been proved finite.

## Discrete Brownian-Bridge Kernel

The finite pairwise-variance identity gives

```text
||[S_K]||_K^2
 =K^(-1)sum_(0<=i<j<K)
   |S_(K,j)-S_(K,i)|^2.                                (MCBVH.24)
```

Writing `x_i:=mu(K+i)`, `1<=i<K`, this is exactly

```text
||[S_K]||_K^2
 =sum_(i,j=1)^(K-1)x_i x_j B_K(i,j),                   (MCBVH.25)

B_K(i,j)
 :=min(i,j)-ij/K
  =min(i,j)(K-max(i,j))/K.
```

`B_K` is the discrete Brownian-bridge Green matrix. Separating its
diagonal and writing `j=i+h` gives

```text
||[S_K]||_K^2=Delta_K+2C_K,                            (MCBVH.26)

Delta_K
 :=sum_(i=1)^(K-1)
   mu(K+i)^2 i(K-i)/K,

C_K
 :=sum_(h=1)^(K-2)
   sum_(i=1)^(K-1-h)
   mu(K+i)mu(K+i+h)i(K-i-h)/K.                         (MCBVH.27)
```

The diagonal is harmless:

```text
0<=Delta_K<=(K^2-1)/6,

sum_(K dyadic)K^(-2-alpha)Delta_K<infinity.            (MCBVH.28)
```

As the bridge energies are nonnegative, the complete nonconstant gate
is therefore

```text
L_alpha<infinity
 iff
sup_J sum_(K dyadic<=2^J)K^(-2-alpha)C_K
 <infinity.                                            (MCBVH.29)
```

This is an anchored but base-averaged, triangularly weighted two-point
Mobius correlation.

## Pre-Collapse Vaughan Handoff

For `1<=h<=K-2`, put `i=n-K` and define

```text
g_(K,h)(n)
 :=mu(n+h)i(K-i-h)/K
   1_(1<=i<=K-1-h).                                    (MCBVH.30)
```

Then

```text
|g_(K,h)(n)|<=K/4,

C_K
 =sum_(h=1)^(K-2)
   sum_(K<n<=2K)mu(n)g_(K,h)(n).                       (MCBVH.31)
```

Choose `1<=U_K,V_K<=K` and define the Vaughan coefficients

```text
a_(K,d)
 :=sum_(bc=d,b<=U_K,c<=V_K)mu(b)mu(c),

b_(K,d)
 :=sum_(c|d,c>V_K)mu(c).                               (MCBVH.32)
```

Green and Tao's finite Vaughan identity gives

```text
TI_K
 :=sum_h sum_(d<=U_KV_K)a_(K,d)
   sum_(K/d<w<=2K/d)g_(K,h)(d*w),

TII_K
 :=sum_h sum_(V_K<d<=2K/U_K)b_(K,d)
   sum_(max(U_K,K/d)<w<=2K/d)
   mu(w)g_(K,h)(d*w),

C_K=-TI_K+TII_K.                                      (MCBVH.33)
```

The exact live Type I/II target is

```text
sup_J sum_(K dyadic<=2^J)
 K^(-2-alpha)(-TI_K+TII_K)
 <infinity.                                            (MCBVH.34)
```

The signed difference must be retained. Replacing it with
`|TI_K|+|TII_K|` is not an equivalent handoff.

The elementary absolute scale is still one full power too large:

```text
|C_K|
 <=sum_(h,i)i(K-i-h)/K
 <=K^3/8.                                             (MCBVH.35)
```

After multiplication by `K^(-2-alpha)`, this gives the divergent scale
`K^(1-alpha)`. The averaged Chowla theorem of Matomaki, Radziwill, and
Tao gives the natural unweighted average `o(K^2)` after collapsing the
redundant common-shift average. Multiplication by the Brownian-bridge
weight `O(K)` still gives only `o(K^3)`, not the
`O_epsilon(K^(2+epsilon))` scale required here. This is a power audit,
not a lower bound on the actual correlation.

## Open Gate

The quotient reduction removes only local constants. It does not control
the affine mean series `M_alpha`. Thus the surviving cofinal gate has
two noninterchangeable pieces:

```text
1. Prove M_alpha<infinity.

2. Prove the signed Brownian-bridge/Vaughan criterion (MCBVH.34).
```

A bound

```text
||[S_K]||_K^2=O_epsilon(K^(2+epsilon))
```

with `epsilon<alpha` is a sufficient calibration for the second piece,
but is not proved here. The exact reduction concerns `R_alpha`, not the
full Burnol energy `Q`. No Type I/II power gain, affine mean estimate,
full-Q bound, RH, PF-infinity, or `Lambda <= 0` follows from this handoff;
all remain open.

Primary sources:

- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1 for the finite Vaughan identity:
  https://doi.org/10.5802/aif.2401
- Matomaki, Radziwill, and Tao, *An Averaged Form of Chowla's
  Conjecture*:
  https://doi.org/10.2140/ant.2015.9.2167
