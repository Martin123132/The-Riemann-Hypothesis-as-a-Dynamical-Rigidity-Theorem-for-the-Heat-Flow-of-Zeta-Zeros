# Exact Weber completion of the modular continuation

Date: 2026-08-13

Status: exact Abel-regularized completion and cancellation; not a proof of a
bound for the surviving common Gaussian tails

Fix any matched continuation label

```text
alpha=A-2k,
n_B=L+k, n_A=k, k>=1.                                (WC1)
```

For endpoint `D` and `n=(D-alpha)/2`, put `p=sqrt(x-z)`.  Including the
principal modular factor `exp(i*pi/4)`, its dual endpoint kernel becomes

```text
K_D=exp(i*pi/4) 2sqrt(x) exp(i*pi*x(D^2/4+n^2))
 int_0^sqrt(x) [2+i*pi*D^2(x-p^2)]
 exp(-i*pi[D^2*p^2/4+n^2*x^2/p^2]) dp.              (WC2)
```

The exact monotone coordinate

```text
u_D=Dp/2-nx/p                                         (WC3)
```

satisfies

```text
x(D^2/4+n^2)-D^2*p^2/4-n^2*x^2/p^2
   =alpha^2*x/4-u_D^2,
u_D(0+)=-infinity,
u_D(sqrt(x))=alpha*sqrt(x)/2.                         (WC4)
```

Thus both endpoints have exactly the same Gaussian phase and interval.  No
stationary-phase approximation has been used.

Now complete the `p` integral in (WC2) to infinity in the Abel sense.  For
`Re(a),Re(b)>0`,

```text
I0(a,b)=int_0^infinity exp(-a*p^2-b/p^2)dp
       =sqrt(pi)/(2sqrt(a)) exp(-2sqrt(a*b)),
I2=-partial_a I0=I0[1/(2a)+sqrt(b/a)].                (WC5)
```

Taking the principal boundary values
`a -> i*pi*D^2/4`, `b -> i*pi*n^2*x^2` gives the exact completed kernel

```text
K_D^[0,infinity]
 =2*i*pi*alpha*x^(3/2)*exp(i*pi*alpha^2*x/4).         (WC6)
```

The right side is independent of `D`.  Hence the B and A completed kernels
cancel exactly for every fixed label in (WC1), including negative `alpha`.
Returning to the original finite interval leaves only

```text
K_B-K_A=-exp(i*pi/4+i*pi*alpha^2*x/4)
 int_(alpha*sqrt(x)/2)^infinity
       [J_B(u)-J_A(u)] exp(-i*pi*u^2) du,             (WC7)

J_D=4sqrt(x)p_D(u)^2[2+i*pi*D^2(x-p_D(u)^2)]
       /[D*p_D(u)^2+2nx].                             (WC8)
```

Equation (WC7) is cancellation preserving: one phase, one interval, and one
amplitude difference.  It also explains why the scalar term in the leading
stationary coefficient was not the final continuation residue; the exact
completion includes every stationary correction and cancels the completed
bulk for each fixed pair.

There is a necessary cutoff guard.  With the same dual cutoff `N>L` on both
endpoint sums, exact finite reindexing gives

```text
sum_(n=0)^N (B_n-A_n)
 =sum_(n=0)^L B_n-A_0
  +sum_(k=1)^(N-L)(B_(L+k)-A_k)
  -sum_(k=N-L+1)^N A_k.                              (WC9)
```

The last shift-defect block has exactly `L=2481422` terms.  It cannot be dropped
merely because every fixed high-index term tends to zero.  Vanishing of this
block in the original cutoff or Abel prescription is an explicit remaining
obligation before (WC7) may be summed through `k=infinity`.

The exceptional `k=0` term is outside this certificate because the A term
has `n_A=0` and is a face contribution, not an interior Weber saddle.  For
`k>=79789` the lower tail endpoint in (WC7) is negative; those terms are
already included in the common tail series and are not a separate bulk.

Pi provenance: every `pi` in (WC2)-(WC8) comes from the equation-(9)
Fourier/Kummer phase and the principal Jacobi modular factor.  Formula (WC5)
is first used in its absolutely convergent half-plane and only then taken to
the stated Abel boundary.

Proof boundary: the exact completion and completed-pair cancellation for
every fixed `k>=1`, together with the finite-cutoff identity (WC9), are
proved.  Infinite-cutoff reindexing, vanishing of the shift defect, a bound
for the resulting common-tail series, the exceptional `k=0` face/fold
splice, the complete source-minus-target estimate, `T_upper`, `Lambda<=0`,
PF-infinity, RH, and any prize-level conclusion remain open.
