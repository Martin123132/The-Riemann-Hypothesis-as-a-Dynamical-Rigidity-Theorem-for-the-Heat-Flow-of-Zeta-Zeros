# Finite geometric RSI contour residual and pole release

Date: 2026-08-27

Status: exact finite identity certified; not a proof of the continuation defect

The A21 infinite interchange is unnecessary for every finite odd-label
prefix.  For any integer `N>=0` and any non-pole `z`, finite geometric algebra
gives

```text
csc(pi z)
 =-2i sum_(k=0)^(N-1) exp(i*pi*(2k+1)*z)
  +exp(2*pi*i*N*z)csc(pi*z).                         (FR1)
```

This is an identity of meromorphic functions, not an asymptotic expansion.
Insert (FR1) directly into the exact Riemann-Siegel contour before any
infinite summation.  If `P_N(s)` is the resulting finite odd-exponential
prefix and

```text
E_N(s)=1/(2i) integral_C
       exp(-i*pi*z^2+2*pi*i*N*z) z^(-s)csc(pi*z) dz, (FR2)
```

then exactly

```text
RSI(s)=P_N(s)+E_N(s).                                (FR3)
```

Translate `u=z-N`.  Since

```text
exp(i*pi*N^2)=(-1)^N,
sin(pi(u+N))=(-1)^N sin(pi*u),                       (FR4)
```

the signs cancel and

```text
E_N(s)=1/(2i) integral_(C-N)
       exp(-i*pi*u^2)(u+N)^(-s)csc(pi*u) du.         (FR5)
```

Move `C-N` back to `C`.  The branch point `u=-N` remains to the left of the
deformation strip.  Exactly the poles `m=1-N,...,0` are crossed, and

```text
Res_(u=m) [exp(-i*pi*u^2)(u+N)^(-s)csc(pi*u)]
 =(m+N)^(-s)/pi.                                     (FR6)
```

Therefore

```text
E_N(s)=sum_(n=1)^N n^(-s)+C_N(s),

C_N(s)=1/(2i) integral_C
       exp(-i*pi*u^2)(u+N)^(-s)csc(pi*u) du.         (FR7)
```

The altered complex-height numerical contour check has maximum discrepancy
`1.7844932832759947e-56` for
`N=1,2,3,5`; it confirms the residue orientation independently of the finite
geometric algebra.

For the actual odd window, put

```text
n_-=(A-1)/2=79788,
n_+=(B+1)/2=2561211.        (FR8)
```

There are `2481423` labels, and subtraction of two
copies of (FR3) gives the exact kernel-level common-contour identity

```text
P_[A,B]=P_(n_+)-P_(n_-)
       =E_(n_-)-E_(n_+)
       =-sum_(n=n_-+1)^(n_+) n^(-s)
        +C_(n_-)-C_(n_+).                            (FR9)
```

Thus the finite window is already an exact difference of two translated
residual contours plus an explicit pole block.  No claim that the A21
infinite remainder vanishes is needed.  The next gate must align the finite
exponential-prefix normalization and branches with the corrected physical
`Q_K`, then combine the released pole block with `G` and `A_transition`
before taking norms.

Pi provenance: every `pi` in (FR1)--(FR7) comes from the original
`sin(pi z)` denominator, Gaussian `exp(-i*pi*z^2)`, and integer contour
translation.  No geometric circle construction or fitted constant is used.

Proof boundary: exact finite meromorphic identity, contour translation,
residue inventory, actual-window index arithmetic, and a surrogate numerical
orientation check only.  The finite exponential integrals have not yet been
matched branch-for-branch to the corrected physical `Q_K`; no `C_K`, `D_K`,
`Delta_KU`, `J_Z`, or non-A enclosure, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
