# Jensen-Window PF Mertens Anchor/Logarithmic Tail-Energy Reduction

Date: 2026-07-23

Status: exact all-cutoff anchor/logarithmic stable-prefix reduction with
one open RH-equivalent arithmetic gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.json
python work/rh_compute/scripts/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py
```

## Adjacent Tails Recover The Anchor

Fix `0<alpha<1` and write

```text
a_n:=mu(n)n^(-alpha),
A_alpha(N):=sum_(n<=N)a_n,

r_N:=sum_(n>N)a_n/n,
r_0=1/zeta(1+alpha).                                  (MALTER.1)
```

Here `A_alpha` is the weighted Mobius prefix. It is not the reciprocal
partial sum called `A_(omega,N)` in the Burnol tail-discrepancy note.

Adjacent tails satisfy

```text
r_(N-1)-r_N=a_N/N.                                    (MALTER.2)
```

Multiplying by `N` and telescoping gives the exact reconstruction

```text
A_alpha(N)
 =sum_(j=0)^(N-1)r_j-N*r_N.                           (MALTER.3)
```

Thus the constant anchor that is invisible in one post-prefix q-block is
not independent once the compatible tail sequence at every cutoff is
retained. The recovery is global in `N`; one tail value still does not
determine the anchor.

The inverse identity is the previously checked Abel transform

```text
r_N
 =-A_alpha(N)/(N+1)
  +sum_(n>N)A_alpha(n)/(n(n+1)).                       (MALTER.4)
```

No zero-free region is needed for either identity.

## Weighted Hardy/Copson Equivalence

Set

```text
E_alpha:=sum_(N>=1)N^(alpha-2)|A_alpha(N)|^2,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2,
b_N:=A_alpha(N)/N,
s_N:=r_(N-1).
```

Then (MALTER.3) says

```text
b_N=H(s)(N)-r_N,
H(s)(N):=N^(-1)sum_(j=1)^N s_j.                       (MALTER.5)
```

On `l2(N^alpha)`, the weighted Hardy and Copson operators admit the
nonsharp Schur bounds

```text
||H||<=h_alpha:=sqrt(2(3-alpha))/(1-alpha),
||C||<=c_alpha:=sqrt(2(3+alpha))/(1+alpha),            (MALTER.6)

(Cb)(N):=sum_(n>=N)b_n/n.
```

The one-step weight shifts and (MALTER.4)-(MALTER.6) give

```text
R_alpha^(1/2)
 <=(1+c_alpha)E_alpha^(1/2),

E_alpha^(1/2)
 <=(1+2^(alpha/2)h_alpha)
   (|r_0|^2+R_alpha)^(1/2).                            (MALTER.7)
```

Since `r_0` is finite for `alpha>0`,

```text
E_alpha<infinity
 iff
R_alpha<infinity.                                     (MALTER.8)
```

This is an exact finiteness equivalence with explicit nonsharp operator
constants. It does not estimate either arithmetic sequence.

## Stable-Prefix Energy Across All Cutoffs

Put `alpha=2omega`. The stable-prefix term in the Burnol
tail-discrepancy split is

```text
P_(omega,N)
 =r_N^2*sum_(k<=N)k^alpha.                             (MALTER.9)
```

The elementary power-sum comparison yields

```text
(1/(1+alpha))*N^alpha*r_N^2
 <=P_(omega,N)/N
 <=N^alpha*r_N^2.                                     (MALTER.10)
```

Therefore

```text
E_alpha<infinity
 iff
R_alpha<infinity
 iff
sum_(N>=1)P_(alpha/2,N)/N<infinity.                    (MALTER.11)
```

This is the precise all-cutoff form of the weighted-prefix anchor. A
uniform bound on each `P_(omega,N)` is one logarithm weaker than the
required series.

Composing (MALTER.11) with the checked weighted-prefix/Mertens criterion
gives, for every fixed cofinal sequence `alpha_j->0`,

```text
RH
 iff
sum_(N>=1)P_(alpha_j/2,N)/N<infinity
for every j.                                           (MALTER.12)
```

This is an RH-equivalent target, not a proof of convergence.

## The Exact Logarithmic Barrier

The pointwise stable-prefix bound gives only

```text
r_N=O(N^(-(1+alpha)/2)).
```

Inserting this envelope into (MALTER.3) gives at best

```text
A_alpha(N)=O(N^((1-alpha)/2)),
```

whose contribution to `E_alpha` is of order `1/N`. The loss is genuinely
logarithmic, as the following scalar model shows.

Let

```text
gamma:=(1+alpha)/2,
r_N:=(N+1)^(-gamma),
a_N:=N[r_(N-1)-r_N].                                  (MALTER.13)
```

Then

```text
r_N=sum_(n>N)a_n/n,
0<a_N<=gamma*N^(-gamma)<=gamma*N^(-alpha).             (MALTER.14)
```

The associated stable-prefix energy obeys

```text
2^(-(1+alpha))/(1+alpha)<=P_N<=1,                      (MALTER.15)
```

so `sum_N P_N/N` diverges. Moreover,

```text
A(N)
 ~[gamma/(1-gamma)]N^(1-gamma),

N^(alpha-2)A(N)^2
 ~[gamma/(1-gamma)]^2/N.                               (MALTER.16)
```

This is a scalar/operator countermodel, not a surrogate for `mu`. It
proves that coefficient size plus the critical pointwise tail rate cannot
close the logarithmic series by generic analysis. It leaves open a
Mobius-specific mean-square cancellation theorem.

## Relation To The Earlier Local Guard

The earlier fixed-cutoff constant-tail countermodel remains correct. At
one cutoff the q-block observation is invariant under a constant prefix
shift. Equations (MALTER.2)-(MALTER.3) add the missing compatibility
between adjacent cutoffs. In short:

```text
one cutoff:       anchor is invisible;
all cutoff tails: anchor is reconstructed;
critical sup P:   reconstruction is one logarithm nonsummable;
sum P_N/N:        exact weighted-prefix/Mertens energy.
```

The obstruction concerns the stable-prefix component `P` alone. A
cofinal uniform theorem for the full Burnol energy `Q` would already
imply RH through the existing three-gate criterion; (MALTER.12) is a
cross-coordinate sharpening, not an extra condition after full `Q`.

## Literature Guard

The primary context checked was:

- Burnol's analytic estimate and strengthened Nyman-Beurling argument:
  https://arxiv.org/abs/math/0202166
- Baez-Duarte's arithmetical Nyman-Beurling reformulation:
  https://arxiv.org/abs/math/0011254
- weighted Hardy/Copson context:
  https://arxiv.org/abs/2508.00388
- weighted Dirichlet-space zero-free criteria:
  https://doi.org/10.1007/s11785-025-01661-2
- the current Hardy-space approximation formulation:
  https://arxiv.org/abs/2606.16097

A focused search found nearby Mobius approximation and zero-free
criteria, but no exact displayed match for (MALTER.11). That is a search
report, not a novelty or priority claim.

## Open Gate

The sharpened arithmetic target is:

```text
for every alpha in one fixed cofinal sequence alpha_j->0,

sum_(N>=1) P_(alpha/2,N)/N < infinity                  (MALTER.17)
```

without RH or a zero-free hypothesis. Equivalently, prove the
weighted-prefix/Mertens energy. The present reduction supplies no Mobius
cancellation and proves neither RH, PF-infinity, nor `Lambda <= 0`.
