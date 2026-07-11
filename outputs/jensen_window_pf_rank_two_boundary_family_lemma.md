# Jensen-Window PF Rank-Two Boundary-Family Lemma

Date: 2026-07-10

Status: exact all-degree boundary-family lemma with mixture countermodel. This is not a proof
of a zeta representation, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_rank_two_boundary_family_lemma`.

```text
work/rh_compute/results/jensen_window_pf_rank_two_boundary_family_lemma.json
python work/rh_compute/scripts/jensen_window_pf_rank_two_boundary_family_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_rank_two_boundary_family_lemma.py
```

Current result:

```text
validated Jensen-window PF rank-two boundary-family lemma: 11 rows, 0 issues, 4 exact identities, 1 all-degree factorization, 1 integer-product closure, 2 exact countermodels, 1 open structural handoff
```

## Exact Family

For `1/2<u<1`, put

```text
a=(2*u-1)/u, b=1-u, c=2*u-1,
A_k(u)=a^(k-1)*(c+k*b)/u.
```

Then `A_0=A_1=1`, and exact binomial algebra gives, for every `d>=1,n>=0`,

```text
J_(d,n)(z)=A_n*(1+a*z)^(d-1)
             *(1+a*(c+n*b+d*b)/(c+n*b)*z).
```

Both root locations are strictly negative. Thus every shifted Jensen window in this
model is hyperbolic. Its contractions and defects are

```text
x_k=1-b^2/(c+k*b)^2,
d_k=1/(k+c/b)^2
   = integral_0^1 t^(k+c/b-1)*(-log t) dt.
```

The cubic critical value `x=5/9,y=21/25,z=45/49` at `u=3/5` is one
member of this all-degree hyperbolic boundary family.

## Mixture Gate

The family is not convex in the required hyperbolicity sense. The equal positive
mixture of `u=3/5` and `u=2/3` gives

```text
J_3(z)=1+3*z+(47/24)*z^2+(41/108)*z^3,
Disc(J_3)=-937/3456<0.
```

Therefore a positive integral or mixture over these boundary sequences is not, by
itself, the missing zeta kernel representation.

## Multiplier Product Route

With `u=(alpha+1)/(alpha+2)`, the family becomes

```text
M_k^(alpha)=(alpha/(alpha+1))^(k-1)*(alpha+k)/(alpha+1),
sum_(k>=0) M_k^(alpha)*z^k/k!
  =(1+z/(alpha+1))*exp(alpha*z/(alpha+1)).
```

This is a classical multiplier sequence. Finite pointwise products remain multiplier
sequences because the corresponding diagonal real-root preservers compose. Their
ratio contractions satisfy

```text
-log x_k=sum_j -log(1-1/(k+alpha_j)^2).
```

The atoms must carry genuine integer multiplicity. Fractional weights are unsafe:
after removing positive geometric scaling, `alpha=1` and exponent `1/2` gives

```text
J_3(z)=1+3*sqrt(2)*z+3*sqrt(3)*z^2+2*z^3,
Disc=-432*sqrt(2)-324*sqrt(3)+378+324*sqrt(6)
    <-27/125<0.
```

Thus a general positive measure representation of `-log x_k` is insufficient; the
live target is a discrete counting-measure factorization or another structure with
an independently proved stability-preserving composition law.

## Handoff

A one-parameter positive multiplier sequence has an exact all-degree/all-shift Jensen factorization with one simple and one repeated negative root. Finite pointwise products give an integer-atomic canonical-product cone. However, both an equal convex mixture and a fractional pointwise power have negative cubic discriminants, so arbitrary positive superposition or measure weights are not the missing kernel theorem.

```text
outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md
outputs/jensen_window_pf_bridge_target.md
```
