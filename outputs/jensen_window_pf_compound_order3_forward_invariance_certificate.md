# Jensen-Window PF Compound Order-Three Forward-Invariance Certificate

Date: 2026-07-12

Status: exact all-shift contiguous order-three propagation through
lambda=0 with noncontiguous and higher compounds open. This is not a
proof of PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order3_forward_invariance_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order3_forward_invariance_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order3_forward_invariance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order3_forward_invariance_certificate.py
```

Current result:

```text
validated Jensen-window PF compound order-three forward-invariance certificate: 10 rows, 0 issues, 2 exact identities, 1 cooperative flow, 1 inward boundary theorem, 1 coefficient-growth lemma, 1 infinite maximum principle, 1 full forward propagation theorem, 1 lambda=0 theorem, 1 open higher-compound handoff
```

## Cooperative Flow

Put `k=n+2` and

```text
a=q_(k-1), b=q_k, c=q_(k+1), d=q_(k+2),
C_n=a*c-b^2+1, C_(n+1)=b*d-c^2+1.
```

The exact reciprocal-defect heat equation is

```text
q_j'=r_j*((2*j-1)*q_j-(2*j+3)*(q_j^3-q_j)/q_(j+1)^2)
```

and direct differentiation factors as

```text
C_n'/r_k=alpha_k*C_(n+1)+beta_k*C_n
alpha_k=(2*k+5)*a*(b*d+c^2-1)/(c*d^2)
beta_k=-N_k/(c^2*(b^2-1))
N_k=(2*k+1)*a^2*c^2+(2*k+1)*a*b^2*c-(2*k+1)*a*c+(4*k+6)*b^4+(-4*k+2)*b^2*c^2-(4*k+6)*b^2
```

Inside the strict ratio cone, `alpha_k>0`. In particular,

```text
at C_n=0, C_n'/r_k=alpha_k*C_(n+1).
```

Thus the local boundary is exactly one-sided and inward whenever the next
compound margin is nonnegative.

## Compact Tail Control

The proved cubic forward theorem gives, on every finite heat interval
`[-100,L]`, a constant `B_L` such that

```text
0<=g_j=q_(j+1)-q_j<=B_L/sqrt(j).
```

The ratio cone gives `q_j>=sqrt(j)`, while summing the increment bound
gives `q_j=O_L(sqrt(j))`. Use the exact scaling

```text
k=h^(-2), b=y/h, a=(y-u*h^2)/h, c=(y+v*h^2)/h, d=(y+(v+w)*h^2)/h
```

where `u,v,w` are the scaled neighboring increments and remain in a
compact box. Exact cancellation gives

```text
lim_(h->0) h^2*alpha_k=4,
lim_(h->0) (alpha_k+beta_k)=2*(u+2*v-3*w)/y.
```

The canceled denominator is nonzero on the compact box because `y>=1`
and `h^2<=1/2`. Therefore

```text
alpha_k=O_L(k), alpha_k+beta_k=O_L(1).
```

## Infinite Maximum Principle

The same tail estimates give `C_n=O_L(1)`. Set

```text
z_n=C_n/(n+1).
```

Then `z_n->0` uniformly on compact heat intervals and

```text
z_n'=r_k*alpha_k*((n+2)/(n+1))*z_(n+1)+r_k*beta_k*z_n.
```

Since `r_k<=r_0`, the effective diagonal coefficient

```text
r_k*(alpha_k+beta_k+alpha_k/(n+1))
```

is bounded above on each compact interval. After one exponential
integrating factor, every negative spatial infimum is attained at a finite
index and has nonnegative upper-right derivative. The standard connected-
component argument excludes a negative component.

The lambda=-100 entry theorem is strict. Variation of constants in the
cooperative equation then preserves strict positivity at every fixed index.
Hence

```text
C_n(lambda)>0 for every n>=0 and finite lambda>=-100,
D_(3,n)(0)<0 for every n>=0.
```

This closes the complete shifted contiguous order-three layer at lambda=0.
Noncontiguous order-three minors and every higher compound order remain
open, so the all-order signed-Hankel/Jensen bridge is not proved.

```text
outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md
outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md
outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
