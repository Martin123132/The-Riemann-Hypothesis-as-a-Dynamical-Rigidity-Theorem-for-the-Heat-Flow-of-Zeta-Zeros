# Jensen-Window PF Heat-Flow Infinite Cone-Invariance Certificate

Date: 2026-07-10

Status: exact infinite ratio-cone propagation theorem and open all-order handoff. This is not a proof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_infinite_cone_invariance_certificate`.

```text
work/rh_compute/results/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.json
python work/rh_compute/scripts/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.py
```

## Exact Defect System

Put `d_k=1-x_k` and `h_k=d_k-d_(k+1)=x_(k+1)-x_k`. The exact
pointwise walls give `0<=d_k<=1/k`, `|h_k|<=1/k`, and `r_k<=r_0`.
Direct algebra gives

```text
d_k'=2*r_k*((2*k+3)*(1-d_k)*d_(k+1)-(2*k-1)*d_k)
h_k'=a_k*(h_(k+1)-h_k)+q_k*h_k+c_k
2*r_k*(2*k+5)*(1-d_(k+1))^2>=0
2*r_k*d_k^2*(6-(2*k+5)*d_k)>=0
```

The adjacent equation can be written

```text
h_k'=a_k*(h_(k+1)-h_k)+q_k*h_k+c_k,
a_k>=0, c_k>=0, q_k<=176*R,
R=sup r_0(lambda)<infinity on each compact heat interval.
```

## Infinite Maximum Principle

Set `z_k(t)=exp(-176*R*t)*h_k(-100+t)`. Since `|h_k|<=1/k`,
`z_k->0` uniformly on compact time intervals. Therefore every negative
spatial infimum is attained at a finite index. Define
`m(t)=min(0,inf_(k>=1) z_k(t))`. Whenever `m(t)<0`, uniform tail decay
leaves a locally finite active minimum set. At each active index,

```text
a_k*(z_(k+1)-z_k)>=0,
(q_k-176*R)*z_k>=0,
exp(-176*R*t)*c_k>=0.
```

Thus `D_+m(t)=min{z_k'(t):z_k(t)=m(t)}>=0` whenever `m(t)<0`. On a
connected component `(alpha,beta)` of `{t:m(t)<0}`, this makes `m`
nondecreasing, contradicting `m(alpha)=0`. Hence a negative component cannot form.
Iterating over compact heat intervals
proves `h_k(lambda)>=0` for every `k>=1` and finite `lambda>=-100`.

## Consequence

Full cone entry at `lambda=-100`, the exact pointwise walls, and this
maximum principle prove the full infinite ratio cone at `lambda=0`.
The remaining problem is the noncircular all-shape Jensen-window or
sign-regular bridge; ratio-cone propagation alone is not PF-infinity.

```text
outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md
outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
