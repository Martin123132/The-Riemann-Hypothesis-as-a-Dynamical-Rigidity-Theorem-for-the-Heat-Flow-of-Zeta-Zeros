# Jensen-Window PF Negative-Lambda -100 Compound Order-Three Entry Certificate

Date: 2026-07-12

Status: all-shift contiguous order-three entry theorem at lambda=-100
with open forward propagation. This is not a proof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda -100 compound order-three entry certificate: 11 rows, 0 issues, 322 positive coefficients, 318 prefix compound margins, 318 prefix defect gaps, 1 exact tail theorem, 1 all-shift entry theorem, 1 open forward handoff
```

## Repaired Prefix

The full repaired source is supplemented by a ten-row precision repair
around `k=210` and a one-row collar at `k=321`. Direct 1024-bit Arb
arithmetic proves

```text
C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0, 0<=n<=317,
d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)>0, 0<=n<=317.
```

The weakest compound margin is `7.70728029027026774090631029652944448849006795880897978491067E-1` at `n=0`.
The weakest defect-gap enclosure is `7.77208729941522315092729901005318452539804324300304408879936E-9` at `n=317`.
All 319 prefix reciprocal-defect increments are positive and below one;
the 318 adjacent increment differences are positive, and the 317
adjacent compound-margin differences are positive. These shape facts are
finite diagnostics, not assumptions in the tail proof.

The splice anchor is

```text
s_319=(639/2)*d_319=[0.502756829007979122036006833538041640 +/- 3.30E-37]
s_319-251/500>7.56829007979122036006833538041639670971857146232348179021168E-4.
```

## Exact Tail

Put `k=n+2`, `m=2*k+1`, and

```text
u_k=q_k-q_(k-1),
v_k=q_(k+1)-q_k.
```

The completed adaptive-defect theorem gives decreasing `d_k` and
increasing `s_k=(2*k+1)*d_k/2`. Hence for every `k>=320`,

```text
s_j>251/500 for j>=319,
q_j^2<250*(2*j+1)/251.
```

Scaled-defect growth also gives

```text
u_k<=q_k*(1-sqrt((2*k-1)/(2*k+1))).
```

For `m>=641`, exact squaring proves

```text
1-sqrt(1-2/m)<1/m+1/m^2
```

because the cleared condition is `m^2-2*m-1>0`. The same estimate at
the next index yields

```text
q_k*u_k<(250/251)*(1+1/m), u_k*v_k<(250/251)*(1/m)*(1+1/m)^2
```

Since

```text
C_(k-2)=1-u_k*v_k+q_k*(v_k-u_k)
        >=1-u_k*v_k-q_k*u_k,
```

we obtain

```text
C_(k-2)>1-(250/251)*(1+2/m+2/m^2+1/m^3)
C_n>57613471/66107054971>0, n>=318.
```

The endpoint comparison is exact; after setting `m=641` the rational
lower margin is approximately `0.0008715177377856874`.

## Theorem And Boundary

Combining the Arb prefix with the analytic tail proves

```text
C_n(-100)>0,
D_(3,n)(-100)=det[A_(n+i+j)]_(i,j=0..2)<0,
for every integer n>=0.
```

This closes all-shift contiguous order-three entry at one heat parameter.
Forward propagation to `lambda=0`, noncontiguous order-three minors, and
every higher compound order remain open.

```text
outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md
outputs/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
