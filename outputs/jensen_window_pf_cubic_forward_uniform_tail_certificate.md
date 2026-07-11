# Jensen-Window PF Cubic Forward-Uniform Tail Certificate

Date: 2026-07-10

Status: exact all-shift degree-3 propagation through lambda=0 with the
higher-degree bridge open. This is not a proof of PF-infinity, RH, or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_cubic_forward_uniform_tail_certificate.json
python work/rh_compute/scripts/jensen_window_pf_cubic_forward_uniform_tail_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_cubic_forward_uniform_tail_certificate.py
```

Current result:

```text
validated Jensen-window PF cubic forward-uniform tail certificate: 10 rows, 0 issues, 3 exact flow identities, 1 weighted source cap, 1 initial weighted tail, 1 forward-uniform tail, 1 full cubic propagation theorem, 1 lambda=0 cubic theorem, 1 open higher-degree handoff
```

## Exact Increment Flow

Set `d_k=1-x_k`, `q_k=d_k^(-1/2)`, and
`g_k=q_(k+1)-q_k`. The defect equation gives

```text
q_k'=r_k*((2*k-1)*q_k-(2*k+3)*(q_k^3-q_k)/q_(k+1)^2)
```

The derivative with respect to the next increment is strictly positive:

```text
2*(2*k + 5)*(g + q - 1)**2*(g + q + 1)**2/((g + q)*(g + j + q)**3)
```

At `g_(k+1)=g_k=g`, exact factorization and deletion of manifestly
negative terms give, for `0<=g<=1` and `q_k>=sqrt(k)`,

```text
sqrt(k)*g_k'/r_k<=6+13/sqrt(k)+5/k+1/k^(3/2)<7
```

The last strict inequality uses `sqrt(319)>17` and the exact remainder
`3843/4913<1`.

## Entry Tail

The lambda=-100 theorem already proves

```text
g_k<=(100000/99999)^2*(5*k+6)^(3/2)/k^2, k>=319
```

After multiplying by `sqrt(k)`, the right side decreases with `k` and
is strictly below 12 at `k=319`. The squared integer margin is

```text
57062151807283727967356093296
```

## Forward Uniformity

For a finite `L`, put `R_L=sup_[-100,L] r_0(lambda)<infinity`. On an
interval where the cubic cone holds, the componentwise barrier

```text
G_k(lambda)=(12+7*R_L*(lambda+100))/sqrt(k)
```

is a supersolution. Use the explicit coercive regularization

```text
G_k^epsilon=G_k+epsilon*k*exp(5*R_L*(lambda+100)).
```

Because `g_k<=1`, any first contact is at a finite index. At contact,
`H_(k+1)-H_k=H_k/k`, while next-increment monotonicity gives the
extra coupling cap `2*R_L*(2*k+5)*H_k/k<5*R_L*H_k=H_k'`.
The weighted source cap excludes contact. Letting `epsilon` vanish proves

```text
sup_(lambda in [-100,L]) g_k(lambda)
 <=(12+7*R_L*(L+100))/sqrt(k)->0.
```

This tail makes the reciprocal-defect first-crossing principle legitimate.
The all-k entry theorem and exact inward frontier condition therefore
propagate `0<=g_k<=1` to every finite `lambda>=-100`.

## Consequence And Boundary

Every shifted degree-3 Jensen polynomial is hyperbolic at `lambda=0`.
No degree-4 invariant or all-degree PF/Jensen theorem is proved here; those
are the remaining obligations before any RH conclusion could be considered.
