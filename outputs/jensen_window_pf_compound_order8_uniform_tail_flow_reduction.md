# Jensen-Window PF Compound Order-Eight Uniform Tail And Flow Reduction

Date: 2026-07-13

Status: exact uniform eventual signed order-eight tail and conditional
forward reduction with one open `lambda=-100` entry. This is not a
proof of all-shift order eight, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_uniform_tail_flow_reduction.py
```

## Universal Signed Tail

The graded-kernel Vandermonde lemma applies at `m=8`, where
`D=binom(8,2)=28` and `epsilon_8=1`. It gives

```text
det K=33664847019245568000*G_2^28*h^28+o(h^28)
Q_8=positive_scale*(33664847019245568000*G_2^28*h^28+o(h^28))
there exists N_8 such that Q_(8,n)(lambda)>0 for every n>=N_8 and -100<=lambda<=0
```

The threshold is finite but non-effective. This is an eventual tail,
not an all-shift order-eight theorem.

## Endpoint Coordinate

Signed Desnanot-Jacobi gives

```text
Q_(8,n)*Q_(6,n+2)=Q_(7,n+1)^2-Q_(7,n)*Q_(7,n+2)
Q_(8,n)>0 iff Q_(7,n+1)^2>Q_(7,n)*Q_(7,n+2)
```

The completed positive `Q_6=-H_6` and `Q_7=-H_7` layers make this
an exact strict log-concavity target for the endpoint `Q_7` sequence.

## Cooperative Heat Flow

The general affine-Hankel and adjacent Plucker identities specialize to

```text
Q_(8,n)'=(4*n+58)*delta(Q_(8,n))
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n=(4*n+58)*Q_(7,n)/Q_(7,n+1)>0, b_n=((4*n+58)/(4*n+54))*(log Q_(7,n+1))'
```

Since `Q_(7,n)>0` on the complete heat interval, `a_n>0`. The uniform
eventual tail and variation of constants therefore prove

```text
[Q_(8,n)(-100)>0 for all n] => [Q_(8,n)(lambda)>0 for all n and -100<=lambda<=0]
completed signed contiguous layers through order eight imply every arbitrary-column signed layer through order eight
```

conditional only on all-shift endpoint entry.

## Countermodel Gate

The exact rational sequence

```text
(1,1,1/2,1/6,1/24,1/120,1/720,1/5040,1/40320,1/362880,1/3628800,1/39916800,1/479001600,1/6227020800,1/87120000000)
```

has every available signed contiguous minor through order seven strict
and positive, but

```text
Q_(7,1)^2-Q_(7,0)Q_(7,2)=-463/10246317406459624992735028091299392073078206541205704645017600000000000000000000000<0,
Q_(8,0)=-463/69210459277322870541707704182767616000000000000000<0.
```

Thus no lower-cone promotion can replace the new endpoint theorem.

## Remaining Endpoint Target

The sole fixed-order-eight input still missing is

```text
Q_(8,n)(-100)>0 for every n>=0, equivalently strict log-concavity of Q_(7,n)(-100)
```

The next decision is whether to build another zeta-specific stable
curvature ladder or replace the layer-by-layer endpoint analysis with
a uniform-in-order theorem.

```text
outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md
outputs/jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
