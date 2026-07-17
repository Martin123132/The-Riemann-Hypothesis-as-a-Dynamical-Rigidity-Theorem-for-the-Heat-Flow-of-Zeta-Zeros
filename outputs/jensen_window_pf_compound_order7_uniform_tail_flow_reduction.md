# Jensen-Window PF Compound Order-Seven Uniform Tail And Flow Reduction

Date: 2026-07-13

Status: exact uniform eventual signed order-seven tail and conditional
forward reduction with one open `lambda=-100` entry. This is not a
proof of all-shift order seven, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_uniform_tail_flow_reduction.py
```

## Universal Signed Tail

The graded-kernel Vandermonde lemma applies at `m=7`, where
`D=binom(7,2)=21` and `epsilon_7=-1`. It gives

```text
det K=-52183852646400*G_2^21*h^21+o(h^21)
Q_7=positive_scale*(52183852646400*G_2^21*h^21+o(h^21))
there exists N_7 such that Q_(7,n)(lambda)>0 for every n>=N_7 and -100<=lambda<=0
```

The threshold is finite but non-effective. This is an eventual tail,
not an all-shift order-seven theorem.

## Endpoint Coordinate

Signed Desnanot-Jacobi gives

```text
Q_(7,n)*Q_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)*Q_(6,n+2)
Q_(7,n)*H_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)*Q_(6,n+2)
Q_(7,n)>0 iff Q_(6,n+1)^2>Q_(6,n)*Q_(6,n+2)
```

The completed positive `Q_5=H_5` and `Q_6=-H_6` layers make this an
exact strict log-concavity target for the endpoint `Q_6` sequence.

## Cooperative Heat Flow

The general affine-Hankel and adjacent Plucker identities specialize to

```text
Q_(7,n)'=(4*n+50)*delta(Q_(7,n))
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n=(4*n+50)*Q_(6,n)/Q_(6,n+1)>0, b_n=((4*n+50)/(4*n+46))*(log Q_(6,n+1))'
```

Since `Q_(6,n)>0` on the complete heat interval, `a_n>0`. The uniform
eventual tail and variation of constants therefore prove

```text
[Q_(7,n)(-100)>0 for all n] => [Q_(7,n)(lambda)>0 for all n and -100<=lambda<=0]
completed signed contiguous layers through order seven imply every arbitrary-column signed layer through order seven
```

conditional only on all-shift endpoint entry.

## Countermodel Gate

The exact rational sequence

```text
(1,1,1/2,1/6,1/24,1/120,1/720,1/5040,1/40320,1/362880,1/3628800,1/39916800,1/480100000)
```

has every available signed contiguous minor through order six strict
and positive, but

```text
Q_(6,1)^2-Q_(6,0)Q_(6,2)=-34201/59436732624881410374275156859658650098073600000000000000000<0,
Q_(7,0)=-34201/88845982207870628726518579200000000000<0.
```

Thus no lower-cone promotion can replace the new endpoint theorem.

## Remaining Endpoint Target

The sole fixed-order-seven input still missing is

```text
Q_(7,n)(-100)>0 for every n>=0, equivalently strict log-concavity of Q_(6,n)(-100)
```

A rigorous coefficient prefix and a cancellation-preserving analytic
tail coordinate are the next two sub-obligations.

```text
outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md
outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
