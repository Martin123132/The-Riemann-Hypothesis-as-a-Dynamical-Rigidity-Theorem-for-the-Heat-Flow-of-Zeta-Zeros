# Jensen-Window PF Compound Order-Nine Uniform Tail And Flow Reduction

Date: 2026-07-13

Status: exact uniform eventual signed order-nine tail and conditional
forward reduction with one open `lambda=-100` entry. This is not a
proof of all-shift order nine, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_uniform_tail_flow_reduction.py
```

## Universal Signed Tail

At `m=9`, `D=36` and `epsilon_9=1`. The all-fixed-order theorem gives

```text
det K=347485857744891213250560000*G_2^36*h^36+o(h^36)
there exists N_9 such that Q_(9,n)(lambda)>0 for every n>=N_9 and -100<=lambda<=0
```

## Endpoint Coordinate And Flow

```text
Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)
Q_(9,n)>0 iff Q_(8,n+1)^2>Q_(8,n)*Q_(8,n+2)
Q_(9,n)'=(4*n+66)*delta(Q_(9,n))
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n=(4*n+66)*Q_(8,n)/Q_(8,n+1)>0, b_n=((4*n+66)/(4*n+62))*(log Q_(8,n+1))'
```

The completed `Q_7,Q_8>0` cone makes the off-diagonal coefficient
strictly positive. The eventual tail and variation of constants prove
the forward and arbitrary-column conclusions conditional only on entry.

## Countermodel Gate

The exact sequence `(1,1,1/2!,...,1/15!,1/20926400000000)` has every
available signed contiguous minor through order eight strictly positive,
but

```text
Q_(8,1)^2-Q_(8,0)Q_(8,2)=-40921/24538241725813248696877825138277835241884876724206315993476458291926306465427019677040640000000000000000000000000000<0,
Q_(9,0)=-40921/1107954002817773794682431864791874987650437949161472000000000000000000<0.
```

Thus lower-cone promotion cannot replace a new zeta-specific theorem.

## Remaining Endpoint Target

```text
Q_(9,n)(-100)>0 for every n>=0, equivalently strict log-concavity of Q_(8,n)(-100)
```

```text
outputs/jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.md
outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
