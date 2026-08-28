# Complete absolute bound for the translated A-face endpoint layer

Date: 2026-08-23

Status: rigorous local interval certificate at `t=10^10`, independently
replayed on a different interval partition; not a proof of `R_after_A` or RH

Let `D_(tr,8)` be the exact eight-cell Fresnel midpoint defect from Section
11.455 and `D_(tail,34),6` the exact rational-log six-current tail.  The
positive midpoint Peano kernel gives

```text
|D_(tr,8)(x)|
 <=sum_(m=39895)^39902 integral_(m-1/2)^(m+1/2)
       K_m(y)|G_(A,yy)(y,x)|dy.                    (EP1)
```

For a real normal coordinate `q0<0`, put

```text
w0=(-q0)sqrt(pi)(1-i)/2,
H(q0)=(1+i)exp(w0^2)erfc(w0)/2.                    (EP2)
```

If an interval for `q` has midpoint `q0` and radius `r`, the exact integrating
factor for `H_q=1-i*pi*qH` and `|exp(i*theta)-1|<=|theta|` give

```text
|H(q)-H(q0)|
 <=r[1+pi|H(q0)|(|q0|+r/2)].                       (EP3)
```

Thus (EP2)--(EP3), followed by the exact ODE for `H_qq`, enclose `G_(A,yy)`
without complex-box evaluation of `erfc`.  Boxes with certified `|q|>=2`
instead use the six-current expansion and its exact remainder; if that
alternating interval expression is indeterminate, the box fails over to
(EP2)--(EP3).

The 34-cell tail is not discarded.  Its exact rational-log defect is enclosed
on every `x` slab and its certified six-current replacement error is attached
after integration.  Consequently

```text
B_end
 <=2(pi/(32T))^(1/4) integral_(x_16)^(1/2)
      x^(-5/4)(1-x)^(-1/4)
      [B_8(x)+|D_(tail,34),6(x)|]dx
    +epsilon_(tail,6).                              (EP4)
```

The builder covers the slightly larger rational-decimal interval
`[0.4999028779462474,1/2]` by `256` `x` slabs and splits each half-cell into
`16` normal boxes.  Arb certifies

```text
eight-cell contribution       <= [0.003583797688220635967531828752917358924911078910548908869597760022053758 +/- 3.25e-73],
34-cell six-current tail      <= [1.780116174555390927677664839309381722377458211022801365518606313071404e-5 +/- 1.79e-5],
tail replacement error       <= [5.407920667053412976178242568568276672982652639121521449041334065558340e-15 +/- 2.09e-85],

B_end                         <= [0.003601598849971597797475658814286630984703421769332119535892067606633513 +/- 1.79e-5]
                              < 0.0037.             (EP5)
```

The independent checker repeats (EP1)--(EP5) at 120 decimal digits on a
`320 x (8*40)` partition and obtains its own bound below `0.0030`; it does not
import the builder or accept sampled floating values.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_peano_interval_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_peano_interval_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_peano_interval_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_peano_interval_gate.py
```

Pi provenance: every `pi` in (EP1)--(EP5) comes from the canonical Fresnel
primitive, its exact ODE, or the inherited equation-(9) normalization.  No
geometric or fitted occurrence of `pi` is inserted.

Proof boundary: (EP5) proves only the complete translated A-face endpoint
layer on `x_16<=x<=1/2` at the saved height.  It does not prove the separate
pre-endpoint Morse integral, the full signed A-face channel, the joined
`R_after_A`, `R_Dir`, or `Q_K-T` bound, an all-height theorem, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
