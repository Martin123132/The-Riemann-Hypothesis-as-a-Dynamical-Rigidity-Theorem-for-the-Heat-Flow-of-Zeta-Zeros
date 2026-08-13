# Event-0 transport through the selector top corridor

Date: 2026-08-13

Status: rigorous one-mode fold/Fresnel transport extension; not a complete
ordinary splice, `T_upper` theorem, or proof of RH

For the first event mode `m=39894`,

```text
tau_0=t*-pi/8,        h=(pi/16)theta.                  (TC1)
```

The exact-minus-beta-minus-four integrand obeys the already certified identity

```text
D_(0,tau_0+h)(z)=exp(i*h*z/beta)D_(0,tau_0)(z).        (TC2)
```

The existing atlas used moments through degree `6` on
`|theta|<=1`.  Reusing those stored interval balls, multiplying the kth
coefficient by `2^k`, and recomputing the exponential and lifted-tail
remainders proves the larger one-sided interval

```text
0<=theta<=2,
tau_0<=t<=t*,
|Delta I_2(0,t)|
 < [2.141000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e-6 +/- 4.41e-107] < 2.15e-6,
physical < [9.640000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e-7 +/- 2.51e-107]. (TC3)
```

The top corridor is exactly `1<=theta<=2`, because its lower face is
`[10000011008.84623879165060443333514568144680285591602134770666267230484984016248766209627646889114177 +/- 5.29e-90]` and `theta=2` is the selector center
`[10000011009.04258833249996651073906089290177178617834443516760678611578387718172618748916453130539395 +/- 1.84e-90]`.  Thus the edge mode remains fold-owned all
the way from its certified event face to the completed selector object; it
need not be compared there with the invalid bare full-saddle carrier.

This does not identify the one-mode transport ball with the complete 399-mode
selector strip.  At `t*`, the mode is only one summand inside the selected
branch projection, whose opposite branch, outer complement, and endpoint
half-current remain grouped.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate.py
```

No remaining-398-mode grouped estimate, finite-integral ordinary amplitude
theorem, all-corridor continuation, complete `Q_K-T` or `T_upper`,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
