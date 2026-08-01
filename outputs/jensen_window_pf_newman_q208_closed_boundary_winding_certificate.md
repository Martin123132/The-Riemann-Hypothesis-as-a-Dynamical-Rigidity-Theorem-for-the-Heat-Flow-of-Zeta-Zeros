# Newman Q208 Closed-Boundary Winding Certificate

Date: 2026-07-26

Status: rigorous finite Q208 zero-winding and no-contact
certificate. This is not a proof of `Lambda<=0`, RH, or the
Clay prize.

```text
work/rh_compute/results/jensen_window_pf_newman_q208_closed_boundary_winding_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_q208_closed_boundary_winding_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_q208_closed_boundary_winding_certificate.py
```

## Exact Boundary

The positively oriented low rectangle is composed as

```text
bottom: t=1/1040, x=0 -> 246
right:  x=246, t=1/1040 -> 1/5
top:    t=1/5, x=246 -> 0
axis:   x=0, t=1/5 -> 1/1040.
```

Every numerical cell encloses a whole path arc and excludes the
origin. The axis image is positive real because
`F_t(0)=16H_t(0)>0`, `F_t'(0)=0`, and `H_t(0)` increases with
`t`. Exact dyadic witnesses lie in every cyclic adjacent-cell
intersection.

```text
bottom cells=492
right cells=20
reversed top cells=492
axis cells=1
total cells=1005
exact witness-polygon winding=0
```

The proxy `F=16(1+x^4)H` is joined to the ordinary first jet by
a positive-determinant triangular homotopy. Zero proxy winding
therefore equals zero `H+iH_x` winding. Since every heat contact
has strictly positive local index in the standard `(x,t)`
orientation, the low rectangle contains
no contact.

The independent `Lambda<=1/5` theorem excludes contacts for
`1/5<t<=1/4`; the explicitly certified top handles equality.
Thus

```text
Q_208 is contact-free and has zero first-jet winding; therefore Q_1 through Q_208 are certified by containment.
```

## Scope

This is one finite successor to Q207. It supplies no stage
`j>=209`, no parameter-uniform termination theorem, no
`Lambda<=0`, and no RH conclusion.
