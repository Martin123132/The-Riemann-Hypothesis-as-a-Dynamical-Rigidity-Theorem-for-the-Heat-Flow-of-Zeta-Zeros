# Quartic Outer-Branch Length-13 Obstruction

Date: 2026-07-25

Status: exact fixed-prefix length-thirteen obstruction. This is
not a proof of a uniform quartic threshold, PF-infinity,
`Lambda <= 0`, or RH.

Artifact kind: `jensen_window_pf_quartic_outer_branch_length13_obstruction`.

```text
work/rh_compute/results/jensen_window_pf_quartic_outer_branch_length13_obstruction.json
python work/rh_compute/scripts/jensen_window_pf_quartic_outer_branch_length13_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_length13_obstruction.py
```

Current result:

```text
validated quartic outer-branch length-13 obstruction: 10 rows, 0 issues, 3 exact corridor parameters, 3 derivative certificates, 1516 derivative Bernstein coefficients, 8 denominator-factor Bernstein coefficients, 1 negative cube-corner maximum, 1 forbidden all-length promotion, 1 repaired handoff
```

## Exact Parameterization

Keep the rational outer-contact prefix `x_2,...,x_9` from the
order-four nonpromotion gate. Put

```text
G_j=d_(j+2)^2-x_(j+2)^2*d_(j+1)*d_(j+3)
C_k=G_(k-4)^2/(x_(k-2)^3*G_(k-5))
```

For each `k=10,11,12`, the strict terminal order-three and
order-four signs are exactly parameterized by

```text
G_(k-3)=y_k*C_k,  0<y_k<1.
```

This is exhaustive: `G_(k-3)>0` is the order-three sign and
`G_(k-3)<C_k` is the order-four sign.
The checker independently reconstructs the same cube by solving
the raw terminal `3x3` and `4x4` Hankel determinant roots at each
step, then reproduces the derivative Bernstein hashes.

## Terminal Compatibility

At the next step define

```text
M_13=d_12^2-x_12^2*d_11*d_12,
Delta_13=C_13-M_13.
```

If `x_13>x_12`, then the new gap `G_10(x_13)` is strictly
larger than `M_13`. The next order-four sign requires
`G_10(x_13)<C_13`. Hence `Delta_13>0` is necessary.

Exact tensor-product Bernstein conversion gives the derivative
numerator sign counts

```text
y_10: +612, -0, 0=0
y_11: +0, -480, 0=64
y_12: +0, -216, 0=144.
```

Together with the independently certified denominator signs,
all three partial derivatives of `Delta_13` are nonnegative on
the open cube. Therefore its supremum is the closed-cube corner
`(1,1,1)`. Exact substitution gives

```text
Delta_13(1,1,1)=-524254295062439831022231455843667625869467186452689465245719504775268755498243950341094754949683301859298469682457827876218799379390895688576267963329930419105243738286554439443726602730795165957041668602994459185114645569013951551839026537101979584509791546986995998894738450218260389107200405986711153542269363281233626701659419093750000000000000000000000000000000000000000000000000/344655609182286642995773369360966254780375865094480599264341686223734721504950962574603697173070444997375009363172701420864420824669985673172055889071113287856553576756107122605282798363163431647862915186523301539756509214311081629497085466980262738678116158104540923683305004724600266303649550993239396090375856179049598973088251326796535248225106780557582258791817078543175255883015503129
                    =-0.00000152109607705576128738324556976<0.
```

Thus no choice inside any of the three signed continuation
corridors can produce an increasing `x_13` while preserving the
next order-four sign.

## Corrected Handoff

The previously displayed `x_10` corridor was real, and the
branch can continue through `x_12`, but it cannot continue
indefinitely inside the increasing-contraction cone. The live
question is no longer whether this particular witness has an
all-length extension. It does not.

This does not promote the quartic threshold. The companion
alternate-tail gate now gives an explicit rational continuation
from the same outer contact that survives through `x_13`, so
this obstruction is demonstrably prefix-specific. The next
theorem target is a uniform finite-tail obstruction beginning
at length 14, or another tail crossing that boundary. Global
Xi, degree-five, and
theta-specific routes remain open.

```text
outputs/jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.md
```
