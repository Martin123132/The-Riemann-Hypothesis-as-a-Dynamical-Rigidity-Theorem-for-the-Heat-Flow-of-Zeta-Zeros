# Event-zero transport recentered at the selector

Date: 2026-08-13

Status: rigorous event-zero exact-minus-beta-four top-corridor enclosure;
not a complete selected-branch or `T_upper` splice

The earlier extension transported seven wide event-centered moment balls over
`0<=theta<=2`.  Its `2.141e-6` majorant was dominated by lost interval
correlation even though the selector-center midpoint was near `10^-13`.

Recenter the exact transport identity at `t*` instead:

```text
D_(0,t*+h)(z)=exp(ihz/beta)D_(0,t*)(z),
h=(pi/16)theta,       -1<=theta<=0.                  (RC1)
```

The `pi/16` width is inherited from the exact selector/event geometry.  A
fresh contour calculation through moment degree `8`, rather
than reparsing the old moment balls, gives

```text
sup_(t*-pi/16<=t<=t*) |Delta I_2(0,t)|
 < [3.25920806144198526443489617522953081171057582794588981743386534007204599693163482941344131e-13 +/- 2.75e-103] <1e-8,
physical < [1.46715909773349559469904100090583557646428682086563789121566081377694687079864846688898606e-13 +/- 5.08e-103]. (RC2)
```

Direct contour integrations at both faces overlap the recentered transport
balls.  This replaces the old top-corridor majorant for event zero; it does
not alter or invalidate the original all-event cells.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_selector_center_recentered_transport_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_selector_center_recentered_transport_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_selector_center_recentered_transport_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_selector_center_recentered_transport_gate.py
```

No remaining-398-mode finite-integral theorem, completed-remainder bound,
all-corridor continuation, complete `Q_K-T` or `T_upper`, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
