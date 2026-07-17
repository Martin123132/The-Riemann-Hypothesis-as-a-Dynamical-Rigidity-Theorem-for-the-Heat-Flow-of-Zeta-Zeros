# Order-Ten Higher-Order Scaling Scout

Date: 2026-07-16

Status: rigorous selected real-interval blocks. This is not a proof;
it is not a
continuous bridge, full-kernel theorem, endpoint tail, or RH proof.

```text
work/rh_compute/results/jensen_window_pf_compound_order10_higher_order_scaling_scout.json
python work/rh_compute/scripts/jensen_window_pf_compound_order10_higher_order_scaling_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_higher_order_scaling_scout.py
```

## Selected Transitions

```text
width 1/8: first selected passing anchor 1350
width 1/4: first selected passing anchor 1500
width 1/2: first selected passing anchor 1800
```

These are strategy transitions, not proved threshold locations.
The blocks between selected anchors remain uncovered.

The scout supports an adaptive design: short blocks near the lower
edge, quarter blocks in the middle transition, and half-unit blocks
from the far lower bridge toward the existing saddle handoff.
