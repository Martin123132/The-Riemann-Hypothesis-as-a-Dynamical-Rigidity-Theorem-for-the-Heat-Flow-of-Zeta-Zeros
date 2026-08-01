# Order-Eleven Compact Adaptive-H23 Certificate

Date: 2026-07-22

Status: **rigorous order-eleven first-summand curvature theorem on 5700<=t<=38020**. This is not a proof of RH or `Lambda <= 0`.

## Theorem

`y_1''(t)<=6000/t^2 for every real 5700<=t<=38020`.

All `2020` segments and `129280` quarter blocks pass.
The largest scaled upper is `2122.86516691011593226179209689464539687227984653146836129212` and the smallest margin is `3877.13483308988406773820790310535460312772015346853163870788`.

## Saddle Overlap

The rigorous upper enclosure `V'(2001/1000)<=3.80196621635193385678801171533426331770390296897595492274576E+4<38020` joins the compact interval to the monotone saddle ray.

## Boundary

This artifact proves the compact first-summand theorem at lambda=-100 only. It does not by itself prove the global first-summand theorem, full-kernel transfer, PF-infinity, Lambda<=0, or RH.

This is not a proof of RH.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.py
```
