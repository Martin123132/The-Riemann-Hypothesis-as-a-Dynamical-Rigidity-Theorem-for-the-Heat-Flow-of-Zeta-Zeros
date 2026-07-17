# Order-ten global first-summand curvature certificate

Date: 2026-07-17

Status: **rigorous global order-ten first-summand curvature theorem on t>=1251**. This certificate is not a proof of PF-infinity, RH, or `Lambda <= 0`.

## Theorem

`z_1''(t)<=4200/t^2 for every real t>=1251`.

The largest source upper bound is `4.09781412038460013699765821854500983262668416719590010539729E+3`.

## Boundary

This artifact composes the continuous first-summand theorem at lambda=-100 only. It does not prove a full-Newman-kernel ceiling, endpoint entry, heat-forward invariance, the Jensen hierarchy, Lambda<=0, or RH.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_compound_order10_first_summand_curvature_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_first_summand_curvature_certificate.py
```
