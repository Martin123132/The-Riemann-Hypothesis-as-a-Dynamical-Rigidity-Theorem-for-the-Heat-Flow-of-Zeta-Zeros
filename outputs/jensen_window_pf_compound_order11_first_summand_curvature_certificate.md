# Order-Eleven Global First-Summand Curvature Certificate

Date: 2026-07-22

Status: **rigorous global order-eleven first-summand curvature theorem on t>=1252**. This is not a proof of RH or `Lambda <= 0`.

## Theorem

`y_1''(t)<=6000/t^2 for every real t>=1252`.

The four exact source ranges are:

- `y_1''(t)<=6000/t^2 for every real 1252<=t<=5700`
- `y_1''(t)<=6000/t^2 for every real 5700<=t<=38020`
- `y_1''(t)<=6000/t^2 for every saddle mode 2001/1000<=u<=20`
- `t^2*y_1''(t)<2000 for every mode u>=20`

The compact-saddle overlap is `V'(2001/1000)<=3.80196621635193385678801171533426331770390296897595492274576E+4<38020`.
The largest scaled upper across the four sources is `2122.86516691011593226179209689464539687227984653146836129212<6000`.

## Boundary

This artifact composes the continuous first-summand theorem at lambda=-100 only. It does not prove full-Newman-kernel entry, heat-forward invariance, all-shift order eleven, PF-infinity, Lambda<=0, or RH.

This is not a proof of RH.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_compound_order11_first_summand_curvature_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_first_summand_curvature_certificate.py
```
