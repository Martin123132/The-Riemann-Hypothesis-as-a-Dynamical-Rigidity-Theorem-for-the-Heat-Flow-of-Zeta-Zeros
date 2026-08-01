# Order-Eleven Lambda=-100 Entry Certificate

Date: 2026-07-22

Status: **rigorous all-shift signed order-eleven entry at lambda=-100**. This is not a proof of RH or `Lambda <= 0`.

## Composition

`z_1''(t)<=4200/t^2 for every real t>=1251`
`y_1''(t)<=6000/t^2 for every real t>=1252`
`Q_(11,n)(-100)>0 for every n>=1243`
`Q_(11,n)(-100)>0 for every 0<=n<=1242`

Together these prove `Q_(11,n)(-100)>0 for every integer n>=0`.

## Boundary

This artifact proves fixed order-eleven contiguous signed Hankel positivity at lambda=-100. It does not prove the heat interval, all orders, PF-infinity, Lambda<=0, or RH.

This is not a proof of RH.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_compound_order11_m100_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_m100_entry_certificate.py
```
