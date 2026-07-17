# Order-ten compact first-summand curvature certificate

Date: 2026-07-17

Status: **rigorous order-ten first-summand curvature theorem on 5700<=t<=38020**. This certificate is not a proof of PF-infinity, RH, or `Lambda <= 0`.

## Theorem

`z_1''(t)<=4200/t^2 for every real 5700<=t<=38020`.

The exact partition uses unit blocks from `5700` through `10000` and width-two blocks from `10000` through `38020`.
All `18310` blocks pass; the largest scaled upper bound is `3.10527904794520576042110924623685242329682444671709504636403E+3`.
The weakest proved `W` floor is `1.63273977810527135997728978132745155444550432524586757993867E-4`.

## Transition

The rigorous upper enclosure `V'(2.001) <= 3.80196621635193385678801171533426331770390296897595492274576E+4` is below `38020`, so this compact theorem overlaps the finite saddle ray.

## Boundary

This certificate concerns the first Newman summand at lambda=-100. It does not yet transfer the theorem to the full Newman kernel and does not prove RH, the Jensen hierarchy, or Lambda<=0.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_compound_order10_nested_curvature_compact_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_nested_curvature_compact_certificate.py
```
