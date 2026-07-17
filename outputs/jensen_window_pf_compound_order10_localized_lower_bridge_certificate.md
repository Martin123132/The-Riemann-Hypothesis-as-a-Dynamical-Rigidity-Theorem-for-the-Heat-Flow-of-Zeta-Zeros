# Order-Ten Localized Lower-Bridge Certificate

Date: 2026-07-16

Status: rigorous first-summand continuous curvature theorem on
`1251<=t<=5700`. This is not yet a full-kernel or RH theorem.

```text
work/rh_compute/results/jensen_window_pf_compound_order10_localized_lower_bridge_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order10_localized_lower_bridge_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_localized_lower_bridge_certificate.py
```

## Theorem

```text
z_1''(t)<=4200/t^2 for every real 1251<=t<=5700
z''=2*w''-s''+phi(W)*W''-chi(W)*(W')^2
z''<=2*w''-s''+phi(W)*max(W'',0)
```

The bridge contains `9,996` contiguous real-interval blocks:

```text
near quarter blocks: 996
middle quarter blocks: 1200
far half blocks: 7800
largest scaled upper: 4.09781412038460013699765821854500983262668416719590010539729E+3
smallest cap margin: 1.40218587961539986300234178145499016737331583280409989460271E+3
```

Every serialized Arb ball was parsed at 512-bit precision. The complete
segment cache and all source caches are SHA-bound by the run contract.

The next analytic task is the upper first-summand saddle range from
`t=5700` onward, followed by the first-summand-to-full-kernel transfer.
