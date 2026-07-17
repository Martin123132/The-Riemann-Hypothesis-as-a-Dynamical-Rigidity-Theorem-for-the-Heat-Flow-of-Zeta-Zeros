# Order-Ten Localized Formula Pilot Certificate

Date: 2026-07-16

Status: rigorous first-summand local curvature theorem. This is not a proof
of RH and does not cover the full half-line, full Newman kernel, or order-ten endpoint tail.

```text
work/rh_compute/results/jensen_window_pf_compound_order10_localized_formula_pilot_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order10_localized_formula_pilot_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_localized_formula_pilot_certificate.py
```

## Theorem

```text
z_1''(t)<=5500/t^2 for every real 1251<=t<=1251.5
```

The interval is split into two exact width-`1/4` blocks, expanded
rightward from `1251` and leftward from `1251.5`. Every block
uses a hundredth-tile `H2-H24` common remainder and an exact
eighth-grid `H0-H8` point jet parsed at 512-bit precision.

```text
z''=2*w''-s''+phi(W)*W''-chi(W)*(W')^2 <=2*w''-s''+phi(W)*max(W'',0)
largest scaled curvature upper: 3.77457452093705165580718107676218208771995796660502457346484E+2
smallest cap margin lower: 5.12254254790629483441928189232378179122800420333949754265351E+3
```

The monolithic sixth-order composed-log remainder is deliberately not
used. Keeping `s`, `w`, and `W` separate preserves the cancellation that
is visible in the exact point jets while retaining a real-interval proof.

The next task is to extend this economical higher-order component bound
through the adaptive lower regimes and into the saddle handoff.
