# Order-Ten Seventh-Nested First-Summand Point Scout

Date: 2026-07-16

Status: rigorous selected-point first-summand scout. This is not a proof;
it is not a
continuous curvature theorem, a full-kernel tail theorem, or a proof
of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order10_first_summand_point_scout.json
python work/rh_compute/scripts/jensen_window_pf_compound_order10_first_summand_point_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_first_summand_point_scout.py
```

## New Stable Layer

```text
W_1(t)=8*B_1(t)-w_1(t-1)+2*w_1(t)-w_1(t+1)
z_1(t)=2*w_1(t)-s_1(t)+log(1-exp(-W_1(t)))
```

At each selected exact cache point
`t in {1251,1300,1500,2000,3000,4000}`, all eight coordinates
`B,J,R,S,T,U,V,W` are rigorously positive and

```text
t^2*z_1''(t)<250<5500.
```

The largest selected-point upper enclosure occurs at `t=4000` and is `1.07995129943424296658936866268643391045197629092152697354229E+2`. The smallest selected `W` lower bound is `1.19444524346138206436133765548537325899020640734554457970615E-3`.

At `t=1251` the rigorous ball is

```text
[29.58 +/- 5.22E-3]
```

## What Remains

The headroom is substantial, but isolated points do not control a real
interval. The next proof object must extend the localized tent hierarchy
through `W,z`, cover `t>=1251`, and add a rigorous higher-summand
transfer whose total stays below the full `5500/k^2` budget.

```text
outputs/jensen_window_pf_compound_order10_m100_tail_curvature_reduction.md
outputs/jensen_window_pf_compound_order9_first_summand_curvature_certificate.md
outputs/jensen_window_pf_endpoint_order10_counterexample.md
outputs/formal_core.md
```
