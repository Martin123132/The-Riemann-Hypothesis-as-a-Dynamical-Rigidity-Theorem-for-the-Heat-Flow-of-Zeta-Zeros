# Jensen-Window PF Negative-Lambda Scaled-Curvature Continuous Bridge

Date: 2026-07-11

Status: interval and analytic all-k scaled-curvature theorem at lambda=-100.
This is not a proof of PF-infinity, the all-order Jensen bridge, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py
```

Current result:

```text
validated Jensen-window PF negative-lambda scaled-curvature continuous bridge: 10 rows, 0 issues, 16074 compact blocks, 318 prefix gaps, 1 analytic ray, 1 all-k scaled-curvature theorem, 0 open requirements
```

## Exact Continuous Bridge

For the first Newman summand, set `H(t)=log Gamma(t+1/2)-log M_t^(1)` and

```text
B(t)=H(t+1)-2*H(t)+H(t-1),
C(t)=(2*t+1)*B(t).
```

The centered second difference has its tent-kernel representation. With
`Q(r)=2*H''(r)+(2*r+1)*H'''(r)`, differentiation gives

```text
C'(t)=integral_[-1,1] (1-|s|)*(Q(t+s)-2*s*H'''(t+s)) ds.
```

Thus `Q(r)-2*abs(H'''(r))>64/(r-3)^5` for `r>=318` proves
strict first-summand growth and leaves enough room for the complete-kernel transfer.

## Compact Certificate

Interval composite Simpson quadrature with a fourth-derivative remainder proves
the buffered inequality continuously on `0.9264<=u<=5`.

```text
compact blocks: 16074
minimum shifted-buffer lower bound: 2.79390590809330483221501780691349081277920947608007082173275E-13
minimum full-kernel transfer margin: 2.79390590385814009594351611157187830860135660899466919665804E-13
regions: {'low': 1074, 'mid_2_3': 5000, 'mid_3_4': 5000, 'high_4_5': 5000}
```

The low region uses the established `|y|<=6`, `V^(8)/a^4<=1/50000` envelope.
A separate 3,000-interval gate proves `V^(8)/a^4<=1/10^10` on the
`|y|<=8` high-mode windows.

## Analytic Ray

For `u>=5`, `q=pi*exp(4u)>=10^9`. Exact polynomial algebra proves

```text
u^2*t*((2*t+1)*V'''/V''^3-2/V'') >= 1/5.
```

After `u=5+v`, `q=10^9+Q`, the cleared numerator has
66 strictly positive coefficients.
The paired ray moment errors, Gamma-series bounds, and `H'''` buffer consume

```text
normalized error <= 2106009471479000117/29250000000000000000
normalized final margin >= 3743990528520999883/29250000000000000000 > 1/10.
```

This also dominates `64/(t-3)^5` on the whole ray.

## Full-Kernel Transfer And Prefix

If `M_k=M_k^(1)*(1+delta_k)` with `0<=delta_k<=2/k^6`, direct expansion gives

```text
abs((C_(k+1)-C_k)-(C1_(k+1)-C1_k)) <= 64/(k-1)^5, k>=319.
```

The repaired Arb source independently proves

```text
C_(k+1)-C_k>0 for k=1..318 (318 rows),
minimum prefix margin=7.86819705122683963600171680885401933766308105326326491770645E-4 at k=318.
```

Therefore the actual lambda=-100 sequence satisfies

```text
C_(k+1)>=C_k for every integer k>=1.
```

Combined with the already proved `B_(k+1)<=B_k` and `B_k>=0`, this closes
the two-sided raw-ratio corridor at lambda=-100. It does not establish any
higher-degree minor cone or the all-order Jensen/PF bridge.

```text
outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md
outputs/signed_hankel_jensen_dependency_graph.md
```
