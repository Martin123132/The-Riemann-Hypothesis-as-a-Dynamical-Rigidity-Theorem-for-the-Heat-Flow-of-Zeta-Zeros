# Jensen-Window PF Negative-Lambda High-Order Taylor Scout

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of the
bounded log-curvature theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_high_order_taylor_scout`.

Proof boundary: this artifact certifies higher local Taylor coefficients
and classifies finite truncation samples. It does not prove an infinite
Taylor-tail remainder theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_high_order_taylor_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_high_order_taylor_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_high_order_taylor_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda high-order Taylor scout: 8 coefficient rows, 35 truncation rows, 9 invalid normalizers, 2 upper-wall violations, 3 overbound rows, 0 ready-to-apply rows, 0 issues
```

## Coefficient Scope

The generator forms `P_j(q)` by exact formal series from:

```text
(2*q^2*exp(9u)-3*q*exp(5u))*exp(-q*(exp(4u)-1))
```

It certifies the even local coefficients through `c14` using the same
finite-sum plus geometric-tail method as the lower-order sign scout.

Coefficient signs:

```text
c0: positive, ratio c0/c0 = [1.000000000000000000000000000000000000000 +/- 1e-44]
c2: negative, ratio c2/c0 = [-37.45380985900664233683260390434643972954 +/- 2.43e-39]
c4: positive, ratio c4/c0 = [606.0506446363092395662753197700082678903 +/- 9.08e-39]
c6: negative, ratio c6/c0 = [-5047.989249594588274331383288623799728021 +/- 3.85e-37]
c8: positive, ratio c8/c0 = [12958.15479367270338908734802047778726347 +/- 4.05e-36]
c10: positive, ratio c10/c0 = [207171.8969088401413552699997367309475546 +/- 1.77e-35]
c12: negative, ratio c12/c0 = [-3304359.712814711536315223533315244487346 +/- 3.86e-35]
c14: positive, ratio c14/c0 = [30122498.10781048441714183429334046716634 +/- 4.82e-33]
```

## Truncation Matrix

Rows classify `k=22` and `T=25,50,100,200,500,1000,2000` for
truncation degrees 6, 8, 10, 12, and 14.

```text
degree=6, T=25: F_k=-3.705176510564669608E+3, status=invalid_normalizer, B/bound=n/a
degree=6, T=50: F_k=-4.108228689777087407E+2, status=invalid_normalizer, B/bound=n/a
degree=6, T=100: F_k=-4.077572511834853796E+1, status=invalid_normalizer, B/bound=n/a
degree=6, T=200: F_k=-3.376514867254146312E+0, status=invalid_normalizer, B/bound=n/a
degree=6, T=500: F_k=7.322730386875998071E-2, status=overbound, B/bound=1.362086182073573375E+1
degree=6, T=1000: F_k=4.123450107885821634E-1, status=target_satisfying_truncation, B/bound=4.583846541015247424E-2
degree=6, T=2000: F_k=6.505827654571352895E-1, status=target_satisfying_truncation, B/bound=3.929354468000889985E-3
degree=8, T=25: F_k=7.253037405300258616E+3, status=target_satisfying_truncation, B/bound=5.730442653658323245E-1
degree=8, T=50: F_k=2.740655007638492733E+2, status=target_satisfying_truncation, B/bound=7.734366042747254382E-1
degree=8, T=100: F_k=2.029797990498837916E+0, status=overbound, B/bound=1.514094382351363726E+1
degree=8, T=200: F_k=-7.011696729511853196E-1, status=invalid_normalizer, B/bound=n/a
degree=8, T=500: F_k=1.417161408429157821E-1, status=target_satisfying_truncation, B/bound=9.195506423257628972E-1
degree=8, T=1000: F_k=4.166255630994669010E-1, status=target_satisfying_truncation, B/bound=1.974853208203007153E-2
degree=8, T=2000: F_k=6.508502999765655856E-1, status=target_satisfying_truncation, B/bound=3.154227475054448076E-3
degree=10, T=25: F_k=1.929621942891629552E+5, status=target_satisfying_truncation, B/bound=5.676012864611807968E-1
degree=10, T=50: F_k=6.077476653384558541E+3, status=target_satisfying_truncation, B/bound=5.775318061363339222E-1
degree=10, T=100: F_k=1.833863965098960025E+2, status=target_satisfying_truncation, B/bound=6.016666587380244242E-1
degree=10, T=200: F_k=4.966224030779976075E+0, status=target_satisfying_truncation, B/bound=6.020892333219957254E-1
degree=10, T=500: F_k=1.997502523691228748E-1, status=upper_wall_violation, B/bound=-7.733885588213444624E-1
degree=10, T=1000: F_k=4.184391290846608726E-1, status=target_satisfying_truncation, B/bound=4.447999641784707720E-3
degree=10, T=2000: F_k=6.509069739136028972E-1, status=target_satisfying_truncation, B/bound=2.907200831239827433E-3
degree=12, T=25: F_k=-3.065273367629861797E+6, status=invalid_normalizer, B/bound=n/a
degree=12, T=50: F_k=-4.483245400160020321E+4, status=invalid_normalizer, B/bound=n/a
degree=12, T=100: F_k=-6.120812699742408998E+2, status=invalid_normalizer, B/bound=n/a
degree=12, T=200: F_k=-7.462958258034663024E+0, status=invalid_normalizer, B/bound=n/a
degree=12, T=500: F_k=1.488403217141381130E-1, status=overbound, B/bound=1.109778718884260427E+0
degree=12, T=1000: F_k=4.176436614181767357E-1, status=target_satisfying_truncation, B/bound=1.343386392858017701E-2
degree=12, T=2000: F_k=6.508945447313140826E-1, status=target_satisfying_truncation, B/bound=2.982498591564098327E-3
degree=14, T=25: F_k=3.079504194598911611E+7, status=target_satisfying_truncation, B/bound=7.559222645144021739E-1
degree=14, T=50: F_k=2.197012593860480617E+5, status=target_satisfying_truncation, B/bound=7.680279120590501359E-1
degree=14, T=100: F_k=1.454588365866761170E+3, status=target_satisfying_truncation, B/bound=7.836668877501897378E-1
degree=14, T=200: F_k=8.682898271973165644E+0, status=target_satisfying_truncation, B/bound=7.755156486980905861E-1
degree=14, T=500: F_k=1.752936930529029395E-1, status=upper_wall_violation, B/bound=-3.820397829277359356E-1
degree=14, T=1000: F_k=4.178503283817608359E-1, status=target_satisfying_truncation, B/bound=1.042608934307475150E-2
degree=14, T=2000: F_k=6.508961593169670834E-1, status=target_satisfying_truncation, B/bound=2.969632014203448549E-3
```

Summary counts:

```text
invalid normalizers: 9
upper-wall violations: 2
overbound rows: 3
target-satisfying truncation rows: 21
```

Interpretation: higher finite Taylor truncation is a theorem-search
diagnostic, not a replacement for an infinite Taylor-tail estimate.
The local/mesoscopic route still needs an analytic remainder theorem
controlling positivity and second/third log-differences.

Integration:

```text
outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md
outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
```

Summary:

High-order local Taylor truncations through degree 14 do not provide a stable finite proof model: the k=22 sample matrix contains invalid normalizers, overbound log-curvature rows, and upper-wall violations, so the proof route still needs an analytic infinite-tail remainder theorem rather than a higher finite truncation.
