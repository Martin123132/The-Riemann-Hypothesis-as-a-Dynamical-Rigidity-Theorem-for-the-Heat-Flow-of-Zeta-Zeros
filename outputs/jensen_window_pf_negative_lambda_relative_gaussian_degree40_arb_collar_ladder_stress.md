# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-40 Arb Collar Ladder Stress

Date: 2026-07-07

Status: finite theorem-search diagnostic. This is not a proof
of an infinite Taylor-tail theorem, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress`.

Proof boundary: this artifact stress-tests finite Arb coefficient-ball
surrogates through degree 40 on the same real-`T` collar. It does not
bound the infinite residual Taylor tail.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress: 8 rows, 0 issues, 13 degree levels, max degree 40, 0 failed Bernstein rows, 0 ready-to-apply rows
```

## Degree Ladder

```text
real T interval: [1156, infinity)
u interval: [0, 1/1156]
degree values: [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]
degree levels: 13
all degree levels certified: True
total failed Bernstein rows: 0
weakest stencil lower: degree 16 companion_product = [0.0630023568050961538569390725332 +/- 2.64e-33]
```

Per-degree stress rows:

```text
degree=16, M=8: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [0.0630023568050961538569390725332 +/- 2.64e-33]
degree=18, M=9: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [22.1926895649987730183430685674 +/- 2.26e-29]
degree=20, M=10: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [15.8977949886685698343390329961 +/- 1.94e-29]
degree=22, M=11: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.5278135161576143829717829687 +/- 1.78e-29]
degree=24, M=12: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.1720202438784095917682475691 +/- 4.79e-29]
degree=26, M=13: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.2292873126090604937716029993 +/- 4.20e-29]
degree=28, M=14: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.2266586021971187019192956599 +/- 4.84e-29]
degree=30, M=15: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.2239315025138966525431637249 +/- 3.34e-29]
degree=32, M=16: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.2254369739258098285094603945 +/- 4.38e-29]
degree=34, M=17: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.2248865820116223134229674160 +/- 4.49e-29]
degree=36, M=18: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.2250578182865955154436111184 +/- 2.41e-29]
degree=38, M=19: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.2250081542294786587911012907 +/- 1.84e-29]
degree=40, M=20: normalizers=4/4, stencils=3/3, failed Bernstein=0, weakest stencil=companion_product [17.2250223756550133023672457518 +/- 2.10e-29]
```

## Degree-40 Endpoint Snapshot

```text
highest ratio c40/c0: [4944469183781380656.77812493750 +/- 3.86e-13] (positive)
B_product endpoint: [37.2033572903651204156102209610 +/- 2.67e-29]
companion_product endpoint: [17.2250223756550133023672457518 +/- 2.10e-29]
weighted_gap_derivative endpoint: [27.4244024269680055723947152492 +/- 3.51e-29]
```

Interpretation:

The degree-16 Arb collar was not an isolated finite-truncation fluke
through the next twelve even truncation levels. The remaining task is
now cleaner but still open: replace finite degree-by-degree stress with
a signed residual-tail bound beyond degree 40 on the same collar, then
remove the fixed-`k` limitation.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The same real-T collar T>=1156 survives Arb/Bernstein finite-surrogate stress for every even Taylor degree 16 through 40 at fixed k=22: all four normalizers and the B, companion, and weighted-gap stencils certify positive at every tested degree, with zero Bernstein failures. This makes the remaining gap sharper, but still finite: the proof still needs a residual-tail bound beyond the last tested degree and then removal of the fixed-k limitation.
