# Jensen-Window PF Negative-Lambda Relative-Gaussian Formal-Tail Obstruction Scout

Date: 2026-07-07

Status: finite theorem-search obstruction. This is not a proof
of an analytic residual estimate, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout`.

Proof boundary: this artifact rejects naive formal-tail summation
templates. It does not reject or prove an actual contour, saddle,
integral-remainder, or optimally truncated asymptotic remainder theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout: 8 rows, 0 issues, 4 profile rows, 4 formal-tail turnaround rows, 0 ready-to-apply rows
```

## Formal-Term Profile

```text
max Taylor degree: 240
tested j range: 21..120
u interval: [0, 1/1156]
value budget A: 5.382819486765314521E-01
derivative budget B: 9.315354075509573936E-03
highest coefficient ratio: c240/c0 = [+/- 8.82e+111]
```

Per-index profile:

```text
F_21: least value j=103 [1.23501699002924070e-33 +/- 4.21e-51], least derivative j=103 [1.10040441594033041e-34 +/- 3.32e-52], first value growth after j>=80: j=103 ratio=[4.38508892036147272 +/- 5.75e-19], first derivative growth after j>=80: j=103 ratio=[4.42766258567677186 +/- 3.96e-18]
  value sums j21..40=[6.75772978733390378e-5 +/- 4.86e-23] (budget fraction [0.000125542567495511743 +/- 2.90e-23]), j41..80=[3.75588531519403356e-15 +/- 4.53e-34] (budget fraction [6.97754276254032287e-15 +/- 3.89e-34]), j81..120=[7.47031824300903847e-24 +/- 1.30e-42] (budget fraction [1.38780768357107956e-23 +/- 4.26e-41])
  derivative sums j21..40=[1.24827307244259553e-6 +/- 2.35e-24] (budget fraction [0.000134001677480446364 +/- 4.42e-23]), j41..80=[1.34951077916362017e-16 +/- 3.20e-34] (budget fraction [1.44869509867750082e-14 +/- 3.24e-32]), j81..120=[7.73085086816123746e-25 +/- 4.76e-43] (budget fraction [8.29904135204687852e-23 +/- 2.87e-41])
F_22: least value j=103 [7.15161000074733337e-33 +/- 1.45e-52], least derivative j=103 [6.37210930224165507e-34 +/- 4.31e-52], first value growth after j>=80: j=103 ratio=[4.42031051072258249 +/- 3.51e-18], first derivative growth after j>=80: j=103 ratio=[4.46322612536023878 +/- 2.41e-18]
  value sums j21..40=[0.000134693747721880453 +/- 3.51e-22] (budget fraction [0.000250228988828346650 +/- 3.72e-22]), j41..80=[1.10118595510691262e-14 +/- 2.18e-32] (budget fraction [2.04574193471355992e-14 +/- 6.66e-33]), j81..120=[4.90371254906031775e-23 +/- 4.10e-41] (budget fraction [9.10993311426665469e-23 +/- 4.11e-41])
  derivative sums j21..40=[2.48931741062047505e-6 +/- 1.22e-24] (budget fraction [0.000267227352867346908 +/- 9.63e-23]), j41..80=[3.95780107786572981e-16 +/- 1.66e-34] (budget fraction [4.24868560634849284e-14 +/- 2.37e-32]), j81..120=[5.07490316995341896e-24 +/- 2.99e-43] (budget fraction [5.44789079278858077e-22 +/- 4.99e-40])
F_23: least value j=103 [3.98900914265024987e-32 +/- 7.27e-51], least derivative j=103 [3.55422095390621751e-33 +/- 3.38e-51], first value growth after j>=80: j=103 ratio=[4.45553211200152908 +/- 2.59e-18], first derivative growth after j>=80: j=103 ratio=[4.49878969722598918 +/- 6.45e-19]
  value sums j21..40=[0.000262589344463092416 +/- 1.43e-22] (budget fraction [0.000487828627931362480 +/- 3.75e-22]), j41..80=[3.13461619777926091e-14 +/- 2.17e-32] (budget fraction [5.82337231535687016e-14 +/- 2.55e-32]), j81..120=[3.09774372235908113e-22 +/- 4.78e-40] (budget fraction [5.75487201451854968e-22 +/- 1.32e-40])
  derivative sums j21..40=[4.85553420391875047e-6 +/- 3.73e-25] (budget fraction [0.000521239897545509047 +/- 4.85e-22]), j41..80=[1.12696025377309577e-15 +/- 1.14e-33] (budget fraction [1.20978788850970019e-13 +/- 7.19e-32]), j81..120=[3.20598971365921075e-23 +/- 3.66e-41] (budget fraction [3.44161873791559000e-21 +/- 1.89e-39])
F_24: least value j=103 [2.14727512952643298e-31 +/- 2.96e-49], least derivative j=103 [1.91322957993488895e-32 +/- 3.18e-50], first value growth after j>=80: j=103 ratio=[4.49075370734479227 +/- 8.55e-19], first derivative growth after j>=80: j=103 ratio=[4.53435324974143468 +/- 2.91e-18]
  value sums j21..40=[0.000501440303600542439 +/- 4.13e-23] (budget fraction [0.000931556974617909435 +/- 2.84e-22]), j41..80=[8.67830153123329798e-14 +/- 2.09e-32] (budget fraction [1.61222228472839426e-13 +/- 3.03e-31]), j81..120=[1.88684846361950440e-21 +/- 1.97e-39] (budget fraction [3.50531625342198492e-21 +/- 3.35e-39])
  derivative sums j21..40=[9.27703715699803406e-6 +/- 6.23e-25] (budget fraction [0.000995886692207194084 +/- 4.68e-22]), j41..80=[3.12098431553720432e-15 +/- 2.32e-33] (budget fraction [3.35036574051693082e-13 +/- 4.85e-31]), j81..120=[1.95284147317557306e-22 +/- 1.35e-41] (budget fraction [2.09636848728022996e-20 +/- 3.78e-38])
```

## Rejected Shortcut

A monotone decreasing or fixed-ratio geometric formal tail from
degree 80 is rejected in the tested window, because all four
index rows grow again from `j=103` to `j=104`.

The fixed-radius Cauchy termwise route is also structurally
insufficient:

```text
A fixed-radius Cauchy coefficient estimate |c_(2j)|<=C*R^(-2j) gives moment terms bounded by C*(a)_j*(u/R^2)^j. For any u/R^2>0, the consecutive ratio (a+j)*(u/R^2) tends to infinity, so this termwise majorant cannot be summed as an infinite absolute tail.
```

This does not close the degree-40 residual budget. It tells us
what the next proof must look like: a genuine asymptotic remainder
estimate for the actual analytic object, not an infinite absolute
sum of the formal Taylor terms.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The degree-240 formal-term scout rejects the naive infinite formal-tail summation route. After the degree-40 residual budget, formal terms continue decreasing through a least-term region near j=103, but every index F_21..F_24 then grows from j=103 to j=104 after j>=80. A fixed-radius Cauchy coefficient bound is also structurally insufficient once multiplied by the moment rising factors. The live route is therefore an actual asymptotic remainder theorem, not an infinite absolute sum of formal Taylor terms.
