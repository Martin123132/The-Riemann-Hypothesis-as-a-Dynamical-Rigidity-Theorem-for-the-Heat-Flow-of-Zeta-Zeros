# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Finite-Part Weighted-Sum Interval Certificate

Date: 2026-07-07

Status: worst-row finite-part weighted-sum interval certificate. This is not a proof
of a quadrature-remainder theorem, a finite-grid interval certificate,
a uniform collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate`.

Proof boundary: this artifact certifies only the finite `n<=30`
weighted Phi/Phi' quadrature sum for `T=10000`, `F_21`, `N=320`.
It does not compose the `n>30` tail source, prove quadrature
remainder, cover all rows, aggregate rounding, or bridge the grid to
a full collar.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval certificate: 6 rows, 0 issues, 320 refined nodes, 320 interval weights, 2 below-one ratios, 0 ready-to-apply rows
```

## Refined Inputs

```text
quadrature order: 320
T: 10000
index: F_21
refinement steps: 120
precision bits: 4096
Phi finite terms: 30
widest refined node width: 3.595897829114163907E-35
maximum relative weight width: 5.764230762018462413243017966627116467124e-34
```

## First-Omitted Comparison

```text
value residual interval: [-7e-34 +/- 5.04e-35]
derivative residual interval: [-1.3806e-32 +/- 8.40e-37]
value scaled abs upper: 6.669422023164208891321813783474733415487e-22
derivative scaled abs upper: 1.380653362177084050322182762738806664470e-24
value ratio upper: 0.9833957992836557769419015895036210773888
derivative ratio upper: 0.9694055674762067320093698741711260875260
both ratios certified below one: True
```

Sample node rows:

```text
root 1: width=0.0000000000000000000000000000000000004566186317413069441729821814890864570575399679580345046157885402493548099300824105739593505859375, rel_weight_width=5.764230762018462413243017966627116467124e-34, Phi/c0=[0.9980938012593711743011320684973196049 +/- 7.52e-38], xPhi'/(2c0)=[-0.00190462765240698040626982056668304800 +/- 4.88e-39]
root 2: width=0.00000000000000000000000000000000000014832761614343335410945522937192146968899117972700300028252896888947276465842151083052158355712890625, rel_weight_width=1.353747411065887054069985401988272801351e-34, Phi/c0=[0.9973642569844503722775592580640961583 +/- 1.84e-38], xPhi'/(2c0)=[-0.00263273829997994072268636054778026935 +/- 1.49e-39]
root 43: width=0.0000000000000000000000000000000000006013459854474487583413526438348406766802673339692265905092292965772315938011161051690578460693359375, rel_weight_width=1.863657069936706284102457198869171719895e-35, Phi/c0=[0.9248549271154701824502856628089643926 +/- 5.61e-38], xPhi'/(2c0)=[-0.07262557169599389666698020241864186291 +/- 4.92e-39]
root 100: width=0.00000000000000000000000000000000000129549225111557898603805534827930446326626431561195444447549818267617638412048108875751495361328125, rel_weight_width=9.043583667091942160666087773689560933891e-36, Phi/c0=[0.7025889714608752691019929003168863863 +/- 2.95e-38], xPhi'/(2c0)=[-0.25355541737907644702050804281315749178 +/- 8.37e-39]
root 200: width=0.0000000000000000000000000000000000029058801026825112742597271063960466694069346116187382519403803460278368220315314829349517822265625, rel_weight_width=5.154064887969916019188198607511703673367e-36, Phi/c0=[0.22796919998577171910525779987721955442 +/- 7.12e-39], xPhi'/(2c0)=[-0.3638740311358604505104443531821494399 +/- 5.36e-38]
root 290: width=0.00000000000000000000000000000000000664409119046345195549758564717349328081883934694704203142816278937488050360116176307201385498046875, rel_weight_width=4.865741127742146292011017009233983195487e-36, Phi/c0=[0.018680982821585348016044217449426163853 +/- 8.83e-40], xPhi'/(2c0)=[-0.08698287041749126596124001343766051303 +/- 3.78e-39]
root 320: width=0.0000000000000000000000000000000000359589782911416390728962280331651698558745318006693101188335492945924443120020441710948944091796875, rel_weight_width=1.803320465338829977162699876611642957881e-35, Phi/c0=[0.001896103534654406590255110301336561003 +/- 7.19e-40], xPhi'/(2c0)=[-0.01463739226993907381988232929811005512 +/- 5.53e-39]
```

Remaining sources:

```text
compose n>30 Phi/Phi'/Phi(0) tail source with this finite-part interval
quadrature-remainder theorem or interval adaptive integration beyond the N=320 sum
rounding and cross-source aggregation
all recorded rows/orders, not just the worst row
finite-grid to full-collar coverage
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The worst-row finite-part weighted-sum interval certificate refines the T=10000, F_21, N=320 Laguerre node brackets to 120 bisection steps, evaluates the finite n<=30 Phi/Phi' contributions with Arb on those node intervals, sums them with certified Christoffel weights, and subtracts the polynomial part through exact Gamma moments. The resulting scaled value and derivative residuals are both certified below one first omitted term. This is still only the worst-row finite-part quadrature sum: n-tail composition, quadrature remainder, rounding aggregation, all-row coverage, and grid-to-collar coverage remain open.
