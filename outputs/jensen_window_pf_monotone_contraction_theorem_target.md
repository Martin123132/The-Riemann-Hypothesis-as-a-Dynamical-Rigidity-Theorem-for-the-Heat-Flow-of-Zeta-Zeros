# Jensen-Window PF Monotone-Contraction Theorem Target

Date: 2026-07-06

Status: open theorem target. This is not a proof of an analytic
monotone-contraction theorem, Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_monotone_contraction_theorem_target`.

Proof boundary: this artifact states a missing analytic theorem target. It
does not prove the theorem, does not construct a Cauchy-Binet kernel or
determinant integral, does not prove `jwpf_06`, and does not prove
`Lambda <= 0`.

## 2026-07-10 Scope Correction

The interval certificate
`outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md`
proves `x_121<x_120` for the actual zeta kernel at `lambda=-1156`.
Consequently, adjacent-`k` monotonicity is not a theorem for every real
`lambda` and every `k`. The live negative-lambda target is now parameter
specific: find one suitable moderate entry parameter, certify a finite
collar there, and prove its eventual-`k` tail. Positive-lambda finite stress
and bounded/frontier applications remain separate diagnostics.

The live negative-lambda instance is now fully localized at `lambda=-100`.
The finite collar reaches `k=318`, and the full-kernel tail differs from its
first-summand wall by at most `16/(k-1)^6`. The remaining theorem is exactly

```text
L_k^(1)>=1/(4*k^2), k>=319.
```

See `outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md`.
The exact sufficient reduction is now

```text
kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318,
```

as recorded in
`outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md`.
The leading term and cubic and fifth-order corrections are certified in
`outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md`;
`outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md`
proves the live seventh-order normalized remainder floor `-79/1000` on
`0.9264<=u<=5`, while
`outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md`
proves `H_t>=-3/250`
on `u>=5`. The lambda=-100 adjacent wall is therefore closed; this broader
monotone-contraction theorem remains open outside that parameter-specific result.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_theorem_target.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_theorem_target.py
```

Current result:

```text
validated Jensen-window PF monotone contraction theorem target: 9 rows, 0 issues, 0 ready-to-apply rows, 2 live routes
```

## Statement

The coefficient normalization used throughout the finite Arb corpus is:

```text
mu_{2k}(lambda) = 2 * integral_0^infty u^(2k) exp(lambda*u^2) Phi(u) du
A_k(lambda) = mu_{2k}(lambda) * k! / (2k)!
```

Define adjacent ratio contractions:

```text
x_k = (A_{k+1}/A_k) / (A_k/A_{k-1})
```

The original theorem target was:

```text
0 < x_k <= 1
x_{k+1} >= x_k
```

for every required `k`, shift, and lambda regime used by the Jensen-window PF
bridge. The `lambda=-1156`, `k=120` counterexample rules out this universal
form. On any surviving parameter-specific regime, the equivalent inequality
is:

```text
A_{k+2}*A_k^3 >= A_{k+1}^3*A_{k-1}
Delta^3 log A_{k-1} >= 0
```

This is the increasing-contraction direction. The older ratio-condition scout
rejected the opposite direction `x1 >= x2 >= x3` by exact countermodel.

## Moment Normalization

The factorial normalization cannot be dropped. With raw even moments
`M_k = mu_{2k}(lambda)`,

```text
x_k(A) = ((2*k-1)/(2*k+1)) * (M_{k+1}*M_{k-1}/M_k^2)
```

and the increasing-contraction target is equivalent to:

```text
M_{k+2}*M_k^3/(M_{k+1}^3*M_{k-1}) >= ((2*k-1)*(2*k+3))/(2*k+1)^2
```

So a proof that stops at ordinary Stieltjes moment log-convexity is too weak
and points at the wrong object.

## Contract

A proof of this target must:

```text
1. prove A_k(lambda)>0 from the Phi/Newman-kernel integral;
2. prove 0 < x_k <= 1 and x_{k+1} >= x_k uniformly in the required regime;
3. keep A_k=mu_{2k} k!/(2k)! explicit;
4. use a noncircular analytic mechanism such as a total-positive kernel,
   determinant identity, or closed heat-flow differential inequality;
5. explain why the rational log-concavity countermodel is outside the
   zeta-specific hypotheses;
6. state whether it covers bounded/frontier column rows, all column rows, or
   all Toeplitz/Jacobi-Trudi shapes; even the full static ratio cone is blocked
   as an all-`m` claim by the exact degree-7 `m=11` counterexample;
7. avoid endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and
   Lambda <= 0 as input assumptions.
```

## Target Rows

```text
mct_01_exact_monotone_contraction_statement:
  open statement; prove A_{k+2}*A_k^3 >= A_{k+1}^3*A_{k-1}

mct_02_factorial_normalized_moment_reduction:
  exact reduction; raw moment log-convexity alone is insufficient

mct_03_total_positivity_kernel_route:
  live route if a zeta-specific TP kernel or determinant identity is built

mct_04_heat_flow_differential_inequality_route:
  live route if a closed noncircular lambda-flow inequality is derived

mct_05_plain_stieltjes_moment_logconvexity:
  insufficient route

mct_06_endpoint_pf_or_laguerre_polya_route:
  circular as an input

mct_07_generic_ratio_conditions:
  rejected by countermodel and constructed-extension scouts, including the
  exact all-m monotone-contraction counterexample

mct_08_finite_stress_support:
  finite evidence only: 2875/2875 checked rows

mct_09_column_frontier_application:
  conditional bounded/frontier column application, not all-m closure
```

## Evidence And Boundaries

The finite evidence currently attached to this target is:

```text
outputs/jensen_window_pf_monotone_contraction_stress.md
validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues
```

The exact frontier application is:

```text
outputs/jensen_window_pf_monotone_contraction_frontier_scout.md
validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues
```

The exact all-m promotion gate is:

```text
outputs/jensen_window_pf_monotone_contraction_all_m_counterexample.md
validated Jensen-window PF monotone-contraction all-m counterexample: degree 7, m=11, exact full-cone witness, 6 lower walls, negative normalized value, 0 issues
```

It proves that the full propagated static ratio cone does not imply all-`m`
column-recurrence positivity. Any all-`m` claim must add heat-flow or
Xi/Phi-specific structure beyond that cone, or use a different positive
construction.

The remaining heat-flow cone-entry problem is now recorded as:

```text
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
validated Jensen-window PF heat-flow cone-entry asymptotic target: 8 rows, 0 issues, 0 ready-to-apply rows, 1 live routes
```

It separates a plausible fixed-k large-negative-lambda asymptotic route from
the harder full infinite or collared finite cone-entry theorem. The heat-flow
route is also sharpened by the exact conditional ratio-cone invariance lemma:

```text
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues
```

It proves inward-pointing boundary algebra for the cone
`(2*k-1)/(2*k+1) <= x_k <= 1`, `x_{k+1}>=x_k`, conditional on a full infinite
or collared finite cone. The heat-flow route also uses the exact
boundary-threshold lemma:

```text
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
```

and by the finite closure scout:

```text
outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md
validated Jensen-window PF heat-flow monotone closure scout: 4 exact rows, 315 threshold rows, 305 flow-bracket rows, 0 issues
```

The threshold subtarget `q >= (2*k-1)/(2*k+5)` is discharged by raw moment
log-convexity. The target remains open because none of these artifacts proves
that the actual zeta heat-flow coefficient sequence enters the full cone at a
suitable lambda, proves the infinite/collared flow legitimacy, or proves the
analytic monotone-contraction theorem for all required coefficients. Even if
that theorem is proved, the all-m counterexample above prevents promoting the
monotone-contraction cone alone into all-column recurrence positivity.

## Integration

This target sharpens:

```text
outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_ratio_condition_scout.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Passing the checker means the theorem target is stated with exact
normalization, finite support, live routes, rejected shortcuts, and proof
boundaries. It does not close the Newman-direction bridge.
