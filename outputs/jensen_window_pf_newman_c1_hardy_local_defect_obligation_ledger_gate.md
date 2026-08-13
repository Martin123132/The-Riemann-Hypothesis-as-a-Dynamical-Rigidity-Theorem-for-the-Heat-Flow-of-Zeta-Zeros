# Newman C1 Hardy Local-Defect Obligation Ledger Gate

Date: 2026-08-06

Status: `exact_defect_reduction_with_componentwise_open_bounds`; this is an exact reduction and source audit, not a proof of the missing external error theorem.

## What equation (81) contains

Primary source: [2607.15310v1](https://arxiv.org/abs/2607.15310).  The paper says its reduction is exact through equation (69).  Approximation begins when each finite integral is replaced by the contour/saddle formulas (71)--(78).  Equation (81) then packages the resulting `W1` correction: an exact digamma base plus diagonal-saddle, lower-endpoint, and upper-endpoint approximations.

The paper says the remaining endpoint-index terms can in principle be summed over all `n`, but uses `P=3` because they appear to decay quickly.  Pages 41--42 then identify `W1` as the usual largest local-error source and leave the higher-order bound unformulated.

## Paper-to-source map

| paper term | equation | Fortran term | role |
|---|---:|---|---|
| W1 | 81 | `t5` | The two ip-truncated endpoint loops and digamma base are paper W1. |
| W2 | 89 | `t2` | The first six-term complex-erf correction is paper W2. |
| W3_or_W4 | 93_or_94 | `t4` | The sign-selected second complex-erf correction is paper W3 or W4. |
| W5 | 97 | `t1` | The non-saddle digamma difference is paper W5. |
| endpoint_half_sum | 96 | `t3` | t3 is the exact half-sum before floating-point evaluation error. |
| CW_assembly_and_recurrence | 120 | `qq_and_csum` | q assembles CW plus the endpoint term and mgausssum applies it once per level. |

Pinned derivative source: `work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90` (`0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d`).

The executable fixes `ip=3`, while the same `q` routine says `ip~20` was found adequate.  It also uses six-term complex-erf approximations, a truncated `PSI` implementation, saddle Newton tolerance `1e-7`, explicitly omits a `psi(2,z)` correction present in a Maple version, and sets the commented `ecor` correction to zero.  These are obligations, not accusations of numerical failure: each needs a bound before the result can be called rigorous.

## Error ledger

| id | scope | unresolved remainder | required certificate |
|---|---|---|---|
| `eld_w1_saddle_phase` | local_MGS | order-three-and-higher phase on the diagonal saddle contour | A contour-domain majorant retaining every Phi_q term. |
| `eld_w1_vertical_phase` | local_MGS | discarded higher phase plus exponential Taylor remainder on both vertical legs | Two explicit integral remainders with endpoint-uniform constants. |
| `eld_w1_index_tail` | local_MGS | all lower and upper endpoint correction indices omitted beyond ip | A summable majorant uniform in Phi, L, and every admissible branch. |
| `eld_w1_saddle_location` | local_MGS | phase error induced by approximate saddle locations | Interval root isolation and a residual-to-location theorem, including multiplicity separation. |
| `eld_digamma_evaluation` | local_MGS_numeric | finite psi series, asymptotic tail, decimal constants, and pole distance | Directed-rounding digamma balls on every source argument domain. |
| `eld_complex_erf_evaluation` | local_MGS_numeric | piecewise six-term power, Taylor, or asymptotic approximation | Region-by-region analytic remainder bounds and selector margins. |
| `eld_other_W_analytic` | local_MGS | non-W1 endpoint contour linearization and sign-selected special case | Explicit bounds for equations (84)-(94), not only W1. |
| `eld_hierarchy_scheme` | local_MGS_model | truncated transformed coefficients and order-raising decision | An exact transformed-level specification plus a tail bound beyond m1. |
| `eld_kernel_model` | local_MGS_model | difference between the finite source kernel model and the exact transformed child quantity | A formal level adapter proving what finite quantity the denominator-weighted kernel represents. |
| `eld_roundoff_accumulation` | local_MGS_numeric | elementary, phase, summation, and recurrence roundoff | Directed-rounding complex balls with a recorded operation graph. |
| `eld_hardy_representation` | whole_Hardy_evaluator | representation error outside the MGS local defect | Explicit Euler-Maclaurin, lower-cutoff, transition, and Riemann-Siegel-tail bounds. |
| `eld_shared_multi_shift` | whole_shifted_batch | shift dependence suppressed by the shared hierarchy parameterization | Differentiate a fixed branch through order six or enclose its analytic surrogate on a disk. |
| `eld_transition_set` | whole_shifted_batch | branch jumps or selector changes across the physical interval/disk | Strict interval margins or an explicit partition with jump-mode vectors. |

## Exact-defect route

Equation (120) already provides a second route that avoids proving every internal cancellation separately:

```text
epsilon_l = S_l - [a_l S_(l+1) + q_l].

e_l = epsilon_l + a_l e_(l+1),

e_0 = sum_k (product_(j<k) a_j) epsilon_k.
```

Enclose the two finite level sums and the source recurrence term directly. This captures W1, W2-W5, coefficient-transform, and numerical defects without requiring their cancellations to be bounded separately.

There is also an algorithmic bypass.  If level `m` is evaluated directly, recurrence begins from that exact level and every descendant defect `epsilon_m,epsilon_(m+1),...` disappears from the root-error identity.  In particular, directly summing the first parent above the kernel removes the very defect the paper says is usually largest.  The physical chain lengths must be measured before claiming this is cheap.

The immediate caveat is structural: The source's denominator-weighted kernel and raised-order coefficients must first be identified with a precise finite target; otherwise a direct defect could compare the wrong two quantities.

## Next target

Add read-only chain telemetry to one accepted low-height fixture. For every level record L, m1, Phi, x, fracL, ip, q components, direct finite level values, and the exact local defect at increasing precision. Then test direct evaluation of the first parent shell. Do not alter the accepted evaluator output path during this scout.

That scout should compare three routes on identical saved chains: termwise `W1` bounds, direct finite local-defect enclosure, and one-parent exact-shell replacement.  The winner is the route that yields a useful uniform constant with the least new machinery, not the route that merely matches samples best.

This gate proves the paper-to-source component map, the exact affine defect accumulation identity, and the exact-shell bypass implication. It does not prove any W1 tail or contour bound, special-function enclosure, transformed-level adapter, evaluator C6 or disk theorem, transition exclusion, outer Hardy representation bound, physical carrier value, interval kernel approximation, determinant sign or bound, complete-current inequality, all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

```text
built Newman C1 Hardy local-defect obligation ledger gate: 16 rows, 6 paper-to-source mappings, 13 error components, 2 exact defect identities, 1 exact-shell bypass, 2 conditional derivative routes and 0 external error bounds
```
