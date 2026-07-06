# Jensen-Window PF Theorem Machinery Fit Matrix

Date: 2026-07-06

Status: theorem-search fit matrix. This is not a proof of Jensen-window
PF-infinity, Jensen hyperbolicity, Laguerre-Polya membership, RH, or
`Lambda <= 0`; it audits which total-positivity and zero-preserver theorem
families could or could not close `jwpf_06_sign_regular_to_jensen_pf_conversion`.

## Purpose

The obligation ledger identifies the central open bridge:

```text
jwpf_06_sign_regular_to_jensen_pf_conversion
```

This matrix asks a sharper question: which known theorem machinery can turn
proved sign-regular information about the actual `A_k(0)` into total
nonnegativity of every binomially weighted Jensen-window Toeplitz matrix?

Machine-readable matrix:

```text
work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
```

Current result:

```text
validated Jensen-window PF theorem machinery fit matrix: 7 rows, 0 issues, 0 ready-to-apply rows
```

## Required Bridge Features

A theorem capable of closing `jwpf_06` must:

```text
output every binomially weighted Jensen-window Toeplitz matrix for all d,n
handle binomial weights binom(d,j)
handle all shifts n uniformly
start from proved noncircular hypotheses about the actual A_k(0)
avoid assuming Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0
```

## Source-Anchored Rows

```text
tm_01_asw_edrei_pf_sequence_characterization:
  endpoint_equivalence_only
  ASW/Edrei tells us what PF-infinity means for a sequence already known to be
  totally positive. It does not convert signed-Hankel data into Jensen-window
  PF-infinity.

tm_02_schoenberg_pf_functions_variation_diminishing:
  possible_if_new_kernel_representation
  promising only if we construct a PF kernel or variation-diminishing operator
  that produces the actual Jensen-window Toeplitz minors.

tm_03_karlin_basic_composition_cauchy_binet:
  possible_if_new_kernel_representation
  promising only if each Jensen-window Toeplitz minor becomes a Cauchy-Binet
  sum/integral of nonnegative determinant products.

tm_04_polya_schur_multiplier_preservers:
  possible_if_preserver_hypotheses_are_proved
  relevant to binomial weights, but it still needs an input family of windows
  already proved to satisfy the appropriate real-root/PF hypothesis.

tm_05_gantmacher_krein_sign_regular_matrices:
  possible_if_preserver_hypotheses_are_proved
  closest in spirit to signed-Hankel evidence, but no transfer theorem from the
  signed/indefinite Hankel pattern to Jensen-window Toeplitz TN is identified.

tm_06_laguerre_polya_jensen_limit:
  conditional_downstream_only
  useful only after all-degree/all-shift Jensen hyperbolicity is proved.

tm_07_finite_grid_or_rh_assuming_shortcuts:
  rejected_circular
  finite grids, local repulsion, or RH-assuming arguments are proof-safety traps.
```

## Primary Source Anchors

The matrix records source anchors for the theorem families, including:

```text
https://www.pnas.org/doi/10.1073/pnas.37.5.303
https://link.springer.com/article/10.1007/BF02786970
https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/proof-of-a-conjecture-of-schoenberg-on-the-generating-function-of-a-totally-positive-sequence/7CC0A2ADBCB2ADF112F737EDE229766C
https://link.springer.com/article/10.1007/BF02790092
https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf
```

These are used as theorem-search anchors, not as claims that their hypotheses
are already satisfied by the zeta heat-flow coefficients.

## Current Best Route

The best current route is still structural:

```text
construct a positive kernel, planar network, production matrix, determinant
integral, or sign-regular-to-Toeplitz transfer theorem that produces every
Jensen-window Toeplitz minor for B^{d,n,0}_j = binom(d,j) A_{n+j}(0)
```

The next proof-search action is therefore not another finite grid by itself.
It is a symbolic/theorem search for a uniform identity that survives the
degree-3 countermodel and handles the binomial weights.

## Boundary

Passing this checker means the theorem-search map does not overstate any known
machinery. It records no `ready_to_apply` row, and it keeps endpoint
equivalences, possible structural routes, downstream limiting arguments, and
rejected shortcuts separated.
