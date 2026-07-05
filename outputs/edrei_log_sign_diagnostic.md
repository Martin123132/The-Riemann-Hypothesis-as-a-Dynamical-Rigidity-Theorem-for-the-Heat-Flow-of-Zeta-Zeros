# Edrei Log-Sign Diagnostic

Date: 2026-07-04

Status: finite necessary-condition diagnostic for the coefficient PF route. This is not a proof of PF-infinity, Laguerre-Polya membership, RH, or `Lambda <= 0`.

## Object

Normalize the ordinary Taylor coefficient sequence:

```text
c_k(lambda) = A_k(lambda) / k!
d_k(lambda) = c_k(lambda) / c_0(lambda)
H_lambda(z) = sum_{k >= 0} d_k(lambda) z^k
```

For an entire restricted Laguerre-Polya / coefficient-PF target with only nonpositive real zeros, one expects:

```text
H(z) = exp(gamma z) product_j (1 + beta_j z)
gamma >= 0
beta_j >= 0
```

Therefore:

```text
log H(z) = sum_{n >= 1} ell_n z^n
q_n = n ell_n
(-1)^(n-1) q_n >= 0
```

Failure of this sign alternation would be a serious obstruction to the entire coefficient-PF route.

## Arb Enclosure Probe

Scripts:

```text
work/rh_compute/scripts/arb_edrei_log_sign_probe.py
work/rh_compute/scripts/check_edrei_log_sign_manifest.py
```

Inputs:

```text
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k33_k48.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k49_k64.jsonl
```

These provide rigorous `c_ball` enclosures for:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
k = 0..64
```

The `k = 49..64` enclosure run used:

```text
dps = 150
abs_tol = 1e-105
rows = 80
max c_ball radius about 1.33e-259
elapsed = 78.74 s
```

## Result

For each lambda:

```text
max_n = 64
dps = 340
rows = 64
failed_or_inconclusive = 0
all_ok = true
```

Manifest validation:

```text
python work/rh_compute/scripts/check_edrei_log_sign_manifest.py
validated 320 finite Edrei-log sign diagnostics
```

Example first terms for `lambda = 0`:

```text
q_1 > 0:  5.7762482788547426972334526075847535e-3
q_2 < 0: -2.3232874553293553853041414044837862e-6
q_3 > 0:  2.2527176781402074952403368139815767e-9
q_4 < 0: -2.5899675009882455854424690701615518e-12
q_5 > 0:  3.1383438970865246250997232420260304e-15
```

All signs above are interval-separated by Arb balls in the row logs.

## Interpretation

This removes one finite obstruction to the coefficient PF route:

```text
c_k = A_k/k!
```

It does not prove that all zeros of `H_lambda` are real and nonpositive, and it does not prove all Toeplitz minors are nonnegative. It is a necessary-condition check that sits between:

```text
finite Toeplitz/PF certificates
```

and the missing all-order theorem:

```text
prove H_0 is restricted Laguerre-Polya
or prove c_k(0) is PF-infinity directly
```

The next useful upgrade would be either:

```text
1. extend the rigorous coefficient enclosures beyond k = 64 and rerun this log-sign probe;
2. derive a representation proving (-1)^(n-1) q_n >= 0 for all n;
3. show that the q_n are actual power sums of a positive zero-parameter measure.
```

## Power-Hankel Necessary Condition

For a representation:

```text
H(z) = exp(gamma z) product_j (1 + beta_j z)
gamma >= 0
beta_j >= 0
```

the signed logarithmic coefficients:

```text
p_n = (-1)^(n-1) q_n
```

should behave like positive power sums. Therefore shifted Hankel determinants:

```text
det(p_{i+j+s})_{i,j=0}^m
```

with `s >= 1` provide another finite necessary-condition diagnostic.

Scripts:

```text
work/rh_compute/scripts/arb_edrei_power_hankel_probe.py
work/rh_compute/scripts/check_edrei_power_hankel_manifest.py
work/rh_compute/scripts/edrei_power_hankel_midpoint_scout.py
work/rh_compute/scripts/check_edrei_power_hankel_frontier_manifest.py
work/rh_compute/scripts/check_edrei_power_hankel_boundary_manifest.py
```

Promoted finite diagnostic:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
m = 0,    s = 1..57
m = 1,    s = 1..55
m = 2,    s = 1..53
m = 3,    s = 1..51
m = 4,    s = 1..49
m = 5,    s = 1..47
m = 6,    s = 1..45
m = 7,    s = 1..43
m = 8,    s = 1..41
m = 9,    s = 1..39
m = 10,   s = 1..37
m = 11,   s = 1..35
m = 12,   s = 1..33
m = 13,   s = 1..31
m = 14,   s = 1..29
m = 15,   s = 1..27
m = 16,   s = 1..25
m = 17,   s = 1..23
m = 18,   s = 1..21
m = 19,   s = 1..19
m = 20,   s = 1..17
m = 21,   s = 1..15
m = 22,   s = 1..13
m = 23,   s = 1..11
m = 24,   s = 1..9
m = 25,   s = 1..7
m = 26,   s = 1..5
m = 27,   s = 1..3
m = 28,   s = 1
validated 4,205 finite Edrei power-Hankel diagnostics
```

Larger scout:

```text
m = 0..8
s = 1..48
needed_max_n = 64
```

found no negative determinants, but has interval-width inconclusives at high shifts. The first inconclusives occur at:

```text
m = 2, s = 47
m = 4, s = 46
m = 5, s = 41
m = 6, s = 37
m = 7, s = 33
m = 8, s = 30
```

The `m = 20, s = 1` edge was interval-width inconclusive from the original
`dps = 340` Edrei rows. A midpoint-only scout was mixed across the lambda
grid, so it was not trusted. Recomputing tighter `k = 0..49` coefficient
enclosures (`dps = 180`, `abs_tol = 1e-120`) and Edrei-log rows through
`n = 49` at `dps = 1000` moves the promoted edge to `m = 24, s = 1` and
fills the lower high-shift wedge through `2m+s <= 49` for `m <= 8`.
The next uniform frontiers use tight coefficient enclosures and promote all
five lambdas through the completed `2m+s <= 57` staircase in the table above.
The `n = 51` layer first needed a `lambda = 1e-6`, `k = 0..51`,
`dps = 220`, `abs_tol = 1e-140` repair. The `n = 53` and later layers needed
the same tighter coefficient treatment through `k = 57` for all five lambdas,
followed by `dps = 2400` Edrei-log reruns.

```text
m = 24, s = 3: [1.985758669646888065740e-2728 +/- 8.57e-2750]
m = 25, s = 1: [3.67828662510165430155e-2732 +/- 5.47e-2753]
```

Boundary history and repair are now checked explicitly:

```text
validated 2 retired inconclusive blocker rows and 3 repaired positive boundary rows
```

The midpoint scout remains recorded as a non-rigorous warning artifact:

```text
validated 5 non-rigorous Edrei midpoint frontier scouts
```

It is a guardrail against promoting center-only signs.

Higher shifts beyond this staircase and the broader high-shift scout remain
diagnostic frontiers, not sign failures.

## Moment-Quadrature Reconstruction Scout

An Arb recurrence reconstruction scout has been added:

```text
work/rh_compute/scripts/edrei_moment_quadrature_scout.py
work/rh_compute/scripts/check_edrei_quadrature_scout_manifest.py
outputs/edrei_moment_quadrature_scout.md
```

It treats the shifted sequence:

```text
a_n = p_{n+1}
```

as candidate moments of `x dnu(x)` and checks finite monic recurrence data using Arb interval arithmetic.

Current result:

```text
orders 2..12: 55/55 Arb recurrence rows positive
orders 2..20: frontier scout with 55 positive, 0 negative, 40 inconclusive rows
```

This remains a finite diagnostic and does not prove the Edrei representation. The frontier beyond order 12 is inconclusive, not negative.
