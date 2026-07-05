# Signed-Hankel Certificate Status

Date: 2026-07-04

Status: finite interval certificate for the reproduced signed-Hankel grid. This is not a proof of RH, not an all-order signed-Hankel theorem, and not a proof of `Lambda <= 0`.

## Certified Object

For the exponential/Jensen coefficient sequence:

```text
A_k(lambda) = mu_{2k}(lambda) k! / (2k)!
```

the tested signed-Hankel determinant is:

```text
(-1)^(m(m+1)/2) det(A_{i+j+s}(lambda))_{i,j=0}^m.
```

## Enclosure-Backed Arb Probe

The certificate uses the rigorous coefficient enclosure rows:

```text
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k33_k48.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k49_k64.jsonl
```

These provide `A_ball` enclosures through:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
k = 0..64
```

The signed-Hankel range:

```text
m = 0..19
s = 0..24
```

uses coefficients only up to:

```text
2*m + s <= 62.
```

## Result

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
```

the enclosure-backed Arb result is:

```text
rows = 500
positive = 500
failed_or_inconclusive = 0
```

Across the grid:

```text
2,500 / 2,500 signed-Hankel determinants are interval sign-separated positive.
```

The checker:

```text
python work/rh_compute/scripts/check_hankel_certificate_manifest.py
```

validates:

```text
2,500 signed-Hankel finite certificates
```

## Why This Changed The Status

A previous uniform-radius propagation using cached decimal centers and:

```text
cache_A +/- 1e-80
```

gave:

```text
lambda = 0
m = 0..12
s = 0..8
108 / 117 sign-separated
9 fragile high-m/high-shift determinants inconclusive
```

That was a propagation-width limitation, not a sign failure. The cache centers are truncated enough that a much smaller uniform cache radius is not justified from the existing cache alone.

The new script:

```text
work/rh_compute/scripts/arb_hankel_enclosure_sign_probe.py
```

uses the actual rigorous `A_ball` fields from the coefficient-enclosure logs instead of wrapping truncated cache centers. With those sharper justified balls, all 2,500 currently promoted signed-Hankel determinants are sign-separated.

## Current Frontier

The larger rectangle:

```text
m = 0..20
s = 0..24
needed_max_k = 64
```

was tested at `dps = 520`, but is not promoted because the far corner contains interval-width inconclusives:

```text
lambda = 0:     m = 20, s = 23 and s = 24 inconclusive
lambda = 1e-6:  m = 20, s = 24 inconclusive
lambda = 1e-4:  m = 20, s = 23 and s = 24 inconclusive
lambda = 1e-2:  m = 20, s = 23 and s = 24 inconclusive
lambda = 1e-1:  m = 20, s = 23 and s = 24 inconclusive
```

A targeted `dps = 1000` rerun of the lambda-zero corner stayed inconclusive, so this is currently an enclosure-width frontier rather than a sign failure.

## Boundary

This proves only the finite signed-Hankel determinant inequalities in the stated grid.

It does not prove:

```text
all signed-Hankel determinants have the pattern;
the signed-Hankel pattern implies Jensen hyperbolicity;
the signed-Hankel pattern implies Laguerre-Polya membership;
RH;
Lambda <= 0.
```

The remaining missing theorem is still an all-order bridge from signed-Hankel/sign-regularity data to Jensen hyperbolicity, Laguerre-Polya membership, or a Newman-boundary obstruction.
