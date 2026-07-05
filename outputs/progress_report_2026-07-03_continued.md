# Progress Report Continued

Date: 2026-07-03

Status: sprint progress report. This is not a proof artifact; it records continued implementation and diagnostic work.

Goal:

Continue extending the RH dynamical-rigidity corpus toward a referee-legible proof programme for `Lambda <= 0`, prioritizing signed Hankel / Jensen / sign-regularity evidence and rigorous guardrails.

## 1. Runner Cache Reuse Added

Updated:

```text
work/rh_compute/scripts/jensen_hankel_runner.py
```

New capabilities:

```text
--coeff-cache PATH
--require-coeff-cache
```

Purpose:

Reuse previously computed `A_k(lambda)` coefficient caches, avoiding repeated high-precision quadrature for Jensen or Hankel reruns.

Verification:

```text
smoke_jensen_cache_reuse:
  lambda = 0
  d = 16
  n = 54
  source = repro_jensen_fragile_d16_n54_coefficients.jsonl
  result = 1/1 passed
  elapsed = about 0.6 seconds
```

Cache reuse also worked for:

```text
smoke_jensen_grid_cache_reuse:
  lambda = 0
  d = 20
  n = 60
  source = repro_jensen_fragile_d16_d20_n50_n60_coefficients.jsonl
  result = 1/1 passed
  elapsed = about 1.6 seconds
```

## 2. Local FLINT/Arb Installed

Installed locally:

```text
work/rh_compute/vendor/python-flint 0.8.0
```

Verified local imports:

```text
flint.arb
flint.acb
flint.arb_mat
```

This does not alter the global Python environment.

## 3. Arb Hankel Sign Probe Added

Added:

```text
work/rh_compute/scripts/arb_hankel_sign_probe.py
```

Purpose:

Given cached decimal coefficients and a chosen absolute coefficient radius, use Arb balls to test whether signed Hankel determinant balls remain separated from zero.

Smoke probe:

```text
coeff cache = repro_hankel_15c_coefficients.jsonl
lambda = 0
m = 0..2
shift = 0..1
absolute coefficient radius = 1e-80
dps = 120
```

Result:

```text
6/6 sign-separated and positive after signature correction
```

Artifact:

```text
work/rh_compute/results/arb_probe_hankel_low_order.jsonl
```

Important caveat:

This is coefficient-ball propagation only. It becomes a finite certificate only after the coefficient balls are proved to enclose the exact transcendental moments.

2026-07-04 update:

The signed-Hankel propagation has now been rerun directly from rigorous `A_ball` coefficient enclosure rows rather than from truncated cache centers with a uniform radius. The resulting manifest:

```text
python work/rh_compute/scripts/check_hankel_certificate_manifest.py
```

validates:

```text
2,500 signed-Hankel finite certificates
```

See:

```text
outputs/signed_hankel_certificate_status.md
```

## 4. Targeted Jensen Fragile Region Expanded

Run:

```text
repro_jensen_fragile_d16_d20_n50_n60
```

Parameters:

```text
lambda = 0
d = 16..20
n = 50..60
dps = 100
rational_digits = 70
n_sum = 120
cutoff = 9
```

Result:

```text
55/55 passed
failures: 0
```

By degree:

```text
d=16: 11/11, positive root count 16
d=17: 11/11, positive root count 17
d=18: 11/11, positive root count 18
d=19: 11/11, positive root count 19
d=20: 11/11, positive root count 20
```

Artifacts:

```text
work/rh_compute/results/repro_jensen_fragile_d16_d20_n50_n60.jsonl
work/rh_compute/results/repro_jensen_fragile_d16_d20_n50_n60.csv
work/rh_compute/results/repro_jensen_fragile_d16_d20_n50_n60_summary.json
work/rh_compute/results/repro_jensen_fragile_d16_d20_n50_n60_coefficients.jsonl
```

## 5. Certification And Countermodel Gates Added

Added:

```text
outputs/certification_and_countermodel_gates.md
```

This defines:

- numerical evidence levels;
- circularity checks;
- local-repulsion countermodel gate;
- signed-Hankel claim gate;
- theorem-search fit/misfit gate;
- safe result language.

The main enforced rule:

```text
No local-repulsion-only argument and no RH-at-lambda=0 assumption may enter
as a proof step toward Lambda <= 0.
```

## Current Position

Stronger now:

- Reproducible Jensen fragile-region evidence expanded from one case to 55 cases.
- Coefficient cache reuse makes future Jensen sweeps practical.
- Local FLINT/Arb is available for certificate development.
- A first Arb determinant sign-probe exists.
- Countermodel and certification gates are explicit.

Still missing:

- rigorous moment enclosures;
- ball-certified Jensen root counts;
- full Jensen reproduction grid;
- theorem connecting signed Hankel structure to all Jensen hyperbolicity, Laguerre-Polya membership, or no positive Newman boundary.

## Recommended Next Work

1. Implement rigorous or semi-rigorous coefficient enclosures for `A_k(lambda)`:

```text
Phi series truncation bound
tail integral bound
quadrature enclosure
```

2. Extract Groechenig's Schoenberg/zeta RH-equivalent condition into a file:

```text
outputs/schoenberg_zeta_condition.md
```

3. Use the Arb probe on increasingly fragile determinants with progressively smaller coefficient radii to estimate what coefficient precision is needed for Level 4 certificates.

4. Add a cache-driven Jensen batch mode for cheap reruns over arbitrary `d,n` windows once coefficients are available.
