# Rigorous Numerics Compute Plan

Date: 2026-07-03

Status: compute plan. This is not a proof artifact; it prioritizes reproducible numerical work and certificate upgrades.

Context:

The machine is now much stronger: 32 GB RAM and roughly 4 TB storage across faster drives. This changes what is practical for the computational side of the RH dynamical-rigidity programme.

The purpose of the compute work is not to "prove RH numerically." The purpose is to:

1. reproduce existing finite claims cleanly;
2. stress-test the signed Hankel/Jensen patterns;
3. find counterexamples quickly if a proposed bridge theorem is false;
4. generate rigorous finite certificates where possible.

## What The Hardware Upgrade Buys

### More RAM

Useful for:

- higher-degree polynomial audits;
- wider lambda grids;
- more simultaneous zero windows;
- high-precision coefficient caches;
- interval or ball arithmetic experiments;
- multiprocessing without constant swapping.

### Faster And Larger Storage

Useful for:

- storing JSONL result logs for every tested degree, shift, lambda, and precision;
- checkpointing long computations;
- keeping coefficient caches rather than recomputing them;
- preserving failed cases for later inspection.

## Compute Philosophy

Every run should produce durable artifacts:

- machine-readable result logs;
- summarized CSV tables;
- exact script version/hash;
- precision settings;
- failure examples, not just pass counts.

No result should exist only as notebook output.

## Tier 0: Machine And Environment Baseline

Deliverable:

```text
work/rh_compute/environment_report.txt
```

Tasks:

- record Python version;
- record installed numeric libraries;
- detect CPU core count;
- detect available memory;
- check whether `mpmath`, `numpy`, `scipy`, `sympy`, and ideally `python-flint` are available;
- run a small benchmark for polynomial coefficient generation and Sturm/root checks.

Success criterion:

We know what tools are actually installed before designing expensive runs.

## Tier 1: Clean Reproduction

Deliverable:

```text
work/rh_compute/results/reproduction_summary.csv
```

Tasks:

- reproduce the existing Jensen/Sturm pass counts from the bundle;
- reproduce the signed Hankel determinant checks already claimed;
- run from clean scripts, not notebooks;
- write every test row to JSONL.

Initial targets from the bundle:

```text
degrees 21..30, n = 0..40
degrees 31..40, n = 0..20
lambda = 0 and small positive lambda values already tested
```

Success criterion:

The current claimed finite evidence can be regenerated on demand.

## Tier 2: Extended Search

Deliverable:

```text
work/rh_compute/results/extended_search_summary.csv
```

Suggested first extension:

```text
Jensen degrees: 21..60
shifts n: 0..200
lambda grid: 0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5
precision: start at 80 decimal digits, rerun suspicious cases at 160+
```

Signed Hankel extension:

```text
matrix sizes m: 1..16
shifts s: 0..50
lambda grid: same as above
```

Success criterion:

Either:

- the signed pattern survives much broader tests, or
- we find the first concrete failure and use it to refine or abandon the bridge route.

## Tier 3: Rigorous Finite Certification

Deliverable:

```text
work/rh_compute/certificates/
```

Candidate methods:

- exact rational enclosure of coefficients when formulas allow it;
- interval arithmetic for coefficient bounds;
- Sturm sequence sign certification with directed rounding;
- Arb/FLINT ball arithmetic via `python-flint` if available;
- independent double-run verification at different precision levels.

Current progress:

```text
work/rh_compute/scripts/coefficient_tail_bounds.py
```

now gives analytic bounds for:

```text
omitted n-series tail beyond n_sum;
omitted u-integral tail beyond cutoff.
```

For the main reproduced Hankel/Toeplitz cache (`n_sum=100`, `cutoff=8`, `k<=32`, lambda up to `0.1`), the worst omitted `A_k` tail is about `1.91e-13597`. For the fragile Jensen cache (`n_sum=120`, `cutoff=9`, `k<=80`), the worst omitted `A_k` tail is about `1.53e-19597`.

Remaining coefficient-enclosure problem:

```text
validated quadrature for the retained finite sum
n = 1..n_sum, u in [0, cutoff].
```

Update:

```text
work/rh_compute/scripts/acb_coefficient_enclosures.py
```

now uses `python-flint`'s rigorous `acb.integral` for selected retained finite integrals, then adds the analytic tail bounds. Current certified coefficient coverage is:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}, k = 0..32
n_sum = 100, cutoff = 8
```

These enclosures justify the cache balls `A_k(cache) +/- 1e-80` for selected structural Toeplitz/PF finite certificates.

Success criterion:

For selected finite claims, produce certificates stronger than floating-point "all roots looked real."

## Tier 4: Counterexample Search

Deliverable:

```text
work/rh_compute/results/counterexample_search.md
```

Tasks:

- build toy entire functions satisfying the same local heat-flow and zero-motion identities;
- test proposed bridge lemmas against these toys;
- preserve minimal failing examples.

Success criterion:

False proof routes fail cheaply in computation before they enter the manuscript.

## Priority Order

1. Build `formal_core.md` first, so the computations know what theorem they are supporting.
2. Build a clean extraction script for Jensen/Hankel definitions.
3. Run Tier 0 environment baseline.
4. Run Tier 1 reproduction.
5. Only then expand to Tier 2.

## Guardrails

- Do not treat broader pass counts as proof.
- Do not bury failures.
- Do not mix physical Newman lambda, forward heat time, and reduced flow parameter.
- Do not run huge jobs until the output format is stable.
- Do not rely on notebooks as the primary record.

## Near-Term Compute Sprint

The first useful compute sprint is:

```text
Sprint C1:
  1. Extract Jensen/Hankel code from notebooks.
  2. Create a clean command-line runner.
  3. Run environment baseline.
  4. Reproduce existing pass counts.
  5. Save JSONL + CSV outputs.
```

After that, the stronger laptop can carry the extended search without turning every run into a waiting-room saga.
