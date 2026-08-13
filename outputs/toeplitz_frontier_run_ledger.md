# Toeplitz Frontier Run Ledger

Date: 2026-07-04

Status: operational ledger for finite Toeplitz/PF frontier work. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.

## Certified Frontier

The promoted certificate manifest currently validates:

```text
95 promoted positive certificate summaries
1 beta negative control
```

The main certified frontiers are:

```text
lambda grid {0, 1e-6, 1e-4, 1e-2, 1e-1}, N <= 30, orders <= 4
lambda grid {0, 1e-6, 1e-4, 1e-2, 1e-1}, N <= 22, orders <= 5
lambda grid {0, 1e-6, 1e-4, 1e-2, 1e-1}, N = 24, orders <= 5
```

These are Level-4 finite certificates because the corresponding coefficient cache balls are covered by the retained-integral/tail enclosure workflow and the Arb determinant intervals are sign-separated.

## Completed Exact Grid Promoted

The full higher-order frontier grid has completed and passed the manifest gate:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 24
orders <= 5
sequence = taylor
abs_radius_on_A = 1e-80
dps = 150
structural_mode = fast
enumeration_mode = nonzero
tests = 1,923,675,220
evaluated_nonstructural_tests = 404,448,400
structural_zero = 1,519,226,820
positive = 404,448,400
negative = 0
inconclusive_contains_zero = 0
unknown = 0
zero = 0
recorded_problem_rows = 0
problem-row file size = 0
stderr size = 0
```

Promoted summaries:

```text
work/rh_compute/results/arb_toeplitz_taylor_lam0_N24_o5_r1e-80_nonzero_summary.json
work/rh_compute/results/arb_toeplitz_taylor_lam1e-6_N24_o5_r1e-80_nonzero_summary.json
work/rh_compute/results/arb_toeplitz_taylor_lam1e-4_N24_o5_r1e-80_nonzero_summary.json
work/rh_compute/results/arb_toeplitz_taylor_lam1em2_N24_o5_r1e-80_nonzero_summary.json
work/rh_compute/results/arb_toeplitz_taylor_lam1em1_N24_o5_r1e-80_nonzero_summary.json
problem-row files = matching *_problem_rows.jsonl files, all zero byte
manifest gate = python work/rh_compute/scripts/check_toeplitz_certificate_manifest.py
```

## Non-Rigorous Prefilter

Before launching the exact run, a bounded CPU Torch screen was run:

```text
label = torch_prefilter_taylor_lam0_N24_o5_500k
lambda = 0
N = 24
orders <= 5
max_evaluations = 500,000
device = cpu
dtype = float64
geometric_scale = auto
```

Result:

```text
evaluated_nonstructural_tests = 500,000
positive = 500,000
negative = 0
near_zero = 0
nonfinite = 0
recorded_problem_rows = 0
truncated = true
```

This is Level-1 evidence only. It is useful as a suspicion screen and does not enter the certificate ledger.

## CUDA Scout Layer

The laptop GPU is now usable through an isolated environment:

```text
environment = work/rh_compute/.venv_cuda
torch = 2.12.1+cu126
device = NVIDIA GeForce GTX 1060
GPU memory = 6 GB
compute capability = (6, 1)
```

CUDA is not used for Arb/interval certification. It is only a Level-1 floating-point scout.

The Torch prefilter now has a windowing option:

```text
--skip-evaluations K
```

It uses exact row-block counts to jump skipped structurally nonzero minors before evaluating a bounded window. This makes late-window scouting practical without changing the certificate pipeline.

After a successful CUDA smoke test, four nonzero-lambda `N = 24`, orders `<= 5` screens were run with `max_evaluations = 20,000,000`:

```text
lambda = 1e-6: positive = 20,000,000, problem_rows = 0, min_abs_det = 9.2881922263158266e-05
lambda = 1e-4: positive = 20,000,000, problem_rows = 0, min_abs_det = 9.2882922181446937e-05
lambda = 1e-2: positive = 20,000,000, problem_rows = 0, min_abs_det = 9.2982975522553835e-05
lambda = 1e-1: positive = 20,000,000, problem_rows = 0, min_abs_det = 9.3898253292749975e-05
```

All four row logs are empty. These are not certificate claims; they only reduce the chance that the next exact nonzero-lambda frontier begins with an obvious floating-point obstruction.

Four later windows were then screened for each nonzero lambda:

```text
window offsets = skip 100M, 200M, 300M, 400M
take per window = 1,000,000 evaluated structurally nonzero minors
window files = 16
window evaluated = 16,000,000
window problem rows = 0
```

Combined nonzero-lambda CUDA scout status:

```text
prefix evaluated = 80,000,000
window evaluated = 16,000,000
total evaluated = 96,000,000
total problem rows = 0
global min_abs_det = 2.96489388972214e-05
```

These windowed scouts make the next Arb choices better informed, but they still do not enter the promoted certificate manifest.

## Next-Frontier GPU Scout

The next larger order-5 matrix frontier is:

```text
N = 26
orders <= 5
evaluated_nonstructural_tests per lambda = 939,484,026
structural_zero = 3,617,893,175
tests = 4,557,377,201
```

A first `lambda = 0` CUDA prefix screen found:

```text
label = torch_prefilter_taylor_lam0_N26_o5_cuda_20m
max_evaluations = 20,000,000
positive = 20,000,000
problem_rows = 0
min_abs_det = 7.3550169057392677e-05
elapsed = 95.962 s
```

This supports the next exact frontier choice now that the `N = 24`, orders `<= 5` jobs have certified cleanly. It is not promoted evidence.

## Next Choices

With the full `N = 24`, orders `<= 5` grid promoted:

```text
1. decide whether to push N = 26, orders <= 5 or N > 30 at orders <= 4;
2. continue using GPU scouts only as suspicion screens before exact Arb runs;
3. keep every promoted frontier behind manifest validation and coefficient-enclosure gates.
```

If a future exact frontier fails:

```text
1. preserve the first problem rows;
2. rerun the flagged minors with higher dps or narrower coefficient balls;
3. decide whether the failure is numerical width, coefficient-enclosure weakness, or genuine sign obstruction.
```
