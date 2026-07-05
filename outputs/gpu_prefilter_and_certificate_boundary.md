# GPU Prefilter and Certificate Boundary

Date: 2026-07-04

Status: implementation note. This is not a proof artifact and not a replacement for Arb/interval certificates.

## Local Hardware Finding

The laptop exposes:

```text
NVIDIA GeForce GTX 1060
6 GB GPU memory
driver-reported CUDA runtime: 13.0
```

The default Python environment has:

```text
torch 2.12.0+cpu
torch.cuda.is_available() = False
```

An isolated CUDA-enabled Torch environment now also exists:

```text
work/rh_compute/.venv_cuda
torch 2.12.1+cu126
torch.cuda.is_available() = True
torch.version.cuda = 12.6
device = NVIDIA GeForce GTX 1060
compute capability = (6, 1)
```

This environment was installed from the official PyTorch CUDA wheel index:

```text
https://download.pytorch.org/whl/cu126
```

The default Arb/certificate environment was left untouched.

## Added Script

```text
work/rh_compute/scripts/torch_toeplitz_prefilter.py
```

Purpose:

```text
fast non-rigorous screening of finite Toeplitz/PF minors
```

It reads the same coefficient cache as the Arb probe and enumerates the same structurally nonzero upper-triangular Toeplitz minors:

```text
c_i >= r_i for every i
```

It then evaluates floating-point Torch determinants in batches.

The script also supports windowed scouting:

```text
--skip-evaluations K
```

This skips the first `K` structurally nonzero minors before evaluating. The skip path uses exact row-block counting, so late windows can be reached without enumerating every skipped minor one by one.

## Sign-Preserving Scaling

The script optionally rescales:

```text
c_k -> c_k p^k, with p > 0.
```

For a Toeplitz minor with row set `R` and column set `C`, this multiplies the determinant by:

```text
p^(sum(C) - sum(R)).
```

This is positive, so signs are preserved. The point is numerical conditioning only: raw Taylor-normalized coefficients can span enough orders of magnitude that direct float determinants become useless.

## Validation Checks

Corrected Taylor sequence, CPU Torch:

```text
lambda = 0
sequence = taylor
N = 10
orders <= 4
tests_accounted = 60,625
evaluated_nonstructural_tests = 19,690
positive = 19,690
structural_zero = 40,935
negative = 0
near_zero = 0
nonfinite = 0
```

Corrected Taylor sequence, CUDA Torch:

```text
lambda = 0
sequence = taylor
N = 10
orders <= 4
tests_accounted = 60,625
evaluated_nonstructural_tests = 19,690
positive = 19,690
structural_zero = 40,935
negative = 0
near_zero = 0
nonfinite = 0
```

Known-bad beta normalization:

```text
lambda = 0
sequence = beta
N = 10
orders <= 3
tests_accounted = 16,525
evaluated_nonstructural_tests = 5,830
negative = 703
```

This matches the Arb negative-control count, so the prefilter can detect the known failed normalization while agreeing with the corrected small Taylor case.

## Performance Note

For the current workload of many very small Toeplitz determinants, the GTX 1060 is not automatically faster than CPU Torch at small caps, but it becomes useful once enough batched determinant work is queued.

Corrected post-cap benchmark:

```text
lambda = 0
sequence = taylor
N = 20
orders <= 4
max_evaluations = 100,000
dtype = float64

CPU Torch:
elapsed = 4.642 s
evaluated_nonstructural_tests = 100,000
positive = 100,000

CUDA Torch:
elapsed = 6.383 s
evaluated_nonstructural_tests = 100,000
positive = 100,000
```

Interpretation:

```text
GPU acceleration is available. CPU Torch remains a good default for small quick probes, but
the CUDA path is now useful for larger Level-1 frontier scouts, especially orders <= 5 at
N = 24 with million-scale caps.
```

## Proof Boundary

The GPU/Torch prefilter is Level 1 evidence only:

```text
non-rigorous floating-point screen
```

It may be used to:

```text
1. find suspicious negative or near-zero minors quickly;
2. choose which larger matrix/order ranges are worth Arb certification;
3. estimate numerical fragility before spending CPU time;
4. stress candidate theorem statements before manuscript wording hardens.
```

It may not be used to claim:

```text
certified finite computation;
PF-infinity;
Laguerre-Polya membership;
RH;
Lambda <= 0.
```

Every GPU-flagged negative, near-zero, nonfinite, or frontier-clean range must be followed by the Arb/interval pipeline before entering the certificate status document.

## Recommended Next Use

1. Use CPU Torch as the default prefilter backend for order `<= 4`.
2. Use the CUDA venv only for explicit experiments, e.g. larger-order determinant batches or custom kernels.
3. Feed suspicious minors or promising frontiers back into `arb_toeplitz_pf_probe.py`.
4. Keep all manuscript-facing claims tied to Arb/interval outputs, not Torch summaries.

## Frontier Screen Added

A bounded CPU Torch screen was run before the `N = 24`, orders `<= 5`, `lambda = 0` exact Arb attempt:

```text
lambda = 0
sequence = taylor
N = 24
orders <= 5
max_evaluations = 500,000
dtype = float64
geometric_scale = auto
positive = 500,000
negative = 0
near_zero = 0
nonfinite = 0
```

This remains non-rigorous and is recorded only as a prefilter. The exact Arb promotion conditions and completed certificate ledger are tracked in:

```text
outputs/toeplitz_frontier_run_ledger.md
```

## CUDA Frontier Scouts Added

After installing the isolated CUDA Torch environment, a real GPU smoke test succeeded:

```text
torch = 2.12.1+cu126
torch.cuda.is_available() = True
device = NVIDIA GeForce GTX 1060
compute capability = (6, 1)
float64 2048x2048 matmul elapsed = 0.4328 s
```

A first `lambda = 1e-6`, `N = 24`, orders `<= 5` CUDA prefilter with `max_evaluations = 1,000,000` found:

```text
positive = 1,000,000
negative = 0
near_zero = 0
nonfinite = 0
recorded_problem_rows = 0
elapsed = 8.343 s
```

A broader nonzero-lambda CUDA scout first used `max_evaluations = 5,000,000` for each lambda, then was extended to `max_evaluations = 20,000,000` for each nonzero lambda:

```text
lambda = 1e-6: positive = 20,000,000, problem_rows = 0, min_abs_det = 9.2881922263158266e-05, elapsed = 47.516 s
lambda = 1e-4: positive = 20,000,000, problem_rows = 0, min_abs_det = 9.2882922181446937e-05, elapsed = 47.576 s
lambda = 1e-2: positive = 20,000,000, problem_rows = 0, min_abs_det = 9.2982975522553835e-05, elapsed = 46.376 s
lambda = 1e-1: positive = 20,000,000, problem_rows = 0, min_abs_det = 9.3898253292749975e-05, elapsed = 45.941 s
```

All four row logs are empty. These results are useful only for frontier triage; none of them is promoted to the certificate manifest.

Windowed scouting was then added and tested. The same `lambda = 1e-6`, `skip = 100,000,000`, `take = 1,000,000` window fell from `84.301 s` with naive skipping to `8.862 s` after row-block skip acceleration, with identical screen output:

```text
positive = 1,000,000
problem_rows = 0
min_abs_det = 2.9648967773718437e-05
```

The nonzero-lambda `N = 24`, orders `<= 5` CUDA scout now consists of:

```text
prefix screens: 4 files x 20,000,000 = 80,000,000 evaluated minors
window screens: 4 lambdas x 4 offsets x 1,000,000 = 16,000,000 evaluated minors
window offsets: skip 100M, 200M, 300M, 400M
total evaluated: 96,000,000
total problem rows: 0
global min_abs_det across these files: 2.96489388972214e-05
```

This is a stronger warning screen for the next exact Arb frontier choices, but it remains Level-1 evidence only.

## Next-Frontier Scout

The next matrix-size frontier is substantially larger:

```text
N = 26
orders <= 5
evaluated_nonstructural_tests per lambda = 939,484,026
structural_zero = 3,617,893,175
total tests = 4,557,377,201
```

A first CUDA prefix scout was run for `lambda = 0`:

```text
sequence = taylor
N = 26
orders <= 5
max_evaluations = 20,000,000
device = cuda
dtype = float64
positive = 20,000,000
problem_rows = 0
min_abs_det = 7.3550169057392677e-05
elapsed = 95.962 s
```

This is only a warning screen for a possible later exact `N = 26` Arb frontier.
