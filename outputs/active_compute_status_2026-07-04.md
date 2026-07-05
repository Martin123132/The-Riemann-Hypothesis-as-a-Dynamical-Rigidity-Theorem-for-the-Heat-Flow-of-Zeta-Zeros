# Active Compute Status - 2026-07-04

Date: 2026-07-04

Status: live operational note. This is not a proof artifact and does not promote any new certificate by itself.

## Exact Arb Jobs

The full `N = 24`, orders `<= 5`, Taylor-normalized lambda grid completed cleanly and is now promoted in the Toeplitz certificate manifest. Current status:

```text
lambda = 0      label = arb_toeplitz_taylor_lam0_N24_o5_r1e-80_nonzero      completed cleanly and promoted
lambda = 1e-6   label = arb_toeplitz_taylor_lam1e-6_N24_o5_r1e-80_nonzero   completed cleanly and promoted
lambda = 1e-4   label = arb_toeplitz_taylor_lam1e-4_N24_o5_r1e-80_nonzero   completed cleanly and promoted
lambda = 1e-2   label = arb_toeplitz_taylor_lam1em2_N24_o5_r1e-80_nonzero   completed cleanly and promoted
lambda = 1e-1   label = arb_toeplitz_taylor_lam1em1_N24_o5_r1e-80_nonzero   completed cleanly and promoted
```

Certified counts per lambda:

```text
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

Promoted grid result:

```text
lambda grid = {0, 1e-6, 1e-4, 1e-2, 1e-1}
summary files present = true
problem-row files are zero byte = true
bad sign/width counts = 0 for all five summaries
manifest validates 95 promoted positive Toeplitz summaries and 1 beta negative control
```

The old nonzero-lambda process ids are no longer running. Their completed summaries are:

```text
work/rh_compute/results/arb_toeplitz_taylor_lam1e-6_N24_o5_r1e-80_nonzero_summary.json
work/rh_compute/results/arb_toeplitz_taylor_lam1e-4_N24_o5_r1e-80_nonzero_summary.json
work/rh_compute/results/arb_toeplitz_taylor_lam1em2_N24_o5_r1e-80_nonzero_summary.json
work/rh_compute/results/arb_toeplitz_taylor_lam1em1_N24_o5_r1e-80_nonzero_summary.json
```

## GPU Scout Layer

CUDA Torch is available in:

```text
work/rh_compute/.venv_cuda
torch = 2.12.1+cu126
device = NVIDIA GeForce GTX 1060
```

Nonzero-lambda `N = 24`, orders `<= 5` scout:

```text
prefix evaluated = 80,000,000
window evaluated = 16,000,000
total evaluated = 96,000,000
total problem rows = 0
global min_abs_det = 2.96489388972214e-05
```

Next-frontier `N = 26`, orders `<= 5`, `lambda = 0` prefix scout:

```text
evaluated = 20,000,000
positive = 20,000,000
problem_rows = 0
min_abs_det = 7.3550169057392677e-05
```

All GPU results are Level-1 warning screens only. They do not enter the promoted certificate manifest.

## Additional Interval Diagnostics

Signed-Hankel certificate manifest:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
m = 0..19
s = 0..24
validated 2,500 signed-Hankel finite certificates
```

Edrei-log necessary-condition diagnostic:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 1..64
validated 320 finite Edrei-log sign diagnostics
```

Edrei power-Hankel necessary-condition diagnostic:

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

The larger scout to `s = 48` found no negative determinants, but high-shift
cases became interval-width inconclusive and are not promoted. The upper
staircase and lower high-shift wedge were promoted only after recomputing
coefficient enclosures through `k = 49` with `dps = 180`,
`abs_tol = 1e-120`, and recomputing Edrei-log rows through `n = 49` at
`dps = 1000`.

The `n = 57` frontier now promotes the full five-lambda staircase
`2m+s <= 57`. The `n = 51` layer first needed a `lambda = 1e-6`,
`k = 0..51`, `dps = 220`, `abs_tol = 1e-140` repair. The `n = 53` and
later layers needed the same tighter coefficient treatment through `k = 57`
for all five lambdas, followed by Edrei-log rows at `dps = 2400`.

Boundary history and repair are executable:

```text
python work/rh_compute/scripts/check_edrei_power_hankel_boundary_manifest.py
validated 2 retired inconclusive blocker rows and 3 repaired positive boundary rows
```

Coefficient enclosures now include:

```text
k = 0..64 for the same five-lambda grid
latest extension k = 49..64:
  rows = 80
  max c_ball radius about 1.33e-259
tight Edrei frontier extension k = 48..49:
  rows = 10
  max c_ball radius about 2.04e-271
tight Edrei frontier extension k = 50..51:
  rows = 10
  max c_ball radius about 7.46e-279
tight Edrei frontier extension k = 52..53:
  rows = 10
  max c_ball radius about 8.92e-287
tight Edrei boundary extension k = 0..53:
  rows = 54 per lambda
  dps = 220
  abs_tol = 1e-140
tight Edrei boundary extension k = 54..55:
  rows = 2 per lambda
  dps = 220
  abs_tol = 1e-140
tight Edrei boundary extension k = 56..57:
  rows = 2 per lambda
  dps = 220
  abs_tol = 1e-140
```

These are finite interval diagnostics. They do not prove the missing all-order bridge theorem.
