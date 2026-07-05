# Progress Report

Date: 2026-07-03

Status: sprint progress report. This is not a proof artifact; it summarizes implementation and evidence-quality work completed in the stated interval.

Goal:

Advance the RH dynamical-rigidity corpus toward a referee-legible proof programme for the missing Newman direction `Lambda <= 0`, with priority on the signed Hankel / Jensen / sign-regularity route.

## Work Completed This Sprint

### 1. Runner Hardened

Updated:

```text
work/rh_compute/scripts/jensen_hankel_runner.py
```

New capabilities:

- `--overwrite` to avoid duplicate appended rows;
- coefficient cache output via `--write-coeff-cache`;
- mixed Jensen/Hankel CSV support;
- near-zero determinant flagging via `--near-zero-log10`;
- summary JSON includes near-zero counts and coefficient-cache path.

Smoke tests passed after patching:

```text
smoke_hankel_v2: 6/6 signed-Hankel cases passed
smoke_both_v2:   4/4 signed-Hankel cases and 4/4 Jensen cases passed
```

## 2. Signed-Hankel Checkpoint 15C Reproduced

Run:

```text
repro_hankel_15c
```

Parameters:

```text
lambda = 0, 1e-6, 1e-4, 1e-2, 1e-1
m = 0..12
shift = 0..8
dps = 90
n_sum = 100
cutoff = 8
```

Result:

```text
585/585 signed determinants positive
failures: 0
```

Artifacts:

```text
work/rh_compute/results/repro_hankel_15c.jsonl
work/rh_compute/results/repro_hankel_15c.csv
work/rh_compute/results/repro_hankel_15c_summary.json
work/rh_compute/results/repro_hankel_15c_coefficients.jsonl
```

Important caveat:

```text
near-zero determinants: 240/585
smallest log10(abs(det)): about -735.09
```

2026-07-04 update:

This was originally strong reproducibility evidence but not an interval certificate. It has now been superseded by the enclosure-backed Arb check recorded in:

```text
outputs/signed_hankel_certificate_status.md
```

which now validates an extended `2,500/2,500` finite signed-Hankel determinant grid using rigorous `A_ball` coefficient enclosures.

## 3. High-Precision Stability Spot Check

Run:

```text
stability_hankel_lambda0_m9_m12
```

Parameters:

```text
lambda = 0
m = 9..12
shift = 0..8
dps = 160
n_sum = 140
cutoff = 9
```

Result:

```text
36/36 signed determinants positive
failures: 0
sign flips versus 90-digit run: 0
```

Artifacts:

```text
work/rh_compute/results/stability_hankel_lambda0_m9_m12.jsonl
work/rh_compute/results/stability_hankel_lambda0_m9_m12.csv
work/rh_compute/results/stability_hankel_lambda0_m9_m12_summary.json
work/rh_compute/results/stability_hankel_lambda0_m9_m12_coefficients.jsonl
```

Interpretation:

The most fragile lambda=0 high-order block is stable under a stricter computation, but still needs interval/ball arithmetic before it can be called certified.

## 4. Jensen Fragile Case Reproduced

Run:

```text
repro_jensen_fragile_d16_n54
```

Parameters:

```text
lambda = 0
d = 16
n = 54
dps = 100
rational_digits = 70
n_sum = 120
cutoff = 9
```

Result:

```text
positive root count of Q_{16,54}: 16
expected: 16
pass: true
```

Artifacts:

```text
work/rh_compute/results/repro_jensen_fragile_d16_n54.jsonl
work/rh_compute/results/repro_jensen_fragile_d16_n54.csv
work/rh_compute/results/repro_jensen_fragile_d16_n54_summary.json
work/rh_compute/results/repro_jensen_fragile_d16_n54_coefficients.jsonl
```

Important caveat:

The Jensen count is exact only after rationalizing high-precision numerical coefficients. It is not yet a proof about the exact transcendental moments.

## 5. Theorem Map Created

Created:

```text
outputs/sign_regularity_theorem_map.md
```

Main conclusion:

The right route is not naive PF-infinity Toeplitz positivity of `A_k`, because the bundle reports robust negative Toeplitz minors. The current best route is:

```text
signed/indefinite Hankel sign-regularity
  -> find correct transformation or kernel-level condition
  -> connect to all Jensen hyperbolicity, Laguerre-Polya membership,
     or no positive Newman boundary.
```

Priority sources identified:

- Groechenig on Schoenberg total positivity and zeta;
- Grussler-Damm on Hankel sign-consistency verification;
- ASW/Edrei/PF-infinity theory as the classical real-rootedness benchmark;
- Farmer as a guardrail against overreading Jensen asymptotics.

Update added 2026-07-03:

This conclusion needs a normalization correction. The negative Toeplitz minors are real for `A_k`, but ASW/Edrei applies to ordinary Taylor coefficients. In the present notation, the ordinary coefficients are `c_k = A_k/k!`. A later exact-rational finite audit found no negative Toeplitz minors for `c_k` in the tested ranges, so the `c_k` PF-infinity route is now active again as a theorem-search target. See `outputs/schoenberg_zeta_condition.md` and `outputs/progress_report_2026-07-03_toeplitz_normalization.md`.

## Current Status

What is stronger now:

- The signed-Hankel finite evidence is no longer just notebook output; it has clean JSONL/CSV reproduction.
- The fragile high-order block has a high-precision stability check.
- The Jensen known-fragile case passes in the clean runner.
- The theorem route is better scoped: signed/indefinite sign-regularity for `A_k`, plus a corrected ordinary-coefficient PF-infinity route for `A_k/k!`.

What remains incomplete:

- 2026-07-04 update: finite signed-Hankel interval certificates now exist for an extended `2,500/2,500` grid; finite Toeplitz/PF certificates also exist for the normalized `c_k = A_k/k!` ranges recorded in `outputs/finite_toeplitz_certificate_status.md`.
- No full Jensen reproduction grid yet.
- No theorem yet connecting signed Hankel signs to all-degree/all-shift Jensen hyperbolicity.
- No theorem yet implying `Lambda <= 0`.

## Next Gates

### Gate 1: Certification Upgrade

Install or otherwise access Arb/FLINT-style ball arithmetic, preferably through:

```text
python-flint
```

Then redo selected near-zero determinants with certified error bounds.

### Gate 2: Theorem Extraction

Read and extract Groechenig's exact Schoenberg/zeta RH-equivalent condition. Compare it against:

```text
Phi_lambda translation minors
A_k signed Hankel minors
Jensen hyperbolicity
```

### Gate 3: Hankel Minor Scope

Use Grussler-Damm to determine whether consecutive shifted Hankel minors can imply broader `k`-sign consistency for the Hankel matrix generated by `A_k`.

### Gate 4: Expand Computation Carefully

Do not run huge Jensen grids until the runner has:

- coefficient cache reuse;
- resumable per-lambda/per-degree mode;
- interval or independent precision comparison for suspicious cases.

Recommended next compute step:

```text
lambda = 0
d = 16..20
n = 50..60
```

This targets the previously fragile region without launching the full grid.
