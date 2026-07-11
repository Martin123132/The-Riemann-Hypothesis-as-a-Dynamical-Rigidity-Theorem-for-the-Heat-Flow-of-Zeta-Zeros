# Arb Jensen-Window Sturm Hyperbolicity Diagnostic

Date: 2026-07-05

Status: finite Arb/Sturm diagnostic. This is not a proof of Jensen-window
PF-infinity, all-degree Jensen hyperbolicity, Laguerre-Polya membership, RH,
or `Lambda <= 0`; it certifies only bounded families of degree-3 through
degree-12 Jensen-window root counts.

## Purpose

This diagnostic is the direct-root companion to the Arb Jensen-window PF
obligation diagnostic. For each lambda, shift, and degree, form:

```text
P_{d,n,lambda}(x)
  = sum_{j=0}^d binom(d,j) A_{n+j}(lambda) x^j

Q_{d,n,lambda}(y) = P_{d,n,lambda}(-y).
```

The finite check verifies that `Q_{d,n,lambda}` has exactly `d` positive real
zeros on the selected grid. Equivalently, the corresponding checked Jensen
polynomials have real nonpositive zeros on this bounded range.

Executable artifacts:

```text
python work/rh_compute/scripts/arb_jensen_window_sturm_hyperbolicity_probe.py
python work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520.jsonl
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520.jsonl
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520_summary.json
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520.jsonl
```

Current degree-3/4 result:

```text
validated 210 Arb Jensen-window Sturm hyperbolicity finite diagnostics with 0 issues
```

Current degree-5 result:

```text
python work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py --summary work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json --rows work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520.jsonl --expected-degrees 5
validated 105 Arb Jensen-window Sturm hyperbolicity finite diagnostics with 0 issues
```

Current degree-6 through degree-12 extension:

```text
python work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py --summary work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520_summary.json --rows work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520.jsonl --expected-degrees 6..12
validated 735 Arb Jensen-window Sturm hyperbolicity finite diagnostics with 0 issues
```

## What Was Checked

The probe uses rigorous `A_ball` coefficient enclosures from:

```text
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl
```

Finite range:

```text
degree d = 3,4 in the base manifest
degree d = 5 in the extension manifest
degree d = 6..12 in the extended terminal-reconnaissance manifest
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n = 0..20
dps = 520
needed max coefficient index = 24 for d=3,4
needed max coefficient index = 25 for d=5
needed max coefficient index = 32 for d=6..12
```

For each row, the probe builds an interval-enclosed Sturm sequence for
`Q_{d,n,lambda}` and counts positive roots by sign variations at `0` and
`+infinity`.

The manifest contains:

```text
210/210 lambda-shift-degree rows positive-root-count certified
105/105 degree-3 rows
105/105 degree-4 rows
105/105 degree-5 extension rows positive-root-count certified
735/735 degree-6 through degree-12 extension rows positive-root-count certified
1050/1050 total degree-3 through degree-12 rows certified
0 failed or inconclusive
```

## Boundary

This is useful pressure on the Jensen-window PF bridge target, but it is still
a finite diagnostic. It does not check:

```text
all degrees d;
all shifts n;
all Jensen-window Toeplitz minors;
Jensen-window PF-infinity;
the limiting Laguerre-Polya/Newman conclusion.
```

Correct use:

```text
finite Arb/Sturm-backed theorem-search evidence for low-degree Jensen-window
hyperbolicity
```

Forbidden use:

```text
proof of all-degree Jensen hyperbolicity, Jensen-window PF-infinity, or
Lambda <= 0
```

## Related Gates

```text
outputs/arb_jensen_window_pf_obligation_diagnostic.md
outputs/jensen_window_sturm_pf_consequence.md
outputs/jensen_window_pf_bridge_target.md
outputs/jensen_window_pf_obligation_algebra.md
python work/rh_compute/scripts/check_arb_jensen_window_pf_obligation_manifest.py
python work/rh_compute/scripts/check_jensen_window_sturm_pf_consequence.py
python work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py
python work/rh_compute/scripts/check_jensen_window_pf_obligation_algebra.py
```
