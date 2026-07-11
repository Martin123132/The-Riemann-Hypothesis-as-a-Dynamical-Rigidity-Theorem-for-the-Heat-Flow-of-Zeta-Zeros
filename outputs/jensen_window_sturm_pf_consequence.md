# Jensen-Window Sturm-To-PF Consequence

Date: 2026-07-05

Status: finite consequence artifact. This is not a proof of all-degree or all-shift Jensen hyperbolicity, not all-minor Jensen-window PF-infinity, not
Laguerre-Polya membership, not RH, and not `Lambda <= 0`; it records the finite
PF consequence of the promoted Arb/Sturm root-count diagnostics.

## Purpose

The finite Arb/Sturm diagnostic verifies positive-root counts for:

```text
Q_{d,n,lambda}(y)=P_{d,n,lambda}(-y)
```

where:

```text
P_{d,n,lambda}(x)
  = sum_{j=0}^d binom(d,j) A_{n+j}(lambda) x^j.
```

If `Q_{d,n,lambda}` has exactly `d` positive real roots and positive constant
term, then `P_{d,n,lambda}` has all roots real and nonpositive. By the finite
Polya-frequency characterization, the checked binomial window:

```text
B^{d,n,lambda}_j = binom(d,j) A_{n+j}(lambda)
```

is a finite PF-infinity sequence for that one checked `d,n,lambda`. This uses
the finite Polya-frequency characterization.

## Checked Range

Executable validation:

```text
python work/rh_compute/scripts/check_jensen_window_sturm_pf_consequence.py
```

Input manifests:

```text
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520.jsonl
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520.jsonl
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520_summary.json
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520.jsonl
```

Finite range:

```text
degree d = 3..12
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n = 0..20
```

Current result:

```text
1050/1050 checked Jensen windows have Sturm-certified finite PF-infinity consequence
validated 1050 finite Sturm-to-PF Jensen-window consequences with 0 issues
```

## Boundary

This consequence is finite and window-by-window. It does not prove:

```text
all degrees d;
all shifts n;
all lambda values;
all-minor Jensen-window PF-infinity as an infinite theorem;
the signed-Hankel/Jensen bridge theorem;
Laguerre-Polya membership;
RH or Lambda <= 0.
```

Correct use:

```text
finite PF-infinity certificate for the checked Jensen windows
```

Forbidden use:

```text
Promoting this finite consequence to the Jensen-window PF bridge.
```

## Related Artifacts

```text
outputs/arb_jensen_window_sturm_hyperbolicity_diagnostic.md
outputs/jensen_window_pf_bridge_target.md
outputs/sign_regularity_theorem_fit_matrix.md
python work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py
python work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py
python work/rh_compute/scripts/check_sign_regularity_theorem_fit_matrix.py
```
