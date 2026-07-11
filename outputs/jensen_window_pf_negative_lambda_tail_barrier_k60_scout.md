# Jensen-Window PF Negative-Lambda Tail-Barrier Scout

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of zeta
cone entry, Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_tail_barrier_scout`.

Proof boundary: this artifact rewrites the missing negative-lambda tail
as defect inequalities and certifies finite-prefix barrier facts only.
It does not prove an all-`k` tail theorem, does not prove a collared
finite flow theorem, does not prove `jwpf_06`, and does not establish
`Lambda <= 0`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k60_scout.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k60.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_tail_barrier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py --target work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k60_scout.json --note outputs/jensen_window_pf_negative_lambda_tail_barrier_k60_scout.md
```

Current result:

```text
validated Jensen-window PF negative-lambda tail-barrier scout: 177 cone-buffer rows, 174 defect-monotone rows, 139 one-third-buffer rows, 174 scaled-defect increase rows, 1 rejected candidate, 0 issues
```

## Defect Form

Write:

```text
d_k = 1 - x_k
x_k = (A_{k+1}/A_k)/(A_k/A_{k-1})
```

Then the ratio cone becomes:

```text
0 <= d_k <= 2/(2*k+1)
d_{k+1} <= d_k
```

This is the exact tail-barrier form that an analytic proof must supply
outside the finite prefix.

## Finite Barrier Diagnostics

```text
lambdas: -100.0, -50.0, -25.0
coefficient range: A_0..A_60
checked contractions: x_1..x_59
cone-buffer rows: 177 / 177
one-third width buffer rows: 139 / 177
defect monotone rows: 174 / 174
scaled-defect increase rows: 174 / 174
```

Minimum and maximum sampled margins:

```text
min upper defect d_k: 4.268538754539301236E-3 at lambda=-100.0, k=59
min lower-wall slack: 9.875541374763750245E-3 at lambda=-25.0, k=59
min one-third buffer margin: -7.907195486822352710E-2 at lambda=-25.0, k=59
min defect monotone gap: 3.050976509934343254E-5 at lambda=-100.0, k=58
max scaled defect: 4.124052882015568604E-1 at lambda=-25.0, k=59
max defect ratio d_{k+1}/d_k: 9.929031354356735843E-1 at lambda=-100.0, k=58
```

## Rejected Shortcut

The scaled defect

```text
s_k = ((2*k+1)/2) * d_k
```

is increasing on every checked adjacent pair. Therefore a proof route
that tries to obtain the lower wall by assuming `s_k` is nonincreasing
is already incompatible with the certified finite prefix.

## Next Upgrade

The current finite-collar contract has active depth `K=57` with collars `x_58` and `x_59`. The next purely finite upgrade to `K=58` needs
`x_60`, hence `A_61`, unless an analytic tail theorem supplies the
missing collar.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```

Summary:

The certified prefix satisfies the defect-form cone barriers 0 <= d_k <= 2/(2k+1), d_{k+1} <= d_k, and the stronger finite one-third width buffer on 139/177 checked rows, with a frontier failure before the end of the prefix. However, the scaled defect is increasing on the checked prefix, so any analytic proof route that assumes scaled-defect monotone decrease is already rejected. The next finite collar needs x_60, hence A_61, unless an analytic tail theorem supplies it.
