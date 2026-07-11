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
work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k80_scout.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k80.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_tail_barrier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py --target work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k80_scout.json --note outputs/jensen_window_pf_negative_lambda_tail_barrier_k80_scout.md
```

Current result:

```text
validated Jensen-window PF negative-lambda tail-barrier scout: 237 cone-buffer rows, 234 defect-monotone rows, 159 one-third-buffer rows, 234 scaled-defect increase rows, 1 rejected candidate, 0 issues
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
coefficient range: A_0..A_80
checked contractions: x_1..x_79
cone-buffer rows: 237 / 237
one-third width buffer rows: 159 / 237
defect monotone rows: 234 / 234
scaled-defect increase rows: 234 / 234
```

Minimum and maximum sampled margins:

```text
min upper defect d_k: 3.745474659438155329E-3 at lambda=-100.0, k=79
min lower-wall slack: 6.974322101168841388E-3 at lambda=-25.0, k=79
min one-third buffer margin: -1.122080596237437763E-1 at lambda=-25.0, k=79
min defect monotone gap: 2.281394370544295027E-5 at lambda=-100.0, k=78
max scaled defect: 4.455413929570771097E-1 at lambda=-25.0, k=79
max defect ratio d_{k+1}/d_k: 9.939458077371220978E-1 at lambda=-100.0, k=78
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

The current finite-collar contract has active depth `K=77` with collars `x_78` and `x_79`. The next purely finite upgrade to `K=78` needs
`x_80`, hence `A_81`, unless an analytic tail theorem supplies the
missing collar.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```

Summary:

The certified prefix satisfies the defect-form cone barriers 0 <= d_k <= 2/(2k+1), d_{k+1} <= d_k, and the stronger finite one-third width buffer on 159/237 checked rows, with a frontier failure before the end of the prefix. However, the scaled defect is increasing on the checked prefix, so any analytic proof route that assumes scaled-defect monotone decrease is already rejected. The next finite collar needs x_80, hence A_81, unless an analytic tail theorem supplies it.
