# Jensen-Window PF Negative-Lambda -100 Full Cone-Entry Certificate

Date: 2026-07-10

Status: interval full cone-entry theorem at lambda=-100 and open flow-legitimacy handoff. This is not a proof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py
```

## Repaired Prefix

The broad lambda=-100 coefficient source is overridden by the dps220
repairs on `k=220..250` and `k=245..320`. Direct Arb arithmetic proves

```text
A_k(-100)>0, 0<=k<=320,
(2*k-1)/(2*k+1)<x_k<1, 1<=k<=319,
x_(k+1)>x_k, 1<=k<=318.
```

The weakest prefix adjacent margin is `3.70342889106601568930111785401330143358127882878709479132476E-6` at `k=318`.

## Global Composition

The exact Cauchy-Schwarz and Mellin log-concavity theorems prove the
pointwise lower and upper walls for every real heat parameter and every
`k>=1`. The paired-remainder wall closure proves the adjacent tail for
every `k>=319`. Hence

```text
(2*k-1)/(2*k+1)<=x_k(-100)<=1,
x_(k+1)(-100)>=x_k(-100),
for every integer k>=1.
```

This is full infinite ratio-cone entry at `lambda=-100`.

## Remaining Handoff

The exact ratio-cone boundary algebra is already available. What remains
is a rigorous infinite-dimensional or collared finite-flow legitimacy
argument before propagating the cone forward in lambda. This certificate
does not promote cone entry directly to the all-order Jensen/Newman result.

```text
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
