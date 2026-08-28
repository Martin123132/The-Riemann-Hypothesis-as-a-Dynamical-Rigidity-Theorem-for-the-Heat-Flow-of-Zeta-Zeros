# Nonoscillatory absolute bound for the translated A-face core

Date: 2026-08-23

Status: rigorous local interval certificate at `t=10^10`, independently
replayed on a different partition; not a proof of `R_after_A` or RH

On `0<x<=x_16`, Section 11.454 gives `D_42=D_6+R_6`.  Apply the positive
midpoint Peano kernel directly to the rational six-current amplitude:

```text
|D_6(x)|
 <=sum_(m=39895)^39936 integral_cell K_m(y)|G_(6,yy)(y,x)|dy.
                                                               (CP1)
```

Termwise differentiation of the exact rational formula factors

```text
G_(6,yy)(y,x)=x M_6(y,x).                           (CP2)
```

No subtraction of the large point sum and strip integral is therefore needed.
With `q=log((1-x)/x)` and `|dx|=x(1-x)dq`, (CP1)--(CP2) give

```text
integral_0^x16 |D_6(x)|x^(-5/4)(1-x)^(-1/4)dx
 <=integral_q16^infinity B_6(x(q))[x(q)(1-x(q))]^(3/4)dq,

B_6(x)=sum_cells integral_cell K_m(y)|M_6(y,x)|dy.  (CP3)
```

The replacement remainder has the still faster regular weight
`x^(23/4)(1-x)^(3/4)`.  The builder encloses the finite interval through
`q=40`, then uses `x(q)<=exp(-q)` and total Peano mass `7/4` for analytic
tails.  At 100 decimal digits its `768`
outer panels and `4` subdivisions per half-cell certify

```text
six-current core contribution <= [1.017349068511713951477088554795874402671360700686272943287947121863281e-5 +/- 3.07e-76],
six-current replacement       <= [1.658535868418944840242605119088391835141506846142221682189364394446745e-15 +/- 3.53e-86],

complete pre-endpoint core    <= [1.017349068677567538318983038820134914580199884200423627902169290082217e-5 +/- 4.09e-75]
                               < 0.000012.              (CP4)
```

The independent checker uses 120 decimal digits, three times the base outer
partition, three subdivisions per half-cell, exact endpoint kernel maxima,
and its own termwise derivative majorant.  It independently obtains a bound
below `0.000012`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_core_absolute_peano_interval_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_core_absolute_peano_interval_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_core_absolute_peano_interval_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_core_absolute_peano_interval_gate.py
```

Pi provenance: every `pi` in (CP1)--(CP4) is inherited from the canonical
Fresnel hierarchy or equation-(9) normalization.  The logistic coordinate
introduces no new occurrence of `pi`.

Proof boundary: (CP4) proves only the translated A-face pre-endpoint core on
`0<x<=x_16` at the saved height.  It does not by itself include the endpoint
layer, assemble the full translated A-face channel, or prove `R_after_A`,
`R_Dir`, `Q_K-T`, an all-height theorem, `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
