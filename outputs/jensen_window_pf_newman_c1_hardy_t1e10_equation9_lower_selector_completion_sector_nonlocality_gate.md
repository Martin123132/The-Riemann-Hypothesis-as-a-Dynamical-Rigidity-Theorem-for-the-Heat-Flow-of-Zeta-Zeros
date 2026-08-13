# Selector completion-sector nonlocality

Date: 2026-08-13

Status: rigorous completion-sector barrier for the selector fold packet;
not a fixed-state source-minus-target estimate

The exact one-cell Poisson completion and paired fold selector identity use
the same source normalization.  With

```text
P_C=Q_half^(C)-Q_half^(C+2),
D_fold=F_C-F_(C+2),
C_out=P_C-D_fold,
```

finite algebra gives

```text
P_C=D_fold+C_out.                                      (CN1)
```

The exact completed common profile satisfies

```text
P_C=kappa_C g_ex(lambda,0),
kappa_C=exp(i*beta^3)2sqrt(2).
```

The beta-minus-four endpoint profile is enclosed on the whole ordinary top
corridor with `lambda` kept as one interval.  Appending the certified
Kummer-ODE exact-profile error and applying odd-square reflection gives

```text
[233.6584800233249421741850821221327660915685896046846649650317494718624 +/- 1.44e-68] < Q_C
 < [233.6877721256657635536148032295636927749201559661766743526828749290710 +/- 3.69e-68].        (CN2)
```

The independent two-packet theorem supplies

```text
[29.53042617185773784003499392716790466468402423517815945484063879176039 +/- 4.07e-69] < Delta Q_fold
 < [32.06339161020397093705532841772994983869768915680920807149458612209941 +/- 3.80e-69].          (CN3)
```

Projecting (CN1) and subtracting interval endpoints therefore proves

```text
[201.5950884131209712371297537044028162528709004478754568935371633497630 +/- 3.11e-68] < Q_out
 < [204.1573459538080257135798093023957881102361317309985148978422361373107 +/- 4.74e-68],
Q_out=2Re[exp(-i*pi/8)C_out].                          (CN4)
```

In particular `Q_out>200`, while `Delta Q_fold<35`.  The fold packet is less
than `[0.1372233167270592747562458949074058667412026684522342839152754460246390 +/- 2.02e-72]` of the completed source
projection, and the outside completion exceeds it by a factor greater than
`[6.287391267396821909342775699524878526433515923484913504939797314889040 +/- 2.92e-70]`.

This is a nonlocality result, not a small-error result.  A fixed-B selector
jump is not closed inside the proposed fold ownership block: its endpoint
half-current, zero and other Fourier sectors carry a mandatory projected
contribution above `200`.  Consequently neither the 200-mode fold packet nor
the 39273-mode positive target roster may replace the completed current one
cell at a time.  The fixed-state endpoint residual must be estimated in the
global cancellation-preserving form of Section 11.350, or in another form
that proves the same completion cancellation explicitly.

Pi provenance: `beta^3=pi*C^2/8`, the Fourier phases `2*pi*n`, and the
half-Kummer phase `pi/8` all descend from the exact equation-(9) transform.
No fitted or geometric constant is introduced.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_selector_completion_sector_nonlocality_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_selector_completion_sector_nonlocality_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_selector_completion_sector_nonlocality_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_selector_completion_sector_nonlocality_gate.py
```

No bound for the global fixed-state endpoint residual, complete `Q_K-T` or
`T_upper`, all-corridor or height-uniform theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
