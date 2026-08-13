# Corrected selected source-projection orientation

Date: 2026-08-13

Status: exact normalization and conditional signed projection theorem; the
embedding into the isolated fold-owned paired-residual block remains open

Let `J_m(y)` denote the raw transformed fold integral with its common source
carrier suppressed.  The common-profile Fourier normalization is

```text
G_m=(1/(2pi)) Integral_0^Y J_m(y)dy.                 (SO1)
```

The exact joint Kummer normal form carries amplitude `sqrt(2)/pi`, so the
original half-domain source mode is

```text
K_m^(half)=(sqrt(2)/pi)Integral_0^Y J_m(y)dy
          =2sqrt(2)G_m.                              (SO2)
```

Restoring the suppressed carrier and applying the exact half-domain
reflection gives

```text
2 Re[e^(-i*pi/8)e^(i*beta^3) 2sqrt(2)W_corr].       (SO3)
```

Here `beta^3=pi*C^2/8`, `C=159577=1 mod 8`, and
`(C^2-1)/8=3183102366` is even.
Therefore the phase in (SO3) is exactly one and

```text
R_corr=4sqrt(2) Re W_corr.                           (SO4)
```

Using the signed theorem of Section 11.402 proves uniformly on the ordinary
top corridor

```text
R_corr<[-0.00012633020694136396822650036607776931238404392877833648447949440001439331656143490 +/- 2.31e-84]<-1.263e-4,
|R_corr|<[0.00023129817306407708082063016174763530876611806139822690924052796930526831295442464 +/- 2.19e-84]<2.314e-4. (SO5)
```

Equation (SO5) fixes the physical scale and sign of the weighted corrected
selected leading-defect statistic on the augmented roster.  It does not
identify that statistic with the unweighted selector strip, or prove that it
embeds unchanged into the isolated `39695..39894` fold-owned paired-residual block:
the common-profile roster is `39696..40094`, and its completed source
remainder reaches the zero, negative, and outer-positive sectors.  An exact
allocation/embedding identity must be proved before (SO5) may be called a
signed component of `Q_K-T`.

Pi provenance is explicit in (SO1)--(SO4): inverse Airy Fourier
normalization, the equation-(9) Kummer amplitude, odd-square reflection, and
`beta^3=pi*C^2/8`.  No fitted constant is used.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.py
```

No weighted-statistic-to-raw-strip identity, selector-to-paired-residual
embedding, complete fold-owned residual, source/initial-data cancellation,
complete `Q_K-T` or `T_upper`, all-corridor
or height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
