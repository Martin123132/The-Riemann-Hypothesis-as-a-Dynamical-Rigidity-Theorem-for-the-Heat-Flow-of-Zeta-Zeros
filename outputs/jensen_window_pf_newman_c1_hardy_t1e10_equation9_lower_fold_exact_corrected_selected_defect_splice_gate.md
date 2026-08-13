# Exact corrected selected-defect splice

Date: 2026-08-13

Status: exact cancellation-preserving splice and a rigorous signed
top-corridor bound; the final source/initial-data orientation remains open

For every mode in the augmented roster, the beta-minus-four branch split and
the exact finite-height correction are

```text
G4_m=G_sel,m+G_opp,m,
Delta G_m=G_ex,m-G4_m.                               (CS1)
```

Therefore exact algebra, without assigning an exact branch decomposition to
`G_ex,m`, gives

```text
G_ex,m-G_opp,m=G_sel,m+Delta G_m,
W_corr=sum_m w_m(G_ex,m-G_opp,m)
      =W_sel,4+Delta W_finite-t.                     (CS2)
```

The opposite beta-minus-four branch remains inside the completed remainder;
the new Kummer-ODE correction is inserted as one whole Fourier profile.  No
modewise finite-integral bound is used.

Combining the two certified estimates gives

```text
sup|W_corr|<[4.0888126662417138234623568368750784037941459236637208937888077887066625567273883e-5 +/- 1.34e-85]<4.088814e-5. (CS3)
```

The exact correction is less than
`[2.9340389942536974707412913447163096342286880480778682647849746168661723664329505e-7 +/- 8.84e-88]` times the previous selected bound.
More importantly, all five grouped real parts are negative.  Appending the
quadrature, height-transport, and exact finite-height errors in the adverse
direction proves throughout the complete ordinary top corridor

```text
Re W_corr<[-2.2332236499234579670996191356006495175340632503244015920485121238255035796888877e-5 +/- 3.48e-85]<-2.233e-5. (CS4)
```

Thus the old scalar-headroom comparison by modulus is not the decisive test:
the corrected term has a certified sign.  Whether that sign is favorable
depends on the still-unproved common source/initial-data orientation in the
final `Q_K-T` or `T_upper` identity.  This gate does not assume that
orientation.

Pi provenance is inherited unchanged from `beta^3=pi*C^2/8`, `hY=2pi`, the
inverse Airy Fourier normalization, and the Kummer-ODE profile theorem.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate.py
```

No source/initial-data sign orientation, local-headroom closure, complete
`Q_K-T` or `T_upper`, all-corridor or height-uniform theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
