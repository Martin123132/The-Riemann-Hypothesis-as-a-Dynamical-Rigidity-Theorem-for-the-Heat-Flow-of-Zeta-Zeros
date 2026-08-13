# Weighted opposite-branch Dirichlet bound

Date: 2026-08-13

Status: rigorous ordinary-top-corridor leading-defect selected/full/opposite
branch decomposition; exact finite-integral remainder remains open

For each beta-minus-four mode write

```text
G4_m=G_sel,m+G_opp,m.
```

With the extended weights from Section 11.393,

```text
W_sel=W_full-W_opp,
W_opp=sum_(m=L)^U w_mG_opp,m.                       (OD1)
```

The opposite branch is `C_+` on modes `39696..39893` and `C_-` on modes
`39895..40094`; event zero has weight zero.  With `s=y/Y` and `hY=2pi`,
factoring a unit-modulus initial frequency leaves two ordinary weighted
Fourier polynomials

```text
P_-(s)=sum_(j=0)^197 c^-_j e^(-2pi i j s),
P_+(s)=sum_(j=0)^199 c^+_j e^(-2pi i j s).          (OD2)
```

For a polynomial `P`, midpoint enclosure on `N_p=16384` equal panels and

```text
sup_s |P'(s)|<=2pi sum_j j|c_j|                     (OD3)
```

give a rigorous whole-kernel `L1` bound.  The saved coefficient balls are
uniform over the complete top corridor.  The independently checked results
are

```text
Integral_0^1|P_-(s)|ds
 <[1.5515123137035420334262011544218455272765247680548991287443430851726963762334298e-5 +/- 3.62e-85],
Integral_0^1|P_+(s)|ds
 <[2.2887702786675737135095542432120434763359567528275985179600580434692067737638988e-5 +/- 1.08e-85],
sum <[3.8402825923711157469357553976338890036124815208824976467044011286419031499973287e-5 +/- 4.70e-85]<4.0e-5.                  (OD4)
```

The full-strip branch envelope from Section 11.394 gives

```text
sup_(lambda,y)|C_+|=sup_(lambda,y)|C_-|
 <[0.35516778152330703687618681565668368564924157986032321433037620681354000356616608 +/- 4.39e-81].              (OD5)
```

Applying (OD4)--(OD5) once to each whole branch proves

```text
|W_opp|<[0.0015869141286729416139853418074259188769992638094770124812458908957377742125707888 +/- 1.95e-83]<0.00166.          (OD6)
```

Together with the auxiliary completed full-roster bound,

```text
|W_sel|<=|W_full|+|W_opp|
 <[0.023164147647532708494214429943579549670672752624080384926972046634488144217071402 +/- 5.23e-82]<0.02324.                (OD7)
```

No modewise finite-integral norm is used.  The estimate is conservative and
does not exploit cancellation between `W_full` and `W_opp`; its role is to
make the completion-to-selected-branch passage mathematically valid.  It
still concerns the certified leading stationary-defect weights only, not the
exact finite-integral minus beta-minus-four amplitude remainder.

Pi provenance: the Fourier phases use the exact period `hY=2pi`, while the
branch amplitude and height interval inherit `beta^3=pi C^2/8`.  No new
geometric normalization is introduced.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.py
```

No exact finite-integral amplitude remainder, completed branch source/initial
data splice, complete `Q_K-T` or `T_upper`, all-corridor theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
