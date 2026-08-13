# Weighted completion partial-sum obligation

Date: 2026-08-13

Status: exact cancellation-preserving Abel reduction and non-identifiability
guard; uniform completed partial-remainder bound remains open

Let `L=39696`, `q=39894`, and `U=40094`.  Extend the 398 selected
leading-defect weights to an auxiliary complete 399-mode beta-minus-four
roster by

```text
w_m=a_m (m!=q),             w_q=0.                  (WC1)
```

For the unweighted beta-minus-four terms define

```text
S_k=sum_(m=L)^k G4_m,
R_k=H4+C4+sum_(m=k+1)^U G4_m.                       (WC2)
```

The completed projection `H4+sum G4+C4=g4_lambda(0)` gives, without
separating `H4` or `C4`,

```text
S_k=g4_lambda(0)-R_k.                               (WC3)
```

Finite Abel summation is exactly

```text
W_full=sum_(m=L)^U w_m G4_m
 =w_U S_U+sum_(k=L)^(U-1)(w_k-w_(k+1))S_k           (WC4)

 =w_L g4_lambda(0)-w_U R_U
  -sum_(k=L)^(U-1)(w_k-w_(k+1))R_k.                 (WC5)
```

An exact rational 399-mode calculation checks both equalities.  Each `R_k`
keeps the endpoint half-current, the zero/negative/outer-positive complement,
and the untouched upper suffix together.  Thus (WC5) does not split bare
endpoint exponentials or take termwise absolute values of the `G4_m`.

The certified weight budget is

```text
sum_(k=L)^(U-1)|w_(k+1)-w_k|
 <1.547051080601302902524677926332508093310025287792086601257324218750000000000000000000000000000000000e-5,

|w_L|+sum|Delta w|+|w_U|
 <3.093179181527883275722471886903797155810025287792086601257324218750000000000000000000000000000000000e-5<3.093180e-5.          (WC6)
```

Consequently a uniform completed-family estimate

```text
max(|g4_lambda(0)|, max_(L<=k<=U)|R_k|) <= B_comp
```

would imply the cancellation-preserving bound

```text
|W_full| < 3.093180e-5 B_comp.                      (WC7)
```

This is not yet the selected-branch weighted sum.  Since
`G4_m=G_sel,m+G_opp,m`, exact algebra requires

```text
W_sel=W_full-W_opp,
W_opp=sum_(m=L)^U w_m G_opp,m.                      (WC7a)
```

The completed identity controls `W_full`; a separate cancellation-preserving
bound for `W_opp` is required before (WC7) can control `W_sel`.

The one global completion is not enough by itself.  The last two selected
weights satisfy

```text
|w_40093-w_40094|
 =[1.917549273811047783055983018130064010620117187500000000000000000000000000000000000000000000000000000e-7 +/- 4.06e-11]>0.                 (WC8)
```

Holding event zero fixed, replace `G4_40093` by `G4_40093+z` and
`G4_40094` by `G4_40094-z`.  The completed total and every other term are
unchanged, while `W_full` changes by `(w_40093-w_40094)z`.  Therefore the global
identity alone cannot identify or cancel the auxiliary weighted full roster.  This is a
logical non-identifiability result for that identity, not a claim that the
actual Fourier coefficients can be varied arbitrarily.

The next theorem targets are a uniform representation or bound for the
completed family `R_k` and a separate weighted opposite-branch estimate;
only their combination through (WC7a) reaches the selected branch.

Pi provenance: this gate introduces no new `pi`; all weights and completed
terms inherit `beta^3=pi C^2/8` and the prior beta-minus-four normalization.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.py
```

No uniform completed partial-remainder or weighted opposite-branch bound,
selected-branch initial-data/source cancellation, exact finite-integral amplitude remainder, grouped 398-mode
splice, complete `Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
