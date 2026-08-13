# Leading-defect weighted Airy Green packets

Date: 2026-08-13

Status: rigorous top-corridor weighted packet interface for the certified
leading stationary defects; no finite-integral amplitude remainder or
completed endpoint splice

Retain the exact leading-defect coefficients from Section 11.385,

```text
a_m=[R_A(s_m,mu)exp(i beta^3F(s_m,mu))-1]/sqrt(m).    (WG1)
```

The saved coherent/leakage balls reconstruct every coefficient on the two
selected branches, with no new midpoint evaluation.  For a contiguous block
`I=[a,b]`, define

```text
V_I=min(|a_a|,|a_b|)+sum_(j=a)^(b-1)|a_(j+1)-a_j|.  (WG2)
```

Finite interval enumeration of every branch block proves

```text
minus branch: V_I <7.705213812399580233921878025249441179767018184065818786621093750000000000000000000000000000000000000e-6,
plus branch:  V_I <7.765296993613448791324901238075639753333234693855047225952148437500000000000000000000000000000000000e-6. (WG3)
```

If `M_j=M(-r_j)`, its positive monotonicity gives, in either orientation,

```text
V_I(aM)/M_max <= V_I(a)+max_(j in I)|a_j|.            (WG4)
```

The final term pays for the entire monotone modulus variation; no lower
bound for `M` and no constant-amplitude replacement is used.  Exhausting all
blocks yields

```text
minus effective variation
 <1.540587155942246786756964771480937770320451818406581878662109375000000000000000000000000000000000000e-5,
plus effective variation
 <1.552592025585636488965507115422859385489573469385504722595214843750000000000000000000000000000000000e-5,
two-branch sum
 <3.093179181527883275722471886903797155810025287792086601257324218750000000000000000000000000000000000e-5.               (WG5)
```

Applying (WG4) to the exact phase packets in Section 11.390 proves, uniformly
over every active contiguous block and the complete top corridor,

```text
sum of the two branchwise weighted K-packet bounds
 <[7.577214657378402093616767903121947819989910359723684903721738138996948948202804065567310294239631462e-7 +/- 4.89e-107]<7.58e-7,

sum of the two branchwise weighted partial_s K-packet bounds
 <[0.001632452081577785218015061517662004244922091890171106705168339022546472367782867117173874291792704844 +/- 1.95e-103]<0.001634.   (WG6)
```

No termwise Green-kernel absolute sum occurs in (WG6).  The first number is
small, but it is not yet an additive residual budget: it must be integrated
against the correctly normalized completed source.  The second number must
still be multiplied by the composed initial-data channel and is not itself
small enough to discard.  Moreover, `a_m` is only the certified leading
stationary defect; the exact finite-integral amplitude and endpoint/Fresnel
remainder is not represented by (WG1).

Pi provenance is inherited from `beta^3=pi C^2/8`, the exact stationary
action bridge, and the Airy Wronskian packet in Section 11.390.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate.py
```

No completed source integral, initial-data cancellation, exact finite-
integral selected/logistic match, grouped 398-mode splice, complete
`Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
