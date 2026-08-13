# Newman C1 Hardy Modewise Derivative Handoff Gate

Date: 2026-08-05

Status: `exact_conditional_handoff_with_open_external_error_theorem`; not a proof of the external error theorem, RH, or a prize-level conclusion.

## Source audit

Primary paper: [2607.15310v1](https://arxiv.org/abs/2607.15310) by David Lewis, Ashley Brereton.

The paper's equations (120)--(123) reduce the recursive MGS error to the largest local iteration error.  On pages 41--42 it then states that a precise upper bound would require an estimate of that term and does not formulate one for the higher-order sums.  Page 44 separately says the hybrid Hardy representation is not exact.  Table 4 is sample evidence; its relative-error statistics omit `|Z(t)|<=1/2`, and the surrounding text says adverse combination of cubic-sum errors is not guaranteed away.

The pinned fifteen-value source is more structurally useful than an independent fifteen-run model: it shares central hierarchy parameters and varies the nearby heights and phases on one `0.01` grid.  But its own comment says there is no hard-and-fast rule behind the multi-value accuracy claim.  Source derivative hash: `0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d`.

## Exact handoff

Put `delta_j=j/100`, `-7<=j<=7`, and let `q_k` be the normalized integer Gram polynomial of degree `k` listed in the artifact.  Assume on the whole shifted interval

```text
sup_(|u-T|<=7/100) |epsilon^(r)(u)| <= D_r(T),  0<=r<=6.
```

Taylor's theorem and exact discrete orthogonality give

```text
|<epsilon_T,q_k>|
 <= sum_(r=k)^5 w_(k,r) D_r(T) + w_(k,6,rem) D_6(T),

||(I-P_5)epsilon_T||_2
 <= [sqrt(sum_j |delta_j|^12)/6!] D_6(T).
```

The weights are exact algebraic numbers `rational/sqrt(integer)`; no fitted decay law enters.  The rough-tail multiplier is `2.5061285157281298155644476165324990346357E-10` and its square is exactly `4069880729/64800000000000000000000000000`.

| mode | nonzero center-derivative weights | sixth-order remainder weight |
|---:|---|---:|
| 0 | D_0: 3.8729833462074168851792653997823996108329E+0, D_2: 3.6147844564602557595006477064635729701107E-3, D_4: 1.0061150070481045197276802782990278100142E-6 | 1.3255644111963977568867256133107094371743E-10 |
| 1 | D_1: 1.6733200530681510959563440515703749787856E-1, D_3: 9.3148149620793744341569818870750873819067E-5, D_5: 1.8408512631431886044919732595906946641617E-8 | 1.9925523343034861439293488000755231830717E-10 |
| 2 | D_2: 3.2114378918692065113922225584599816582644E-3, D_4: 1.2654594550103658991319353176788737248637E-6 | 1.9441251075863297256468392404753386342106E-10 |
| 3 | D_3: 3.9889848332627187365556697850526262676024E-5, D_5: 1.2077759634045453952349111293631562865796E-8 | 1.6914164195820790827515101713269090495067E-10 |
| 4 | D_4: 3.6327478776894416977547279546939807080724E-7 | 1.6870563341594626378434719958639512299356E-10 |
| 5 | D_5: 2.5816814995423748202781130360123577861209E-9 | 1.5988237015914757154807206924994209230318E-10 |

A stronger alternative is one fixed-branch holomorphic disk certificate: if `|epsilon(z)|<=M_rho(T)` on `|z-T|<=rho` with `rho>0.07`, Cauchy's estimate supplies `D_r<=r! M_rho/rho^r`.  If any combinatorial or analytic transition crosses the physical interval, the interval must instead be partitioned and its jump terms entered as additional exact mode vectors.

## Result

The scalar two-height route remains retired.  The missing input is now exact: prove one branch-stable `C^6` envelope, one stronger complex-disk envelope, or a transition-aware piecewise substitute for the external fifteen-value error. Only then may the saved physical coefficient vectors be applied.

This gate proves exact finite derivative-to-mode algebra and audits the published/source-code error claims. It does not prove any evaluator derivative bound, complex-disk bound, transition exclusion, external absolute error constant, physical carrier value, interval enclosure, retained observation, determinant sign or bound, complete-current inequality, all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

```text
built Newman C1 Hardy modewise derivative handoff gate: 12 rows, 6 exact low modes, 7 derivative orders, 1 exact rough-tail constant, 4 primary-paper locations, 7 source markers, 2 conditional routes and 1 open external error theorem
```
