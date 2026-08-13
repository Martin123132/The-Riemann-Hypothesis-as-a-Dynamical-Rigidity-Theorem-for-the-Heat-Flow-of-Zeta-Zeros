# Hardy block-20 coefficient-cell joint recurrence identity

Date: 2026-08-09
Status: joint_corrected_recurrence_reduces_exactly_to_cell_endpoint_error_on_all_374_physical_cell_slices; corrected finite model only, not a proof of the source evaluator or RH

## Correlated Identity

Let `P` be the exact parent sum, `I69` the exact equation-(69) saddle family,
`M*C` any transformed-child model, and

```text
W69(M,C)=I69-M*C.
```

Endpoint-complete finite Poisson summation gives `P=I69+Q_exact`.  Therefore

```text
M*C+W69(M,C)+Q_paper = I69+Q_paper,
P-[M*C+W69(M,C)+Q_paper] = Q_exact-Q_paper.
```

This is invariant under arbitrary changes of `M` or `C` when `W69` is changed
with them.  The broad parent, child, multiplier, and saddle Lipschitz columns
must therefore not be added independently in the corrected joint model.

## Cell Guards

On the physical slice `F'(N)=L+fracL`, all 374 occupied
cells preserve `jbot`, `L`, cubic sign, positive curvature, and positive
upper/lower nonsaddle gaps.  The enclosing product boxes are deliberately
wider in their four derived coordinates; their independent coefficient bounds
certify the endpoint majorant.  The 374 saved equation-(69) integrations
overlap across precision, and all 374 independent point Poisson identities
contain zero.

```text
minimum dual lower selector margin      6.05009695398906329292898658110051081293858261500100E-5
minimum dual upper selector margin      1.46352251113189036331552675157089327410956445155411E-3
minimum upper nonsaddle gap             1.41348685639084811636811128023101257876867003904791E-3
minimum lower nonsaddle gap             5.05260792376543342849797417350049128979911950708127E-1
minimum curvature floor                 1.18452670813453492307342823259646530885652350154945E-2
maximum corrected recurrence budget     2.19652337832631474465367310836148849889919479338472E-3
remaining 0.005 slack                   2.80347662167368525534632689163851150110080520661528E-3
outputs below 0.005                     15 / 15
```

The corrected joint recurrence budget is exactly the already certified
coefficient-cell endpoint budget.  No Legendre-tail column is added because a
change in `M*C` is cancelled by the defining exact W1 correction.

## Boundary

The unmodified source does not compute exact `W69` or exact `Q_paper` by
definition.  Its W1 replacement, source-only `t2` term, special functions,
binary normalization, and floating arithmetic must now be transported as a
separate source-to-corrected-model column.  Other blocks, uniform height,
cross-block accumulation, the outer Hardy representation, `Lambda<=0`,
PF-infinity, RH, and a prize-level theorem remain open.
