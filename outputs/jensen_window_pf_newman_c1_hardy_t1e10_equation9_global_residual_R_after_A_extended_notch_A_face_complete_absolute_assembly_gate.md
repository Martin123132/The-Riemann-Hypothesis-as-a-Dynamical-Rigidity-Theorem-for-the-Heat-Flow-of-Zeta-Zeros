# Complete translated A-face absolute assembly

Date: 2026-08-23

Status: rigorous saved-height local assembly; not a proof of `R_after_A` or RH

The exact handoff `x_16` partitions the translated 42-mode A-face channel as

```text
I_(A,42)=I_(A,42)^core+I_(A,42)^end,

core: 0<x<=x_16,       endpoint: x_16<=x<=1/2.      (AA1)
```

Both dependency gates use the same equation-(9) physical normalization and
their certified `x_16` balls overlap.  Their absolute bounds are

```text
|I_(A,42)^core| <= [1.0173490686775675383189830388201349145801998842004236279021692900822170000000000e-5 +/- 4.10e-75],
|I_(A,42)^end|  <= [0.0036015988499715977974756588142866309847034217693321195358920676066335130000000000 +/- 1.80e-5].  (AA2)
```

Therefore the triangle inequality, with no oscillatory cancellation, gives

```text
|I_(A,42)| <= [0.0036117723406583734728588486446748323338492237681741237721710892995343351700000000 +/- 1.80e-5]
             < 0.00364.                (AA3)
```

Against the inherited joined target `R_after_A<0.0368147039947`, this
leaves a certified allowance greater than `0.03317` for the other post-A
channels.  That subtraction is budget bookkeeping only; it is not a proof
that those channels fit the allowance.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_complete_absolute_assembly_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_complete_absolute_assembly_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_complete_absolute_assembly_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_complete_absolute_assembly_gate.py
```

Pi provenance: no new `pi` is introduced in (AA1)--(AA3).  Both summands
inherit their occurrences from the canonical Fresnel hierarchy and the same
equation-(9) normalization.

Proof boundary: (AA3) proves only the complete translated 42-mode A-face
channel at `t=10^10`.  The zero sector, lower edge, B-owned terms, negative
bulk, and every other post-A remainder remain outside this assembly.  No
joined `R_after_A`, `R_Dir`, `Q_K-T`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
