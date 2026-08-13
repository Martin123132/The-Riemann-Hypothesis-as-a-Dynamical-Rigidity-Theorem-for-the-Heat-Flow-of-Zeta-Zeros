# Exact endpoint residual versus source telemetry: nonidentification guard

Date: 2026-08-13
Status: exact guard certificate; not a proof of the endpoint residual bound

Use three different upper objects:

```text
Q_K = exact finite equation-(9) Kummer roster,
Q_P = published asymptotic/generalized-Gauss hybrid,
Q_S = implemented source hybrid,
T   = exact classical upper main.
```

They obey the exact bookkeeping identity

```text
E_source=Q_S-T
        =(Q_K-T)+(Q_P-Q_K)+(Q_S-Q_P).                 (NI1)
```

The symmetric endpoint calculation and the double-Weber involution reconstruct
`Q_K`, not `Q_P` or `Q_S`.  The paper audit independently records that the
exact Kummer integral precedes the asymptotic grouping, that the sine member
of each equation-(10) pair is dropped before equation (51), and that the
published hybrid is not exact.  Consequently the saved source discrepancy

```text
Q_S-T = [-0.0069830180084699455319725803829133704069288214989573701554341205283354305539931196430885535174224228367531020376 +/- 9.47e-96]               (NI2)
```

cannot be transferred to the exact endpoint residual.  Its absolute size is
more than `[811.9788381941796497942198107680223176374 +/- 4.91e-38]` times the
reference `8.6e-6` scale, but that fact belongs only to `Q_S-T`.

The Gamma calculation does close a different bridge.  With `G` denoting the
projected exact Gamma bulk,

```text
R_end^Gamma=Q_K-G=(Q_K-T)-(G-T),                      (NI3)
|G-T| <= [2.1850771683060808893917330959836276744161995433604234851443178245796785783235933093415905860820286790000000000e-19 +/- 1.59e-88].
```

Thus the Gamma replacement changes the exact Kummer-to-classical residual by
less than `2.186e-19`; it does not identify that residual with source telemetry.
The missing quantities in (NI1) are the analytic Kummer-to-published defect
`Q_P-Q_K` and the complete published-to-implemented bridge `Q_S-Q_P` in the
same global normalization.  Partial evaluator certificates must not be
silently substituted for either complete quantity.

Pi provenance: `pi` comes from the equation-(9) Kummer/Fourier phase, the
Riemann-Siegel theta normalization, and the common `pi/8` reflection.  No
fitted or geometric parameter is introduced.

Proof boundary: this gate rejects a false equality and certifies the tiny
Gamma-to-classical carrier defect at `t=10^10`.  It proves no numerical bound
for `Q_K-T`, no complete B face estimate, no A-fold splice, no complete
source-hybrid error theorem, no height-uniform theorem, no `Lambda<=0`, no
PF-infinity, no RH, and no prize-level conclusion.
