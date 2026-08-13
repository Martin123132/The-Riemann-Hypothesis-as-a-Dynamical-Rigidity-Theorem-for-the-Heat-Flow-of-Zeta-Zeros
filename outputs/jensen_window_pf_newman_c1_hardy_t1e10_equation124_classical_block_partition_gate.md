# Equation-(124) classical block-partition atlas

Date: 2026-08-09

Status: finite interval diagnostic validated; not a proof of RH

The source checkpoint journal preserves cumulative `ZP` values after block
zero and every block through 35.  Paper equation (124) supplies the coordinate

```text
alpha(n)=2n+t/(pi*n).
```

The executable fixes the discrete convention: it assigns `raecutoff=RN1`,
evaluates the lower root `CE` in equation (124), sets `NC=AINT(CE)`, and sums
the direct Riemann--Siegel part over `n=1,...,NC`.  Because `CE>0`, `AINT` is
the floor operation here, so the complementary classical tail starts at
`N_k=NC+1`.  No `RN1+1` cell-edge shift occurs in the source.

For each saved alpha endpoint `RN1_k`, the gate certifies the unique integer
`N_k` satisfying

```text
alpha(N_k)<RN1_k<alpha(N_k-1).
```

It then defines an exact classical comparison target, independently of the
non-exact hybrid formula,

```text
T_k(t)=2 sum_(n=N_k)^39894 cos(theta(t)-t log n)/sqrt(n).
```

The final cutoff is `N_35=622`, so `T_35` is exactly the full classical upper
main sum used in Section 11.298.  The 36 cutoffs partition all 39,273 terms
without overlap or omission.  Every cutoff inequality and all
540 cumulative comparisons overlap at 70 and 110 decimal
digits; an independent checker repeats them at higher precision.

Block zero is already discrepant on all fifteen outputs.  Its signed source
minus exact-target residual is strictly negative and exceeds `0.005` in
absolute value on all fifteen:

```text
minimum |E_0| = [0.016120329380248886474667277959726967994504125015726608953449790204645683872187886845711332553940626527420734472 +/- 1.06e-98]
maximum |E_0| = [0.016159333143851197781251020273789022130732088096949811669214527227749982595334811872695784053702997072879715180 +/- 1.06e-98]
```

The net error added by blocks 1--35 is positive on all fifteen outputs, so it
partially cancels block zero rather than creating the final discrepancy.  The
saved final residual still matches the independently certified classical
upper-main residual on every output.

This atlas identifies the earliest block-aligned defect but not its analytic
cause.  Block zero combines the transition value, direct `ter` terms, and the
final `sqrt(8/a)` normalization.  Those columns must be separated against the
exact `n=N_0,...,39894` target.  Equation (124) is used only to define exact
integer boundaries; the paper explicitly says the hybrid representation is
not exact.  The atlas is finite at the saved heights and has no RH implication.
