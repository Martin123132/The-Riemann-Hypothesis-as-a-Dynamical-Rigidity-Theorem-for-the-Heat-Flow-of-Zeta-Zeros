# Hardy Fastcode Diagnostic Adapter

This directory contains only the RH-owned build and validation wrapper for the
external Hardy-function implementation. The third-party source remains in the
separate local checkout:

```text
<local-third-party-root>/Hardy-function-fastcodes
```

Pinned upstream state:

```text
repository: https://github.com/dml2391/Hardy-function-fastcodes
commit: 2e16dac3206b707052c3ac4cacdf3d1a2325e636
license: GPL-3.0
```

The upstream program evaluates the ordinary Hardy function and documents an
asymptotic relative-error floor. It is a candidate numerical diagnostic, not a
proof oracle. No output from it may be described as an interval enclosure or as
an `h^2` theorem without a separate validated error analysis.

The intended bridge keeps the upstream source unmodified. Ordinary Hardy
samples at nearby heights are first stripped of their Riemann--Siegel endpoint
correction. A certified real exponential approximation on the fixed logarithmic
interval then reconstructs the heat-weighted polynomial carrier observations.
The approximation error, Hardy-code error, endpoint error, and normalization
error must remain separate in every result.

The pinned source declares the MPI reduction operator `mysum` without assigning
it. OpenMPI therefore rejects the first `MPI_Allreduce`. The local build applies
`patches/0001-initialize-mpi-sum.patch`, which adds only `mysum = MPI_SUM` after
MPI initialization.

The pinned source also executes `STOP` before its later `MPI_Finalize` call,
making finalization unreachable and causing OpenMPI to return a failure after a
numerically completed run. The local build therefore also applies
`patches/0002-finalize-mpi-before-stop.patch`, which moves the existing finalize
call immediately before `STOP`. The upstream checkout remains byte-for-byte
unchanged.

The same unreachable-finalize defect occurs in the 15-value
`zeta14cubicmult.f90` program. The multi-value diagnostic build applies the
parallel minimal patch
`patches/0003-finalize-multi-mpi-before-stop.patch`. It does not alter the
upstream shift grid, asymptotic expansion, arithmetic, or output format.

The pinned multi-value source additionally tries to register an undefined
callback `qsum`, although every reduction actually uses `MPI_SUM`, and passes
the eight-byte `nchalf` value to an OpenMPI interface requiring a default
integer count. `patches/0004-fix-multi-openmpi-interface.patch` removes only
that dead registration and converts the two reduction counts with
`int(nchalf)`. These are build-interface repairs, not arithmetic changes.

For conditioning experiments only, a separately named build may also apply
`patches/0005-multi-shift-step-004.patch`. It replaces the hard-coded `0.01`
spacing by the declared constant `shift_step=0.04` in all four phase/height
updates and in output labels. This is an intentional sample-grid variant and
must pass independent direct-sum calibration before use. It is not part of the
minimal upstream portability repair.

`patches/0006-multi-shift-step-002.patch` is the corresponding intermediate
`shift_step=0.02` experiment, applied instead of patch 0005. Each spacing
variant is a distinct executable and output provenance chain; neither may be
silently substituted for the upstream `0.01` grid.
