# Hardy Chain Telemetry Derivative

Status: read-only finite diagnostic; not an interval evaluator or RH proof.

## Source Boundary

`zeta14cubicmult_telemetry.f90` is a copied derivative of the accepted
resumable evaluator. The accepted source remains at:

`work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90`

Its pinned SHA-256 is:

`0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d`

Every derivative-only line is enclosed in a named `RH_TELEMETRY` block.
`check_hardy_chain_telemetry_source_contract.py` verifies the accepted hash,
the exact body of all observer blocks, and exact equality of all nonblank
source lines outside those blocks.

The observer is inert unless both of these variables are set:

- `RH_HARDY_CHAIN_TELEMETRY`: append-only JSONL output path.
- `RH_HARDY_TELEMETRY_MAX_CHAINS`: positive bounded chain count.

`RH_HARDY_TELEMETRY_MAX_DIRECT_TERMS` optionally caps independent direct
summation. Each JSONL record is flushed before the evaluator continues.

## Recorded Quantities

For each admitted chain the observer records:

- block, sum index, conjugate branch, `ip`, `MIT`, and the actual `tpm`;
- the original generalized-sum coefficients passed into `computarrs`;
- every transformed level's inclusive length, degree, coefficients, `xr`,
  `fracL`, conjugation flag, and subtraction flag;
- Fortran `t1` through `t5`, endpoint term, and `qq`;
- recurrence multiplier and source state before and after every transform;
- direct finite parent and child sums, adapted sums, and local defects.

Every proof-relevant real scalar is recorded twice: as a readable 40-place
scientific decimal and as its exact 32-hex-digit `real(kind=16)` payload. The
hex payload lets a ball-arithmetic checker reconstruct the evaluator input as
an exact dyadic rational, without assuming that decimal formatting was
correctly rounded.

The direct diagnostic uses the negative phase convention appearing in `q`
and the kernel. It then applies the source's level conjugation and subtraction
decisions before comparison. The denominator-weighted kernel remains a
separate adapter error; it is not silently identified with an exact finite
child sum.

## Phase Constant Provenance

No unexplained value of pi is inserted. The accepted evaluator defines:

```fortran
C=1.0
p=4*ATAN(C)
tpp=2*p
tpm=-tpp
```

The telemetry chain header records that computed `tpm`. The independent
Python reconstruction consumes the recorded value rather than substituting a
separate library constant.

## Bounded Fixture

The one-worker fixture command is:

```text
python work/rh_compute/scripts/run_hardy_chain_telemetry_fixture.py --max-chains 64
```

It compiles a separate binary with one Docker CPU and below-normal priority,
runs telemetry-disabled and telemetry-enabled copies, and requires:

- exact enabled/disabled full-precision checkpoint journals;
- exact accepted-reference checkpoint state after removing only `run_id`;
- exact enabled/disabled and accepted-reference displayed values;
- valid bounded telemetry and exact logged `qq` decomposition.

The route analyzer then re-sums every logged finite level at 50, 80, and 120
decimal digits and compares termwise W1, direct-local-defect, and exact-shell
routes:

```text
python work/rh_compute/scripts/analyze_hardy_chain_telemetry_routes.py
```

## Proof Boundary

This machinery can establish non-interference and finite pointwise identities
on the saved low-height chains. It does not establish interval coefficient
control, a transformed-level adapter theorem, special-function remainders,
selector or hierarchy-transition margins, the outer Hardy representation
bound, physical-height complexity, Lambda <= 0, PF-infinity, RH, or a
prize-level conclusion.
