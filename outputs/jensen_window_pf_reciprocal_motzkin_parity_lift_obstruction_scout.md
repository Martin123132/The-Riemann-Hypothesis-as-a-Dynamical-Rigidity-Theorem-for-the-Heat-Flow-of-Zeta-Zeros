# Jensen-Window PF Reciprocal Motzkin Parity-Lift Obstruction Scout

Date: 2026-07-06

Status: global parity/sign-lift obstruction diagnostic. This is not a proof against state-space doubled,
modified signed continued-fraction, oscillatory, or production-matrix models,
not Schur positivity, not Jensen-window PF-infinity, not RH, and not
`Lambda <= 0`.

## Purpose

The raw ordinary Motzkin scout showed that the ordinary J-fraction path model
has negative path weights. This note tests a cheap rescue attempt: multiply
all paths of length `m` by a global length-parity sign, or use diagonal sign
conjugation to change edge signs.

The obstruction is that closed excursions of the same length already have
opposite signs:

```text
H0^m > 0
U D H0^(m-2) < 0
```

Both paths start and end at height `0` and have the same length. So a global
length-parity sign cannot make both path weights nonnegative. A diagonal sign
conjugation also cannot change the sign of any closed excursion, because the
endpoint factors telescope.

## Outputs

Machine-readable outputs:

```text
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_lamgrid_n0_n20_d2_d8_m2_m8_dps520.jsonl
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.py
```

Current result:

```text
validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues
```

## Finite Grid

The Arb grid checks:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n=0..20
degrees d=2..8
path lengths m=2..8
dps=520
```

It validates:

```text
5,145 / 5,145 rows have same-length mixed-sign witnesses
positive witness: H0^m
negative witness: U D H0^(m-2)
```

This sharpens:

```text
rp_04_companion_or_production_matrix_total_positivity
rp_09_signed_or_modified_continued_fraction
```

and is finite evidence only.

## Interpretation

The ordinary Motzkin model cannot be repaired by a sign depending only on path
length. The obstruction is not numerical delicacy; it is structural:
same-length closed excursions already have both signs. The same diagonal sign conjugation
cannot repair it either because closed excursion signs are
invariant.

This does not reject a truly modified route. A surviving path proof would need
to change the state space, the weights, or the representation itself, for
example through a state-space doubled model, a new positive network, or an
Xi/Phi-specific path formula.

## Boundary

Passing this checker rejects only global length-parity signs and diagonal sign
conjugation for the raw ordinary Motzkin model. It does not rule out
state-space doubled models, modified signed continued fractions, oscillatory
matrix theorems, production matrices with different weights, Schur positivity,
Jensen-window PF-infinity, or `Lambda <= 0`.
