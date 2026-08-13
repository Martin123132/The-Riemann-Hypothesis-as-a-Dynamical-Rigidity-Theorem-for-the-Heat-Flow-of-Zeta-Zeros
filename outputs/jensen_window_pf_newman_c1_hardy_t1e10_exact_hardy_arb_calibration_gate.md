# Exact Hardy-Z Arb calibration at t=10^10

Date: 2026-08-09

Status: finite interval certificate validated; not a proof of RH

The fifteen preserved binary128 source outputs at
`t=10^10-0.07,...,10^10+0.07` were compared directly with

```text
Z(t)=exp(i theta(t)) zeta(1/2+i t),
theta(t)=Im log Gamma(1/4+i t/2)-t log(pi)/2.
```

Arb evaluated each target independently at 70 and
110 decimal digits.  Every pair of Hardy enclosures overlaps,
every imaginary enclosure contains zero, and every signed
`source-Hardy` interval is strictly negative.

All 15 source outputs are rigorously farther than `0.005`
from the exact Hardy target.  The smallest certified absolute gap occurs at
output 1 and is enclosed by

```text
[0.0066363141932320866944312171077425062526848171437596925364837220689307200591473902958222704059785270817580938599 +/- 9.33e-96]
```

The largest occurs at output
15 and is enclosed by

```text
[0.0072473879139780114569208124954228834383962645695155431663441329268340541446785973841202840027861237261124131350 +/- 8.70e-96]
```

This falsifies the absolute `0.005` accuracy target for this saved legacy
run.  It does not falsify RH, and it does not prove that every implementation
or every parameter choice in the paper has the same error.

The existing blocks-20--35 internal corrected-model column is too small to
explain the discrepancy by itself: even an adjustment with its most helpful
possible sign leaves a reverse-triangle residual greater than `0.005` on all
fifteen outputs.  At least one still-open outer column from Section 11.295 is
therefore essential.  The next proof-facing task is to replace the non-exact
hybrid outer formula by an exact identity with explicit remainders; fitting or
retuning the source output would not close that theorem obligation.
