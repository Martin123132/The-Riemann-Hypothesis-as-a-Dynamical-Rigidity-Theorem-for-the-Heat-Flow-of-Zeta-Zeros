# Jensen-Window PF Monotone-Contraction Frontier Scout

Date: 2026-07-06

Status: symbolic frontier scout. This is not a proof of Jensen-window
PF-infinity, all-shape Schur positivity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_monotone_contraction_frontier_scout`.

Proof boundary: this artifact gives exact sufficient positivity certificates
for the two first hard column-frontier polynomials under an extra
ratio-contraction hypothesis. It does not prove that the zeta heat-flow
coefficients satisfy that hypothesis for all shifts or lambda values, does not
construct a Cauchy-Binet kernel identity, and does not close `jwpf_06`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_frontier_scout.json
```

Builder:

```text
python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_frontier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_frontier_scout.py
```

Current result:

```text
validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues
```

## Question

The Cauchy-Binet cone frontier matrix left two live hard rows:

```text
d3_column_recurrence_m8
d4_column_recurrence_m6
```

Adjacent log-concavity gives only the ratio-contraction box
`0 <= x_i <= 1`, where

```text
x_i = (A_{n+i+1}/A_{n+i}) / (A_{n+i}/A_{n+i-1})
```

The exact rational countermodel lies in that box, so the box is too weak. This
scout tests the sharper region:

```text
x1 <= x2 <= x3
```

truncated to the contractions available in the degree.

## Exact Certificates

For `d3_column_recurrence_m8`, substituting

```text
x1 = s
x2 = s + u*(1-s)
0 <= s,u <= 1
```

puts the normalized hard polynomial on the monotone-contraction square. Its
Bernstein certificate has multidegree `[7, 2]`, 24 coefficients, and
Bernstein minimum 45.

For `d4_column_recurrence_m6`, substituting

```text
x1 = s
x2 = s + u*(1-s)
x3 = x2 + v*(1-x2)
0 <= s,u,v <= 1
```

puts the normalized hard polynomial on the monotone-contraction cube. Its
Bernstein certificate has multidegree `[7, 3, 1]`, 64 coefficients, and
Bernstein minimum 84.

These are exact symbolic sufficient conditions for the two listed frontier
rows only.

## Countermodel Split

The exact rational countermodel has ratio contractions:

```text
x1 = 65/66
x2 = 44/65
x3 = 50/77
```

It violates monotone contractions by `x2 < x1` and `x3 < x2`. Its normalized
hard-row values remain negative:

```text
d3_m8 = -122505653/6324912
d4_m6 = -19180303/754677
```

So the monotone-contraction condition escapes the known rational
log-concavity countermodel without assuming endpoint PF, Laguerre-Polya
membership, RH, or `Lambda <= 0`.

## Finite Zeta Diagnostics

Using the existing Arb coefficient enclosures on five lambda values and shifts
`n=0..20`, the checked zeta windows satisfy the sharper condition on the two
frontier degrees:

```text
d=3: 105 monotone-contraction rows, 105 positive hard-value rows
d=4: 105 monotone-contraction rows, 105 positive hard-value rows
```

Across these 210 finite zeta rows, the sampled ratio contractions stay close
to one. The smallest sampled monotone gap is positive:

```text
d=3 min gap sample: 5.678678882796781911E-4
d=4 min gap sample: 5.293868421064408833E-4
```

The smallest sampled normalized hard values are also positive:

```text
d=3 min normalized hard value sample: 5.700813557379444575E+1
d=4 min normalized hard value sample: 9.811781725060354254E+1
```

This is finite diagnostic evidence only.

## Consequence

The live column-frontier route has sharpened from "find any determinant
integral" to a more concrete theorem target:

```text
Prove a noncircular zeta-specific monotone-contraction theorem, or a stronger
positive kernel identity, that implies x1 <= x2 <= x3 for the required
Jensen-window coefficient ratios.
```

That would cover the first hard column-frontier rows. It would still need
extension to all `m,d,n`, all relevant lambda values, and eventually the
all-shape Toeplitz/Jacobi-Trudi contract.

## Links

```text
outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_log_concavity_frontier_scout.md
work/rh_compute/results/jensen_window_pf_cauchy_binet_cone_frontier_matrix.json
work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json
```
