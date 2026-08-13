# Hardy t=1e10 later complete-recurrence stress-cell pilot

Date: 2026-08-09

Status: four_later_stress_cells_close_complete_tail_and_pinned_binary128_rounding_pilot; four-cell finite pilot only, not a proof of the complete recurrence or RH

The maximum-inflation cell (chain 626), narrowest cell (chain 653), and both
longest block-28 cells (chains 3406 and 3423) preserve the complete factored
cubic Legendre remainder at 70 and 110 decimal digits.  The common transformed
child phase is removed at unit modulus before absolute values.

Every parent, reconstructed child, anchored multiplier, child state, and
`tpm=-2*pi` round-to-nearest binary128 preimage lies inside its corresponding
analytic interval.  The emitted source complex multiply-add also lies inside
the explicit `u=2^-113` operation envelope on all four saved calls.

```text
stress cells                                      4
precision overlaps                                4 / 4
maximum local complete tail                       1.63664601327194522170888907502650825178431619309559E-3
maximum cell/point tail inflation                 1.00000000637430190747318927666098581757905045240105E+0
maximum parent rounding-radius/cell-radius ratio  3.23379091813020018243039201214410599062930200464260E-22
maximum stress-augmented output                   4.56639665774875056167646386825368796568692156096651E-4
minimum remaining 0.005 margin                    4.54336033422512494383235361317463120343130784390335E-3
outputs below 0.005                               15 / 15
```

The exact corrected recurrence still uses the covariance
`W69(M,C)=I69-M*C`; the child and multiplier are not charged twice.

## Boundary

The source-derived probe pins round-to-nearest ties-to-even binary128 in the
accepted Docker/compiler environment.  Source `q` and special-function variation over
the cells, the other 1036 later tails, the outer Hardy representation, height uniformity,
`Lambda<=0`, PF-infinity, RH, and a prize-level conclusion remain open.
