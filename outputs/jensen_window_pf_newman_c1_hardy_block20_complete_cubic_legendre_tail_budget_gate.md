# Hardy block-20 complete cubic Legendre-tail budget

Date: 2026-08-09
Status: complete_saved_point_cubic_legendre_tail_and_endpoint_budget_close; finite saved-point theorem only, not a proof of the complete evaluator or RH

## Exact Joint Remainder

For `u=k-a1` and `z=12*a3*u/y^2`, the exact cubic stationary map gives

```text
F''(X)=y*sqrt(1+z),
A(z)=(1+z)^(1/4),
D2(z)=1+z/4-3*z^2/32.
```

Thus the omitted transformed-child term has both a phase part `H-H3` and an
amplitude part `1/A-1/D2`.  The gate evaluates the exact radicals rather than
stopping at the first omitted Taylor coefficient.  It preserves the source's
binary `tpm` and fixed integer/parity/orientation normalization, so no exact
`2*pi` periodicity is assumed.

## Certified Roster

```text
recursive calls                         374
saved child indices                     753
precision overlaps                      374 / 374
maximum |z|                             1.29757399980552990362988882312330031833776619363891515685099214115505336429660659E-3
minimum 1+z                             9.98702426000194470096370111176876699681662233806361084843149007858844946635703393E-1
maximum |H-H3|                          8.42281119291988804337324269552774205384748176632539682429059080783897176805870467E-6
maximum local recurrence tail           4.36531639513907282303320676432780228539195403522825722467785282518843721344812488E-4
maximum transported tail                1.50953619023699380186406155621874321071796894326424470014008956400111392130294699E-4
maximum endpoint plus tail              2.34747699735001412484007926398336281997099168771114447001400895640011139213029470E-3
outputs below 0.005                     15 / 15
```

The endpoint and Legendre-tail columns remain separately recorded.  The
combined number uses absolute call-by-call transport and does not use the
signed output diagnostic.

## Boundary

This closes the all-orders cubic Legendre phase-and-amplitude truncation only
at the exact saved block-20 child indices.  Tail transport over coefficient
subcells, source arithmetic, other blocks, cross-block accumulation, the outer
Hardy representation, `Lambda<=0`, PF-infinity, RH, and a prize-level theorem
remain open.
